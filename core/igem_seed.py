"""iGEM Registry parçalarını katalog kaydına çevirir.

Veri `data/igem_parts.json` içinde durur ve `tools/fetch_igem_parts.py` ile
çekilir. Registry kayıtları Creative Commons Attribution 4.0 (CC BY 4.0)
ile dağıtılır; betik her parçanın lisansını API'den okuyup kayda yazar, bu
modül de cc-by olmayan bir kayıt görürse onu atlar. Her parçanın
`accession` alanı kaynağı ve lisansı taşır, atıf böylece ekranda görünür.

Registry'nin `role` alanı tutarsızdır - aynı ailenin üyeleri bir yerde
"Promoter", başka yerde "Regulatory" işaretlidir. Bu yüzden kategori,
rol yerine parça kimliğinin kalıbından belirlenir.
"""

from __future__ import annotations

import functools
import html
import json
import re
from pathlib import Path
from typing import Any

DATA = Path(__file__).resolve().parent.parent / "data" / "igem_parts.json"

# BioBrick iskeletleri parça değil tam plazmittir; onlar referans vektör
# olarak `data/reference_vectors/` altına yazılır ve buradan atlanır.
BACKBONES = {"pSB1C3", "pSB1A3", "pSB1K3", "pSB4C5"}

# ---------------------------------------------------------------------------
# Kategori eşlemesi
# ---------------------------------------------------------------------------

# Kimlik kalıbı -> kategori. Sıra önemlidir, ilk eşleşen kazanır.
_PATTERN_CATEGORY: list[tuple[str, str]] = [
    (r"^BBa_J231\d\d$", "promoter"),      # Anderson sabit promotör ailesi
    (r"^BBa_R00\d\d$", "promoter"),       # düzenlenebilir promotörler
    (r"^BBa_I719005$", "promoter"),       # T7 promotörü
    (r"^BBa_B00(29|3\d)$", "rbs"),
    (r"^BBa_B(001\d|100\d)$", "terminator"),
    (r"^BBa_M005\d$", "tag"),             # SsrA bozunma etiketleri
    (r"^BBa_C00(12|40|51|62|79)$", "regulator"),
    (r"^BBa_C0(061|161)$", "cds"),        # LuxI sentetaz
    (r"^BBa_E(0040|1010|0020|0030)$", "reporter"),
    (r"^BBa_(I13504|E0240|J04450|F2620)$", "device"),
]

_ROLE_FALLBACK = {
    "Promoter": "promoter",
    "Regulatory": "promoter",
    "Ribosome Entry Site": "rbs",
    "Terminator": "terminator",
    "Tag": "tag",
    "Reporter": "reporter",
    "CDS": "cds",
    "Engineered Region": "device",
    "Signalling": "device",
}


def _category(rec: dict[str, Any]) -> str:
    for pattern, cat in _PATTERN_CATEGORY:
        if re.match(pattern, rec["igem_id"]):
            return cat
    return _ROLE_FALLBACK.get(rec.get("role", ""), "cds")


# ---------------------------------------------------------------------------
# sigma70 konsensüsüne uzaklık
# ---------------------------------------------------------------------------
# Anderson ailesinin 20 üyesi yalnızca -35 ve -10 kutularında ayrışır;
# geri kalan 23 nükleotid aynıdır. Kutuların yeri bu dosyadaki dizilerden
# ölçülerek bulundu: toplam konsensüs uzaklığını en küçük yapan hizalama
# -35 icin 1..6, -10 icin 24..29 ve bu hizalamada BBa_J23119 konsensüsün
# tam kendisi çıkıyor - ailenin bilinen tasarımı da budur.
#
# ÖNEMLİ: buradaki sayı HESAPLANMIŞ bir benzerliktir, ölçülmüş promotör
# gücü değildir. İkisi aynı yönde gider ama birebir değildir; ölçüm
# değerleri Registry'nin kendi parça sayfalarındadır.

SIGMA70_35 = "TTGACA"
SIGMA70_10 = "TATAAT"
BOX_35 = (1, 6)    # 1 tabanlı, kapsayıcı
BOX_10 = (24, 29)


def consensus_distance(sequence: str) -> dict[str, Any] | None:
    """35 bp Anderson promotörü için -35/-10 kutularını ve uzaklığı döndür."""
    s = (sequence or "").upper()
    if len(s) != 35:
        return None
    b35 = s[BOX_35[0] - 1:BOX_35[1]]
    b10 = s[BOX_10[0] - 1:BOX_10[1]]
    d35 = sum(1 for a, b in zip(b35, SIGMA70_35) if a != b)
    d10 = sum(1 for a, b in zip(b10, SIGMA70_10) if a != b)
    return {"box35": b35, "box10": b10, "d35": d35, "d10": d10,
            "total": d35 + d10}


# ---------------------------------------------------------------------------
# Türkçe açıklamalar
# ---------------------------------------------------------------------------
# Registry'deki açıklamaların bir kısmı boş ya da "Replace later" gibi yer
# tutucudur (Anderson ailesinin tamamı böyle). Aşağıdaki notlar bu yüzden
# aileye göre yazıldı; parçaya özel olan `_PART_NOTE` içinde.

_FAMILY_NOTE: list[tuple[str, str]] = [
    (r"^BBa_J231", (
        "Anderson sabit promotör ailesinden. İndükleyici gerektirmez, "
        "konağın kendi sigma70 RNA polimerazıyla sürekli çalışır. Ailenin "
        "20 üyesi yalnızca -35 ve -10 kutularında ayrışır; konsensüse "
        "yaklaştıkça güç artar. Sentetik devrelerde ifade düzeyini "
        "promotörü değiştirerek ayarlamanın en pratik yolu budur."
    )),
    (r"^BBa_B00(29|3\d)$", (
        "Weiss serisinden ribozom bağlanma bölgesi. Aynı promotör ve aynı "
        "gen ile yalnızca RBS'i değiştirerek çeviri düzeyini kademeli "
        "ayarlamak için tasarlanmıştır; güçleri birbirine göre tanımlıdır. "
        "Başlangıç kodonuna uzaklık (spacing) burada da belirleyicidir - "
        "Dizi analizi sayfası bu aralığı ölçer."
    )),
    (r"^BBa_B(001\d|100\d)$", (
        "Transkripsiyon terminatörü. Kasetin sonuna konur; yoksa RNA "
        "polimeraz aşağı akıştaki dizileri de okur ve komşu genlerin "
        "ifadesi bozulur. Çift terminatörler (B0014, B0015) tek olanlara "
        "göre daha güvenli kapatır."
    )),
    (r"^BBa_R00\d\d$", (
        "Düzenlenebilir promotör: çalışması için kendi baskılayıcısının ya "
        "da aktivatörünün hücrede bulunması gerekir. Düzenleyiciyi aynı "
        "plazmite ya da konağa eklemeyi atlarsanız promotör ya hep açık ya "
        "hep kapalı kalır."
    )),
    (r"^BBa_C00(12|40|51|62|79)$", (
        "Düzenleyici protein kodlayan dizi (ribozom bağlanma bölgesi "
        "İÇERMEZ - önüne bir RBS eklemeniz gerekir). Kendi hedef "
        "promotörüyle birlikte kullanılır."
    )),
    (r"^BBa_M005\d$", (
        "SsrA bozunma etiketi: proteinin C ucuna eklenir ve ClpXP "
        "proteazına hedef gösterir. Devrenin düşüş süresini kısaltmak, "
        "yani sinyal kesildiğinde proteinin hızla yok olmasını sağlamak "
        "için kullanılır. Kararlı protein üretmek istiyorsanız bunu "
        "EKLEMEYİN."
    )),
    (r"^BBa_E(0040|1010|0020|0030)$", (
        "Floresan reporter kodlayan dizi. Ribozom bağlanma bölgesi "
        "içermez; önüne bir RBS, sonuna terminatör gerekir."
    )),
    (r"^BBa_(I13504|E0240|J04450|F2620)$", (
        "Bileşik parça: içinde promotör, RBS, kodlayan dizi ve terminatör "
        "gibi birden fazla eleman birlikte gelir. Tek başına çalışır, bu "
        "yüzden ayrıca promotör ya da RBS eklemeniz gerekmez - eklerseniz "
        "iki kaset üst üste biner."
    )),
]

_PART_NOTE: dict[str, str] = {
    "BBa_J23119": (
        "Ailenin konsensüs üyesi: -35 ve -10 kutuları sigma70 "
        "konsensüsünün tam kendisi. Aileyi karşılaştırırken üst sınır "
        "olarak bunu alın."
    ),
    "BBa_B0015": (
        "Çift terminatör (B0010 + B0012). iGEM yapılarında kasetin sonunda "
        "fiilen standart haline gelmiştir."
    ),
    "BBa_R0010": (
        "pLac: LacI tarafından baskılanır, IPTG ile açılır. Ayrıca CAP "
        "bağlanma bölgesi taşıdığı için ortamda glukoz varken katabolit "
        "baskılanmasıyla zayıflar."
    ),
    "BBa_R0040": (
        "pTet: TetR tarafından baskılanır, anhidrotetrasiklin (aTc) ile "
        "açılır. aTc antibiyotik etkisi neredeyse olmadan baskıyı "
        "kaldırdığı için devrelerde tetrasikline tercih edilir."
    ),
    "BBa_I719005": (
        "T7 promotörü. Konakta T7 RNA polimeraz yoksa (DE3 lizojeni "
        "değilse) hiç transkripsiyon olmaz."
    ),
}


def _note_for(rec: dict[str, Any]) -> str:
    note = ""
    for pattern, text in _FAMILY_NOTE:
        if re.match(pattern, rec["igem_id"]):
            note = text
            break
    extra = _PART_NOTE.get(rec["igem_id"])
    if extra:
        note = f"{note} {extra}".strip()
    own = _clean_text(rec.get("description", ""))
    title = _clean_text(rec.get("title", ""))
    # Registry'nin kendi açıklaması yer tutucuysa taşımanın anlamı yok.
    if own and own.lower().strip(" .") in ("replace later", "later", "none"):
        own = ""
    kaynak = f"Registry kaydı: “{title}”" if title else ""
    if own:
        kaynak = f"{kaynak} - {own}" if kaynak else own
    return f"{note}\n\n{kaynak}".strip() if kaynak else note


_TAG_RE = re.compile(r"<[^>]+>")


def _clean_text(raw: str) -> str:
    """Registry açıklamalarındaki HTML kalıntılarını temizle."""
    s = raw or ""
    s = s.replace("\\n", " ").replace("\\r", " ")
    s = _TAG_RE.sub(" ", s)
    s = html.unescape(s)
    s = re.sub(r"\s+", " ", s).strip()
    return s[:600]


# ---------------------------------------------------------------------------
# Yükleme
# ---------------------------------------------------------------------------

# Hangi düzenleyicinin hangi promotöre lazım olduğu - tasarım denetiminin
# `requires` alanını besler.
_REQUIRES: dict[str, list[str]] = {
    "BBa_R0010": ["LacI baskılayıcısı (konak lacI ya da BBa_C0012)"],
    "BBa_R0011": ["LacI baskılayıcısı (konak lacI ya da BBa_C0012)"],
    "BBa_R0040": ["TetR baskılayıcısı (BBa_C0040)"],
    "BBa_R0051": ["Lambda cI baskılayıcısı (BBa_C0051)"],
    "BBa_R0062": ["LuxR aktivatörü (BBa_C0062) + 3OC6HSL"],
    "BBa_R0063": ["LuxR aktivatörü (BBa_C0062) + 3OC6HSL"],
    "BBa_R0071": ["RhlR aktivatörü + C4-HSL"],
    "BBa_R0079": ["LasR aktivatörü (BBa_C0079) + 3OC12HSL"],
    "BBa_R0080": ["AraC düzenleyicisi + arabinoz"],
    "BBa_R0082": ["OmpR / EnvZ iki bileşenli sistem"],
    "BBa_I719005": ["T7 RNA polimeraz (konak genomunda, örn. DE3)"],
}

_INDUCTION: dict[str, str] = {
    "BBa_R0010": "IPTG",
    "BBa_R0011": "IPTG",
    "BBa_R0040": "anhidrotetrasiklin (aTc)",
    "BBa_R0062": "3OC6HSL (AHL)",
    "BBa_R0063": "3OC6HSL (AHL)",
    "BBa_R0071": "C4-HSL",
    "BBa_R0080": "arabinoz",
    "BBa_I719005": "IPTG (konaktaki T7 kaseti üzerinden)",
}


@functools.lru_cache(maxsize=1)
def _raw() -> dict[str, Any]:
    if not DATA.exists():
        return {"parts": []}
    return json.loads(DATA.read_text(encoding="utf-8"))


def part_id_for(igem_id: str) -> str:
    """BBa_J23100 -> IGEM_J23100, pSB1C3 -> IGEM_PSB1C3."""
    core = igem_id[4:] if igem_id.startswith("BBa_") else igem_id
    return "IGEM_" + core.upper()


@functools.lru_cache(maxsize=1)
def igem_parts() -> tuple[dict[str, Any], ...]:
    """Katalog kaydı biçiminde iGEM parçaları (iskeletler hariç)."""
    out: list[dict[str, Any]] = []
    for rec in _raw().get("parts", []):
        if rec["igem_id"] in BACKBONES:
            continue
        if not str(rec.get("license", "")).startswith("cc-by"):
            continue
        seq = (rec.get("sequence") or "").upper()
        if not seq:
            continue
        out.append({
            "part_id": part_id_for(rec["igem_id"]),
            "name": f"{rec['igem_id']} - {_clean_text(rec.get('title', ''))[:60]}",
            "category": _category(rec),
            "sequence": seq,
            "seq_status": "full",
            "length": len(seq),
            "strength": None,
            "induction": _INDUCTION.get(rec["igem_id"]),
            "requires": _REQUIRES.get(rec["igem_id"], []),
            "incompatible_with": [],
            "copy_number": None,
            "inc_group": None,
            "notes": _note_for(rec),
            "accession": (f"iGEM Registry {rec['igem_id']} "
                          f"({rec.get('license', 'cc-by-4.0').upper()})"),
            "pmid": "",
            "source": "igem",
            "igem_id": rec["igem_id"],
            "igem_url": rec.get("url", ""),
        })
    return tuple(out)


def anderson_family() -> list[dict[str, Any]]:
    """Anderson promotör ailesi, konsensüs uzaklığına göre sıralı."""
    rows = []
    for p in igem_parts():
        if not p["igem_id"].startswith("BBa_J231"):
            continue
        d = consensus_distance(p["sequence"])
        if not d:
            continue
        rows.append({"igem_id": p["igem_id"], "part_id": p["part_id"],
                     "sequence": p["sequence"], "url": p["igem_url"], **d})
    return sorted(rows, key=lambda r: (r["total"], r["igem_id"]))


def license_note() -> str:
    return _raw().get("license_note", "")


def source_name() -> str:
    return _raw().get("source", "iGEM Registry")
