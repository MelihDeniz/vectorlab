"""Referans bakteriyel vektörler: gerçek, tam dizili plazmitler.

Diziler ve GenBank öznitelikleri `data/reference_vectors/*.json` içinde
durur; o dosyalar `tools/fetch_reference_vectors.py` ile NCBI'dan çekilir.
Uygulama çalışırken ağa çıkmaz.

Bu modüldeki `CURATED` sözlüğü editöryel katmandır: hangi plazmit ne işe
yarar, kopya sayısı, uyumsuzluk grubu, dikkat edilecekler. Dizi verisi
NCBI'dan, yorum buradan gelir - ikisi bilinçli olarak ayrı tutulur.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Any

from core.seq_utils import (
    Feature,
    annotate,
    annotate_similar,
    gc_percent,
    revcomp,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "reference_vectors"

# ---------------------------------------------------------------------------
# Editöryel katman
# ---------------------------------------------------------------------------
# `lit_query`: Literatür sayfasının Europe PMC'de arayacağı sorgu. Bilerek
# PMID yazmıyoruz - ezberden atıf numarası uydurmak sessiz hata üretir.

CURATED: dict[str, dict[str, Any]] = {
    "PUC19": {
        "name": "pUC19",
        "ozet": "Yüksek kopyalı standart klonlama vektörü; mavi/beyaz tarama.",
        "kullanim": "klonlama",
        "kopya": "yüksek (~500-700)",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac (lacUV5)",
        "mcs": "HindIII-SphI-PstI-SalI-XbaI-BamHI-SmaI-KpnI-SacI-EcoRI",
        "notlar": (
            "pBR322 iskeleti + M13mp19'dan gelen lac bölgesi. MCS, lacZ-alfa "
            "parçasının içine yerleştirilmiştir: insert girdiğinde alfa "
            "komplementasyon bozulur ve koloni X-gal üzerinde beyaz kalır. "
            "Bu taramanın çalışması için konağın lacZ-delta-M15 taşıması "
            "gerekir (DH5-alfa, TOP10). Kopya sayısının yüksekliği "
            "orijindeki nokta mutasyonundan gelir; rop/rom geni yoktur."
        ),
        "lit_query": "pUC19 cloning vector sequence Yanisch-Perron",
    },
    "PUC18": {
        "name": "pUC18",
        "ozet": "pUC19'un MCS'i ters yönde olan ikizi.",
        "kullanim": "klonlama",
        "kopya": "yüksek (~500-700)",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac (lacUV5)",
        "mcs": "EcoRI-SacI-KpnI-SmaI-BamHI-XbaI-SalI-PstI-SphI-HindIII",
        "notlar": (
            "Dizisi pUC19 ile neredeyse aynıdır, fark MCS'in yönüdür. "
            "Insertin lac promotörüne göre hangi yönde gireceğini seçmek "
            "istediğinizde ikisi arasından seçim yapılır."
        ),
        "lit_query": "pUC18 pUC19 vector Yanisch-Perron Vieira Messing",
    },
    "PBR322": {
        "name": "pBR322",
        "ozet": "Klonlamanın atası; iki ayrı direnç markerı taşır.",
        "kullanim": "klonlama",
        "kopya": "düşük (~15-20)",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)", "TetR (tetA)"],
        "promotor": "yok (kendi genlerinin promotörleri var)",
        "mcs": "yok - tek kesen bölgeler gen içlerine dağılmış",
        "notlar": (
            "Iki markerı olması insersiyonel inaktivasyona izin verir: "
            "BamHI ya da SalI ile tetA içine klonlarsanız koloni ampisiline "
            "dirençli ama tetrasikline duyarlı kalır. rop (rom) geni kopya "
            "sayısını baskılar - pUC19'un yüksek kopyalı olmasının bir "
            "nedeni de bu genin oradan çıkarılmış olmasıdır."
        ),
        "lit_query": "pBR322 construction characterization Bolivar Rodriguez",
    },
    "PACYC184": {
        "name": "pACYC184",
        "ozet": "p15A orijinli; ColE1 vektörlerle aynı hücrede durur.",
        "kullanim": "klonlama",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["CmR (cat)", "TetR (tetA)"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "Iki plazmitli kurulumların klasik ikinci üyesi: orijini ColE1 "
            "ile aynı uyumsuzluk grubunda olmadığı için pUC/pET ile birlikte "
            "tutulabilir. Rosetta ve pLysS suşlarının taşıdığı yardımcı "
            "plazmitler de p15A orijinlidir - üç plazmitli kurulumda bu "
            "çakışmaya dikkat edin."
        ),
        "lit_query": "pACYC184 p15A replicon Chang Cohen plasmid",
    },
    "PACYC177": {
        "name": "pACYC177",
        "ozet": "p15A orijinli, ampisilin + kanamisin dirençli.",
        "kullanim": "klonlama",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["AmpR (bla)", "KanR"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "pACYC184'ün kardeşi; marker çifti farklıdır. GenBank kaydı "
            "doğrusal olarak yatırılmıştır ama molekül daireseldir - "
            "uygulama bu kaydı dairesel olarak işler."
        ),
        "lit_query": "pACYC177 p15A cloning vector Chang Cohen",
    },
    "PBAD24": {
        "name": "pBAD24",
        "ozet": "Arabinozla indüklenen, sıkı kontrollü ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "araBAD (P-BAD), araC düzenleyicisiyle",
        "mcs": "NheI-KpnI-SmaI-XbaI-SalI-PstI-SphI-HindIII",
        "notlar": (
            "T7 sistemine göre asıl avantajı doz yanıtlı olmasıdır: "
            "arabinoz derişimini değiştirerek ifade düzeyini kademeli "
            "ayarlayabilirsiniz, IPTG'de bu kadar düzgün bir kademe yoktur. "
            "Konak arabinozu yıkmamalı (araBAD mutantı, örn. TOP10) ve "
            "glukoz katabolit baskılanması yaptığı için besiyerinde "
            "bulunmamalıdır."
        ),
        "lit_query": "pBAD24 tight regulation arabinose promoter Guzman",
    },
    "PBAD18": {
        "name": "pBAD18",
        "ozet": "pBAD24'ün MCS'li, çok yönlü kardeşi.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "araBAD (P-BAD)",
        "mcs": "pUC tipi genişletilmiş MCS",
        "notlar": (
            "pBAD24 ile aynı düzenleme mantığı; farkı klonlama bölgesidir. "
            "Başlangıç kodonunun promotöre uzaklığı pBAD24'te daha kısadır, "
            "pBAD18'de araya MCS girer."
        ),
        "lit_query": "pBAD18 arabinose expression vector Guzman Beckwith",
    },
    "PBAD33": {
        "name": "pBAD33",
        "ozet": "Arabinozla indüklenen, p15A orijinli ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["CmR (cat)"],
        "promotor": "araBAD (P-BAD)",
        "mcs": "pBAD24 ile aynı dizilim",
        "notlar": (
            "Iki proteini ayrı plazmitlerden ifade edeceğiniz kurulumların "
            "arabinoz kolu: ColE1 orijinli bir pET ya da pUC ile birlikte "
            "durabilir ve markerı kloramfenikoldür, yani ampisilin ya da "
            "kanamisin ile çakışmaz."
        ),
        "lit_query": "pBAD33 arabinose vector p15A chloramphenicol Guzman",
    },
    "PTRC99A": {
        "name": "pTrc99A",
        "ozet": "IPTG ile indüklenen trc promotörlü ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "trc (trp -35 + lacUV5 -10)",
        "mcs": "NcoI-EcoRI-SmaI-BamHI-XbaI-SalI-PstI-HindIII",
        "notlar": (
            "T7 gerektirmez: trc promotörü konağın kendi RNA polimerazıyla "
            "çalışır, yani DE3 lizojeni olmayan suşlarda da ifade alırsınız. "
            "Plazmit kendi lacI genini taşır, bu yüzden yüksek kopyada bile "
            "baskılama ayakta kalır. Aşağı akışta rrnB T1T2 terminatörü var."
        ),
        "lit_query": "pTrc99A expression vector Amann Brosius trc promoter",
    },
    "PBSIIKS": {
        "name": "pBluescript II KS(+)",
        "ozet": "MCS'i T3 ve T7 promotörleri arasında; in vitro transkripsiyon.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac + MCS'i kuşatan T3 / T7",
        "mcs": "pUC tipi MCS, iki yanında T7 ve T3 promotörleri",
        "notlar": (
            "Inserti iki yönden de RNA'ya çevirebilirsiniz: bir yandan T7, "
            "diğer yandan T3 promotörü. f1 (M13) orijini sayesinde tek "
            "zincirli DNA üretilebilir - dizileme ve mutageneze yönelik "
            "tasarlanmıştır. lacZ-alfa taraması pUC19'da olduğu gibi çalışır."
        ),
        "lit_query": "pBluescript II phagemid vector Short Alting-Mees",
    },
    "PSC101": {
        "name": "pSC101",
        "ozet": "Doğal düşük kopyalı plazmit; üçüncü uyumsuzluk grubu.",
        "kullanim": "klonlama",
        "kopya": "çok düşük (~5)",
        "inc_group": "pSC101",
        "markerlar": ["TetR (tet)"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "Klonlama vektörü değil, doğal bir plazmit: moleküler biyolojide "
            "klonlanan ilk plazmit olduğu için ve orijini ColE1 ile p15A'nın "
            "ikisinden de ayrı olduğu için burada. Üç plazmitli kurulumlarda "
            "üçüncü orijin olarak bu replikon kullanılır. Kopya sayısı çok "
            "düşük olduğundan toksik ürünler için elverişlidir; buna karşılık "
            "plazmit izolasyonu az DNA verir."
        ),
        "lit_query": "pSC101 replicon copy number incompatibility Cohen Chang",
    },
    # --- iGEM BioBrick iskeletleri (dizi: iGEM Registry, CC BY 4.0) --------
    # Aşağıdaki içerik tahmin değil: her iskelet katalog parçalarıyla
    # taranarak doğrulandı. pSB1C3/1A3/1K3'te pUC orijini 212..866
    # konumunda, pSB4C5'te pSC101 orijini 689..2019 konumunda çıkıyor -
    # yani addaki "1" yüksek kopyayı, "4" pSC101'i gösteriyor.
    "PSB1C3": {
        "name": "pSB1C3",
        "ozet": "iGEM'in standart yüksek kopyalı montaj iskeleti (CmR).",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["CmR (cat)"],
        "promotor": "yok - parçanın kendi promotörü gelir",
        "mcs": "BioBrick ön eki (EcoRI-XbaI) ve arka eki (SpeI-PstI)",
        "notlar": (
            "iGEM dağıtımının varsayılan iskeleti: parçalar buraya ön ek ve "
            "arka ek arasına takılır. Insert bölgesi numaralandırmanın "
            "başlangıcını kuşatır - EcoRI ve XbaI dizinin sonunda, SpeI ve "
            "PstI başında görünür. RFC10 montajında ön ek EcoRI-XbaI, arka "
            "ek SpeI-PstI ile kesilir; XbaI ve SpeI uçları ligasyonda "
            "birleşerek üründe 6 bp'lik bir iz (scar) bırakır, bu yüzden "
            "protein füzyonları için uygun değildir."
        ),
        "lit_query": "BioBrick standard assembly RFC10 Knight registry",
    },
    "PSB1A3": {
        "name": "pSB1A3",
        "ozet": "pSB1C3'ün ampisilin dirençli eşi.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "pSB1C3 ile aynı iskelet, markerı ampisilindir. İki BioBrick "
            "yapısını aynı hücrede tutmak istediğinizde marker farkı işe "
            "yarar - ama orijinleri aynı olduğu için uyumsuzdur, uzun "
            "sürede biri kaybolur."
        ),
        "lit_query": "BioBrick standard assembly plasmid backbone registry",
    },
    "PSB1K3": {
        "name": "pSB1K3",
        "ozet": "pSB1C3'ün kanamisin dirençli eşi.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["KanR"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "Kanamisin seçilimi ampisiline göre daha temizdir: beta-laktamaz "
            "besiyerine salınıp antibiyotiği bozduğu için ampisilin "
            "plakalarında satellit koloniler çıkar, kanamisinde çıkmaz. "
            "Haritada direnç geni otomatik etiketlenmez: bu iskeletin "
            "kanamisin geni katalogdaki nptII varyantından yeterince "
            "farklı, benzerlik eşiğinin altında kalıyor."
        ),
        "lit_query": "BioBrick kanamycin backbone pSB1K3 registry",
    },
    "PSB4C5": {
        "name": "pSB4C5",
        "ozet": "Düşük kopyalı BioBrick iskeleti (pSC101 orijini, CmR).",
        "kullanim": "klonlama",
        "kopya": "düşük (pSC101)",
        "inc_group": "pSC101",
        "markerlar": ["CmR (cat)"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "Yüksek kopyalı iskelet taşıyıcıya yük bindirdiğinde ya da "
            "ürün toksikse buna geçilir. Orijini pSC101 olduğu için "
            "pSB1 ailesiyle (ColE1) aynı hücrede durabilir - iki ayrı "
            "devreyi iki plazmitte taşımanın yolu bu."
        ),
        "lit_query": "pSB4C5 low copy BioBrick pSC101 backbone",
    },
}

ORDER = ["PUC19", "PUC18", "PBSIIKS", "PBR322", "PACYC184", "PACYC177",
         "PTRC99A", "PBAD24", "PBAD18", "PBAD33", "PSC101",
         "PSB1C3", "PSB1A3", "PSB1K3", "PSB4C5"]


# ---------------------------------------------------------------------------
# Yükleme
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=32)
def load(vector_id: str) -> dict[str, Any] | None:
    """Tam kaydı oku (dizi dahil). Bulunamazsa None."""
    path = DATA / f"{vector_id}.json"
    if not path.exists():
        return None
    rec = json.loads(path.read_text(encoding="utf-8"))
    rec.update(CURATED.get(vector_id, {}))
    # Dosya adındaki kimlik her zaman kazanır.
    rec["vector_id"] = vector_id
    rec.setdefault("name", vector_id)
    rec["gc"] = gc_percent(rec["sequence"])
    return rec


def available() -> list[str]:
    """Diskte verisi olan vektör kimlikleri, editöryel sırayla."""
    have = {p.stem for p in DATA.glob("*.json")}
    return [v for v in ORDER if v in have] + sorted(have - set(ORDER))


def list_vectors() -> list[dict[str, Any]]:
    """Listeleme için özet kayıtlar (dizi taşımaz)."""
    out = []
    for vid in available():
        rec = load(vid)
        if not rec:
            continue
        out.append({k: v for k, v in rec.items()
                    if k not in ("sequence", "features")})
    return out


# ---------------------------------------------------------------------------
# Öznitelik katmanları
# ---------------------------------------------------------------------------


def _spans(f: Feature, length: int) -> list[tuple[int, int]]:
    if f.start <= f.end:
        return [(f.start, f.end)]
    return [(f.start, length), (1, f.end)]


def _overlaps(a: Feature, b: Feature, length: int) -> bool:
    """Iki öznitelik büyük ölçüde üst üste mi biniyor (dairesel)."""
    shared = 0
    for s1, e1 in _spans(a, length):
        for s2, e2 in _spans(b, length):
            shared += max(0, min(e1, e2) - max(s1, s2) + 1)
    return shared >= 0.6 * min(a.length, b.length)


def catalog_features(
    vector_id: str, parts: list[dict[str, Any]], min_identity: float = 90.0
) -> list[Feature]:
    """Katalog parçalarını vektörde bul; çakışanlardan en iyisini tut.

    `pUC ori` ile `pBR322 ori` arasında tek nükleotid fark vardır, bu
    yüzden ikisi de aynı bölgeye düşer. Aynı kategoriden çakışan
    eşleşmelerde önce tam eşleşme, sonra en yüksek benzerlik kazanır -
    yoksa harita aynı bölgeyi iki kez etiketler.
    """
    rec = load(vector_id)
    if not rec:
        return []
    L = rec["length"]
    exact = annotate(rec["sequence"], parts, circular=True)
    fuzzy = annotate_similar(rec["sequence"], parts, circular=True,
                             min_identity=min_identity)

    keep: list[Feature] = list(exact)
    for f in sorted(fuzzy, key=lambda x: (-x.identity, -x.length)):
        if any(k.part_id == f.part_id and k.start == f.start for k in keep):
            continue
        if any(k.category == f.category and _overlaps(k, f, L) for k in keep):
            continue
        keep.append(f)
    return sorted(keep, key=lambda f: (f.start, -f.length))


def genbank_features(vector_id: str) -> list[dict[str, Any]]:
    """Kaydın kendi öznitelik tablosu (NCBI'dan geldiği gibi)."""
    rec = load(vector_id)
    return list(rec.get("features") or []) if rec else []


# ---------------------------------------------------------------------------
# Bölge alma ve dışa aktarım
# ---------------------------------------------------------------------------


def region(vector_id: str, start: int, end: int, strand: int = 1) -> str:
    """1 tabanlı, kapsayıcı bölgeyi döndür; dairesel sarmayı destekler."""
    rec = load(vector_id)
    if not rec:
        return ""
    seq, L = rec["sequence"], rec["length"]
    start = max(1, min(L, int(start)))
    end = max(1, min(L, int(end)))
    sub = seq[start - 1:end] if start <= end else seq[start - 1:] + seq[:end]
    return revcomp(sub) if strand == -1 else sub


def to_fasta(vector_id: str) -> str:
    rec = load(vector_id)
    if not rec:
        return ""
    body = "\n".join(rec["sequence"][i:i + 70]
                     for i in range(0, rec["length"], 70))
    return (f">{rec['name']} | {rec['accession']} | {rec['length']} bp | "
            f"dairesel\n{body}\n")
