"""Var olan bir plazmiti denetim motoruna sokar.

Uygulamanın 42 denetim kuralı şimdiye kadar yalnızca *parçalardan
kurduğunuz* yapılara uygulanabiliyordu. Oysa elinizde hocanın verdiği
ya da bir kayıttan indirdiğiniz bir plazmit varsa, asıl sormak
istediğiniz soru odur: "bu yapıda bir sorun var mı?"

Buradaki iş, dizi üzerinde bulunan şeyleri denetim motorunun
anlayacağı biçime çevirmek. İki kaynak birleştiriliyor:

  1. Katalog taraması - parçalarımız dizide birebir ya da yakın
     eşleşmeyle bulunuyor mu (`annotate` + `annotate_similar`).
  2. Dosyanın kendi öznitelik tablosu - GenBank kaydı "burada bir
     rep_origin var" diyorsa, katalogda karşılığı olmasa bile bunu
     biliyoruz sayılır.

İkincisi olmazsa denetim yanıltır: katalogda olmayan bir orijin
taşıyan plazmit için "replikasyon orijini yok" hatası verilir ki bu
yanlıştır. Bulunamayan bölgeler ayrıca raporlanır - kullanıcı neyin
tanınmadığını görmeli.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from io import BytesIO, StringIO
from typing import Any

from core.assembly import AssembledVector
from core.seq_utils import (
    Feature, annotate, annotate_similar, clean_seq,
    collapse_overlapping,
)

# GenBank öznitelik tipi -> katalog kategorisi. Yalnızca kategorisi
# kesin olanlar; belirsizler denetime sokulmaz.
GB_TYPE_CATEGORY = {
    "rep_origin": "ori",
    "terminator": "terminator",
    "promoter": "promoter",
    "protein_bind": "operator",
    "RBS": "rbs",
}

# Desenler tam sözcük arar. Önceki sürüm `str` önekiyle eşleştiği için
# "L-strand Y effector site" gibi bir özniteliği direnç markerı sanmıştı;
# kısa önekler bu yüzden tam adlara çevrildi.
_MARKER_RE = re.compile(
    r"\b(bla|cat|catA|neo|nptII|aadA|aadD|tetA|tetR|ermC|erm|acc3|aph|"
    r"bleO|sul1|strA|ampicillin|chloramphenicol|kanamycin|tetracycline|"
    r"spectinomycin|gentamicin|erythromycin|resistance)\b", re.I)
_REPORTER_RE = re.compile(
    r"(gfp|cfp|yfp|rfp|mcherry|dsred|lacZ|lux[ABCDE]|luc|phoA)", re.I)
_REGULATOR_RE = re.compile(
    r"\b(lacI|araC|xylS|tetR|cI\b|luxR|lasR|rhaR|trfA|rep[ABC]?|rop|rom)",
    re.I)
_PROMOTER_RE = re.compile(r"(promoter|^P[A-Z_]|lacUV5|araBAD|T7|trc|tac)", re.I)
_TERMINATOR_RE = re.compile(r"(terminator|^T[01]\b|rrnB|tL3)", re.I)
_ORI_RE = re.compile(r"(\bori\b|origin|rep_origin|pMB1|ColE1|p15A|pSC101|R6K)",
                     re.I)


@dataclass
class Inspection:
    name: str
    sequence: str
    circular: bool
    length: int
    source: str                       # "fasta" | "genbank" | "metin"
    catalog_features: list[Feature] = field(default_factory=list)
    gb_features: list[dict[str, Any]] = field(default_factory=list)
    derived_parts: list[dict[str, Any]] = field(default_factory=list)
    unidentified: list[dict[str, Any]] = field(default_factory=list)
    coverage: float = 0.0             # dizinin tanınan oranı (%)
    notes: list[dict[str, str]] = field(default_factory=list)


def _note(kind: str, title: str, detail: str) -> dict[str, str]:
    return {"kind": kind, "title": title, "detail": detail}


# ---------------------------------------------------------------------------
# Dosya okuma
# ---------------------------------------------------------------------------


def parse_upload(data: bytes, filename: str = "") -> dict[str, Any]:
    """FASTA, GenBank ya da düz diziden plazmit oku."""
    metin = data.decode("utf-8", errors="ignore")
    ad = (filename or "plazmit").rsplit(".", 1)[0]

    if metin.lstrip().startswith("LOCUS"):
        from Bio import SeqIO
        try:
            rec = SeqIO.read(StringIO(metin), "genbank")
        except Exception:
            return {"error": "GenBank dosyası okunamadı."}
        feats = []
        for f in rec.features:
            if f.type in ("source", "gene"):
                continue
            q = f.qualifiers
            etiket = ""
            for anahtar in ("gene", "product", "standard_name", "label",
                            "note"):
                if q.get(anahtar):
                    etiket = str(q[anahtar][0])
                    break
            feats.append({
                "type": f.type,
                "start": int(f.location.start) + 1,
                "end": int(f.location.end),
                "strand": int(f.location.strand or 1),
                "label": etiket[:80],
            })
        return {
            "name": rec.name or ad,
            "sequence": str(rec.seq).upper(),
            "circular": rec.annotations.get("topology", "") == "circular",
            "features": feats,
            "source": "genbank",
        }

    if metin.lstrip().startswith(">"):
        ilk = metin.lstrip().splitlines()[0]
        return {"name": ilk[1:].strip()[:40] or ad,
                "sequence": clean_seq(metin), "circular": True,
                "features": [], "source": "fasta"}

    return {"name": ad, "sequence": clean_seq(metin), "circular": True,
            "features": [], "source": "metin"}


# ---------------------------------------------------------------------------
# Öznitelikten kategoriye
# ---------------------------------------------------------------------------


def _category_from_gb(feat: dict[str, Any]) -> str | None:
    tip = feat.get("type", "")
    etiket = (feat.get("label") or "").strip()
    if tip in GB_TYPE_CATEGORY:
        return GB_TYPE_CATEGORY[tip]
    if not etiket:
        return None
    if tip == "regulatory":
        if _TERMINATOR_RE.search(etiket):
            return "terminator"
        if _PROMOTER_RE.search(etiket):
            return "promoter"
        return None
    if tip in ("CDS", "misc_feature"):
        # Orijin YALNIZCA rep_origin tipinden turetilir. misc_feature
        # icinde "pSC101" gecmesi orada orijin oldugu anlamina gelmez;
        # pBR322 kaydindaki o not bir benzerlik aciklamasidir.
        if _MARKER_RE.search(etiket):
            return "marker"
        if _REPORTER_RE.search(etiket):
            return "reporter"
        if _REGULATOR_RE.search(etiket):
            return "regulator"
    return None


# ---------------------------------------------------------------------------
# İnceleme
# ---------------------------------------------------------------------------


def inspect(
    data: dict[str, Any], parts: list[dict[str, Any]],
    min_identity: float = 92.0,
) -> Inspection:
    """Plazmiti tara: ne tanındı, ne tanınmadı, denetime ne girecek."""
    seq = clean_seq(data.get("sequence", ""))
    ins = Inspection(
        name=data.get("name") or "plazmit", sequence=seq,
        circular=bool(data.get("circular", True)), length=len(seq),
        source=data.get("source", "metin"),
        gb_features=list(data.get("features") or []),
    )
    if not seq:
        ins.notes.append(_note("hata", "Dizi okunamadı",
                               "Dosyadan dizi çıkarılamadı."))
        return ins

    tam = annotate(seq, parts, circular=ins.circular)
    yakin = annotate_similar(seq, parts, circular=ins.circular,
                             min_identity=min_identity)
    gorulen = {(f.part_id, f.start) for f in tam}
    birlesik = list(tam)
    for f in yakin:
        if (f.part_id, f.start) not in gorulen:
            birlesik.append(f)
    # Aynı bölgeye düşen kayıtları daralt: katalogda `bla` birden çok
    # kez geçiyor ve üçü de aynı yere oturuyor. Daraltmazsak denetim
    # "birden fazla orijin / marker" diye yanlış hata verir.
    ins.catalog_features = collapse_overlapping(birlesik, len(seq))

    # --- Denetime girecek parça listesi
    parts_by_id = {p["part_id"]: p for p in parts}
    eklenen: set[str] = set()
    for f in ins.catalog_features:
        rec = parts_by_id.get(f.part_id)
        if not rec or f.part_id in eklenen:
            continue
        kayit = dict(rec)
        kayit["_identity"] = f.identity
        kayit["_span"] = (f.start, f.end)
        ins.derived_parts.append(kayit)
        eklenen.add(f.part_id)

    # --- Dosyanın kendi öznitelikleri: katalogda karşılığı olmayanlar
    kapsanan = _covered_positions(ins.catalog_features, len(seq))
    for feat in ins.gb_features:
        kat = _category_from_gb(feat)
        if not kat:
            continue
        bas, son = feat["start"], feat["end"]
        if son < bas or son > len(seq):
            continue
        orta = (bas + son) // 2
        if orta in kapsanan:
            continue          # katalog zaten buldu
        sahte_id = f"GB_{kat.upper()}_{bas}"
        if sahte_id in eklenen:
            continue
        ins.derived_parts.append({
            "part_id": sahte_id,
            "name": f"{feat['label'] or kat} (dosyadan)",
            "category": kat,
            "sequence": seq[bas - 1:son],
            "seq_status": "full",
            "length": son - bas + 1,
            "strength": None, "induction": None,
            "requires": [], "incompatible_with": [],
            "copy_number": None, "inc_group": None,
            "notes": "Kaydın kendi öznitelik tablosundan alındı; "
                     "katalogda karşılığı bulunamadı.",
            "accession": "", "pmid": "", "source": "dosya",
            "_identity": None, "_span": (bas, son),
        })
        eklenen.add(sahte_id)

    # --- Tanınmayan bölgeler
    ins.unidentified = _gaps(ins.derived_parts, len(seq))
    taninan = len(seq) - sum(g["length"] for g in ins.unidentified)
    ins.coverage = round(100 * taninan / max(len(seq), 1), 1)
    ins.notes = _interpret(ins)
    return ins


def _covered_positions(feats: list[Feature], length: int) -> set[int]:
    kapsanan: set[int] = set()
    for f in feats:
        if f.start <= f.end:
            kapsanan.update(range(f.start, f.end + 1))
        else:
            kapsanan.update(range(f.start, length + 1))
            kapsanan.update(range(1, f.end + 1))
    return kapsanan


def _gaps(parcalar: list[dict[str, Any]], length: int,
          min_gap: int = 120) -> list[dict[str, Any]]:
    kapsanan = set()
    for p in parcalar:
        bas, son = p.get("_span", (0, 0))
        if not bas:
            continue
        if bas <= son:
            kapsanan.update(range(bas, son + 1))
        else:
            kapsanan.update(range(bas, length + 1))
            kapsanan.update(range(1, son + 1))
    bosluklar = []
    bas = None
    for x in range(1, length + 1):
        if x not in kapsanan and bas is None:
            bas = x
        elif x in kapsanan and bas is not None:
            if x - bas >= min_gap:
                bosluklar.append({"start": bas, "end": x - 1,
                                  "length": x - bas})
            bas = None
    if bas is not None and length - bas + 1 >= min_gap:
        bosluklar.append({"start": bas, "end": length,
                          "length": length - bas + 1})
    return bosluklar


def to_vector(ins: Inspection) -> AssembledVector:
    """Denetim motorunun anlayacağı yapıyı kur.

    Önemli ayrım: `sequence` yüklenen dizinin KENDİSİDİR, parçaların
    birleşimi değil. Yani uzunluk, GC ve enzim haritası gerçek
    plazmitin değerleridir; `part_order` ise yalnızca denetim
    kurallarının "bu kategoriden eleman var mı" sorusuna cevap verir.
    """
    return AssembledVector(
        name=ins.name, sequence=ins.sequence, features=ins.catalog_features,
        circular=ins.circular, part_order=list(ins.derived_parts),
    )


def _interpret(ins: Inspection) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    katalog = len(ins.catalog_features)
    dosyadan = sum(1 for p in ins.derived_parts
                   if p.get("source") == "dosya")

    out.append(_note(
        "bilgi" if ins.coverage > 50 else "uyari", "Tanıma oranı",
        f"Dizinin %{ins.coverage}'i tanımlı bir bölgeye denk geldi. "
        f"Katalogdan {katalog} eşleşme"
        + (f", dosyanın kendi tablosundan {dosyadan} öznitelik"
           if dosyadan else "")
        + f" kullanıldı. Denetim yalnızca TANINAN elemanlara bakar - "
          f"tanınmayan bir bölgede sorun varsa göremez."))

    if ins.source in ("fasta", "metin") and katalog < 3:
        out.append(_note(
            "uyari", "Öznitelik tablosu yok",
            "Düz dizi yüklediniz ve katalog taraması az eşleşme buldu. "
            "GenBank biçiminde yüklerseniz dosyanın kendi öznitelikleri "
            "de kullanılır ve denetim çok daha doğru olur."))

    if ins.unidentified:
        en_buyuk = max(ins.unidentified, key=lambda g: g["length"])
        toplam = sum(g["length"] for g in ins.unidentified)
        out.append(_note(
            "bilgi", "Tanınmayan bölgeler",
            f"{len(ins.unidentified)} bölge tanınmadı, toplam {toplam} bp. "
            f"En büyüğü {en_buyuk['start']}–{en_buyuk['end']} "
            f"({en_buyuk['length']} bp). Bunlar büyük ihtimalle insertiniz "
            f"ya da katalogda karşılığı olmayan elemanlardır; Dizi analizi "
            f"sayfasında açık okuma çerçevesi arayarak bakabilirsiniz."))

    varyant = [f for f in ins.catalog_features if f.identity < 100]
    if varyant:
        out.append(_note(
            "bilgi", "Varyant eşleşmeler",
            f"{len(varyant)} bölge katalogdakine yakın ama birebir değil "
            f"(en düşüğü %{min(f.identity for f in varyant)}). Bu "
            f"beklenen bir durumdur: aynı gen farklı plazmit soylarında "
            f"birkaç nükleotid farkla taşınır."))
    return out
