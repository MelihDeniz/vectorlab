"""Referans vektörlerin öznitelik tablolarından katalog parçası çıkarır.

Elde 126 tam dizili plazmit ve her birinin GenBank öznitelik tablosu
var. Bu tablolardaki orijinler, terminatörler, promotörler, direnç
genleri ve muhabirler zaten tanımlı bölgelerdir - onları katalog
parçasına çevirmek, elle parça yazmaktan hem daha hızlı hem daha
güvenilirdir, çünkü çıkan her parça **aksesyon + koordinat** taşır.

Ölçütler bilinçli olarak dar: yalnızca kategorisi kesin olan öznitelik
tipleri alınır, dizisi belirsiz (N içeren) ya da aşırı kısa/uzun olanlar
elenir ve katalogda zaten birebir bulunan diziler atlanır. Amaç sayıyı
şişirmek değil, kanıtlı parça eklemek.

    python tools/extract_parts_from_vectors.py
"""

from __future__ import annotations

import datetime as _dt
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core import db, ref_vectors as rv              # noqa: E402
from core.seq_utils import revcomp                  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "data" / "extracted_parts.json"

MIN_LEN, MAX_LEN = 25, 4000

# Öznitelik tipi -> kategori. Yalnızca kategorisi kesin olanlar.
TYPE_CATEGORY = {
    "rep_origin": "ori",
    "terminator": "terminator",
    "promoter": "promoter",
    "protein_bind": "operator",
    "RBS": "rbs",
}

# CDS'ler etikete göre ayrılır; tanınmayan CDS alınmaz.
MARKER_RE = re.compile(
    r"^(bla|neo|cat|catA|aadA|aadD|tetA|tet|acc3|aph|erm|ermC|kan|"
    r"bleO|sul|str)\b", re.I)
REPORTER_RE = re.compile(
    r"(gfp|cfp|yfp|rfp|mcherry|dsred|lacZ|lux[ABCDE]|luc|phoA)", re.I)
REGULATOR_RE = re.compile(
    r"^(lacI|araC|xylS|tetR|cI|luxR|lasR|rhaR|rhaS|trfA|repA|repB|repC|"
    r"rep|mobA|mobB|mobC|mobV|traD|korB|incC)\b", re.I)

# Etiketi bunlara benzeyen oznitelikler atlanir: kategorisi belirsiz ya
# da parca olarak anlamsiz.
# Bu kayitlarda terminatorler `regulatory` tipinde ve kisa adlarla
# geciyor (T0, T1, rrnB, tL3) - tip tek basina yetmiyor, etikete de
# bakmak gerekiyor.
TERMINATOR_RE = re.compile(
    r"^(T[01]\b|T7\s*term|rrnB|tL3|TE\b|term)", re.I)
OPERATOR_RE = re.compile(r"operator", re.I)

SKIP_RE = re.compile(
    r"(priming site|conflict|misc_difference|unknown|hypothetical|"
    r"inverted repeat|fragment|^bp \d|complement of)", re.I)


def category_for(feat: dict) -> str | None:
    tip = feat.get("type", "")
    etiket = (feat.get("label") or "").strip()
    if not etiket or SKIP_RE.search(etiket):
        return None
    if tip in TYPE_CATEGORY:
        return TYPE_CATEGORY[tip]
    if tip == "regulatory":
        if (re.search(r"terminator", etiket, re.I)
                or TERMINATOR_RE.match(etiket)):
            return "terminator"
        if re.search(r"promoter|^P[a-z]|promot", etiket, re.I):
            return "promoter"
        return None
    if tip == "misc_feature" and OPERATOR_RE.search(etiket):
        return "operator"
    if tip == "CDS":
        if MARKER_RE.match(etiket):
            return "marker"
        if REPORTER_RE.search(etiket):
            return "reporter"
        if REGULATOR_RE.match(etiket):
            return "regulator"
    return None


def slug(etiket: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", etiket).strip("_").upper()
    return s[:22] or "PARCA"


def main() -> int:
    db.init_db()
    mevcut = {p["sequence"]: p["part_id"] for p in db.list_parts()
              if p["seq_status"] == "full" and p["sequence"]}

    bulunan: dict[str, dict] = {}      # dizi -> kayit
    atlanan = {"kategori_yok": 0, "uzunluk": 0, "belirsiz": 0,
               "katalogda_var": 0}

    for vid in rv.available():
        rec = rv.load(vid)
        if not rec:
            continue
        for feat in rec.get("features", []):
            kat = category_for(feat)
            if not kat:
                atlanan["kategori_yok"] += 1
                continue
            bas, son, yon = feat["start"], feat["end"], feat["strand"]
            if son < bas or (son - bas + 1) < MIN_LEN or \
                    (son - bas + 1) > MAX_LEN:
                atlanan["uzunluk"] += 1
                continue
            dizi = rec["sequence"][bas - 1:son]
            if yon < 0:
                dizi = revcomp(dizi)
            if set(dizi) - set("ACGT"):
                atlanan["belirsiz"] += 1
                continue
            if dizi in mevcut:
                atlanan["katalogda_var"] += 1
                continue
            if dizi in bulunan:
                bulunan[dizi]["also_in"].append(rec["name"])
                continue

            etiket = feat["label"].strip()
            bulunan[dizi] = {
                "part_id": f"REF_{slug(etiket)}_{len(dizi)}",
                "name": f"{etiket} ({rec['name']})",
                "category": kat,
                "sequence": dizi,
                "length": len(dizi),
                "accession": (f"{rec['accession']}:{bas}..{son}"
                              + ("" if yon > 0 else " (ters zincir)")),
                "source_vector": rec["name"],
                "feature_type": feat["type"],
                "note": (feat.get("note") or "")[:300],
                "also_in": [],
            }

    # Kimlik çakışmalarını tekilleştir
    sayac: dict[str, int] = {}
    for kayit in bulunan.values():
        pid = kayit["part_id"]
        sayac[pid] = sayac.get(pid, 0) + 1
        if sayac[pid] > 1:
            kayit["part_id"] = f"{pid}_{sayac[pid]}"

    kayitlar = sorted(bulunan.values(),
                      key=lambda k: (k["category"], -len(k["also_in"]),
                                     k["part_id"]))
    OUT.write_text(json.dumps({
        "source": "VectorLab referans vektörlerinin GenBank öznitelikleri",
        "generated": _dt.date.today().isoformat(),
        "note": ("Her parça bir referans vektörün öznitelik tablosundan "
                 "çıkarıldı; accession alanı kaynağı ve koordinatı "
                 "taşır. Katalogda zaten birebir bulunan diziler "
                 "alınmadı."),
        "parts": kayitlar,
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    import collections
    dagilim = collections.Counter(k["category"] for k in kayitlar)
    print(f"{len(kayitlar)} parça çıkarıldı -> {OUT}")
    for kat, n in sorted(dagilim.items(), key=lambda x: -x[1]):
        print(f"   {kat:<12} {n}")
    print(f"\natlananlar: {atlanan}")
    paylasilan = [k for k in kayitlar if k["also_in"]]
    print(f"birden çok vektörde geçen: {len(paylasilan)}")
    for k in paylasilan[:5]:
        print(f"   {k['name'][:46]:<46} +{len(k['also_in'])} vektör")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
