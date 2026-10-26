"""NCBI'dan referans bakteriyel vektörleri çeker.

Kaynak:  NCBI Nucleotide (E-utilities). GenBank kayıtları serbestçe
         kullanılabilir; NCBI veritabanı içeriğine telif kısıtı koymaz.
Çıktı:   data/reference_vectors/<ID>.json  (dizi + GenBank öznitelikleri)

Uygulama çalışırken ağa çıkmaz; bu betik yalnızca veriyi yenilemek için
elle çalıştırılır:
    python tools/fetch_reference_vectors.py
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import requests
from Bio import SeqIO

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli vektor katalogu)"}
OUT = Path(__file__).resolve().parent.parent / "data" / "reference_vectors"
RAW = OUT / "_genbank"

# vector_id -> (aksesyon, beklenen uzunluk). Uzunluk bir sağlama: NCBI farklı
# bir sürüm döndürürse sessizce yanlış plazmitle çalışmayalım.
WANTED: dict[str, tuple[str, int]] = {
    "PUC19":    ("L09137", 2686),
    "PUC18":    ("L09136", 2686),
    "PBR322":   ("J01749", 4361),
    "PACYC184": ("X06403", 4245),
    "PACYC177": ("X06402", 3941),
    "PBAD24":   ("X81837", 4542),
    "PBAD18":   ("X81838", 4613),
    "PBAD33":   ("LC760188", 5352),
    "PTRC99A":  ("U13872", 4176),
    "PBSIIKS":  ("X52327", 2961),
    "PSC101":   ("NC_002056", 9263),
}

KEEP_QUALIFIERS = ("gene", "product", "note", "label", "standard_name",
                   "bound_moiety", "regulatory_class", "mobile_element_type")
SKIP_TYPES = {"source", "gene"}  # gene, CDS ile birebir tekrar ediyor


def fetch_genbank(acc: str) -> str:
    r = requests.get(f"{EUTILS}/efetch.fcgi", headers=HEADERS, timeout=90,
                     params={"db": "nuccore", "id": acc,
                             "rettype": "gbwithparts", "retmode": "text",
                             "tool": "vectorlab"})
    r.raise_for_status()
    if not r.text.startswith("LOCUS"):
        raise RuntimeError(f"{acc}: GenBank kaydı gelmedi: {r.text[:120]!r}")
    return r.text


def to_record(vector_id: str, acc: str, gb_text: str) -> dict:
    from io import StringIO
    rec = SeqIO.read(StringIO(gb_text), "genbank")
    feats = []
    for f in rec.features:
        if f.type in SKIP_TYPES:
            continue
        q = {k: v[0] for k, v in f.qualifiers.items() if k in KEEP_QUALIFIERS}
        label = (q.get("gene") or q.get("product") or q.get("standard_name")
                 or q.get("label") or q.get("note") or "")
        feats.append({
            "type": f.type,
            "start": int(f.location.start) + 1,   # 1 tabanlı, kapsayıcı
            "end": int(f.location.end),
            "strand": int(f.location.strand or 1),
            "label": label,
            "note": q.get("note", ""),
        })
    return {
        "vector_id": vector_id,
        "accession": acc,
        "version": rec.id,
        "definition": (rec.description or "").rstrip("."),
        "organism": rec.annotations.get("organism", ""),
        "genbank_topology": rec.annotations.get("topology", ""),
        "length": len(rec.seq),
        "sequence": str(rec.seq).upper(),
        "features": feats,
        "fetched": _dt.date.today().isoformat(),
        "source": f"NCBI Nucleotide {rec.id}",
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    bad = []
    for vid, (acc, expect) in WANTED.items():
        raw = RAW / f"{acc}.gb"
        if not raw.exists() or raw.stat().st_size < 500:
            raw.write_text(fetch_genbank(acc), encoding="utf-8")
            time.sleep(0.5)
        rec = to_record(vid, acc, raw.read_text(encoding="utf-8"))
        if rec["length"] != expect:
            bad.append(f"{vid}: {rec['length']} bp geldi, {expect} bekleniyordu")
            continue
        (OUT / f"{vid}.json").write_text(
            json.dumps(rec, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"  {vid:<10} {acc:<10} {rec['length']:>6} bp  "
              f"{len(rec['features']):>3} öznitelik  {rec['definition'][:45]}")
        sys.stdout.flush()
    print(f"\n{len(list(OUT.glob('*.json')))} vektör yazıldı -> {OUT}")
    if bad:
        print("UYUŞMAYAN UZUNLUK:\n  " + "\n  ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
