"""E. coli K-12 kodon kullanım tablosunu genomdan hesaplar.

Tabloyu ezberden yazmak yerine kaynaktan üretiyoruz: NC_000913.3
(E. coli K-12 MG1655) kaydındaki bütün kodlayan diziler indirilir,
kodonlar sayılır ve CAI için gereken göreli uyarlanabilirlik (w)
değerleri hesaplanır.

CAI, kodon kullanımını "çok ifade edilen genler neye benziyor"
sorusuna göre ölçer. Bu yüzden iki tablo üretilir:

  genome : bütün genlerin kodon kullanımı (karşılaştırma tabanı)
  ref    : çok ifade edilen genlerin kullanımı (CAI'nin ölçütü)

Referans küme gen adlarından seçilir - ribozomal proteinler,
uzama faktörleri, glikoliz enzimleri ve dış zar porinleri. Bu
seçim çıktının içine de yazılır, böylece hangi genlerin ölçüt
alındığı sonradan görülebilir.

    python tools/fetch_codon_usage.py
"""

from __future__ import annotations

import datetime as _dt
import json
import re
import sys
from collections import Counter
from pathlib import Path

import requests
from Bio.Data import CodonTable

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ACCESSION = "NC_000913.3"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli kodon tablosu)"}
OUT = Path(__file__).resolve().parent.parent / "data" / "ecoli_codon_usage.json"
CACHE = Path(__file__).resolve().parent.parent / "data" / "_cache"

# Çok ifade edilen genler: CAI'nin ölçütü bunların kodon tercihidir.
HIGH_EXPRESSION = re.compile(
    r"^(rps[A-Z]|rpl[A-Z]|rpm[A-Z]|tuf[AB]|fusA|tsf|"
    r"eno|gapA|pgk|fbaA|tpiA|gpmA|pykF|"
    r"ompA|ompC|ompF|lpp|"
    r"groL|groS|dnaK|tig|ahpC|"
    r"atpA|atpD|nusA|nusG)$"
)

_TABLE11 = CodonTable.unambiguous_dna_by_id[11]
STOPS = set(_TABLE11.stop_codons)
STARTS = ("ATG", "GTG", "TTG")


def codon_to_aa(codon: str) -> str:
    if codon in STOPS:
        return "*"
    return _TABLE11.forward_table.get(codon, "")


def fetch_cds() -> str:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{ACCESSION}_cds.fna"
    if path.exists() and path.stat().st_size > 1_000_000:
        print(f"önbellekten okunuyor: {path.name}")
        return path.read_text(encoding="utf-8")
    print(f"{ACCESSION} kodlayan dizileri indiriliyor (birkaç MB)…")
    r = requests.get(f"{EUTILS}/efetch.fcgi", headers=HEADERS, timeout=300,
                     params={"db": "nuccore", "id": ACCESSION,
                             "rettype": "fasta_cds_na", "retmode": "text",
                             "tool": "vectorlab"})
    r.raise_for_status()
    if not r.text.startswith(">"):
        raise RuntimeError(f"FASTA gelmedi: {r.text[:120]!r}")
    path.write_text(r.text, encoding="utf-8")
    return r.text


def parse(fasta: str) -> list[tuple[str, str]]:
    """(gen adı, dizi) çiftleri. Bozuk kayıtlar elenir."""
    out: list[tuple[str, str]] = []
    gene, parts = "", []

    def flush():
        if not parts:
            return
        s = "".join(parts).upper()
        if (len(s) >= 90 and len(s) % 3 == 0 and s[:3] in STARTS
                and s[-3:] in STOPS and set(s) <= set("ACGT")):
            out.append((gene, s))

    for line in fasta.splitlines():
        if line.startswith(">"):
            flush()
            parts = []
            m = re.search(r"\[gene=([^\]]+)\]", line)
            gene = m.group(1) if m else ""
        elif line:
            parts.append(line.strip())
    flush()
    return out


def count_codons(seqs: list[str]) -> Counter:
    c: Counter = Counter()
    for s in seqs:
        # Başlama ve durdurma kodonları tercih istatistiğini bozar.
        for i in range(3, len(s) - 3, 3):
            c[s[i:i + 3]] += 1
    return c


def relative_adaptiveness(counts: Counter) -> dict[str, float]:
    """Her amino asit içinde en sık kodona göre oran (CAI'nin w değeri)."""
    by_aa: dict[str, list[str]] = {}
    for codon in _TABLE11.forward_table:
        by_aa.setdefault(codon_to_aa(codon), []).append(codon)
    w: dict[str, float] = {}
    for aa, codons in by_aa.items():
        en_cok = max((counts.get(c, 0) for c in codons), default=0)
        for c in codons:
            # Hiç görülmeyen kodona sıfır vermek CAI'yi çökertir; çok
            # küçük bir taban değer verilir (Sharp & Li'nin yaklaşımı).
            w[c] = (counts.get(c, 0) / en_cok) if en_cok else 0.0
            if w[c] == 0.0:
                w[c] = 0.005
    return w


def main() -> int:
    genes = parse(fetch_cds())
    print(f"{len(genes)} geçerli kodlayan dizi okundu.")
    if len(genes) < 1000:
        print("Beklenenden az gen geldi, çıktı yazılmadı.")
        return 1

    ref = [(g, s) for g, s in genes if HIGH_EXPRESSION.match(g or "")]
    print(f"referans küme: {len(ref)} gen "
          f"({', '.join(sorted(g for g, _ in ref)[:12])}…)")
    if len(ref) < 40:
        print("Referans küme çok küçük, çıktı yazılmadı.")
        return 1

    genome_counts = count_codons([s for _, s in genes])
    ref_counts = count_codons([s for _, s in ref])
    w = relative_adaptiveness(ref_counts)

    toplam_genom = sum(genome_counts.values())
    by_aa: dict[str, list[str]] = {}
    for codon in list(_TABLE11.forward_table) + list(STOPS):
        by_aa.setdefault(codon_to_aa(codon), []).append(codon)

    tablo: dict[str, dict] = {}
    for aa, codons in by_aa.items():
        aa_toplam = sum(genome_counts.get(c, 0) for c in codons) or 1
        for c in codons:
            tablo[c] = {
                "aa": aa,
                "count": genome_counts.get(c, 0),
                "per_1000": round(1000 * genome_counts.get(c, 0)
                                  / max(toplam_genom, 1), 2),
                "fraction": round(genome_counts.get(c, 0) / aa_toplam, 4),
                "w": round(w.get(c, 0.0), 4),
            }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "organism": "Escherichia coli K-12 MG1655",
        "accession": ACCESSION,
        "source": "NCBI Nucleotide, fasta_cds_na",
        "fetched": _dt.date.today().isoformat(),
        "gene_count": len(genes),
        "reference_set": sorted(g for g, _ in ref),
        "reference_note": ("CAI'nin w değerleri bu çok ifade edilen "
                           "genlerden hesaplandı; fraction ve per_1000 "
                           "ise genomun tamamından."),
        "codons": tablo,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"\n{len(tablo)} kodon yazıldı -> {OUT}")

    ornek = sorted(tablo.items(), key=lambda kv: -kv[1]["per_1000"])[:5]
    for c, d in ornek:
        print(f"  {c} {d['aa']}  {d['per_1000']}/1000  w={d['w']}")
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
