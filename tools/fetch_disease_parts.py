"""Hastalık vektörlerinin terapötik genlerini UniProt ve NCBI'dan çeker.

Yol, sinyal peptitlerinde işe yarayan yolun aynısı: UniProt'tan
incelenmiş kaydı bul, EMBL çapraz başvurusundaki protein kimliğiyle
NCBI'dan kodlayan diziyi al, **çevirisinin UniProt'un bildirdiği
proteine birebir eşit olduğunu sına**. Eşit değilse kayıt yazılmaz.

Ezberden tek baz yazılmaz. Bir terapötik gende çerçeve kayması ya da
yanlış başlangıç, haritada görünmeyen ama proteini tamamen bozan bir
hatadır; bu yüzden doğrulama isteğe bağlı değil.

    python tools/fetch_disease_parts.py
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import requests
from Bio.Data import CodonTable

OUT = Path(__file__).resolve().parent.parent / "data" / "disease_parts.json"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli parca katalogu)"}

_T11 = CodonTable.unambiguous_dna_by_id[11]

# (parça kimliği, gen, organizma, kategori, görünen ad)
WANTED = [
    ("CDS_ARGA", "argA", "Escherichia coli", "cds",
     "argA (N-asetilglutamat sentaz)"),
    ("CDS_OXC", "oxc", "Oxalobacter formigenes", "cds",
     "oxc (oksalil-KoA dekarboksilaz)"),
    ("CDS_FRC", "frc", "Oxalobacter formigenes", "cds",
     "frc (formil-KoA transferaz)"),
    ("CDS_MDEA", "mdeA", "Pseudomonas putida", "cds",
     "mdeA / MGL (metiyonin gama-liyaz)"),
    ("CDS_IL10", "IL10", "Homo sapiens", "cds",
     "IL-10 (insan interlökin-10)"),
    ("MRK_THYA_LL", "thyA", "Lactococcus lactis", "marker",
     "thyA (timidilat sentaz, tamamlama markeri)"),
    ("CDS_LASR", "lasR", "Pseudomonas aeruginosa", "regulator",
     "LasR (3OC12-HSL algılayıcı düzenleyici)"),
]

# Protein adıyla aranacaklar (gen adı belirsiz olanlar).
# UniProt'ta "incelenmis" (Swiss-Prot) karsiligi olmayanlar. TrEMBL
# kaydi otomatik olusturulur ve elle denetlenmemistir; dizi yine de
# kaydin kendisinden gelir ve ceviri dogrulamasi uygulanir, ama bu
# fark parcanin notuna yazilir.
WANTED_UNREVIEWED = [
    ("CDS_PYOS5", "pyoS5", "Pseudomonas aeruginosa", "cds",
     "Pyocin S5 (P. aeruginosa'ya özgü bakteriyosin)"),
]


def translate(dna: str) -> str:
    """Bakteriyel çeviri; ilk kodon başlangıç kodonudur (ATG/GTG/TTG)."""
    out = []
    for i in range(0, len(dna) - 2, 3):
        k = dna[i:i + 3]
        if i == 0 and k in ("ATG", "GTG", "TTG"):
            out.append("M")
            continue
        out.append("*" if k in _T11.stop_codons
                   else _T11.forward_table.get(k, "X"))
    return "".join(out)


def uniprot(query: str, reviewed: bool = True) -> dict | None:
    if reviewed:
        query += " AND reviewed:true"
    r = requests.get(UNIPROT, headers=HEADERS, timeout=60,
                     params={"query": query, "format": "json", "size": 5})
    r.raise_for_status()
    for rec in r.json().get("results", []):
        pid = embl = None
        for x in rec.get("uniProtKBCrossReferences", []):
            if x["database"] != "EMBL":
                continue
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            if props.get("ProteinId", "-") != "-":
                pid, embl = props["ProteinId"], x["id"]
                break
        if not pid:
            continue
        ad = (rec.get("proteinDescription", {})
              .get("recommendedName", {}).get("fullName", {}).get("value", ""))
        return {"accession": rec["primaryAccession"],
                "entry": rec.get("uniProtkbId", ""),
                "reviewed": rec.get("entryType", "").endswith("(Swiss-Prot)"),
                "protein": rec["sequence"]["value"],
                "length_aa": rec["sequence"]["length"],
                "full_name": ad, "protein_id": pid, "embl": embl}
    return None


def cds(protein_id: str) -> str | None:
    r = requests.get(f"{EUTILS}/efetch.fcgi", headers=HEADERS, timeout=120,
                     params={"db": "protein", "id": protein_id,
                             "rettype": "fasta_cds_na", "retmode": "text",
                             "tool": "vectorlab"})
    if not r.ok or not r.text.startswith(">"):
        return None
    return "".join(l.strip() for l in r.text.splitlines()[1:]
                   if l and not l.startswith(">")).upper()


def main() -> int:
    kayitlar: list[dict] = []
    basarisiz: list[str] = []

    istekler = [(pid, f"gene:{g}", org, kat, ad, True)
                for pid, g, org, kat, ad in WANTED]
    istekler += [(pid, f"gene:{g}", org, kat, ad, False)
                 for pid, g, org, kat, ad in WANTED_UNREVIEWED]

    for part_id, sorgu, org, kategori, ad, incelenmis in istekler:
        print(f"{ad}…")
        try:
            up = uniprot(f'{sorgu} AND organism_name:"{org}"', incelenmis)
        except Exception as exc:
            basarisiz.append(f"{part_id}: UniProt hatası {exc}")
            continue
        if not up:
            basarisiz.append(f"{part_id}: UniProt'ta uygun kayıt yok")
            continue

        dna = cds(up["protein_id"])
        time.sleep(0.4)
        if not dna:
            basarisiz.append(f"{part_id}: CDS indirilemedi "
                             f"({up['protein_id']})")
            continue

        cevrilen = translate(dna).rstrip("*")
        beklenen = up["protein"]
        if cevrilen != beklenen:
            basarisiz.append(
                f"{part_id}: ÇEVİRİ TUTMADI ({len(cevrilen)} aa vs "
                f"{len(beklenen)} aa)")
            continue

        kayitlar.append({
            "part_id": part_id, "name": ad, "category": kategori,
            "sequence": dna, "length": len(dna),
            "length_aa": up["length_aa"], "organism": org,
            "uniprot": up["accession"], "uniprot_entry": up["entry"],
            "full_name": up["full_name"],
            "reviewed": up.get("reviewed", False),
            "nucleotide": f"{up['embl']} / {up['protein_id']}",
        })
        print(f"   {up['accession']}  {up['length_aa']} aa / {len(dna)} bp "
              f"— çeviri doğrulandı"
              + ("" if up.get("reviewed") else "  (TrEMBL, incelenmemiş)"))
        sys.stdout.flush()

    OUT.write_text(json.dumps({
        "source": "UniProt (protein) + NCBI (kodlayan dizi)",
        "generated": _dt.date.today().isoformat(),
        "note": ("Her kaydın nükleotid dizisi çevrilip UniProt'un "
                 "bildirdiği proteine eşit olduğu doğrulandı; eşit "
                 "olmayan kayıt yazılmadı."),
        "parts": kayitlar,
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"\n{len(kayitlar)} gen yazıldı -> {OUT}")
    if basarisiz:
        print("Alınamayanlar:")
        for b in basarisiz:
            print(f"   {b}")
    return 0 if kayitlar else 1


if __name__ == "__main__":
    raise SystemExit(main())
