"""Sinyal peptidi dizilerini UniProt ve NCBI'dan çeker.

Sinyal peptidi, proteini hücre içinde nereye gideceğine yönlendiren ve
sonra kesilip atılan N-ucu dizisidir. Önemli olan ayrıntı, **nerede
kesildiğidir**: bir baz fazla ya da eksik alırsanız okuma çerçevesi
kayar ya da olgun proteinin başına birkaç yabancı amino asit eklenir.

Bu yüzden kesim noktası ezberden yazılmıyor: UniProt'un `Signal`
öznitelik aralığı okunuyor (deneysel kanıt ya da SignalP tahmini ile
işaretlenmiş), karşılık gelen nükleotid dizisi NCBI'dan alınıyor ve
**çevirisi UniProt'un bildirdiği peptide birebir eşit mi diye
sınanıyor**. Eşit değilse kayıt yazılmıyor.

    python tools/fetch_signal_peptides.py
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import requests
from Bio.Data import CodonTable

OUT = Path(__file__).resolve().parent.parent / "data" / "signal_peptides.json"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli parca katalogu)"}

_T11 = CodonTable.unambiguous_dna_by_id[11]

# (gen, organizma, kısa ad, yol, not)
WANTED = [
    ("pelB", "Dickeya chrysanthemi", "pelB",
     "Sec (post-translasyonel)",
     "Periplazmik hedeflemenin en yaygın sinyali; pET-22b bu sinyali "
     "taşır. Protein sitoplazmada katlanmadan, açık hâlde zara "
     "taşınır."),
    ("ompA", "Escherichia coli", "OmpA",
     "Sec (post-translasyonel)",
     "E. coli'nin kendi dış zar proteininin sinyali. pelB'ye "
     "alternatiftir ve çoğu hedefte benzer verim verir."),
    ("phoA", "Escherichia coli", "PhoA",
     "Sec (post-translasyonel)",
     "Alkalen fosfataz sinyali. Klasik periplazmik hedefleme "
     "çalışmalarının çoğu bununla yapılmıştır."),
    ("malE", "Escherichia coli", "MalE",
     "Sec (post-translasyonel)",
     "Maltoz bağlayan proteinin kendi sinyali. MBP füzyonunu "
     "periplazmaya göndermek isterseniz bu sinyal kullanılır."),
    ("dsbA", "Escherichia coli", "DsbA",
     "SRP (ko-translasyonel)",
     "Diğerlerinden farklı yol kullanır: protein daha ribozomdan "
     "çıkarken zara yönlendirilir. Sitoplazmada toplanıp çökelen "
     "hedeflerde bu yüzden daha iyi sonuç verebilir."),
    ("torA", "Escherichia coli", "TorA",
     "Tat (twin-arginine)",
     "Tat yolu, KATLANMIŞ proteini taşır - ötekiler açık zincir "
     "taşır. Kofaktör bağlamış ya da sitoplazmada katlanması gereken "
     "proteinler için tek seçenektir. Buna karşılık kapasitesi "
     "düşüktür, yüksek verim beklemeyin."),
    ("sacB", "Bacillus subtilis", "SacB",
     "Sec (Gram-pozitif)",
     "B. subtilis levansukraz sinyali. Gram-pozitiflerde dış zar "
     "olmadığı için protein periplazmada sıkışmaz, doğrudan "
     "besiyerine salgılanır."),
    ("amyE", "Bacillus subtilis", "AmyE",
     "Sec (Gram-pozitif)",
     "B. subtilis alfa-amilaz sinyali; endüstriyel enzim üretiminde "
     "en çok kullanılan salgılama sinyallerinden biridir."),
]


def translate(dna: str) -> str:
    """Bakteriyel çeviri: ilk kodon başlangıç kodonudur.

    Bakterilerde genler ATG dışında GTG ve TTG ile de başlar ve o
    konumda hepsi formil-metiyonin olarak okunur - `phoA` örneğin GTG
    ile başlar. Baştaki kodonu normal kodon gibi çevirirsek valin
    çıkar ve UniProt'un bildirdiği peptide eşleşmez.
    """
    out = []
    for i in range(0, len(dna) - 2, 3):
        k = dna[i:i + 3]
        if i == 0 and k in ("ATG", "GTG", "TTG"):
            out.append("M")
            continue
        out.append("*" if k in _T11.stop_codons
                   else _T11.forward_table.get(k, "X"))
    return "".join(out)


def uniprot_record(gene: str, organism: str) -> dict | None:
    q = (f'gene:{gene} AND organism_name:"{organism}" AND reviewed:true')
    r = requests.get(f"{UNIPROT}/search", headers=HEADERS, timeout=60,
                     params={"query": q, "format": "json", "size": 5})
    r.raise_for_status()
    for rec in r.json().get("results", []):
        sinyal = [f for f in rec.get("features", [])
                  if f.get("type") == "Signal"]
        if not sinyal:
            continue
        loc = sinyal[0]["location"]
        bas, son = loc["start"]["value"], loc["end"]["value"]
        if not bas or not son:
            continue
        protein = rec["sequence"]["value"]
        pid = None
        embl = None
        for x in rec.get("uniProtKBCrossReferences", []):
            if x["database"] != "EMBL":
                continue
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            if props.get("ProteinId", "-") != "-":
                pid, embl = props["ProteinId"], x["id"]
                break
        return {
            "accession": rec["primaryAccession"],
            "entry": rec.get("uniProtkbId", ""),
            "signal_start": bas, "signal_end": son,
            "peptide": protein[bas - 1:son],
            "protein_id": pid, "embl": embl,
            "evidence": sinyal[0].get("evidences", [{}])[0]
                        .get("evidenceCode", ""),
        }
    return None


def cds_nucleotide(protein_id: str) -> str | None:
    r = requests.get(f"{EUTILS}/efetch.fcgi", headers=HEADERS, timeout=90,
                     params={"db": "protein", "id": protein_id,
                             "rettype": "fasta_cds_na", "retmode": "text",
                             "tool": "vectorlab"})
    if not r.ok or not r.text.startswith(">"):
        return None
    satirlar = r.text.splitlines()
    return "".join(l.strip() for l in satirlar[1:] if l and not
                   l.startswith(">")).upper()


def main() -> int:
    kayitlar = []
    basarisiz = []
    for gen, org, kisa, yol, aciklama in WANTED:
        print(f"{kisa} ({gen}, {org})…")
        try:
            up = uniprot_record(gen, org)
        except Exception as exc:
            basarisiz.append(f"{kisa}: UniProt hatası {exc}")
            continue
        if not up:
            basarisiz.append(f"{kisa}: UniProt'ta sinyal peptidi yok")
            continue
        if not up["protein_id"]:
            basarisiz.append(f"{kisa}: nükleotid çapraz başvurusu yok")
            continue

        cds = cds_nucleotide(up["protein_id"])
        time.sleep(0.4)
        if not cds:
            basarisiz.append(f"{kisa}: CDS indirilemedi "
                             f"({up['protein_id']})")
            continue

        n = len(up["peptide"]) * 3
        dna = cds[:n]
        cevrilen = translate(dna)
        if cevrilen != up["peptide"]:
            basarisiz.append(
                f"{kisa}: ÇEVİRİ TUTMADI - DNA {cevrilen!r} vs "
                f"UniProt {up['peptide']!r}")
            continue

        kayitlar.append({
            "part_id": f"SIG_{kisa.upper()}",
            "name": f"{kisa} sinyal peptidi",
            "peptide": up["peptide"],
            "sequence": dna,
            "length": len(dna),
            "length_aa": len(up["peptide"]),
            "pathway": yol,
            "organism": org,
            "gene": gen,
            "uniprot": up["accession"],
            "uniprot_entry": up["entry"],
            "nucleotide": f"{up['embl']} / {up['protein_id']}",
            "note": aciklama,
        })
        print(f"   {up['accession']} {up['peptide']} "
              f"({len(up['peptide'])} aa) — çeviri doğrulandı")
        sys.stdout.flush()

    OUT.write_text(json.dumps({
        "source": "UniProt (sinyal peptidi aralığı) + NCBI (nükleotid)",
        "generated": _dt.date.today().isoformat(),
        "note": ("Her kaydın nükleotid dizisi çevrilip UniProt'un "
                 "bildirdiği peptide eşit olduğu doğrulandı; eşit "
                 "olmayan kayıt yazılmadı."),
        "peptides": kayitlar,
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"\n{len(kayitlar)} sinyal peptidi yazıldı -> {OUT}")
    if basarisiz:
        print("Alınamayanlar:")
        for b in basarisiz:
            print(f"   {b}")
    return 0 if kayitlar else 1


if __name__ == "__main__":
    raise SystemExit(main())
