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
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb"
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
    ("CDS_SODA", "sodA", "Lactococcus lactis", "cds",
     "sodA (manganez superoksit dismutaz)"),
    ("CDS_KATE", "katE", "Escherichia coli", "cds",
     "katE (katalaz HPII)"),
    ("CDS_LASR", "lasR", "Pseudomonas aeruginosa", "regulator",
     "LasR (3OC12-HSL algılayıcı düzenleyici)"),
    ("CDS_CQSA", "cqsA", "Vibrio cholerae", "cds",
     "cqsA (CAI-1 otoindükleyici sentazı)"),
    ("CDS_LACZ_LD", "lacZ", "Lactobacillus delbrueckii", "cds",
     "lacZ / beta-galaktosidaz (gıda sınıfı LAB kaynaklı)"),
    ("CDS_BLAC_BL", "penP", "Bacillus licheniformis", "cds",
     "P1A beta-laktamaz (bağırsakta antibiyotik yıkan terapötik)"),
    ("CDS_LSRK", "lsrK", "Escherichia coli", "cds",
     "lsrK (AI-2 kinaz, quorum söndürme)"),
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

# Proteinin TAMAMI değil, UniProt'un kendi anotasyonunda işaretlenmiş bir
# BÖLGESİ kullanılacak genler. Terapötik yapıların çoğu olgun enzimi ya
# da tek bir peptidi ifade eder; ön-pro proteinin bütününü koymak yanlış
# olur. Bölge sınırları UniProt'un öznitelik tablosundan OKUNUR,
# ezberden yazılmaz: (öznitelik türü, açıklamada aranan metin).
# Tam CDS'in çevirisi önce her zamanki gibi doğrulanır; kesme ondan
# sonra yapılır, yani kesilen parçanın doğruluğu tamamın doğruluğuna
# dayanır.
WANTED_REGION = [
    ("CDS_LSS", "gene:lss", "Staphylococcus simulans", "cds",
     "Lizostafin (olgun enzim, S. aureus'a özgü)",
     ("Chain", "Lysostaphin")),
    ("CDS_GLP1", "gene:GCG", "Homo sapiens", "cds",
     "GLP-1 (1-37) (insan proglukagonundan)",
     ("Peptide", "Glucagon-like peptide 1")),
    ("CDS_ELAFIN", "gene:PI3", "Homo sapiens", "cds",
     "Elafin (olgun, elastaz inhibitörü)",
     ("Chain", "Elafin")),
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


def regions(accession: str) -> list[dict[str, Any]]:
    """UniProt kaydının öznitelik tablosundaki bölgeler (1 tabanlı)."""
    r = requests.get(f"{UNIPROT_ENTRY}/{accession}.json", headers=HEADERS,
                     timeout=60)
    r.raise_for_status()
    out = []
    for f in r.json().get("features", []):
        if f["type"] not in ("Signal", "Propeptide", "Chain", "Peptide"):
            continue
        out.append({"type": f["type"],
                    "start": f["location"]["start"]["value"],
                    "end": f["location"]["end"]["value"],
                    "desc": f.get("description", "")})
    return out


def pick_region(bolgeler: list[dict[str, Any]], tur: str,
                desc: str) -> dict[str, Any] | None:
    """İstenen öznitelik; aynı ada birden fazla uyarsa EN KISA olanı.

    UniProt'ta "Glucagon-like peptide 1" hem 92-128 (asıl peptit) hem
    "Glucagon-like peptide 1(7-37)" gibi türevler olarak geçer. Tam
    eşleşmeyi tercih eder, yoksa açıklaması bu metinle BAŞLAYAN en kısa
    kaydı alır - böylece seçim rastlantıya bırakılmaz.
    """
    adaylar = [b for b in bolgeler if b["type"] == tur]
    tam = [b for b in adaylar if b["desc"] == desc]
    if tam:
        return tam[0]
    baslayan = sorted((b for b in adaylar if b["desc"].startswith(desc)),
                      key=lambda b: b["end"] - b["start"])
    return baslayan[0] if baslayan else None


def cut_region(dna: str, protein: str, bolge: dict[str, Any],
               ) -> tuple[str, str, bool]:
    """Doğrulanmış CDS'ten bir amino asit bölgesini kodon sınırından kes.

    Döner (dizi, bolge_proteini, baslangic_eklendi). Kesilen parça tek
    başına bir ORF değildir: başlangıç kodonu yoksa ATG eklenir ve
    sonuna stop kodonu konur. Bu EKLEME bizimdir, kayıttan gelmez; o
    yüzden ayrıca bildirilir ve parçanın notuna yazılır.
    """
    lo, hi = bolge["start"], bolge["end"]
    parca_aa = protein[lo - 1:hi]
    parca_dna = dna[(lo - 1) * 3:hi * 3]
    eklendi = parca_aa[0] != "M"
    if eklendi:
        parca_dna = "ATG" + parca_dna
        parca_aa = "M" + parca_aa
    return parca_dna + "TAA", parca_aa, eklendi


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

    istekler = [(pid, f"gene:{g}", org, kat, ad, True, None)
                for pid, g, org, kat, ad in WANTED]
    istekler += [(pid, f"gene:{g}", org, kat, ad, False, None)
                 for pid, g, org, kat, ad in WANTED_UNREVIEWED]
    istekler += [(pid, sorgu, org, kat, ad, True, bolge)
                 for pid, sorgu, org, kat, ad, bolge in WANTED_REGION]

    for part_id, sorgu, org, kategori, ad, incelenmis, bolge_spec in istekler:
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

        kayit: dict[str, Any] = {
            "part_id": part_id, "name": ad, "category": kategori,
            "sequence": dna, "length": len(dna),
            "length_aa": up["length_aa"], "organism": org,
            "uniprot": up["accession"], "uniprot_entry": up["entry"],
            "full_name": up["full_name"],
            "reviewed": up.get("reviewed", False),
            "nucleotide": f"{up['embl']} / {up['protein_id']}",
        }

        if bolge_spec:
            tur, desc = bolge_spec
            try:
                bolgeler = regions(up["accession"])
            except Exception as exc:
                basarisiz.append(f"{part_id}: öznitelik tablosu alınamadı "
                                 f"{exc}")
                continue
            b = pick_region(bolgeler, tur, desc)
            if not b:
                basarisiz.append(
                    f"{part_id}: UniProt'ta '{tur}: {desc}' bölgesi yok "
                    f"(var olanlar: "
                    + "; ".join(f"{x['type']} {x['start']}-{x['end']} "
                                f"{x['desc']}" for x in bolgeler) + ")")
                continue
            kesilen, parca_aa, eklendi = cut_region(dna, beklenen, b)
            # Kesilen parça da kendi başına doğrulanır.
            if translate(kesilen).rstrip("*") != parca_aa:
                basarisiz.append(f"{part_id}: kesilen bölgenin çevirisi "
                                 f"tutmadı")
                continue
            kayit.update({
                "sequence": kesilen, "length": len(kesilen),
                "length_aa": len(parca_aa),
                "region": f"{tur} {b['start']}-{b['end']}"
                          + (f" ({b['desc']})" if b["desc"] else ""),
                "region_full_aa": up["length_aa"],
                "start_codon_added": eklendi,
            })
            print(f"   {up['accession']}  bölge {b['start']}-{b['end']} "
                  f"= {len(parca_aa)} aa / {len(kesilen)} bp — "
                  f"çeviri doğrulandı"
                  + ("  (ATG eklendi)" if eklendi else ""))
        else:
            print(f"   {up['accession']}  {up['length_aa']} aa / "
                  f"{len(dna)} bp — çeviri doğrulandı"
                  + ("" if up.get("reviewed") else
                     "  (TrEMBL, incelenmemiş)"))

        kayitlar.append(kayit)
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
