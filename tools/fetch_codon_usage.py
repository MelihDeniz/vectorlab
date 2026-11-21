"""Konak genomlarından kodon kullanım tablosu hesaplar.

Tabloları ezberden yazmak yerine kaynaktan üretiyoruz: her konağın
RefSeq genom kaydındaki bütün kodlayan diziler indirilir, kodonlar
sayılır ve CAI için gereken göreli uyarlanabilirlik (w) değerleri
hesaplanır.

Katalogdaki beş organizma için ayrı tablo üretilir. Bu tek tablonun
yetmemesinin sebebi basit: kodon tercihi organizmaya göre değişir, ve
uygulama hem E. coli'de hem *L. cremoris*'te terapötik tasarlıyor. Bir
*L. cremoris* yapısını E. coli tablosuyla optimize etmek, optimize
etmemekten daha iyi değildir - yanlış yöne çeker.

CAI, kodon kullanımını "çok ifade edilen genler neye benziyor"
sorusuna göre ölçer. Bu yüzden her organizma için iki tablo üretilir:

  genome : bütün genlerin kodon kullanımı (karşılaştırma tabanı)
  ref    : çok ifade edilen genlerin kullanımı (CAI'nin ölçütü)

**Referans küme nasıl seçiliyor.** Gen adlarına göre seçmek beş
organizmada işe yaramaz: ribozomal protein genleri her yerde `rpl/rps`
diye adlanır ama glikoliz genleri adlanmaz (`gapA` / `gap` / `gapB`),
ve bazı kayıtlarda gen adı hiç yoktur, yalnızca locus_tag bulunur. Bu
yüzden seçim **proteinin tanımına** göre yapılıyor: RefSeq'in
`[protein=...]` alanı beş genomda da standart bir sözcük kullanıyor
("50S ribosomal protein L2", "elongation factor Tu"). Böylece beş
organizmanın hepsine AYNI kural uygulanıyor; kural organizmaya göre
elle ayarlanmış olsaydı tablolar birbiriyle karşılaştırılamazdı.

Kümeye giren protein sınıfları, klasik CAI referans kümesindekiler
(Sharp & Li): ribozomal proteinler, uzama faktörleri, şaperonlar,
glikoliz enzimleri ve ATP sentaz altbirimleri. Hangi kalıbın kaç gen
getirdiği çıktının içine yazılır, yani seçim sonradan denetlenebilir.

    python tools/fetch_codon_usage.py            # beş organizmanın hepsi
    python tools/fetch_codon_usage.py ecoli      # yalnızca biri
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
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli kodon tablosu)"}
OUT_DIR = Path(__file__).resolve().parent.parent / "data" / "codon_usage"
CACHE = Path(__file__).resolve().parent.parent / "data" / "_cache"

#: Çok ifade edilen protein sınıfları, proteinin TANIMINA göre.
#:
#: Kalıplar tanımın tamamına bağlanmıştır (`^...$`). Serbest bırakılsaydı
#: "ribosomal protein L11 methyltransferase" gibi düzenleyici enzimler de
#: kümeye girerdi; onlar çok ifade edilen gen değildir ve ölçütü
#: bozarlardı. Aynı şekilde "transcription elongation factor GreB", bir
#: uzama faktörü adı taşımasına rağmen kümeye girmez.
#:
#: Alternatifler beş genomun kayıtlarından OKUNARAK yazıldı, tahmin
#: edilmedi: aynı protein beş kayıtta farklı adlanıyor. Ribozomal
#: protein L4 örneğin E. coli'de "50S ribosomal subunit protein L4",
#: B. subtilis'te "ribosomal protein L4", P. putida'da "50S ribosomal
#: protein L4", C. glutamicum'da "ribosomal protein bL12" (yeni
#: b/u önekli adlandırma) biçiminde geçiyor.
REFERENCE_CLASSES = {
    "ribozomal protein":
        r"(?:30S |50S )?ribosomal (?:subunit )?protein [bu]?[LS]\d+[A-Za-z]?",
    "uzama faktörü":
        r"(?:translation )?elongation factor (?:Tu|G|Ts|P)(?: \d)?",
    "şaperon":
        r"(?:co)?chaperonin(?: (?:GroEL|GroES|large subunit|small subunit))?"
        r"|(?:molecular chaperone |chaperone protein )DnaK"
        r"|trigger factor|prolyl isomerase \(trigger factor\)",
    "glikoliz":
        r"enolase|phosphoglycerate kinase|triose-?phosphate isomerase"
        r"|fructose-bisphosphate aldolase(?: class (?:I|II))?"
        r"|(?:type I )?glyceraldehyde-3-phosphate dehydrogenase"
        r"|pyruvate kinase",
    "ATP sentaz":
        r"(?:F0F1 )?ATP synthase[^$]*subunit \(?(?:alpha|beta)\)?",
}

_CLASS_RE = {ad: re.compile(r"^(?:" + kalip + r")$", re.IGNORECASE)
             for ad, kalip in REFERENCE_CLASSES.items()}

#: Korunmuş gen adları. Protein tanımı tek başına yetmiyor çünkü bazı
#: kayıtlarda tanım organizmaya özgü yazılmış; gen adı tek başına da
#: yetmiyor çünkü bu beş kayıtta gen adı oranı %14 ile %100 arasında
#: değişiyor (P. putida 1372/5616, L. cremoris 597/2502). İki kural
#: BİRLEŞTİRİLİYOR: ribozomal protein ve uzama faktörü genleri, adı olan
#: genlerin arasında neredeyse her zaman bulunur, çünkü adlandırması
#: bakterilerde standarttır (rplA-rplY, rpsA-rpsU, rpmA-rpmJ).
GENE_PATTERN = re.compile(
    r"^(?:rps[A-Z]\d?|rpl[A-Z]\d?|rpm[A-Z]|tuf[AB]?|fusA|tsf"
    r"|eno|gap[AB]?|pgk|tpiA?|fba[AB]?|gpmA|pgm|pyk[AF]?"
    r"|groE?[LS]|dnaK|tig|atpA|atpD)$", re.IGNORECASE)

#: Katalogdaki konakların organizmaları. `accession` NCBI'ya sorularak
#: teyit edildi (bkz. `organism_ncbi`): ezberden yazılmadı.
GENOMES = [
    {
        "key": "ecoli",
        "accession": "NC_000913.3",
        "organism": "Escherichia coli K-12 MG1655",
        "organism_ncbi": "Escherichia coli str. K-12 substr. MG1655",
        "catalog_organism": "E. coli",
        "min_genes": 3000,
        "taxonomy_note": "",
    },
    {
        "key": "bsubtilis",
        "accession": "NC_000964.3",
        "organism": "Bacillus subtilis subsp. subtilis str. 168",
        "organism_ncbi": "Bacillus subtilis subsp. subtilis str. 168",
        "catalog_organism": "Bacillus subtilis",
        "min_genes": 3000,
        "taxonomy_note": "",
    },
    {
        "key": "pputida",
        "accession": "NC_002947.4",
        "organism": "Pseudomonas putida KT2440",
        "organism_ncbi": "Pseudomonas putida KT2440",
        "catalog_organism": "Pseudomonas putida",
        "min_genes": 4000,
        "taxonomy_note": "",
    },
    {
        "key": "llactis",
        "accession": "NC_009004.1",
        "organism": "Lactococcus cremoris subsp. cremoris MG1363",
        "organism_ncbi": "Lactococcus cremoris subsp. cremoris MG1363",
        "catalog_organism": "Lactococcus lactis",
        # 2,5 Mb'lik küçük bir genom; E. coli'nin yarısı kadar gen taşır.
        "min_genes": 2000,
        "taxonomy_note": (
            "Katalogdaki konak 'Lactococcus lactis NZ9000' olarak "
            "adlanır; NZ9000 MG1363 türevidir. NCBI bu suşu artık "
            "*Lactococcus cremoris* subsp. cremoris olarak kaydediyor: "
            "L. lactis subsp. cremoris 2021'de ayrı bir tür mertebesine "
            "çıkarıldı. Tablo MG1363'ün kendi genomundan hesaplandığı "
            "için adlandırmadaki bu değişiklik sayıları etkilemez, ama "
            "kaynağı ararken karşınıza iki ad çıkacağı için yazıldı."
        ),
    },
    {
        "key": "cglutamicum",
        "accession": "NC_003450.3",
        "organism": "Corynebacterium glutamicum ATCC 13032",
        "organism_ncbi": "Corynebacterium glutamicum ATCC 13032",
        "catalog_organism": "Corynebacterium glutamicum",
        "min_genes": 2500,
        "taxonomy_note": "",
    },
]

_TABLE11 = CodonTable.unambiguous_dna_by_id[11]
STOPS = set(_TABLE11.stop_codons)
STARTS = ("ATG", "GTG", "TTG")


def codon_to_aa(codon: str) -> str:
    if codon in STOPS:
        return "*"
    return _TABLE11.forward_table.get(codon, "")


def fetch_cds(accession: str) -> str:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{accession}_cds.fna"
    if path.exists() and path.stat().st_size > 500_000:
        print(f"  önbellekten okunuyor: {path.name}")
        return path.read_text(encoding="utf-8")
    print(f"  {accession} kodlayan dizileri indiriliyor (birkaç MB)…")
    r = requests.get(f"{EUTILS}/efetch.fcgi", headers=HEADERS, timeout=300,
                     params={"db": "nuccore", "id": accession,
                             "rettype": "fasta_cds_na", "retmode": "text",
                             "tool": "vectorlab"})
    r.raise_for_status()
    if not r.text.startswith(">"):
        raise RuntimeError(f"FASTA gelmedi: {r.text[:120]!r}")
    path.write_text(r.text, encoding="utf-8")
    return r.text


_GENE_RE = re.compile(r"\[gene=([^\]]+)\]")
_PROTEIN_RE = re.compile(r"\[protein=([^\]]+)\]")
_LOCUS_RE = re.compile(r"\[locus_tag=([^\]]+)\]")


def parse(fasta: str) -> list[dict[str, str]]:
    """Her kayıt için ad, protein tanımı ve dizi. Bozuk kayıtlar elenir."""
    out: list[dict[str, str]] = []
    bilgi: dict[str, str] = {}
    parts: list[str] = []

    def flush():
        if not parts or not bilgi:
            return
        s = "".join(parts).upper()
        if (len(s) >= 90 and len(s) % 3 == 0 and s[:3] in STARTS
                and s[-3:] in STOPS and set(s) <= set("ACGT")):
            out.append({**bilgi, "seq": s})

    for line in fasta.splitlines():
        if line.startswith(">"):
            flush()
            parts = []
            g = _GENE_RE.search(line)
            p = _PROTEIN_RE.search(line)
            lt = _LOCUS_RE.search(line)
            bilgi = {
                "gene": g.group(1) if g else "",
                "protein": p.group(1) if p else "",
                "locus": lt.group(1) if lt else "",
            }
        elif line:
            parts.append(line.strip())
    flush()
    return out


def pick_reference(genes: list[dict[str, str]]
                   ) -> tuple[list[dict[str, str]], dict[str, int]]:
    """Çok ifade edilen genleri seç: protein tanımı VEYA korunmuş gen adı.

    Dönen sayaç hangi kuralın kaç gen getirdiğini söyler; çıktıya
    yazılır ki seçim sonradan denetlenebilsin. İki kuralın birleşimi
    alınır, kesişimi değil: bir gen ya tanımından ya adından tanınırsa
    kümeye girer.
    """
    secilen: list[dict[str, str]] = []
    sayac: dict[str, int] = {ad: 0 for ad in _CLASS_RE}
    sayac["gen adı (korunmuş)"] = 0
    for g in genes:
        tanim = (g.get("protein") or "").strip()
        ad = (g.get("gene") or "").strip()
        bulundu = False
        for sinif, rx in _CLASS_RE.items():
            if rx.match(tanim):
                sayac[sinif] += 1
                bulundu = True
                break
        if not bulundu and ad and GENE_PATTERN.match(ad):
            sayac["gen adı (korunmuş)"] += 1
            bulundu = True
        if bulundu:
            secilen.append(g)
    return secilen, sayac


def cai_of(seq: str, w: dict[str, float],
           syn: dict[str, list[str]]) -> float | None:
    """CAI, `core.codon.cai` ile AYNI kuralla.

    Başlama ve durdurma kodonları ile tek seçenekli amino asitler
    dışarıda bırakılır. Burada ayrıca hesaplanmasının sebebi, "iyi bir
    CAI kaçtır" sorusunun cevabının organizmaya göre değişmesi: ölçek
    her tablonun kendi referans genlerinden çıkarılıyor.
    """
    import math

    logs = []
    for i in range(3, len(seq), 3):
        c = seq[i:i + 3]
        if len(c) < 3 or c in STOPS:
            continue
        aa = _TABLE11.forward_table.get(c)
        if not aa or len(syn.get(aa, [])) < 2:
            continue
        logs.append(math.log(max(w.get(c, 0.005), 1e-4)))
    if not logs:
        return None
    return math.exp(sum(logs) / len(logs))


def reference_scale(ref: list[dict[str, str]],
                    w: dict[str, float]) -> dict[str, float]:
    """Referans genlerin CAI dağılımı: "iyi bir CAI kaçtır"ın cevabı.

    Eskiden bu iki sayı (E. coli `rpsA` ve `tufA`) koda gömülüydü. Beş
    organizma olunca gömülü sayı yanlış oldu: her tablonun ölçeği
    kendi genlerinden çıkmak zorunda.
    """
    syn: dict[str, list[str]] = {}
    for codon, aa in _TABLE11.forward_table.items():
        syn.setdefault(aa, []).append(codon)

    degerler = sorted(x for x in (cai_of(g["seq"], w, syn) for g in ref)
                      if x is not None)
    if not degerler:
        return {}

    def yuzdelik(p: float) -> float:
        i = min(int(p * (len(degerler) - 1)), len(degerler) - 1)
        return round(degerler[i], 4)

    return {
        "gene_count": len(degerler),
        "median": yuzdelik(0.5),
        "p10": yuzdelik(0.10),
        "p90": yuzdelik(0.90),
    }


def count_codons(seqs: list[str]) -> Counter:
    c: Counter = Counter()
    for s in seqs:
        # Başlama ve durdurma kodonları tercih istatistiğini bozar.
        for i in range(3, len(s) - 3, 3):
            c[s[i:i + 3]] += 1
    return c


def count_terminal_stops(seqs: list[str]) -> Counter:
    """Genlerin GERÇEK durdurma kodonu.

    `count_codons` ilk ve son kodonu atlar, yani durdurma kodonlarını
    hiç saymaz. Buna rağmen tabloda sıfırdan farklı durdurma sayıları
    görünüyordu: onlar iç çerçevede durdurma kodonu taşıyan genlerden
    (selenoproteinler, okumadan geçme) geliyordu ve gerçek kullanımı
    temsil etmiyordu - E. coli tablosu bu yüzden "TGA, durdurma
    kodonlarının %100'ü" diyordu, oysa E. coli'de baskın durdurma
    kodonu TAA'dır.

    Hangi durdurma kodonunun kullanılacağı gerçek bir tasarım kararı
    olduğu için burada doğru sayılıyor: her genin SON kodonu.
    """
    return Counter(s[-3:] for s in seqs if s[-3:] in STOPS)


def count_internal_stops(seqs: list[str]) -> int:
    """İç çerçevede durdurma kodonu taşıyan gen sayısı (veri kalitesi)."""
    n = 0
    for s in seqs:
        if any(s[i:i + 3] in STOPS for i in range(3, len(s) - 3, 3)):
            n += 1
    return n


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


def build(cfg: dict) -> int:
    print(f"\n=== {cfg['organism']} ({cfg['accession']}) ===")
    genes = parse(fetch_cds(cfg["accession"]))
    adli = sum(1 for g in genes if g["gene"])
    tanimli = sum(1 for g in genes if g["protein"])
    print(f"  {len(genes)} geçerli kodlayan dizi "
          f"({adli} tanesinde gen adı, {tanimli} tanesinde protein tanımı)")
    if len(genes) < cfg["min_genes"]:
        print(f"  Beklenenden az gen geldi (< {cfg['min_genes']}), "
              f"çıktı yazılmadı.")
        return 1

    ref, sayac = pick_reference(genes)
    print(f"  referans küme: {len(ref)} gen")
    for ad, n in sayac.items():
        print(f"    {ad:20} {n:>4}")
    if len(ref) < 40:
        print("  Referans küme çok küçük, çıktı yazılmadı.")
        return 1

    genome_counts = count_codons([g["seq"] for g in genes])
    ref_counts = count_codons([g["seq"] for g in ref])
    w = relative_adaptiveness(ref_counts)

    toplam_genom = sum(genome_counts.values())
    by_aa: dict[str, list[str]] = {}
    for codon in list(_TABLE11.forward_table) + list(STOPS):
        by_aa.setdefault(codon_to_aa(codon), []).append(codon)

    # Durdurma kodonları ayrı sayılır: genin son kodonu.
    stop_counts = count_terminal_stops([g["seq"] for g in genes])
    ic_durdurma = count_internal_stops([g["seq"] for g in genes])

    tablo: dict[str, dict] = {}
    for aa, codons in by_aa.items():
        durdurma = aa == "*"
        sayim = stop_counts if durdurma else genome_counts
        aa_toplam = sum(sayim.get(c, 0) for c in codons) or 1
        for c in codons:
            tablo[c] = {
                "aa": aa,
                "count": sayim.get(c, 0),
                # Durdurma kodonları anlamlı kodonların 1000'de-bir
                # istatistiğine girmez; aynı denominatörle karışmasın.
                "per_1000": 0.0 if durdurma else round(
                    1000 * genome_counts.get(c, 0) / max(toplam_genom, 1), 2),
                "fraction": round(sayim.get(c, 0) / aa_toplam, 4),
                # CAI'de durdurma kodonunun seçimi yoktur.
                "w": 0.0 if durdurma else round(w.get(c, 0.0), 4),
            }

    baskin_stop = max(STOPS, key=lambda c: stop_counts.get(c, 0))
    print(f"  durdurma kodonu: " + ", ".join(
        f"{c} %{round(100 * stop_counts.get(c, 0) / max(sum(stop_counts.values()), 1))}"
        for c in sorted(STOPS, key=lambda x: -stop_counts.get(x, 0)))
        + f"  (baskın: {baskin_stop})")

    # Referans kümedeki genlerin adı; adı olmayanlar locus_tag ile yazılır.
    ref_adlari = sorted({(g["gene"] or g["locus"] or "?") for g in ref})

    olcek = reference_scale(ref, w)
    if olcek:
        print(f"  referans CAI: ortanca {olcek['median']} "
              f"(%10-%90: {olcek['p10']}-{olcek['p90']})")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f"{cfg['key']}.json"
    out.write_text(json.dumps({
        "key": cfg["key"],
        "organism": cfg["organism"],
        "organism_ncbi": cfg["organism_ncbi"],
        "catalog_organism": cfg["catalog_organism"],
        "accession": cfg["accession"],
        "source": "NCBI Nucleotide, fasta_cds_na",
        "fetched": _dt.date.today().isoformat(),
        "gene_count": len(genes),
        "taxonomy_note": cfg["taxonomy_note"],
        "reference_rule": (
            "Çok ifade edilen genler iki kuralın BİRLEŞİMİYLE seçildi ve "
            "beş organizmanın hepsine aynı kurallar uygulandı: (1) RefSeq'in "
            "[protein=...] tanımı şu sınıflardan birine tam olarak uyuyorsa: "
            + ", ".join(REFERENCE_CLASSES) + "; (2) ya da gen adı korunmuş "
            "ribozomal protein / uzama faktörü / glikoliz / şaperon "
            "adlandırmasına uyuyorsa. İkinci kural gerekli, çünkü bu beş "
            "kayıtta gen adı taşıyan CDS oranı %14 ile %100 arasında "
            "değişiyor; birinci kural gerekli, çünkü gen adı olmayan "
            "kayıtlarda tanınacak tek şey tanımdır."
        ),
        "reference_classes": sayac,
        "stop_codon_note": (
            "Durdurma kodonlarının `count` ve `fraction` değerleri her "
            "genin SON kodonundan sayıldı; anlamlı kodonların "
            "istatistiğine girmedikleri için `per_1000` ve `w` sıfırdır. "
            f"Bu kayıtta {ic_durdurma} gen iç çerçevede durdurma kodonu "
            "taşıyor (selenoproteinler ve okumadan geçme); onlar "
            "sayılmadı, çünkü gerçek durdurma kullanımını temsil "
            "etmiyorlar."
        ),
        "internal_stop_genes": ic_durdurma,
        "reference_cai": olcek,
        "reference_cai_note": (
            "Referans kümedeki genlerin kendi tablosuyla hesaplanan CAI "
            "dağılımı. 'İyi bir CAI kaçtır' sorusunun cevabı budur ve "
            "organizmaya göre değişir; optimizasyonun 1,0'a gitmemesi "
            "gerektiğini gösteren ölçek de bu."
        ),
        "reference_set": ref_adlari,
        "reference_note": ("CAI'nin w değerleri bu çok ifade edilen "
                           "genlerden hesaplandı; fraction ve per_1000 "
                           "ise genomun tamamından."),
        "codons": tablo,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"  {len(tablo)} kodon yazıldı -> {out.name}")

    ornek = sorted(tablo.items(), key=lambda kv: -kv[1]["per_1000"])[:3]
    for c, d in ornek:
        print(f"    {c} {d['aa']}  {d['per_1000']}/1000  w={d['w']}")
    return 0


def main(argv: list[str]) -> int:
    istenen = [a for a in argv if not a.startswith("-")]
    seçilmiş = [g for g in GENOMES if not istenen or g["key"] in istenen]
    if not seçilmiş:
        print(f"Bilinmeyen anahtar. Seçenekler: "
              f"{', '.join(g['key'] for g in GENOMES)}")
        return 1
    hata = 0
    for cfg in seçilmiş:
        hata |= build(cfg)
    sys.stdout.flush()
    return hata


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
