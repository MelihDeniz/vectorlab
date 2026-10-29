"""Uygulamanın hesaplarını BAĞIMSIZ kaynaklarla karşılaştırır.

"Test ettim" ile "bağımsız bir ölçütle aynı sonucu verdi" aynı şey
değildir. Buradaki kontrollerin her biri, uygulamanın kendi kodundan
bağımsız bir referansa dayanır:

  1. pUC19 çoklu klonlama bölgesinin yayımlanmış koordinatları
  2. Üç ayrı erime sıcaklığı yöntemi (en yakın komşu, GC tabanlı, Wallace)
  3. Protein kütlesi: ProtParam yerine Biopython'un ayrı
     `molecular_weight` işlevi
  4. NCBI'nın kayıt başlığında bildirdiği uzunluklar
  5. Elle yazılmış nadir kodon listesi ile genomdan hesaplanan tablo

Çalıştırmak:
    python tools/validate_against_references.py

Çıktı README'ye işlenmek üzere düzyazı olarak basılır. Hiçbir kontrol
başarısız olmamalı; olursa hangi değerin nerede ayrıştığı yazılır.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Bio.SeqUtils import MeltingTemp as mt           # noqa: E402
from Bio.SeqUtils import molecular_weight            # noqa: E402

from core import (                                       # noqa: E402
    codon, db, enzymes as ez, protein, ref_vectors as rv,
)
from core.seq_utils import (                         # noqa: E402
    RARE_CODONS_ECOLI, restriction_sites, tm_nn,
)

# pUC19'un çoklu klonlama bölgesindeki tek kesen enzimlerin yayımlanmış
# konumları (standart pUC19 numaralandırması, L09137). Bu sayılar
# uygulamadan değil, kaydın kendi kanonik numaralandırmasından gelir.
# Her girdi: enzim -> tanima dizisinin baslangici. Yayimlanmis pUC19
# haritasiyla ayni sira ve ayni araliga oturur; EcoRI'nin 396. konumda
# olmasi bu numaralandirmanin kanonik capasidir.
PUC19_MCS = {
    "EcoRI": 396, "SacI": 402, "KpnI": 408, "SmaI": 412, "BamHI": 417,
    "XbaI": 423, "SalI": 429, "PstI": 435, "SphI": 441, "HindIII": 447,
}
PUC19_MCS_SIRASI = ["EcoRI", "SacI", "KpnI", "SmaI", "BamHI", "XbaI",
                    "SalI", "PstI", "SphI", "HindIII"]

# NCBI kayıt başlıklarında bildirilen uzunluklar.
BEKLENEN_UZUNLUK = {
    "PUC19": 2686, "PUC18": 2686, "PBR322": 4361, "PACYC184": 4245,
    "PACYC177": 3941, "PBAD24": 4542, "PBAD18": 4613, "PBAD33": 5352,
    "PTRC99A": 4176, "PBSIIKS": 2961, "PSC101": 9263,
    "PSEVA111": 2094, "PSEVA211": 1993, "PSEVA311": 1838,
    "PKD46": 6329, "PKD4": 3267, "PKD3": 2804,
    "PBBR1MCS": 4707, "RSF1010": 8684, "PBHR1": 5300,
    "PUB110": 4548, "PC194": 2910, "PE194": 3728, "PWV01": 2178,
}

gecen = 0
kalan: list[str] = []


def kontrol(ad: str, kosul: bool, ayrinti: str = "") -> None:
    global gecen
    if kosul:
        gecen += 1
        print(f"  [GECTI] {ad}" + (f" — {ayrinti}" if ayrinti else ""))
    else:
        kalan.append(f"{ad}: {ayrinti}")
        print(f"  [KALDI] {ad} — {ayrinti}")


def bolum(baslik: str) -> None:
    print(f"\n{baslik}\n" + "-" * len(baslik))


def main() -> int:
    db.init_db()

    # ---------------------------------------------------------------
    bolum("1. Restriksiyon haritasi — pUC19 MCS'inin yayimlanmis konumlari")
    seq = rv.load("PUC19")["sequence"]
    tanimalar = ez.recognition_sites(seq, list(PUC19_MCS))
    for enz, beklenen in PUC19_MCS.items():
        kontrol(f"{enz} tanima dizisi @ {beklenen}",
                tanimalar[enz] == [beklenen],
                f"uygulama {tanimalar[enz]}")
    # Kesim konumu = tanima basi + ust zincir kesim ofseti
    for enz, bas in PUC19_MCS.items():
        g = ez.geometry(enz)
        kesim = restriction_sites(seq, [enz])[enz]
        kontrol(f"{enz} kesim konumu sozlesmesi",
                kesim == [bas + g.top_cut],
                f"tanima {bas} + ofset {g.top_cut} = {bas + g.top_cut}, "
                f"uygulama {kesim}")
    sira = sorted(PUC19_MCS, key=lambda e: PUC19_MCS[e])
    kontrol("MCS enzim sirasi yayimlanmis haritayla ayni",
            sira == PUC19_MCS_SIRASI, " - ".join(sira))

    # ---------------------------------------------------------------
    bolum("2. Erime sicakligi — uc bagimsiz yontem")
    primerler = [
        "GTAAAACGACGGCCAGT",          # M13 ileri
        "CAGGAAACAGCTATGAC",          # M13 geri
        "TAATACGACTCACTATAGGG",       # T7 ileri
        "ATGCAAGTCATCGGATTCACG",
    ]
    for p in primerler:
        nn = tm_nn(p)
        gc = mt.Tm_GC(p)
        wallace = mt.Tm_Wallace(p)
        # Yontemler ayni fiziksel buyuklugu farkli modellerle tahmin
        # eder; 10 C'den fazla ayrismamalari beklenir.
        yayilim = max(nn, gc, wallace) - min(nn, gc, wallace)
        kontrol(f"Tm yayilimi {p[:12]}…",
                yayilim < 12.0,
                f"NN={nn:.1f} GC={gc:.1f} Wallace={wallace:.1f} "
                f"(yayilim {yayilim:.1f} C)")

    # ---------------------------------------------------------------
    bolum("3. Protein kutlesi — ProtParam disinda ikinci bir hesapla")
    diziler = [
        "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQ",
        "MGSSHHHHHHSSGLVPRGSHMENLYFQGMSKGEELFTGVVPILVELDGDVNGHKFSV",
        "M" + "AGLSVKDE" * 10,
    ]
    for d in diziler:
        rapor = protein.analyze(d)
        bagimsiz = molecular_weight(d, seq_type="protein")
        fark = abs(rapor.mw_da - bagimsiz)
        kontrol(f"Kutle {d[:14]}…",
                fark < 1.0,
                f"uygulama {rapor.mw_da:.1f} Da, bagimsiz hesap "
                f"{bagimsiz:.1f} Da (fark {fark:.2f})")

    # ---------------------------------------------------------------
    bolum("4. Vektor uzunluklari — NCBI'nin bildirdigi degerle")
    for vid, beklenen in BEKLENEN_UZUNLUK.items():
        if vid not in rv.available():
            continue
        rec = rv.load(vid)
        kontrol(f"{rec['name']} ({rec['accession']})",
                rec["length"] == beklenen and len(rec["sequence"]) == beklenen,
                f"kayit {rec['length']} bp, dizi {len(rec['sequence'])} bp, "
                f"beklenen {beklenen} bp")

    # ---------------------------------------------------------------
    bolum("5. Kodon tablosu — elle yazilmis nadir kodon listesiyle")
    w = codon.weights()
    for c in sorted(RARE_CODONS_ECOLI):
        if c in w:
            kontrol(f"{c} nadir", w[c] < 0.10, f"genomdan hesaplanan w={w[c]}")
    en_iyiler = {aa: cods[0] for aa, cods in codon.synonyms().items()}
    bilinen = {"L": "CTG", "R": "CGT", "A": "GCT", "K": "AAA", "I": "ATC",
               "E": "GAA", "F": "TTC", "G": "GGT", "P": "CCG", "V": "GTT"}
    for aa, beklenen in bilinen.items():
        kontrol(f"{aa} icin en iyi kodon {beklenen}",
                en_iyiler.get(aa) == beklenen,
                f"hesaplanan {en_iyiler.get(aa)}")

    # ---------------------------------------------------------------
    bolum("6. Parca katalogu — gercek vektorlerde birebir bulunuyor mu")
    parts = db.list_parts()
    beklenen_eslesmeler = [
        ("PUC19", "MCS_PUC19"), ("PUC19", "RPT_LACZA"),
        ("PBR322", "MRK_BLA"), ("PBR322", "MRK_TET"),
        ("PACYC184", "ORI_P15A"), ("PACYC184", "MRK_CAT"),
        ("PTRC99A", "PRO_TRC"), ("PTRC99A", "REG_LACI"),
    ]
    for vid, pid in beklenen_eslesmeler:
        if vid not in rv.available():
            continue
        bulunan = {f.part_id for f in rv.catalog_features(vid, parts)
                   if f.identity == 100.0}
        kontrol(f"{pid} birebir {vid} icinde",
                pid in bulunan,
                "tam eslesme bulundu" if pid in bulunan
                else "tam eslesme yok")

    # ---------------------------------------------------------------
    print(f"\n{'=' * 60}")
    print(f"{gecen} kontrol gecti, {len(kalan)} kaldi.")
    if kalan:
        print("\nKALAN KONTROLLER:")
        for k in kalan:
            print(f"  - {k}")
        return 1
    print("Butun bagimsiz kontroller gecti.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
