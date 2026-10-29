"""Dizi motorunun temel işlemleri.

Bu testlerin amacı "çalışıyor mu" değil, "doğru mu". Her biri bağımsız
olarak doğrulanabilir bir gerçeğe dayanır: ters tümleyenin iki kez
alınması diziyi geri verir, pUC19'un EcoRI bölgesi 396. konumdadır,
REBASE'in bildirdiği çıkıntı gerçek diziden çıkarılanla aynıdır.
"""

from __future__ import annotations

import pytest

from core import enzymes as ez
from core.seq_utils import (
    clean_seq,
    find_orfs,
    gc_percent,
    primer_dimer,
    primer_hairpin,
    restriction_sites,
    revcomp,
    tm_nn,
    translate_cds,
)

GEN = ("ATGAAACTGAAATTTGCAGAAGTTCTGAAACCGGGTTCTAGCGGTAAAATTGATCTG"
       "GAAGTTCTGAAACCGGGTTCTAGCGGTAAAATTTAA")


def test_revcomp_tersi_kendisidir():
    for s in ("ATGC", GEN, "AAAA", "GCGCGC"):
        assert revcomp(revcomp(s)) == s


def test_revcomp_bilinen_deger():
    assert revcomp("GAATTC") == "GAATTC"      # EcoRI palindromik
    assert revcomp("ATGC") == "GCAT"


def test_gc_hesabi():
    assert gc_percent("GGCC") == 100.0
    assert gc_percent("ATAT") == 0.0
    assert gc_percent("ATGC") == 50.0


def test_clean_seq_fasta_ve_bosluk_temizler():
    assert clean_seq(">başlık\nATGC\nGGTT\n") == "ATGCGGTT"
    assert clean_seq("atg c\tg") == "ATGCG"


def test_tm_uzun_dizide_daha_yuksek():
    kisa = tm_nn("ATGCATGCATGCATGCAT")
    uzun = tm_nn("ATGCATGCATGCATGCATGCATGCATGC")
    assert kisa is not None and uzun is not None
    assert uzun > kisa


def test_gc_zengin_primer_daha_yuksek_tm():
    at = tm_nn("ATATATATATATATATATAT")
    gc = tm_nn("GCGCGCGCGCGCGCGCGCGC")
    assert gc > at


# ---------------------------------------------------------------------------
# ORF ve çeviri
# ---------------------------------------------------------------------------


def test_orf_bulunur_ve_koordinatlari_tutar():
    orfs = find_orfs(GEN, min_aa=20, circular=False)
    assert orfs, "açık okuma çerçevesi bulunamadı"
    en_uzun = max(orfs, key=lambda o: o.length_aa)
    assert en_uzun.start == 1
    assert GEN[en_uzun.start - 1:en_uzun.start + 2] == "ATG"


def test_orf_ters_zincirde_de_bulunur():
    orfs = find_orfs(revcomp(GEN), min_aa=20, circular=False)
    assert any(o.strand == -1 for o in orfs)


def test_ceviri_durdurma_kodonunu_isaretler():
    """Durdurma kodonu `*` olarak korunur - çeviriyi okuyan kişi
    çerçevenin gerçekten kapandığını görmelidir."""
    r = translate_cds(GEN)
    assert r["protein"].startswith("MKLKFAEVLKPGSSGKID")
    assert r["protein"].endswith("*")
    assert r["protein"].count("*") == 1


# ---------------------------------------------------------------------------
# Restriksiyon
# ---------------------------------------------------------------------------


def test_ecori_bolgesi_bulunur():
    s = "AAAA" + "GAATTC" + "TTTT"
    sites = restriction_sites(s, ["EcoRI"], circular=False)
    assert len(sites["EcoRI"]) == 1


def test_bulunmayan_enzim_bos_liste_verir():
    sites = restriction_sites("AAAAAAAAAAAA", ["EcoRI"], circular=False)
    assert sites["EcoRI"] == []


@pytest.mark.parametrize("name", [
    "EcoRI", "BamHI", "HindIII", "PstI", "KpnI", "SacI",
    "SmaI", "EcoRV", "NotI", "XhoI", "NdeI", "SalI",
])
def test_cikinti_rebase_ile_tutarli(name):
    """Gerçek diziden çıkarılan çıkıntı, REBASE'in bildirdiğiyle aynı mı.

    Bu, kesim geometrisinin doğruluğunu bağımsız olarak sınar: biz
    koordinatlardan çıkıntıyı diziden kesip alıyoruz, REBASE ise
    çıkıntıyı kendi alanında bildiriyor. İkisi tutmazsa koordinat
    hesabı yanlıştır.
    """
    d = ez.details(name)
    assert d, f"{name} REBASE'de bulunamadı"
    bolge = d["site"]
    seq = "GGGGGGGGGG" + bolge + "CCCCCCCCCC"
    konumlar = ez.sites_in(seq, [name], circular=False)[name]
    assert konumlar, f"{name} kendi tanıma dizisini kesmedi"
    ctx = ez.cut_context(seq, name, konumlar[0], flank=10)
    assert ctx is not None
    beklenen = (d.get("overhang_seq") or "").upper()
    if beklenen and "N" not in beklenen:
        assert ctx["overhang_seq"].upper() == beklenen, (
            f"{name}: dizi {ctx['overhang_seq']}, REBASE {beklenen}")


def test_kor_uc_enzimleri_cikinti_birakmaz():
    for name in ("SmaI", "EcoRV", "PvuII"):
        d = ez.details(name)
        assert d["blunt"] is True
        assert not d["overhang_seq"]


def test_tip_iis_tanima_disindan_keser():
    for name in ("BsaI", "BsmBI", "BbsI"):
        g = ez.geometry(name)
        assert g is not None, f"{name} geometrisi okunamadı"
        assert g.outside, f"{name} Tip IIS olarak işaretlenmedi"


# ---------------------------------------------------------------------------
# Primer dimer
# ---------------------------------------------------------------------------


def test_tamamen_tumleyen_primerler_dimer_verir():
    a = "GCATCGATCGGGCCCAAATTT"
    d = primer_dimer(a, revcomp(a))
    assert d["risk"] == "yüksek"
    assert d["end_match"] >= 3


def test_temiz_primer_cifti_dimer_vermez():
    d = primer_dimer("ATGCAAGTCATCGGATTCACG", "TTCAGGAAGTTGCAATCACCT")
    assert d["risk"] == "yok"


def test_palindromik_primer_sac_tokasi_verir():
    h = primer_hairpin("GCGCGCATTTTGCGCGC")
    assert h["risk"] != "yok"


def test_duz_primer_sac_tokasi_vermez():
    h = primer_hairpin("ATGCAAGTCATCGGATTCACG")
    assert h["risk"] == "yok"
