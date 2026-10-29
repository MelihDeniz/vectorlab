"""Kodon, mRNA yapısı, montaj, klonlama ve doğrulama modülleri.

Buradaki testlerin çoğu bir *değişmezi* sınar: optimizasyon proteini
değiştirmemeli, hizalama kendi kendine %100 benzerlik vermeli, kesim
ürünleri toplamı kalıbın uzunluğuna eşit olmalı. Değişmezler, beklenen
değer yazmaktan daha güçlüdür çünkü girdi değişse de geçerli kalırlar.
"""

from __future__ import annotations

import pytest

from core import (
    assembly_methods as am,
    cloning,
    codon,
    db,
    gel,
    protein,
    protocol as pr,
    ref_vectors as rv,
    structure as stx,
    verification as vf,
)
from core.assembly import assemble, sort_canonical
from core.seq_utils import revcomp
from core.validate import validate

GEN = ("ATG" + "AAACTGAAATTTGCAGAAGTTCTGAAACCGGGTTCTAGCGGTAAAATT" * 4
       + "TAA")
RBS = "AAAGAGGAGAAATACTAG"


@pytest.fixture(scope="module", autouse=True)
def _veritabani():
    db.init_db()


# ---------------------------------------------------------------------------
# Kodon
# ---------------------------------------------------------------------------


def test_kodon_tablosu_yuklendi():
    t = codon.table()
    assert t, "kodon tablosu yok - tools/fetch_codon_usage.py çalıştırın"
    assert t["gene_count"] > 3000
    assert len(t["codons"]) == 64


def test_her_amino_asidin_bir_en_iyi_kodonu_var():
    w = codon.weights()
    syn = codon.synonyms()
    for aa, codons in syn.items():
        assert max(w[c] for c in codons) == pytest.approx(1.0)


def test_bilinen_nadir_kodonlar_dusuk_w_alir():
    """Elle yazılmış nadir kodon listesi ile genomdan hesaplanan tablo
    bağımsız olarak aynı sonucu vermeli."""
    from core.seq_utils import RARE_CODONS_ECOLI
    w = codon.weights()
    for c in RARE_CODONS_ECOLI:
        if c in w:
            assert w[c] < 0.10, f"{c} nadir sayılıyor ama w={w[c]}"


def test_optimizasyon_proteini_degistirmez():
    for kip in ("nadir", "tam"):
        r = codon.optimize(GEN, mode=kip)
        assert codon._translate(r.sequence) == codon._translate(GEN)


def test_optimizasyon_cai_yi_1_e_cikarmaz():
    """1,0 doğada karşılığı olmayan bir hedeftir; gerçek genler 0,8
    civarındadır. Optimizasyon oraya gitmemeli."""
    r = codon.optimize(GEN, mode="tam")
    assert r.cai_after < 0.97, (
        f"CAI {r.cai_after} - tek kodona yüklenmiş olabilir")
    assert r.cai_after > r.cai_before


def test_optimizasyon_nadir_kodonlari_azaltir():
    kotu = "ATG" + ("AGG" "CTA" "ATA" "CCC" "CGA") * 12 + "TAA"
    r = codon.optimize(kotu, mode="nadir")
    assert r.rare_after < r.rare_before
    assert r.rare_after == 0


def test_optimizasyon_yasakli_bolgeyi_kaldirir():
    hedef = "ATG" + "GAATTCAAACTGAAAGGATCCGGTTCT" * 5 + "TAA"
    r = codon.optimize(hedef, mode="tam",
                       forbidden_sites=["GAATTC", "GGATCC"])
    from core.seq_utils import restriction_sites
    kalan = restriction_sites(r.sequence, ["EcoRI", "BamHI"], circular=False)
    assert sum(len(v) for v in kalan.values()) == 0
    assert codon._translate(r.sequence) == codon._translate(hedef)


def test_kaldirilamayan_bolge_durustce_raporlanir():
    """NdeI (CATATG) başlama kodonunun üstündeyse kaldırılamaz, çünkü
    Met'in tek kodonu vardır. Sessizce geçmemeli."""
    m = "ATG" + "CATATGAAACTGAAATAA"
    r = codon.optimize(m, mode="tam", forbidden_sites=["CATATG"])
    assert r.sites_remaining, "kaldırılamayan bölge bildirilmedi"


def test_tekrar_tespiti_calisir():
    tekrarli = "ATG" + "CTGAAACTGAAACTGAAAGAAGTTCTGAAA" * 4 + "TAA"
    tekrarlar = codon.find_repeats(tekrarli, min_len=20)
    assert tekrarlar and tekrarlar[0]["count"] >= 2


def test_gercek_e_coli_genleri_yuksek_cai_alir():
    """rpsA gibi çok ifade edilen genler yüksek, yabancı benzeri
    diziler düşük CAI almalı."""
    iyi = "ATG" + "CTGGAAAAAGCGACCGTTATCGAAGGTCTGGAA" * 4 + "TAA"
    kotu = "ATG" + "CTAGAGAAGGCGACAGTTATAGAGGGTCTAGAG" * 4 + "TAA"
    assert codon.cai(iyi) > codon.cai(kotu)


# ---------------------------------------------------------------------------
# mRNA yapısı
# ---------------------------------------------------------------------------


def test_istif_enerjisi_simetriktir():
    """TC/AG tabloda yok ama aynı fiziksel istif olan GA/CT var."""
    assert stx._stack_dg("TC", "AG") == pytest.approx(
        stx._stack_dg("GA", "CT"))
    assert stx._stack_dg("CC", "GG") == pytest.approx(
        stx._stack_dg("GG", "CC"))


def test_gc_istifi_at_istifinden_kararli():
    assert stx._stack_dg("GC", "CG") < stx._stack_dg("AA", "TT")


def test_baslama_kodonu_promotor_icindeki_atg_yi_secmez():
    """Promotörün -10 kutusundan sonra gelen ATG tuzaktır."""
    dizi = ("TTTACAGCTAGCTCAGTCCTAGGTATAATGCTAGC" + RBS
            + "ATGAAACTGAAATTTGCAGAAGTTCTGAAA")
    pos = stx.guess_start(dizi)
    assert pos == dizi.index(RBS) + len(RBS) + 1


def test_sac_tokasi_baslama_kodonunu_kapatiyorsa_isaretlenir():
    sap = "GCCGGCGGCCGCCGGCA"
    dizi = RBS + "ATG" + sap + "TTTT" + revcomp(sap) + "AAACTGAAATAA"
    rap = stx.analyze_start_region(dizi)
    assert rap.blocking_start is not None
    # Şiddet yapının kararlılığına bağlı (-8 kcal/mol altı hata); burada
    # aranan, bölgenin işaretlenmiş olması.
    assert any(n["kind"] in ("hata", "uyari")
               and "başlama kodonunu kapatıyor" in n["title"]
               for n in rap.notes)


def test_acik_baslangic_bolgesinde_uyari_yok():
    dizi = RBS + "ATGAAACTGAAATTTGCAGAAGTTCTGAAACCGGGTTCTAGC"
    rap = stx.analyze_start_region(dizi)
    assert rap.blocking_start is None
    assert rap.blocking_rbs is None


# ---------------------------------------------------------------------------
# Klonlama ve montaj
# ---------------------------------------------------------------------------


def test_ligasyon_molar_orani():
    """ng insert = ng vektör × (insert bp / vektör bp) × oran"""
    r = pr.ligation_amounts(5000, 1000, 50, 3.0)
    assert r["insert_ng"] == pytest.approx(30.0)
    r2 = pr.ligation_amounts(3000, 3000, 100, 1.0)
    assert r2["insert_ng"] == pytest.approx(100.0)


def test_fragman_uzunluklari_toplami_kalibi_verir():
    seq = rv.load("PUC19")["sequence"]
    for enzimler in (["EcoRI", "BamHI"], ["HaeIII"], ["PstI"]):
        frag = cloning.fragments(seq, enzimler)
        assert sum(f["length"] for f in frag) == len(seq)


def test_tek_enzimle_kesim_yonsuzluk_uyarisi_verir():
    seq = rv.load("PUC19")["sequence"]
    p = cloning.plan(seq, "pUC19", ["EcoRI"], "ATG" * 60)
    assert any("yönsüz" in n.title for n in p.notes)


def test_insert_icinde_kesen_enzim_hata_verir():
    seq = rv.load("PUC19")["sequence"]
    p = cloning.plan(seq, "pUC19", ["EcoRI", "BamHI"],
                     "ATGGAATTCAAA" * 10)
    assert not p.ok


def test_golden_gate_ayni_cikintiyi_yakalar():
    parcalar = [("A", "ATGGCCGGTTTAAACCCGGGAAATTTCCCGGG"),
                ("B", "ATGGTTTTAACCCGGGAAATTTCCCGGGAAAT")]
    p = am.golden_gate(parcalar, "BsaI")
    assert not p.ok


def test_golden_gate_ic_bolgeyi_yakalar():
    parcalar = [("X", "ATGGGTCTCAAACCCGGGAAATTTCCCGGGAAA"),
                ("Y", "TTGACCGGTTTAAACCCGGG")]
    p = am.golden_gate(parcalar, "BsaI")
    assert any("içinde kesiyor" in n["title"] for n in p.notes)


def test_gibson_kisa_parcayi_reddeder():
    p = am.gibson([("kisa", "ATGCATGC"), ("uzun", "ATGC" * 30)], overlap=20)
    assert not p.ok


def test_mutagenez_degisikligi_ortada_tasir():
    p = am.mutagenesis("ATGAAACTGAAATTTGCAGAAGTTCTGAAACCGGGTTAA", 10, "GAT")
    assert len(p.primers) == 2
    f, r = p.primers
    assert r.sequence == revcomp(f.sequence)
    assert "GAT" in f.sequence


# ---------------------------------------------------------------------------
# Jel
# ---------------------------------------------------------------------------


def test_jel_goc_sirasi_dogru():
    """Büyük fragman yukarıda, küçük aşağıda olmalı."""
    y_buyuk = gel._y_for(5000, 1.0, 0.0, 100.0)
    y_orta = gel._y_for(1000, 1.0, 0.0, 100.0)
    y_kucuk = gel._y_for(500, 1.0, 0.0, 100.0)
    assert y_buyuk < y_orta < y_kucuk


def test_yakin_boyutlar_tek_bant_olur():
    ln = gel.build_lane("t", [2000, 2050, 500], 1.0)
    birlesen = [b for b in ln.bands if b.merged_from]
    assert birlesen, "yakın boyutlar birleştirilmedi"


def test_kesilmemis_plazmit_uc_bant_verir():
    ln = gel.uncut_lane("kesilmemiş", 3000, 1.0)
    assert len(ln.bands) == 3
    etiketler = " ".join(b.label for b in ln.bands)
    assert "süperkıvrımlı" in etiketler


def test_superkivrimli_daha_hizli_gocer():
    ln = gel.uncut_lane("k", 3000, 1.0)
    boyutlar = {b.label.split()[0]: b.size for b in ln.bands}
    assert boyutlar["süperkıvrımlı"] < boyutlar["doğrusal"]
    assert boyutlar["doğrusal"] < boyutlar["açık"]


def test_agaroz_onerisi_kucuk_fragmanda_yuksek_yuzde():
    assert gel.recommend_agarose([150, 200, 300]) > \
        gel.recommend_agarose([3000, 5000, 8000])


# ---------------------------------------------------------------------------
# Doğrulama
# ---------------------------------------------------------------------------


def test_kendi_kendine_hizalama_tam_benzerlik():
    r = vf.compare_read(GEN, GEN)
    assert r.identity == pytest.approx(100.0)
    assert not r.differences


def test_ters_zincir_okumasi_cevrilir():
    r = vf.compare_read(GEN, revcomp(GEN))
    assert r.reversed_read
    assert r.identity == pytest.approx(100.0)


def test_tek_baz_degisimi_bulunur():
    bozuk = list(GEN)
    bozuk[30] = "T" if bozuk[30] != "T" else "A"
    r = vf.compare_read(GEN, "".join(bozuk))
    assert len(r.differences) == 1
    assert r.differences[0].position == 31
    assert isinstance(r.differences[0].position, int)


def test_silme_cerceve_kaymasi_olarak_isaretlenir():
    silik = GEN[:40] + GEN[41:]
    r = vf.compare_read(GEN, silik)
    assert any(f.kind != "degisim" for f in r.differences)


def test_erken_durdurma_hata_verir():
    dur = list(GEN)
    dur[3:6] = list("TAA")
    r = vf.compare_read(GEN, "".join(dur), cds_start=1)
    assert any(n["kind"] == "hata" for n in r.notes)


def test_dizileme_plani_kapsama_hesaplar():
    seq = rv.load("PUC19")["sequence"]
    plan = vf.sequencing_plan(seq)
    assert plan.primers
    assert plan.covered > plan.length * 0.9
    for p in plan.primers:
        assert seq.count(p.sequence) >= 1


# ---------------------------------------------------------------------------
# Protein
# ---------------------------------------------------------------------------


def test_protein_kutlesi_uzunlukla_orantili():
    kisa = protein.analyze("M" + "A" * 50)
    uzun = protein.analyze("M" + "A" * 100)
    assert uzun.mw_da > kisa.mw_da


def test_trp_tyr_yoksa_a280_olculemez():
    r = protein.analyze("M" + "AGLSVKDEAGLSVKDE" * 6)
    assert r.extinction_reduced == 0
    assert r.a280_for_1mg_ml is None
    assert any("A280" in n["title"] for n in r.notes)


def test_his_etiketi_bulunur():
    r = protein.analyze("MGSSHHHHHHSSGLVPRGSH" + "AGLSVKDE" * 8)
    assert any(t[0] == "His6" for t in r.tags)


def test_cok_kisa_dizi_none_doner():
    assert protein.analyze("MKV") is None


# ---------------------------------------------------------------------------
# Kromatogram kalite kırpması
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ad,kalite,beklenen", [
    ("iki ucu kötü", [5] * 40 + [38] * 500 + [6] * 60, (40, 540)),
    ("tamamı iyi", [40] * 300, (0, 300)),
    ("tamamı kötü", [3] * 300, (0, 0)),
    ("ortada tek kötü baz", [38] * 100 + [2] + [38] * 100, (0, 201)),
    ("başı kötü", [4] * 80 + [36] * 400, (80, 480)),
])
def test_kalite_kirpmasi(ad, kalite, beklenen):
    """Okumanın başı ve sonu neredeyse her zaman güvenilmezdir;
    kırpma o bölgeleri atmalı ama ortadaki tek bir kötü bazda
    okumayı ikiye bölmemeli."""
    assert vf._mott_trim(kalite, 20) == beklenen, ad


def test_bozuk_ab1_dosyasi_none_doner():
    assert vf.read_chromatogram(b"bu bir ab1 dosyasi degil") is None


# ---------------------------------------------------------------------------
# Çıkıntı kümesi
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [2, 4, 8, 12, 20])
def test_onerilen_cikintilar_kurallara_uyar(n):
    """Her çift arasında en az iki baz fark olmalı, hiçbiri palindromik
    ya da bir diğerinin tersi olmamalı."""
    kume = am.suggest_overhangs(n)
    assert len(kume) == n
    for i, a in enumerate(kume):
        assert a != revcomp(a), f"{a} palindromik"
        for b in kume[i + 1:]:
            assert a != b
            assert revcomp(a) != b
            assert am._hamming(a, b) >= 2, f"{a} ile {b} tek baz farklı"


def test_bozuk_cikinti_kumesi_yakalanir():
    sorunlar = am.check_overhang_set(["AATG", "AATG", "GGTA"])
    assert any(n["kind"] == "hata" for n in sorunlar)


def test_ters_cikinti_cifti_yakalanir():
    a = "AACA"
    sorunlar = am.check_overhang_set([a, revcomp(a)])
    assert any("tersi" in n["title"] for n in sorunlar)


def test_temiz_kume_sorun_vermez():
    sorunlar = am.check_overhang_set(am.suggest_overhangs(6))
    assert all(n["kind"] == "bilgi" for n in sorunlar)


def test_hizalama_gorunumu_farki_isaretler():
    bozuk = list(GEN)
    bozuk[30] = "T" if bozuk[30] != "T" else "A"
    gorunum = vf.alignment_view(GEN, "".join(bozuk))
    assert "*" in gorunum, "fark işaretlenmedi"
    assert "benzerlik" in gorunum


def test_hizalama_gorunumu_tam_eslesmede_yildiz_yok():
    gorunum = vf.alignment_view(GEN, GEN)
    govde = "\n".join(l for l in gorunum.splitlines()
                      if l.startswith(" ") or l[:7].strip().isdigit())
    assert "*" not in govde


# ---------------------------------------------------------------------------
# Yüklenen plazmiti denetleme
# ---------------------------------------------------------------------------


def test_genbank_yuklemesi_okunur():
    from pathlib import Path

    from core import inspect_plasmid as ip
    gb = Path("data/reference_vectors/_genbank/J01749.gb")
    if not gb.exists():
        pytest.skip("ham GenBank önbelleği yok")
    d = ip.parse_upload(gb.read_bytes(), "pBR322.gb")
    assert d["source"] == "genbank"
    assert len(d["sequence"]) == 4361
    assert d["circular"]
    assert len(d["features"]) > 20


def test_fasta_yuklemesi_okunur():
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PUC19").encode(), "pUC19.fasta")
    assert d["source"] == "fasta"
    assert len(d["sequence"]) == 2686


def test_puc19_incelemesi_dogru_elemanlari_bulur():
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PUC19").encode(), "pUC19.fasta")
    ins = ip.inspect(d, db.list_parts())
    kategoriler = {p["category"] for p in ins.derived_parts}
    assert {"ori", "marker", "promoter"} <= kategoriler
    assert ins.coverage > 50


def test_inceleme_ayni_bolgeyi_iki_kez_saymaz():
    """Katalogda `bla` birden çok kayıt; üçü de aynı yere oturur.
    Daraltılmazsa 'birden fazla marker/orijin' yanlış hatası çıkar."""
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PBR322").encode(), "pBR322.fasta")
    ins = ip.inspect(d, db.list_parts())
    oris = [p for p in ins.derived_parts if p["category"] == "ori"]
    assert len(oris) == 1, [p["name"] for p in oris]


def test_yuklenen_plazmitte_eksiklik_hataya_donmez():
    """Anotasyon eksikliği, elemanın yokluğu demek değildir."""
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PSEVA231").encode(), "x.fasta")
    ins = ip.inspect(d, db.list_parts())
    vec = ip.to_vector(ins)
    kati = validate(vec, host=db.get_host("DH5A"), purpose="klonlama")
    yumusak = validate(vec, host=db.get_host("DH5A"), purpose="klonlama",
                       annotation_complete=False)
    assert any(x.severity == "hata" for x in kati)
    assert not any(x.rule_id.startswith("CORE") and x.severity == "hata"
                   for x in yumusak)


def test_f1_orijini_coklu_orijin_hatasi_vermez():
    """f1/M13 tek zincir orijinidir, ColE1 ile çatışmaz."""
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PBAD24").encode(), "pBAD24.fasta")
    ins = ip.inspect(d, db.list_parts())
    f = validate(ip.to_vector(ins), host=db.get_host("TOP10"),
                 purpose="ifade", annotation_complete=False)
    assert not any(x.rule_id == "ORI-01" for x in f)


def test_inceleme_dizisi_yuklenen_dizinin_kendisi():
    """Yapının dizisi parçaların birleşimi değil, yüklenen dizi olmalı."""
    from core import inspect_plasmid as ip
    d = ip.parse_upload(rv.to_fasta("PUC19").encode(), "x.fasta")
    ins = ip.inspect(d, db.list_parts())
    assert ip.to_vector(ins).sequence == rv.load("PUC19")["sequence"]


# ---------------------------------------------------------------------------
# Kombinatoryal kütüphane
# ---------------------------------------------------------------------------


def test_kutuphane_butun_birlesimleri_uretir():
    from core import library as lib
    v, _ = lib.build(
        fixed_part_ids=["ORI_P15A", "MRK_CAT", "TER_B0015"],
        variable={"promoter": ["IGEM_J23119", "IGEM_J23105"],
                  "rbs": ["RBS_B0034", "RBS_B0031"]},
        host=db.get_host("DH5A"), purpose="ifade")
    assert len(v) == 4
    assert len({x.name for x in v}) == 4


def test_kutuphane_ifadeye_gore_sirali():
    from core import library as lib
    v, _ = lib.build(
        fixed_part_ids=["ORI_P15A", "MRK_CAT"],
        variable={"promoter": ["IGEM_J23119", "IGEM_J23113"],
                  "rbs": ["RBS_B0034", "RBS_B0031"]},
        host=db.get_host("DH5A"), purpose="ifade")
    skorlar = [x.expression for x in v if x.expression is not None]
    assert skorlar == sorted(skorlar, reverse=True)


def test_anderson_gucu_konsensusten_hesaplanir():
    from core import library as lib
    guclu = db.get_part("IGEM_J23119")
    zayif = db.get_part("IGEM_J23113")
    s1, k1 = lib.strength_score(guclu)
    s2, k2 = lib.strength_score(zayif)
    assert s1 > s2
    assert "konsensüs" in k1 and "konsensüs" in k2


def test_katalog_etiketi_yedek_olarak_kullanilir():
    from core import library as lib
    skor, kaynak = lib.strength_score(db.get_part("RBS_B0031"))
    assert skor is not None
    assert "katalog" in kaynak


def test_kutuphane_sinir_asilinca_uyarir():
    from core import library as lib
    promotorler = [p["part_id"] for p in db.list_parts("promoter")][:20]
    rbsler = [p["part_id"] for p in db.list_parts("rbs")][:10]
    v, uyari = lib.build(
        fixed_part_ids=["ORI_P15A", "MRK_CAT"],
        variable={"promoter": promotorler, "rbs": rbsler},
        host=db.get_host("DH5A"), max_variants=25)
    assert len(v) == 25
    assert uyari


def test_kutuphane_araligi_hesaplanir():
    from core import library as lib
    v, _ = lib.build(
        fixed_part_ids=["ORI_P15A", "MRK_CAT"],
        variable={"promoter": ["IGEM_J23119", "IGEM_J23113"],
                  "rbs": ["RBS_B0034", "RBS_B0031"]},
        host=db.get_host("DH5A"))
    a = lib.span(v)
    assert a["fold"] and a["fold"] > 1
