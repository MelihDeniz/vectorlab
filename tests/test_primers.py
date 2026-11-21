"""Vektöre özel primer tasarımı.

Buradaki testlerin çoğu bir *değişmezi* sınar: koloni PCR primerleri
insertin dışında kalmalı, iki bandın farkı insert uzunluğuna eşit
olmalı, 3' ucu oturmayan primer bağlanmış sayılmamalı.
"""

from __future__ import annotations

import pytest

from core import db, primers as pm
from core.assembly import assemble, sort_canonical
from core.seq_utils import clean_seq, revcomp

PARCALAR = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"]
INSERT = "ATG" + "AAACTGAAAGTTCTGGAAGCGACCGTT" * 12 + "TAA"


@pytest.fixture(scope="module")
def vektor():
    db.init_db()
    return assemble("pTestPrimer", sort_canonical(db.get_parts(PARCALAR)),
                    insert_seq=INSERT, insert_name="hedef")


@pytest.fixture(scope="module")
def insert_araligi(vektor):
    f = [x for x in vektor.features if x.part_id == "INSERT"]
    assert f, "insert özniteliği bulunamadı"
    return f[0].start, f[0].end


# ---------------------------------------------------------------------------
# Özgüllük
# ---------------------------------------------------------------------------


def test_vektorden_kesilen_dizi_kendi_yerinde_bulunuyor(vektor):
    parca = vektor.sequence[200:224]
    yerler = pm.binding_sites(parca, vektor.sequence)
    tam = [y for y in yerler if y.mismatches == 0 and y.strand == 1]
    assert tam, "diziden kesilen parça kendi kalıbında bulunamadı"
    assert any(y.position == 201 for y in tam), (
        f"beklenen konum 201, bulunan: {[y.position for y in tam]}")


def test_ters_tumleyen_alt_zincirde_bulunuyor(vektor):
    parca = revcomp(vektor.sequence[300:324])
    yerler = pm.binding_sites(parca, vektor.sequence)
    assert any(y.strand == -1 and y.mismatches == 0 for y in yerler), (
        "ters tümleyen alt zincirde bulunamadı")


def test_uc_prime_uyumsuzlugu_baglanmayi_eliyor(vektor):
    """3' ucu oturmayan primer bağlanmış sayılmamalı.

    Gerekçe biyolojik: polimeraz uzamayı 3' ucundan başlatır. Bu
    ayrım olmasaydı tarama, laboratuvarda hiç ürün vermeyecek
    primerleri "bağlanıyor" diye gösterirdi.
    """
    parca = vektor.sequence[400:424]
    # 5' ucunda iki baz boz: bağlanma SAYILMALI.
    bes_bozuk = ("T" if parca[0] != "T" else "A") + \
                ("T" if parca[1] != "T" else "A") + parca[2:]
    assert pm.binding_sites(bes_bozuk, vektor.sequence), (
        "5' ucundaki uyumsuzluk bağlanmayı elememeli")

    # 3' ucunun SON bazını boz: bağlanma SAYILMAMALI.
    son = parca[-1]
    uc_bozuk = parca[:-1] + ("T" if son != "T" else "A")
    tam_yer = [y for y in pm.binding_sites(uc_bozuk, vektor.sequence)
               if y.position == 401 and y.strand == 1]
    assert not tam_yer, "3' ucu oturmayan primer bağlanmış sayıldı"


def test_tekrarlanan_dizi_birden_fazla_yer_veriyor():
    """Aynı dizi iki kez geçiyorsa özgüllük uyarı vermeli."""
    birim = "ATGCCGTTACGGATCAGTTC"
    kalip = "AAAA" + birim + "GGGGCCCCTTTTAAAACCCC" + birim + "TTTT"
    oz = pm.specificity(birim, kalip, circular=False)
    assert not oz["unique"]
    assert len(oz["sites"]) >= 2
    assert oz["notes"][0]["kind"] == "uyari"


def test_dairesel_kalipta_baslangictan_tasan_primer_bulunuyor():
    """Dairesel plazmitte primer 1. konumun üstünden geçebilir."""
    kalip = "ATGCCGTTACGGATCAGTTCAAAGGGCCCTTTACGTACGTACGT"
    # Sonun son 10 bazı + başın ilk 10 bazı.
    tasan = kalip[-10:] + kalip[:10]
    assert pm.binding_sites(tasan, kalip, circular=True), (
        "dairesel kalıpta sarılan primer bulunamadı")
    assert not pm.binding_sites(tasan, kalip, circular=False), (
        "doğrusal kalıpta sarılan primer bulunmamalıydı")


# ---------------------------------------------------------------------------
# Koloni PCR
# ---------------------------------------------------------------------------


def test_koloni_primerleri_insertin_disinda(vektor, insert_araligi):
    """Koloni PCR'ın bütün anlamı budur.

    Insertin içine bağlanan primer boş vektörde hiç ürün vermez, yani
    boş koloni ile tutmamış PCR ayırt edilemez. Primerler iskelette
    olmak ZORUNDA.
    """
    bas, son = insert_araligi
    r = pm.colony_pcr(vektor.sequence, bas, son)
    insert_dizi = vektor.sequence[bas - 1:son]
    for p in (r.forward, r.reverse):
        assert p.sequence not in insert_dizi, (
            f"{p.name} insertin içinden alınmış")
        assert revcomp(p.sequence) not in insert_dizi, (
            f"{p.name} (ters) insertin içinden alınmış")


def test_bantlar_arasindaki_fark_insert_uzunlugu(vektor, insert_araligi):
    """Dolu ve boş bant arasındaki fark tam olarak insert kadar olmalı."""
    bas, son = insert_araligi
    r = pm.colony_pcr(vektor.sequence, bas, son)
    assert r.fark == son - bas + 1, (
        f"fark {r.fark}, insert {son - bas + 1} bp")
    assert r.with_insert > r.without_insert > 0


def test_koloni_primerleri_vektorde_tek_yere_baglaniyor(vektor,
                                                        insert_araligi):
    bas, son = insert_araligi
    r = pm.colony_pcr(vektor.sequence, bas, son)
    for p in (r.forward, r.reverse):
        oz = pm.specificity(p.sequence, vektor.sequence)
        assert oz["unique"], (
            f"{p.name} vektörde {len(oz['sites'])} yere bağlanıyor")


def test_kisa_insert_ayirt_edilemez_diye_uyariyor(vektor):
    """İki bant jelde ayrılamayacaksa bunu söylemeli."""
    r = pm.colony_pcr(vektor.sequence, 100, 140)   # 41 bp insert
    assert r.fark < 100
    assert any("ayrılmayabilir" in n["title"] for n in r.notes), (
        "kısa insert için uyarı verilmedi")


def test_gecersiz_koordinat_hata_veriyor(vektor):
    r = pm.colony_pcr(vektor.sequence, 500, 100)
    assert r.notes and r.notes[0]["kind"] == "hata"


# ---------------------------------------------------------------------------
# Dizileme çapaları
# ---------------------------------------------------------------------------


def test_capalar_katalog_dizisinden_geliyor(vektor):
    """Önerilen primer uydurulmuş değil, vektörün kendi dizisinden.

    Evrensel dizileme primerlerini ezberden yazmak yerine katalogda
    doğrulanmış parçalardan kesmenin sebebi bu: her primerin kaynağı
    gösterilebilir olmalı.
    """
    capalar = pm.catalog_anchors(vektor.sequence, db.list_parts())
    assert capalar, "hiç çapa bulunamadı"
    for a in capalar:
        assert a["primer"] in vektor.sequence, (
            f"{a['part_id']} primeri vektörde birebir bulunmuyor")
        parca = db.get_parts([a["part_id"]])[0]
        assert a["primer"] in clean_seq(parca["sequence"]), (
            f"{a['part_id']} primeri katalog parçasının içinde değil")


def test_t7_capasi_bulunuyor(vektor):
    """T7 promotörü taşıyan vektörde T7 çapası çıkmalı."""
    capalar = pm.catalog_anchors(vektor.sequence, db.list_parts())
    assert any("T7" in a["part_id"] for a in capalar), (
        f"T7 çapası yok: {[a['part_id'] for a in capalar]}")


def test_dizisi_eksik_parca_capa_olmuyor(vektor):
    """N bloğu taşıyan parçadan primer çıkarılmamalı."""
    sahte = [{"part_id": "SAHTE", "name": "eksik",
              "category": "promoter", "sequence": "N" * 40}]
    assert pm.catalog_anchors(vektor.sequence, sahte) == []


# ---------------------------------------------------------------------------
# Primer inceleme
# ---------------------------------------------------------------------------


def test_incele_olculeri_dolduruyor(vektor):
    p = vektor.sequence[150:172]
    r = pm.analyze(p, vektor.sequence)
    assert "error" not in r
    assert r["metrics"]["length"] == len(p)
    assert r["metrics"]["tm"] is not None
    assert r["specificity"] is not None


def test_cok_kisa_primer_reddediliyor():
    assert "error" in pm.analyze("ATGC")


def test_kotu_primer_uyari_aliyor():
    """GC'si çok düşük, 3' ucu A olan primer uyarı almalı."""
    kotu = "ATATATATATATAAATATATATA"
    r = pm.analyze(kotu)
    basliklar = " ".join(n["title"] for n in r["notes"])
    assert "GC" in basliklar
    assert any(n["kind"] == "uyari" for n in r["notes"])


def test_kalipsiz_incelemede_ozgulluk_yok():
    r = pm.analyze("ATGCCGTTACGGATCAGTTCGG")
    assert r["specificity"] is None
