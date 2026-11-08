"""core/refbuild.py: referans vektörü kes-yapıştır hesabı (Biopython gerekmez)."""

import json
from pathlib import Path

import pytest

from core import refbuild as rb

DATA = Path(__file__).resolve().parent.parent / "data" / "reference_vectors"


@pytest.fixture(scope="module")
def pacyc():
    return json.loads((DATA / "PACYC184.json").read_text())["sequence"].upper()


def test_tek_kesen_enzimler(pacyc):
    assert rb.find_sites(pacyc, "HindIII") == [1524]
    assert rb.find_sites(pacyc, "SalI") == [2146]
    assert rb.find_sites(pacyc, "BclI") == [3542]
    assert rb.find_sites(pacyc, "EcoRI") == [1]


def test_iki_asama_boyutlari(pacyc):
    r = rb.build(pacyc, "N" * 3202, "N" * 804)
    s1, s2 = r.steps
    assert s1.length == 6831 and s1.removed == 616
    assert s2.length == 6937 and s2.removed == 698
    assert r.final_length == 6937


def test_enzim_yerleri_korunur(pacyc):
    r = rb.build(pacyc, "N" * 3202, "N" * 804)
    for e in ("HindIII", "SalI"):
        assert len(rb.find_sites(r.steps[0].sequence, e)) == 1
    for e in ("BclI", "EcoRI"):
        assert len(rb.find_sites(r.final, e)) == 1


def test_koordinat_esleme(pacyc):
    s1 = rb.build(pacyc, "N" * 3202, "N" * 804).steps[0]
    assert rb.map_position(s1, 1493) == 1493          # ori sonu aynı
    assert rb.map_position(s1, 1530) is None          # atılan bölge
    assert rb.map_position(s1, 2771) == 5357          # tet sonu kaydı
    assert rb.map_position(s1, 4245) == 6831


def test_sindirim_boyutlari(pacyc):
    assert rb.digest_sizes(pacyc, ["HindIII", "SalI"]) == [3623, 622]
    r = rb.build(pacyc, "N" * 3202, "N" * 804)
    assert sum(rb.digest_sizes(r.final, ["BclI", "EcoRI"])) == 6937


def test_tek_kesen_olmayan_hata(pacyc):
    with pytest.raises(ValueError):
        rb.single_cutter(pacyc + "AAGCTT", "HindIII")


# ---------------------------------------------------------------------------
# map_feature: bölgenin bütün korunup korunmadığı
# ---------------------------------------------------------------------------


def _al(seq: str, a: int, b: int) -> str:
    """1 tabanlı, kapsayıcı; a > b ise dizinin başını sarar."""
    return seq[a - 1:b] if a <= b else seq[a - 1:] + seq[:b]


def test_map_feature_butun_korunan(pacyc):
    s1 = rb.build(pacyc, "N" * 3202, "N" * 804).steps[0]
    # HindIII (1524) öncesindeki p15A ori bölgesi olduğu gibi kalır.
    assert rb.map_feature(s1, 581, 1493, len(pacyc)) == (581, 1493, "tam")


def test_map_feature_tamamen_atilan(pacyc):
    s1 = rb.build(pacyc, "N" * 3202, "N" * 804).steps[0]
    # HindIII ile SalI arasında kalan bölge atılır.
    assert rb.map_feature(s1, 1540, 2100, len(pacyc)) is None


def test_map_feature_kesilen_bolge(pacyc):
    """Kesim bölgesini İÇİNE alan gen bozulur, gerilmez.

    map_position uçları ayrı taşıdığı için böyle bir geni insertin
    üstüne binmiş halde sapasağlam gösterirdi; bu testin varlık sebebi
    o hatayı bir daha yapmamak.
    """
    s1 = rb.build(pacyc, "N" * 3202, "N" * 804).steps[0]
    L = len(pacyc)
    tetA = (1581, 2771)                       # HindIII/SalI bölgesini kapsar
    a, b, durum = rb.map_feature(s1, *tetA, L)
    assert durum == "kesildi"
    kalan = b - a + 1
    assert kalan < tetA[1] - tetA[0] + 1
    # Kalan parça gerçekten orijinal genin içinden bir parça olmalı.
    assert _al(s1.sequence, a, b) in _al(pacyc, *tetA)


def test_map_feature_saran_bolge(pacyc):
    """Dizinin başını saran gen (pACYC184'te cat) doğru taşınır."""
    s1 = rb.build(pacyc, "N" * 3202, "N" * 804).steps[0]
    L = len(pacyc)
    a, b, durum = rb.map_feature(s1, 3805, 219, L)
    assert durum == "tam"
    assert _al(s1.sequence, a, b) == _al(pacyc, 3805, 219)


# ---------------------------------------------------------------------------
# Tarifler
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tarif_sonuclari():
    """Her tarifi kur; sayfa açılmadan önce kurulabildiklerini sına."""
    from core import db, ref_vectors
    from core.assembly import assemble

    out = {}
    for tarif in rb.RECIPES:
        rec = ref_vectors.load(tarif.vector_id)
        assert rec, f"{tarif.vector_id} kaydı yok"
        kasetler = [assemble(f"k{i}", db.get_parts(s.insert_parts))
                    for i, s in enumerate(tarif.steps)]
        res = rb.build_recipe(
            rec["sequence"],
            [(s.left, s.right, k.sequence, s.title)
             for s, k in zip(tarif.steps, kasetler)])
        out[tarif.key] = (tarif, rec, kasetler, res)
    return out


def test_tarif_anahtarlari_tekil():
    anahtarlar = [r.key for r in rb.RECIPES]
    assert len(anahtarlar) == len(set(anahtarlar))
    assert all(rb.recipe_by_key(k) for k in anahtarlar)
    assert rb.recipe_by_key("olmayan") is None


def test_her_tarif_kuruluyor(tarif_sonuclari):
    """Ürün, eklenen kaseti gerçekten taşımalı.

    Ölçüt "ürün orijinalden büyük" DEĞİL: pNZ8048 kaydı boş vektör
    değildir, içinde kargo vardır ve o kargoyu atıp yerine daha küçük
    bir kaset koymak ürünü küçültür. Doğru değişmez, kasetin ürünün
    içinde bulunmasıdır.
    """
    for key, (tarif, rec, kasetler, res) in tarif_sonuclari.items():
        son = kasetler[-1].sequence
        assert son in res.final + res.final[:len(son)], key
        assert res.final_length > sum(len(k.sequence) for k in kasetler), key


def test_tarifler_gercekten_kurulabilir(tarif_sonuclari):
    """Kullanılan enzim, EKLENECEK kasetin içinde kesmemeli.

    Bu yalnızca bir incelik değil: inserti o iki enzimle hazırlarsınız.
    Kasetin içinde aynı enzimin yeri varsa insert ortadan bölünür ve
    tarif kağıt üzerinde kalır. Üç tarif bu sınavı ilk yazıldığında
    geçemedi (pACYC184/PKU'da HindIII, pSEVA231'de EcoRI, pNZ8048'de
    hem EcoRI hem HindIII kasetin içindeydi) ve enzim çiftleri
    değiştirildi.
    """
    for key, (tarif, _rec, kasetler, _res) in tarif_sonuclari.items():
        for adim, kaset in zip(tarif.steps, kasetler):
            for e in (adim.left, adim.right):
                n = len(rb.find_sites(kaset.sequence, e))
                assert n == 0, f"{key}: {e} kasetin içinde {n} yerde"


def test_tarif_uclari_farkli(tarif_sonuclari):
    """Çıkıntılar aynıysa insert ters de girer.

    Enzimlerin farklı olması yetmez: SalI/XhoI ikisi de 5' TCGA,
    XbaI/SpeI ikisi de 5' CTAG, BamHI/BclI ikisi de 5' GATC bırakır.
    Böyle bir çiftle yönlü klonlama yapamazsınız.
    """
    for key, (tarif, _rec, _k, res) in tarif_sonuclari.items():
        for adim in tarif.steps:
            assert (rb.ENZYMES[adim.left]["ovh"]
                    != rb.ENZYMES[adim.right]["ovh"]), f"{key}/{adim.title}"
        for basamak in res.steps:
            uyari = [n for n in basamak.notes if n.startswith("UYARI")]
            assert not uyari, f"{key}: {uyari}"


def test_uyumlu_uc_uyarisi_calisiyor(pacyc):
    """Uyumlu uç kontrolü gerçekten uyarıyor mu?

    Kontrolün kendisi sınanmazsa, tarifler onu geçiyor diye içinin boş
    olduğunu anlamayız.
    """
    adim = rb.cut_and_paste(pacyc, "SalI", "XhoI", "A" * 60)         if len(rb.find_sites(pacyc, "XhoI")) == 1 else None
    if adim is None:
        # pACYC184'te XhoI tek kesmiyorsa kontrolü doğrudan sına.
        assert rb.ENZYMES["SalI"]["ovh"] == rb.ENZYMES["XhoI"]["ovh"]
        return
    assert any("AYNI çıkıntıyı" in n for n in adim.notes)


def test_tarif_vektorleri_farkli():
    """Her tarif ayrı bir referans vektör kullanmalı.

    Aynı vektörü iki kez kullanmak tarif sayısını artırır ama
    öğretilen şeyi artırmaz.
    """
    vektorler = [r.vector_id for r in rb.RECIPES]
    assert len(vektorler) == len(set(vektorler)), vektorler


def test_tarif_anahtari_vektoru_yansitiyor():
    """Anahtar, kullanılan vektörle uyuşmalı (pbad24_kolera → PBAD24)."""
    for r in rb.RECIPES:
        onek = r.key.split("_")[0]
        assert r.vector_id.lower().startswith(onek[:5]),             f"{r.key} ile {r.vector_id} uyuşmuyor"


def test_her_tarifin_dersi_ve_gerekcesi_var():
    for r in rb.RECIPES:
        assert len(r.why_vector) > 80, r.key
        assert len(r.teaches) > 80, r.key
        assert r.source, r.key
        assert r.disease, r.key
        assert r.steps, r.key


def test_tarif_enzimleri_tek_kesiyor(tarif_sonuclari):
    """Her aşamada iki enzim de o aşamanın girdisini tek yerden kesmeli."""
    for key, (tarif, rec, _k, res) in tarif_sonuclari.items():
        girdi = rec["sequence"]
        for adim, basamak in zip(tarif.steps, res.steps):
            for e in (adim.left, adim.right):
                assert len(rb.find_sites(girdi, e)) == 1, f"{key}/{e}"
            girdi = basamak.sequence


def test_tarif_enzimleri_farkli():
    """Aynı enzim iki uçta kullanılırsa insert ters girebilir."""
    for tarif in rb.RECIPES:
        for adim in tarif.steps:
            assert adim.left != adim.right, tarif.key


def test_tarif_bolge_tasimasi_dizi_ile_dogrulanir(tarif_sonuclari):
    """Taşınan her bölge, ürün üzerinde gerçekten o diziyi göstermeli."""
    from core import db, ref_vectors

    parts = db.list_parts()
    for key, (tarif, rec, _k, res) in tarif_sonuclari.items():
        L0 = rec["length"]
        for f in ref_vectors.catalog_features(tarif.vector_id, parts):
            beklenen = _al(rec["sequence"], f.start, f.end)
            a, b, durum, L = f.start, f.end, "tam", L0
            atildi = False
            for basamak in res.steps:
                m = rb.map_feature(basamak, a, b, L)
                if m is None:
                    atildi = True
                    break
                a, b, d = m
                L = basamak.length
                if d == "kesildi":
                    durum = "kesildi"
            if atildi:
                continue
            gercek = _al(res.final, a, b)
            if durum == "tam":
                assert gercek == beklenen, f"{key}/{f.name}"
            else:
                assert gercek in beklenen, f"{key}/{f.name}"
