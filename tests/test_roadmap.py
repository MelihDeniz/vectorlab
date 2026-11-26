"""core/roadmap.py: gerçekleştirme yol haritası.

Bu modülün yaptığı iş, bir listeyi ekrana basmak değil; bir **sıra**
iddiasında bulunmak. Yol haritasının değeri de riski de buradan
gelir: sıra yanlışsa, kullanıcı etik kurul onayı almadan hayvan
çalışmasına başlar ya da içerme ölçmeden canlı GDO'yu hayvana verir.
Buradaki sınavların çoğu bu yüzden **sıra** ve **koşul** sınavıdır,
metin sınavı değil.

İkinci risk: yol haritasının statik bir metne dönüşmesi. Tasarımdan
türemiyorsa hiçbir değeri yok - herkes genel bir "önce klonla sonra
ifade et" listesi yazabilir. `test_gorevler_tasarimdan_turuyor` ve
kardeşleri bunu koruyor.
"""

from __future__ import annotations

import pytest

from core import db, disease_seed as ds, roadmap as rm
from core.assembly import assemble, sort_canonical
from core.validate import validate

IBD_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
             "CDS_IL10", "TER_B0015", "MRK_THYA_LL"]
LAB_KASET = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_PUC", "MRK_KAN"]


@pytest.fixture(scope="module", autouse=True)
def _db():
    db.init_db()


@pytest.fixture(scope="module")
def terapotik() -> rm.Roadmap:
    h = ds.by_key("ibd")
    vec = assemble("pIBD", sort_canonical(db.get_parts(IBD_KASET)))
    b = validate(vec, host=db.get_host("LLNZ9000"), purpose=h.purpose)
    return rm.build("pIBD", IBD_KASET, host_id="LLNZ9000", disease=h,
                    findings=b)


@pytest.fixture(scope="module")
def laboratuvar() -> rm.Roadmap:
    return rm.build("pLab", LAB_KASET, host_id="BL21DE3")


def _gorevler(yol: rm.Roadmap) -> list[rm.Task]:
    return [t for p in yol.phases for t in p.tasks]


def _sira(yol: rm.Roadmap, key: str) -> int:
    """Görevin yol haritasındaki sırası; yoksa -1."""
    for i, g in enumerate(_gorevler(yol)):
        if g.key == key:
            return i
    return -1


# ---------------------------------------------------------------------------
# Sıra: bu modülün asıl iddiası
# ---------------------------------------------------------------------------


def test_etik_onayi_her_hayvan_isinden_once(terapotik):
    """Onaysız yapılan hayvan çalışması yayımlanamaz.

    Etik kurul başvurusu yol haritasında hayvanla ilgili HER işten
    önce gelmeli. Bu sınav metne değil **sıraya** bakar: başlıkta
    "önce etik kurul" yazıp görevi sona koymak işe yaramaz.
    """
    etik = _sira(terapotik, "etik")
    assert etik >= 0, "etik kurul görevi yok"
    for anahtar in ("model", "uygula", "olc"):
        s = _sira(terapotik, anahtar)
        assert s > etik, f"{anahtar} etik onayından önce geliyor"


def test_icerme_olculmeden_hayvana_gecilmiyor(terapotik):
    """Canlı GDO'yu hayvana vermeden önce kaçak oranı ölçülmeli."""
    icerme = _sira(terapotik, "icerme")
    uygula = _sira(terapotik, "uygula")
    assert icerme >= 0 and uygula >= 0
    assert icerme < uygula


def test_ifade_gosterilmeden_islev_aranmiyor(terapotik):
    """Protein üretilmeden etkinlik ölçmek sonucu yorumlanamaz kılar."""
    assert _sira(terapotik, "ifade") < _sira(terapotik, "etkinlik")


def test_dizileme_her_seyden_once(terapotik):
    """Yanlış bir yapıyla yapılan altı aylık iş boşa gider."""
    dizi = _sira(terapotik, "dizile")
    assert dizi >= 0
    for g in _gorevler(terapotik)[dizi + 1:]:
        assert g.key not in ("klonla",)
    # Dizilemeden sonra gelen her şey onun doğruluğuna dayanır.
    assert _sira(terapotik, "ifade") > dizi


def test_asamalarin_kapisi_yazili(terapotik):
    """Bir aşamadan ötekine geçiş koşulu tanımlı olmalı.

    Kapısı olmayan aşama, ne zaman bittiği belli olmayan aşamadır.
    """
    for p in terapotik.phases[:-1]:
        assert p.kapi, f"Aşama {p.no} ({p.title}): geçiş koşulu yok"


# ---------------------------------------------------------------------------
# Tasarımdan türeme
# ---------------------------------------------------------------------------


def test_sinyal_peptidi_salgilanma_gorevi_ekliyor(terapotik):
    """Kasette sinyal varsa salgılanma AYRICA doğrulanmalı."""
    g = next((x for x in _gorevler(terapotik) if x.key == "salgi"), None)
    assert g is not None, "sinyal peptidi var ama salgılanma görevi yok"
    assert "SIG_USP45" in g.kosul
    # Lizis kontrolü olmadan salgılanma iddiası kurulamaz.
    assert "liz" in g.nasil.lower()


def test_sinyalsiz_kasette_salgilanma_gorevi_yok(laboratuvar):
    assert _sira(laboratuvar, "salgi") == -1


def test_gram_pozitif_konak_donusum_gorevi_ekliyor(terapotik):
    g = next((x for x in _gorevler(terapotik) if x.key == "aktar"), None)
    assert g is not None
    assert "Lactococcus" in g.kosul or "Gram-pozitif" in g.kosul
    assert "elektroporasyon" in g.nasil.lower()


def test_terapotik_olmayan_tasarimda_in_vivo_yok(laboratuvar):
    """Sıradan bir klonlama vektörü için hayvan modeli anlamsız."""
    assert laboratuvar.phase("invivo") is None
    assert _sira(laboratuvar, "etik") == -1
    assert _sira(laboratuvar, "icerme") == -1
    # Neden eklenmediği SÖYLENMELİ; sessizce atlamak yanıltır.
    assert any("terapötik bir bağlamda işaretlenmedi" in n
               for n in laboratuvar.notes)


def test_in_vivo_kapatilinca_asama_ve_izin_dusuyor():
    h = ds.by_key("ibd")
    acik = rm.build("p", IBD_KASET, host_id="LLNZ9000", disease=h)
    kapali = rm.build("p", IBD_KASET, host_id="LLNZ9000", disease=h,
                      in_vivo=False)
    assert acik.phase("invivo") is not None
    assert kapali.phase("invivo") is None
    assert kapali.gun < acik.gun
    assert "hadyek" not in {a.key for a in kapali.approvals}


def test_hastalik_biyobelirteci_in_vivo_okumasi_oluyor(terapotik):
    """Başarı ölçütü uydurulmaz: kayıt defterinden gelir."""
    h = ds.by_key("ibd")
    assert terapotik.readout == h.biomarker
    g = next(x for x in _gorevler(terapotik) if x.key == "olc")
    assert h.biomarker in g.nasil


def test_tezgah_protokolu_tekrarlanmiyor(terapotik):
    """Kesim/ligasyon/transformasyon `core/protocol.py` içinde.

    İki yerde iki ayrı protokol tutmak, biri değiştiğinde öbürünün
    sessizce yanlış kalması demektir. Birinci aşama o protokole
    YÖNLENDİRİR, adımlarını yeniden yazmaz.
    """
    g = next(x for x in _gorevler(terapotik) if x.key == "klonla")
    assert "Laboratuvar sayfas" in g.nasil
    # Reaksiyon kurulumu burada OLMAMALI.
    metin = " ".join(x.nasil for x in _gorevler(terapotik)).lower()
    for ayrinti in ("µl", "ng/µl", "16 °c", "37 °c bir gece"):
        assert ayrinti not in metin, f"tezgâh ayrıntısı sızmış: {ayrinti}"


# ---------------------------------------------------------------------------
# Riskler
# ---------------------------------------------------------------------------


def test_denetim_bulgulari_riske_cevriliyor():
    """Denetimde görünen eksik, başvuruda risk olarak durmalı."""
    kaset = ["PRO_T7LAC", "RBS_PET", "CDS_IL10"]       # orijin ve marker yok
    vec = assemble("x", sort_canonical(db.get_parts(kaset)))
    b = validate(vec, host=db.get_host("BL21DE3"))
    hatalar = {f.rule_id for f in b if f.severity in ("hata", "uyari")}
    assert hatalar, "sınav için hata üreten bir kaset gerekiyordu"
    yol = rm.build("x", kaset, host_id="BL21DE3", findings=b)
    risk_kimlikleri = {k.key for k in yol.risks}
    assert hatalar <= risk_kimlikleri, hatalar - risk_kimlikleri


def test_antibiyotik_markeri_risk_aciyor(laboratuvar):
    k = next((x for x in laboratuvar.risks if x.key == "ANTIBIYOTIK"), None)
    assert k is not None
    assert "MRK_KAN" in k.etki
    assert "yatay gen aktarım" in k.etki
    # Önlem somut olmalı: ne yapılacağı ve ne zaman.
    assert "antibiyotiksiz" in k.onlem


def test_insan_geni_kodon_riski_aciyor(terapotik):
    k = next((x for x in terapotik.risks if x.key == "INSAN_CDS"), None)
    assert k is not None
    assert "CDS_IL10" in k.etki


def test_bilinen_acik_riske_tasiniyor(terapotik):
    h = ds.by_key("ibd")
    k = next((x for x in terapotik.risks if x.key == "BILINEN_ACIK"), None)
    assert k is not None
    assert k.etki == h.gap


def test_her_riskin_onlemi_var(terapotik, laboratuvar):
    """Önlemsiz risk, başvuruda işe yaramaz."""
    for yol in (terapotik, laboratuvar):
        for k in yol.risks:
            assert k.onlem and len(k.onlem) > 30, k.title


# ---------------------------------------------------------------------------
# İzinler
# ---------------------------------------------------------------------------


def test_izinlerin_yasal_dayanagi_kunyeli():
    """Dayanak künyesiz yazılırsa başvuruda kullanılamaz."""
    for a in rm.APPROVALS.values():
        assert a.body and a.basis and a.when, a.key
        assert len(a.basis) > 25, a.key


def test_mevzuatin_degistigi_sayfada_yaziyor():
    """Biyogüvenlik Kurulu 2018'de kaldırıldı.

    Bu, ezberden yazılan bir başvurunun var olmayan bir kurula
    başvurmasına yol açar. Uyarı metinde DURMALI; silinirse bu sınav
    düşer.
    """
    gdo = rm.APPROVALS["gdo"]
    assert "2018" in gdo.caution
    assert "703" in gdo.caution
    assert "teyit" in gdo.caution.lower()


def test_hadyek_dayanagi_resmi_gazete_kunyesiyle():
    h = rm.APPROVALS["hadyek"]
    assert "28914" in h.basis and "2014" in h.basis
    assert "5199" in h.basis


def test_izinler_yalnizca_gerekenler(terapotik, laboratuvar):
    """Sıradan bir klonlama vektörüne klinik izni koymak yanıltır."""
    t_izin = {a.key for a in terapotik.approvals}
    l_izin = {a.key for a in laboratuvar.approvals}
    assert {"hadyek", "gdo", "klinik"} <= t_izin
    assert "hadyek" not in l_izin
    assert "gdo" not in l_izin
    assert "klinik" not in l_izin
    assert "ibg" in l_izin          # laboratuvar çalışması yine de izinli


# ---------------------------------------------------------------------------
# Görev bütünlüğü ve çıktı
# ---------------------------------------------------------------------------


def test_her_gorevin_nasil_gereken_cikti_alanlari_dolu(terapotik,
                                                       laboratuvar):
    """Çıktısı tanımlanmayan görev, ne zaman bittiği belli olmayan görev."""
    for yol in (terapotik, laboratuvar):
        for g in _gorevler(yol):
            assert g.nasil and len(g.nasil) > 40, g.key
            assert g.gereken, g.key
            assert g.cikti and len(g.cikti) > 10, g.key
            assert g.gun > 0, g.key


def test_gorev_anahtarlari_tekil(terapotik):
    anahtarlar = [g.key for g in _gorevler(terapotik)]
    assert len(anahtarlar) == len(set(anahtarlar))


def test_sure_toplami_tutuyor(terapotik):
    assert terapotik.gun == sum(p.gun for p in terapotik.phases)
    assert terapotik.gun == sum(g.gun for g in _gorevler(terapotik))
    assert terapotik.ay == round(terapotik.gun / 21.0, 1)


def test_sureler_disaridan_degistirilebiliyor():
    """Tahmin olduğu için ayarlanabilir olmalı."""
    normal = rm.build("x", LAB_KASET, host_id="BL21DE3")
    hizli = rm.build("x", LAB_KASET, host_id="BL21DE3",
                     phase_days={"klonlama": 5.0})
    assert hizli.gun < normal.gun


def test_sureler_tahmin_oldugunu_soyluyor(terapotik):
    metin = " ".join(terapotik.notes)
    assert "tahmin" in metin.lower()
    assert "paralelleştirme" in metin


def test_markdown_butun_asamalari_basiyor(terapotik):
    md = rm.to_markdown(terapotik)
    for p in terapotik.phases:
        assert p.title in md
    for g in _gorevler(terapotik):
        assert g.title in md
    assert "İzinler ve yasal dayanaklar" in md
    assert "Risk yönetimi" in md


def test_basvuru_taslagi_kullanicinin_yazacagini_isaretliyor(laboratuvar):
    """Aracın yazdığı yerle sizin yazacağınız yer karışmamalı."""
    metin = rm.to_proposal(laboratuvar)
    assert "[" in metin and "]" in metin, "doldurulacak yer işaretlenmemiş"
    assert "Özet" in metin and "Risk yönetimi" in metin
    assert "İş paketleri" in metin or "iş paketleri" in metin
    # Danışmanlık olmadığı taslağın kendisinde yazılı olmalı.
    assert "danışmanlık" in metin


def test_basvuru_taslagi_hastaliktan_icerik_aliyor(terapotik):
    h = ds.by_key("ibd")
    metin = rm.to_proposal(terapotik, disease=h)
    assert h.name in metin
    assert h.logic in metin
    # Kaynakça uydurulmaz: kayıt defterindeki atıflar.
    for makale in h.papers:
        assert makale in metin


def test_basvuru_taslagi_is_paketlerini_asamalardan_kuruyor(terapotik):
    metin = rm.to_proposal(terapotik)
    for p in terapotik.phases:
        assert f"İP{p.no}" in metin
    assert f"{terapotik.gun:g}" in metin
