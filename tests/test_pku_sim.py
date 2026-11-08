"""pku_sim modelinin sınanması.

Model bir tahmindir; burada sınanan şey tahminin doğruluğu değil,
modelin kendi iç tutarlılığıdır: kalibrasyon hedefi tutturuluyor mu,
aktivite arttıkça Phe düşüyor mu, tedavisiz durum dengede mi.
"""

from core import pku_sim as S


def test_tedavisiz_denge():
    """Tedavi yokken günlük ortalama, varsayılan tabana oturmalı."""
    p = S.Params()
    r = S.simulate(p, 0.0)
    assert abs(r.son_gun_ortalama - S.BASELINE_PHE) < 5.0


def test_kalibrasyon_hedefi_tutturuyor():
    """Aktivite 1.0, klinik olarak bildirilen %40 düşüşü vermeli."""
    p = S.calibrate()
    taban = S.simulate(p, 0.0).aclik_phe
    tedavi = S.simulate(p, 1.0).aclik_phe
    dusus = (taban - tedavi) / taban
    assert abs(dusus - S.CLINICAL_DROP) < 0.01


def test_aktivite_arttikca_phe_duser():
    """Monotonluk: daha fazla aktivite, daha düşük açlık Phe'si."""
    p, runs = S.scenarios()
    degerler = [r.aclik_phe for r in runs]
    assert degerler == sorted(degerler, reverse=True)


def test_aclik_sabah_yediden_okunuyor():
    """Açlık değeri ilk öğünden (08:00) önceki saatten alınmalı."""
    p = S.Params()
    r = S.simulate(p, 0.0)
    # 07:00 değeri, aynı günün öğün sonrası tepesinden düşük olmalı
    assert r.aclik_phe < r.tepe


def test_bakteri_penceresi_ogunlerle_ortusuyor():
    """Bakteri yalnızca öğün sonrası pencerede çalışmalı."""
    p = S.Params()
    assert S._bakteri_acik(8.0 + p.bakteri_gecikme + 1.0, p)
    assert not S._bakteri_acik(6.0, p)       # sabah aç
    assert not S._bakteri_acik(3.0, p)       # gece


def test_vmax_sifirsa_tedavi_etkisiz():
    p = S.Params(vmax=0.0)
    assert abs(S.simulate(p, 1.0).aclik_phe - S.simulate(p, 0.0).aclik_phe) < 1e-9


def test_klinik_sabitler_kaynakla_uyumlu():
    """Vockley 2023'ten alınan sayılar bozulmadan dursun."""
    assert S.CLINICAL_DROP == 0.40
    assert S.CLINICAL_DROP_CI == (0.24, 0.52)
    assert S.SYNB1934_AUC_DROP > S.SYNB1618_AUC_DROP
