"""Sayfaların açılıp açıldığını sınar (duman testi).

Buradaki testler bir hesabın doğruluğunu değil, **uygulamanın
çalıştığını** sınar: her sayfa istisna atmadan açılıyor mu, bir sayfaya
girilen değer diğerine doğru geçiyor mu, seçiciler gerçekten bağlı mı.

Neden ayrı bir dosya: diğer testlerin hepsi `core/` içindeki modülleri
doğrudan çağırıyor ve sayfaların kendisine hiç dokunmuyordu. Bir
sayfada yazım hatası ya da yanlış bir `st.` çağrısı olsa test paketi
bunu görmüyordu - ancak uygulamayı elle açınca ortaya çıkıyordu. Bir
kullanıcı (ya da hoca) o sayfayı açtığında çıkan hata, testin
bulamadığı hatadır.

`streamlit.testing.v1.AppTest` sayfayı gerçekten çalıştırır, tarayıcı
olmadan.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

#: AppTest göreli yolu kendi dosyasına göre çözer, bu yüzden mutlak.
PAGES = Path(__file__).resolve().parent.parent / "app_pages"

#: `streamlit_app.py` içindeki paylaşılan durum varsayılanları.
#: Sayfalar tek başına çalıştırıldığında bunlar kurulu olmaz.
DEFAULTS = {
    "design_name": "pVL-01",
    "design_parts": [],
    "design_host": "BL21DE3",
    "design_insert": "",
    "design_insert_name": "hedefGen",
    "design_purpose": "ifade",
    "design_enzymes": [],
    "analyze_seq": "",
    "lit_results": [],
}

SAYFALAR = [
    "wizard", "designer", "disease", "refbuild", "vectors",
    "inspect", "analyze", "primers", "enzyme_lab", "labwork",
    "parts", "projects",
    "literature",
]

# Kısa ama gerçek bir kodlayan dizi: 189 bp, çerçevede, stop ile biter.
CDS = ("ATG" + "AAACTGAAATTTGCAGAAGTTCTGAAACCGGGTTCTAGCGGTAAAATT" * 3
       + "CTGGAAAAAGCGACCGTTATCGAAGGTCTGGAA" + "TAA")


def calistir(sayfa: str, durum: dict | None = None,
             timeout: float = 60) -> AppTest:
    at = AppTest.from_file(str(PAGES / f"{sayfa}.py"),
                           default_timeout=timeout)
    for k, v in DEFAULTS.items():
        at.session_state[k] = v
    for k, v in (durum or {}).items():
        at.session_state[k] = v
    return at.run()


@pytest.mark.parametrize("sayfa", SAYFALAR)
def test_sayfa_istisna_atmadan_aciliyor(sayfa):
    at = calistir(sayfa)
    assert not at.exception, (
        f"{sayfa} sayfası açılırken hata verdi: "
        + "; ".join(str(e.value) for e in at.exception))


@pytest.mark.parametrize("sayfa", SAYFALAR)
def test_sayfa_bos_ekran_dondurmuyor(sayfa):
    """Her sayfa en azından bir şey yazmalı.

    İstisna atmayan ama hiçbir şey göstermeyen bir sayfa da bozuktur.
    """
    at = calistir(sayfa)
    icerik = (len(at.markdown) + len(at.dataframe) + len(at.metric)
              + len(at.info) + len(at.warning) + len(at.error)
              + len(at.caption) + len(at.subheader) + len(at.tabs))
    assert icerik > 0, f"{sayfa} sayfası hiçbir şey göstermedi"


# ---------------------------------------------------------------------------
# Kodon organizma seçicisi
# ---------------------------------------------------------------------------


def test_kodon_organizma_secicisi_varsayilani_konaktan_aliyor():
    """Seçici, tasarımın konağının organizmasına açılmalı.

    Konak BL21(DE3) ise E. coli, L. lactis konağı seçiliyse llactis.
    Bu bağ kopuksa kullanıcı yanlış tabloyla optimize eder.
    """
    at = calistir("analyze", {"analyze_seq": CDS})
    assert not at.exception
    assert at.session_state["kodon_org"] == "ecoli"

    at2 = calistir("analyze", {"analyze_seq": CDS,
                               "design_host": "LLNZ9000"})
    assert not at2.exception
    assert at2.session_state["kodon_org"] == "llactis", (
        "L. lactis konağı seçiliyken seçici hâlâ E. coli'de")


def test_kodon_organizmasi_degisince_olcum_degisiyor():
    """Seçiciyi değiştirmek gerçekten başka bir tablo kullanmalı.

    Sayfa aynı sayıyı veriyorsa seçici bağlı değil demektir - hata tam
    olarak buydu: tek bir E. coli tablosu vardı ve seçim yoktu.
    """
    e = calistir("analyze", {"analyze_seq": CDS, "kodon_org": "ecoli"})
    ll = calistir("analyze", {"analyze_seq": CDS, "kodon_org": "llactis"})
    assert not e.exception and not ll.exception

    def cai_degeri(at: AppTest) -> float:
        for m in at.metric:
            if m.label == "CAI":
                return float(m.value)
        raise AssertionError("CAI ölçümü sayfada yok")

    assert cai_degeri(e) != cai_degeri(ll), (
        f"E. coli ve L. cremoris aynı CAI'yi verdi: {cai_degeri(e)}")


# ---------------------------------------------------------------------------
# Projeler sayfası
# ---------------------------------------------------------------------------


@pytest.fixture
def kayitli_yapi():
    """Geçici bir kayıtlı tasarım; testten sonra silinir."""
    from core import db

    db.init_db()
    ad = "pTestProje"
    db.delete_construct(ad)
    db.save_construct(ad, "BL21DE3",
                      ["PRO_T7LAC", "RBS_PET", "ORI_COLE1", "MRK_KAN"],
                      CDS, "testGen", "ilk kayıt")
    yield ad
    db.delete_construct(ad)


def test_projeler_sayfasi_kayitli_yapiyi_gosteriyor(kayitli_yapi):
    at = calistir("projects", {"proje_secim": kayitli_yapi})
    assert not at.exception
    metin = " ".join(d.value for d in at.caption) + " ".join(
        m.value for m in at.markdown)
    assert kayitli_yapi in str(at.session_state["proje_secim"])
    # Yapı kurulduysa uzunluk ölçümü dolu olmalı.
    uzunluk = [m.value for m in at.metric if m.label == "Uzunluk"]
    assert uzunluk and uzunluk[0] != "—", (
        f"yapı kurulamadı: {uzunluk}")


def test_projeler_bos_durumda_yol_gosteriyor(tmp_path, monkeypatch):
    """Hiç kayıt yoksa sayfa boş ekran değil, yönlendirme vermeli.

    Geçici bir veritabanı kullanılır. Kullanıcının kendi kayıtlarını
    silip "artık boş" demek, testin veri kaybettirmesi olurdu -
    `DB_PATH` `connect()` içinde okunduğu için yolu değiştirmek yeterli.
    """
    from core import db

    monkeypatch.setattr(db, "DB_PATH", tmp_path / "bos.db")
    db.init_db()
    assert db.list_constructs() == []

    at = calistir("projects")
    assert not at.exception
    metin = " ".join(i.value for i in at.info)
    assert "Tasarımcı" in metin, "boş durumda yol gösterilmedi"


def test_surum_gecmisi_her_degisikligi_yaziyor(kayitli_yapi):
    """Değişiklik yeni sürüm bırakmalı, aynı kayıt bırakmamalı."""
    from core import db

    assert len(db.list_construct_versions(kayitli_yapi)) == 1
    # Aynısını tekrar kaydetmek sürüm açmamalı.
    assert db.save_construct(
        kayitli_yapi, "BL21DE3",
        ["PRO_T7LAC", "RBS_PET", "ORI_COLE1", "MRK_KAN"],
        CDS, "testGen", "ilk kayıt") is False
    assert len(db.list_construct_versions(kayitli_yapi)) == 1

    # Parça eklemek sürüm açmalı.
    assert db.save_construct(
        kayitli_yapi, "BL21DE3",
        ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"],
        CDS, "testGen", "terminatör eklendi") is True
    surumler = db.list_construct_versions(kayitli_yapi)
    assert len(surumler) == 2
    # Her sürüm kendi damgasını taşımalı: tekrarlanabilirliğin dayanağı.
    assert all(s["provenance"] for s in surumler)

    at = calistir("projects", {"proje_secim": kayitli_yapi})
    assert not at.exception
    fark = " ".join(m.value for m in at.markdown)
    assert "TER_T7" in fark, "son iki sürüm arasındaki fark yazılmadı"


def test_kayitli_yapi_tasarimciya_geri_yuklenebiliyor(kayitli_yapi):
    """Geri yükleme Tasarımcı'nın widget anahtarlarını doldurmalı.

    `design_parts` yazmak yetmez: key verilmiş bir multiselect kendi
    durumunu korur ve `default` parametresini yok sayar. Bu test o
    ayrımı sabitliyor.
    """
    from core import db
    from core.parts_seed import DESIGNER_CATEGORIES

    kayit = db.get_construct(kayitli_yapi)
    parts_by_id = {p["part_id"]: p for p in db.list_parts()}

    # Projeler sayfasının yaptığı eşlemenin aynısı.
    durum = {}
    for cat in DESIGNER_CATEGORIES:
        durum[f"ms_{cat}"] = [
            pid for pid in kayit["part_ids"]
            if parts_by_id.get(pid, {}).get("category") == cat
        ]
    # Kayıttaki her parça bir kategoriye düşmeli; düşmeyen parça
    # sessizce kaybolur.
    dusen = {pid for v in durum.values() for pid in v}
    assert dusen == set(kayit["part_ids"]), (
        f"kategoriye düşmeyen parça: {set(kayit['part_ids']) - dusen}")

    at = calistir("designer", {**durum,
                               "design_name": kayit["name"],
                               "design_host": kayit["host_id"]})
    assert not at.exception
    assert set(at.session_state["design_parts"]) == set(kayit["part_ids"])


def test_konak_tablodan_farkliysa_sayfa_soyluyor():
    """Konak ile seçilen tablo farklıysa kullanıcı bunu görmeli.

    Sessizce başka bir organizmaya göre optimize etmek kabul edilemez.
    """
    at = calistir("analyze", {"analyze_seq": CDS,
                              "design_host": "LLNZ9000",
                              "kodon_org": "ecoli"})
    assert not at.exception
    metinler = " ".join(i.value for i in at.info)
    assert "Escherichia coli" in metinler or "E. coli" in metinler, (
        "konak-tablo farkı bildirilmedi")
