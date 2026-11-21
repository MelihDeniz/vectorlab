"""Arayüz dili: sözlük kodla tutarlı mı, çeviri gerçekten uygulanıyor mu.

Kaynak metni anahtar yapmanın bilinen bedeli şudur: Türkçe bir etiket
değişirse çeviri bağı **sessizce** kopar ve İngilizce kipte o satır
Türkçe görünmeye başlar. Buradaki testler o sessizliği kaldırıyor.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from core import i18n

KOK = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(KOK / "tools"))

import i18n_scan  # noqa: E402


@pytest.fixture(autouse=True)
def turkceye_don():
    """Her testten sonra dili varsayılana döndür.

    `i18n` etkin dili modül düzeyinde tutuyor; bir test İngilizceye
    geçip bırakırsa sonraki testler yanlış dilde çalışır.
    """
    yield
    i18n.set_language(i18n.DEFAULT_LANG)


def test_sozlukte_eksik_ceviri_yok():
    metinler = i18n_scan.kaynak_metinler()
    assert metinler, "koddan hiç çevrilebilir metin çıkmadı"
    sozluk = i18n.catalog("en")
    eksik = [m for m in metinler if m not in sozluk]
    assert not eksik, (
        f"{len(eksik)} metnin İngilizcesi yok. "
        f"`python tools/i18n_scan.py --yaz` ile ekleyin. "
        f"İlk üç: {eksik[:3]}")


def test_sozlukte_artik_kullanilmayan_kayit_yok():
    """Türkçe metin değişmişse sözlükte öksüz kayıt kalır.

    Öksüz kayıt zararsız görünür ama kalan şey, o satırın artık
    çevrilmediğidir: ekranda Türkçe görünür ve kimse fark etmez.
    """
    kodda = set(i18n_scan.kaynak_metinler())
    sozluk = i18n.catalog("en")
    oksuz = sorted(k for k in sozluk if k not in kodda)
    assert not oksuz, (
        f"{len(oksuz)} sözlük kaydı kodda bulunamıyor - Türkçe metin "
        f"değişmiş olabilir. `python tools/i18n_scan.py --temizle` "
        f"ile silin. İlk üç: {oksuz[:3]}")


def test_ceviri_gercekten_uygulaniyor():
    i18n.set_language("en")
    assert i18n.t("Yapı adı") == "Construct name"
    assert i18n.t("Sürüm geçmişi") == "Version history"
    assert i18n.t("Parça kataloğu") == "Part catalog"


def test_turkce_kipte_metin_degismiyor():
    i18n.set_language("tr")
    for m in ("Yapı adı", "Sürüm geçmişi", "bilinmeyen metin"):
        assert i18n.t(m) == m


def test_karsiligi_olmayan_metin_turkce_kaliyor():
    """Eksik çeviri, bozuk ekran değil çalışan Türkçe bırakmalı.

    Soyut anahtar kullanılsaydı burada anahtarın kendisi görünürdü
    (`ui.designer.name`); kaynak metni anahtar yapmanın asıl sebebi bu.
    """
    i18n.set_language("en")
    assert i18n.t("Hiç çevrilmemiş bir cümle.") == "Hiç çevrilmemiş bir cümle."


def test_bilinmeyen_dil_varsayilana_duser():
    i18n.set_language("klingon")
    assert i18n.get_language() == i18n.DEFAULT_LANG


def test_bos_ceviriler_sozluge_girmiyor():
    """Boş karşılık 'henüz çevrilmedi' demektir, boş ekran değil."""
    ham = json.loads(
        (KOK / "data" / "i18n" / "en.json").read_text(encoding="utf-8"))
    sozluk = i18n.catalog("en")
    for k, v in ham.items():
        if not v.strip():
            assert k not in sozluk


def test_ceviriler_bicim_yeri_tasimiyor():
    """Sarılan metinlerde `{...}` olmamalı.

    Biçim yeri içeren metinler elle çevrilir; sözlüğe girmiş olsalardı
    `t()` onları olduğu gibi döndürür ve `.format()` çağrılmadığı için
    ekranda süslü parantez görünürdü.
    """
    for kaynak, ceviri in i18n.catalog("en").items():
        assert "{" not in kaynak, f"kaynakta biçim yeri: {kaynak[:50]}"
        assert "{" not in ceviri, f"çeviride biçim yeri: {ceviri[:50]}"


def test_kapsam_hesabi_dogru():
    metinler = i18n_scan.kaynak_metinler()
    cevrilen, toplam = i18n.coverage("en", metinler)
    assert toplam == len(metinler)
    assert cevrilen == toplam, (
        f"İngilizce kapsam %{100 * cevrilen / toplam:.1f}")
    # Türkçe kaynak dildir: her zaman tam.
    assert i18n.coverage("tr", metinler) == (len(metinler), len(metinler))
