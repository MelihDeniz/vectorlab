"""Arayüz dili: Türkçe ve İngilizce.

**Anahtar Türkçe metnin kendisidir.** `t("Yapı adı")` çağrısı Türkçe
kipte aynı dizeyi döndürür, İngilizce kipte sözlükten karşılığını
arar, bulamazsa yine Türkçeyi döndürür.

Neden soyut anahtar (`designer.name_label`) kullanılmadı: bu uygulama
Türkçe yazıldı ve Türkçe birincil dil olarak kalıyor. Soyut anahtara
geçmek 2000'den fazla çağrı yerini okunmaz hale getirirdi
(`t("ui.designer.field.name")` kodu okuyan için hiçbir şey söylemez) ve
çeviri eksikse ekranda anahtarın kendisi görünürdü. Kaynak metni
anahtar yapmak, çeviri eksik olduğunda **çalışan Türkçe** bırakır -
bozuk bir ekran değil.

Bedeli şudur: Türkçe metin değişirse çeviri bağı kopar. `tools/
i18n_scan.py` bunu yakalar - sözlükte olup kodda artık bulunmayan
kayıtları "artık kullanılmıyor" diye bildirir.

Dil, uygulamanın her çalışmasında `streamlit_app.py` tarafından
ayarlanır; bu modül Streamlit'e bağımlı değildir, çünkü `core/`
içindeki hiçbir şey arayüze bağımlı olmamalı (testler ve çekme
betikleri de bu modülleri kullanıyor).
"""

from __future__ import annotations

import functools
import json
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "i18n"

#: Desteklenen diller: kod -> o dildeki kendi adı.
LANGS = {"tr": "Türkçe", "en": "English"}

DEFAULT_LANG = "tr"

_current = DEFAULT_LANG


def set_language(lang: str) -> None:
    """Etkin dili ayarla. Tanınmayan kod varsayılana döner."""
    global _current
    _current = lang if lang in LANGS else DEFAULT_LANG


def get_language() -> str:
    return _current


@functools.lru_cache(maxsize=None)
def catalog(lang: str) -> dict[str, str]:
    """Bir dilin sözlüğü: Türkçe kaynak metin -> çeviri.

    Türkçe için sözlük yoktur (kaynak dil odur). Dosya yoksa boş
    sözlük döner ve her şey Türkçe görünür - uygulama yine çalışır.
    """
    if lang == DEFAULT_LANG:
        return {}
    yol = DATA_DIR / f"{lang}.json"
    if not yol.exists():
        return {}
    ham = json.loads(yol.read_text(encoding="utf-8"))
    # Boş karşılıklar "henüz çevrilmedi" demektir; sözlüğe girmezler ki
    # `t()` kaynağa düşsün, ekrana boşluk basmasın.
    return {k: v for k, v in ham.items() if isinstance(v, str) and v.strip()}


def t(metin: str) -> str:
    """Metni etkin dile çevir; karşılığı yoksa olduğu gibi bırak."""
    if _current == DEFAULT_LANG:
        return metin
    return catalog(_current).get(metin, metin)


def coverage(lang: str, kaynaklar: list[str]) -> tuple[int, int]:
    """(çevrilen, toplam) - arayüzün ne kadarının çevrildiği.

    `kaynaklar`, koddan çıkarılmış çevrilebilir metinlerin listesidir
    (`tools/i18n_scan.py` üretir). Sayı uydurulmaz: sözlükte karşılığı
    olan kaç kaynak metin varsa o sayılır.
    """
    if lang == DEFAULT_LANG:
        return len(kaynaklar), len(kaynaklar)
    sozluk = catalog(lang)
    return sum(1 for k in kaynaklar if k in sozluk), len(kaynaklar)
