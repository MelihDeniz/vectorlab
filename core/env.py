"""Çalışma ortamı tespiti.

Uygulama iki şekilde çalışabilir:

  yerel    - tek kullanıcı, kendi bilgisayarında. Kaydedilen yapılar
             kalıcıdır ve kimseyle paylaşılmaz.
  barındırılan - Streamlit Community Cloud gibi bir sunucuda. Orada TEK
             bir kap ve TEK bir SQLite dosyası vardır: bütün ziyaretçiler
             aynı veritabanını paylaşır ve kap yeniden başlayınca
             (yeni sürüm yüklendiğinde ya da uygulama uykuya geçip
             uyandığında) dosya sıfırlanır.

Katalog her açılışta `parts_seed.py` içinden yeniden kurulduğu için
barındırılan ortamda da eksiksiz çalışır; kaybolabilecek tek şey
kullanıcının kaydettiği yapılar ve kaynakçadır. Arayüz bunu kaydetme
noktalarında söyler - sessizce veri kaybettirmek kabul edilemez.
"""

from __future__ import annotations

import os
from pathlib import Path

#: Streamlit Community Cloud depoyu buraya bağlar.
_CLOUD_MOUNT = Path("/mount/src")


def is_shared_deployment() -> bool:
    """Uygulama paylaşılan bir sunucuda mı çalışıyor?

    `VECTORLAB_SHARED` ortam değişkeni her şeyi geçersiz kılar; Streamlit
    Cloud'da Settings > Secrets altından ayarlanabilir. Otomatik tespit
    Streamlit Cloud'un bağlama noktasına bakar.
    """
    override = os.environ.get("VECTORLAB_SHARED", "").strip().lower()
    if override in ("1", "true", "evet", "yes"):
        return True
    if override in ("0", "false", "hayir", "hayır", "no"):
        return False
    return _CLOUD_MOUNT.exists()


SHARED_NOTICE = (
    "Bu uygulama paylaşılan bir sunucuda çalışıyor: kaydettiğiniz yapılar "
    "ve kaynakça **bu linki açan herkesle ortaktır** ve uygulama yeniden "
    "başladığında silinir. Tasarımınızı saklamak için GenBank ya da FASTA "
    "olarak indirin. Parça kataloğu her açılışta yeniden kurulduğu için "
    "ondan bir şey kaybolmaz."
)
