# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller tarifi: tek klasörlük Windows kurulumu.

    pip install -e ".[desktop,build]"
    pyinstaller vectorlab.spec

Çıktı `dist/VectorLab/VectorLab.exe` olur. Python kurulu olmayan bir
makinede çalışır; müşterinin uygulamayı deneyebilmesi için gereken şey
budur.

**Neden onedir, onefile değil.** Tek dosyalık .exe her açılışta kendini
geçici bir klasöre açar; bu uygulamada `data/` altında 136 veri dosyası
var ve açılış belirgin biçimde yavaşlıyor. Ayrıca uygulama kendi
`data/` klasörüne SQLite veritabanı yazıyor: tek dosyalık kipte o yol
her çalıştırmada değişen geçici bir klasör olur ve kullanıcının
kaydettiği yapılar kaybolur. Klasörlük kipte veritabanı .exe'nin
yanında durur ve kalıcıdır.

**Streamlit neden elle toplanıyor.** Streamlit, sayfalarını ve statik
dosyalarını çalışma anında kendi paket klasöründen okur; PyInstaller
bunları kendiliğinden bulamaz, `collect_all` ile açıkça toplanır.
"""

from PyInstaller.utils.hooks import collect_all, collect_data_files

datas = []
binaries = []
hiddenimports = []

# Streamlit'in statik dosyaları (frontend) ve metadata'sı şart.
for paket in ("streamlit", "altair", "pydeck"):
    p_datas, p_binaries, p_hidden = collect_all(paket)
    datas += p_datas
    binaries += p_binaries
    hiddenimports += p_hidden

# Biopython'un veri dosyaları (REBASE enzim tabloları burada).
datas += collect_data_files("Bio")

# Uygulamanın kendi dosyaları. `data/` içindeki JSON'lar olmadan
# katalog kurulamaz; `.streamlit/config.toml` tema ve sunucu ayarı.
datas += [
    ("streamlit_app.py", "."),
    ("app_pages", "app_pages"),
    ("core", "core"),
    ("data", "data"),
    (".streamlit", ".streamlit"),
    ("README.md", "."),
    ("LICENSE", "."),
    ("DISCLAIMER.md", "."),
    ("THIRD-PARTY-NOTICES.md", "."),
]

# Streamlit sayfaları ve çekme betikleri dinamik yüklendiği için
# statik çözümleyici göremez.
hiddenimports += [
    "streamlit.web.cli",
    "streamlit.runtime.scriptrunner.magic_funcs",
    "webview",
    "core.notice",
]

a = Analysis(
    ["desktop.py"],
    pathex=[],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=["matplotlib", "tkinter", "PyQt5", "PySide6"],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="VectorLab",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    # console=False: kullanıcı siyah konsol penceresi görmesin.
    # Hata olursa desktop.py zaten MessageBox açıyor ve
    # data/desktop.log dosyasına yazıyor.
    console=False,
    icon="assets/vectorlab.ico" if __import__("pathlib").Path(
        "assets/vectorlab.ico").exists() else None,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="VectorLab",
)
