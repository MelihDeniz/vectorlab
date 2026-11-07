"""VectorLab'i masaüstü uygulaması olarak çalıştırır.

Streamlit bir web sunucusudur ve normalde varsayılan tarayıcıda açılır.
Bu betik onu arka planda, görünmez bir şekilde başlatır ve kendi
penceresinde gösterir: adres çubuğu yok, sekme yok, tarayıcı menüsü
yok. Kullanıcı açısından sıradan bir Windows uygulamasıdır.

Pencere, Windows'un WebView2 bileşenini kullanır - bu, Edge
tarayıcısının kendisi değil, uygulamaların içine gömülen görüntüleme
motorudur (Teams ve VS Code gibi uygulamalar da aynısını kullanır).
Windows 11'de hazır gelir, ayrıca kurulum gerekmez.

Çalıştırmak:
    pythonw desktop.py      (konsol penceresi açılmaz)
    python desktop.py       (hata ayıklarken; günlük konsola da basar)
"""

from __future__ import annotations

import ctypes
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

KOK = Path(__file__).resolve().parent
UYGULAMA = KOK / "streamlit_app.py"
GUNLUK = KOK / "data" / "desktop.log"

BASLIK = "VectorLab"
GENISLIK, YUKSEKLIK = 1440, 900
MIN_GENISLIK, MIN_YUKSEKLIK = 1024, 680
BASLATMA_SURESI = 90          # saniye; ilk açılışta veritabanı kurulur


def gunlukle(mesaj: str) -> None:
    satir = f"{time.strftime('%Y-%m-%d %H:%M:%S')}  {mesaj}"
    print(satir)
    try:
        GUNLUK.parent.mkdir(parents=True, exist_ok=True)
        with GUNLUK.open("a", encoding="utf-8") as f:
            f.write(satir + "\n")
    except OSError:
        pass


def hata_kutusu(mesaj: str) -> None:
    """pythonw ile çalışırken konsol olmadığı için hata görünmez kalır."""
    gunlukle(f"HATA: {mesaj}")
    try:
        ctypes.windll.user32.MessageBoxW(
            None, f"{mesaj}\n\nAyrıntı: {GUNLUK}", f"{BASLIK} - hata", 0x10)
    except Exception:
        pass


def bos_port() -> int:
    """İşletim sisteminden kullanılmayan bir port iste.

    Sabit 8501 kullanmıyoruz: uygulama zaten açıkken ikinci kez
    başlatılırsa ya da başka bir Streamlit çalışıyorsa çakışır ve
    pencere boş gelir.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def sunucuyu_baslat(port: int) -> subprocess.Popen:
    komut = [
        sys.executable, "-m", "streamlit", "run", str(UYGULAMA),
        "--server.port", str(port),
        "--server.address", "127.0.0.1",
        "--server.headless", "true",        # tarayıcıyı kendisi açmasın
        "--browser.gatherUsageStats", "false",
        "--server.fileWatcherType", "none",  # masaüstünde yeniden yükleme yok
        "--global.developmentMode", "false",
    ]
    gunlukle(f"sunucu baslatiliyor: port {port}")
    bayraklar = 0
    if sys.platform == "win32":
        # Alt süreç kendi konsol penceresini açmasın.
        bayraklar = getattr(subprocess, "CREATE_NO_WINDOW", 0)
    return subprocess.Popen(
        komut, cwd=str(KOK), creationflags=bayraklar,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        encoding="utf-8", errors="replace",
    )


def hazir_mi(port: int, surec: subprocess.Popen) -> bool:
    """Sunucu yanıt verene kadar bekle; süreç ölürse erken çık."""
    adres = f"http://127.0.0.1:{port}/_stcore/health"
    bitis = time.time() + BASLATMA_SURESI
    while time.time() < bitis:
        if surec.poll() is not None:
            cikti = (surec.stdout.read() if surec.stdout else "") or ""
            gunlukle(f"sunucu beklenmedik sekilde kapandi:\n{cikti[-1500:]}")
            return False
        try:
            with urllib.request.urlopen(adres, timeout=2) as y:
                if y.status == 200:
                    gunlukle("sunucu hazir")
                    return True
        except (urllib.error.URLError, OSError):
            pass
        time.sleep(0.4)
    gunlukle(f"sunucu {BASLATMA_SURESI} saniyede hazir olmadi")
    return False


def sureci_durdur(surec: subprocess.Popen) -> None:
    """Streamlit alt süreçler açar; ağacın tamamı kapatılmalı."""
    if surec.poll() is not None:
        return
    gunlukle("sunucu kapatiliyor")
    if sys.platform == "win32":
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(surec.pid)],
                capture_output=True, timeout=15,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            return
        except Exception as exc:
            gunlukle(f"taskkill basarisiz: {exc}")
    surec.terminate()
    try:
        surec.wait(timeout=10)
    except subprocess.TimeoutExpired:
        surec.kill()


def main() -> int:
    if not UYGULAMA.exists():
        hata_kutusu(f"Uygulama dosyasi bulunamadi:\n{UYGULAMA}")
        return 1
    try:
        import webview
    except ImportError:
        hata_kutusu(
            "pywebview kurulu degil. Komut isteminde su satiri calistirin:\n\n"
            "pip install -r requirements.txt")
        return 1

    port = bos_port()
    surec = sunucuyu_baslat(port)

    if not hazir_mi(port, surec):
        sureci_durdur(surec)
        hata_kutusu(
            "VectorLab sunucusu baslatilamadi.\n\n"
            "En sik sebep eksik bagimliliktir; su satiri calistirin:\n"
            "pip install -r requirements.txt")
        return 1

    pencere = webview.create_window(
        BASLIK, f"http://127.0.0.1:{port}",
        width=GENISLIK, height=YUKSEKLIK,
        min_size=(MIN_GENISLIK, MIN_YUKSEKLIK),
        confirm_close=False, text_select=True,
    )

    def kapanista() -> None:
        sureci_durdur(surec)

    pencere.events.closed += kapanista

    # Pencere kapanmadan süreç ölürse (çökme) kullanıcı boş ekranla
    # kalmasın; izleyip pencereyi kapatıyoruz.
    def izle() -> None:
        while surec.poll() is None:
            time.sleep(1.0)
        gunlukle("sunucu durdu, pencere kapatiliyor")
        try:
            pencere.destroy()
        except Exception:
            pass

    threading.Thread(target=izle, daemon=True).start()

    gunlukle("pencere aciliyor")
    simge = KOK / "assets" / "vectorlab.ico"
    try:
        if simge.exists():
            webview.start(icon=str(simge))
        else:
            webview.start()
    finally:
        sureci_durdur(surec)
    gunlukle("cikildi")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
