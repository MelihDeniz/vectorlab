"""Çeviri sözlüğünü kodla eşitler ve kapsamı bildirir.

`t("...")` çağrılarındaki Türkçe metinleri toplar, sözlükle
karşılaştırır ve üç şeyi söyler:

  eksik       kodda var, sözlükte yok  -> çevrilmeyi bekliyor
  kullanılmıyor  sözlükte var, kodda yok  -> Türkçe metin değişmiş
  kapsam      yüzde kaçı çevrilmiş

İkincisi, kaynak metni anahtar yapmanın bedelidir: Türkçe bir etiket
değişirse çeviri bağı sessizce kopar. Bu tarayıcı onu sessiz
bırakmıyor.

    python tools/i18n_scan.py              # rapor
    python tools/i18n_scan.py --yaz        # eksikleri bos deger olarak ekle
    python tools/i18n_scan.py --temizle    # kullanilmayanlari sil
"""

from __future__ import annotations

import argparse
import ast
import json
import pathlib
import sys

KOK = pathlib.Path(__file__).resolve().parent.parent
SAYFALAR = KOK / "app_pages"
I18N = KOK / "data" / "i18n"


def taranan_dosyalar() -> list[pathlib.Path]:
    """Sayfalar ve uygulama girişi.

    `streamlit_app.py` de taranmalı: gezinme başlıkları orada ve onlar
    da çevriliyor. Yalnızca `app_pages/` taransaydı başlıklar sözlüğe
    hiç girmez, çevrilmediği de görülmezdi.
    """
    return sorted(SAYFALAR.glob("*.py")) + [KOK / "streamlit_app.py"]


def kaynak_metinler() -> list[str]:
    """Taranan dosyalardaki `t("...")` argümanları, sırayla, tekrarsız."""
    out: list[str] = []
    gorulen: set[str] = set()
    for yol in taranan_dosyalar():
        agac = ast.parse(yol.read_text(encoding="utf-8"))
        for node in ast.walk(agac):
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "t"
                    and node.args
                    and isinstance(node.args[0], ast.Constant)
                    and isinstance(node.args[0].value, str)):
                s = node.args[0].value
                if s not in gorulen:
                    gorulen.add(s)
                    out.append(s)
    return out


def f_string_sayisi() -> int:
    """Elle çevrilmesi gereken f-string sayısı (sarılamayanlar)."""
    n = 0
    for yol in sorted(SAYFALAR.glob("*.py")):
        agac = ast.parse(yol.read_text(encoding="utf-8"))
        for node in ast.walk(agac):
            if not isinstance(node, ast.Call):
                continue
            ad = (node.func.attr if isinstance(node.func, ast.Attribute)
                  else None)
            if ad in {"markdown", "caption", "write", "info", "warning",
                      "error", "success", "metric", "subheader", "button"}:
                for a in node.args:
                    if isinstance(a, ast.JoinedStr):
                        n += 1
    return n


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dil", default="en")
    ap.add_argument("--yaz", action="store_true",
                    help="eksik kayitlari bos deger olarak sozluge ekle")
    ap.add_argument("--temizle", action="store_true",
                    help="kodda bulunmayan kayitlari sozlukten sil")
    ns = ap.parse_args(argv)

    metinler = kaynak_metinler()
    yol = I18N / f"{ns.dil}.json"
    sozluk: dict[str, str] = {}
    if yol.exists():
        sozluk = json.loads(yol.read_text(encoding="utf-8"))

    kodda = set(metinler)
    eksik = [m for m in metinler if not sozluk.get(m, "").strip()]
    kullanilmiyor = sorted(k for k in sozluk if k not in kodda)
    cevrilen = len(metinler) - len(eksik)

    print(f"Dil            : {ns.dil}")
    print(f"Koddaki metin  : {len(metinler)}")
    print(f"Cevrilen       : {cevrilen}")
    print(f"Eksik          : {len(eksik)}")
    print(f"Kullanilmiyor  : {len(kullanilmiyor)}")
    oran = 100 * cevrilen / max(len(metinler), 1)
    print(f"Kapsam         : %{oran:.1f}")
    print(f"Elle cevrilecek f-string: {f_string_sayisi()}")

    if eksik[:10]:
        print("\nIlk eksikler:")
        for m in eksik[:10]:
            print(f"  - {m[:70]!r}")
    if kullanilmiyor[:10]:
        print("\nKodda artik bulunmayanlar (Turkce metin degismis olabilir):")
        for m in kullanilmiyor[:10]:
            print(f"  - {m[:70]!r}")

    degisti = False
    if ns.yaz:
        for m in metinler:
            sozluk.setdefault(m, "")
        degisti = True
    if ns.temizle:
        for m in kullanilmiyor:
            sozluk.pop(m, None)
        degisti = True
    if degisti:
        # Kod sirasiyla yaz: sozluk okunurken sayfa sayfa ilerlesin.
        sirali = {m: sozluk.get(m, "") for m in metinler}
        for k, v in sozluk.items():
            sirali.setdefault(k, v)
        I18N.mkdir(parents=True, exist_ok=True)
        yol.write_text(json.dumps(sirali, ensure_ascii=False, indent=1),
                       encoding="utf-8")
        print(f"\n{len(sirali)} kayit yazildi -> {yol}")

    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
