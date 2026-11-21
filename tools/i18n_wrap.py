"""Sayfalardaki arayüz metinlerini `t(...)` ile sarar (tek seferlik).

Bu bir kod dönüştürücüdür. 2000'den fazla çağrı yerini elle sarmak
hem uzun sürer hem de hata yapılır; mekanik olan mekanik yapılmalı.

**Neden izin listesi var.** Çağrıdaki her dizeyi sarmak yanlış olur:
`"utf-8".decode(...)`, `dict.get("ad")`, `df.index("Kategori")` gibi
çağrılar da dize alır ve onlar arayüz metni değildir. Bu yüzden
yalnızca Streamlit'in metin gösteren çağrıları ve yalnızca bilinen
argüman konumları sarılır.

**Neden f-string'ler sarılmıyor.** `f"{n} parça"` gibi bir metni
sarmak için önce biçim yerlerinin çıkarılması gerekir
(`t("{n} parça").format(n=n)`); bunu mekanik yapmak, değişken adlarını
ve ifadeleri bozma riski taşır. Onlar elle çevrilir; tarayıcı
(`tools/i18n_scan.py`) kaç tane kaldığını söyler.

**Neden `ast.unparse` kullanılmıyor.** Ağacı yeniden yazdırmak en kolay
yoldu ama kabul edilemez: `ast.unparse` YORUMLARI SİLER. Bu depoda
yorumlar kodun yarısı kadar yer tutuyor ve kararların gerekçesi orada
yazılı; onları kaybetmek, uygulamanın anlatma işlevini kaybetmesi
demek. Bunun yerine her dizenin kaynak dosyadaki tam konumu
(`lineno`/`col_offset`) alınıp yalnızca o aralık değiştiriliyor.
Dosyanın geri kalanına - boşluğuna, yorumuna, satır sonuna -
dokunulmuyor. Değişiklikler sondan başa uygulanıyor, yoksa ilk
değişiklik sonraki konumları kaydırırdı.

    python tools/i18n_wrap.py --kontrol     # ne değişecek, yazmadan
    python tools/i18n_wrap.py               # uygula
"""

from __future__ import annotations

import argparse
import ast
import pathlib
import sys

KOK = pathlib.Path(__file__).resolve().parent.parent
SAYFALAR = KOK / "app_pages"

#: Metin gösteren Streamlit çağrıları. `st.` ya da bir kapsayıcı
#: üzerinden (`c.metric(...)`) çağrılabildikleri için yalnızca
#: öznitelik adına bakılır.
UI_CALLS = {
    "button", "caption", "checkbox", "code", "column_config",
    "date_input", "download_button", "error", "expander",
    "file_uploader", "form_submit_button", "header", "info",
    "link_button", "markdown", "metric", "multiselect", "number_input",
    "page_link", "pills", "radio", "segmented_control", "select_slider",
    "selectbox", "slider", "subheader", "success", "tabs", "text_area",
    "text_input", "title", "toggle", "warning", "write", "popover",
    "status", "toast", "badge", "dataframe", "data_editor",
    "NumberColumn", "ProgressColumn", "TextColumn", "LinkColumn",
}

#: Sarılacak anahtar argümanlar. Konumsal argümanlarda yalnızca ilki
#: (etiket) sarılır - `st.metric("Uzunluk", deger)` çağrısında ikinci
#: argüman veridir, metin değil.
KW_ARGS = {"label", "help", "placeholder", "body", "text"}

#: Bunlar arayüz metni değil: biçim adı, MIME türü, dosya adı, renk.
KW_ATLA = {"mime", "file_name", "language", "key", "icon", "type",
           "format", "width", "alt", "delta_color", "label_visibility"}


def cevrilebilir(s: str) -> bool:
    if len(s.strip()) < 2:
        return False
    if not any(c.isalpha() for c in s):
        return False
    if s.startswith(":material/"):
        return False
    # Biçim yeri içeren metinler elle çevrilir.
    if "{" in s or "}" in s:
        return False
    # Saf ASCII ve tek sözcükse büyük ihtimalle teknik bir değer
    # ("utf-8", "stretch", "text/plain", "ifade" gibi kipler hariç
    # tutulamaz, o yüzden yalnızca belirgin teknik kalıplar elenir).
    dusuk = s.lower()
    if dusuk in {"utf-8", "stretch", "text/plain", "text/html", "ltr",
                 "rtl", "auto", "wide", "centered"}:
        return False
    if "/" in s and " " not in s:          # text/plain, application/json
        return False
    return True


def hedefleri_bul(agac: ast.AST) -> list[ast.Constant]:
    """Sarılacak dize düğümleri."""
    out: list[ast.Constant] = []

    def uygun(n: ast.AST) -> bool:
        return (isinstance(n, ast.Constant) and isinstance(n.value, str)
                and cevrilebilir(n.value))

    for node in ast.walk(agac):
        if not isinstance(node, ast.Call):
            continue
        ad = node.func.attr if isinstance(node.func, ast.Attribute) else None
        if ad not in UI_CALLS:
            continue

        # Konumsal: yalnızca ilk argüman (etiket / gövde).
        # `st.metric("Uzunluk", deger)` çağrısında ikincisi veridir.
        if node.args:
            ilk = node.args[0]
            if uygun(ilk):
                out.append(ilk)
            # `st.tabs([...])`, `st.radio(..., [...])`: seçenek
            # etiketleri de arayüz metnidir.
            elif isinstance(ilk, (ast.List, ast.Tuple)):
                out.extend(e for e in ilk.elts if uygun(e))

        for kw in node.keywords:
            if kw.arg in KW_ATLA or kw.arg not in KW_ARGS:
                continue
            if uygun(kw.value):
                out.append(kw.value)
    return out


def _satir_baslari(ham: bytes) -> list[int]:
    """Her satırın dosya içindeki başlangıç konumu, BAYT olarak.

    Bayt olması şart: `col_offset` karakter değil UTF-8 bayt konumudur.
    Türkçe metinde 'ı', 'ş', 'ğ' iki bayt tuttuğu için karakterle
    hesaplamak konumları kaydırır - ilk denemede bütün dosyalar tam
    olarak bu yüzden derlenmedi.
    """
    bas = [0, 0]
    for satir in ham.splitlines(keepends=True):
        bas.append(bas[-1] + len(satir))
    return bas


def islem(yol: pathlib.Path, yaz: bool) -> int:
    kaynak = yol.read_text(encoding="utf-8")
    ham = kaynak.encode("utf-8")
    agac = ast.parse(kaynak)
    hedefler = hedefleri_bul(agac)
    if not hedefler:
        print(f"  {yol.name:16} değişiklik yok")
        return 0

    bas = _satir_baslari(ham)
    # Sondan başa: her değişiklik kendinden öncekilerin konumunu
    # bozmasın.
    araliklar = sorted(
        ((bas[n.lineno] + n.col_offset,
          bas[n.end_lineno] + n.end_col_offset, n.value)
         for n in hedefler),
        key=lambda x: -x[0])

    yeni_ham = ham
    ornekler: list[str] = []
    for a, b, deger in araliklar:
        yeni_ham = yeni_ham[:a] + b"t(" + yeni_ham[a:b] + b")" + yeni_ham[b:]
        if len(ornekler) < 3:
            ornekler.append(deger[:48])
    yeni = yeni_ham.decode("utf-8")

    # `t` içe aktarılmadıysa ekle. Sarmalar hep çağrı yerlerinde, yani
    # içe aktarma bloğundan SONRA olduğu ve satır eklemediği için
    # buradaki satır numaraları hâlâ geçerli.
    if "from core.i18n import t" not in yeni:
        son_import = max(
            (n.end_lineno for n in agac.body
             if isinstance(n, (ast.Import, ast.ImportFrom))), default=0)
        satirlar = yeni.splitlines(keepends=True)
        satirlar.insert(son_import, "from core.i18n import t\n")
        yeni = "".join(satirlar)

    # Yazmadan önce sonucun hâlâ geçerli Python olduğunu doğrula.
    try:
        ast.parse(yeni)
    except SyntaxError as e:
        print(f"  {yol.name:16} ATLANDI - sonuç derlenmiyor: {e}")
        return 0

    print(f"  {yol.name:16} {len(hedefler):>4} metin  "
          f"(ör. {', '.join(repr(o) for o in reversed(ornekler))})")
    if yaz:
        yol.write_text(yeni, encoding="utf-8")
    return len(hedefler)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kontrol", action="store_true",
                    help="dosyaları değiştirme, yalnızca say")
    ns = ap.parse_args(argv)
    toplam = 0
    for yol in sorted(SAYFALAR.glob("*.py")):
        toplam += islem(yol, yaz=not ns.kontrol)
    print(f"\n{toplam} metin sarıldı"
          f"{' (yazılmadı, --kontrol)' if ns.kontrol else ''}.")
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
