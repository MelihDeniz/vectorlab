"""f-string'leri çevrilebilir hale getirir: f"..{x}.." -> t("..{0}..").format(x)

`tools/i18n_wrap.py` düz dizeleri `t(...)` içine alıyor ama f-string'e
dokunamıyor: f-string çalışma anında değerlendirilir, yani `t()` onu
görmeden birleşmiş bir metin olur ve sözlükte anahtarı bulunmaz. Bu
yüzden İngilizce arayüzde 158 metin Türkçe kalmıştı.

Çözüm, f-string'i bir **şablona** çevirmektir: değişken kısımlar
`{0}`, `{1}` olur, şablon `t()` içine girer ve sonuç `.format()` ile
doldurulur. Şablon artık sabit bir dizedir, yani sözlükte anahtar
olabilir.

Neden `ast.unparse` kullanılmıyor: yorumları siler. Bu depodaki
gerekçe yorumları kodun yarısı kadar ve onları kaybetmek kabul
edilemez. Bu yüzden değişiklik **kaynak metnin kendisinde**, bayt
aralığı değiştirilerek yapılır.

Neden BAYT: `col_offset` bir UTF-8 bayt sapmasıdır, karakter sapması
değil. Türkçe harfler iki bayt tuttuğu için karakter sayarak
hesaplamak dosyayı bozar - bu hata bir kez yapıldı, bu yüzden burada
her şey `bytes` üzerinde yürür.

Hangi f-string'e dokunulur: `ui_metni()` karar verir ve kasıtlı olarak
İHTİYATLIDIR. Dosya adı, FASTA başlığı, markdown tablo satırı ve
sabit metni kısa olan f-string'ler atlanır; bunlar çeviri gerektirmez
ve yanlış dokunmak çıktıyı bozar. Atlanan her f-string Türkçe kalır
ve bu, bozulmasından iyidir.

Kullanım:
    python tools/i18n_fstring.py            # kuru çalıştırma, rapor
    python tools/i18n_fstring.py --uygula   # değişikliği yaz
"""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HEDEF = ROOT / "app_pages"

#: Bu anahtar sözcüklere verilen f-string çevrilmez: kullanıcıya
#: gösterilen metin değil, makineye verilen değerdir.
ATLANAN_KW = {"file_name", "key", "mime", "language", "url", "name"}

TURKCE = set("çğıöşüÇĞİÖŞÜ")


def satir_baslari(ham: bytes) -> list[int]:
    """Her satırın dosyadaki bayt sapması."""
    out, n = [0], 0
    for satir in ham.splitlines(keepends=True):
        n += len(satir)
        out.append(n)
    return out


def aralik(baslar: list[int], node: ast.AST) -> tuple[int, int]:
    return (baslar[node.lineno - 1] + node.col_offset,
            baslar[node.end_lineno - 1] + node.end_col_offset)


def kaynak(ham: bytes, baslar: list[int], node: ast.AST) -> str:
    a, b = aralik(baslar, node)
    return ham[a:b].decode("utf-8")


#: Türkçeye özgü harf içermeyen ama Türkçe olan kısa parçalar. Birim
#: ve simge (bp, %, kDa, °C, nt, kcal/mol) İngilizcede de aynı
#: olduğu için listede YOK: onlara dokunmak gereksiz anahtar üretir.
TURKCE_SOZ = {"okundu", "genotip", "kodon", "adet", "kesim", "tane",
              "bulundu", "secildi", "kayit", "satir", "sonra", "once"}


def ui_metni(literal: str) -> bool:
    """Bu f-string'in sabit kısmı kullanıcıya gösterilen metin mi?

    Üç ölçüt var ve sırası önemli:

    1. Türkçeye özgü bir harf (ç, ğ, ı, ö, ş, ü) varsa metin kesin
       arayüz metnidir - kısa olsa bile. `**Promotör:** {0}` ya da
       `güç: {0}` böyle yakalanır; bunlar boşluk içermediği ya da
       kısa olduğu için ilk sürümde atlanmış ve ekranda Türkçe
       kalmıştı.
    2. Türkçe ama yalnızca ASCII harf kullanan kısa parçalar
       (`{0} nt okundu.`) bir sözcük listesiyle yakalanır.
    3. Geri kalanda ihtiyatlı davranılır: boşluk ve en az 10 harf
       istenir, böylece dosya adı, markdown tablo satırı ve kimlik
       dizeleri atlanır.

    Birimler (` bp`, `%`, ` kDa`, ` °C`) bilerek atlanır: İngilizcede
    de aynı yazılır, çevirmek yalnızca sözlüğü şişirir.
    """
    s = literal.strip()
    if not s or s.startswith(">"):          # FASTA başlığı
        return False
    harf = sum(ch.isalpha() for ch in s)
    if harf < 3:
        return False
    if set(s) & TURKCE:
        return True
    kucuk = s.lower()
    if any(w in kucuk for w in TURKCE_SOZ):
        return True
    if " " not in s:                        # dosya adı, kimlik, kod
        return False
    return harf >= 10


def sablon(node: ast.JoinedStr, ham: bytes, baslar: list[int],
           ) -> tuple[str, list[str], str]:
    """JoinedStr'den (şablon, argüman kaynakları, sabit metin) üret.

    Dönüşüm (`!r`) ve biçim belirteci (`:.1f`) korunur; kaybolursa
    çıktı değişir. Süslü parantezler `.format()` için yeniden
    kaçırılır: f-string'de `{{` yazılmış bir parantez AST'de tek
    parantez olarak gelir ve şablona tek konursa `.format()` onu yer
    tutucu sanar.
    """
    parcalar: list[str] = []
    args: list[str] = []
    sabit: list[str] = []
    for v in node.values:
        if isinstance(v, ast.Constant) and isinstance(v.value, str):
            sabit.append(v.value)
            parcalar.append(v.value.replace("{", "{{").replace("}", "}}"))
            continue
        if not isinstance(v, ast.FormattedValue):
            raise ValueError("beklenmeyen düğüm")
        i = len(args)
        args.append(kaynak(ham, baslar, v.value))
        yt = f"{{{i}"
        if v.conversion and v.conversion != -1:
            yt += "!" + chr(v.conversion)
        if v.format_spec is not None:
            spec = "".join(
                c.value for c in v.format_spec.values
                if isinstance(c, ast.Constant) and isinstance(c.value, str))
            # İçinde kendi yer tutucusu olan biçim belirteci (`{:{w}}`)
            # desteklenmez: nadir ve yanlış dönüştürmek kolay.
            if len(v.format_spec.values) != len(
                    [c for c in v.format_spec.values
                     if isinstance(c, ast.Constant)]):
                raise ValueError("iç içe biçim belirteci")
            yt += ":" + spec
        parcalar.append(yt + "}")
    return "".join(parcalar), args, "".join(sabit)


def adaylar(yol: Path) -> list[tuple[ast.JoinedStr, str, list[str], str]]:
    ham = yol.read_bytes()
    baslar = satir_baslari(ham)
    agac = ast.parse(ham)

    # `t(...)` içindeki her düğümü ve atlanacak anahtar sözcük
    # değerlerini işaretle.
    yasak: set[int] = set()
    for n in ast.walk(agac):
        if isinstance(n, ast.Call):
            if getattr(n.func, "id", None) == "t":
                for c in ast.walk(n):
                    yasak.add(id(c))
            for kw in n.keywords:
                if kw.arg in ATLANAN_KW:
                    for c in ast.walk(kw.value):
                        yasak.add(id(c))

    out = []
    for n in ast.walk(agac):
        if not isinstance(n, ast.JoinedStr) or id(n) in yasak:
            continue
        try:
            sab, args, metin = sablon(n, ham, baslar)
        except ValueError:
            continue
        if not ui_metni(metin):
            continue
        # Çok satırlı argüman kaynağı tek satıra sığmaz; atla.
        if any("\n" in a for a in args):
            continue
        out.append((n, sab, args, metin))
    return out


def yeni_kaynak(sab: str, args: list[str]) -> str:
    gosterim = json.dumps(sab, ensure_ascii=False)
    if not args:
        return f"t({gosterim})"
    return f"t({gosterim}).format({', '.join(args)})"


def dosyayi_isle(yol: Path, uygula: bool) -> list[tuple[str, str]]:
    ham = yol.read_bytes()
    baslar = satir_baslari(ham)
    bulunan = adaylar(yol)
    if not bulunan:
        return []

    # Sondan başa doğru değiştir: baştan değiştirilse sonraki
    # sapmalar kayar.
    degisiklikler = []
    parcalar = sorted(((aralik(baslar, n), sab, args)
                       for n, sab, args, _ in bulunan),
                      key=lambda x: x[0][0], reverse=True)
    yeni = ham
    for (a, b), sab, args in parcalar:
        eski = ham[a:b].decode("utf-8")
        yerine = yeni_kaynak(sab, args)
        degisiklikler.append((eski, yerine))
        yeni = yeni[:a] + yerine.encode("utf-8") + yeni[b:]

    # Yazmadan önce ayrıştırılabilirliği DENETLE. Bozuk Python
    # yazmaktansa hiç yazmamak iyidir.
    try:
        ast.parse(yeni)
    except SyntaxError as exc:
        raise SystemExit(f"{yol.name}: sonuç ayrıştırılamadı ({exc}); "
                         f"dosya değiştirilmedi")
    if uygula:
        yol.write_bytes(yeni)
    return degisiklikler


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--uygula", action="store_true",
                    help="değişikliği dosyaya yaz (yoksa kuru çalıştırma)")
    ap.add_argument("--goster", type=int, default=6,
                    help="dosya başına kaç örnek gösterilsin")
    a = ap.parse_args()

    toplam = 0
    for yol in sorted(HEDEF.glob("*.py")):
        d = dosyayi_isle(yol, a.uygula)
        if not d:
            continue
        toplam += len(d)
        print(f"\n{yol.name}: {len(d)} f-string")
        for eski, yeni in d[:a.goster]:
            print(f"  - {eski[:96]}")
            print(f"  + {yeni[:96]}")
    print(f"\n{'Yazildi' if a.uygula else 'Kuru calistirma'}: "
          f"{toplam} f-string")
    if not a.uygula:
        print("Uygulamak icin: python tools/i18n_fstring.py --uygula")
        print("Sonra: python tools/i18n_scan.py --yaz  (yeni anahtarlar)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
