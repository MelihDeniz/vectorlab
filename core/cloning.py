"""Restriksiyonla klonlama planı: vektörü kes, inserti tak, kontrol et.

Burada hesaplanan şeyler deneyde gerçekten başarısızlığa yol açan
noktalar: enzim vektörü birden fazla yerden kesiyor mu, insertin içinde
kesiyor mu, iki uç birbirine uyumlu mu (uyumluysa insert ters de girer),
vektör kendi üzerine kapanabilir mi.

Modül Streamlit'ten bağımsızdır, bu yüzden ayrıca sınanabilir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core import enzymes as ez
from core.seq_utils import clean_seq, compatible_ends, restriction_sites


@dataclass
class Note:
    severity: str   # "hata" | "uyari" | "bilgi"
    title: str
    detail: str


@dataclass
class Plan:
    vector_name: str
    vector_length: int
    enzymes: list[str]
    cuts: dict[str, list[int]] = field(default_factory=dict)
    insert_cuts: dict[str, list[int]] = field(default_factory=dict)
    backbone_length: int = 0      # insertin takılacağı iskelet
    removed_length: int = 0       # kesilip atılan parça
    product_length: int = 0       # beklenen son yapı
    notes: list[Note] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(n.severity == "hata" for n in self.notes)


def _span(a: int, b: int, length: int) -> int:
    """a'dan b'ye dairesel mesafe (a dahil değil, b dahil değil)."""
    return (b - a) % length or length


def plan(
    vector_seq: str,
    vector_name: str,
    enzyme_names: list[str],
    insert_seq: str = "",
    insert_name: str = "insert",
    circular: bool = True,
) -> Plan:
    """İki (ya da bir) enzimle klonlama planını hesapla ve denetle."""
    v = clean_seq(vector_seq)
    ins = clean_seq(insert_seq)
    names = [n for n in dict.fromkeys(enzyme_names) if n]
    p = Plan(vector_name=vector_name, vector_length=len(v), enzymes=names)

    if not v:
        p.notes.append(Note("hata", "Vektör dizisi yok",
                            "Önce bir vektör seçin."))
        return p
    if not names:
        p.notes.append(Note("bilgi", "Enzim seçilmedi",
                            "Kesim için en az bir enzim seçin."))
        return p

    p.cuts = {n: restriction_sites(v, [n], circular).get(n, []) for n in names}
    if ins:
        p.insert_cuts = {
            n: restriction_sites(ins, [n], circular=False).get(n, [])
            for n in names
        }

    # --- Vektörde kaç kesim var
    for n, pos in p.cuts.items():
        if not pos:
            p.notes.append(Note(
                "hata", f"{n} vektörü kesmiyor",
                f"{vector_name} içinde {n} tanıma dizisi yok, bu enzimle "
                f"vektör açılamaz. Dizide kesen bir enzim seçin."))
        elif len(pos) > 1:
            p.notes.append(Note(
                "uyari", f"{n} vektörü {len(pos)} yerden kesiyor",
                f"Kesim konumları: {', '.join(str(x) for x in pos)}. "
                f"Tek kesen bir enzim seçmezseniz iskelet parçalanır ve "
                f"doğru fragmanı jelden ayırmanız gerekir."))

    # --- İnsertin içinde kesiyor mu
    for n, pos in p.insert_cuts.items():
        if pos:
            p.notes.append(Note(
                "hata", f"{n} insertin içinde kesiyor",
                f"{insert_name} dizisinde {len(pos)} kesim var "
                f"({', '.join(str(x) for x in pos)}). Bu enzimle "
                f"keserseniz insert parçalanır. Başka bir enzim seçin ya da "
                f"inserti sessiz mutasyonla bu bölgeden arındırın."))

    # --- Uç geometrisi
    geoms = {n: ez.geometry(n) for n in names}
    for n, g in geoms.items():
        if g is None:
            continue
        if g.overhang == 0:
            p.notes.append(Note(
                "uyari", f"{n} kör uç bırakıyor",
                "Kör uç ligasyonunda yön kontrolü yoktur: insert iki yönde "
                "de girer ve verim yapışkan uca göre düşüktür. Mümkünse "
                "yapışkan uçlu bir enzim seçin."))

    if len(names) == 1:
        p.notes.append(Note(
            "uyari", "Tek enzimle kesim yönsüzdür",
            "Vektörün iki ucu aynı olduğu için insert iki yönde de girer ve "
            "vektör kendi üzerine kapanabilir (boş koloni). İki farklı "
            "enzimle yönlü klonlama yapın; tek enzimle çalışacaksanız "
            "kesilmiş vektörü fosfatazla muamele edin."))
    elif len(names) >= 2:
        a, b = names[0], names[1]
        if compatible_ends(a, b):
            p.notes.append(Note(
                "uyari", f"{a} ve {b} uyumlu uç bırakıyor",
                "İki uç birbirine ligate olabildiği için insert ters de "
                "girer ve vektör kendi üzerine kapanır. Yönlü klonlama "
                "istiyorsanız uyumsuz uç bırakan bir çift seçin - Enzimler "
                "sayfasındaki Uyumlu uçlar sekmesi bu çiftleri listeler."))

    # --- Metilasyon
    for n in names:
        blocked = ez.METHYLATION_BLOCKED.get(n)
        if blocked:
            p.notes.append(Note(
                "bilgi", f"{n} metilasyona duyarlı",
                f"{blocked}. E. coli'den saflaştırılmış plazmitte bu bölge "
                f"kesilmeyebilir; gerekirse "
                f"{ez.DAM_DCM_FREE_STRAINS[0]} gibi dam-/dcm- bir suşta "
                f"çoğaltın."))

    # --- Boyut hesabı: iki tek kesen enzim varsa iskelet ile atılan parça
    singles = [n for n in names if len(p.cuts.get(n, [])) == 1]
    if len(singles) >= 2:
        c1 = p.cuts[singles[0]][0]
        c2 = p.cuts[singles[1]][0]
        # İki kesim arasındaki iki yay: kısası atılır varsayımı yanlış
        # olabilir, bu yüzden ikisini de hesaplayıp küçüğünü "atılan"
        # sayıyoruz ve bunu kullanıcıya söylüyoruz.
        arc1 = _span(c1, c2, len(v))
        arc2 = len(v) - arc1
        p.removed_length = min(arc1, arc2)
        p.backbone_length = max(arc1, arc2)
        p.product_length = p.backbone_length + len(ins)
        p.notes.append(Note(
            "bilgi", "Beklenen boyutlar",
            f"Kesimler {c1} ve {c2}. İki parça çıkar: "
            f"{p.backbone_length} bp ve {p.removed_length} bp. Büyük olanı "
            f"iskelet kabul edildi; insert ({len(ins)} bp) takıldığında "
            f"son yapı {p.product_length} bp olur. MCS'in iki yakın bölgesi "
            f"arasında kesiyorsanız bu doğrudur; aksi halde hangi parçanın "
            f"iskelet olduğunu haritadan doğrulayın."))
    elif len(singles) == 1 and len(names) == 1:
        p.backbone_length = len(v)
        p.product_length = len(v) + len(ins)
        p.notes.append(Note(
            "bilgi", "Beklenen boyut",
            f"Tek kesimle vektör açılır (hiçbir parça atılmaz): son yapı "
            f"{p.product_length} bp olur."))

    return p


def fragments(vector_seq: str, enzyme_names: list[str],
              circular: bool = True) -> list[dict[str, Any]]:
    """Jel benzeri fragman listesi (büyükten küçüğe)."""
    from core.seq_utils import digest
    return digest(clean_seq(vector_seq), enzyme_names, circular)
