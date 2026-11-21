"""Tasarımdan malzeme listesi: ne sipariş edilecek, ne kadar, kaç günde.

Tasarım bitince iş bitmiyor. Sıradaki soru "ne satın alacağım" ve bu
soru pratikte en çok para ve en çok zaman kaybettiren yer: eksik
sipariş edilen bir oligo bir haftaya, yanlış tampon bir deneye mal
olur. Bu modül tasarımdan doğrudan bir malzeme listesi (BOM) çıkarır.

**Bu modül FİYAT İÇERMEZ ve KATALOG NUMARASI İÇERMEZ.** İkisi de
uydurmaya en açık veri türüdür: fiyat ülkeye, kuruma, sözleşmeye ve
yıla göre değişir; katalog numarası firmaya göre değişir ve biri
yanlış yazılırsa yanlış ürün gelir. Burada hesaplanan şey
**TÜKETİM**'dir - kaç baz oligo, kaç ünite enzim, kaç reaksiyon, kaç
iş günü. Birim fiyatı kullanıcı kendi laboratuvarının fiyat
listesinden girer ve `cost()` çarpımı yapar. Hesabın doğruluğu
bizden, fiyat sizden.

Enzimi hangi firmaların sattığı **REBASE'den** okunur (Biopython'un
`supplier_list()` çağrısıyla), ezberden yazılmaz. Firma adı
listelenmemiş bir enzim için liste boş döner ve bu gizlenmez.

Süre tahminleri ise ne kaynaktan ne hesaptan gelir: laboratuvarda
yerleşik iş akışının alışıldık süreleridir. Bu yüzden hepsi
`PHASE_DAYS` içinde tek yerde durur, `build()` çağrısında
değiştirilebilir ve çıktıda "varsayım" diye etiketlenir. Dizileme
süresi özellikle hizmet sağlayıcıya bağlıdır; varsayılanı kendi
sağlayıcınızla değiştirin.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass, field
from typing import Any

from core import protocol
from core.seq_utils import Primer

# ---------------------------------------------------------------------------
# Varsayımlar - hepsi tek yerde, hepsi değiştirilebilir
# ---------------------------------------------------------------------------

#: Oligo sipariş ölçeği: (üst sınır nt, ölçek, saflaştırma, gerekçe).
#: Sıra önemli: ilk uyan kullanılır. Klonlama primerleri için standart
#: tuz giderme yeterlidir; uzun oligoda hatalı sentez oranı
#: uzunlukla arttığı için saflaştırma önerilir.
OLIGO_SCALE: list[tuple[int, str, str, str]] = [
    (60, "25 nmol", "Standart (tuz giderme)",
     "Klonlama ve koloni PCR primerleri için yeterli; bir sipariş "
     "yüzlerce reaksiyona yeter."),
    (90, "100 nmol", "Standart (tuz giderme)",
     "60 nt üstünde verim düşer, daha büyük ölçek istemek gerekir. "
     "Uzun kuyruklu Gibson primerleri burada çıkar."),
    (200, "100 nmol", "HPLC ya da PAGE",
     "90 nt üstünde hatalı sentez (kısa ürün) oranı ciddileşir; "
     "saflaştırma olmadan karışımın bir kısmı eksik oligodur."),
]

#: 200 nt üstü için oligo değil sentez önerilir.
SYNTHESIS_THRESHOLD = 200

#: Bir sindirim reaksiyonunda kullanılan enzim ünitesi. 1 µg DNA'yı
#: 1 saatte kesmek için 5-10 U alışıldık aralıktır; üst sınır alındı,
#: çünkü eksik kesim en sık klonlama hatasıdır ve enzim pahalı değildir.
UNITS_PER_DIGEST = 10.0

#: Reaksiyon sayısı varsayımları: her adım kaç kez kurulur. Deneme
#: yanılma payı dahildir - ilk seferde tutmayan adım tekrarlanır.
DEFAULT_REACTIONS: dict[str, int] = {
    "pcr": 2,              # insert PCR'ı, bir tekrar payıyla
    "digest": 2,           # vektör ve insert ayrı ayrı
    "ligation": 2,         # bağlanma ve negatif kontrol (vektör-tek)
    "transformation": 2,   # bağlanma ve kontrol
    "colony_pcr": 8,       # bir plakta sekiz koloni taranır
    "miniprep": 2,         # doğru çıkan iki koloni
    "sequencing": 2,       # iki yönden okuma
}

#: Aşama süreleri (iş günü). Kaynaktan gelmez; yerleşik iş akışının
#: alışıldık süreleridir ve `build(phase_days=...)` ile değiştirilir.
PHASE_DAYS: dict[str, float] = {
    "oligo_teslim": 2.0,
    "pcr": 0.5,
    "digest_ligation": 1.0,
    "transformation": 1.0,     # dönüşüm + bir gece inkübasyon
    "colony_pcr": 0.5,
    "miniprep": 0.5,
    "sequencing": 3.0,         # hizmet sağlayıcıya bağlı
    "expression": 2.0,
}


# ---------------------------------------------------------------------------
# Veri yapıları
# ---------------------------------------------------------------------------


@dataclass
class Item:
    """Sipariş listesindeki tek satır.

    `quantity` ve `unit` birlikte okunur: 2 "reaksiyon", 20 "ünite",
    46 "baz". Fiyat burada YOK; `cost()` kullanıcının birim fiyatıyla
    çarpar.
    """

    kind: str                  # oligo | enzim | sus | kit | kimyasal
    name: str
    quantity: float
    unit: str
    note: str = ""
    suppliers: list[str] = field(default_factory=list)
    sequence: str = ""
    scale: str = ""
    purification: str = ""

    @property
    def supplier_text(self) -> str:
        """Satan firmalar - yalnızca enzimler için anlamlı.

        REBASE bir restriksiyon enzimi veritabanıdır; bir miniprep
        kiti ya da yetkin hücre orada hiç bulunmaz. Her kalemde
        "REBASE'de listelenmemiş" yazmak, o kalemin aranıp
        bulunamadığı izlenimi verir - oysa aranmadı bile. Bu yüzden
        enzim dışındaki kalemde sütun boş kalır.
        """
        if self.suppliers:
            return ", ".join(self.suppliers)
        if self.kind == "enzim":
            return "REBASE'de firma listelenmemiş"
        return ""


@dataclass
class Phase:
    name: str
    days: float
    note: str = ""


@dataclass
class Bom:
    """Malzeme listesi (bill of materials) ve takvim."""

    title: str
    items: list[Item] = field(default_factory=list)
    phases: list[Phase] = field(default_factory=list)
    assumptions: list[str] = field(default_factory=list)

    @property
    def total_days(self) -> float:
        return sum(p.days for p in self.phases)

    def of_kind(self, kind: str) -> list[Item]:
        return [i for i in self.items if i.kind == kind]

    @property
    def oligo_bases(self) -> int:
        """Toplam sipariş edilecek baz sayısı.

        Oligo fiyatlandırması baz başınadır, bu yüzden toplam baz
        tek başına anlamlı bir sayıdır.
        """
        return int(sum(i.quantity for i in self.of_kind("oligo")))


# ---------------------------------------------------------------------------
# Oligolar
# ---------------------------------------------------------------------------


def oligo_scale(length: int) -> tuple[str, str, str]:
    """Uzunluğa göre (ölçek, saflaştırma, gerekçe) döndür."""
    for ust, olcek, saf, gerekce in OLIGO_SCALE:
        if length <= ust:
            return olcek, saf, gerekce
    return ("sentez", "—",
            f"{SYNTHESIS_THRESHOLD} nt üstünde oligo sentezi güvenilir "
            "değildir; bu uzunlukta parçayı oligo olarak değil "
            "**gen sentezi** olarak sipariş edin.")


def as_primer(x: Any, fallback_name: str = "primer") -> Primer:
    """Primer benzeri her şeyi `Primer`a çevir.

    Bu modülün girdisi tek bir yerden gelmiyor: `primers.colony_pcr`
    gerçek `Primer` nesnesi döner, `primers.catalog_anchors` ise sözlük
    içinde düz dize döner, kullanıcı da elle bir dizi yapıştırabilir.
    Çeviriyi burada bir kez yapmak, çağıran her yerin aynı dönüşümü
    kendi başına yapmasından iyidir - o yolda biri eksik yapar.

    Tm ve GC verilmemişse HESAPLANIR, boş bırakılmaz: sipariş
    sayfasındaki boş Tm sütunu, primerin Tm'si yok demek değil, bizim
    yazmadığımız demektir ve ikisi karışır.
    """
    if isinstance(x, Primer):
        return x
    from core.seq_utils import clean_seq, gc_percent, tm_nn
    if isinstance(x, dict):
        dizi = clean_seq(str(x.get("primer") or x.get("sequence") or ""))
        ad = str(x.get("name") or x.get("part_id") or fallback_name)
    else:
        dizi = clean_seq(str(x))
        ad = fallback_name
    return Primer(name=ad, sequence=dizi, tm=tm_nn(dizi) if dizi else None,
                  gc=gc_percent(dizi) if dizi else 0.0, length=len(dizi))


def oligo_items(primers: list[Any]) -> list[Item]:
    """Primer listesini sipariş satırlarına çevir.

    Aynı dizi iki kez geçerse (örneğin koloni PCR'ı ile dizileme
    çapası aynı primeri kullanıyorsa) TEK satır yazılır: iki kez
    sipariş etmek çıplak gözle görülmeyen ama parayla ölçülen bir
    hatadır. Birleştirilen satırın notunda hangi adlarla geçtiği
    yazılır.
    """
    bulunan: dict[str, Item] = {}
    adlar: dict[str, list[str]] = {}
    for sira, ham in enumerate(primers, start=1):
        p = as_primer(ham, f"oligo-{sira}")
        if not p.sequence:
            continue
        dizi = p.sequence.upper()
        adlar.setdefault(dizi, []).append(p.name)
        if dizi in bulunan:
            continue
        olcek, saf, gerekce = oligo_scale(len(dizi))
        bulunan[dizi] = Item(
            kind="oligo", name=p.name, quantity=float(len(dizi)),
            unit="baz", sequence=dizi, scale=olcek, purification=saf,
            note=gerekce)
    for dizi, item in bulunan.items():
        if len(adlar[dizi]) > 1:
            item.note = (f"Aynı dizi {len(adlar[dizi])} ayrı adla geçiyor "
                         f"({', '.join(adlar[dizi])}); bir kez sipariş "
                         f"edilir. " + item.note)
    return list(bulunan.values())


def oligo_order_csv(primers: list[Any]) -> str:
    """Oligo sipariş sayfası (CSV).

    Sütun adları oligo sağlayıcılarının toplu yükleme sayfalarında
    yerleşik olan dörtlüdür (ad, dizi, ölçek, saflaştırma). Sağlayıcı
    farklı bir sütun düzeni isterse bu dosya yine okunur, çünkü düz
    CSV'dir ve başlık satırı açıktır.

    Dizi BÜYÜK HARF yazılır ve kuyruk ile çakışma bölgesi ayrı
    sütunlarda da verilir: bir Gibson primerinde hangi kısmın şablona
    oturduğunu sipariş sayfasından okumak, sonradan sorun ararken
    işe yarar.
    """
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    w.writerow(["Name", "Sequence", "Scale", "Purification",
                "Length", "Tail", "Anneal", "Tm"])
    yazilan: set[str] = set()
    for sira, ham in enumerate(primers, start=1):
        p = as_primer(ham, f"oligo-{sira}")
        dizi = p.sequence.upper()
        if not dizi or dizi in yazilan:
            continue
        yazilan.add(dizi)
        olcek, saf, _ = oligo_scale(len(dizi))
        w.writerow([p.name, dizi, olcek, saf, len(dizi),
                    p.tail.upper(), p.anneal.upper(),
                    "" if p.tm is None else f"{p.tm:.1f}"])
    return out.getvalue()


# ---------------------------------------------------------------------------
# Enzimler
# ---------------------------------------------------------------------------


def enzyme_suppliers(name: str) -> list[str]:
    """Enzimi satan firmalar - REBASE'den, Biopython üzerinden.

    Bulunamayan ya da firma listesi boş olan enzim için boş liste
    döner. Boş listeyi "firma yok" diye değil "REBASE'de
    listelenmemiş" diye okumak gerekir; `Item.supplier_text` bunu
    böyle yazar.
    """
    try:
        from Bio import Restriction
    except ImportError:                                 # pragma: no cover
        return []
    enz = getattr(Restriction, name, None)
    if enz is None:
        return []
    try:
        return list(enz.supplier_list())
    except Exception:                                   # pragma: no cover
        return []


def enzyme_items(enzymes: list[str], reactions: int | None = None,
                 units_per_digest: float = UNITS_PER_DIGEST) -> list[Item]:
    """Kesim enzimlerini sipariş satırlarına çevir."""
    n = DEFAULT_REACTIONS["digest"] if reactions is None else reactions
    out: list[Item] = []
    for ad in dict.fromkeys(enzymes):                   # sırayı koru, tekille
        out.append(Item(
            kind="enzim", name=ad,
            quantity=n * units_per_digest, unit="ünite",
            suppliers=enzyme_suppliers(ad),
            note=(f"{n} sindirim reaksiyonu x {units_per_digest:g} U. "
                  f"1 µg DNA'yı 1 saatte kesmek için 5-10 U alışıldık "
                  f"aralıktır; üst sınır alındı çünkü eksik kesim en sık "
                  f"klonlama hatasıdır ve enzim bu listede en pahalı "
                  f"kalem değildir.")))
    return out


# ---------------------------------------------------------------------------
# Suş, kit ve kimyasallar
# ---------------------------------------------------------------------------


def strain_items(host_name: str, transformations: int | None = None,
                 ) -> list[Item]:
    n = (DEFAULT_REACTIONS["transformation"] if transformations is None
         else transformations)
    if not host_name:
        return []
    return [Item(
        kind="sus", name=f"{host_name} yetkin hücre", quantity=float(n),
        unit="dönüşüm",
        note=("Biri bağlanma karışımı, biri vektör-tek negatif kontrol "
              "için. Negatif kontrolü atlamak, kesilmemiş vektörden "
              "gelen kolonileri doğru klon sanmaya yol açar."))]


def kit_items(sequencing: int | None = None,
              miniprep: int | None = None) -> list[Item]:
    ns = DEFAULT_REACTIONS["sequencing"] if sequencing is None else sequencing
    nm = DEFAULT_REACTIONS["miniprep"] if miniprep is None else miniprep
    out = [
        Item(kind="kit", name="Plazmit miniprep", quantity=float(nm),
             unit="kolon",
             note="Doğru çıkan koloniden plazmit saflaştırma."),
        Item(kind="kit", name="PCR / jel saflaştırma", quantity=float(nm + 2),
             unit="kolon",
             note="Kesilmiş vektör ve insertin jelden geri alınması."),
        Item(kind="kit", name="Sanger dizileme reaksiyonu",
             quantity=float(ns), unit="okuma",
             note=("İki yönden okuma. Uzun bir kaseti tek okumayla "
                   "kapatamazsınız: Sanger okuması 700-900 bazda "
                   "güvenilirliğini kaybeder, bu yüzden Doğrulama "
                   "sayfası okuma sayısını kaset uzunluğundan hesaplar.")),
    ]
    return out


def chemical_items(part_ids: list[str]) -> list[Item]:
    """Antibiyotik ve indükleyiciyi tasarımın parçalarından çıkar.

    Derişimler `core/protocol.py` içindeki tablodan gelir; iki yerde
    ayrı ayrı yazılsa biri değiştiğinde öbürü sessizce yanlış olurdu.
    """
    out: list[Item] = []
    ab = protocol._antibiotic(part_ids)
    if ab:
        ad, derisim, not_ = ab
        out.append(Item(kind="kimyasal", name=ad, quantity=1.0, unit="şişe",
                        note=f"Çalışma derişimi {derisim}. {not_}"))
    ind = protocol._induction(part_ids)
    if ind:
        ad, derisim, not_ = ind
        out.append(Item(kind="kimyasal", name=ad, quantity=1.0, unit="şişe",
                        note=f"Çalışma derişimi {derisim}. {not_}"))
    return out


# ---------------------------------------------------------------------------
# Takvim
# ---------------------------------------------------------------------------


def timeline(needs_pcr: bool = True, needs_sequencing: bool = True,
             needs_expression: bool = True,
             phase_days: dict[str, float] | None = None) -> list[Phase]:
    """Aşama aşama iş günü.

    Toplam, adımların toplamıdır - paralelleştirme varsayılmaz. Bir
    laboratuvarda iki adım aynı güne sığabilir, ama planlarken üst
    sınırı bilmek işe yarar.
    """
    g = dict(PHASE_DAYS)
    if phase_days:
        g.update(phase_days)
    out = [Phase("Oligo teslimi", g["oligo_teslim"],
                 "Sipariş verildikten sonra beklenen süre; "
                 "sağlayıcınıza göre değişir.")]
    if needs_pcr:
        out.append(Phase("PCR", g["pcr"],
                         "Insert çoğaltma ve jelde doğrulama."))
    out += [
        Phase("Kesim ve bağlanma", g["digest_ligation"],
              "Vektör ve insert sindirimi, jelden geri alma, ligasyon."),
        Phase("Dönüşüm ve plak", g["transformation"],
              "Dönüşüm ve bir gece inkübasyon."),
        Phase("Koloni PCR", g["colony_pcr"],
              "Plaktan koloni tarama."),
        Phase("Miniprep", g["miniprep"], "Doğru kolonilerden plazmit."),
    ]
    if needs_sequencing:
        out.append(Phase("Dizileme", g["sequencing"],
                         "Hizmet sağlayıcıya bağlı; en değişken kalem "
                         "budur."))
    if needs_expression:
        out.append(Phase("İfade denemesi", g["expression"],
                         "İndükleme ve ürünün jelde görülmesi."))
    return out


# ---------------------------------------------------------------------------
# Bütünü kur
# ---------------------------------------------------------------------------


def build(title: str, part_ids: list[str], enzymes: list[str] | None = None,
          primers: list[Any] | None = None, host_name: str = "",
          needs_pcr: bool = True, needs_sequencing: bool = True,
          needs_expression: bool = True,
          reactions: dict[str, int] | None = None,
          phase_days: dict[str, float] | None = None) -> Bom:
    """Tasarımdan malzeme listesi ve takvim üret."""
    r = dict(DEFAULT_REACTIONS)
    if reactions:
        r.update(reactions)

    items: list[Item] = []
    items += oligo_items(primers or [])
    items += enzyme_items(enzymes or [], reactions=r["digest"])
    items += strain_items(host_name, transformations=r["transformation"])
    items += kit_items(sequencing=r["sequencing"], miniprep=r["miniprep"])
    items += chemical_items(part_ids)

    varsayimlar = [
        "**Fiyat yok, katalog numarası yok.** İkisi de ülkeye, kuruma "
        "ve yıla göre değişir; uydurulmuş bir fiyat yanlış bütçeye, "
        "uydurulmuş bir katalog numarası yanlış ürüne yol açar. "
        "Buradaki sayılar TÜKETİM'dir; birim fiyatı kendi fiyat "
        "listenizden girin.",
        f"Enzim ünitesi: reaksiyon başına {UNITS_PER_DIGEST:g} U "
        f"(1 µg DNA / 1 saat için alışıldık aralık 5-10 U, üst sınır "
        f"alındı).",
        f"Reaksiyon sayıları tekrar payı içerir: koloni PCR "
        f"{r['colony_pcr']} koloni, dizileme {r['sequencing']} okuma, "
        f"dönüşüm {r['transformation']} (biri negatif kontrol).",
        "Süreler kaynaktan gelmez: yerleşik iş akışının alışıldık "
        "süreleridir ve toplam, adımların toplamıdır "
        "(paralelleştirme varsayılmaz). Dizileme süresi hizmet "
        "sağlayıcınıza göre değişir.",
        "Enzimi satan firmalar REBASE'den okunur (Biopython "
        "`supplier_list()`); liste boşsa enzim satılmıyor demek "
        "değil, REBASE'de firma listelenmemiş demektir.",
    ]

    return Bom(title=title, items=items,
               phases=timeline(needs_pcr, needs_sequencing,
                               needs_expression, phase_days),
               assumptions=varsayimlar)


# ---------------------------------------------------------------------------
# Fiyat - birim fiyat KULLANICIDAN gelir
# ---------------------------------------------------------------------------


def cost(bom: Bom, unit_prices: dict[str, float],
         currency: str = "TL") -> dict[str, Any]:
    """Birim fiyatlarla toplam hesapla.

    `unit_prices` anahtarı kalemin `kind` değeri ya da tam adıdır; tam
    ad önce denenir, böylece "oligo" için genel bir baz fiyatı verip
    tek bir enzime ayrı fiyat yazmak mümkün olur.

    Fiyatı girilmeyen kalem toplama KATILMAZ ve ayrı listelenir.
    Eksik fiyatı sıfır saymak, bütçeyi sessizce küçük gösterir -
    bütçe hatalarının en sık sebebi budur.
    """
    satirlar: list[dict[str, Any]] = []
    eksik: list[str] = []
    toplam = 0.0
    for i in bom.items:
        birim = unit_prices.get(i.name, unit_prices.get(i.kind))
        if birim is None:
            eksik.append(f"{i.name} ({i.quantity:g} {i.unit})")
            continue
        tutar = i.quantity * float(birim)
        toplam += tutar
        satirlar.append({"kalem": i.name, "tur": i.kind,
                         "miktar": i.quantity, "birim": i.unit,
                         "birim_fiyat": float(birim), "tutar": tutar})
    return {"satirlar": satirlar, "toplam": toplam, "para_birimi": currency,
            "fiyati_girilmeyen": eksik,
            "not": ("Fiyatı girilmeyen kalemler toplama katılmadı; "
                    "toplam bu yüzden ALT SINIR'dır."
                    if eksik else
                    "Bütün kalemlerin birim fiyatı girildi.")}


# ---------------------------------------------------------------------------
# Çıktı
# ---------------------------------------------------------------------------


def to_csv(bom: Bom) -> str:
    """Malzeme listesini CSV olarak ver (fiyat sütunu boş bırakılır).

    Boş `Birim fiyat` sütunu kasıtlıdır: dosyayı tabloda açıp kendi
    fiyatınızı yazmanız içindir. Sütunu doldurulmuş sanmaktansa boş
    görmek iyidir.
    """
    out = io.StringIO()
    w = csv.writer(out, lineterminator="\n")
    w.writerow(["Tur", "Kalem", "Miktar", "Birim", "Olcek", "Saflastirma",
                "Dizi", "Satan firmalar", "Birim fiyat", "Not"])
    for i in bom.items:
        w.writerow([i.kind, i.name, f"{i.quantity:g}", i.unit, i.scale,
                    i.purification, i.sequence, i.supplier_text, "",
                    i.note.replace("\n", " ")])
    return out.getvalue()


def to_markdown(bom: Bom) -> str:
    """İnsan okuyacak sürüm."""
    s = [f"# {bom.title} — malzeme listesi", ""]
    s.append(f"**Toplam süre:** {bom.total_days:g} iş günü "
             f"(adımların toplamı; paralelleştirme varsayılmadı)")
    if bom.oligo_bases:
        s.append(f"**Toplam oligo:** {bom.oligo_bases} baz")
    s.append("")

    TUR_ADI = {"oligo": "Oligolar", "enzim": "Enzimler", "sus": "Suşlar",
               "kit": "Kitler", "kimyasal": "Kimyasallar"}
    for tur, baslik in TUR_ADI.items():
        kalemler = bom.of_kind(tur)
        if not kalemler:
            continue
        s += [f"## {baslik}", ""]
        if tur == "oligo":
            s.append("| Ad | Dizi | nt | Ölçek | Saflaştırma |")
            s.append("|---|---|---|---|---|")
            for i in kalemler:
                s.append(f"| {i.name} | `{i.sequence}` | {i.quantity:g} | "
                         f"{i.scale} | {i.purification} |")
        elif tur == "enzim":
            s.append("| Enzim | Ünite | Satan firmalar (REBASE) |")
            s.append("|---|---|---|")
            for i in kalemler:
                s.append(f"| {i.name} | {i.quantity:g} | {i.supplier_text} |")
        else:
            s.append("| Kalem | Miktar | Not |")
            s.append("|---|---|---|")
            for i in kalemler:
                s.append(f"| {i.name} | {i.quantity:g} {i.unit} | {i.note} |")
        s.append("")

    s += ["## Takvim", "", "| Aşama | İş günü | Not |", "|---|---|---|"]
    for p in bom.phases:
        s.append(f"| {p.name} | {p.days:g} | {p.note} |")
    s.append(f"| **Toplam** | **{bom.total_days:g}** | |")
    s += ["", "## Varsayımlar", ""]
    for a in bom.assumptions:
        s.append(f"* {a}")
    return "\n".join(s) + "\n"
