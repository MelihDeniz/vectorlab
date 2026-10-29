"""Tıkandığında ne yapacağın: çift kesim uyumluluğu ve sorun giderme.

İki ayrı iş yapar.

**Çift kesim denetimi.** İki enzimi aynı tüpte kullanmak her zaman
mümkün değildir. Burada doğrulanabilir olanlar kontrol edilir: optimum
sıcaklıklar tutuyor mu, biri metilasyona duyarlı mı, kesim sonrası uçlar
birbirine uyumlu mu. Tampon seçimi bilinçli olarak DIŞARIDA bırakıldı -
hangi enzimin hangi tamponda yüzde kaç çalıştığı üreticiye özgü bir
veridir ve uydurulamaz; kod bunu uydurmak yerine nereye bakılacağını
söyler.

**Sorun giderme.** Belirtiden yola çıkıp olası sebepleri sıralar. Bir
kısmı yalnızca metin, ama bir kısmı tasarımın kendisine bakarak
"sizin yapınızda bu sorun VAR" diyebilir - asıl değeri burada.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core import enzymes as ez

SUPPLIER_NOTE = (
    "Tampon seçimini üreticinin kendi uyumluluk aracından doğrulayın "
    "(NEB Double Digest Finder, Thermo DoubleDigest). Hangi enzimin "
    "hangi tamponda yüzde kaç aktif olduğu üreticiye özgü ölçüm "
    "verisidir; bu uygulama onu tahmin etmez."
)


@dataclass
class Check:
    kind: str          # "hata" | "uyari" | "bilgi"
    title: str
    detail: str


@dataclass
class Cause:
    title: str
    detail: str
    check: str = ""            # kullanıcının bakacağı şey
    applies: bool | None = None   # tasarımdan otomatik belirlenebildiyse


@dataclass
class Symptom:
    key: str
    label: str
    causes: list[Cause] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Çift kesim
# ---------------------------------------------------------------------------


def double_digest(enzyme_names: list[str]) -> list[Check]:
    """İki (ya da daha fazla) enzimi aynı tüpte kullanmak uygun mu."""
    adlar = [n for n in dict.fromkeys(enzyme_names) if n]
    out: list[Check] = []
    if len(adlar) < 2:
        return [Check("bilgi", "Tek enzim",
                      "Çift kesim denetimi için en az iki enzim seçin.")]

    detaylar = {n: ez.details(n) for n in adlar}
    eksik = [n for n, d in detaylar.items() if not d]
    if eksik:
        out.append(Check("hata", "Tanınmayan enzim",
                         f"{', '.join(eksik)} REBASE'de bulunamadı."))
        return out

    # --- Sıcaklık
    sicakliklar = {n: d.get("opt_temp") for n, d in detaylar.items()}
    bilinen = {n: t for n, t in sicakliklar.items() if t}
    if len(set(bilinen.values())) > 1:
        liste = ", ".join(f"{n}: {t} °C" for n, t in bilinen.items())
        out.append(Check(
            "uyari", "Optimum sıcaklıklar farklı",
            f"{liste}. Aynı tüpte kesecekseniz düşük olanda çalışın ve "
            f"süreyi uzatın; ya da sırayla kesin: önce düşük sıcaklıktaki "
            f"enzimle, sonra tamponu ayarlayıp yüksek olanla. Sıralı "
            f"kesimde aradaki kolon saflaştırması çoğu zaman gereklidir."))
    elif bilinen:
        out.append(Check(
            "bilgi", "Sıcaklıklar uyumlu",
            f"İkisi de {list(bilinen.values())[0]} °C'de çalışıyor, aynı "
            f"tüpte kesilebilir."))

    # --- Metilasyon
    for n, d in detaylar.items():
        if d.get("methylation"):
            out.append(Check(
                "uyari", f"{n} metilasyona duyarlı",
                f"{d['methylation']}. E. coli'den saflaştırdığınız "
                f"plazmitte bu bölge kesilmeyebilir ve 'enzim çalışmadı' "
                f"diye yorumlanır. Plazmiti "
                f"{', '.join(ez.DAM_DCM_FREE_STRAINS[:2])} gibi dam⁻/dcm⁻ "
                f"bir suşta çoğaltın."))

    # --- Uç uyumu
    for i, a in enumerate(adlar):
        for b in adlar[i + 1:]:
            if ez.compatible_with(a) and b in ez.compatible_with(a):
                out.append(Check(
                    "uyari", f"{a} ve {b} uyumlu uç bırakıyor",
                    "Kesilmiş vektörün iki ucu birbirine ligate "
                    "olabildiği için insert ters de girer ve vektör "
                    "kendi üzerine kapanır. Yönlü klonlama istiyorsanız "
                    "uyumsuz uç bırakan bir çift seçin."))

    # --- Kör uç
    kor = [n for n, d in detaylar.items() if d.get("blunt")]
    if len(kor) == len(adlar):
        out.append(Check(
            "uyari", "İki uç da kör",
            "Kör uç ligasyonunda yön kontrolü yoktur ve verim düşüktür. "
            "Mümkünse en az bir yapışkan uç kullanın."))

    out.append(Check("bilgi", "Tampon", SUPPLIER_NOTE))
    out.append(Check(
        "bilgi", "Yıldız aktivitesi",
        "Enzim hacmi toplam reaksiyonun %10'unu geçmesin. Enzim stoğu "
        "gliserol içerir; yoğunlaşınca enzim tanıma dizisine benzeyen "
        "yerlerden de keser (yıldız aktivitesi) ve jelde beklenmedik "
        "bantlar çıkar. Uzun inkübasyon ve yüksek enzim aynı etkiyi "
        "yapar."))
    return out


# ---------------------------------------------------------------------------
# Sorun giderme ağacı
# ---------------------------------------------------------------------------


def symptoms() -> list[Symptom]:
    return [
        Symptom("koloni_yok", "Plakta hiç koloni çıkmadı", [
            Cause("Ligasyon olmadı",
                  "Vektör ve insert uçları birbirine uymuyor olabilir, ya "
                  "da kesim tam değildi. Kesilmiş vektörü jelde kontrol "
                  "edin: tek, keskin bir bant görmelisiniz; kesilmemiş "
                  "plazmit üç bant verir.",
                  "Kesilmiş vektörü jelde yürütün"),
            Cause("Yetkin hücre ölü ya da verimsiz",
                  "Yetkin hücreler donma-çözülmeye çok duyarlıdır. Bir "
                  "kontrol olarak yanına kesilmemiş, sağlam bir plazmit "
                  "(örneğin saf pUC19) transforme edin: ondan da koloni "
                  "çıkmıyorsa sorun hücrededir, ligasyonda değil.",
                  "Pozitif kontrol transformasyonu"),
            Cause("Toparlanma adımı atlandı",
                  "Isı şokundan sonra antibiyotiksiz besiyerinde 37 °C'de "
                  "1 saat beklemeden plakaya yayarsanız, direnç proteini "
                  "henüz üretilmemişken hücre antibiyotikle karşılaşır ve "
                  "ölür. Kanamisin ve kloramfenikolde bu adım zorunludur.",
                  "Protokoldeki toparlanma süresi"),
            Cause("Antibiyotik derişimi yüksek ya da plak eski",
                  "Plakları dökerken antibiyotiği besiyeri 50 °C'nin "
                  "altına indikten sonra ekleyin; sıcakken eklenen "
                  "ampisilin bozulmaz ama bazı antibiyotikler bozulur. "
                  "Aylarca buzdolabında kalan plaklarda da derişim düşer.",
                  "Plak tarihi ve derişim"),
            Cause("Ürün konağa toksik",
                  "Sabit promotörlü ya da kaçak ifade yapan bir yapı, "
                  "protein toksikse hücreyi öldürür. Koloni çıkmaması ya "
                  "da yalnızca mutasyona uğramış yapı taşıyan koloni "
                  "çıkması buna işaret eder.",
                  "Promotör sıkı kontrollü mü"),
        ]),
        Symptom("hepsi_bos", "Koloni çıktı ama hepsi boş vektör", [
            Cause("Vektör kendi üzerine kapandı",
                  "Tek enzimle kestiyseniz vektörün iki ucu birbirine "
                  "uyumludur ve insert olmadan da kapanır. Kesilmiş "
                  "vektörü alkalen fosfatazla (rSAP/CIP) muamele edin: "
                  "5' fosfatı alınan uç ligaz tarafından birleştirilemez.",
                  "Tek enzim mi kullanıldı"),
            Cause("Kesim tam değildi",
                  "Kısmen kesilmiş vektör dairesel kalır ve yüksek verimle "
                  "transforme olur. Kesim süresini uzatın, enzimi tazeleyin "
                  "ve kesilmiş vektörü MUTLAKA jelden saflaştırın.",
                  "Jelden saflaştırma yapıldı mı"),
            Cause("İnsert eksik ya da yanlış molar oranda",
                  "Kütlece eşit koymak yanlıştır; molar oran gerekir. "
                  "Laboratuvar sekmesindeki hesap bunu veriyor.",
                  "Ligasyondaki insert miktarı"),
            Cause("Kontrol kurulmadı",
                  "İnsert yerine su koyduğunuz bir kontrol tüpü, boş "
                  "koloni arka planını doğrudan gösterir. Kontrolde de "
                  "bol koloni varsa sorun vektördedir.",
                  "Su kontrolü"),
        ]),
        Symptom("protein_yok", "Klon doğru ama protein yok", [
            Cause("Konakta T7 polimeraz yok",
                  "T7 ya da T7lac promotörlü yapı, DE3 lizojeni olmayan "
                  "suşta hiç okunmaz. DH5-alfa ve TOP10'da sessiz kalır - "
                  "bu klonlama için avantajdır ama ifade için konağı "
                  "değiştirmeniz gerekir.",
                  "Konak suşu DE3 mü"),
            Cause("mRNA başlangıcı kapalı",
                  "Shine-Dalgarno dizisini ya da başlama kodonunu kapatan "
                  "bir saç tokası, promotör ne kadar güçlü olursa olsun "
                  "çeviriyi başlatmaz. Dizi analizi → mRNA yapısı sekmesi "
                  "bunu tarar.",
                  "mRNA yapısı sekmesi"),
            Cause("Nadir kodon kümesi",
                  "Özellikle genin ilk 30-50 kodonunda yoğunlaşan nadir "
                  "kodonlar ribozomu durdurur. Kodon optimizasyonu "
                  "sekmesi kümeleri gösterir; Rosetta gibi bir suş da "
                  "aynı sorunu başka yoldan çözer.",
                  "Kodon optimizasyonu sekmesi"),
            Cause("Okuma çerçevesi kaymış",
                  "Etiket ile gen arasındaki birleşme noktası üçün katı "
                  "değilse protein ya hiç çıkmaz ya da saçma bir dizi "
                  "olarak çıkar. Tasarım denetimindeki FRM kuralları bunu "
                  "yakalar.",
                  "Denetimdeki FRM uyarıları"),
            Cause("Protein çözünmez halde",
                  "Protein üretiliyor olabilir ama çökelti (inclusion "
                  "body) halindedir ve çözünür kısma bakıyorsanız "
                  "göremezsiniz. Liziz sonrası üst sıvı ile dipteki "
                  "çökeltiyi AYRI AYRI jele yükleyin.",
                  "Çökelti jelde kontrol edildi mi"),
        ]),
        Symptom("yanlis_bant", "Jelde beklenmedik bant boyutu", [
            Cause("Kesilmemiş plazmit",
                  "Dairesel plazmit üç bant verir: süperkıvrımlı (kendi "
                  "boyutundan hızlı göçer), açık dairesel ve doğrusal. "
                  "Kesilmemiş kontrolü yan kuyuya koyarsanız bu üçlüyü "
                  "tanırsınız.",
                  "Kesilmemiş kontrol kuyusu"),
            Cause("Yıldız aktivitesi",
                  "Fazla enzim, fazla gliserol ya da çok uzun inkübasyon, "
                  "enzimin benzer dizilerden de kesmesine yol açar. Enzim "
                  "hacmini toplam reaksiyonun %10'unun altında tutun.",
                  "Enzim hacmi / süre"),
            Cause("Yakın boyutlu fragmanlar tek bant olmuş",
                  "Birbirine %5-10 yakın fragmanlar jelde ayrışmaz. "
                  "Laboratuvar → Jel simülasyonu bunu önceden gösterir.",
                  "Jel simülasyonu"),
            Cause("Jel yüzdesi uygun değil",
                  "Çok büyük fragmanlar yüksek yüzdeli jelde kuyunun "
                  "altında sıkışır, çok küçükler düşük yüzdelide boya "
                  "cephesiyle çıkar.",
                  "Agaroz yüzdesi"),
        ]),
        Symptom("plazmit_kayboluyor", "Plazmit kararsız, kayboluyor", [
            Cause("İki plazmit aynı uyumsuzluk grubunda",
                  "ColE1 ile pUC aynı gruptandır; aynı hücrede uzun süre "
                  "duramazlar. İkinci plazmit için p15A ya da pSC101 "
                  "orijini seçin. Tasarım denetimindeki ORI kuralları "
                  "bunu kontrol eder.",
                  "Denetimdeki ORI uyarıları"),
            Cause("Metabolik yük",
                  "Yüksek kopyalı orijin + güçlü sabit promotör, hücreyi "
                  "yorar ve plazmiti kaybeden ya da promotörü mutasyona "
                  "uğratan hücreler seçilir. Düşük kopyalı orijine ya da "
                  "indüklenebilir promotöre geçin.",
                  "Orijin kopya sayısı"),
            Cause("Tekrar dizileri rekombinasyona uğruyor",
                  "Uzun doğrudan tekrarlar recA⁺ suşta birbirine "
                  "rekombine olur ve aradaki bölge silinir. Klonlamayı "
                  "recA1 taşıyan bir suşta (DH5-alfa, TOP10) yapın.",
                  "Konak recA durumu"),
        ]),
    ]


def evaluate(
    symptom_key: str,
    part_ids: list[str] | None = None,
    host: dict[str, Any] | None = None,
    enzymes: list[str] | None = None,
    purpose: str = "ifade",
) -> Symptom:
    """Belirtinin sebeplerini sırala ve mümkün olanları tasarımdan işaretle.

    `applies` alanı True ise sebep bu tasarımda GERÇEKTEN var demektir;
    None ise uygulama buna bakamıyor, kullanıcının kontrol etmesi gerek.
    """
    pids = set(part_ids or [])
    enz = [e for e in (enzymes or []) if e]
    hedef = next((s for s in symptoms() if s.key == symptom_key), None)
    if hedef is None:
        return Symptom(symptom_key, symptom_key, [])

    for c in hedef.causes:
        if c.title == "Vektör kendi üzerine kapandı":
            c.applies = len(enz) < 2 if enz else None
        elif c.title == "Konakta T7 polimeraz yok":
            t7_parca = any(p in pids for p in ("PRO_T7", "PRO_T7LAC"))
            if t7_parca and host is not None:
                c.applies = not host.get("t7_polymerase", False)
            elif not t7_parca:
                c.applies = False
        elif c.title == "İki plazmit aynı uyumsuzluk grubunda":
            c.applies = None     # ikinci plazmiti uygulama bilmiyor
        elif c.title == "Metabolik yük":
            yuksek = "ORI_PUC" in pids
            sabit = any(p.startswith(("PRO_J23", "IGEM_J231")) for p in pids)
            c.applies = bool(yuksek and sabit)
        elif c.title == "Tekrar dizileri rekombinasyona uğruyor":
            if host is not None:
                c.applies = bool(host.get("recA", True))
        elif c.title == "Ürün konağa toksik":
            sabit = any(p.startswith(("PRO_J23", "IGEM_J231")) for p in pids)
            c.applies = True if (sabit and purpose == "ifade") else None
        elif c.title == "Nadir kodon kümesi":
            c.applies = None
    return hedef
