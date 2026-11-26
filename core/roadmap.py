"""Gerçekleştirme yol haritası: tasarımdan tedavi adayına.

Tasarım bitince "artık tedavi var" olmuyor. Arada yıllar süren,
basamakları belli bir yol var ve bu yolun en büyük sorunu şudur:
basamakların **sırası** keyfi değildir. Hayvan deneyine etik kurul
onayından önce başlanmaz; içerme testi yapılmadan canlı GDO ile
hayvan çalışılmaz; ifade gösterilmeden işlev aranmaz. Sırayı bilmemek
bir projeyi, yanlış bir promotör seçmekten daha kesin batırır.

Bu modül o yolu **tasarımın kendisinden** üretir. Statik bir metin
değil: kasette sinyal peptidi varsa salgılanma doğrulaması eklenir,
konak Gram-pozitifse dönüşüm optimizasyonu gelir, markerde antibiyotik
direnci varsa "in vivo öncesi değiştirilmeli" riski açılır, tasarım
hastalık kayıt defterinden geldiyse in vivo okuma o hastalığın
biyobelirteci olur.

**Tezgâh protokolü burada tekrarlanmaz.** Kesme, ligasyon,
transformasyon, koloni taraması ve dizileme `core/protocol.py`
içindedir ve Laboratuvar sayfası onu basar. Yol haritası oradan
SONRASINI anlatır; birinci aşamada yalnızca o protokole yönlendirir.
İki yerde iki ayrı protokol tutmak, biri değiştiğinde öbürünün
sessizce yanlış kalması demektir.

**Süreler tahmindir.** Kaynaktan gelmez; yerleşik akademik iş
akışının alışıldık süreleridir, `PHASE_DAYS` içinde tek yerde durur
ve `build(...)` çağrısında değiştirilir. Bir doktora öğrencisinin tek
başına çalıştığı varsayılır; ekip büyüdükçe kısalır, ilk denemede
tutmayan adım uzatır.

**Bu bir düzenleyici danışmanlık değildir.** Yasal dayanaklar
aşağıda kaynağıyla yazılıdır ama mevzuat değişir - nitekim
Biyogüvenlik Kurulu 2018'de kaldırılmıştır ve bunu ezberden yazan bir
başvuru, var olmayan bir kurula başvurmuş olur. Her izin için güncel
mevzuat ve kurum kendiniz doğrulanmalıdır.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core import db

# ---------------------------------------------------------------------------
# Yasal dayanaklar - kaynağıyla, ve "doğrulayın" uyarısıyla
# ---------------------------------------------------------------------------


@dataclass
class Approval:
    """Alınması gereken bir izin ve yasal dayanağı."""

    key: str
    title: str
    body: str                  # hangi kuruma
    basis: str                 # yasal dayanak
    when: str                  # hangi aşamadan ÖNCE
    note: str = ""
    caution: str = ""


#: İzinler. Her birinin dayanağı Resmî Gazete künyesiyle yazılı.
#: Künyeler arama ile doğrulandı; buna rağmen "güncel halini kontrol
#: edin" uyarısı her birinde duruyor, çünkü mevzuat değişir ve bu
#: alanda bir kez değişmiştir.
APPROVALS: dict[str, Approval] = {
    "hadyek": Approval(
        key="hadyek",
        title="Hayvan deneyi etik kurul onayı (HADYEK)",
        body="Kurumunuzun Hayvan Deneyleri Yerel Etik Kurulu",
        basis=("Hayvan Deneyleri Etik Kurullarının Çalışma Usul ve "
               "Esaslarına Dair Yönetmelik (R.G. 15.02.2014, 28914); "
               "dayanağı 5199 sayılı Hayvanları Koruma Kanunu"),
        when="Hayvan çalışmasına başlamadan ÖNCE",
        note=("Deney hayvanı tanımı insan dışındaki bütün omurgalıları "
              "ve canlı kafadanbacaklıları kapsar. TÜBİTAK projelerinde "
              "deney hayvanı kullanılacaksa onay belgesi başvuru "
              "ekinde istenir - yani proje başvurusundan önce "
              "başlatılması gereken bir süreçtir.\n\n"
              "Başvuru formu proje adını, yürütücüyü, uygulama yerini, "
              "hayvan türü ve sayısını, ve **neden o sayıda hayvan "
              "gerektiğini** sorar. Güç analizi olmadan yazılan hayvan "
              "sayısı en sık geri dönen başvuru sebebidir."),
        caution=("Kurul takvimleri aylıktır: onay beklemek haftalar "
                 "alır ve bu süre projenin takvimine YAZILMALIDIR.")),
    "gdo": Approval(
        key="gdo",
        title="Genetiği değiştirilmiş organizma (GDO) izni",
        body="Tarım ve Orman Bakanlığı",
        basis=("5977 sayılı Biyogüvenlik Kanunu "
               "(kabul 18.03.2010, R.G. 26.03.2010, 27533) ve "
               "uygulama yönetmelikleri"),
        when="GDO'yu kapalı alan dışında kullanmadan ÖNCE",
        note=("Kanun, yurt içinde geliştirilen GDO ve ürünleri için "
              "gerçek ve tüzel kişilerin Bakanlığa başvurmasını "
              "öngörür. Araştırma amaçlı kapalı alan çalışmaları ile "
              "çevreye salım ve piyasaya sürme farklı rejimlere "
              "tabidir; hangi rejimde olduğunuzu baştan belirlemek "
              "gerekir."),
        caution=("**Kurum adını ezberden yazmayın.** Biyogüvenlik "
                 "Kurulu 09.07.2018 tarihli 703 sayılı KHK ile "
                 "kaldırılmış ve görevleri Bakanlık içine "
                 "aktarılmıştır. Başvuruyu yapacağınız güncel birimi "
                 "ve yürürlükteki yönetmeliği başvuru öncesinde "
                 "Bakanlıktan teyit edin - bu bilgi bir kez "
                 "değişmiştir, yine değişebilir.")),
    "ibg": Approval(
        key="ibg",
        title="Kurumsal biyogüvenlik onayı ve laboratuvar seviyesi",
        body="Kurumunuzun biyogüvenlik kurulu / sorumlusu",
        basis="Kurumsal biyogüvenlik yönergesi",
        when="Laboratuvar çalışmasına başlamadan ÖNCE",
        note=("Çalışılacak organizmanın risk grubu, laboratuvarın "
              "biyogüvenlik seviyesini (BSL) belirler. *E. coli* K-12 "
              "türevleri ve gıda sınıfı laktik asit bakterileri "
              "genellikle en düşük seviyede çalışılır; patojen "
              "kaynaklı bir gen (örneğin bir bakteriyosin ya da "
              "virülans düzenleyicisi) taşıyan yapı için kurumunuzun "
              "değerlendirmesi gerekir.\n\n"
              "Atık imhası da bu onayın parçasıdır: GDO taşıyan katı "
              "ve sıvı atık otoklavlanmadan çöpe ya da lavaboya "
              "verilmez."),
        caution=""),
    "klinik": Approval(
        key="klinik",
        title="Klinik araştırma izni (bu projenin kapsamı dışında)",
        body="Türkiye İlaç ve Tıbbi Cihaz Kurumu (TİTCK)",
        basis="İlaç ve biyolojik ürün klinik araştırma mevzuatı",
        when="İnsanda herhangi bir uygulamadan ÖNCE",
        note=("Canlı bakteri terapötiği insanda denenecekse klinik "
              "araştırma iznine, ondan önce de klinik öncesi "
              "güvenlik paketine ve GMP koşullarında üretime ihtiyaç "
              "vardır. Bir lisans ya da yüksek lisans projesinin "
              "kapsamı burayı içermez; yol haritasında durmasının "
              "nedeni, projenin **nereye kadar** gittiğini "
              "başvuruda doğru anlatabilmenizdir."),
        caution=("Başvuru metninde 'tedavi geliştiriyoruz' demek ile "
                 "'tedavi adayının klinik öncesi ilk basamağını "
                 "kuruyoruz' demek arasındaki fark, hakemin "
                 "ciddiye alıp almamasını belirler.")),
}


# ---------------------------------------------------------------------------
# Görev ve aşama
# ---------------------------------------------------------------------------


@dataclass
class Task:
    """Tek bir iş: ne, nasıl, neyle, çıktısı ne.

    `nasil` alanı yöntemin adını verir (SDS-PAGE, Western blot, CFU
    sayımı...). `cikti` o işin BİTTİĞİNİ kanıtlayan şeydir - bir
    görev çıktısı tanımlanmadan takvime yazılırsa ne zaman bittiği
    hiç belli olmaz.
    """

    key: str
    title: str
    nasil: str
    gereken: str
    cikti: str
    gun: float
    neden: str = ""
    kosul: str = ""            # bu tasarımda NEDEN var
    approval: str = ""         # önce alınması gereken izin


@dataclass
class Phase:
    key: str
    no: int
    title: str
    amac: str
    tasks: list[Task] = field(default_factory=list)
    kapi: str = ""             # bir sonrakine geçmek için ne olmalı

    @property
    def gun(self) -> float:
        return sum(t.gun for t in self.tasks)


@dataclass
class Risk:
    """Projeyi durdurabilecek şey ve ne yapılacağı."""

    key: str
    title: str
    etki: str
    onlem: str
    kaynak: str = ""           # hangi denetim kuralından geldi


@dataclass
class Roadmap:
    design_name: str
    host_name: str
    phases: list[Phase] = field(default_factory=list)
    risks: list[Risk] = field(default_factory=list)
    approvals: list[Approval] = field(default_factory=list)
    readout: str = ""          # başarıyı ölçecek biyobelirteç
    notes: list[str] = field(default_factory=list)

    @property
    def gun(self) -> float:
        return sum(p.gun for p in self.phases)

    @property
    def ay(self) -> float:
        """Takvim ayı: 21 iş günü/ay."""
        return round(self.gun / 21.0, 1)

    def phase(self, key: str) -> Phase | None:
        return next((p for p in self.phases if p.key == key), None)


# ---------------------------------------------------------------------------
# Süre varsayımları - tek yerde
# ---------------------------------------------------------------------------

#: İş günü. Tek kişinin çalıştığı varsayılır. Kaynaktan gelmez.
PHASE_DAYS: dict[str, float] = {
    "klonlama": 15.0,
    "dizileme": 5.0,
    "suşa_aktarim": 5.0,
    "ifade": 10.0,
    "salgilanma": 8.0,
    "etkinlik": 12.0,
    "kararlilik": 15.0,
    "yuk": 5.0,
    "sindirim": 10.0,
    "hucre": 20.0,
    "sitotoksisite": 8.0,
    "icerme": 15.0,
    "etik": 30.0,
    "hayvan_kurulum": 20.0,
    "hayvan_uygulama": 30.0,
    "hayvan_olcum": 20.0,
    "gdo": 40.0,
    "yayim": 20.0,
}


#: İnsan kaynaklı terapötik yükler. Bakteri vektörü kapsamının
#: belgelenmiş istisnası; kodon uyumu ve katlanma ek iş demektir.
INSAN_CDS = {"CDS_IL10", "CDS_ELAFIN", "CDS_GLP1", "CDS_TFF1"}

#: Terapötik bağlamda çalışan konaklar.
TERAPOTIK_KONAK = {"ECN", "ECN_DAPA", "LLNZ9000"}

#: Antibiyotik direnci taşıyan markerlar.
ANTIBIYOTIK_MARKER = {"MRK_BLA", "MRK_KAN", "MRK_CAT", "MRK_TET",
                      "MRK_AADA", "REF_ERM_C_735"}

#: Antibiyotiksiz idame markerları.
ICERME_MARKER = {"MRK_THYA_LL", "MRK_DAPA", "MRK_CCDA", "MRK_CCDB"}


def _g(ad: str, ozel: dict[str, float] | None = None) -> float:
    return (ozel or {}).get(ad, PHASE_DAYS[ad])


# ---------------------------------------------------------------------------
# Yol haritasını kur
# ---------------------------------------------------------------------------


def build(design_name: str, part_ids: list[str], host_id: str = "",
          disease: Any = None, findings: list[Any] | None = None,
          in_vivo: bool = True,
          phase_days: dict[str, float] | None = None) -> Roadmap:
    """Tasarımdan gerçekleştirme yol haritası üret.

    `disease`, `core/disease_seed.Disease` nesnesidir (ya da None):
    verilirse in vivo okuma o hastalığın biyobelirteci olur ve
    "bilinen açık" notu riske dönüşür.

    `findings`, `core/validate.validate()` çıktısıdır: her hata ve
    uyarı bir proje riskine çevrilir. Denetimde görünen bir eksiğin
    başvuruda risk olarak yazılmaması, hakemin ilk yakalayacağı
    şeydir.
    """
    parts = {p["part_id"]: p for p in db.get_parts(part_ids)}
    kategoriler = {p["category"] for p in parts.values()}
    konak = db.get_host(host_id) if host_id else None
    konak_ad = (konak or {}).get("name", host_id or "—")
    organizma = (konak or {}).get("organism", "")
    gram_poz = organizma in ("Lactococcus lactis", "Bacillus subtilis",
                             "Corynebacterium glutamicum")
    terapotik = bool(disease) or host_id in TERAPOTIK_KONAK

    sinyal = [p for p in part_ids if parts.get(p, {}).get("category")
              == "signal"]
    insan = [p for p in part_ids if p in INSAN_CDS]
    antibiyotik = [p for p in part_ids if p in ANTIBIYOTIK_MARKER]
    icerme = [p for p in part_ids if p in ICERME_MARKER]
    cds = [p for p in part_ids if parts.get(p, {}).get("category") == "cds"]

    yol = Roadmap(design_name=design_name, host_name=konak_ad)

    # ---------------------------------------------------------------
    # Aşama 1 - yapıyı kur
    # ---------------------------------------------------------------
    a1 = Phase(
        key="kurulum", no=1, title="Yapıyı kur ve dizisini doğrula",
        amac=("Kâğıttaki tasarımı tüpteki plazmite çevirmek. Buradan "
              "sonraki her aşama bu plazmitin DOĞRU olduğuna "
              "dayanır; yanlış bir yapıyla yapılan altı aylık iş "
              "boşa gider."),
        kapi=("Dizileme okuması tasarımla birebir tutmalı. Tek baz "
              "farkı bile 'küçük fark' değildir: çerçeve içinde bir "
              "insersiyon proteini tamamen bozar ve haritada "
              "görünmez."))
    a1.tasks.append(Task(
        "klonla", "Klonlama ve koloni doğrulaması",
        nasil=("Tezgâh protokolü **Laboratuvar sayfasında** hazır: "
               "kesim, jelden saflaştırma, ligasyon, transformasyon, "
               "koloni PCR. Reaksiyon miktarları ve antibiyotik "
               "derişimleri tasarımınızdan hesaplanır."),
        gereken=("Sipariş ve denetim izi sayfasındaki malzeme "
                 "listesi: oligolar, kesim enzimleri, yetkin hücre, "
                 "kitler."),
        cikti="Beklenen bant boyunu veren koloni PCR jeli",
        gun=_g("klonlama", phase_days),
        neden=("Koloni PCR'ı atlayıp doğrudan dizilemeye gitmek "
               "pahalıdır: on koloninin dokuzu yanlışsa dokuz "
               "dizileme parası ödersiniz."),
        approval="ibg"))
    a1.tasks.append(Task(
        "dizile", "Sanger dizilemesiyle doğrulama",
        nasil=("Doğrulama sayfası kaset uzunluğundan kaç okuma "
               "gerektiğini hesaplar ve okunan diziyi tasarımla "
               "hizalayıp farkları gösterir. Tek okuma 700-900 bazdan "
               "sonra güvenilirliğini kaybeder."),
        gereken="Dizileme primerleri (Primer tasarımı sayfası), plazmit",
        cikti="Tasarımla %100 örtüşen, birleştirilmiş okuma",
        gun=_g("dizileme", phase_days),
        neden=("Ticari dizileme hizmetinin verdiği ham okumayı gözle "
               "karşılaştırmak hata kaynağıdır; hizalama aracı kullanın.")))

    if host_id and (gram_poz or terapotik):
        a1.tasks.append(Task(
            "aktar", f"Yapıyı {konak_ad} suşuna aktar",
            nasil=("Elektroporasyon. Gram-pozitif konaklarda hücre "
                   "duvarı kalın olduğu için dönüşüm verimi "
                   "E. coli'nin çok altındadır: büyütme sırasında "
                   "glisin ya da DL-treonin ile duvar zayıflatması ve "
                   "yüksek sükroz içeren yıkama tamponu yerleşik "
                   "yaklaşımdır. Verim yine de düşük kalırsa "
                   "konjugasyonla aktarım denenir."
                   if gram_poz else
                   "Elektroporasyon ya da kimyasal yetkinlik; "
                   "terapötik suş laboratuvar suşundan daha zor "
                   "dönüşebilir."),
            gereken=("Elektroporatör, suşa uygun yetkin hücre "
                     "hazırlama protokolü"),
            cikti="Seçilim plağında koloni ve plazmitin geri alınması",
            gun=_g("suşa_aktarim", phase_days),
            neden=("Yapı E. coli'de kurulur ama orada kalmaz. Bu "
                   "basamak projelerin takvimde en çok hafife aldığı "
                   "yerdir: bir suşa ilk kez plazmit sokmak haftalar "
                   "alabilir."),
            kosul=(f"{organizma} Gram-pozitif, dönüşümü zor"
                   if gram_poz else
                   "Terapötik suşa aktarım gerekiyor")))
    yol.phases.append(a1)

    # ---------------------------------------------------------------
    # Aşama 2 - ifade ve işlev
    # ---------------------------------------------------------------
    a2 = Phase(
        key="ifade", no=2, title="İfadeyi ve işlevi göster (in vitro)",
        amac=("Plazmitin doğru olması, proteinin üretildiği anlamına "
              "gelmez; proteinin üretilmesi de ÇALIŞTIĞI anlamına "
              "gelmez. Bu aşama ikisini ayrı ayrı kanıtlar."),
        kapi=("Hem ifade (bant) hem etkinlik (ölçülebilir dönüşüm) "
              "gösterilmeden sonraki aşamaya geçilmez. Etkinliği "
              "olmayan bir proteini hayvana vermek, sonucu baştan "
              "belli bir deneydir."))
    a2.tasks.append(Task(
        "ifade", "Protein ifadesini göster",
        nasil=("İndüklenmiş ve indüklenmemiş kültürden eşit hücre "
               "alıp SDS-PAGE; beklenen boyda bant aranır. Bant "
               "belirsizse etikete ya da antiserumla Western blot. "
               "Çözünürlük için hücreyi parçalayıp çözünen ve "
               "çökelti kısımlarını AYRI yürütün: inklüzyon "
               "cisimciğinde biriken protein jelde bol görünür ama "
               "işlevsizdir."),
        gereken="SDS-PAGE düzeneği, protein boyası, (gerekirse) antikor",
        cikti=("İndüklenmişte olup indüklenmemişte olmayan, beklenen "
               "boyda bant; çözünen kısımda olduğu gösterilmiş"),
        gun=_g("ifade", phase_days),
        neden=("İfade edilmeyen bir kasetin sebebi promotör, RBS, "
               "kodon kullanımı ya da toksisite olabilir - hangisi "
               "olduğunu ancak bu basamakta anlarsınız.")))

    if sinyal:
        a2.tasks.append(Task(
            "salgi", "Salgılanmayı doğrula",
            nasil=("Kültür süpernatantını hücreden ayırıp ayrı "
                   "yürütün. **Kritik kontrol:** hücre içi bir "
                   "proteini (örneğin sitoplazmik bir enzim) aynı "
                   "süpernatantta aramak. O da çıkıyorsa gördüğünüz "
                   "salgılanma değil, hücre lizisidir - bu, salgı "
                   "çalışmalarının en sık yanılgısıdır."),
            gereken="Santrifüj, süpernatant yoğunlaştırma, kontrol antiserumu",
            cikti=("Süpernatantta terapötik protein VAR, lizis "
                   "göstergesi protein YOK"),
            gun=_g("salgilanma", phase_days),
            neden=("Sinyal peptidi koymak salgılandığını kanıtlamaz; "
                   "sinyal işlenmemiş olabilir ya da protein zarda "
                   "sıkışabilir."),
            kosul=f"Kasette sinyal peptidi var ({', '.join(sinyal)})"))

    if cds:
        a2.tasks.append(Task(
            "etkinlik", "Terapötik etkinliği ölç",
            nasil=("Substrat tüketimi ya da ürün oluşumu üzerinden "
                   "enzim tayini; sitokin ya da peptit terapötiğinde "
                   "ELISA ile miktar, ardından hücre temelli etkinlik "
                   "testi. Ölçüm BİRİM ile verilmelidir (ör. birim "
                   "enzim/mg toplam protein, ya da ng/mL/OD)."),
            gereken="Spektrofotometre ya da plaka okuyucu, substrat/ELISA kiti",
            cikti=("Boş vektör taşıyan kontrol suşuna göre anlamlı "
                   "fark gösteren, birimi yazılı etkinlik değeri"),
            gun=_g("etkinlik", phase_days),
            neden=("Boş vektör kontrolü olmadan ölçülen etkinlik, "
                   "konağın kendi arka planından ayırt edilemez.")))

    a2.tasks.append(Task(
        "yuk", "Metabolik yükü ölç",
        nasil=("Yapıyı taşıyan ve taşımayan suşun büyüme eğrilerini "
               "aynı koşulda karşılaştır (OD600, en az üç tekrar). "
               "Üreme hızı ve son yoğunluk ayrı ayrı bakılır."),
        gereken="Plaka okuyucu ya da spektrofotometre",
        cikti="Büyüme eğrisi karşılaştırması ve üreme hızı farkı",
        gun=_g("yuk", phase_days),
        neden=("Ağır bir kaset taşıyan hücre, bağırsakta taşımayan "
               "rakipleri tarafından birkaç nesilde elenebilir. "
               "Yükü ölçmemek, in vivo'da sebebi anlaşılmayan "
               "başarısızlık demektir.")))
    yol.phases.append(a2)

    # ---------------------------------------------------------------
    # Aşama 3 - kararlılık, içerme ve model ortam
    # ---------------------------------------------------------------
    a3 = Phase(
        key="kararlilik", no=3,
        title="Kararlılık, içerme ve benzetilmiş koşullar",
        amac=("Laboratuvar şişesinde çalışan bir suşun bağırsakta "
              "çalışacağı varsayılamaz. Bu aşama, canlı bir ilacın "
              "laboratuvar suşundan ayrıldığı yerdir."),
        kapi=("Plazmit antibiyotiksiz pasajda tutulmalı ve içerme "
              "mekanizmasının kaçak oranı ölçülmüş olmalı. İkisi de "
              "hayvan çalışmasının ön koşuludur."))
    a3.tasks.append(Task(
        "kararlilik", "Plazmit kararlılığını antibiyotiksiz ölç",
        nasil=("Suşu antibiyotiksiz besiyerinde seri pasajla (tipik "
               "olarak 50-100 nesil) büyüt; belirli aralıklarla "
               "seyreltip seçilimli ve seçilimsiz plağa ek, "
               "koloni oranını say. Kayıp varsa nesil başına kayıp "
               "oranını hesapla."),
        gereken="Seri pasaj kültürü, iki tip plak, koloni sayımı",
        cikti="Nesil sayısına karşı plazmit taşıyan hücre yüzdesi eğrisi",
        gun=_g("kararlilik", phase_days),
        neden=("Antibiyotikle büyütülen bir kültürde plazmit elbette "
               "durur; soru, antibiyotik YOKKEN ne olduğudur. "
               "Hastaya antibiyotik verilmeyecek.")))

    if terapotik:
        a3.tasks.append(Task(
            "sindirim", "Benzetilmiş mide-bağırsak koşullarına dayanım",
            nasil=("Düşük pH'lı benzetilmiş mide sıvısında ve safra "
                   "tuzu içeren benzetilmiş bağırsak sıvısında "
                   "belirli sürelerle bekletip canlı hücre sayımı "
                   "(CFU/mL). Başlangıç ve son sayımın oranı "
                   "verilir."),
            gereken="pH ayarlı tamponlar, safra tuzu, pepsin/pankreatin",
            cikti="Her koşul için log azalma değeri",
            gun=_g("sindirim", phase_days),
            neden=("Ağızdan verilen bir bakteri mideyi geçemezse "
                   "bağırsakta hiç bulunmaz; kasetin ne kadar iyi "
                   "olduğunun önemi kalmaz."),
            kosul="Tasarım bağırsakta çalışacak bir terapötik"))
        a3.tasks.append(Task(
            "icerme", "Biyolojik içermeyi sına",
            nasil=("İçerme mekanizmasının izin verdiği koşulu "
                   "kaldırıp (timin ya da DAP içermeyen besiyeri) "
                   "büyük hücre sayısını plağa ekerek **kaçak "
                   "oranını** ölç: kaçak = büyüyen koloni / ekilen "
                   "hücre. Kaçan kolonileri ayrıca dizileyip nasıl "
                   "kaçtıklarını belirleyin."),
            gereken="Eksiltilmiş besiyeri, yüksek yoğunluklu ekim",
            cikti="Hücre başına kaçak oranı (bilimsel gösterimle)",
            gun=_g("icerme", phase_days),
            neden=("İçerme mekanizmasının varlığı değil, KAÇAK ORANI "
                   "raporlanır. Bir mekanizma koymak onun çalıştığını "
                   "göstermez; kaçak oranı ölçülmeden hiçbir "
                   "düzenleyici çerçeveye girilmez."),
            kosul=(f"Tasarımda içerme markeri var "
                   f"({', '.join(icerme)})" if icerme else
                   "Terapötik canlı GDO: içerme ölçülmeden "
                   "hayvana verilemez"),
            approval="gdo"))
        a3.tasks.append(Task(
            "hucre", "Bağırsak epiteli hücre modelinde ko-kültür",
            nasil=("Farklılaşmış bağırsak epiteli hücre hattı "
                   "(yaygın olarak Caco-2 ya da HT-29) üzerinde "
                   "ko-kültür. Bariyer bütünlüğü transepitelyal "
                   "elektriksel direnç (TEER) ile, iltihap yanıtı "
                   "salgılanan sitokinlerle, terapötik etki hastalığa "
                   "özgü okuma ile izlenir."),
            gereken="Hücre kültürü altyapısı, TEER ölçer, ELISA",
            cikti=("Kontrol suşuna göre bariyer ve/veya iltihap "
                   "göstergelerinde fark"),
            gun=_g("hucre", phase_days),
            neden=("Hayvan deneyinden önce hücre modelinde bir etki "
                   "göstermek hem etik kurul başvurusunu "
                   "güçlendirir hem de hayvan sayısını azaltır."),
            approval="ibg"))
        a3.tasks.append(Task(
            "sitotoksisite", "Sitotoksisite ve güvenlik taraması",
            nasil=("Hücre canlılığı testi (MTT/XTT ya da LDH "
                   "salınımı) ile suşun ve süpernatantının "
                   "toksisitesi; terapötik protein saflaştırılabiliyorsa "
                   "ayrıca saf protein üzerinde."),
            gereken="Plaka okuyucu, canlılık testi kiti",
            cikti="Doz-yanıt eğrisi ve toksik olmayan aralık",
            gun=_g("sitotoksisite", phase_days),
            neden=("Terapötik olduğu varsayılan bir molekülün yüksek "
                   "derişimde toksik çıkması sık görülür.")))
    yol.phases.append(a3)

    # ---------------------------------------------------------------
    # Aşama 4 - hayvan modeli
    # ---------------------------------------------------------------
    if in_vivo and terapotik:
        a4 = Phase(
            key="invivo", no=4, title="Hayvan modelinde etkinlik",
            amac=("Tedavi iddiasının sınandığı yer burasıdır. "
                  "Öncesindeki her aşama, bu deneyi anlamlı kılmak "
                  "içindi."),
            kapi=("Önceden belirlenmiş birincil sonlanım noktasında "
                  "istatistiksel olarak anlamlı fark. Sonlanım "
                  "noktası deneyden SONRA seçilirse sonuç "
                  "yayımlanabilir olmaz."))
        a4.tasks.append(Task(
            "etik", "Etik kurul onayını al",
            nasil=("HADYEK başvuru formu: hayvan türü, sayısı, "
                   "deney düzeni, ağrı kategorisi, insancıl sonlanım "
                   "ölçütleri. Hayvan sayısı **güç analiziyle** "
                   "gerekçelendirilir."),
            gereken="Başvuru formu, güç analizi, deney protokolü",
            cikti="HADYEK onay belgesi (proje başvurusunun eki)",
            gun=_g("etik", phase_days),
            neden=("Onaysız yapılan hayvan çalışması yayımlanamaz ve "
                   "hukuki sorumluluk doğurur. Kurul takvimi aylık "
                   "olduğu için bu süre projeye YAZILMALIDIR."),
            approval="hadyek"))
        a4.tasks.append(Task(
            "model", "Hastalık modelini kur ve doğrula",
            nasil=("Model önce TEDAVİSİZ kurulur ve hastalığın "
                   "kendisinin oluştuğu gösterilir. Sağlıklı kontrol "
                   "ile model grubu arasında beklenen farkı "
                   "göremiyorsanız tedaviyi sınayacak bir zemin de "
                   "yok demektir."),
            gereken="Hayvan, model ajanı, barınak koşulları",
            cikti="Model ve sağlıklı grup arasında doğrulanmış fark",
            gun=_g("hayvan_kurulum", phase_days),
            neden=("Modeli tedaviyle aynı anda kurmak, sonuç "
                   "çıkmadığında sebebin model mi tedavi mi olduğunu "
                   "belirsiz bırakır."),
            approval="hadyek"))
        a4.tasks.append(Task(
            "uygula", "Suşu uygula ve kolonizasyonu izle",
            nasil=("Ağızdan gavaj ya da içme suyuyla uygulama. "
                   "Dışkıdan seri CFU sayımı ve plazmite özgü qPCR "
                   "ile suşun ne kadar kaldığı izlenir. "
                   "**Gruplar:** sağlıklı kontrol, hastalık + boş "
                   "vektörlü suş, hastalık + terapötik suş. Boş "
                   "vektör grubu olmadan, görülen etkinin bakterinin "
                   "kendisinden mi kasetten mi geldiği ayrılamaz."),
            gereken="Gavaj düzeneği, seçici besiyeri, qPCR",
            cikti="Zamana karşı CFU/g dışkı eğrisi; grup başına n",
            gun=_g("hayvan_uygulama", phase_days),
            neden=("Nissle gibi suşlar insan bağırsağını kalıcı "
                   "kolonize etmez; doz tekrarı gerekir ve bu, "
                   "tedavi tasarımının parçasıdır."),
            approval="hadyek"))
        okuma = (disease.biomarker if disease and
                 getattr(disease, "biomarker", "") else
                 "Hastalığa özgü, önceden tanımlanmış birincil "
                 "sonlanım noktası")
        a4.tasks.append(Task(
            "olc", "Birincil sonlanım noktasını ölç",
            nasil=(f"Birincil okuma: {okuma}. Ölçüm körleştirilerek "
                   f"yapılır ve istatistiksel plan deneyden ÖNCE "
                   f"yazılır."),
            gereken="Biyokimya/histoloji altyapısı, istatistik planı",
            cikti=("Gruplar arası fark, etki büyüklüğü ve güven "
                   "aralığıyla"),
            gun=_g("hayvan_olcum", phase_days),
            neden=("p değeri tek başına yetmez; etki büyüklüğü "
                   "klinik anlamlılığın göstergesidir."),
            approval="hadyek"))
        a4.tasks.append(Task(
            "cevre", "Çevresel kaçışı ve yatay gen aktarımını ölç",
            nasil=("Dışkı ve kafes altlığı örneklerinde suşun "
                   "çevrede ne kadar canlı kaldığını izle; "
                   "bağırsak mikrobiyotasına plazmit aktarımını "
                   "plazmite özgü qPCR ve seçilimli kültürle ara."),
            gereken="Örnekleme planı, qPCR",
            cikti="Kaçış ve aktarım verisi (negatif sonuç da sonuçtur)",
            gun=_g("hayvan_olcum", phase_days) / 2,
            neden=("GDO izin dosyasının çevresel risk bölümü bu "
                   "veriyle doldurulur. Sonradan toplanamaz: "
                   "örnekler hayvan çalışması sırasında alınır."),
            approval="gdo"))
        yol.phases.append(a4)

    # ---------------------------------------------------------------
    # Aşama 5 - yaygınlaştırma ve sonrası
    # ---------------------------------------------------------------
    a5 = Phase(
        key="yayim", no=len(yol.phases) + 1,
        title="Yaygınlaştırma ve sonraki basamak",
        amac=("Projenin çıktısının başkaları tarafından "
              "kullanılabilir hale gelmesi ve bir sonraki aşamanın "
              "ne olduğunun dürüstçe yazılması."),
        kapi="")
    a5.tasks.append(Task(
        "veri", "Veri ve dizileri erişime aç",
        nasil=("Plazmit dizisini aksesyon alacak şekilde bir dizi "
               "veritabanına yatır; ham veriyi ve analiz kodunu "
               "kalıcı tanımlayıcısı olan bir depoda yayımla. "
               "Tasarımın denetim izi dosyası (Sipariş ve denetim "
               "izi sayfası) hangi veri sürümüyle üretildiğini "
               "gösterir."),
        gereken="Dizi kaydı, veri deposu",
        cikti="Aksesyon numarası ve kalıcı veri bağlantısı",
        gun=_g("yayim", phase_days) / 2,
        neden=("Yeniden üretilebilirlik, proje değerlendirmesinde "
               "giderek daha çok sorulan bir ölçüt.")))
    a5.tasks.append(Task(
        "yayin", "Bulguları yayımla ve sonraki basamağı tanımla",
        nasil=("Olumsuz sonuçlar da yazılır. Bir sonraki basamak "
               "(klinik öncesi güvenlik paketi, GMP üretim, klinik "
               "araştırma izni) projenin KAPSAMI DIŞINDADIR ve "
               "başvuruda böyle söylenmelidir."),
        gereken="—",
        cikti="Yayın ve gerçekçi bir sonraki adım tanımı",
        gun=_g("yayim", phase_days) / 2,
        neden=("'Tedavi geliştiriyoruz' ile 'tedavi adayının klinik "
               "öncesi ilk basamağını kuruyoruz' arasındaki fark, "
               "hakemin başvuruyu ciddiye alıp almamasını "
               "belirler."),
        # Klinik yol yalnızca terapötik bağlamda anılır. Sıradan bir
        # klonlama vektörünün yol haritasına klinik araştırma izni
        # koymak, okuyana yanlış bir ölçek duygusu verir.
        approval="klinik" if terapotik else ""))
    yol.phases.append(a5)

    # ---------------------------------------------------------------
    # İzinler, riskler, notlar
    # ---------------------------------------------------------------
    gereken_izin = {t.approval for p in yol.phases for t in p.tasks
                    if t.approval}
    yol.approvals = [APPROVALS[k] for k in
                     ("ibg", "hadyek", "gdo", "klinik")
                     if k in gereken_izin]

    yol.risks = _riskler(findings, part_ids, antibiyotik, insan,
                         gram_poz, disease, organizma)
    yol.readout = (disease.biomarker if disease and
                   getattr(disease, "biomarker", "") else "")

    yol.notes = [
        "**Süreler tahmindir.** Kaynaktan gelmez; tek kişinin "
        "çalıştığı varsayılan, yerleşik akademik iş akışının "
        "alışıldık süreleridir. Ekip büyüdükçe kısalır, ilk denemede "
        "tutmayan adım uzatır. Takvim adımların toplamıdır; "
        "paralelleştirme varsayılmamıştır.",
        "**Basamakların sırası keyfi değildir.** Her aşamanın bir "
        "kapısı var ve o kapı açılmadan sonrakine geçilmez: ifade "
        "gösterilmeden işlev aranmaz, içerme ölçülmeden canlı GDO "
        "hayvana verilmez, etik kurul onayı alınmadan hayvan "
        "çalışması başlamaz.",
        "**Bu bir düzenleyici danışmanlık değildir.** Yasal "
        "dayanaklar künyesiyle yazılıdır ama mevzuat değişir; her "
        "izin için güncel durumu ilgili kurumdan kendiniz teyit "
        "edin.",
    ]
    if not terapotik:
        yol.notes.insert(0, (
            "Bu tasarım terapötik bir bağlamda işaretlenmedi "
            "(konak terapötik suş değil ve hastalık kaydı "
            "seçilmedi), bu yüzden içerme, hayvan modeli ve GDO "
            "basamakları eklenmedi. Terapötik bir yol haritası "
            "istiyorsanız tasarımı Hastalık vektörü sayfasından "
            "başlatın ya da terapötik bir konak seçin."))
    return yol


def _riskler(findings: list[Any] | None, part_ids: list[str],
             antibiyotik: list[str], insan: list[str], gram_poz: bool,
             disease: Any, organizma: str) -> list[Risk]:
    """Denetim bulgularını ve tasarım seçimlerini proje riskine çevir.

    Denetimde görünen bir eksiğin başvuruda risk olarak yazılmaması,
    hakemin ilk yakalayacağı şeydir. Tersi de doğru: bilinen bir
    açığı risk olarak yazıp önlemini söylemek, başvuruyu zayıf değil
    GÜÇLÜ gösterir.
    """
    out: list[Risk] = []

    for f in findings or []:
        onem = getattr(f, "severity", "")
        if onem not in ("hata", "uyari"):
            continue
        out.append(Risk(
            key=getattr(f, "rule_id", ""),
            title=getattr(f, "title", ""),
            etki=getattr(f, "detail", ""),
            onlem=(getattr(f, "fix", "")
                   or "Denetim bu bulgu için bir düzeltme önermiyor; "
                      "önlemi başvuruda kendiniz tanımlayın."),
            kaynak=f"Denetim kuralı {getattr(f, 'rule_id', '')}"))

    if antibiyotik:
        out.append(Risk(
            key="ANTIBIYOTIK",
            title="Antibiyotik direnç markeri in vivo kabul görmez",
            etki=("Tasarımda antibiyotik direnç geni var "
                  f"({', '.join(antibiyotik)}). Terapötik bir "
                  "bakteride direnç geni, bağırsak mikrobiyotasına "
                  "yatay gen aktarımıyla geçebilir; düzenleyici "
                  "çerçevelerde kabul görmez."),
            onlem=("Laboratuvar prototipi bu markerla kurulur, ama "
                   "hayvan çalışmasına geçmeden önce antibiyotiksiz "
                   "idameye (thyA/dapA tamamlaması ya da "
                   "toksin-antitoksin) geçilmelidir. Hastalık "
                   "vektörü sayfasındaki idame merdiveni bu geçişi "
                   "gösterir ve takvime ayrı bir klonlama turu "
                   "olarak YAZILMALIDIR."),
            kaynak="Tasarımdaki marker"))

    if insan:
        out.append(Risk(
            key="INSAN_CDS",
            title="İnsan kaynaklı terapötik yük: kodon ve katlanma",
            etki=(f"{', '.join(insan)} insan genidir. Bakteride kodon "
                  "kullanımı uyumsuzdur ve disülfür bağı taşıyan "
                  "proteinler sitoplazmada yanlış katlanır."),
            onlem=("Dizi analizi sayfasındaki kodon optimizasyonunu "
                   "HEDEF KONAĞIN tablosuyla çalıştırın; disülfür "
                   "bağı varsa periplazmaya ya da hücre dışına "
                   "gönderin. Her iki durumda da çözünen kısımda "
                   "olduğunu jelde göstermek zorunludur."),
            kaynak="Tasarımdaki insan kaynaklı CDS"))

    if gram_poz:
        out.append(Risk(
            key="GRAM_POZ_DONUSUM",
            title=f"{organizma} dönüşüm verimi düşük",
            etki=("Gram-pozitif konaklarda elektroporasyon verimi "
                  "E. coli'nin çok altındadır; yapı kurulsa bile "
                  "suşa sokulamayabilir."),
            onlem=("Duvar zayıflatıcı büyütme (glisin) ve yüksek "
                   "sükrozlu yıkama ile optimizasyon için takvimde "
                   "ayrı süre ayırın; olmazsa konjugasyonla aktarımı "
                   "B planı olarak başvuruya yazın."),
            kaynak="Konak seçimi"))

    if disease is not None and getattr(disease, "gap", ""):
        out.append(Risk(
            key="BILINEN_ACIK",
            title="Tasarımın kayıt defterinde yazılı bilinen açığı",
            etki=disease.gap,
            onlem=("Bu açık tasarımın kendisinden gelir ve "
                   "başvuruda saklanmamalıdır. Hakem zaten görür; "
                   "yazıp önlemini söylemek başvuruyu güçlendirir."),
            kaynak="Hastalık kayıt defteri"))

    return out


# ---------------------------------------------------------------------------
# Çıktı
# ---------------------------------------------------------------------------


def to_markdown(r: Roadmap) -> str:
    """İnsan okuyacak yol haritası."""
    s = [f"# {r.design_name} — gerçekleştirme yol haritası", "",
         f"**Konak:** {r.host_name}  ",
         f"**Toplam:** {r.gun:g} iş günü ≈ {r.ay:g} ay "
         f"(adımların toplamı, paralelleştirme varsayılmadı)", ""]
    if r.readout:
        s += [f"**Başarıyı ölçecek okuma:** {r.readout}", ""]

    for p in r.phases:
        s += [f"## Aşama {p.no} — {p.title}", "",
              f"*{p.amac}*", "",
              f"**Bu aşama {p.gun:g} iş günü.**", ""]
        for t in p.tasks:
            s += [f"### {t.title}  ({t.gun:g} gün)", ""]
            if t.kosul:
                s += [f"> Bu adım şu yüzden var: {t.kosul}", ""]
            s += [f"**Nasıl:** {t.nasil}", "",
                  f"**Gereken:** {t.gereken}", "",
                  f"**Çıktı:** {t.cikti}", ""]
            if t.neden:
                s += [f"**Neden önemli:** {t.neden}", ""]
            if t.approval:
                s += [f"**Önce izin:** {APPROVALS[t.approval].title}", ""]
        if p.kapi:
            s += [f"**Sonraki aşamaya geçiş koşulu:** {p.kapi}", ""]

    if r.approvals:
        s += ["## İzinler ve yasal dayanaklar", ""]
        for a in r.approvals:
            s += [f"### {a.title}", "",
                  f"* **Kime:** {a.body}",
                  f"* **Dayanak:** {a.basis}",
                  f"* **Ne zaman:** {a.when}", ""]
            if a.note:
                s += [a.note, ""]
            if a.caution:
                s += [f"> {a.caution}", ""]

    if r.risks:
        s += ["## Risk yönetimi", "",
              "| Risk | Etki | Önlem |", "|---|---|---|"]
        for k in r.risks:
            s.append(f"| {k.title} | {k.etki.replace(chr(10), ' ')} | "
                     f"{k.onlem.replace(chr(10), ' ')} |")
        s.append("")

    s += ["## Notlar", ""]
    for n in r.notes:
        s.append(f"* {n}")
    return "\n".join(s) + "\n"


def to_proposal(r: Roadmap, disease: Any = None,
                bom: Any = None) -> str:
    """Proje başvurusu iskeleti.

    Başvuru formlarının istediği bölümleri, yol haritasından türeyen
    içerikle doldurur: iş paketleri aşamalardan, takvim sürelerden,
    risk yönetimi denetim bulgularından, bütçe kalemleri (verilirse)
    malzeme listesinden, kaynakça hastalık kaydından.

    Boş bırakılan bölümler **köşeli parantezle** işaretlidir: bir
    başvuru metni, aracın yazdığı yerlerle sizin yazmanız gereken
    yerleri karıştırmamalıdır.
    """
    s = [f"# Proje başvurusu taslağı — {r.design_name}", "",
         "> Bu taslak VectorLab tarafından tasarımdan üretildi. "
         "Köşeli parantezli bölümler sizin yazmanız gerekenlerdir; "
         "geri kalanı tasarımın kendisinden gelir ve üstünde "
         "çalışılması beklenir.", "",
         "## 1. Özet", "",
         "[Projenin bir paragraflık özeti.]", "",
         "## 2. Amaç ve gerekçe", ""]
    if disease is not None:
        s += [f"**Hedef hastalık:** {getattr(disease, 'name', '')}", "",
              getattr(disease, "mechanism", ""), "",
              f"**Terapötik mantık:** {getattr(disease, 'logic', '')}", ""]
    else:
        s += ["[Hangi soruna, neden bu yaklaşımla.]", ""]

    s += ["## 3. Yöntem ve iş paketleri", "",
          "| İP | Başlık | İş günü | Çıktı |", "|---|---|---|---|"]
    for p in r.phases:
        ciktilar = "; ".join(t.cikti for t in p.tasks[:2])
        s.append(f"| İP{p.no} | {p.title} | {p.gun:g} | {ciktilar} |")
    s += [f"| | **Toplam** | **{r.gun:g}** (≈{r.ay:g} ay) | |", ""]

    for p in r.phases:
        s += [f"### İP{p.no} — {p.title}", "", p.amac, ""]
        for t in p.tasks:
            s.append(f"* **{t.title}** ({t.gun:g} gün) — {t.nasil} "
                     f"*Çıktı:* {t.cikti}")
        if p.kapi:
            s += ["", f"*Geçiş koşulu:* {p.kapi}"]
        s.append("")

    s += ["## 4. Başarı ölçütleri", ""]
    if r.readout:
        s += [f"Birincil okuma: **{r.readout}**", ""]
    for p in r.phases:
        if p.kapi:
            s.append(f"* İP{p.no}: {p.kapi}")
    s.append("")

    s += ["## 5. Risk yönetimi ve B planı", ""]
    if r.risks:
        s += ["| Risk | Etki | Önlem |", "|---|---|---|"]
        for k in r.risks:
            s.append(f"| {k.title} | {k.etki.replace(chr(10), ' ')} | "
                     f"{k.onlem.replace(chr(10), ' ')} |")
    else:
        s.append("[Denetim bulgusu çıkmadı; riskleri kendiniz yazın.]")
    s.append("")

    s += ["## 6. Etik ve yasal izinler", ""]
    for a in r.approvals:
        s.append(f"* **{a.title}** — {a.body}. Dayanak: {a.basis}. "
                 f"{a.when}.")
    if not r.approvals:
        s.append("[Bu tasarım için ayrı bir izin çıkmadı.]")
    s.append("")

    s += ["## 7. Bütçe", ""]
    if bom is not None:
        s += ["| Tür | Kalem | Miktar |", "|---|---|---|"]
        for i in getattr(bom, "items", []):
            s.append(f"| {i.kind} | {i.name} | {i.quantity:g} {i.unit} |")
        s += ["", "> Birim fiyatlar kurumunuzun fiyat listesinden "
                  "girilir; VectorLab fiyat üretmez.", ""]
    else:
        s += ["[Sipariş ve denetim izi sayfasından malzeme listesini "
              "indirip buraya ekleyin.]", ""]

    s += ["## 8. Kaynakça", ""]
    if disease is not None and getattr(disease, "papers", None):
        for m in disease.papers:
            s.append(f"* {m}")
    else:
        s.append("[Kaynaklar.]")
    s += ["", "---", "",
          "*Bu taslak bir danışmanlık ya da uygunluk belgesi "
          "değildir. Mevzuat değişir; her izin için güncel durumu "
          "ilgili kurumdan teyit edin.*", ""]
    return "\n".join(s) + "\n"
