"""Yöntem kütüphanesi: numaralı, uygulanabilir tezgâh protokolleri.

Yol haritası "ifadeyi göster" diyor; bu modül **nasıl** gösterileceğini
basamak basamak yazıyor. Hedef, bu işi hiç yapmamış birinin okuyup
uygulayabilmesi.

Üç kural var ve üçü de bu dosyanın tasarımını belirliyor:

1. **Her yöntemin doğrulanmış bir kaynağı var.** Atıflar Europe
   PMC'den PMID ve DOI'siyle alındı; ezberden künye yazılmadı. Bir
   protokolü kaynaksız vermek, onu uygulayacak kişiyi ortada bırakır:
   bir şey tutmadığında bakacağı yer olmaz.

2. **Sayılar BAŞLANGIÇ NOKTASIDIR, emir değil.** Jel yüzdesi,
   voltaj, inkübasyon süresi; hepsi yaygın kullanılan değerlerdir ama
   kesin değerler kaynağın kendisinden ve kullandığınız kitin
   talimatından gelir. Her kaydın `uyari` alanı bunu söyler ve bir
   test o alanın boş olmamasını şart koşar.

3. **Her yöntemin KONTROLÜ ve en sık yapılan HATASI yazılı.**
   Kontrolsüz bir deney sonuç vermez, yorum verir. Yeni başlayan
   birinin düşeceği çukur genellikle protokolün kendisinde değil,
   atladığı kontroldedir - bu yüzden `kontrol` ve `tuzak` alanları
   zorunlu.

Bu protokoller eğitim amaçlıdır ve laboratuvarınızın kendi
onaylanmış prosedürünün yerine geçmez.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Method:
    """Tek bir yöntem: ne için, nasıl, neyle, hangi kontrol."""

    key: str
    title: str
    amac: str
    steps: list[str] = field(default_factory=list)
    malzeme: list[str] = field(default_factory=list)
    sure: str = ""
    kontrol: str = ""
    tuzak: str = ""
    papers: list[str] = field(default_factory=list)
    uyari: str = ("Buradaki sayılar yaygın kullanılan başlangıç "
                  "değerleridir. Kesin değerler kaynağın kendisinden "
                  "ve kullandığınız kitin talimatından gelir; "
                  "laboratuvarınızın onaylanmış prosedürü varsa o "
                  "geçerlidir.")


#: Doğrulanmış atıflar - hepsi Europe PMC'den PMID ve DOI'siyle alındı.
_LAEMMLI = ("Laemmli UK. Cleavage of structural proteins during the "
            "assembly of the head of bacteriophage T4. Nature "
            "227:680-685 (1970). doi:10.1038/227680a0 · PMID 5432063")
_TOWBIN = ("Towbin H, Staehelin T, Gordon J. Electrophoretic transfer "
           "of proteins from polyacrylamide gels to nitrocellulose "
           "sheets. Proceedings of the National Academy of Sciences "
           "76:4350-4354 (1979). doi:10.1073/pnas.76.9.4350 · "
           "PMID 388439")
_BRADFORD = ("Bradford MM. A rapid and sensitive method for the "
             "quantitation of microgram quantities of protein "
             "utilizing the principle of protein-dye binding. "
             "Analytical Biochemistry 72:248-254 (1976). "
             "doi:10.1016/0003-2697(76)90527-3 · PMID 942051")
_MOSMANN = ("Mosmann T. Rapid colorimetric assay for cellular growth "
            "and survival: application to proliferation and "
            "cytotoxicity assays. Journal of Immunological Methods "
            "65:55-63 (1983). doi:10.1016/0022-1759(83)90303-4 · "
            "PMID 6606682")
_HOLO = ("Holo H, Nes IF. High-Frequency Transformation, by "
         "Electroporation, of Lactococcus lactis subsp. cremoris Grown "
         "with Glycine in Osmotically Stabilized Media. Applied and "
         "Environmental Microbiology 55:3119-3123 (1989). "
         "doi:10.1128/aem.55.12.3119-3123.1989 · PMID 16348073")
_SPIZIZEN = ("Anagnostopoulos C, Spizizen J. Requirements for "
             "transformation in Bacillus subtilis. Journal of "
             "Bacteriology 81:741-746 (1961). "
             "doi:10.1128/jb.81.5.741-746.1961 · PMID 16561900")
_INFOGEST = ("Brodkorb A ve ark. INFOGEST static in vitro simulation "
             "of gastrointestinal food digestion. Nature Protocols "
             "14:991-1014 (2019). doi:10.1038/s41596-018-0119-1 · "
             "PMID 30886367")
_STEIDLER = ("Steidler L ve ark. Biological containment of genetically "
             "modified Lactococcus lactis for intestinal delivery of "
             "human interleukin 10. Nature Biotechnology 21:785-789 "
             "(2003). doi:10.1038/nbt840 · PMID 12808464")
_TEER = ("Srinivasan B ve ark. TEER measurement techniques for in "
         "vitro barrier model systems. Journal of Laboratory "
         "Automation 20:107-126 (2015). doi:10.1177/2211068214561025 · "
         "PMID 25586998")
_DROP = ("Herigstad B, Hamilton M, Heersink J. How to optimize the "
         "drop plate method for enumerating bacteria. Journal of "
         "Microbiological Methods 44:121-129 (2001). "
         "doi:10.1016/s0167-7012(00)00241-4 · PMID 11165341")
_ARRIVE = ("Percie du Sert N ve ark. The ARRIVE guidelines 2.0: "
           "updated guidelines for reporting animal research. Journal "
           "of Cerebral Blood Flow and Metabolism 40:1769-1777 (2020). "
           "doi:10.1177/0271678x20943823 · PMID 32663096")


METHODS: dict[str, Method] = {

    # -------------------------------------------------------------- SDS-PAGE
    "sdspage": Method(
        key="sdspage",
        title="SDS-PAGE ile protein ifadesini görme",
        amac=("Hücrenin beklediğiniz proteini üretip üretmediğini "
              "görmek. Jelde beklenen boyda, indüklenmiş örnekte olup "
              "indüklenmemişte olmayan bir bant arıyorsunuz."),
        steps=[
            "Bir gece önce kültürü kur. Ertesi gün taze besiyerine "
            "1/100 seyrelt ve OD600 yaklaşık 0,5-0,6 olana kadar "
            "37 °C'de çalkalayarak büyüt (tipik olarak 2-3 saat).",
            "**İndüksiyondan HEMEN ÖNCE 1 mL örnek al.** Bu senin "
            "'0. saat' örneğin; bunu almayı unutmak, jelde "
            "karşılaştıracak bir şeyin olmaması demektir.",
            "Kültürü ikiye böl. Birine indükleyiciyi ekle (tasarımın "
            "promotörüne göre IPTG ya da arabinoz; derişimi "
            "Laboratuvar sayfası tasarımından hesaplıyor), ötekine "
            "hiçbir şey ekleme.",
            "3-4 saat daha büyüt. Her saat başı iki koldan da 1 mL "
            "örnek al.",
            "Her örneği santrifüjle (yaklaşık 13.000 g, 1 dakika), "
            "üstteki sıvıyı at, çökeltiyi sakla. Örnekleri -20 °C'de "
            "biriktirebilirsin.",
            "**Hücre sayısını eşitle.** Her çökeltiyi, alındığı "
            "andaki OD600 ile orantılı hacimde yükleme tamponunda "
            "çöz - örneğin OD 0,5 için 50 µL, OD 1,0 için 100 µL. "
            "Eşitlemezsen jeldeki fark ifadeden değil hücre "
            "sayısından gelir.",
            "95 °C'de 5 dakika kaynat, kısa santrifüjle.",
            "Jele yükle. Protein boyuna göre jel yüzdesi seç: 10-20 "
            "kDa için %15, 20-60 kDa için %12, 60 kDa üstü için %10 "
            "yaygın seçimlerdir. Bir kuyuya mutlaka **boy belirteci** "
            "yükle.",
            "Yığma jelde düşük voltajda (yaklaşık 80-100 V), ayırma "
            "jeline geçince daha yüksek voltajda (yaklaşık 120-150 V) "
            "yürüt; boya cephesi jelin altına yaklaşınca durdur.",
            "Coomassie ile boya, sonra renk açıcıyla arka planı aç. "
            "Bant zayıfsa gümüş boyama daha duyarlıdır.",
        ],
        malzeme=["Dikey jel düzeneği ve güç kaynağı",
                 "Akrilamid/bisakrilamid çözeltisi, APS, TEMED",
                 "Yükleme tamponu (SDS + indirgeyici ajan içeren)",
                 "Protein boy belirteci", "Coomassie boyası",
                 "Masa üstü santrifüj, ısı bloğu"],
        sure="Kültür dahil 2 gün; jelin kendisi yarım gün",
        kontrol=("**İndüklenmemiş kol zorunludur.** Yalnızca "
                 "indüklenmiş örneği yürütmek hiçbir şey söylemez: "
                 "o bant zaten orada olabilir. Ayrıca boş vektör "
                 "taşıyan bir suş yürütmek, bandın kasetten geldiğini "
                 "kanıtlar."),
        tuzak=("Bandı görüp 'oldu' demek. Protein **çözünür** mü, "
               "yoksa inklüzyon cisimciğinde mi? Hücreyi parçalayıp "
               "santrifüjleyin ve çözünen (üst sıvı) ile çökelti "
               "kısmını AYRI kuyulara yükleyin. İnklüzyon "
               "cisimciğindeki protein jelde bol görünür ama "
               "işlevsizdir - bu, yeni başlayanın en sık yanıldığı "
               "yerdir."),
        papers=[_LAEMMLI, _BRADFORD]),

    # ---------------------------------------------------------------- Western
    "western": Method(
        key="western",
        title="Western blot ile bandın kimliğini doğrulama",
        amac=("SDS-PAGE'de gördüğünüz bandın gerçekten aradığınız "
              "protein olduğunu kanıtlamak. Doğru boyda bir bant tek "
              "başına kanıt değildir; o boyda başka proteinler de "
              "vardır."),
        steps=[
            "Önce normal bir SDS-PAGE yürüt (yukarıdaki yöntem), ama "
            "jeli boyama.",
            "Jeli ve membranı (nitroselüloz ya da PVDF) transfer "
            "tamponunda ıslat. PVDF kullanıyorsan önce metanolle "
            "etkinleştir.",
            "Transfer sandviçini kur: sünger - kâğıt - **jel** - "
            "**membran** - kâğıt - sünger. Jel ile membran arasında "
            "hava kabarcığı bırakma; kabarcık olan yerde transfer "
            "olmaz ve jelde beyaz leke olarak görünür.",
            "Membranın negatif kutba değil, jelin karşısına geldiğinden "
            "emin ol: protein negatif yüklüdür ve membrana doğru "
            "göçer. Ters kurarsan protein tampona kaçar.",
            "Transfer et (yaş transferde yaklaşık 100 V, 1 saat, "
            "buzda; yarı kuru düzeneklerde süre daha kısadır).",
            "Membranı bloklama çözeltisinde (yaygın olarak %5 yağsız "
            "süt tozu ya da BSA) oda sıcaklığında 1 saat çalkala.",
            "Birincil antikorla inkübe et (çoğu protokol +4 °C'de bir "
            "gece önerir). Seyreltme oranı antikorun kendi "
            "talimatındadır.",
            "Yıka (yaklaşık 3 x 5 dakika), işaretli ikincil antikorla "
            "1 saat inkübe et, tekrar yıka.",
            "Görüntüle (kemilüminesans ya da renk oluşturan substrat) "
            "ve boy belirtecinin yerine göre bandı konumlandır.",
        ],
        malzeme=["Transfer düzeneği", "Nitroselüloz ya da PVDF membran",
                 "Bloklama çözeltisi", "Birincil ve işaretli ikincil "
                 "antikor", "Görüntüleme sistemi"],
        sure="2 gün (bir gece antikor inkübasyonu dahil)",
        kontrol=("Birincil antikoru atlayan bir membran şeridi "
                 "yürütün: orada bant çıkıyorsa gördüğünüz şey ikincil "
                 "antikorun özgül olmayan bağlanmasıdır. Ayrıca boş "
                 "vektörlü suşun örneği negatif kontroldür."),
        tuzak=("Antikor bulamamak. Terapötik proteininize karşı ticari "
               "antikor yoksa Western yapamazsınız. Çözüm tasarım "
               "aşamasındadır: kasete bir etiket (His, FLAG) ekleyip "
               "etikete karşı antikor kullanmak. Bunu proteinin "
               "işlevini bozmayacak uçtan yapmak gerekir - Tasarımcı "
               "sayfasındaki çerçeve denetimi etiketin okuma "
               "çerçevesinde olup olmadığını söyler."),
        papers=[_TOWBIN]),

    # -------------------------------------------------------------- salgılanma
    "salgi": Method(
        key="salgi",
        title="Salgılanmayı doğrulama (lizis kontrolüyle)",
        amac=("Proteinin hücre dışına GERÇEKTEN salgılandığını "
              "göstermek. Süpernatantta protein bulmak tek başına "
              "salgılanma kanıtı değildir: hücreler patlamış da "
              "olabilir."),
        steps=[
            "Kültürü büyüt ve indükle (SDS-PAGE yöntemindeki gibi).",
            "Kültürü santrifüjle (yaklaşık 4.000-5.000 g, 10 dakika, "
            "+4 °C). Üst sıvı = **süpernatant**, çökelti = **hücre**.",
            "Süpernatantı 0,22 µm filtreden geçir: kalan hücreler "
            "sonucu bozar.",
            "Süpernatanttaki protein seyreltiktir; yoğunlaştır "
            "(trikloroasetik asitle çöktürme ya da santrifüj "
            "filtresiyle). Yoğunlaştırma oranını yaz, çünkü jeldeki "
            "yoğunluğu yorumlarken gerekecek.",
            "Hücre çökeltisini tamponla çöz ve parçala (sonikasyon ya "
            "da lizozim). Bu senin **hücre içi** örneğin.",
            "**Kritik adım:** aynı süpernatantta hücre içinde kalması "
            "gereken bir proteini ara. Sitoplazmik bir enzim (örneğin "
            "bir dehidrogenaz) ya da etiketli bir sitoplazmik "
            "kontrol kullan.",
            "Süpernatant ve hücre içi örneklerini **aynı jele** yan "
            "yana yükle, SDS-PAGE ve gerekirse Western yürüt.",
            "Yorumla: terapötik protein süpernatantta VAR + "
            "sitoplazmik kontrol süpernatantta YOK ⇒ salgılanma. "
            "İkisi de varsa gördüğünüz şey lizistir.",
        ],
        malzeme=["Soğutmalı santrifüj", "0,22 µm filtre",
                 "Protein yoğunlaştırma aracı",
                 "Sitoplazmik kontrol için antikor ya da enzim tayini"],
        sure="2-3 gün",
        kontrol=("Sitoplazmik lizis göstergesi. Bu kontrol olmadan "
                 "yapılan salgılanma iddiası, hakem tarafından "
                 "sorulacak ilk şeydir."),
        tuzak=("Geç faz kültürle çalışmak. Kültür ne kadar uzun "
               "büyürse o kadar çok hücre ölür ve süpernatanta "
               "sitoplazmik protein sızar. Örneği erken alın ve "
               "lizis kontrolünü her zaman birlikte yürütün."),
        papers=[_LAEMMLI, _TOWBIN]),

    # ---------------------------------------------------------------- büyüme
    "buyume": Method(
        key="buyume",
        title="Büyüme eğrisi ile metabolik yükü ölçme",
        amac=("Kasetin hücreye ne kadar yük bindirdiğini sayıyla "
              "göstermek. Ağır bir kaset taşıyan hücre, bağırsakta "
              "taşımayan rakipleri tarafından birkaç nesilde elenir."),
        steps=[
            "Üç suş hazırla: (a) plazmitsiz konak, (b) boş vektör "
            "taşıyan konak, (c) terapötik kaseti taşıyan konak.",
            "Her birinden bir gece kültürü kur.",
            "Ertesi gün üçünü de aynı başlangıç OD600'üne (tipik "
            "olarak 0,05) taze besiyerine seyrelt. **Aynı besiyeri, "
            "aynı sıcaklık, aynı çalkalama hızı** - yoksa "
            "karşılaştırma geçersizdir.",
            "Her suş için en az üç biyolojik tekrar kur (ayrı "
            "kolonilerden gelen ayrı kültürler; aynı kültürden üç "
            "kuyu teknik tekrardır ve yerine geçmez).",
            "OD600'ü düzenli aralıklarla ölç (plaka okuyucuda 15-30 "
            "dakikada bir otomatik; elle ölçüyorsan saat başı), "
            "en az 8-10 saat.",
            "Üssel fazda eğrinin doğrusal kısmından **üreme hızını** "
            "hesapla ve ikilenme süresine çevir.",
            "Üç şeyi ayrı ayrı raporla: gecikme (lag) süresi, "
            "ikilenme süresi, ulaşılan son yoğunluk. Bunlar farklı "
            "şeyler söyler - kaset yalnızca gecikmeyi uzatıyorsa "
            "sorun farklıdır, son yoğunluğu düşürüyorsa farklı.",
        ],
        malzeme=["Plaka okuyucu ya da spektrofotometre",
                 "Çalkalamalı inkübatör", "96 kuyulu plaka ya da şişe"],
        sure="2 gün",
        kontrol=("Plazmitsiz konak ve boş vektör birlikte. Boş vektör "
                 "olmadan, gördüğünüz yükün plazmitin kendisinden mi "
                 "kasetten mi geldiğini ayıramazsınız."),
        tuzak=("Tek kuyu okumak ve eğriyi gözle yorumlamak. Üreme "
               "hızı hesaplanmadan 'biraz yavaş büyüyor' demek "
               "başvuruda veri sayılmaz."),
        papers=[]),

    # ------------------------------------------------------------ kararlılık
    "kararlilik": Method(
        key="kararlilik",
        title="Plazmit kararlılığı (antibiyotiksiz seri pasaj)",
        amac=("Antibiyotik YOKKEN plazmitin tutulup tutulmadığını "
              "ölçmek. Hastaya antibiyotik verilmeyecek, dolayısıyla "
              "asıl soru budur."),
        steps=[
            "Tek koloniden, seçilimli besiyerinde bir gece kültürü "
            "kur (başlangıçta plazmitin hepsinde olduğundan emin ol).",
            "**Antibiyotiksiz** taze besiyerine 1/1000 seyrelt. Her "
            "1/1000 seyreltme yaklaşık 10 nesildir.",
            "Her gün aynı seyreltmeyi tekrarla. 50-100 nesil hedefle "
            "(yani 5-10 pasaj).",
            "Her pasajda kültürden bir örnek al, uygun seyreltmeyi "
            "yap ve **iki plağa** ek: biri antibiyotiksiz (toplam "
            "hücre), biri antibiyotikli (plazmit taşıyan hücre).",
            "Koloni sayımı için damla plak yöntemi, yayma plağa göre "
            "daha az iş ve daha az malzeme ister.",
            "Her pasaj için yüzdeyi hesapla: (antibiyotikli plak "
            "kolonisi / antibiyotiksiz plak kolonisi) x 100.",
            "Nesil sayısına karşı yüzdeyi grafiğe dök. Düşüş varsa "
            "eğimden nesil başına kayıp oranını hesapla.",
            "İçerme markeri (thyA, dapA) kullanıyorsan antibiyotik "
            "yerine o maddenin yokluğunu seçilim olarak kullan.",
        ],
        malzeme=["Seçilimli ve seçilimsiz plaklar",
                 "Seri seyreltme malzemesi", "İnkübatör"],
        sure="10-15 gün (pasaj sayısına göre)",
        kontrol=("Başlangıç (0. nesil) ölçümü. %100 ile başlamıyorsa "
                 "kültürünüz zaten karışıktır ve geri kalan veri "
                 "yorumlanamaz."),
        tuzak=("Yalnızca antibiyotikli plağa ekmek. O plakta koloni "
               "görmek 'plazmit duruyor' demek değildir; kaç hücrenin "
               "kaybettiğini ancak iki plağın ORANI söyler."),
        papers=[_DROP]),

    # --------------------------------------------------------------- sindirim
    "sindirim": Method(
        key="sindirim",
        title="Benzetilmiş mide-bağırsak koşullarına dayanım",
        amac=("Ağızdan verilen bakterinin mideyi ve safrayı geçip "
              "bağırsağa canlı ulaşıp ulaşmadığını ölçmek. Geçemezse "
              "kasetin ne kadar iyi olduğunun önemi kalmaz."),
        steps=[
            "Protokolü INFOGEST uzlaşı yönteminden al: ağız, mide ve "
            "ince bağırsak aşamalarının tampon bileşimi, enzim "
            "etkinlikleri ve pH değerleri orada tanımlıdır. Kendi "
            "tamponunuzu uydurmayın - sonuç başka çalışmalarla "
            "karşılaştırılamaz hale gelir.",
            "Suşu durağan faza kadar büyüt, santrifüjle topla ve "
            "tamponla yıka.",
            "Başlangıç canlı hücre sayısını belirle (CFU/mL). Bu "
            "senin referans noktan.",
            "Mide aşamasını uygula: örneği benzetilmiş mide sıvısına "
            "koy, 37 °C'de çalkalayarak INFOGEST'in belirttiği süre "
            "boyunca bekle.",
            "Belirli zaman noktalarında (örneğin 0, 30, 60, 120 "
            "dakika) örnek al, hemen seyrelt ve plağa ek.",
            "Bağırsak aşamasına geç: safra tuzu ve pankreatin içeren "
            "benzetilmiş bağırsak sıvısı, INFOGEST'teki pH'da.",
            "Aynı şekilde zaman noktalarında örnekle ve say.",
            "Sonucu **log azalma** olarak ver: log10(başlangıç CFU) - "
            "log10(son CFU). 'Hayatta kaldı' demek yerine sayı verin.",
        ],
        malzeme=["pH metre", "Benzetilmiş sıvı bileşenleri "
                 "(pepsin, safra tuzu, pankreatin)",
                 "37 °C çalkalamalı su banyosu", "Plaklar"],
        sure="8-10 gün (tekrarlarla birlikte)",
        kontrol=("Bilinen bir probiyotik suşu aynı koşuldan geçirin: "
                 "o da ölüyorsa sorun sizin suşunuzda değil "
                 "düzeneğinizdedir."),
        tuzak=("Örneği seyreltmeden plağa ekmek. Mide aşamasından "
               "çıkan örnek hâlâ asidiktir ve plakta büyümeyi "
               "baskılar; hemen nötr tamponla seyreltin."),
        papers=[_INFOGEST, _DROP]),

    # ----------------------------------------------------------------- içerme
    "icerme": Method(
        key="icerme",
        title="Biyolojik içerme: kaçak oranını ölçme",
        amac=("İçerme mekanizmasının ne kadar sızdırdığını sayıyla "
              "vermek. Mekanizmanın VARLIĞI değil, KAÇAK ORANI "
              "raporlanır; düzenleyici çerçeveye girilen yer burasıdır."),
        steps=[
            "Suşu izin verici koşulda büyüt (timin ya da DAP içeren "
            "besiyeri) ve durağan faza getir.",
            "Hücreleri topla, izin verici maddeyi İYİCE uzaklaştır: "
            "en az üç kez yıka. Yıkamadan kalan eser madde, kaçak "
            "gibi görünen sahte koloni üretir.",
            "Başlangıç canlı hücre sayısını belirle (CFU/mL).",
            "**Çok yüksek hücre sayısı** ek: kaçak oranı küçük bir "
            "sayı olduğu için az hücre ekerseniz hiç koloni görmez "
            "ve yanlışlıkla 'kaçak yok' dersiniz. 10⁹-10¹⁰ hücre "
            "ekmek yaygın bir yaklaşımdır.",
            "Eksiltilmiş plakları uzun süre inkübe et (48-72 saat): "
            "kaçan hücreler yavaş büyür ve bir gecede görünmeyebilir.",
            "Kaçak oranını hesapla: (eksiltilmiş plaktaki koloni) / "
            "(ekilen toplam canlı hücre). Sonucu bilimsel gösterimle "
            "ver.",
            "**Kaçan kolonileri topla ve incele.** Nasıl kaçtılar - "
            "besiyerinden bulaşma mı, kromozomal geri dönüş mü, "
            "plazmitin kromozoma girmesi mi? Bunu dizileme ile "
            "belirleyin; düzenleyici dosyanın sorduğu soru tam "
            "olarak budur.",
        ],
        malzeme=["Eksiltilmiş besiyeri (timinsiz ya da DAP'sız)",
                 "Yıkama tamponu", "Büyük plaklar", "İnkübatör"],
        sure="12-15 gün",
        kontrol=("İzin verici plak: orada koloni çıkmıyorsa "
                 "hücreleriniz zaten ölmüştür ve kaçak ölçümü "
                 "anlamsızdır. Ayrıca boş plak (steril kontrol)."),
        tuzak=("Az hücre ekip 'kaçak yok' demek. Kaçak oranı ancak "
               "ekilen hücre sayısı kadar hassas ölçülebilir; "
               "raporda 'kaçak < 1/ekilen hücre sayısı' diye "
               "yazılır, '0' diye değil."),
        papers=[_STEIDLER, _DROP]),

    # ------------------------------------------------------------------ TEER
    "teer": Method(
        key="teer",
        title="Bağırsak epiteli modeli ve bariyer bütünlüğü (TEER)",
        amac=("Suşun bağırsak epiteline ne yaptığını hayvandan önce "
              "görmek: bariyeri koruyor mu, bozuyor mu, iltihabı "
              "azaltıyor mu."),
        steps=[
            "Bağırsak epiteli hücre hattını geçirgen zarlı kültür "
            "kaplarına (Transwell tipi) ek ve **farklılaşana kadar** "
            "büyüt. Caco-2 için bu tipik olarak 21 gün sürer ve "
            "kısaltılamaz: farklılaşmamış hücrede bariyer yoktur.",
            "Her 2-3 günde bir besiyerini değiştir.",
            "TEER'i düzenli ölç. Ölçümden önce plakayı oda "
            "sıcaklığına getir: TEER sıcaklığa duyarlıdır ve "
            "inkübatörden çıkar çıkmaz ölçmek sahte düşüş gösterir.",
            "Değerleri boş (hücresiz) kuyunun direncini çıkararak ve "
            "zar alanıyla çarparak düzelt; sonucu Ω·cm² olarak ver.",
            "TEER plato çizdiğinde tabaka hazırdır. Bakteriyi üst "
            "bölmeye ekle (antibiyotiksiz besiyerinde).",
            "Gruplar: yalnız hücre, hücre + boş vektörlü suş, hücre + "
            "terapötik suş, ve varsa hücre + iltihap uyarıcı.",
            "Belirli saatlerde TEER'i tekrar ölç ve alt bölmeden "
            "besiyeri alarak salgılanan sitokinleri ELISA ile ölç.",
            "Ko-kültür süresini kısa tut (tipik olarak 4-8 saat): "
            "bakteri uzun sürede hücreleri besiyeri rekabetiyle "
            "öldürür ve bu, terapötik etki sanılabilir.",
        ],
        malzeme=["Hücre kültürü altyapısı ve CO2 inkübatörü",
                 "Geçirgen zarlı kültür kapları", "TEER ölçer (voltohmmetre)",
                 "ELISA kitleri"],
        sure="25-30 gün (farklılaşma dahil)",
        kontrol=("Hücresiz kuyu (zar direnci) ve boş vektörlü suş. "
                 "Ayrıca bakterisiz hücre, zaman içinde TEER'in "
                 "kendiliğinden nasıl değiştiğini gösterir."),
        tuzak=("Farklılaşmayı beklememek. 7 günlük Caco-2 tabakası "
               "bariyer kurmamıştır; ölçtüğünüz TEER değeri gerçek "
               "bağırsak epitelini temsil etmez."),
        papers=[_TEER]),

    # ---------------------------------------------------------- sitotoksisite
    "mtt": Method(
        key="mtt",
        title="Sitotoksisite (MTT/XTT ile hücre canlılığı)",
        amac=("Terapötik olduğu varsayılan şeyin hücreye zarar verip "
              "vermediğini ölçmek. Yüksek derişimde toksik çıkması "
              "sık görülür."),
        steps=[
            "Hücreleri 96 kuyulu plakaya, kuyu başına eşit sayıda ek "
            "ve bir gece yapışmalarını bekle.",
            "Test edeceğin şeyin seri seyreltmesini hazırla "
            "(süpernatant, saflaştırılmış protein ya da bakteri "
            "süspansiyonu). En az 5-6 derişim basamağı kullan - tek "
            "derişim doz-yanıt eğrisi vermez.",
            "Belirlenen süre inkübe et (24 ya da 48 saat yaygındır).",
            "MTT çözeltisini ekle ve 2-4 saat inkübe et. Canlı "
            "hücrelerin enzimleri sarı tetrazolyum tuzunu mor "
            "formazan kristallerine çevirir.",
            "Besiyerini dikkatle çek (kristalleri emme!) ve "
            "çözücü ekleyerek kristalleri çöz.",
            "570 nm civarında absorbansı oku.",
            "Canlılığı işlenmemiş kontrole göre yüzde olarak hesapla "
            "ve doz-yanıt eğrisini çiz.",
        ],
        malzeme=["96 kuyulu plaka", "MTT ya da XTT kiti",
                 "Plaka okuyucu", "Çözücü (DMSO vb.)"],
        sure="4-5 gün",
        kontrol=("İşlenmemiş hücre (%100 canlılık) ve bilinen bir "
                 "toksik madde (pozitif kontrol). Pozitif kontrol "
                 "olmadan 'toksik değil' sonucunu, testin "
                 "çalışmadığından ayıramazsınız."),
        tuzak=("Bakteri süspansiyonuyla çalışırken bakterinin kendi "
               "metabolizmasının da MTT'yi indirgemesi. Bakteri "
               "içeren kuyularda bakteri-tek kontrolü mutlaka "
               "yürütün, yoksa canlılığı olduğundan yüksek "
               "ölçersiniz."),
        papers=[_MOSMANN]),

    # ---------------------------------------------------------- elektroporasyon
    "elektroporasyon": Method(
        key="elektroporasyon",
        title="Gram-pozitif konağa elektroporasyonla aktarım",
        amac=("Yapıyı E. coli'den gıda sınıfı ya da endüstriyel "
              "Gram-pozitif konağa sokmak. Kalın hücre duvarı yüzünden "
              "verim E. coli'nin çok altındadır."),
        steps=[
            "Konağı **glisin içeren** besiyerinde büyüt. Glisin hücre "
            "duvarının çapraz bağlanmasını zayıflatır ve dönüşüm "
            "verimini belirgin artırır; kullanılacak derişim suşa "
            "göre değişir, kaynaktaki aralıktan başlayıp kendi "
            "suşunuz için en iyisini bulun.",
            "Üssel fazın ortasında hücreleri topla (soğutmalı "
            "santrifüj, +4 °C). Bu noktadan sonra her şey buzda.",
            "Ozmotik olarak dengeleyen bir tamponla (yüksek sükroz "
            "içeren) en az üç kez yıka. Tuzu tamamen uzaklaştırmak "
            "şarttır: kalan tuz elektroporasyonda ark oluşturur ve "
            "örneği yakar.",
            "Hücreleri küçük hacimde çöz ve soğuk küvetlere dağıt. "
            "Hemen kullanmayacaksan -80 °C'de sakla.",
            "Plazmit DNA'sını ekle (tuzsuz, saf; etanol çöktürmesi "
            "yapılmış olması iyidir), karıştır, 1-2 dakika buzda "
            "bekle.",
            "Elektroporasyonu uygula. Darbe parametreleri cihaza ve "
            "küvet aralığına bağlıdır; kaynaktaki değerlerden "
            "başlayın. **Zaman sabiti (time constant) cihazda "
            "görünür** - beklenenden çok düşükse örnekte tuz var "
            "demektir.",
            "Hemen ozmotik olarak dengeleyen, besleyici katkılı "
            "besiyeri ekle ve seçilimsiz olarak 1-3 saat dinlendir. "
            "Bu adım atlanırsa direnç geni ifade edilmeye fırsat "
            "bulamaz ve koloni çıkmaz.",
            "Seçilimli plağa ek ve konağa uygun sıcaklıkta inkübe et. "
            "Koloniler E. coli'den yavaş çıkabilir, 48 saat bekle.",
        ],
        malzeme=["Elektroporatör ve soğuk küvetler",
                 "Glisin", "Yüksek sükrozlu yıkama tamponu",
                 "Soğutmalı santrifüj", "Saf, tuzsuz plazmit DNA'sı"],
        sure="4-5 gün (yetkin hücre hazırlığı dahil)",
        kontrol=("DNA'sız küvet (sahte dönüşüm): orada koloni "
                 "çıkıyorsa seçiliminiz sızdırıyordur. Ayrıca o "
                 "konakta çalıştığı bilinen bir plazmitle pozitif "
                 "kontrol - koloni çıkmazsa sorun yapınızda değil "
                 "yetkin hücrenizdedir."),
        tuzak=("Ark olunca denemeyi sürdürmek. Ark, tuzun "
               "kalıntısıdır ve o örnekten koloni çıkmaz; yıkamayı "
               "tekrarlayıp baştan başlamak zaman kazandırır. "
               "Verim yine düşükse konjugasyonla aktarımı B planı "
               "olarak düşünün."),
        papers=[_HOLO, _SPIZIZEN]),

    # -------------------------------------------------------------- in vivo
    "gavaj": Method(
        key="gavaj",
        title="Hayvanda uygulama ve kolonizasyon izleme",
        amac=("Suşun bağırsakta ne kadar kaldığını ve terapötik "
              "etkiyi ölçmek. Yol haritasının en pahalı ve en çok "
              "izin isteyen basamağı."),
        steps=[
            "**Önce HADYEK onayını al.** Onaysız başlanan çalışma "
            "yayımlanamaz. Başvuruya deney düzenini, hayvan sayısını "
            "ve güç analizini ekleyin.",
            "Deney düzenini ARRIVE 2.0 kılavuzuna göre önceden yazın: "
            "gruplar, randomizasyon, körleştirme, dışlama ölçütleri "
            "ve birincil sonlanım noktası. Bunları deneyden SONRA "
            "seçmek sonucu yayımlanamaz kılar.",
            "Grupları kur. En az: sağlıklı kontrol, hastalık + taşıyıcı "
            "(bakterisiz), hastalık + **boş vektörlü** suş, hastalık + "
            "terapötik suş. Boş vektör grubu olmadan etkinin "
            "bakteriden mi kasetten mi geldiği ayrılamaz.",
            "Hayvanları rastgele gruplara ayır ve kafesleri karıştır: "
            "aynı kafesteki hayvanlar birbirinin dışkısını yer "
            "(koprofaji) ve bu, gruplar arası bulaşma demektir.",
            "Suşu taze hazırla, uygulamadan hemen önce CFU'sunu "
            "doğrula - verdiğiniz dozu bilmeden sonuç yorumlanamaz.",
            "Gavajla ya da içme suyuyla uygula. Doz ve sıklık "
            "kolonizasyon verinize göre belirlenir; Nissle gibi "
            "suşlar kalıcı kolonize etmediği için tekrar gerekir.",
            "Dışkıyı düzenli topla ve iki yöntemle izle: seçici "
            "besiyerinde CFU sayımı (canlı hücre) ve plazmite özgü "
            "qPCR (toplam hedef DNA). İkisi birlikte, suşun canlı mı "
            "yoksa yalnızca DNA olarak mı bulunduğunu ayırır.",
            "Birincil sonlanım noktasını körleştirilerek ölçün.",
            "Dışkı ve altlık örneklerinde çevresel kaçışı ve "
            "mikrobiyotaya yatay gen aktarımını da örnekleyin - bu "
            "veri sonradan toplanamaz.",
        ],
        malzeme=["Etik kurul onayı", "Gavaj iğnesi",
                 "Seçici besiyeri", "qPCR düzeneği ve özgül primerler",
                 "Hastalığa özgü ölçüm altyapısı"],
        sure="60-80 gün (onay süresi hariç)",
        kontrol=("Boş vektörlü suş grubu ve bakterisiz taşıyıcı "
                 "grubu. İkisi de olmadan sonuç yorumlanamaz."),
        tuzak=("Hayvan sayısını güç analizi olmadan yazmak. HADYEK'in "
               "en sık geri çevirme sebebi budur ve haklıdır: az "
               "hayvanla yapılan deney sonuç vermez, fazlasıyla "
               "yapılan deney gereksiz hayvan kullanır."),
        papers=[_ARRIVE, _DROP]),
}


def get(key: str) -> Method | None:
    return METHODS.get(key)


def to_markdown(m: Method) -> str:
    """Tek bir yöntemi okunur protokol olarak yaz."""
    s = [f"### {m.title}", "", f"*{m.amac}*", "",
         f"**Süre:** {m.sure}", "", "**Basamaklar**", ""]
    for i, adim in enumerate(m.steps, start=1):
        s.append(f"{i}. {adim}")
    s += ["", "**Gereken malzeme**", ""]
    for x in m.malzeme:
        s.append(f"* {x}")
    s += ["", f"**Şart olan kontrol:** {m.kontrol}", "",
          f"**En sık yapılan hata:** {m.tuzak}", ""]
    if m.papers:
        s += ["**Kaynak**", ""]
        for p in m.papers:
            s.append(f"* {p}")
        s.append("")
    s += [f"> {m.uyari}", ""]
    return "\n".join(s)
