# VectorLab

Bakteriyel vektör tasarım, analiz ve literatür aracı. Addgene + SnapGene
ikilisinin ders düzeyinde yerine geçecek, tek dosyalık veritabanıyla
çalışan yerel bir uygulama.

**Kapsam bilinçli olarak sadece bakteriyel sistemlerdir** (ağırlıklı
*E. coli*). Memeli, bitki, böcek, maya **vektör elemanı** — orijin,
promotör, marker, iskelet — yoktur ve eklenmeyecektir.

İstisna kargodadır, elemanda değil: Hastalık vektörü sayfasında üç
insan geni vardır — IL-10 (IBD, Steidler 2000), elafin (IBD, Motta
2012) ve GLP-1 (diyabet, Duan 2015) — çünkü bu üç tasarımda bakterinin
*salgıladığı* terapötik protein insan proteinidir ve kaynak makaleler
de bunları kullanır. Yani bakteri bir insan proteini üretir; vektörün
kendisi (orijin, promotör, marker, iskelet) baştan sona bakteriyeldir
ve insan geni hiçbir yerde vektör elemanı olarak kullanılmaz.

## Masaüstü uygulaması olarak çalıştırma

Masaüstündeki **VectorLab** kısayolu uygulamayı kendi penceresinde
açar: adres çubuğu yok, sekme yok, tarayıcı menüsü yok. Tarayıcı
açılmaz.

```bash
pythonw desktop.py
```

Pencere Windows'un **WebView2** bileşeni üzerine kuruludur. Bu, Edge
tarayıcısının kendisi değil; uygulamaların içine gömülen görüntüleme
motorudur (Teams ve VS Code da aynısını kullanır) ve Windows 11'de
hazır gelir. Doğrulandı: uygulama açıldığında `msedge.exe` ve
`chrome.exe` pencere sayısı değişmiyor, açılan pencerenin sınıfı
`WindowsForms10.Window` (tarayıcı olsaydı `Chrome_WidgetWin_1`
olurdu).

Üç ayrıntı bilinçli:

- **Port dinamik alınıyor.** Sabit 8501 kullanılsaydı uygulama ikinci
  kez açıldığında ya da başka bir Streamlit çalışırken çakışır ve
  pencere boş gelirdi.
- **Kapanış iki yönlü.** Pencereyi kapatınca sunucu süreç ağacıyla
  birlikte kapanır; sunucu önce çökerse pencere de kapanır, boş ekranla
  kalmazsınız.
- **Konsol hiç açılmaz.** `pythonw` kullanıldığı için pencere
  parlaması olmaz; sorun çıkarsa yerel bir hata kutusu belirir ve
  ayrıntı `data/desktop.log` içine yazılır.

Hata ayıklarken `baslat-konsollu.bat` günlüğü konsolda da gösterir.

## Kurulum ve çalıştırma

```bash
pip install -r requirements.txt
```

```bash
streamlit run streamlit_app.py
```

Tarayıcıda `http://localhost:8501` açılır. Veritabanı ilk açılışta
`data/vectorlab.db` olarak kendiliğinden kurulur ve 122 parça + 8 konak
suşuyla doldurulur (60'ı kendi çekirdek kataloğu, 62'si iGEM
Registry'den). Ayrıca 15 referans vektörün tam dizisi depoda hazır
gelir; uygulama çalışırken hiçbir ağ isteği yapmaz (tek istisna
Literatür sayfasıdır).

## Sayfalar

### Sihirbaz
Birkaç soru sorar (amaç, toksisite, gen kaynağı, saflaştırma, disülfit,
ikinci plazmit, ölçek) ve çalışan bir vektör önerir. Önceki yanıtlara
göre anlamsız kalan sorular gizlenir — sadece klonlama seçilirse 8 soru
yerine 3 soru sorulur.

**Her parça seçiminin gerekçesi yazılır.** Motor "T7lac kullan" demez,
"T7lac kullan çünkü …" der; bu, ders bağlamında kararın kendisinden daha
öğreticidir. Öneri Tasarımcı sayfasına aktarılabilir, orada istediğiniz
parçayı değiştirebilirsiniz.

Karar zincirinden örnekler:

| Yanıt | Sonuç | Neden |
|---|---|---|
| Ürün toksik olabilir | T7lac yerine P-BAD + AraC, konak TOP10 | Arabinoz dozu kademeli ayarlanır; T7 ya kapalı ya tam açıktır. TOP10 arabinozu tüketmez, indüksiyon sabit kalır |
| Gen ökaryot kaynaklı | Konak Rosetta(DE3) | pRARE nadir kodonların tRNA'larını sağlar |
| Disülfit bağı var | Konak SHuffle T7 | Sitoplazma oksitleyici, bağlar periplazma olmadan kurulur |
| Çözünürlük sorunu | His yerine His-SUMO | Çözünürlüğü artırır ve kesimde artık bırakmaz |
| İkinci plazmit var | p15A orijin | ColE1 ailesiyle uyumlu |
| Ökaryot **ve** ikinci plazmit | p15A değil ColE1 | Konaktaki pRARE zaten p15A taşır; üçüncüsü için pSC101 önerilir |
| Büyük ölçek | Ampisilin yerine kanamisin | Beta-laktamaz besiyerindeki ampisilini parçalar |

Motor kural tabanlıdır: sizin proteininizin katlanma davranışını ya da
gerçek toksisitesini bilemez, bilinen tuzaklara düşmeyen bir başlangıç
noktası verir. Son söz kullanıcıdadır.

Doğrulama: 1296 soru kombinasyonunun tamamı üretilip denetimden
geçirilir ve hiçbiri `hata` seviyesinde bulgu vermez. Bu test
`core/recommend.py` değiştirildiğinde yeniden çalıştırılmalıdır —
ORI-04 çakışması (Rosetta + p15A) tam olarak böyle yakalandı.

### Tasarımcı
Katalogdan parça seçerek vektör kurar, dairesel ve doğrusal haritasını
çizer, tasarımı denetler, GenBank/FASTA olarak dışa aktarır.

Dört hazır şablon: pET tipi ifade, pUC tipi klonlama, arabinozla
indüklenen (P-BAD), sabit ifadeli sentetik devre.

### Dizi analizi
Yapıştırılan ya da yüklenen diziyi inceler: katalog parçalarına karşı
anotasyon, restriksiyon haritası (tek/çift/çok kesen ve kesmeyen
enzimler), sindirim ve fragman simülasyonu, ORF taraması, çeviri,
nadir kodon analizi, Shine-Dalgarno boşluk kontrolü, primer tasarımı.

### Enzimler
REBASE'den gelen **1088 restriksiyon enzimi** (623'ü ticari olarak
satılıyor). Veri Biopython ile birlikte dağıtılır; kaynak
Roberts ve ark., rebase.neb.com.

- **Kütüphane** — ada ya da tanıma dizisine göre arama; tanıma dizisi
  uzunluğu, uç tipi (5'/3' çıkıntı, kör), Tip IIS ve ticari erişilebilirlik
  filtreleri. Her enzim için optimum sıcaklık, beklenen kesim sıklığı,
  izoşizomerler ve metilasyon duyarlılığı.
- **Nereden keser** — tanıma dizisi çift zincir hâlinde, bazlar dikeyde
  hizalı çizilir; `▼` üst zincirdeki, `▲` alt zincirdeki kesimi gösterir.
  İki ok aynı sütundaysa kör uç, kaymışsa yapışkan uç oluşur. Altında
  kesimden sonra oluşan uçlar ayrıca gösterilir.
- **Dizide kesim** — Tasarımcıdaki yapının ya da yapıştırılan bir dizinin
  her kesim noktası, çevresindeki 12 bazla birlikte çift zincir olarak
  çizilir; çıkıntı dizisi gerçek DNA'dan okunur.
- **Uyumlu uçlar** — aynı çıkıntıyı bırakan enzimler. Yönlü klonlamada
  bunlardan kaçınmak, uyumlu uç ararken bulmak gerekir.

Tip IIS enzimleri (BsaI, BsmBI, SapI…) tanıma dizisinin **dışından**
keser, bu yüzden diyagramda çıkıntı bölgesi `N` görünür: çıkıntı dizisi
tanıma bölgesinde değil çevredeki DNA'dadır. Golden Gate montajının
çalışma mantığı budur ve gerçek bir dizide kesim gösterilirken o bölge
DNA'dan okunur.

Doğrulama: kesim geometrisi, gerçek dizilerden çıkarılan çıkıntının
REBASE'in bildirdiği çıkıntıyla birebir eşleşmesiyle sınanır
(5' çıkıntı, 3' çıkıntı ve kör uç enzimleri dahil).

### Referans vektörler
Gerçek plazmitlerin **tam dizileri** üzerinde çalışır — 41 vektör, beş
organizma grubunda. Hedef organizmaya göre süzülür.

| Grup | Vektörler |
| --- | --- |
| E. coli klonlama / ifade | pUC19, pUC18, pBluescript II KS(+), pBR322, pACYC184, pACYC177, pTrc99A, pBAD24, pBAD18, pBAD33, pSC101 |
| BioBrick iskeletleri | pSB1C3, pSB1A3, pSB1K3, pSB4C5 |
| SEVA standart koleksiyonu | pSEVA111/211/311/411/511/611 ve 221/231/241/251/321/331/341/351 |
| Geniş konak aralığı | pBBR1MCS, pBHR1, RSF1010 |
| Genom düzenleme | pKD46, pKD4, pKD3, pE194 |
| Gram-pozitif ve diğer | pUB110, pC194, pNZ8048, pWV01, pCG1 |

Beş sekme:

- **Harita** — katalog parçalarının bu dizide bulunduğu yerler, dairesel
  ve doğrusal harita, seçtiğiniz tek kesen enzim bölgeleriyle.
- **Öznitelikler** — katalog eşleşmeleri (benzerlik yüzdesiyle) ve
  kaydın kendi GenBank öznitelik tablosu, ayrı ayrı.
- **Enzim haritası** — 623 ticari enzimin kaç kez kestiği; tek kesen bir
  enzim seçtiğinizde kesim çift zincirli diyagramla gösterilir.
- **Klonlama planı** — enzim(ler) + insert verince: vektörü kaç yerden
  kestiği, insertin içinde kesip kesmediği, uçların yönlü klonlamaya
  uygun olup olmadığı, iskelet/atılan parça/son yapı boyutları ve
  jeldeki fragman listesi.
- **Parça çıkar** — vektörün bir bölgesini (etiketli bir bölgeyi ya da
  elle verdiğiniz koordinatı) katalog parçası olarak kaydeder; kayıt
  aksesyon ve koordinatı taşır.

### Kodon optimizasyonu ve mRNA yapısı
Dizi analizi sayfasına eklenen iki sekme, ifadenin gerçekten tıkandığı
iki noktaya bakar.

**Kodon optimizasyonu.** CAI ölçer, nadir kodonları ve kümelendikleri
yerleri gösterir, istenirse düzeltir. İki kip var: "sadece nadir
kodonlar" yalnızca sorunlu kodonlara dokunur ve dizinin geri kalanını
korur; "tam optimizasyon" her kodonu kısıtların izin verdiği en iyi
seçenekle değiştirir. Seçtiğiniz enzimlerin tanıma dizileri insertin
içinden kaldırılır - bu, klonlamada kullanacağınız enzimin inserti
kesmesini önler.

**mRNA yapısı.** Shine-Dalgarno dizisini ve başlama kodonunu bulur,
çevrelerindeki saç tokalarını tarar. RBS'i ya da başlama kodonunu
kapatan bir yapı varsa uyarır: promotör güçlü, kodonlar doğru olduğu
halde protein alınamamasının en sinsi sebebi budur ve haritada
görünmez.

### Plazmit denetle
42 denetim kuralı şimdiye kadar yalnızca *parçalardan kurduğunuz*
yapılara uygulanabiliyordu. Bu sayfa elinizdeki bir plazmiti - hocanın
verdiği, bir kayıttan indirdiğiniz ya da daha önce kurduğunuz - aynı
kurallardan geçirir. GenBank, FASTA ya da düz dizi kabul eder.

İki kaynak birleştirilir: katalog taraması (parçalarımız dizide birebir
ya da yakın eşleşmeyle bulunuyor mu) ve dosyanın kendi öznitelik
tablosu. İkincisi olmazsa denetim yanıltır — katalogda karşılığı olmayan
bir orijin taşıyan plazmit için "replikasyon orijini yok" hatası verilir
ki bu yanlıştır.

İki ayar bilinçli:

- **Aynı bölgeye düşen eşleşmeler daraltılır.** Katalogda `bla` birden
  çok kayıt olarak var ve üçü de aynı yere oturuyor; daraltılmazsa
  "aynı plazmitte birden fazla marker" diye yanlış hata çıkıyordu.
- **"Eleman yok" bulguları hataya değil uyarıya çevrilir.** Dışarıdan
  yüklenmiş bir plazmitte eleman *vardır da tanınmamış* olabilir.
  Parçalardan kurduğunuz yapıda bu bilgi kesindir, yüklenende değil.

Tanınmayan bölgeler ayrıca listelenir — en büyüğünü tek tıkla Dizi
analizine gönderip açık okuma çerçevesi arayabilirsiniz.

### Laboratuvar
Tasarımı tezgâha taşıyan sayfa. Dört sekme:

- **Protokol** — tasarımın kendisinden üretilen, sayılarla dolu bir plan:
  hangi antibiyotik kaç µg/mL, ligasyona kaç ng insert, hangi
  indükleyici hangi derişimde. Ligasyon kütlesi molar orandan hesaplanır
  (`ng insert = ng vektör × insert bp / vektör bp × oran`) — eşit kütle
  koymak ligasyonun en sık sessiz hatasıdır. Markdown olarak indirilir.
- **Jel simülasyonu** — sindirimin agaroz jelde nasıl görüneceği.
  Fragman listesi iki şeyi gizler: yakın boyutlar tek bant olur ve kısa
  fragmanlar az boya bağladığı için sönük kalır. İkisi de burada
  çizilir, agaroz yüzdesi fragmanlara göre önerilir.
- **Protein** — insertin açık okuma çerçevesi çevrilir ve kütle, pI,
  sönüm katsayısı, GRAVY, kararsızlık indeksi hesaplanır. Her sayının
  pratik karşılığı yazılır: SDS-PAGE'de bandı nerede arayacağınız,
  A280'den derişimi nasıl çevireceğiniz, inclusion body riski, tampon
  pH'ının pI'ya göre nerede durması gerektiği.
- **Montaj yöntemleri** — Golden Gate, Gibson ve bölge-yönelimli
  mutagenez. Her yöntemin kendi tuzağı denetlenir: Golden Gate'te iki
  birleşme noktası aynı (ya da birbirinin tersi) çıkıntı alırsa parçalar
  yanlış sırada bağlanır ve ürün jelde doğru boyutta göründüğü için hata
  dizilemeye kadar fark edilmez; Gibson'da örtüşmenin Tm'i ve GC'si
  tutmazsa birleşme olmaz. Sipariş edilecek primerler kuyruklarıyla
  birlikte üretilir ve CSV olarak indirilir.
- **Dizileme doğrulaması** — iki bölüm. *Dizileme planı* yapıyı baştan
  sona okutmak için kaç primer gerektiğini ve nereye konacaklarını
  hesaplar; Sanger okumasının ~900 nt olduğunu ve primerin hemen
  ardındaki ~50 bazın güvenilmez olduğunu hesaba katar, kapanmayan
  bölgeleri açıkça gösterir. *Okumayı karşılaştır* dizileme sonucunu
  tasarımla hizalar (Biopython `PairwiseAligner`), her farkı çıkarır ve
  kodlayan bölgedekilerin amino asit karşılığını hesaplar: sessiz mi,
  değiştiren mi, erken durdurma mı. Ters zincirden okunan diziyi kendi
  bulup çevirir.
- **Sorun giderme** — iki bölüm. *Çift kesim uyumluluğu* iki enzimin
  aynı tüpte kullanılıp kullanılamayacağını denetler: sıcaklıklar,
  metilasyon, uç uyumu. Tampon seçimi bilinçli olarak dışarıda bırakıldı
  - hangi enzimin hangi tamponda yüzde kaç çalıştığı üreticiye özgü
  ölçüm verisidir, uygulama onu tahmin etmek yerine nereye bakılacağını
  söyler. *Belirtiden yola çık* beş yaygın belirti için sebepleri
  sıralar ve tasarımınıza bakarak hangisinin sizde gerçekten var
  olduğunu işaretler: T7 promotörlü yapı DH5-alfa'ya konmuşsa "sizde
  var", BL21(DE3)'e konmuşsa o sebep elenir.
- **Rapor** — harita, parça listesi, denetim sonucu, protein özellikleri
  ve protokol tek bir HTML dosyasında. Harita ve jel dosyanın içine
  gömülür, ayrı görsel gerekmez; tarayıcıda açılıp PDF'e basılabilir.

### Hastalık vektörü
Literatürden bir terapötik canlı bakteri vektörünü alıp revize etmek ve
hastalığa göre işe yarayıp yaramayacağını denetlemek için ayrı bir
sayfa. Katalogda **12 hastalık ve 45 strateji** var; 18 atfın hepsi
Europe PMC'den PMID ve DOI'siyle **doğrulanmıştır**. Liste uzunluk olsun
diye uzun değil: her hastalık farklı bir tasarım desenini gösterir.

| Hastalık | Tasarım deseni | Kaynak |
|---|---|---|
| Fenilketonüri (PKU) | Substratı tüketen enzim (PAL) + taşıyıcı (PheP) | Isabella 2018, *Nat Biotechnol* 36:857-864, PMID 30102294 |
| Hiperammonemi | Azotu biyosenteze çeken yol (argA) | Kurtz 2019, *Sci Transl Med* 11:eaau7975, PMID 30651324 |
| Hiperoksalüri | İki enzimli döngü (oxc + frc) | Lubkowicz 2022, *Mol Syst Biol* 18:e10539, PMID 35253995 |
| Homosistinüri | Substratı tüketen enzim, taşıyıcısız (MGL) | Perreault 2024, *Cell Host Microbe* 32:382-395.e10, PMID 38309259 |
| İnflamatuar bağırsak hastalığı | Salgılanan protein terapötiği (IL-10 **ve** elafin) | Steidler 2000, *Science* 289:1352-1355, PMID 10958782 · Motta 2012, *Sci Transl Med* 4:158ra144, PMID 23115353 |
| Laktoz intoleransı | Hücre içi enzim, substrat alımına bağlı (beta-galaktosidaz) | Li 2012, *Appl Microbiol Biotechnol* 96:1499-1506, PMID 22395907 |
| Tip 2 diyabet | Salgılanan peptit terapötiği (GLP-1) | Duan 2015, *Diabetes* 64:1794-1803, PMID 25626737 |
| *P. aeruginosa* enfeksiyonu | Algılayıcı-öldürücü devre (LasR → pyocin S5) | Saeidi 2011, *Mol Syst Biol* 7:521, PMID 21847113 |
| MRSA kolonizasyonu | Hücre duvarını kesen enzim (lizostafin) | Guan 2022, *PLoS One* 17:e0276795, PMID 36520793 |
| Kolera | Sinyal taklidi — virülansı kapatır (cqsA) | Duan 2010, *PNAS* 107:11260-11264, PMID 20534565 |
| AI-2 quorum söndürme | Sinyali ortamdan çeken enzim (LsrK) | Xavier 2005, *Nature* 437:750-753, PMID 16193054 |
| Antibiyotiğe bağlı disbiyoz | Bağırsakta antibiyotik parçalayan enzim (P1A) | Kokai-Kun 2019, *Lancet Infect Dis* 19:487-496, PMID 30885591 |

Testler atıfları de denetler: PMID'si ya da DOI'si olmayan bir atıf
eklenirse `test_atiflarda_pmid_ve_doi_var` düşer.

Sayfa her hastalık için aynı akışı yürütür: makaledeki vektör ile bizim
revizyonumuzu yan yana koyar, her satırı **Literatür** (makaleden) ya da
**Öneri** (bizim kararımız) diye etiketler, mekanizmayı ve metabolik
yolu çizer, sonra idame/içerme stratejilerini aynı denetimden geçirir.

**İdame merdiveni.** Stratejiler rastgele değil: aynı dört soru her
hastalıkta aynı sırada sorulur, böylece karşılaştırılabilir olurlar.
E. coli konaklı hastalıklarda (A) toksin-antitoksin idamesi, (B) dapA
tamamlama, (C) antibiyotik markerlı prototip, (D) hiç idame yok.
Gram-pozitif konaklı olanlarda (A) gıda sınıfı konak, (B) E. coli'de ön
üretim, (C) ön üretim + ccdAB. Denetim her basamağa farklı şey söyler ve
bunu görmek merdivenin asıl amacıdır.

**Terapötik genler.** On beş gen UniProt + NCBI'dan çekilmiştir ve her
birinin nükleotid dizisi çevrilip UniProt'un bildirdiği proteine
**birebir eşit olduğu doğrulanmıştır** — eşit olmayan kayıt dosyaya
yazılmaz (`tools/fetch_disease_parts.py`):

| Parça | Gen | Organizma | UniProt | Uzunluk |
|---|---|---|---|---|
| CDS_ARGA | argA | *E. coli* | P0A6C5 | 1332 bp |
| CDS_OXC | oxc | *O. formigenes* | P40149 | 1707 bp |
| CDS_FRC | frc | *O. formigenes* | O06644 | 1287 bp |
| CDS_MDEA | mdeA / MGL | *P. putida* | P13254 | 1197 bp |
| CDS_IL10 | IL10 | *H. sapiens* | P22301 | 537 bp |
| MRK_THYA_LL | thyA | *L. lactis* | Q02Y43 | 840 bp |
| CDS_LASR | lasR | *P. aeruginosa* | P25084 | 720 bp |
| CDS_PYOS5 | pyoS5 | *P. aeruginosa* | Q9I4Y4 | 1497 bp |
| CDS_CQSA | cqsA | *V. cholerae* | Q9KM65 | 1170 bp |
| CDS_LACZ_LD | lacZ | *L. delbrueckii* | P0C1Y0 | 3024 bp |
| CDS_BLAC_BL | penP (P1A) | *B. licheniformis* | P00808 | 924 bp |
| CDS_LSRK | lsrK | *E. coli* | P77432 | 1593 bp |
| CDS_LSS | lss (olgun) | *S. simulans* | P10547 | 744 bp |
| CDS_GLP1 | GCG → GLP-1(1-37) | *H. sapiens* | P01275 | 117 bp |
| CDS_ELAFIN | PI3 (olgun) | *H. sapiens* | P19957 | 177 bp |

**Bölgesi kesilmiş genler.** Son üçü proteinin tamamı değildir.
Terapötik yapıların çoğu olgun enzimi ya da tek bir peptidi ifade eder;
ön-pro proteinin bütününü koymak yanlış olur. Bölge sınırları
**UniProt'un kendi öznitelik tablosundan okunur**, ezberden yazılmaz:
lizostafin `Chain 248-493`, GLP-1 `Peptide 92-128`, elafin
`Chain 61-117`. Tam kaydın çevirisi önce doğrulanır, kesme ondan sonra
yapılır ve kesilen parçanın çevirisi de ayrıca sınanır. Üçünün de
başlangıç kodonu (ATG) ve stop kodonu **bizim eklediğimizdir** — bu
parçanın notunda açıkça yazılıdır ve bir test notun orada durduğunu
şart koşar. Bağımsız doğrulama bölge uzunluklarını UniProt
koordinatlarından yeniden hesaplar.

Pyocin S5 kaydı **TrEMBL** durumundadır, yani otomatik üretilmiş ve
elle incelenmemiştir; bu fark gizlenmez, parçanın notuna yazılır.

Denetim kuralları:

| Kural | Yakaladığı hata |
|---|---|
| THR-01 | Phe taşıyıcısı (pheP) var ama onu kullanan enzim (PAL) yok |
| THR-02 | PAL var, taşıyıcı yok (verim düşük kalabilir) |
| ANA-01 | Anaerobik promotör: aerobik şişede "çalışmıyor" sonucuna varmayın |
| SEL-01 | dapA markeri, kendi dapA'sı olan konakta seçilim sağlamaz |
| SEL-02 | DAP bağımlı konak: besiyerine DAP eklenmeli |
| BIO-01 | dapA'yı plazmite koymak biyolojik içermeyi kaldırır |
| BIO-02 | Terapötik bakteride antibiyotik direnç geni |
| BIO-03 | Nissle'ın kriptik plazmitleri, orijin uyumu doğrulanmalı |

**Bilinen açıklar, kapatılmadan yazılır.** 45 stratejinin 13'ü
denetimden **bilerek** hata alır ve ikisi de aracın bozukluğu değil,
tasarımın sınırıdır:

* **CORE-01 (7 strateji)** — Gram-pozitif konaklı tasarımlarda
  replikasyon orijini yok, çünkü katalogda L. lactis ya da B. subtilis'te
  çoğalan bir replikon bulunmuyor. Sebebi somut: bu plazmitler yuvarlanan
  halka ile çoğalır ve çoğalma başlangıcı (dso) GenBank kayıtlarında ayrı
  bir öznitelik olarak işaretlenmez, rep geninin içinde durur. Sınırı
  uydurmadan parça kataloğuna eklenemez. İskelet "Referanstan
  kes-yapıştır" sayfasındaki pNZ8048 ve pUB110 tariflerinden gelir.
* **CORE-02 (6 strateji)** — "idame yok" karşılaştırma basamağında
  marker yoktur, yani transformasyondan sonra doğru koloniyi seçmenin
  yolu bile yok. Basamak, merdivenin diğer üçünün ne işe yaradığını
  göstermek için duruyor.

Test paketi bu açıkları **koruma altına alır**:
`test_orijinsiz_yapi_denetimde_yakalaniyor` orijini olmayan her yapının
hata almasını, `test_acigi_olan_hastalik_bunu_yaziyor` hata veren her
hastalığın açığını yazmasını, `test_beklenen_denetim_hatalari_aciklanmis`
ise hata veren her stratejinin o kural kimliğini **kendi metninde**
söylemesini şart koşar. "Her yapıda orijin var" diye bir sınav
kurulmadı, çünkü böyle bir sınav açığı kapatmak yerine testi
gevşetmeye iter.

### Referanstan kes-yapıştır
Vektör sıfırdan kurulmaz: gerçek bir klonlama vektörü alınır, iki farklı
restriksiyon enzimiyle kesilir, araya terapötik kaset yapıştırılır.
**12 tarif** var, **12 ayrı referans vektör** kullanır ve her biri ayrı
bir şey öğretir:

| Tarif | Enzimler | Sonuç | Öğrettiği |
|---|---|---|---|
| pACYC184 → PKU | BamHI/SalI, sonra NcoI/EcoRI | 4245 → 7714 bp | İki aşamalı kurulum; ikinci aşamanın enzimi *birinci aşamanın ürününde* seçilir |
| pUC19 → hiperammonemi | EcoRI/BamHI | 2686 → 4244 bp | Mavi/beyaz tarama: kaset lacZ-α'yı bozar, doğru klon **beyaz** koloni verir |
| pBR322 → hiperoksalüri | EcoRI/BamHI | 4361 → 7238 bp | İnsersiyonel inaktivasyon: doğru klon Amp<sup>R</sup> / Tet<sup>S</sup> olur |
| pSEVA231 → *P. aeruginosa* | HindIII/SpeI | 3123 → 5681 bp | Modüler iskelet; T0/T1 devreyi iskeletten yalıtır |
| pNZ8048 → IBD | XhoI/BclI | 4332 → 3355 bp | Konak iskeleti belirler; ayrıca ürün orijinalden **küçük** olabilir |
| pBluescript II KS(+) → homosistinüri | EcoRI/XhoI | 2961 → 4372 bp | **Uçların uyumluluğu**, enzimlerin farklı olmasıyla aynı şey değil |
| pTrc99A → laktoz intoleransı | BamHI/SalI | 4176 → 7424 bp | Enzimi vektör değil **insert** seçer |
| pBAD24 → kolera | EcoRI/SalI | 4542 → 6165 bp | Arabinozla kademeli ayarlanabilir ifade |
| pSB1C3 → AI-2 söndürme | EcoRI/SpeI | 2070 → 3787 bp | iGEM BioBrick standardı (RFC10) ve **skar** |
| pBBR1MCS → MRSA | HindIII/SalI | 4707 → 5699 bp | Geniş konak aralıklı replikon tasarımın sınırını çizer |
| pUB110 → C. difficile riski | XbaI/BamHI | 4548 → 5411 bp | Yuvarlanan halka replikasyonu ve bedeli |
| pKD46 → diyabet | SalI/BclI | 6329 → 6669 bp | Plazmiti **kaybettirmek** de bir tasarım aracıdır (30/37 °C) |

Her tarif için sayfa altı şeyi gösterir: referans vektörün haritası ve
neden o vektör seçildiği, enzim tablosu (tanıma dizisi, bıraktığı uç, o
aşamada kaç yerden kestiği, Dam metilasyonuna duyarlı mı), aşama aşama
ne çıktı ne girdi, orijinal ile ürünün haritaları yan yana, beklenen jel
deseni ve sınırlılıklar.

**Tarifler gerçekten kurulabilir mi?** İki sınav var ve ikisi de
otomatik:

1. **Kullanılan enzim, eklenecek kasetin içinde kesmemeli.** İnserti o
   iki enzimle hazırlarsınız; kasetin içinde aynı enzimin yeri varsa
   insert ortadan bölünür ve tarif kağıt üzerinde kalır. Üç tarif bu
   sınavı ilk yazıldığında geçemedi — pACYC184/PKU'da HindIII,
   pSEVA231'de EcoRI, pNZ8048'de hem EcoRI hem HindIII kasetin
   içindeydi — ve enzim çiftleri değiştirildi.
2. **Çıkıntılar farklı olmalı.** Enzimlerin farklı olması yetmez: SalI
   ile XhoI ikisi de 5′ TCGA, XbaI ile SpeI ikisi de 5′ CTAG, BamHI ile
   BclI ikisi de 5′ GATC bırakır. Böyle bir çiftte uçlar birbirine uyar,
   insert ters de girer ve klonların yarısı yanlış olur. `cut_and_paste`
   bunu her aşamada denetler ve uyarır.

Buna ek olarak enzim **sırası** da denetlenir. Bölge sol enzimden sağ
enzime doğru saat yönünde alınır: pKD46 tarifinde `BclI+SalI` yazmak
plazmitin neredeyse tamamını atarken `SalI+BclI` yalnızca 10 bp'yi atar.
Bu hata da ilk yazımda yapıldı ve ürün uzunluğuna bakan kontrol yakaladı
(6329 bp'lik vektörden 372 bp'lik bir "ürün" çıkmıştı).

**Orijinal bölgelere ne oldu?** Bu tablonun sayıları diziden ölçülür,
elle yazılmaz. Her bölge adım adım koordinat eşlemesinden geçirilir ve
ürün üzerinde ne kadarının kaldığına bakılır: pUC19 tarifinde lacZ-α
324 bp'den 256 bp'ye düşer (mavi/beyaz taramanın çalışma sebebi tam
olarak bu), pBR322'de tetA 1191 bp'den 902 bp'ye iner. **Kesilen gen
bozulmuş gendir**: kalan parça çerçeve içinde bile olsa tam protein
vermez.

Bunu doğru göstermek ayrı bir işlev gerektirdi. `map_position` tek bir
konumu taşır; bir *bölgenin* iki ucunu ayrı ayrı taşımak, kesim
bölgesini içine alan geni sapasağlam ve insertin üstüne binmiş halde
gösteriyordu — yani lacZ'yi bozmak üzerine kurulu bir tarifte lacZ
haritada bütün görünüyordu. `map_feature` bölgeyi bütün olarak taşır ve
üç sonuç döndürür: `tam`, `kesildi` (yalnızca kalan en büyük parça) ya
da `None` (tamamen atıldı). Dizinin başını saran genleri de (pACYC184'te
cat, 3805..219) doğru taşır. Testler bunu **dizinin kendisiyle** sınar:
taşınan her bölge, ürün üzerinde gerçekten o diziyi göstermek
zorundadır.

**Gerçek dizi, N bloğundan farklı davranır.** pACYC184 tarifi ilk
kurulduğunda ikinci aşama BclI kullanıyordu ve her şey yolundaydı —
çünkü PKU kasetinin dizisi henüz çekilmemişti ve kaset haritada N bloğu
olarak duruyordu. Dizi gerçekleşince BclI'nın birinci aşamanın ürününde
**dört yerden** kestiği görüldü ve yerine NcoI geldi. Ders şu: bir enzim
vektörde tek kesiyor diye *ürününde* de tek keseceği varsayılamaz.
Uygulama bunu her aşamada yeniden hesaplar.

Ürün FASTA olarak indirilebilir ya da tek tuşla **Plazmit denetle**
sayfasına aktarılıp bütün kurallardan geçirilebilir.

### Parça kataloğu
Parçalara göz atar, dizisi eksik olanları doldurur, kendi parçanızı
ekler. Her parçaya bağlanan literatür kaynakları burada görünür.
Kaynağa göre süzülebilir (çekirdek katalog / iGEM / sizin
eklediğiniz).

**iGEM serileri** sekmesi iki güç serisini karşılaştırır: Anderson
sabit promotör ailesi (20 üye) ve Weiss ribozom bağlanma bölgesi
serisi.

### Literatür
Europe PMC üzerinden arama (API anahtarı gerekmez), tasarımdaki her
parça için otomatik kaynak taraması, kaynakça üretimi (numaralı liste /
BibTeX / Markdown), metin ve isimlendirme düzeltmesi, makale metninden
parça adayı çıkarma.

## Doğrulama

Bir tasarım aracının "çalışıyor" olması yetmez, "doğru" olması gerekir.
İkisi için iki ayrı mekanizma var.

### Test paketi

```bash
python -m pytest
```

223 test, yedi dosyada: dizi motoru, tasarım modülleri, veri
bütünlüğü, PKU tasarımı ve simülasyonu, hastalık kayıt defteri ve
referanstan kes-yapıştır. Çoğu bir *değişmezi* sınar — optimizasyon
proteini değiştirmemeli, hizalama kendi kendine %100 benzerlik
vermeli, kesim ürünlerinin toplamı kalıbın uzunluğuna eşit olmalı.
Değişmezler beklenen değer yazmaktan güçlüdür, çünkü girdi değişse de
geçerli kalırlar.

Paket yazıldığından beri **yedi gerçek hata** buldu:

1. Ürettiğimiz GenBank'ın LOCUS satırı sütun kaymasıyla bozuktu
   (Biopython "malformed locus line" uyarısı veriyordu).
2. Kromatogram kalite kırpmasında işaret tersti — algoritma yüksek
   kaliteyi değil düşüğü arıyordu.
3. Kesim bölgesini içine alan gen, haritada bütün ve insertin üstüne
   binmiş görünüyordu; `map_feature` buna karşılık yazıldı. Bu, lacZ'yi
   bozmak üzerine kurulu bir tarifte lacZ'yi sağlam göstermek demekti.
4. pACYC184 tarifinde BclI, gerçek PKU dizisi gelince ürünü dört
   yerden kesiyordu; N bloğuyla sınandığı sürece görünmüyordu.
5. Üç tarifin enzimi, eklenecek kasetin **içinde** kesiyordu: insert
   ortadan bölünürdü, yani tarif kağıt üzerinde kalmıştı.
6. pKD46 tarifinde enzim sırası ters yazılmıştı ve 6329 bp'lik
   vektörden 372 bp'lik bir "ürün" çıkıyordu.
7. İki ön üretim stratejisi T7 promotörünü T7 polimerazı olmayan bir
   konakla eşleştiriyordu — uygulamanın kendi PRO-03 kuralı yakaladı.
   Aracın kendi kuralının, aracı yazanın hatasını bulması iyi bir
   işarettir.

Yedisi de düzeltildi ve yedisinin de testi duruyor.

### Bağımsız doğrulama

```bash
python tools/validate_against_references.py
```

Testler kodun kendi içinde tutarlı olduğunu gösterir; bu betik ise
hesapları **uygulamanın dışındaki** ölçütlerle karşılaştırır. 135
kontrol, sekiz başlıkta:

| Ne | Bağımsız ölçüt |
| --- | --- |
| Restriksiyon haritası | pUC19 MCS'inin yayımlanmış koordinatları ve enzim sırası |
| Erime sıcaklığı | Üç ayrı yöntem (en yakın komşu, GC tabanlı, Wallace) 12 °C içinde uyuşmalı |
| Protein kütlesi | ProtParam yerine Biopython'un ayrı `molecular_weight` işlevi |
| Vektör uzunlukları | NCBI'nın kayıt başlığında bildirdiği değerler |
| Kodon tablosu | Elle yazılmış nadir kodon listesi ve bilinen en iyi kodonlar |
| Parça kataloğu | Parçaların gerçek vektörlerde birebir bulunması |
| Terapötik genler | Biopython'un `translate(cds=True)` işlevi: başlangıç, stop, çerçeve ve erken stop denetimini o yapar |
| Bölge sınırları | UniProt öznitelik tablosundaki koordinatlardan yeniden hesap |
| Tarifler | Ürün uzunluğu atılan/eklenenden yeniden hesaplanır; kaset ürün dizisinde birebir aranır |

Bu betik de bir şey yakaladı — ama uygulamada değil, karşılaştırma
tablosunda. Aşağıya bakın.

### Koordinat sözleşmesi

Uygulama **kesim** konumunu verir; yayımlanmış plazmit haritaları ve
üretici katalogları genellikle **tanıma dizisinin başlangıcını** verir.
İkisi aynı değildir ve fark enzime göre değişir: EcoRI (G^AATTC) için
1, SacI (GAGCT^C) için 5, SmaI (CCC^GGG) için 3.

Bunu fark etmemizin yolu bağımsız doğrulama oldu: on enzimin onunda da
sayılar tutmadı. Farklar incelenince rastgele olmadıkları, her birinin
o enzimin kesim ofsetine eşit olduğu görüldü. Uygulama doğruydu;
karşılaştırma tablosundaki iki sayı (PstI ve SphI) ise yanlıştı.

Sonuç olarak arayüz artık **iki konumu birden** gösteriyor ve aradaki
farkı açıklıyor. pUC19 MCS'i şöyle doğrulanıyor:

| Enzim | Tanıma | Kesim | Enzim | Tanıma | Kesim |
| --- | --- | --- | --- | --- | --- |
| EcoRI | 396 | 397 | XbaI | 423 | 424 |
| SacI | 402 | 407 | SalI | 429 | 430 |
| KpnI | 408 | 413 | PstI | 435 | 440 |
| SmaI | 412 | 415 | SphI | 441 | 446 |
| BamHI | 417 | 418 | HindIII | 447 | 448 |

### Sürüm ve veri kimliği

NCBI kayıtları revize edilir, iGEM Registry güncellenir. Altı ay sonra
aynı tasarım farklı bir sonuç verirse, sebebin kodda mı veride mi
olduğunu bilmek gerekir. Her veri dosyasının SHA-256 özeti
`data/CHECKSUMS.txt` içinde durur, testler tutarlılığı sınar ve
kaydedilen her yapı hangi sürümle, hangi veri kimliğiyle üretildiğini
taşır. Uygulamanın altında ve raporun başlığında bu kimlik yazılıdır.

Veri çekme betiklerinden birini çalıştırdıktan sonra:

```bash
python tools/update_checksums.py
```

## Veri kaynakları ve lisanslar

Uygulamadaki dizilerin hiçbiri ezberden yazılmadı; hepsi açık bir
kaynaktan çekildi ve kaydı depoda duruyor. Addgene ve SnapGene içeriği
kullanılmadı — o veritabanlarının kayıtları kendi kullanım koşullarına
tabidir ve kopyalanamaz.

| Kaynak | Ne geliyor | Lisans / durum | Çeken betik |
| --- | --- | --- | --- |
| NCBI Nucleotide | 11 referans vektörün tam dizisi ve GenBank öznitelikleri | NCBI veritabanı içeriğine telif kısıtı koymaz | `tools/fetch_reference_vectors.py` |
| iGEM Registry | 62 parça + 4 BioBrick iskeleti | Creative Commons Attribution 4.0 (CC BY 4.0) | `tools/fetch_igem_parts.py` |
| UniProt | 8 sinyal peptidinin kesim noktası | CC BY 4.0 | `tools/fetch_signal_peptides.py` |
| REBASE (Biopython ile dağıtılır) | 1088 restriksiyon enzimi, 623'ü ticari | Biopython paketiyle birlikte gelir | — |
| Europe PMC | makale üst verisi ve kaynakça | açık API, anahtar gerekmez | çalışma anında, Literatür sayfası |

iGEM kayıtları CC BY 4.0 olduğu için her parçanın `accession` alanına
Registry kimliği ve lisansı yazılır, Parça kataloğu sayfasında da
kaynak satırı gösterilir. `tools/fetch_igem_parts.py` her parçanın
lisansını API'den okur ve `cc-by` ile başlamayan bir kayıt gelirse onu
atlar — lisans varsayılmaz, doğrulanır.

Lisans ayrıntıları: kodun kendisi [MIT](LICENSE), veri dosyalarının
lisansları [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md) içinde
dosya dosya yazılı.

Betikler yalnızca veriyi yenilemek için elle çalıştırılır:

```bash
python tools/fetch_reference_vectors.py
```

```bash
python tools/fetch_igem_parts.py
```

İkisi de yeniden çalıştırılabilir: var olan kayıtları atlar, iGEM
betiği hız sınırına (429) takılırsa bekleyip kaldığı yerden devam eder.
Referans vektör betiği indirdiği her kaydın uzunluğunu beklenen
uzunlukla karşılaştırır; tutmazsa o vektörü yazmaz ve hata verir.

## Tasarım kararları

### Dizi "tam" ve "eksik" ayrımı
Katalogdaki 60 parçanın tamamı artık gerçek dizi taşır (`seq_status`
= `full`). Başlangıçta uzun CDS'ler ve orijinler bilinçli olarak dizisiz
bırakılmıştı: ezberden uzun gen dizisi yazmak sessiz hatalar üretir ve
sessiz hata, bu aracın önlemeye çalıştığı şeyin ta kendisidir.

Çözüm, dizileri uydurmak değil **kaynağından çekmek** oldu. 29 parçanın
dizisi NCBI'dan, her biri için açıkça yazılmış bir erişim numarası +
koordinat çiftinden alındı (aşağıdaki tablo). Hiçbir baz elle yazılmadı;
her parça adlı bir kayıttaki adlı bir aralıktan gelir ve
yeniden üretilebilir.

`stub` mekanizması yerinde duruyor, çünkü kendi parçanızı dizisiz
ekleyebilirsiniz: eksik parçalar montajda doğru uzunlukta `N` bloğu
olarak yer alır, harita ve uzunluk hesapları doğru kalır, arayüz bunu
her yerde işaretler ve `PHY-01` bulgusu çıkar.

### Dizi kaynakları ve koordinatlar

Koordinatlar 1-bazlı ve kapalı aralıktır (GenBank üslubu). Her dizi
çekildikten sonra doğrulandı: uzunluk, 3'ün katı olma, başlama kodonu,
iç stop kodonu yokluğu ve çevirinin referans protein kaydıyla
birebir aynı olması.

| Parça | Ad | Kaynak | Uzunluk (bp) |
|---|---|---|---|
| `PRO_ARABAD` | araBAD promotörü | NC_000913.3: complement(70049..70386) | 338 (eski tahmin 130) |
| `PRO_RHABAD` | rhaBAD promotörü | NC_000913.3: complement(4097449..4097735) | 287 (eski tahmin 160) |
| `PRO_PL` | lambda P-L | NC_001416: complement(35575..35654) | 80 |
| `TAG_GST` | GST etiketi | M14654 (Sj26): 17..670 | 654 (eski tahmin 669) |
| `TAG_MBP` | MBP füzyonu | NC_000913.3 malE: ATG + complement(4245232..4246341) | 1113 |
| `TAG_SUMO` | SUMO füzyonu | NM_001180818 (SMT3): 1..294 | 294 (eski tahmin 300) |
| `TAG_TRXA` | TrxA füzyonu | NC_000913.3 trxA: 3965761..3966087 | 327 |
| `TER_RRNB` | rrnB T1T2 | NC_000913.3: 4171757..4171974 | 218 (eski tahmin 360) |
| `ORI_PUC` | pUC ori | L09137 (pUC19): complement(805..1459) | 655 (eski tahmin 589) |
| `ORI_PBR322` | pBR322 ori | J01749 (pBR322): complement(2472..3126) | 655 (eski tahmin 589) |
| `ORI_COLE1` | ColE1 ori | J01566 (ColE1): 601..1257 | 657 (eski tahmin 600) |
| `ORI_P15A` | p15A ori | X06403 (pACYC184): 581..1493 | 913 (eski tahmin 570) |
| `ORI_PSC101` | pSC101 ori + repA | K00042: 370..1700 | 1331 (eski tahmin 1000) |
| `ORI_R6K` | R6K gamma ori | M65025: 1..442 | 442 (eski tahmin 400) |
| `ORI_F1` | f1 / M13 ori | V00604 (M13): 5501..6005 | 505 (eski tahmin 460) |
| `MRK_BLA` | bla (AmpR) | J01749 bla: complement(3293..4153) | 861 |
| `MRK_KAN` | nptII (KanR) | U00004 (Tn5) nptII: 1551..2345 | 795 (eski tahmin 816) |
| `MRK_CAT` | cat (CmR) | X06403 cat: complement(join(3805..4245,1..219)) | 660 |
| `MRK_TET` | tetA (TetR) | J01749 tetA: 86..1276 | 1191 |
| `MRK_AADA` | aadA (SpecR) | U57972 aadA: 1..792 | 792 |
| `MRK_CCDB` | ccdB | AP001918 (F plazmiti) ccdB/letD: 46560..46865 | 306 |
| `REG_LACI` | lacI | NC_000913.3 lacI: complement(366428..367510) | 1083 |
| `REG_ARAC` | araC | NC_000913.3 araC: 70387..71265 | 879 |
| `REG_TETR` | tetR | J01830 (Tn10) tetR: complement(903..1526) | 624 |
| `REG_ROP` | rop / rom | J01749 rop: 1915..2106 | 192 |
| `RPT_SFGFP` | sfGFP | MW132720: 1..717 | 717 (eski tahmin 714) |
| `RPT_MCHERRY` | mCherry | AY678264 (Shaner 2004): 1..711 | 711 |
| `RPT_LACZ` | lacZ | NC_000913.3 lacZ: complement(363231..366305) | 3075 |
| `RPT_LACZA` | lacZ-alfa | L09137 (pUC19): complement(146..469) | 324 (eski tahmin 450) |

13 parça katalogdaki tahminle birebir tuttu; 16'sında tahmin
düzeltildi. Her parçanın `notes` alanı dizinin nasıl türetildiğini ve
tahminin neden değiştiğini tek paragrafta anlatır.

Birkaçı ayrıca anlatmaya değer, çünkü iki bağımsız kayıt birbirini
doğruladı:

- **pUC / pBR322 / ColE1 orijinleri.** Hiçbir kayıtta "ori" diye
  sınırları çizilmiş bir özellik yok. ColE1 kaydının (J01566) kendi
  anotasyonlarından (RNA II başlangıcı 641, RNA I 645..751, replikasyon
  orijini 1197) yola çıkıldı; bu bölge pBR322'ye %93.8 kimlikle
  hizalandığında ColE1'in ori noktası pBR322'nin kendi `rep_origin`
  anotasyonuna (2535) 3 bp içinde oturdu. Aynı bölge pUC19'da
  `complement(805..1459)`. Sonuçta çıkan iki dizi 655 bp boyunca
  **tek bir bazda** ayrılıyor (152. baz, pBR322'de C, pUC19'da T) -
  yani katalog notlarında anlatılan pUC yüksek kopya mutasyonunun ta
  kendisi, bağımsız olarak doğrulanmış halde.
- **lambda P-L.** Kayıtta promotör özelliği yok, sadece oL1/oL2/oL3
  operatörleri var. Bölge pL yönünde okunduğunda -35 = `TTGACA` ve
  tam 17 bp aralıkla -10 = `GATACT` çıkıyor, +1 ise lambda 35582'ye
  denk geliyor; yani promotör konumu operatör anotasyonlarından
  bağımsız olarak türetildi.
- **cat.** pACYC184'te *cat* CDS'i kaydın numaralandırma başlangıcını
  aştığı için kayıtta CDS olarak işaretli değil. İki zincirde 660
  nt'lik tek açık okuma çerçevesi arandı; çeviri UniProt P62577 (Tn9
  CAT, 219 aa) ile birebir aynı çıktı.
- **araBAD promotörü.** 3' ucu araB'nin ATG'sinden hemen önce bitiyor
  ve pBAD vektörlerinin bilinen `ACTGTTTCTCCATACCCGTTTTTTTGG` ucuyla
  birebir uyuşuyor.
- **sfGFP.** NCBI'da "superfolder GFP" adıyla duran kayıtlar birbirini
  tutmuyor. avGFP (P42212) ile karşılaştırma yapıldı: `MW132720` tam
  olarak 11 kanonik sfGFP mutasyonunu taşıyor (S30R, Y39N, F64L, S65T,
  F99S, N105T, Y145F, M153T, V163A, I171V, A206V), `JQ341914` ise bu
  imzayı taşımıyor (Y39I, N105K ve 10 ek fark). İlki alındı. Bir
  kaydın adının doğru olması dizisinin doğru olduğunu göstermiyor.

### Yaklaşık anotasyon (`annotate_similar`)
Gerçek plazmitlerde tam eşleşme yetmiyor. Aynı gen farklı vektör
soylarında birkaç nükleotid farkla taşınır: `bla` pBR322 ve pACYC177'de
katalogdaki diziyle birebir eşleşirken pUC19, pBAD24, pTrc99A ve
pBluescript'te %99,8 benzer çıkıyor. `pUC ori` ile `pBR322 ori`
arasındaki fark ise tek nükleotid — pUC'yi yüksek kopyalı yapan
mutasyon tam olarak bu.

Sadece tam eşleşme arasaydık pUC19 haritasında `bla` hiç görünmezdi. Bu
yüzden ikinci bir katman var: parçanın başından, ortasından ve sonundan
alınan 20 bp'lik tohumlarla aday konumlar bulunur, sonra parçanın
tamamı boyunca uyuşmayan baz sayılır. Eşik arayüzden ayarlanabilir
(varsayılan %90).

Sınırı: **hizalama yapılmaz.** Yalnızca aynı uzunluktaki nokta farkları
yakalanır; insersiyon ya da delesyon taşıyan varyantlar bulunamaz.
Tablo her satırın benzerlik yüzdesini gösterir, %100 olmayanlar
varyanttır.

Aynı bölgeye iki parça birden düşerse (ori örneğinde olduğu gibi) aynı
kategoriden çakışanlardan yalnızca en iyisi tutulur — yoksa harita aynı
yeri iki kez etiketler.

### Bileşik parçalar (`device` kategorisi)
iGEM'in bazı kayıtları tek parça değil, kendi promotörünü, ribozom
bağlanma bölgesini, kodlayan dizisini ve terminatörünü birlikte taşıyan
bir kasettir (BBa_J04450, BBa_I13504, BBa_E0240, BBa_F2620). Bunlara
`device` kategorisi verildi ve montaj sırasında kasetin başına
konuyorlar.

Denetimde iki kural buna bağlı: bileşik parça varken "promotör yok",
"RBS yok", "terminatör yok" kuralları uygulanmaz (DEV-01) — eleman
eksik değil, parçanın içinde. Buna karşılık ayrıca promotör ya da RBS
seçilmişse uyarı çıkar (DEV-02), çünkü montajda iki kaset üst üste
biner.

### Sinyal peptitleri
Protein modülü "inclusion body riski yüksek" ya da "disülfit için
periplazmik hedefleme gerekir" dediğinde, katalogda o öneriyi
uygulayacak parça yoktu. Uygulama kendi açtığı döngüyü kapatamıyordu.

Sekiz sinyal peptidi eklendi; her biri kendi kategorisinde
(`signal`), çünkü afinite etiketinden işlevsel olarak farklıdır:
N ucunda OLMAK ZORUNDADIR, taşıma sırasında kesilip atılır ve hangi
yoldan taşınacağını belirler. Montaj sırasında etiketten de önce
yerleşir — etiket araya girerse sinyal tanınmaz.

| Yol | Sinyal | Ne zaman |
| --- | --- | --- |
| Sec (post-translasyonel) | pelB, OmpA, PhoA, MalE | Varsayılan periplazmik hedefleme |
| SRP (ko-translasyonel) | DsbA | Sitoplazmada çökelen hedefler |
| Tat (twin-arginine) | TorA | **Katlanmış** protein taşır; kofaktörlü hedefler |
| Sec (Gram-pozitif) | SacB, AmyE | *B. subtilis*; doğrudan besiyerine salgılama |

**Kesim noktası ezberden yazılmadı.** Önemli olan nerede kesildiğidir:
bir baz fazla ya da eksik alırsanız okuma çerçevesi kayar.
`tools/fetch_signal_peptides.py` UniProt'un `Signal` öznitelik
aralığını okur, karşılık gelen nükleotid dizisini NCBI'dan alır ve
**çevirisinin UniProt'un bildirdiği peptide birebir eşit olduğunu
sınar**. Eşit değilse kayıt yazılmaz — nitekim `phoA` ilk denemede
elendi, çünkü GTG ile başlıyor ve çevirici onu valin okuyordu;
bakteriyel başlangıç kodonları düzeltilince geçti.

Beş denetim kuralı eklendi: birden fazla sinyal (hata), Gram
uyuşmazlığı, Tat kapasitesi, floresan muhabirlerin periplazmada
katlanmaması ve klonlama amacında sinyalin anlamsız olması.

### Kombinatoryal kütüphane
Tasarımcıdaki **Kütüphane kur** sekmesi, tek tasarım yerine bir seri
üretir: sabit tuttuğunuz parçaların üstüne, değişken bıraktığınız
kategorilerde bütün birleşimleri kurar, her birini denetimden geçirir ve
beklenen ifade düzeyine göre sıralar.

Sebebi şu: ifade düzeyini önceden hesaplayamazsınız, ama bir aralığı
tarayabilirsiniz. Gerçek laboratuvarda da böyle yapılır.

**Beklenen ifade bir ölçüm değildir** ve arayüzde de böyle yazar. Üç
kaynaktan türetilmiş sıralamadır, her satırda hangisinden geldiği
belirtilir:

1. Anderson promotörlerinde sigma70 konsensüsüne uzaklık — diziden
   hesaplanır, en nesnel olanı budur.
2. Katalogdaki güç etiketi ("çok güçlü", "orta", "zayıf").
3. iGEM Registry başlığındaki strong / medium / weak nitelemesi.

Promotör ve RBS skorları çarpılır, çünkü transkripsiyon ve çeviri
ardışık adımlardır; biri sıfıra yaklaşırsa diğeri telafi edemez.
Mekanizmaya uygun ama yine de kaba bir yaklaşım.

Beş promotör × dört RBS'lik tipik bir seri yaklaşık 11 kat aralık
tarıyor. Seçtiğiniz varyantı tek tıkla Tasarımcıya yükleyebilir, seriyi
CSV olarak indirebilirsiniz.

### Jel simülasyonunun modeli
Agarozda göç mesafesi kabaca `log(boyut)` ile doğrusaldır; jel yüzdesi
hangi aralığın ayrıştığını belirler. İki ayrıntı bilinçli olarak
modellendi:

- **Bant birleşmesi.** Boyutları birbirine %6'dan yakın fragmanlar tek
  bant olarak çizilir, çünkü jelde de öyle görünürler. Tabloda iki
  satır, jelde bir bant.
- **Parlaklık kütleyle orantılı.** İnterkalan boya uzunlukla orantılı
  bağlandığı için, eşit molde olan 100 bp'lik fragman 3000 bp'linin
  otuzda biri kadar sinyal verir. Sönük kalacak bantlar ayrıca uyarı
  olarak yazılır - "bu bandı doğrulama ölçütü saymayın".

Agaroz yüzdesi, en çok fragmanı ayırma aralığına alan değer seçilerek
önerilir; 50 bp'nin altındaki fragmanlar bu seçime katılmaz, çünkü
hiçbir jelde güvenilir biçimde görünmezler.

### Protein hesapları nereden geliyor
Kütle, pI, sönüm katsayısı, GRAVY, kararsızlık indeksi ve ikincil yapı
oranları Biopython'un `Bio.SeqUtils.ProtParam` modülünden gelir (ExPASy
ProtParam ile aynı yöntemler). Bu değerler **tahmindir, ölçüm
değildir**: kararsızlık indeksi diziden hesaplanan bir göstergedir,
gerçek yarı ömür değil; SDS-PAGE'de göç kütleyle birebir orantılı
değildir.

Uygulamanın kattığı şey yorum: hangi sayı hangi deney kararını
değiştiriyor. Örneğin triptofan ve tirozin yoksa sönüm katsayısı sıfır
çıkar ve Nanodrop sıfıra yakın okur - bu protein olmadığı anlamına
gelmez, A280'in o protein için kullanılamayacağı anlamına gelir.

### Protokoldeki sayılar
Antibiyotik ve indükleyici derişimleri standart E. coli çalışma
değerleridir ve tasarımın parçalarından seçilir: `MRK_KAN` varsa
kanamisin 50 µg/mL, `PRO_ARABAD` varsa arabinoz. Laboratuvarınızın
kendi protokolü farklıysa onunki önceliklidir - bu sayılar başlangıç
noktası, reçete değil.

### SEVA adlandırması tahmin edilmedi
pSEVA[marker][orijin][kargo] şeması, kayıtların kendi GenBank öznitelik
tablolarından okundu: pSEVA111'de `bla` ve `R6K`, pSEVA211'de `neo` ve
`R6K`, pSEVA221'de `neo` ve `trfA` görünüyor. Haneler böyle eşleşti,
hiçbir hane ezberden yorumlanmadı.

| Hane | 1 | 2 | 3 | 4 | 5 | 6 |
| --- | --- | --- | --- | --- | --- | --- |
| Marker | bla | neo | cat | aadA | tetA | acc3 |
| Orijin | R6K | RK2 | pBBR1 | pRO1600/ColE1 | RSF1010 | — |

Katalog anotasyonu bunu bağımsız doğruladı: pSEVA341'in dizisinde `cat`
ve pUC orijini, pSEVA111'inkinde R6K gama orijini bulunuyor.

### Konak organizması denetimi
Konak listesinde artık 22 suş var ve beşi E. coli dışı (*B. subtilis*,
*P. putida*, *L. lactis*, *C. glutamicum*). Bu, denetime dört yeni kural
getirdi, çünkü bir E. coli plazmiti Bacillus'ta çoğalmaz ama haritası
gayet doğru görünür:

- **ORI-05** — R6K orijini yalnızca *pir*⁺ konakta çoğalır. Katalogdaki
  pSEVA x1x iskeletleri ile pKD3/pKD4 bu orijini taşıdığı için kural
  gerçekten gerekli; PIR1 suşu bu yüzden eklendi.
- **ORG-01** — orijin konak organizmasında çalışmıyorsa hata.
- **ORG-02** — promotör E. coli'ye özgüyse, başka organizmada sigma
  faktörü ve ribozom bağlanma bölgesi farklı olduğu için uyarı.
- **ORG-04** — E. coli konakta Gram-pozitif replikon varsa mekik vektör
  önerisi.

### Kodon tablosu nereden geliyor
`data/ecoli_codon_usage.json` ezberden yazılmadı:
`tools/fetch_codon_usage.py` E. coli K-12 MG1655 genomunu (NC_000913.3)
indirir, 4226 geçerli kodlayan diziyi okur ve kodonları sayar. CAI'nin
ölçütü olan göreli uyarlanabilirlik (w) değerleri ise çok ifade edilen
78 genden hesaplanır - ribozomal proteinler, uzama faktörleri, glikoliz
enzimleri ve dış zar porinleri. Hangi genlerin ölçüt alındığı çıktının
içine de yazılır.

Bu hesap, uygulamada zaten elle yazılmış olan "nadir kodon" listesini
bağımsız olarak doğruladı: o listedeki dokuz kodonun hepsi genomdan
hesaplanan tabloda w < 0,05 çıkıyor. En iyi kodonlar da bilinen E. coli
tercihleriyle örtüşüyor (Leu için CTG, Arg için CGT, Ala için GCT).

### Saç tokası hesabının sınırı
Sap bölgesinin serbest enerjisi Biopython'un RNA en yakın komşu
tablosundan (`RNA_NN3`, Xia ve ark. 1998) gelir; ilmek için standart
biçimde yaklaşık bir ceza eklenir.

Bu **tam bir katlanma hesabı değildir.** ViennaRNA ya da mfold gibi
programlar bütün olası yapıları birlikte değerlendirip en düşük
enerjili katlanmayı bulur; burada yalnızca tek tek saç tokaları
taranır. Çıkan sayı mutlak bir ΔG olarak değil, "bu bölgede güçlü bir
yapı var mı" sorusunun göreli yanıtı olarak okunmalıdır. Arayüzde de
böyle yazıyor.

Bir ayrıntı: en yakın komşu tabloları her istifi tek bir kanonik
biçimde tutar (`TC/AG` yok ama aynı fiziksel istif olan `GA/CT` var).
Simetrik eşdeğer denenmezse istiflerin yarısı sessizce kaybolur ve
yapılar olduğundan kararsız görünür - kod bunu açıkça ele alıyor.

### Denetim kuralları
Her kural gerçek bir laboratuvar hatasına karşılık gelir ve üç seviyeye
ayrılır: `hata` (yapı çalışmaz), `uyari` (çalışır ama kötü sonuç verir),
`bilgi` (dikkat edilmeli).

Kapsanan başlıca hatalar:

| Kural | Yakaladığı hata |
|---|---|
| PRO-03 | T7 promotörünün T7 polimeraz taşımayan konakta kullanılması |
| REG-01 | lacI / araC / tetR eksikliği, kaçak ifade |
| ORI-04 | Konaktaki pLysS / pRARE ile uyumsuzluk grubu çakışması |
| MRK-04 | İfade plazmiti ile yardımcı plazmitin aynı direnci kullanması |
| MRK-03 | tetA markerı ile P-tet promotörünün aynı plazmitte olması |
| MRK-05 | ccdB taşıyan vektörün normal konakta çoğaltılamaması |
| RPT-01 | lacZ-dM15 taşımayan konakta mavi/beyaz tarama |
| FRM-02 | Insertin stop kodonunun C-terminal etiketi kapatması |
| FRM-05 | pET-28a'da NcoI/NdeI ve stop kodonu kaynaklı etiket kaybı |
| INS-05 | Klonlama enziminin insertin içinde kesmesi |
| ENZ-02 | Dam/Dcm metilasyonunun enzimi bloke etmesi |
| INS-01 | Nadir kodon yükü, Rosetta gereksinimi |
| PRO-04 | Yüksek kopya orijin + çok güçlü promotör, metabolik yük |

### Türkçe metin düzeltme
Gen ve tür adlarının italiği, enzim ve suş yazımı, birim biçimleri
otomatik düzeltilir. **Türkçe ek alan terimler otomatik düzeltilmez**:
ünlü uyumu ve kaynaştırma harfi yüzünden kökü değiştirip eki olduğu gibi
bırakmak bozuk sonuç verir ("ekspresyonu" → "ifadeu"). Bunlar öneri
olarak işaretlenir, çekimi yazar yapar.

## Dosya düzeni

```
streamlit_app.py       gezinme ve paylaşılan durum
app_pages/
  wizard.py            gerekçeli vektör önerisi
  designer.py          vektör kurma, harita, denetim
  vectors.py           referans vektörler (harita, enzim, klonlama)
  analyze.py           dizi analizi
  enzyme_lab.py        enzim kütüphanesi, kesim diyagramları
  inspect.py           var olan plazmiti denetle
  disease.py           hastalık vektörü: 12 hastalık, revizyon, denetim
  refbuild.py          referanstan kes-yapıştır: 12 tarif
  labwork.py           protokol, jel, protein, rapor
  parts.py             katalog yönetimi, iGEM serileri
  literature.py        literatür modülü
core/
  parts_seed.py        60 bakteriyel parça + 8 konak suşu
  pku_seed.py          PKU vektörü parçaları ve Nissle konakları
  pku_sim.py           PKU tasarımının davranış simülasyonu
  disease_seed.py      hastalık kayıt defteri (12 hastalık, 45 strateji)
  refbuild.py          kes-yapıştır hesabı ve 12 tarif
  igem_seed.py         iGEM kayıtlarını katalog parçasına çevirir
  signal_seed.py       sinyal peptitleri
  ref_vectors.py       referans vektörler + editöryel katman
  cloning.py           restriksiyonla klonlama planı ve denetimi
  inspect_plasmid.py   yüklenen plazmiti denetime hazırla
  library.py           kombinatoryal kütüphane
  assembly_methods.py  Golden Gate, Gibson, mutagenez
  troubleshoot.py      çift kesim denetimi, sorun giderme ağacı
  verification.py      dizileme planı ve okuma karşılaştırması
  codon.py             CAI ve kodon optimizasyonu
  structure.py         mRNA saç tokası taraması
  protocol.py          tasarımdan laboratuvar protokolü
  gel.py               agaroz jel simülasyonu
  protein.py           protein özellikleri ve yorumu
  report.py            tek dosyalık HTML rapor
  enzymes.py           REBASE enzim kütüphanesi
  recommend.py         sihirbazın kural motoru
  db.py                SQLite katmanı
  seq_utils.py         dizi motoru (Biopython üzerine)
  assembly.py          montaj, GenBank/FASTA çıktısı
  validate.py          denetim kuralları
  plasmid_map.py       SVG harita çizici
  literature.py        Europe PMC, kaynakça, metin düzeltme
  env.py               paylaşılan sunucu algılama
  provenance.py        sürüm ve veri kimliği
tools/
  fetch_reference_vectors.py   NCBI'dan vektör dizileri
  fetch_igem_parts.py          iGEM Registry'den parçalar
  fetch_codon_usage.py         E. coli kodon kullanım tablosu
  fetch_pku_parts.py           stlA, PfnrS, ccdAB dizileri (NCBI)
  fetch_signal_peptides.py     9 sinyal peptidi (UniProt)
  fetch_disease_parts.py       15 terapötik gen (UniProt + NCBI)
  extract_parts_from_vectors.py  vektörlerden parça çıkarma
  validate_against_references.py  bağımsız doğrulama
  update_checksums.py          veri özetlerini yenile
tests/
  test_sequences.py    dizi motoru
  test_design.py       tasarım modülleri
  test_data.py         veri bütünlüğü ve denetim kuralları
  test_pku.py          PKU parçaları, konakları ve kuralları
  test_pku_sim.py      PKU simülasyonu
  test_disease.py      hastalık kayıt defteri ve terapötik genler
  test_refbuild.py     kes-yapıştır hesabı, map_feature, tarifler
data/
  reference_vectors/   126 vektörün dizisi ve öznitelikleri (JSON)
  igem_parts.json      iGEM parça kayıtları
  extracted_parts.json vektörlerden çıkarılmış parçalar
  signal_peptides.json 9 sinyal peptidi (UniProt)
  pku_parts.json       PKU genlerinin çekilmiş dizileri
  disease_parts.json   15 terapötik gen (UniProt + NCBI)
  ecoli_codon_usage.json  genomdan hesaplanan kodon tablosu
  CHECKSUMS.txt        izlenen veri dosyalarının SHA-256 özetleri
  vectorlab.db         veritabanı (ilk çalıştırmada kurulur)
```

## Bilinen sınırlar

### Dizi sınırlarının seçimi
Katalogda dizisi eksik parça kalmadı, ama birkaç yerde **sınırı nereye
çekeceğiniz bir karar**, kayıttan okunan bir gerçek değil. Bunlar
bilinçli seçimlerdir ve kullanmadan önce bilmeniz gerekir:

- **Orijinlerin sınırları konvansiyondur.** pUC / pBR322 / ColE1 için
  sınır "RNA II başlangıcından 40 bp önce, replikasyon orijininden 60 bp
  sonra" kuralıyla çizildi. Bu bölge işlevsel orijini (iki RNA'nın
  promotörü + ori) kesinlikle içeriyor, ancak SnapGene'in aynı parçaya
  verdiği 589 bp'lik sınırla birebir aynı değil. Klonlamada gerçek
  sınır sizin seçtiğiniz restriksiyon bölgelerinin nerede olduğuyla
  belirlenir.
- **`ORI_PSC101` orijin + repA'yı birlikte içerir** (1331 bp). pSC101
  orijini RepA proteini olmadan replike olmaz, dolayısıyla tek başına
  "ori" parçası yanıltıcı olurdu. Buna karşılık **`ORI_R6K`,
  pi (pir) proteinini içermez**: gamma orijininin koşulluluğu zaten
  pir'in konaktan gelmesine dayanır.
- **`PRO_ARABAD` ve `PRO_RHABAD` genler arası bölgenin tamamıdır.**
  Yani düzenleyicinin ters yönlü kendi promotörünü (P-C, P-rhaSR) de
  içerirler - pBAD ve pRHA vektörlerinde olduğu gibi. Tek yönlü bir
  promotör parçası bekliyorsanız bu beklentinizi karşılamaz.
- **`TAG_MBP`'nin başındaki ATG kayıttan gelmez.** malE CDS'inden
  sinyal peptidi (kodon 1-26) çıkarıldı ve pMAL vektörlerindeki olgun
  MBP gibi başa ATG eklendi. Eklenen tek şey bu üç bazdır.
- **`TAG_SUMO` Smt3'ün 1-98. kodonlarıdır**, tam CDS değil: Ulp1
  Gly97-Gly98'den sonra kestiği için 99-101 kuyruğu çıkarıldı.
- **`TER_RRNB` 5S geninin bitiminden sonrasıdır** (218 bp). Vektörlerde
  "rrnB T1T2" diye geçen parçalar genelde daha uzundur, çünkü 5S geninin
  içinden başlayan bir restriksiyon fragmanıdır; buradaki parça T1 +
  ara bölge + T2'nin kendisidir.
- Doldurma sırasında ortaya çıkan ve işe yarayan bir ayrıntı:
  araBAD promotörünün içinde **gerçek bir BamHI bölgesi** var
  (`ATTAGCGGATCC...`). Dizi eksikken bu görünmüyordu; artık BamHI ile
  klonlama yapan bir tasarımda restriksiyon haritasında çıkıyor.

### Diğer
- Anotasyon yalnızca dizisi tam olan parçaları bulabilir; bu, kendi
  eklediğiniz `stub` parçalar için geçerli olmaya devam eder.
- Diziler çekildiği andaki NCBI sürümünü yansıtır. Kayıtlar revize
  edilebilir; erişim numarası + koordinat yazılı olduğu için
  karşılaştırma her zaman yapılabilir.
- Golden Gate çıkıntıları varsayılan olarak parçaların birleştiği
  yerdeki doğal diziden alınır (izsiz montaj). Çakışma çıkarsa uygulama
  hesaplanmış bir küme önerir: 256 olası dört-baz çıkıntıdan,
  palindromik olmayan ve birbirinden en az iki bazla ayrılan en büyük
  ayrık küme. Yayımlanmış bir MoClo tablosu kopyalanmadı - küme
  hesaplanıyor ve ölçütleri kodda açıkça duruyor. MoClo'nun düzey
  mantığı (level 0/1/2 kasetleri) henüz yok.
- pET serisi referans vektörler arasında yok: Novagen kayıtları
  GenBank'ta yayımlanmamış (arama boş dönüyor), Addgene'deki kopyalar
  ise kopyalanamaz. Katalogda pET'in parçaları (T7lac promotörü,
  pET-28a MCS'i, T7 terminatörü) tek tek var, tam plazmit dizisi yok.
- Anderson promotör tablosundaki "konsensüse benzerlik" sütunu
  hesaplanmış bir değerdir, ölçülmüş promotör gücü değil. İkisi aynı
  yönde gider ama birebir değildir; ölçüm değerleri Registry'nin kendi
  parça sayfalarındadır.
- Jel simülasyonunda kesilmemiş plazmit tek bant olarak çizilir.
  Gerçekte dairesel plazmit üç bant verir (süperkıvrımlı, açık
  dairesel, doğrusal) ve süperkıvrımlı olan kendi boyutundan daha hızlı
  göçer. Kesilmiş fragmanlar için model doğrudur, kesilmemiş kontrol
  için basitleştirmedir.
- Hiçbir tasarım deneysel olarak sınanmadı. Uygulama "bu tasarımda şu
  hata var" diyebilir, "bu tasarım şu kadar protein verir" diyemez;
  ifade düzeyini, katlanmayı ve gerçek verimi hiçbir tasarım aracı
  hesaplayamaz.
- Editöryel bilgiler (kopya sayıları, suş özellikleri, protokol
  derişimleri) kaynak gösterilmeden yazıldı. Uydurma atıf numarası
  yazmamak bilinçli bir karardı; sonucu, bu değerlerin ders kitabı
  düzeyinde doğru ama referanssız olmasıdır. Her kaydın `lit_query`
  alanı, Literatür sayfasının o konuyu aramasını sağlar.
- Metilasyon duyarlılığı listesi elle küratörlüdür (15 enzim); REBASE
  bu alanı taşımıyor.
- Protein özellikleri dizi üzerinden yapılan tahminlerdir. Katlanma,
  çeviri sonrası değişiklikler ve gerçek çözünürlük bu sayılardan
  okunamaz; GRAVY yüksek diye protein mutlaka çöker, düşük diye mutlaka
  çözünür demek değildir.
- Çift kesim denetimi tampon uyumluluğunu söylemez. Sıcaklık, metilasyon
  ve uç uyumu REBASE verisinden doğrulanabilir; tampon aktivite yüzdeleri
  üreticiye özgüdür ve uydurulmadı.
- Klonlama planı iki tek kesen enzimde büyük parçayı iskelet sayar. MCS
  içinden kesiyorsanız bu doğrudur; uzak iki bölgeden kesiyorsanız
  hangi parçanın iskelet olduğunu haritadan doğrulayın.
- Europe PMC araması internet bağlantısı ister; bağlantı yoksa diğer
  tüm sayfalar çalışmaya devam eder.
- Harita, çok sayıda kısa parça aynı bölgede toplandığında sıkışır;
  o durumda doğrusal harita sekmesi daha okunaklıdır.

## Haritanın `st.markdown` ile basılması

`st.html` bu Streamlit sürümünde `<svg>` içeriğini tamamen temizliyor
(eleman kabı boş geliyor). `st.image` çalışır ama SVG'yi base64 görsele
çevirdiği için açık/koyu tema uyumu (`currentColor`) ve fare üzerine
gelince çıkan `<title>` ipuçları kaybolur. Bu yüzden harita
`st.markdown(..., unsafe_allow_html=True)` ile basılır; haritadaki tüm
metinler `html.escape` ile kaçırıldığından parça adları ham HTML
taşıyamaz.
