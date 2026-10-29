# VectorLab

Bakteriyel vektör tasarım, analiz ve literatür aracı. Addgene + SnapGene
ikilisinin ders düzeyinde yerine geçecek, tek dosyalık veritabanıyla
çalışan yerel bir uygulama.

**Kapsam bilinçli olarak sadece bakteriyel sistemlerdir** (ağırlıklı
*E. coli*). Memeli, bitki, böcek, maya elemanları yoktur.

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

114 test, üç dosyada: dizi motoru, tasarım modülleri ve veri bütünlüğü.
Çoğu bir *değişmezi* sınar — optimizasyon proteini değiştirmemeli,
hizalama kendi kendine %100 benzerlik vermeli, kesim ürünlerinin
toplamı kalıbın uzunluğuna eşit olmalı. Değişmezler beklenen değer
yazmaktan güçlüdür, çünkü girdi değişse de geçerli kalırlar.

Paket ilk çalıştırmada iki gerçek hata buldu: ürettiğimiz GenBank'ın
LOCUS satırı sütun kaymasıyla bozuktu (Biopython "malformed locus line"
uyarısı veriyordu) ve kromatogram kalite kırpmasında işaret tersti —
algoritma yüksek kaliteyi değil düşüğü arıyordu. İkisi de düzeltildi.

### Bağımsız doğrulama

```bash
python tools/validate_against_references.py
```

Testler kodun kendi içinde tutarlı olduğunu gösterir; bu betik ise
hesapları **uygulamanın dışındaki** ölçütlerle karşılaştırır. 79
kontrol, altı başlıkta:

| Ne | Bağımsız ölçüt |
| --- | --- |
| Restriksiyon haritası | pUC19 MCS'inin yayımlanmış koordinatları ve enzim sırası |
| Erime sıcaklığı | Üç ayrı yöntem (en yakın komşu, GC tabanlı, Wallace) 12 °C içinde uyuşmalı |
| Protein kütlesi | ProtParam yerine Biopython'un ayrı `molecular_weight` işlevi |
| Vektör uzunlukları | NCBI'nın kayıt başlığında bildirdiği değerler |
| Kodon tablosu | Elle yazılmış nadir kodon listesi ve bilinen en iyi kodonlar |
| Parça kataloğu | Parçaların gerçek vektörlerde birebir bulunması |

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
  labwork.py           protokol, jel, protein, rapor
  parts.py             katalog yönetimi, iGEM serileri
  literature.py        literatür modülü
core/
  parts_seed.py        60 bakteriyel parça + 8 konak suşu
  igem_seed.py         iGEM kayıtlarını katalog parçasına çevirir
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
  validate_against_references.py  bağımsız doğrulama
  update_checksums.py          veri özetlerini yenile
tests/
  test_sequences.py    dizi motoru
  test_design.py       tasarım modülleri
  test_data.py         veri bütünlüğü ve denetim kuralları
data/
  reference_vectors/   15 vektörün dizisi ve öznitelikleri (JSON)
  igem_parts.json      iGEM parça kayıtları
  ecoli_codon_usage.json  genomdan hesaplanan kodon tablosu
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
