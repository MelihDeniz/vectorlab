# Kullanım amacı ve sınırları

## Yalnızca araştırma ve eğitim amaçlıdır

VectorLab bir **tasarım ve öğrenme aracıdır**. Yalnızca araştırma ve
eğitim amaçlı kullanım için yapılmıştır.

**Tanı, tedavi ya da klinik karar için kullanılamaz.** Bir tıbbi cihaz
değildir; hiçbir düzenleyici kurum (TİTCK, EMA, FDA) tarafından
değerlendirilmemiş ve onaylanmamıştır. İnsan ya da hayvan üzerinde
uygulanacak hiçbir işlemin dayanağı olarak kullanılmamalıdır.

## Uygulamanın ürettiği şey bir öneridir, sonuç değildir

Uygulama dizi hesabı yapar, plazmit haritası çıkarır, kesim ürünlerini
ve jeli benzetir, kodon kullanımını ölçer ve literatür arar. Bunların
hepsi **hesaptır ve modeldir**; laboratuvarda doğrulanması gerekir.

Özellikle:

- **Jel simülasyonu** bir ayrışma modelidir, gerçek bir jel görüntüsü
  değildir. Bant konumları gerçek bir elektroforezde farklı çıkabilir.
- **mRNA saç tokası hesabı** komşu baz etkileşimlerini içeren tam bir
  katlanma hesabı değildir; `README.md` içinde sınırı yazılıdır.
- **Kodon optimizasyonu** ifade düzeyini garanti etmez. CAI tek başına
  ifadeyi öngörmez: promotör, mRNA yapısı, protein katlanması ve konak
  fizyolojisi de belirleyicidir.
- **Denetim kuralları** (CORE, PRO, HOST ve diğerleri) yakaladıkları
  hatayı bildirir; hata vermemesi tasarımın çalışacağı anlamına
  gelmez. Uygulama yalnızca kendisine öğretilen kuralları bilir.
- **Yaklaşık anotasyon** benzerliğe dayanır; bir parçayı "tanıdım"
  demesi o parçanın gerçekten orada olduğunun kanıtı değildir.

Her sayfada hesabın nereden geldiği ve nerede durduğu yazılıdır.
Sayıyı kullanmadan önce o açıklamayı okuyun.

## Biyogüvenlik ve yasal sorumluluk

Uygulama, canlı bakteriyel sistemler için tasarım yapmaya yarar; bazı
sayfalar terapötik amaçlı yapılar (örneğin bağırsakta çalışan canlı
biyoterapötikler) üzerine örnekler içerir. Bunlar **öğretici
örneklerdir**.

Genetiği değiştirilmiş organizmalarla çalışmak Türkiye'de ve her
ülkede yasal izne ve kurumsal onaya bağlıdır. Kullanıcı:

- kendi kurumunun biyogüvenlik kurulundan gerekli onayları almakla,
- yürürlükteki ulusal mevzuata (Biyogüvenlik Kanunu ve ilgili
  yönetmelikler) uymakla,
- çalıştığı organizma ve yapı için uygun biyogüvenlik düzeyinde
  (BSL) çalışmakla,

kendisi yükümlüdür. VectorLab bu yükümlülüklerin hiçbirini
karşılamaz, azaltmaz ve bunlar hakkında hukuki görüş vermez.

## Veri kaynakları ve doğruluk

Diziler ve atıflar kaynaktan çekilir, ezberden yazılmaz; her kaydın
erişim numarası, sürümü ve indirilme tarihi saklanır
(`THIRD-PARTY-NOTICES.md` ve `data/CHECKSUMS.txt`). Buna rağmen:

- Kaynak veritabanları (NCBI, UniProt, iGEM Registry, Europe PMC)
  kayıtlarını **revize edebilir**. Elinizdeki tablo, indirildiği
  tarihteki haldir.
- Bir kaydın kendisi yanlış olabilir. Uygulama kaynağı doğru
  aktarmayı garanti eder, kaynağın doğruluğunu garanti etmez.
- Uygulamanın kendi hesapları `tools/validate_against_references.py`
  ile uygulama dışındaki ölçütlerle karşılaştırılır; geçen kontrol
  sayısı `README.md` içinde yazılıdır. Bu, hatasızlık kanıtı değil,
  hatanın aranmış olduğunun kanıtıdır.

## Garanti yoktur

Kod [MIT lisansıyla](LICENSE) dağıtılır. Lisansın kendi metninde
yazdığı gibi, yazılım **"olduğu gibi"** sunulur; açık ya da örtük
hiçbir garanti verilmez ve yazar, kullanımından doğan hiçbir
zarardan sorumlu tutulamaz.
