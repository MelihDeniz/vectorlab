"""Kullanım amacı beyanı ve üçüncü taraf atıfları.

Metinler tek yerde durur, çünkü aynı beyanın arayüzde, dışa aktarılan
dosyalarda ve basılan raporda birbirinden farklı yazılması, hangisinin
geçerli olduğunu belirsiz bırakır.

Uzun hali `DISCLAIMER.md` dosyasındadır; buradakiler onun arayüze sığan
kısaltmalarıdır ve ona işaret eder.
"""

from __future__ import annotations

#: Her sayfanın altında duran tek satır.
FOOTER = (
    "Yalnızca araştırma ve eğitim amaçlıdır; tanı ya da tedavi için "
    "kullanılamaz. Uygulamanın ürettiği her sonuç laboratuvarda "
    "doğrulanmalıdır."
)

#: Terapötik tasarım yapılan sayfalarda gösterilen uzun hali. Orada
#: kısa bir satır yetmez: o sayfalar insan hastalıklarına yönelik
#: canlı bakteriyel yapılar tasarlıyor.
THERAPEUTIC = (
    "**Bu sayfa öğretici örnekler içerir.** Buradaki yapılar, canlı "
    "biyoterapötik tasarımının nasıl düşünüldüğünü göstermek için "
    "yayımlanmış çalışmalardan türetilmiştir. Hiçbiri klinik kullanıma "
    "hazır bir tasarım değildir ve uygulama bir tıbbi cihaz değildir: "
    "hiçbir düzenleyici kurum tarafından değerlendirilmemiştir.\n\n"
    "Genetiği değiştirilmiş organizmalarla çalışmak yasal izne ve "
    "kurumsal biyogüvenlik onayına bağlıdır. Bu izinleri almak "
    "kullanıcının sorumluluğundadır."
)

#: Dışa aktarılan GenBank/FASTA dosyalarına ve rapora giren satır.
#: Dosya uygulamadan çıkıp başka bir yere gittiğinde beyanın onunla
#: birlikte gitmesi gerekir.
EXPORT = (
    "VectorLab ile uretildi - yalnizca arastirma ve egitim amaclidir; "
    "tani ya da tedavi icin kullanilamaz. Laboratuvarda dogrulanmadan "
    "kullanilmamalidir."
)

#: REBASE'in kendi istedigi atif. rebase.neb.com ana sayfasinda
#: "If you include or otherwise utilize REBASE information in or for
#: your work, please cite the latest review" diye yaziyor; istenen
#: kunye bu. Restriksiyon enzimi verisi Biopython'un Bio.Restriction
#: modulu uzerinden REBASE'ten geldigi icin atif burada duruyor.
REBASE_CITATION = (
    "Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for "
    "DNA restriction and modification: enzymes, genes and genomes. "
    "Nucleic Acids Res. 51: D629-D630 (2023). "
    "doi:10.1093/nar/gkac975 (PMID 36318248)."
)

#: REBASE'in kullanim kosulu: acik veritabani, ticari kullanimi
#: yasaklamiyor, ama turetilen dagitimin da acik kalmasini bekliyor.
#: Kaynak: rebase.neb.com ana sayfasi.
REBASE_TERMS = (
    "REBASE acik bir veritabanidir: \"REBASE is an open database "
    "provided to benefit the scientific community. Adaptation or "
    "distribution of information from REBASE should also remain open "
    "and accessible to the public.\" Bu depo MIT lisansiyla acik "
    "oldugu icin kosul karsilaniyor."
)
