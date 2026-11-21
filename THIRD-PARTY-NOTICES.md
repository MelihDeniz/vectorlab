# Üçüncü taraf veriler

Bu depoda VectorLab'in kendi kodunun yanında, açık kaynaklardan çekilmiş
veri dosyaları da bulunur. Kodun lisansı [MIT](LICENSE)'tır; aşağıdaki
dosyalar kendi lisanslarıyla gelir.

## iGEM Registry — CC BY 4.0

- `data/igem_parts.json`
- `data/reference_vectors/PSB1A3.json`
- `data/reference_vectors/PSB1C3.json`
- `data/reference_vectors/PSB1K3.json`
- `data/reference_vectors/PSB4C5.json`

Kaynak: iGEM Registry of Standard Biological Parts, iGEM Foundation —
<https://registry.igem.org>

Lisans: Creative Commons Attribution 4.0 International (CC BY 4.0) —
<https://creativecommons.org/licenses/by/4.0/>

Her parçanın lisansı `tools/fetch_igem_parts.py` tarafından Registry
API'sinden okunur ve kayda yazılır; `cc-by` ile başlamayan bir lisans
görülürse o parça alınmaz. Uygulamada her parçanın kaynak satırında
Registry kimliği ve lisansı gösterilir.

## NCBI Nucleotide (GenBank) — kısıtsız

- `data/reference_vectors/` altındaki diğer tüm dosyalar

Kaynak: NCBI Nucleotide — <https://www.ncbi.nlm.nih.gov/nuccore>

NCBI, veritabanı içeriğinin kullanımına ya da dağıtımına kısıt koymaz.
Her kayıt kendi erişim numarasını, sürümünü ve indirilme tarihini
taşır.

## REBASE — çalışma anında Biopython ile

Restriksiyon enzimi verisi bu depoda saklanmaz; çalışma anında
Biopython'un `Bio.Restriction` modülünden okunur. Veri REBASE'ten
gelir: <http://rebase.neb.com>

REBASE kendi ana sayfasında kullanım koşulunu şöyle yazıyor:

> REBASE is an open database provided to benefit the scientific
> community. Adaptation or distribution of information from REBASE
> should also remain open and accessible to the public.

Yani kullanımı akademik çalışmayla **sınırlı değildir**, ama
REBASE'ten türetilen bilginin dağıtımının da açık kalması bekleniyor.
Bu depo MIT lisansıyla açık olduğu için koşul karşılanıyor. Aynı sayfa
künye de istiyor ("please cite the latest review"); istenen künye:

> Roberts RJ, Vincze T, Posfai J, Macelis D. REBASE: a database for
> DNA restriction and modification: enzymes, genes and genomes.
> *Nucleic Acids Res.* 51: D629-D630 (2023).
> doi:[10.1093/nar/gkac975](https://doi.org/10.1093/nar/gkac975)
> (PMID 36318248)

Künye `core/notice.py` içinde tek yerde durur ve Enzimler sayfasında,
basılan protokol raporunda ve bu dosyada aynı metinle görünür.

## Kod bağımlılıklarının lisansları

Hepsi ticari kullanıma açık, permissive lisanslar:

| Paket | Lisans |
| --- | --- |
| Streamlit | Apache License 2.0 |
| Biopython | Biopython License Agreement (BSD benzeri) |
| Requests | Apache License 2.0 |
| pywebview | BSD 3-Clause |
| pytest (yalnızca geliştirme) | MIT |

## Kullanılmayan kaynaklar

Addgene ve SnapGene içeriği kullanılmamıştır. Bu veritabanlarının
kayıtları kendi kullanım koşullarına tabidir ve buraya kopyalanamaz.
