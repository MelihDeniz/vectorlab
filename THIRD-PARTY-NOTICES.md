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

## Kullanılmayan kaynaklar

Addgene ve SnapGene içeriği kullanılmamıştır. Bu veritabanlarının
kayıtları kendi kullanım koşullarına tabidir ve buraya kopyalanamaz.
