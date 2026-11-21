"""Hastalık vektörleri: terapötik genler ve hastalık kayıt defteri.

İki şey tutar:

  1. `data/disease_parts.json` içindeki terapötik genleri katalog
     parçasına çevirir. Diziler `tools/fetch_disease_parts.py` ile
     UniProt + NCBI'dan çekilir ve **çevirisi UniProt'un bildirdiği
     proteine eşit olduğu doğrulanmıştır**; ezberden tek baz
     yazılmamıştır.

  2. `DISEASES` kayıt defteri: her hastalık için mekanizma, kaynak
     makale, terapötik kaset ve karşılaştırılacak strateji seçenekleri.

Kayıt defterindeki her atıf Europe PMC'den PMID ve DOI'siyle
doğrulanmıştır. Bu, sunumda sorulacak ilk sorudur ve uydurma atıf
yazmak bu aracın önlemeye çalıştığı hatanın ta kendisidir.

PKU kaydı `core/pku_seed.py` ile birlikte çalışır: o dosya PKU'ya
özgü parçaları (stlA, pheP, dapA, PfnrS, ccdAB) sağlar, burası
hastalığın kendisini ve karşılaştırmayı tarif eder.
"""

from __future__ import annotations

import functools
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

DATA = Path(__file__).resolve().parent.parent / "data" / "disease_parts.json"


# ---------------------------------------------------------------------------
# Terapötik genler
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def _raw() -> dict[str, Any]:
    if not DATA.exists():
        return {"parts": []}
    return json.loads(DATA.read_text(encoding="utf-8"))


#: Her gen için ne işe yaradığı. Dizi kaynaktan gelir, yorum buradan.
GENE_NOTES: dict[str, str] = {
    "CDS_ARGA": (
        "Arginin biyosentezinin ilk ve hız sınırlayıcı adımı. Terapötik "
        "mantık dolaylıdır: bakteri arginin üretirken amonyağı tüketir, "
        "yani kandaki fazla amonyak bağırsakta bakteriyel biyokütleye "
        "çevrilir.\n\n"
        "**Önemli fark:** makaledeki suş yabanıl argA değil, geri "
        "besleme direnci kazandırılmış bir varyant (argA-fbr) kullanır. "
        "Yabanıl enzim kendi ürünü arginin tarafından baskılanır ve "
        "yeterince amonyak tüketmez. Buradaki dizi yabanıl E. coli "
        "argA'sıdır; terapötik sürüm için nokta mutasyonu gerekir "
        "(Tasarımcı'da bölge-yönelimli mutagenez aracı bunun için var)."),
    "CDS_OXC": (
        "Oksalil-KoA'yı formil-KoA ve CO2'ye parçalar. Oksalat "
        "yıkımının iki enziminden biri; tek başına işe yaramaz, frc ile "
        "birlikte bir döngü kurar."),
    "CDS_FRC": (
        "KoA grubunu formattan oksalata aktarır ve oksalil-KoA üretir, "
        "yani oxc'nin substratını hazırlar. İkisi birlikte oksalatı "
        "formata ve CO2'ye çevirir. Oxalobacter formigenes bu çifti "
        "doğal olarak taşır; bağırsakta bu bakteri yoksa oksalat emilimi "
        "artar ve böbrek taşı riski yükselir."),
    "CDS_MDEA": (
        "Metiyonini metantiyol, alfa-ketobütirat ve amonyağa parçalar. "
        "Metiyonin kısıtlaması gereken metabolik hastalıklarda bağırsakta "
        "metiyonin tüketmek için kullanılır."),
    "CDS_IL10": (
        "Anti-inflamatuar sitokin. Bağırsak iltihabında sistemik olarak "
        "verildiğinde yan etkisi yüksektir; fikir, onu bağırsakta "
        "yerinde ve az miktarda üretmektir.\n\n"
        "**İnsan geni olduğuna dikkat edin:** bakteride ifade etmek için "
        "kodon kullanımı uyumsuzdur (Dizi analizi > Kodon optimizasyonu "
        "bunu ölçer) ve salgılanması için bir sinyal peptidi gerekir."),
    "MRK_THYA_LL": (
        "Timidilat sentaz. Terapötik kullanımda markerden çok biyolojik "
        "içerme aracıdır: konağın thyA'sı silinir, gen plazmite konur. "
        "Timin bulunmayan ortamda (bağırsak dışı) plazmitsiz hücre "
        "ölür ve kaçan bakteri çevrede çoğalamaz. Antibiyotik direnci "
        "gerekmez."),
    "CDS_LASR": (
        "Pseudomonas aeruginosa'nın quorum sinyalini (3OC12-HSL) algılayan "
        "düzenleyici. Sinyali bağlayınca hedef promotörü açar. "
        "Terapötik devrenin 'algılayıcı' yarısıdır: ilaç yalnızca patojen "
        "oradayken üretilir."),
    "CDS_PYOS5": (
        "P. aeruginosa'ya özgü bakteriyosin; hedef hücrenin zarında "
        "gözenek açar. Dar spektrumludur, yani bağırsak florasının geri "
        "kalanına dokunmaz - geniş spektrumlu antibiyotiğe göre asıl "
        "üstünlüğü budur.\n\n"
        "Kayıt UniProt'ta TrEMBL (otomatik, elle incelenmemiş) "
        "durumundadır; dizi PAO1 referans genomundan gelir ve çevirisi "
        "doğrulanmıştır."),
    "CDS_CQSA": (
        "Vibrio cholerae'nın CAI-1 quorum sinyalini üreten sentaz. "
        "Terapötik mantık terstir: V. cholerae'da **yüksek** hücre "
        "yoğunluğu sinyali virülansı KAPATIR. Bol CAI-1 üreten zararsız "
        "bir bakteri, patojeni ortamı kalabalık sanmaya zorlar ve "
        "kolera toksini üretimi açılmaz.\n\n"
        "Patojeni öldürmez, yalnızca susturur: seçilim baskısı "
        "yaratmaması bunun üstünlüğü, etkinin geçici olması zayıflığıdır."),
    "CDS_LACZ_LD": (
        "Beta-galaktosidaz: laktozu glukoz ve galaktoza ayırır. Dizi "
        "Lactobacillus delbrueckii'den alındı, E. coli'den değil - "
        "gıda sınıfı bir bakteride kullanılacak genin gıda sınıfı bir "
        "kaynaktan gelmesi düzenleyici açıdan önemlidir.\n\n"
        "1007 amino asit: katalogdaki en büyük terapötik gen. Bu boyut "
        "salgılanmasını zorlaştırır, o yüzden tasarım hücre içi "
        "ifadeye dayanır ve laktozun bakteriye girmesini gerektirir."),
    "CDS_BLAC_BL": (
        "P1A beta-laktamaz. Bağırsağa safrayla sızan beta-laktam "
        "antibiyotiği orada parçalar; kandaki tedaviye dokunmaz, çünkü "
        "enzim bağırsak lümeninde kalır.\n\n"
        "**Dikkat: bu enzim bir direnç genidir.** Terapötik niyetle "
        "kullanılması moleküler gerçeği değiştirmez. Yatay gen "
        "aktarımıyla bir patojene geçerse önlemeye çalıştığı sorunu "
        "yaratır. Kaynak çalışmada bakteri değil, doğrudan enzimin "
        "kendisi veriliyor olması tam olarak bu yüzdendir."),
    "CDS_LSS": (
        "Lizostafin: Staphylococcus aureus hücre duvarındaki "
        "pentaglisin köprülerini kesen glisil-glisin endopeptidaz. "
        "Hedefi çok dardır, bu yüzden cilt ve bağırsak florasının geri "
        "kalanına dokunmaz.\n\n"
        "Dizi ön-pro proteinin tamamı değildir: UniProt'un "
        "anotasyonundaki **olgun enzim bölgesi** (Chain 248-493) "
        "kesilmiştir. Ön-pro bölge S. simulans'ın kendi salgı yoluna "
        "aittir; başka bir bakteride kendi sinyal peptidinizi "
        "kullanırsınız."),
    "CDS_GLP1": (
        "GLP-1(1-37): yemekten sonra salgılanan ve insülin salımını "
        "artıran inkretin peptidi. İnsan proglukagonunun "
        "(UniProt P01275) içindeki bir parçadır; UniProt'un "
        "anotasyonundaki 'Glucagon-like peptide 1' bölgesi (92-128) "
        "kesilerek alınmıştır.\n\n"
        "37 amino asit - katalogdaki en kısa terapötik gen. Peptit "
        "terapötiklerinde asıl sorun ifade değil **kararlılıktır**: "
        "bakteri içinde proteazlara, salgılandıktan sonra DPP-4'e "
        "açıktır."),
    "CDS_ELAFIN": (
        "Elafin: nötrofil elastazını ve proteinaz-3'ü bağlayan küçük "
        "bir inhibitör. İltihaplı bağırsakta doku yıkımının önemli "
        "bir kısmını bu proteazlar yapar; elafin onları durdurur. "
        "IL-10'dan farklı bir mantıktır: bağışıklığı baskılamaz, "
        "hasarı yapan enzimi engeller.\n\n"
        "Olgun bölge (Chain 61-117) kesilmiştir; ön-pro kısmı "
        "insan hücresinin salgı yoluna aittir."),
    "CDS_LSRK": (
        "AI-2 kinaz. AI-2 türler arası ortak bir quorum sinyalidir; "
        "Salmonella, V. cholerae ve EHEC gibi patojenler virülans ve "
        "biyofilm kararlarını buna göre verir. LsrK sinyali fosforile "
        "eder, fosforile AI-2 hücreye çekilir ve ortamdan silinir. "
        "Fazla LsrK üreten bir bakteri komşu patojenin sinyalini "
        "keser.\n\n"
        "**Özgül değildir:** AI-2 ortak bir dildir, onu tüketmek "
        "yararlı komşuların haberleşmesini de bozar."),
    "CDS_PUCL": (
        "Ürat yıkımının ÜÇ basamağından birinci ve üçüncüsünü birden "
        "yapar: ürik asidi 5-hidroksiizoürata (HIU) yükseltger "
        "(N-uç ürat oksidaz alanı) ve zincirin sonundaki OHCU'yu "
        "(S)-allantoine dekarboksiller (C-uç alanı). UniProt kaydın "
        "adını zaten böyle veriyor: *Uric acid degradation "
        "bifunctional protein PucL*.\n\n"
        "**Terapötik mantığın özü:** insanda ürat oksidaz geni "
        "işlevsizdir, o yüzden pürin yıkımı ürik asitte durur ve "
        "çözünürlüğü düşük olan bu molekül birikir. Enzimi bağırsağa "
        "koymak, zincirin insanın kaybettiği basamağını geri "
        "koymaktır.\n\n"
        "**Makaledeki dizi bununla aynı değil.** Kaynak çalışma kodon "
        "optimize edilmiş bir pucL kullanır ve üstüne iki nokta "
        "mutasyonu koyar: Asp-44 yerine Val, Gln-268 yerine Arg "
        "(pucL-M adını verdikleri sürüm). Buradaki dizi yabanıl "
        "*B. subtilis* 168 pucL'idir.\n\n"
        "Bu mutasyonları koordinatına bakıp körü körüne taşımayın - "
        "denetledik: bu dizinin 44. kalıntısı gerçekten Asp, ama "
        "268. kalıntısı zaten Arg, Gln değil (262-274 arası sırayla "
        "L H Y V A H **R** F L D T Y S). Demek ki makalenin "
        "numaralaması bu kayıtla örtüşmüyor ve ikinci mutasyonun yeri "
        "bu dizide başka bir yerdedir. Bir makaleden mutasyon "
        "koordinatı alırken o konumdaki kalıntının beklenen kalıntı "
        "olduğunu doğrulamak zorunludur; Tasarımcı'daki "
        "bölge-yönelimli mutagenez aracı tam bu yüzden değiştirdiği "
        "kodonu ve kalıntıyı size gösterir."),
    "CDS_PUCM": (
        "5-hidroksiizoüratı (HIU) hidrolize edip OHCU'ya çevirir "
        "(EC 3.5.2.17): üç basamaklı ürat yıkım zincirinin "
        "ORTASINDAKİ basamak. En çok atlanan basamak da budur, çünkü "
        "ürikaz koyduk denince genellikle yalnızca pucL'in birinci "
        "işlevi düşünülür.\n\n"
        "HIU kendiliğinden de bozunur - kaynak makale ikinci ve "
        "üçüncü basamağın kendiliğinden yürüyebildiğini, ama "
        "birincisinin vazgeçilmez olduğunu açıkça yazar. Enzimle "
        "yapıldığında zincir hızlı ve tek ürüne gider; enzimsiz "
        "bırakıldığında ara ürünler birikir. 114 amino asitlik küçük "
        "bir proteindir, yani kasete eklenmesinin metabolik bedeli de "
        "düşüktür."),
    "CDS_UACT": (
        "Ürik asidi hücre içine taşıyan, protona bağımlı yüksek "
        "kapasiteli taşıyıcı (UniProt: *Uric acid transporter "
        "UacT*).\n\n"
        "**Katalogda uacT yazıyor, makalede ygfU yazıyor - aynı "
        "gendir:** UniProt Q46821'in gen adı uacT, eş adı ygfU'dur. "
        "Kaynak makaleyi açıp bizim listede ygfU yok diye "
        "düşünmeyin.\n\n"
        "Taşıyıcı olmadan kaset çalışmaz ve bu, enzim listesine "
        "bakarken görünmeyen bir eksiktir: ürat yıkım enzimleri "
        "sitoplazmada çalışır, ürik asit ise bağırsak lümeninde "
        "durur. Zar geçişi sağlanmazsa substrat enzime hiç ulaşmaz. "
        "Terapötik bir kasette enzimi koymak ile enzimin substratına "
        "erişmesi iki ayrı sorundur."),
    "CDS_KATG": (
        "Katalaz-peroksidaz HPI (EC 1.11.1.21). Ürat kasetinde iki iş "
        "yapar ve ikincisi şaşırtıcı olanıdır.\n\n"
        "**Birincisi, zehir temizliği:** ürat oksidaz her ürik asit "
        "için bir hidrojen peroksit (H2O2) üretir. Yıkımı "
        "hızlandırdıkça ürettiğiniz peroksit de artar; "
        "temizlenmezse tasarımın kendisi hücreye zarar verir. Kaynak "
        "makale katalaz eklenince hücrenin ürettiği reaktif oksijen "
        "türlerinin azaldığını ölçmüştür.\n\n"
        "**İkincisi, oksijen geri dönüşümü:** katalazın tepkimesi "
        "2 H2O2 -> O2 + 2 H2O'dur, yani harcanan oksijenin bir "
        "kısmını GERİ VERİR. Kalın bağırsak neredeyse anaerobiktir ve "
        "ürat oksidazın substratlarından biri O2'dir - yani orada "
        "sınırlayıcı olan enzim miktarı değil, oksijendir. Makalenin "
        "oksijen geri dönüşüm sistemi dediği şeyin yarısı budur, "
        "öbür yarısı Vhb'dir."),
    "CDS_VHB": (
        "*Vitreoscilla stercoraria*'nın hemoglobini (Vhb), 146 amino "
        "asitlik küçük bir hem proteini. Biyoteknolojideki klasik "
        "kullanımı şudur: düşük oksijenli kültürde hücreye ifade "
        "edildiğinde oksijen kullanımını iyileştirir. Kaynak makale "
        "onu tam bu nedenle kasete koyar - suşun az oksijende de "
        "ürat yıkabilmesi için.\n\n"
        "**Dizi denetimi:** makalenin kullandığı kayıt NCBI "
        "WP_019959060.1'dir; buradaki UniProt P04252 dizisinin "
        "çevirisi o kayda birebir eşittir (146 aa, karşılaştırıldı). "
        "Yani bu parça makaledeki proteinin aynısıdır.\n\n"
        "**Dürüst olmak gerekirse:** UniProt bu proteinin işlevini "
        "terminal oksidaz olarak çalışır diye anotlar; bu, oksijen "
        "deposu demekle aynı şey değildir. Vhb'nin hücre içi oksijen "
        "erişimini nasıl iyileştirdiği literatürde hâlâ tartışılır. "
        "Onu kasete makale koyduğu için koyuyoruz, mekanizması "
        "kesinleşmiş olduğu için değil."),
    "CDS_PEP_MX": (
        "*Myxococcus xanthus*'un prolil endopeptidazı (prolil "
        "oligopeptidaz, EC 3.4.21.26): peptit bağını içteki bir "
        "prolinin karboksil tarafından keser.\n\n"
        "**Neden tam olarak bu kesim:** glutenin gliadin kısmı prolin "
        "ve glutamin bakımından çok zengindir, insanın mide, pankreas "
        "ve fırça kenarı proteazları ise prolinden sonra verimli "
        "kesemez. Bu yüzden sindirimden uzun, prolinden zengin "
        "peptitler sağ çıkar; en iyi çalışılmışı alfa-2 gliadinin "
        "33 amino asitlik parçasıdır (33-mer). Prolil endopeptidaz "
        "insanın kesemediği bağı keser, yani eksik olan enzim "
        "etkinliğini yerine koyar.\n\n"
        "Kayıt TrEMBL durumundadır (otomatik üretilmiş, elle "
        "incelenmemiş). Dizinin çevirisinin UniProt'un bildirdiği "
        "proteine eşit olduğu doğrulanmıştır, ama kaydın küratörlüğü "
        "yapılmamıştır."),
    "CDS_TFF1": (
        "İnsan trefoil faktör 1'inin (hTFF1) olgun kısmı. Trefoil "
        "peptitleri mukus jelinin içinde duran, proteaza dirençli "
        "küçük peptitlerdir; UniProt işlevi şöyle özetler: sindirim "
        "sistemi mukozasının üstündeki mukus jelini kararlı kılan ve "
        "zarar verici etkenlere karşı fiziksel bariyer sağlayan "
        "protein.\n\n"
        "Yara iyileşmesindeki rolü restitüsyondur: epitel hücreleri "
        "yaranın üstüne göç eder ve bariyer, hücre bölünmesini "
        "beklemeden kapanır. Işın ya da kemoterapi sonrası bölünmenin "
        "durduğu bir dokuda işe yarayabilecek mekanizma tam olarak "
        "budur.\n\n"
        "**Katlanma uyarısı:** 61 amino asitlik bu parça ÜÇ disülfür "
        "bağı taşır (UniProt: 31-57, 41-56, 51-68 - tam kayıt "
        "numaralamasıyla). Disülfür bağı sitoplazmanın indirgen "
        "ortamında kurulmaz; bu yüzden E. coli'de ön üretim "
        "yapılacaksa periplazmaya gönderilmesi gerekir ve kaset "
        "DsbA sinyaliyle kurulur. Elafinde karşılaşılan sorunun "
        "aynısıdır.\n\n"
        "**İnsan geni olduğuna dikkat edin:** bu araç bakteri "
        "vektörleri için yapıldı; insan dizileri yalnızca "
        "salgılanacak terapötik yük olarak bulunur (IL-10, elafin, "
        "GLP-1, hTFF1) ve bu istisna belgelenmiştir. Vektörün hiçbir "
        "parçası - orijin, promotör, marker, iskelet - memeli "
        "kaynaklı değildir."),
}


@functools.lru_cache(maxsize=1)
def disease_parts() -> tuple[dict[str, Any], ...]:
    out: list[dict[str, Any]] = []
    for rec in _raw().get("parts", []):
        incelenmis = rec.get("reviewed", True)
        notlar = [GENE_NOTES.get(rec["part_id"], "")]
        bolge = rec.get("region")
        if bolge:
            notlar.append(
                f"Dizi {rec['organism']} kaynaklıdır: UniProt "
                f"{rec['uniprot']} ({rec['uniprot_entry']}), nükleotid "
                f"{rec['nucleotide']}. Proteinin TAMAMI değil, "
                f"UniProt'un öznitelik tablosundaki **{bolge}** bölgesi "
                f"kullanıldı: {rec['region_full_aa']} amino asitlik "
                f"kayıttan {rec['length_aa']} amino asit / "
                f"{rec['length']} bp. Tam kaydın çevirisi önce "
                f"doğrulandı, kesilen parçanın çevirisi de ayrıca "
                f"sınandı."
                + (" Kesilen bölge metiyoninle başlamadığı için "
                   "**başlangıç kodonu (ATG) bizim eklediğimizdir**; "
                   "sonundaki stop kodonu da öyle. Bölge sınırları "
                   "kayıttan, bu iki kodon bizden."
                   if rec.get("start_codon_added") else
                   " Sonundaki stop kodonu bizim eklediğimizdir.")
                + ("" if incelenmis else
                   " Kayıt TrEMBL durumundadır, yani otomatik üretilmiş "
                   "ve elle incelenmemiştir."))
        else:
            notlar.append(
                f"Dizi {rec['organism']} kaynaklıdır: UniProt "
                f"{rec['uniprot']} ({rec['uniprot_entry']}), nükleotid "
                f"{rec['nucleotide']}. {rec['length_aa']} amino asit / "
                f"{rec['length']} bp; çevirisinin UniProt'taki proteine "
                f"eşit olduğu doğrulandı."
                + ("" if incelenmis else
                   " Kayıt TrEMBL durumundadır, yani otomatik üretilmiş "
                   "ve elle incelenmemiştir."))
        out.append({
            "part_id": rec["part_id"],
            "name": rec["name"],
            "category": rec["category"],
            "sequence": rec["sequence"],
            "seq_status": "full",
            "length": rec["length"],
            "strength": None, "induction": None,
            "requires": [], "incompatible_with": [],
            "copy_number": None, "inc_group": None,
            "notes": "\n\n".join(x for x in notlar if x),
            "accession": (f"UniProt {rec['uniprot']} · "
                          f"{rec['nucleotide']}"),
            "pmid": "", "source": "disease",
        })
    return tuple(out)


# ---------------------------------------------------------------------------
# Hastalık kayıt defteri
# ---------------------------------------------------------------------------


@dataclass
class Strategy:
    key: str
    title: str
    parts: list[str]
    host: str
    summary: str


@dataclass
class Disease:
    key: str
    name: str
    short: str
    mechanism: str
    logic: str                 # terapötik fikrin tek cümlelik özeti
    papers: list[str]
    pathway: str               # graphviz kaynağı
    pathway_caption: str
    cassette: list[str]
    strategies: list[Strategy]
    comparison: list[dict[str, str]] = field(default_factory=list)
    gap: str = ""
    biomarker: str = ""
    purpose: str = "ifade"


def _gv(*satirlar: str) -> str:
    govde = "\n        ".join(satirlar)
    return ("digraph {\n        rankdir=LR;\n"
            '        node [shape=box, style="rounded,filled", '
            'fillcolor="#F4F4F5", fontname="Helvetica"];\n        '
            + govde + "\n    }")


PKU_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_STLA", "RBS_B0032", "CDS_PHEP",
             "TER_B0015"]
AMONYAK_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_ARGA", "TER_B0015"]
OKSALAT_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_OXC", "RBS_B0032", "CDS_FRC",
                 "TER_B0015"]
IBD_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
             "CDS_IL10", "TER_B0015"]
PSA_KASET = ["PRO_J23105", "RBS_B0032", "CDS_LASR", "TER_B0015",
             "IGEM_R0062", "RBS_B0034", "CDS_PYOS5", "IGEM_B0015"]


DISEASES: list[Disease] = [
    # -------------------------------------------------------------- PKU
    Disease(
        key="pku",
        name="Fenilketonüri (PKU)",
        short="PKU",
        logic="Bağırsakta fenilalanini tüketen bir probiyotik.",
        mechanism=(
            "Fenilalanin hidroksilaz (PAH) eksikliği yüzünden fenilalanin "
            "(Phe) vücutta birikir ve sinir sistemine zarar verir. Standart "
            "tedavi Phe'den kısıtlı diyettir ve uyumu zordur. Makaledeki "
            "fikir: bağırsakta Phe'yi tüketen bir probiyotik, yani "
            "**bağırsakta çalışan yaşayan bir ilaç**."),
        papers=[
            "Isabella VM ve ark. Development of a synthetic live bacterial "
            "therapeutic for the human metabolic disease phenylketonuria. "
            "Nature Biotechnology 36:857-864 (2018). doi:10.1038/nbt.4222 "
            "· PMID 30102294",
            "Adolfsen KJ ve ark. Improvement of a synthetic live bacterial "
            "therapeutic for phenylketonuria with biosensor-enabled "
            "enzyme engineering. Nature Communications 12:6215 (2021). "
            "doi:10.1038/s41467-021-26524-0 · PMID 34711827",
        ],
        pathway=_gv(
            'a [label="Diyetten gelen\\nfenilalanin (Phe)"];',
            'b [label="PheP taşıyıcısı\\nPhe\'yi hücreye alır"];',
            'c [label="PAL (stlA)\\nPhe → trans-sinnamat"];',
            'd [label="Karaciğerde\\nhippurata dönüşür"];',
            'e [label="İdrarla atılır"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "Tüm genler PfnrS promotörü altındadır, yani yalnızca "
            "bağırsağın oksijensiz ortamında açılır."),
        cassette=PKU_KASET,
        biomarker="İdrarda hippurat, kanda trans-sinnamat",
        gap=("Nissle insan bağırsağını kalıcı kolonize etmez, bu yüzden "
             "tedavi kalıcı değildir ve doz tekrarlanır."),
        strategies=[
            Strategy("D", "Toksin-antitoksin (ccdAB) idamesi, içerme korunur",
                     PKU_KASET + ["ORI_P15A", "PRO_CCDAB", "MRK_CCDA",
                                  "MRK_CCDB"], "ECN_DAPA",
                     "Plazmiti kaybeden hücrede kararsız antitoksin "
                     "yıkılır, toksin kalır ve hücre ölür; antibiyotik "
                     "gerekmez. Konak ΔdapA kaldığı için makaledeki DAP "
                     "bağımlılığı bozulmaz: iki ayrı güvenlik katmanı."),
            Strategy("A", "dapA tamamlama (içerme bedelli)",
                     PKU_KASET + ["ORI_P15A", "MRK_DAPA"], "ECN_DAPA",
                     "Plazmit dapA'yı geri verir; plazmitsiz hücre ölür. "
                     "Bedeli: plazmitli hücre DAP'sız da büyür, yani "
                     "biyolojik içerme kalkar (BIO-01)."),
            Strategy("B", "Kanamisin markerlı laboratuvar prototipi",
                     PKU_KASET + ["ORI_P15A", "MRK_KAN"], "ECN",
                     "Klonlaması kolay, ama hastaya verilecek bakteride "
                     "direnç geni kabul görmez (BIO-02)."),
            Strategy("C", "İdame yok (yalnızca DAP bağımlılığı)",
                     PKU_KASET + ["ORI_P15A"], "ECN_DAPA",
                     "İçerme korunur ama plazmit idamesi için yöntem yok: "
                     "antibiyotiksiz pasajlarda plazmit kaybolur.\n\n"
                     "**Denetim burada CORE-02 hatası verir ve bu "
                     "beklenen:** markersiz bir plazmitte transformasyondan "
                     "sonra doğru koloniyi seçmenin yolu bile yoktur. "
                     "Basamak, öbür üçünün ne işe yaradığını "
                     "göstermek için duruyor."),
        ],
        comparison=[
            {"Özellik": "Kasa organizma",
             "Makale": "E. coli Nissle 1917",
             "Bizim tasarım": "Aynı (Nissle 1917 ΔdapA)",
             "Etiket": "Literatür"},
            {"Özellik": "Genetik yapı",
             "Makale": "Kromozoma entegre: 2 kopya pheP, 3 kopya stlA",
             "Bizim tasarım": "Plazmitte tek operon + ccdAB idame kaseti",
             "Etiket": "Öneri"},
            {"Özellik": "Enzim ve taşıyıcı",
             "Makale": "PAL (stlA) + LAAD (pma) + PheP",
             "Bizim tasarım": "PAL + PheP; LAAD kapsam dışı",
             "Etiket": "Öneri"},
            {"Özellik": "RBS tekrarı",
             "Makale": "—",
             "Bizim tasarım": "İki gen farklı RBS kullanır (B0034/B0032): "
                              "aynı RBS iki kez geçerse tekrar dizisi "
                              "rekombinasyona açık olur",
             "Etiket": "Öneri"},
        ],
    ),
    # ----------------------------------------------------- hiperammonemi
    Disease(
        key="amonyak",
        name="Hiperammonemi (üre döngüsü bozukluğu / karaciğer yetmezliği)",
        short="Hiperammonemi",
        logic="Amonyağı bağırsakta arginine çeviren bir probiyotik.",
        mechanism=(
            "Üre döngüsü enzimlerinden biri eksik olduğunda ya da karaciğer "
            "yetmezliğinde amonyak kanda birikir; beyne geçer ve "
            "ensefalopatiye yol açar. Fikir dolaylıdır: bakteriye arginin "
            "ürettirmek. Arginin biyosentezi azot tüketir, yani bakteri "
            "çoğalırken bağırsaktaki amonyağı kendi biyokütlesine çevirir "
            "ve dışkıyla atılır."),
        papers=[
            "Kurtz CB ve ark. An engineered E. coli Nissle improves "
            "hyperammonemia and survival in mice and shows dose-dependent "
            "responses in healthy humans. Science Translational Medicine "
            "11:eaau7975 (2019). doi:10.1126/scitranslmed.aau7975 "
            "· PMID 30651324",
        ],
        pathway=_gv(
            'a [label="Bağırsakta\\nfazla amonyak (NH3)"];',
            'b [label="ArgA (geri besleme\\ndirençli)"];',
            'c [label="Arginin biyosentezi\\nazotu tüketir"];',
            'd [label="Amonyak bakteri\\nbiyokütlesine girer"];',
            'e [label="Dışkıyla atılır"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "Buradaki terapötik etki bir enzimin substratı yok etmesi "
            "değil, bir biyosentez yolunun azotu çekmesidir - PKU'dan "
            "mantıkça farklıdır."),
        cassette=AMONYAK_KASET,
        biomarker="Kanda amonyak, idrarda azot bileşikleri",
        gap=("Katalogdaki argA YABANIL tiptir. Makaledeki suş geri besleme "
             "direnci kazandırılmış bir varyant (argA-fbr) kullanır; "
             "yabanıl enzim kendi ürünü arginin tarafından baskılandığı "
             "için yeterince amonyak tüketmez. Bu fark tasarımın en "
             "kritik noktasıdır ve Tasarımcı'daki mutagenez aracıyla "
             "kapatılmalıdır."),
        strategies=[
            Strategy("A", "ccdAB idamesi, ΔdapA içerme korunur",
                     AMONYAK_KASET + ["ORI_P15A", "PRO_CCDAB", "MRK_CCDA",
                                      "MRK_CCDB"], "ECN_DAPA",
                     "PKU'daki D seçeneğinin aynısı: antibiyotiksiz "
                     "idame + DAP bağımlılığı. Terapötik kaset tek gen "
                     "olduğu için yapı PKU'dan daha küçüktür."),
            Strategy("B", "dapA tamamlama",
                     AMONYAK_KASET + ["ORI_P15A", "MRK_DAPA"], "ECN_DAPA",
                     "Antibiyotiksiz ama içerme bedelli (BIO-01)."),
            Strategy("C", "Kanamisin markerlı prototip",
                     AMONYAK_KASET + ["ORI_P15A", "MRK_KAN"], "ECN",
                     "Yalnızca laboratuvar karşılaştırması için."),
        ],
        comparison=[
            {"Özellik": "Kasa organizma", "Makale": "E. coli Nissle 1917",
             "Bizim tasarım": "Aynı (ΔdapA)", "Etiket": "Literatür"},
            {"Özellik": "Terapötik gen",
             "Makale": "argA-fbr (geri besleme dirençli)",
             "Bizim tasarım": "argA yabanıl — mutagenez gerekiyor",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Promotör", "Makale": "PfnrS (anaerobik)",
             "Bizim tasarım": "PfnrS", "Etiket": "Literatür"},
            {"Özellik": "Yerleşim",
             "Makale": "Kromozoma entegre",
             "Bizim tasarım": "Plazmit (ayar kolaylığı için)",
             "Etiket": "Öneri"},
        ],
    ),
    # ------------------------------------------------------ hiperoksalüri
    Disease(
        key="oksalat",
        name="Hiperoksalüri (böbrek taşı)",
        short="Hiperoksalüri",
        logic="Bağırsakta oksalatı parçalayan iki enzimli bir döngü.",
        mechanism=(
            "Oksalat idrarda kalsiyumla birleşip böbrek taşı oluşturur. "
            "Sağlıklı bağırsakta Oxalobacter formigenes oksalatı parçalar; "
            "bu bakteri antibiyotik kullanımıyla kaybolunca oksalat emilimi "
            "artar. Fikir: oksalat yıkım enzimlerini taşıyan bir probiyotik "
            "vermek, böylece oksalat bağırsakta parçalansın ve kana "
            "geçmesin."),
        papers=[
            "Lubkowicz D ve ark. An engineered bacterial therapeutic lowers "
            "urinary oxalate in preclinical models and in healthy "
            "volunteers. Molecular Systems Biology 18:e10539 (2022). "
            "doi:10.15252/msb.202110539 · PMID 35253995",
        ],
        pathway=_gv(
            'a [label="Diyetten gelen\\noksalat"];',
            'b [label="Frc\\nformil-KoA + oksalat\\n→ oksalil-KoA"];',
            'c [label="Oxc\\noksalil-KoA →\\nformil-KoA + CO2"];',
            'd [label="Formatla\\natılır"];',
            "a -> b -> c -> d;",
            "c -> b [label=\"KoA geri döner\", style=dashed];"),
        pathway_caption=(
            "İki enzim bir DÖNGÜ kurar: Oxc'nin ürettiği formil-KoA, "
            "Frc'nin substratı olarak geri döner. Bu yüzden ikisinden "
            "biri eksikse yol hiç çalışmaz - tek gen yeterli değildir."),
        cassette=OKSALAT_KASET,
        biomarker="İdrarda oksalat, dışkıda format",
        gap=("Oksalatın bakteriye girmesi için bir taşıyıcı gerekir. "
             "Oxalobacter'in kendi oksalat/format antiporteri (OxlT) "
             "katalogda yok; onsuz enzimler hücre içinde substrat "
             "bulamayabilir. PKU'da PheP'in oynadığı rolün karşılığıdır "
             "ve eklenmesi gerekir."),
        strategies=[
            Strategy("A", "İki enzim tek operonda, ccdAB idamesi",
                     OKSALAT_KASET + ["ORI_P15A", "PRO_CCDAB", "MRK_CCDA",
                                      "MRK_CCDB"], "ECN_DAPA",
                     "oxc ve frc tek promotörden, farklı RBS'lerle "
                     "okunur. İki genin oranı önemlidir: döngü olduğu "
                     "için biri diğerini sınırlar."),
            Strategy("B", "dapA tamamlama",
                     OKSALAT_KASET + ["ORI_P15A", "MRK_DAPA"], "ECN_DAPA",
                     "Antibiyotiksiz, içerme bedelli."),
            Strategy("C", "Yalnızca oxc (eksik tasarım - karşılaştırma)",
                     ["PRO_FNRS", "RBS_B0034", "CDS_OXC", "TER_B0015",
                      "ORI_P15A", "MRK_KAN"], "ECN",
                     "Bilerek eksik bırakılmış seçenek: frc olmadan "
                     "oksalil-KoA üretilemez, yani oxc substratsız kalır. "
                     "Sunumda 'neden iki gen gerekiyor' sorusunun "
                     "cevabıdır."),
        ],
        comparison=[
            {"Özellik": "Kasa organizma", "Makale": "E. coli Nissle 1917",
             "Bizim tasarım": "Aynı", "Etiket": "Literatür"},
            {"Özellik": "Terapötik genler",
             "Makale": "Oxalobacter formigenes oxc + frc",
             "Bizim tasarım": "Aynı ikili", "Etiket": "Literatür"},
            {"Özellik": "Oksalat taşıyıcısı",
             "Makale": "Oksalat alımı için taşıyıcı",
             "Bizim tasarım": "Katalogda yok — açık",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Gen oranı",
             "Makale": "—",
             "Bizim tasarım": "Farklı RBS'lerle oran ayarlanabilir "
                              "(Kütüphane kur ile taranabilir)",
             "Etiket": "Öneri"},
        ],
    ),
    # -------------------------------------------------------------- IBD
    Disease(
        key="ibd",
        name="İnflamatuar bağırsak hastalığı (kolit)",
        short="IBD / kolit",
        logic="Bağırsakta yerinde IL-10 salgılayan gıda sınıfı bir bakteri.",
        mechanism=(
            "IL-10 iltihabı baskılayan bir sitokindir, ama damardan "
            "verildiğinde sistemik yan etkileri ağırdır. Fikir: onu "
            "iltihabın olduğu yerde, bağırsakta, az miktarda üretmek. "
            "Steidler ve ark. bunu *Lactococcus lactis* ile yaptı - gıda "
            "sınıfı (GRAS) bir bakteri, proteini doğrudan besiyerine "
            "salgılar.\n\n"
            "Bu, kataloğun geri kalanından **mimari olarak farklı** bir "
            "tasarımdır: kasa Gram-pozitif, ürün bir enzim değil "
            "salgılanan bir insan proteini, içerme antibiyotikle değil "
            "timin bağımlılığıyla sağlanır."),
        papers=[
            "Steidler L ve ark. Treatment of murine colitis by Lactococcus "
            "lactis secreting interleukin-10. Science 289:1352-1355 (2000). "
            "doi:10.1126/science.289.5483.1352 · PMID 10958782",
            "Steidler L ve ark. Biological containment of genetically "
            "modified Lactococcus lactis for intestinal delivery of human "
            "interleukin 10. Nature Biotechnology 21:785-789 (2003). "
            "doi:10.1038/nbt840 · PMID 12808464",
        ],
        pathway=_gv(
            'a [label="L. lactis\\nbağırsakta"];',
            'b [label="Usp45 sinyali\\nIL-10\'u dışarı taşır"];',
            'c [label="İltihaplı dokuda\\nyerel IL-10"];',
            'd [label="İltihap baskılanır"];',
            'e [label="thyA eksikliği:\\ntiminsiz ortamda ölür", '
            'fillcolor="#FEE2E2"];',
            "a -> b -> c -> d;",
            "a -> e [style=dashed];"),
        pathway_caption=(
            "İki ayrı mekanizma: Usp45 sinyali ürünü dışarı taşır, thyA "
            "delesyonu ise bakterinin bağırsak dışında yaşamasını "
            "engeller. İkincisi 2003 makalesinin katkısıdır."),
        cassette=IBD_KASET,
        purpose="ifade",
        biomarker="Dışkıda kalprotektin, doku histolojisi",
        gap=("Katalogda Gram-pozitif bir replikasyon orijini parça olarak "
             "yok. Yapı bu hâliyle E. coli'de çoğalır ama L. lactis'te "
             "çoğalmaz - denetim bunu ORG-01 ile yakalar. Doğru yol, "
             "yapıyı pNZ8048 ya da pWV01 iskeletine takmaktır: ikisi de "
             "Referans vektörler sayfasında var ve Referanstan "
             "kes-yapıştır sayfasında bu iş yapılabilir."),
        strategies=[
            Strategy("A", "thyA içermesi (makaledeki yol)",
                     IBD_KASET + ["MRK_THYA_LL"], "LLNZ9000",
                     "Konağın thyA'sı silinir, gen plazmite konur. "
                     "Timin bulunmayan ortamda plazmitsiz hücre ölür ve "
                     "kaçan bakteri çevrede çoğalamaz. Antibiyotik "
                     "direnci taşımaz - terapötik kullanımın şartı budur.\n\n"
                     "**Denetim burada CORE-01 hatası verir ve bu "
                     "beklenen:** katalogda L. lactis'te çoğalan bir "
                     "replikon yok. Eksik olan iskelettir, kaset değil; "
                     "Referanstan kes-yapıştır sayfasındaki pNZ8048 "
                     "tarifi onu verir."),
            Strategy("B", "Kloramfenikol markerlı prototip",
                     IBD_KASET + ["REF_CAT_651"], "LLNZ9000",
                     "pNZ8048'den çıkarılmış cat markeri. Laboratuvarda "
                     "kurması kolay, hastada kullanılamaz.\n\n"
                     "**Denetim burada CORE-01 hatası verir ve bu "
                     "beklenen:** katalogda L. lactis'te çoğalan bir "
                     "replikon yok. Eksik olan iskelettir, kaset değil; "
                     "Referanstan kes-yapıştır sayfasındaki pNZ8048 "
                     "tarifi onu verir."),
            Strategy("C", "E. coli'de ön üretim (mekik yaklaşımı)",
                     IBD_KASET + ["ORI_P15A", "MRK_CAT"], "DH5A",
                     "Klonlama ve doğrulama E. coli'de yapılır, sonra "
                     "Gram-pozitif replikon eklenip L. lactis'e aktarılır. "
                     "Denetim burada promotörün E. coli'ye uymadığını "
                     "söyleyecek - bu beklenen bir uyarıdır."),
        ],
        comparison=[
            {"Özellik": "Kasa organizma",
             "Makale": "Lactococcus lactis (gıda sınıfı)",
             "Bizim tasarım": "Aynı (NZ9000)", "Etiket": "Literatür"},
            {"Özellik": "Ürün", "Makale": "İnsan IL-10, salgılanan",
             "Bizim tasarım": "Aynı; Usp45 sinyaliyle",
             "Etiket": "Literatür"},
            {"Özellik": "İçerme",
             "Makale": "thyA delesyonu (2003 makalesi)",
             "Bizim tasarım": "thyA tamamlama plazmitte",
             "Etiket": "Literatür"},
            {"Özellik": "Kodon kullanımı",
             "Makale": "—",
             "Bizim tasarım": "IL-10 insan genidir; bakteride kodon "
                              "uyumsuzluğu ölçülmeli (Kodon "
                              "optimizasyonu sekmesi)",
             "Etiket": "Öneri"},
        ],
    ),
    # ---------------------------------------------- P. aeruginosa enfeksiyonu
    Disease(
        key="psa",
        name="Pseudomonas aeruginosa enfeksiyonu",
        short="P. aeruginosa",
        logic="Patojeni kendi sinyalinden tanıyıp yalnızca o varken öldüren "
              "bir devre.",
        mechanism=(
            "P. aeruginosa hastane enfeksiyonlarının en inatçı "
            "etkenlerinden biridir ve çok ilaca dirençlidir. Bakteriler "
            "birbirleriyle küçük moleküllerle haberleşir (quorum sensing); "
            "P. aeruginosa 3OC12-HSL salgılar.\n\n"
            "Fikir: bu sinyali duyan bir bakteri yapmak. LasR algılayıcısı "
            "sinyali bağlayınca aşağı akıştaki promotörü açar ve hücre "
            "yalnızca o zaman öldürücü protein üretir. Geniş spektrumlu "
            "antibiyotikten farkı budur: ilaç patojen oradayken ve "
            "yalnızca ona karşı üretilir."),
        papers=[
            "Saeidi N ve ark. Engineering microbes to sense and eradicate "
            "Pseudomonas aeruginosa, a human pathogen. Molecular Systems "
            "Biology 7:521 (2011). doi:10.1038/msb.2011.55 · PMID 21847113",
        ],
        pathway=_gv(
            'a [label="P. aeruginosa\\n3OC12-HSL salgılar", '
            'fillcolor="#FEE2E2"];',
            'b [label="LasR algılayıcı\\nsinyali bağlar"];',
            'c [label="Hedef promotör\\naçılır"];',
            'd [label="Pyocin S5\\nüretilir"];',
            'e [label="P. aeruginosa\\nzarında gözenek", '
            'fillcolor="#FEE2E2"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "Bu bir ALGILAYICI devredir, metabolik bir yol değil: iki "
            "kaset var, birincisi algılayıcıyı sürekli üretir, ikincisi "
            "yalnızca sinyal varken açılır. Katalogdaki öteki "
            "hastalıklardan tasarım deseni olarak ayrılır."),
        cassette=PSA_KASET,
        biomarker="Kültürde P. aeruginosa canlılığı, öldürme eğrisi",
        gap=("Öldürücü proteinin hücre dışına çıkması gerekir. Makaledeki "
             "devre bunun için ayrıca bir liziz proteini (E7) kullanır; "
             "katalogda yok. Ayrıca LasR'nin açtığı promotör olarak "
             "pLux (BBa_R0062) kullanılmıştır: LasR ile LuxR'nin hedef "
             "promotörleri çapraz tanınır, ama bu bir YAKLAŞIKLIKTIR ve "
             "deneyle doğrulanmalıdır."),
        strategies=[
            Strategy("A", "Algılayıcı + öldürücü, kanamisin markerlı",
                     PSA_KASET + ["ORI_P15A", "MRK_KAN"], "DH5A",
                     "İki kasetli devre: LasR sabit promotörden sürekli "
                     "üretilir, pyocin ise yalnızca sinyal varken. "
                     "Laboratuvar prototipi."),
            Strategy("B", "Yalnızca algılayıcı (raportörlü ön test)",
                     ["PRO_J23105", "RBS_B0032", "CDS_LASR", "TER_B0015",
                      "IGEM_R0062", "RBS_B0034", "RPT_SFGFP", "IGEM_B0015",
                      "ORI_P15A", "MRK_KAN"], "DH5A",
                     "Öldürücü yerine GFP konur. Devrenin sinyale "
                     "gerçekten yanıt verip vermediğini ölçmenin ilk "
                     "adımı budur; öldürme denemesinden önce yapılmalı."),
            Strategy("C", "Sürekli öldürücü (kontrol - kötü tasarım)",
                     ["PRO_J23105", "RBS_B0034", "CDS_PYOS5", "TER_B0015",
                      "ORI_P15A", "MRK_KAN"], "DH5A",
                     "Algılayıcı yok, pyocin sürekli üretilir. "
                     "Karşılaştırma için: devrenin asıl değeri "
                     "seçiciliktir, bu seçenek onu kaybeder."),
        ],
        comparison=[
            {"Özellik": "Kasa organizma", "Makale": "E. coli",
             "Bizim tasarım": "E. coli (DH5-alfa prototip)",
             "Etiket": "Literatür"},
            {"Özellik": "Algılayıcı", "Makale": "LasR + 3OC12-HSL",
             "Bizim tasarım": "LasR", "Etiket": "Literatür"},
            {"Özellik": "Öldürücü", "Makale": "Pyocin S5 + E7 liziz proteini",
             "Bizim tasarım": "Pyocin S5; liziz proteini kapsam dışı",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Hedef promotör",
             "Makale": "LasR'nin hedef promotörü",
             "Bizim tasarım": "pLux (BBa_R0062) — çapraz tanıma "
                              "varsayımı, doğrulanmalı",
             "Etiket": "Öneri"},
        ],
    ),
]


def by_key(key: str) -> Disease | None:
    return next((d for d in DISEASES if d.key == key), None)


def keys() -> list[str]:
    return [d.key for d in DISEASES]


HCU_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_MDEA", "TER_B0015"]
KOLERA_KASET = ["PRO_J23105", "RBS_B0034", "CDS_CQSA", "TER_B0015"]
CDIFF_KASET = ["PRO_J23100", "RBS_B0034", "SIG_AMYE", "CDS_BLAC_BL",
               "TER_B0015"]
#: Oksidatif bagirsak hasari: SOD superoksidi H2O2'ye, katalaz
#: H2O2'yi suya cevirir. Ikisi arka arkaya calisir - yalnizca SOD
#: koymak superoksidi H2O2'de BIRAKIR, ki o da zarar vericidir.
SOD_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034",
             "CDS_SODA", "CDS_KATE", "TER_B0015"]

LAKTOZ_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034",
                "CDS_LACZ_LD", "TER_B0015"]
MRSA_KASET = ["PRO_J23105", "RBS_B0034", "SIG_USP45", "CDS_LSS",
              "TER_B0015"]
GLP1_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
              "CDS_GLP1", "TER_B0015"]
AI2_KASET = ["PRO_J23100", "RBS_B0034", "CDS_LSRK", "TER_B0015"]
ELAFIN_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
                "CDS_ELAFIN", "TER_B0015"]

# İdame ve içerme merdiveni. Aynı dört soru her hastalıkta sorulur:
# plazmit nasıl tutulur, içerme korunur mu, antibiyotik var mı, hiç
# idame olmasa ne olur. Karşılaştırılabilir olması için kasıtlı olarak
# aynı sırada durur; denetim her birine farklı şey söyler.
_CCDAB = ["PRO_CCDAB", "MRK_CCDA", "MRK_CCDB"]


def _ecoli_merdiven(kaset: list[str], ori: str, ne: str) -> list[Strategy]:
    """E. coli konaklı hastalıklar için dört basamaklı merdiven."""
    return [
        Strategy("A", "Toksin-antitoksin (ccdAB) idamesi, içerme korunur",
                 kaset + [ori] + _CCDAB, "ECN_DAPA",
                 f"{ne} Plazmiti kaybeden hücrede kararsız antitoksin "
                 "yıkılır, toksin kalır ve hücre ölür; antibiyotiğe "
                 "gerek kalmaz. Konak ΔdapA kaldığı için makalelerdeki "
                 "DAP bağımlılığı da bozulmaz: iki ayrı güvenlik "
                 "katmanı üst üste durur."),
        Strategy("B", "dapA tamamlama (içerme bedelli)",
                 kaset + [ori, "MRK_DAPA"], "ECN_DAPA",
                 "Plazmit konağın eksik dapA'sını tamamlar, yani "
                 "plazmitini kaybeden hücre DAP'sız ortamda çoğalamaz. "
                 "Bedeli şudur: dapA artık plazmitte olduğu için "
                 "bakteri bağırsak dışında da yaşayabilir hale gelir "
                 "ve biyolojik içerme kalkar. Denetim bunu BIO-01 ile "
                 "söyler."),
        Strategy("C", "Kanamisin markerlı laboratuvar prototipi",
                 kaset + [ori, "MRK_KAN"], "ECN",
                 "Klonlaması ve doğrulaması en kolay sürüm; "
                 "laboratuvarda ilk kuracağınız budur. Hastaya "
                 "verilecek bir bakteride direnç geni kabul görmez, "
                 "çünkü bağırsak mikrobiyotasına yatay gen "
                 "aktarımıyla geçebilir. Denetim BIO-02 ile uyarır."),
        Strategy("D", "İdame yok (yalnızca konak bağımlılığı)",
                 kaset + [ori], "ECN_DAPA",
                 "Karşılaştırma için: plazmiti tutan hiçbir seçilim "
                 "yok. Antibiyotiksiz pasajlarda plazmit kaybolur ve "
                 "popülasyon birkaç nesilde etkisiz hale gelir. "
                 "Konağın ΔdapA bağımlılığı bakterinin bağırsak "
                 "dışında yaşamasını engeller ama plazmiti tutmaz - "
                 "ikisi farklı şeydir ve karıştırılması sık bir "
                 "hatadır.\n\n"
                 "**Denetim burada CORE-02 hatası verir ve bu "
                 "beklenen:** markersiz bir plazmitte transformasyondan "
                 "sonra doğru koloniyi seçmenin yolu bile yoktur. "
                 "Basamak, merdivenin diğer üçünün ne işe yaradığını "
                 "göstermek için duruyor."),
    ]


#: Konak kimliğinin okunur karşılığı (strateji başlıklarında kullanılır).
_KONAK_ADI = {"LLNZ9000": "L. lactis NZ9000",
              "BSUB168": "B. subtilis 168",
              "BL21DE3": "BL21(DE3)"}


#: Gram-pozitif konakta kullanılabilecek idame markeri. L. lactis'te
#: thyA tamamlaması gerçek bir yoldur (Steidler 2003, timin bağımlı
#: içerme); B. subtilis için katalogda thyA eksik bir konak yok, o
#: yüzden orada antibiyotik markeri kullanılır ve bu fark gizlenmez.
_GP_MARKER = {
    "LLNZ9000": ("MRK_THYA_LL", "thyA içermesi",
                 "Konağın thyA'sı silinir ve gen plazmite konur: timin "
                 "bulunmayan ortamda (bağırsak dışı) plazmitsiz hücre "
                 "ölür. Antibiyotik gerekmez - Steidler 2003'ün "
                 "L. lactis'te kurduğu içerme yolu budur."),
    "BSUB168": ("REF_ERM_C_735", "eritromisin markeri",
                "pE194'ten çıkarılmış erm(C): Gram-pozitif bakterilerde "
                "çalışan gerçek bir seçilim markeri. Burada thyA "
                "içermesi kullanılamıyor çünkü katalogda thyA'sı "
                "silinmiş bir B. subtilis suşu yok - antibiyotiksiz "
                "idame bu konakta henüz kurulamıyor ve bu bir eksiktir, "
                "tercih değil."),
}


def _gram_poz_merdiven(kaset: list[str], ekoli_kaset: list[str],
                       konak: str, ne: str) -> list[Strategy]:
    """Gram-pozitif konaklı hastalıklar için üç basamaklı merdiven.

    A basamağı katalogdaki açığı bilerek gösterir: Gram-pozitif
    replikon yok, o yüzden denetim CORE-01 hatası verir. Bu hata
    gizlenmedi; iskelet "Referanstan kes-yapıştır" sayfasındaki
    pNZ8048 ve pUB110 tariflerinden gelir ve açık orada kapanır.
    """
    mrk, mrk_ad, mrk_aciklama = _GP_MARKER[konak]
    return [
        Strategy("A", f"Gram-pozitif konak, {mrk_ad} "
                 f"({_KONAK_ADI.get(konak, konak)})",
                 kaset + [mrk], konak,
                 f"{ne} {mrk_aciklama}\n\n"
                 "**Denetim burada CORE-01 hatası verir ve bu "
                 "beklenen:** katalogda bu bakteride çoğalan bir "
                 "replikon yok. Eksik olan iskelettir, kaset değil; "
                 "iskeleti Referanstan kes-yapıştır sayfasındaki "
                 "pNZ8048 tarifi verir."),
        Strategy("B", "E. coli'de ön üretim (mekik yaklaşımı)",
                 ekoli_kaset + ["ORI_PUC", "MRK_KAN"], "BL21DE3",
                 "Kaseti önce E. coli'de kurup dizileyip doğrularsınız; "
                 "gıda sınıfı konağa ancak doğru olduğundan emin olunca "
                 "geçersiniz. Pratikte herkes böyle çalışır. Promotör "
                 "ve sinyal peptidi E. coli'ye uygun olanla "
                 "değiştirildi: Gram-pozitif bir promotör E. coli'de "
                 "okunmayabilir, Usp45 sinyali de E. coli'nin salgı "
                 "yoluna ait değildir.\n\n"
                 "Konak **BL21(DE3)**, DH5-alfa değil: kaset T7 "
                 "promotörü taşıyor ve T7 polimerazı yalnızca DE3 "
                 "lizojenleri üretir. Bu tasarım ilk yazıldığında konak "
                 "DH5-alfa'ydı ve denetim PRO-03 ile yakaladı - en sık "
                 "yapılan ifade hatası tam olarak budur."),
        Strategy("C", "E. coli ön üretim + ccdAB idamesi",
                 ekoli_kaset + ["ORI_PBR322"] + _CCDAB, "BL21DE3",
                 "Aynı ön üretim, antibiyotik yerine toksin-antitoksin "
                 "idamesiyle. Düşük kopyalı pBR322 orijini seçildi: "
                 "yüksek kopya, büyük bir terapötik geni taşırken "
                 "metabolik yükü gereksiz yükseltir ve plazmit "
                 "kararsızlaşır."),
    ]


DISEASES += [
    # ------------------------------------------------------ Homosistinüri
    Disease(
        key="homosistinuri",
        name="Homosistinüri (sistationin beta-sentaz eksikliği)",
        short="HCU",
        logic="Bağırsakta metiyonini tüketen bir probiyotik.",
        mechanism=(
            "Sistationin beta-sentaz (CBS) eksikliğinde homosistein "
            "yıkılamaz ve birikir: damar tıkanıklığı, göz merceğinin "
            "yerinden kayması ve iskelet bozuklukları yapar. Homosistein "
            "metiyoninden gelir, bu yüzden standart tedavi metiyoninden "
            "kısıtlı diyettir - ve uyumu çok zordur, çünkü metiyonin "
            "neredeyse bütün proteinlerde vardır.\n\n"
            "Makaledeki fikir PKU'dakinin aynısıdır ama substrat "
            "farklıdır: bağırsakta metiyonini tüketen bir bakteri, "
            "diyetin yapacağı işi diyet olmadan yapar."),
        papers=[
            "Perreault M ve ark. The live biotherapeutic SYNB1353 "
            "decreases plasma methionine via directed degradation in "
            "animal models and healthy volunteers. Cell Host & Microbe "
            "32:382-395.e10 (2024). doi:10.1016/j.chom.2024.01.005 "
            "· PMID 38309259",
        ],
        pathway=_gv(
            'a [label="Diyetten gelen\\nmetiyonin"];',
            'b [label="Bakteri metiyonini\\niçine alır"];',
            'c [label="MGL (mdeA)\\nmetiyonin → metantiyol\\n+ '
            'alfa-ketobütirat + NH3"];',
            'd [label="Homosistein havuzu\\nküçülür"];',
            'e [label="Plazma metiyonini\\ndüşer"];',
            "a -> b -> c -> d; c -> e;"),
        pathway_caption=(
            "PfnrS altındadır, yani enzim yalnızca bağırsağın oksijensiz "
            "ortamında üretilir. Metantiyol kokulu bir bileşiktir; "
            "makalede yan etki olarak izlenir."),
        cassette=HCU_KASET,
        biomarker=("Plazma metiyonin (çalışmanın birincil ölçütü); dolaylı "
                   "olarak total homosistein"),
        gap=("Enzim hücre içindedir, yani metiyoninin bakteriye girmesi "
             "gerekir. PKU tasarımında bu işi PheP taşıyıcısı yapar; "
             "katalogda metiyonin taşıyıcısı (metNIQ / metP) yok, bu "
             "yüzden buradaki tasarım taşıyıcısızdır ve verimi düşük "
             "kalabilir. Denetim THR-02'nin metiyonin karşılığını "
             "bilmediği için bunu söylemez - bu, aracın bir sınırıdır."),
        strategies=_ecoli_merdiven(
            HCU_KASET, "ORI_P15A",
            "Tek enzimli en basit tasarım: taşıyıcı yok, yalnızca MGL."),
        comparison=[
            {"Özellik": "Konak", "Makale": "E. coli Nissle 1917 türevi",
             "Bizim tasarım": "E. coli Nissle 1917 ΔdapA",
             "Etiket": "Literatür"},
            {"Özellik": "Enzim", "Makale": "Metiyonin gama-liyaz (MGL)",
             "Bizim tasarım": "mdeA / MGL (P. putida, UniProt P13254)",
             "Etiket": "Literatür"},
            {"Özellik": "Metiyonin alımı",
             "Makale": "Alım yolu da güçlendirilir",
             "Bizim tasarım": "Taşıyıcı yok - katalogda metNIQ bulunmuyor",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Promotör", "Makale": "Bağırsakta açılan promotör",
             "Bizim tasarım": "PfnrS (anaerobik)", "Etiket": "Öneri"},
            {"Özellik": "İdame", "Makale": "Kromozoma entegre",
             "Bizim tasarım": "Plazmit + ccdAB (A stratejisi)",
             "Etiket": "Öneri"},
        ],
    ),
    # ------------------------------------------------------------- Kolera
    Disease(
        key="kolera",
        name="Kolera (Vibrio cholerae enfeksiyonu)",
        short="Kolera",
        logic="Patojenin haberleşmesini taklit ederek virülansını kapatır.",
        mechanism=(
            "V. cholerae'nın virülans kararı quorum sensing'e bağlıdır ve "
            "mantık TERSTİR: hücre yoğunluğu **düşükken** ToxR/TcpP yolu "
            "açılır ve kolera toksini üretilir; yoğunluk **yüksekken** "
            "CAI-1 sinyali birikir ve virülans kapanır. Patojen, "
            "bağırsağa yeni girdiğinde (seyrekken) saldırır.\n\n"
            "Makaledeki fikir bu mantığı kullanır: bağırsakta bol CAI-1 "
            "üreten zararsız bir E. coli, patojeni ortamı kalabalık "
            "sanmaya zorlar. Toksin üretimi hiç açılmaz. Patojen "
            "öldürülmediği için seçilim baskısı da yaratılmaz."),
        papers=[
            "Duan F, March JC. Engineered bacterial communication prevents "
            "Vibrio cholerae virulence in an infant mouse model. "
            "Proceedings of the National Academy of Sciences 107:11260-"
            "11264 (2010). doi:10.1073/pnas.1001294107 · PMID 20534565",
        ],
        pathway=_gv(
            'a [label="Bizim bakteri\\ncqsA ifade eder"];',
            'b [label="CAI-1 sinyali\\nbağırsakta birikir"];',
            'c [label="V. cholerae CqsS\\nsinyali algılar"];',
            'd [label="HapR aktif\\n→ ToxT baskılanır"];',
            'e [label="Kolera toksini\\nüretilmez"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "Sabit (J23105) promotör seçildi: sinyal patojen gelmeden "
            "ÖNCE ortamda birikmiş olmalı, yoksa koruma gecikir."),
        cassette=KOLERA_KASET,
        biomarker=("Yavru fare modelinde sağ kalım ve bağırsakta kolera "
                   "toksini; in vitro CqsS raportör aktivitesi"),
        gap=("CAI-1 üretimi konağın kendi yağ asidi havuzuna bağlıdır: "
             "cqsA tek başına yeterli sinyal vermeyebilir. Ayrıca "
             "virülansı baskılamak öldürmek değildir - patojen bağırsakta "
             "kalır ve koruma yalnızca sinyal sürdüğü kadar sürer."),
        strategies=_ecoli_merdiven(
            KOLERA_KASET, "ORI_P15A",
            "Sinyal molekülü üreten tek enzimli tasarım: ürün protein "
            "değil, küçük bir molekül."),
        comparison=[
            {"Özellik": "Konak", "Makale": "E. coli Nissle 1917",
             "Bizim tasarım": "E. coli Nissle 1917 ΔdapA",
             "Etiket": "Literatür"},
            {"Özellik": "Gen", "Makale": "cqsA (V. cholerae)",
             "Bizim tasarım": "cqsA (UniProt Q9KM65)",
             "Etiket": "Literatür"},
            {"Özellik": "Etki mekanizması",
             "Makale": "CAI-1 taklidi, virülans baskılanır",
             "Bizim tasarım": "Aynı", "Etiket": "Literatür"},
            {"Özellik": "Promotör", "Makale": "Sabit ifade",
             "Bizim tasarım": "J23105 (sabit, orta güç)",
             "Etiket": "Öneri"},
            {"Özellik": "Öldürücü etki", "Makale": "Yok (bilinçli)",
             "Bizim tasarım": "Yok", "Etiket": "Literatür"},
        ],
    ),
    # --------------------------------------------- Antibiyotiğe bağlı disbiyoz
    Disease(
        key="cdiff",
        name="Antibiyotiğe bağlı disbiyoz ve C. difficile riski",
        short="C. difficile",
        logic="Bağırsağa sızan antibiyotiği orada parçalayıp florayı korur.",
        mechanism=(
            "Damardan verilen seftriakson gibi beta-laktamlar safrayla "
            "bağırsağa geçer. Orada hedefi olmayan bir antibiyotik "
            "mikrobiyotayı kırar; boşalan yeri Clostridioides difficile "
            "doldurur ve ağır bir kolit yapar.\n\n"
            "SYN-004 (ribaxamase) ağızdan verilen bir beta-laktamazdır: "
            "bağırsak lümeninde antibiyotiği parçalar, kana geçmediği "
            "için damardaki tedaviye dokunmaz. Faz 2b çalışmasında "
            "C. difficile enfeksiyonunu azalttı."),
        papers=[
            "Kokai-Kun JF ve ark. Use of ribaxamase (SYN-004), a "
            "beta-lactamase, to prevent Clostridium difficile infection "
            "in beta-lactam-treated patients: a double-blind, phase 2b, "
            "randomised placebo-controlled trial. The Lancet Infectious "
            "Diseases 19:487-496 (2019). "
            "doi:10.1016/S1473-3099(18)30731-X · PMID 30885591",
            "Kaleko M ve ark. Development of SYN-004, an oral "
            "beta-lactamase treatment to protect the gut microbiome from "
            "antibiotic-mediated damage and prevent Clostridium difficile "
            "infection. Anaerobe 41:58-67 (2016). "
            "doi:10.1016/j.anaerobe.2016.05.015 · PMID 27262694",
        ],
        pathway=_gv(
            'a [label="Damardan\\nseftriakson"];',
            'b [label="Safrayla bağırsağa\\nsızar"];',
            'c [label="P1A beta-laktamaz\\nhalkayı açar"];',
            'd [label="Mikrobiyota\\nkorunur"];',
            'e [label="C. difficile\\nçoğalamaz"];',
            'f [label="Kandaki tedavi\\netkilenmez"];',
            "a -> b -> c -> d -> e; a -> f;"),
        pathway_caption=(
            "Enzimin bağırsakta kalması şarttır: kana geçerse asıl "
            "tedaviyi de bozar. Salgılanan enzim lümende çalışır."),
        cassette=CDIFF_KASET,
        biomarker="Dışkıda seftriakson düzeyi ve mikrobiyota çeşitliliği",
        gap=("**Bu tasarım bir `Öneri`dir, literatürden alınmamıştır.** "
             "Makalelerde enzimin kendisi ağızdan veriliyor; canlı "
             "bakteri kullanılmıyor. Sebebi ciddi: terapötik "
             "beta-laktamaz, direnç geniyle AYNI enzimdir. Yatay gen "
             "aktarımıyla bir patojene geçerse önlemeye çalıştığı sorunu "
             "yaratır. Canlı bakteriyle yapılacaksa gen plazmitte değil "
             "kromozomda, mobil olmayan bir bağlamda durmalıdır - bu "
             "sayfadaki plazmit tasarımları o yüzden yalnızca birer "
             "düşünce denemesidir."),
        strategies=_gram_poz_merdiven(
            CDIFF_KASET,
            ["PRO_T7", "RBS_PET", "SIG_PELB", "CDS_BLAC_BL", "TER_T7"],
            "BSUB168",
            "Bacillus subtilis salgılama konağı: AmyE sinyali enzimi "
            "besiyerine çıkarır."),
        comparison=[
            {"Özellik": "Verilme şekli",
             "Makale": "Ağızdan enzim (canlı bakteri YOK)",
             "Bizim tasarım": "Enzimi bağırsakta üreten bakteri",
             "Etiket": "Öneri"},
            {"Özellik": "Enzim", "Makale": "P1A beta-laktamaz (ribaxamase)",
             "Bizim tasarım": "penP / P1A (B. licheniformis, P00808)",
             "Etiket": "Literatür"},
            {"Özellik": "Etki yeri", "Makale": "Bağırsak lümeni",
             "Bizim tasarım": "Aynı - salgılanan enzim",
             "Etiket": "Literatür"},
            {"Özellik": "Direnç geni riski",
             "Makale": "Yok (gen taşıyan bakteri verilmiyor)",
             "Bizim tasarım": "VAR - tasarımın temel zaafı",
             "Etiket": "Öneri (riskli)"},
        ],
    ),
    # ---------------------------------------------------- Laktoz intoleransı
    Disease(
        key="laktoz",
        name="Laktoz intoleransı (laktaz eksikliği)",
        short="Laktoz",
        logic="Bağırsakta laktozu parçalayan gıda sınıfı bir bakteri.",
        mechanism=(
            "Erişkinlerin çoğunda ince bağırsak laktazı (LCT) azalır. "
            "Sindirilmeyen laktoz kolona geçer, orada bakteriler onu "
            "fermente eder: gaz, şişkinlik, kramp ve ishal. Hastalık "
            "değil, bir enzim yetmezliğidir - ama yaşam kalitesini "
            "belirgin düşürür.\n\n"
            "Fikir: laktozu parçalayan beta-galaktosidazı bağırsakta "
            "üreten gıda sınıfı bir laktik asit bakterisi. Hap halinde "
            "verilen enzime göre üstünlüğü, enzimin yerinde ve "
            "sürekli üretilmesidir."),
        papers=[
            "Li J ve ark. Lactococcus lactis expressing food-grade "
            "beta-galactosidase alleviates lactose intolerance symptoms "
            "in post-weaning Balb/c mice. Applied Microbiology and "
            "Biotechnology 96:1499-1506 (2012). "
            "doi:10.1007/s00253-012-3977-4 · PMID 22395907",
        ],
        pathway=_gv(
            'a [label="Sütteki laktoz"];',
            'b [label="İnce bağırsakta\\nlaktaz yetersiz"];',
            'c [label="Bakteri laktozu\\niçine alır"];',
            'd [label="Beta-galaktosidaz\\nlaktoz → glukoz + galaktoz"];',
            'e [label="Kolona geçen\\nlaktoz azalır"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "Enzim hücre içindedir: laktozun bakteriye girmesi gerekir. "
            "Gıda sınıfı konakta kendi laktoz taşıyıcısı (lacS) vardır."),
        cassette=LAKTOZ_KASET,
        biomarker="Laktoz solunum hidrojen testi; dışkı pH'ı ve semptom skoru",
        gap=("Beta-galaktosidaz 1007 amino asitlik büyük bir proteindir: "
             "salgılanması kolay değildir, bu yüzden tasarım hücre içi "
             "ifadeye dayanır. Bu da etkiyi bakterinin laktoz alım hızına "
             "bağlar - lümendeki laktoza doğrudan dokunulmaz. Hap "
             "enzimine karşı asıl üstünlük ancak salgılanan bir sürümle "
             "kurulabilir."),
        strategies=_gram_poz_merdiven(
            LAKTOZ_KASET,
            ["PRO_LAC", "RBS_B0034", "CDS_LACZ_LD", "TER_RRNB"],
            "LLNZ9000",
            "Gıda sınıfı Lactococcus lactis: düzenleyici açıdan en "
            "uygun konak, çünkü yoğurt ve peynirde zaten kullanılır."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Lactococcus lactis",
             "Bizim tasarım": "Lactococcus lactis NZ9000",
             "Etiket": "Literatür"},
            {"Özellik": "Enzim kaynağı",
             "Makale": "Gıda sınıfı beta-galaktosidaz",
             "Bizim tasarım": "lacZ (L. delbrueckii, P0C1Y0) - LAB "
                              "kaynaklı",
             "Etiket": "Öneri"},
            {"Özellik": "Konum", "Makale": "Hücre içi",
             "Bizim tasarım": "Hücre içi", "Etiket": "Literatür"},
            {"Özellik": "Promotör", "Makale": "Gıda sınıfı sabit promotör",
             "Bizim tasarım": "L. reuteri EF-Tu promotörü",
             "Etiket": "Öneri"},
            {"Özellik": "İskelet", "Makale": "L. lactis plazmiti",
             "Bizim tasarım": "Katalogda Gram-pozitif replikon yok; "
                              "pNZ8048 tarifinden gelir",
             "Etiket": "Öneri (eksik)"},
        ],
    ),
    # --------------------------------------------------------------- MRSA
    Disease(
        key="mrsa",
        name="MRSA (metisiline dirençli S. aureus) kolonizasyonu",
        short="MRSA",
        logic="S. aureus'un hücre duvarını kesen lizostafin salgılar.",
        mechanism=(
            "S. aureus'un deride ve burunda taşınması enfeksiyonun "
            "kaynağıdır; metisiline dirençli suşlarda (MRSA) tedavi "
            "seçeneği azdır. Antibiyotikle dekolonizasyon hem direnci "
            "seçer hem de cilt florasının tamamını vurur.\n\n"
            "Lizostafin, S. aureus hücre duvarındaki pentaglisin "
            "köprülerini kesen bir glisil-glisin endopeptidazdır. Bu "
            "köprü başka bakterilerde yoktur: hedef son derece dardır. "
            "Fikir, enzimi deride yerinde üreten bir probiyotiktir."),
        papers=[
            "Guan C ve ark. Engineering a \"detect and destroy\" skin "
            "probiotic to combat methicillin-resistant Staphylococcus "
            "aureus. PLoS One 17:e0276795 (2022). "
            "doi:10.1371/journal.pone.0276795 · PMID 36520793",
            "Turner MS ve ark. Antimicrobial activity of lysostaphin and "
            "a Listeria monocytogenes bacteriophage endolysin produced "
            "and secreted by lactic acid bacteria. Systematic and Applied "
            "Microbiology 30:58-67 (2007). "
            "doi:10.1016/j.syapm.2006.01.013 · PMID 16490333",
        ],
        pathway=_gv(
            'a [label="Probiyotik deride\\nlizostafin salgılar"];',
            'b [label="Enzim S. aureus\\nhücre duvarına bağlanır"];',
            'c [label="Pentaglisin köprüleri\\nkesilir"];',
            'd [label="Hücre duvarı\\ndayanamaz, liziz"];',
            'e [label="Cilt florasının geri\\nkalanı etkilenmez"];',
            "a -> b -> c -> d; a -> e;"),
        pathway_caption=(
            "Usp45 sinyali enzimi hücre dışına çıkarır: lizostafin "
            "hedefine ancak salgılanırsa ulaşır, hücre içinde işe yaramaz."),
        cassette=MRSA_KASET,
        biomarker=("Deri sürüntüsünde S. aureus kolonisi sayımı; in vitro "
                   "liziz halkası (zone of clearance)"),
        gap=("Makaledeki devrenin **'algıla' yarısı burada yok.** Özgün "
             "tasarım S. aureus'un AIP sinyalini okuyan AgrC/AgrA "
             "çiftini kullanır ve enzimi yalnızca patojen varken üretir; "
             "katalogda bu çift bulunmuyor, bu yüzden buradaki tasarım "
             "SÜREKLİ üretir. Sürekli üretim, hedefi dar olsa bile "
             "seçilim baskısı yaratır ve lizostafine direnç gelişimini "
             "hızlandırır. Algılamalı bir sürüm için P. aeruginosa "
             "sayfasındaki LasR devresine bakın: aynı desenin "
             "katalogdaki çalışan örneği odur."),
        strategies=_gram_poz_merdiven(
            MRSA_KASET,
            ["PRO_T7", "RBS_PET", "SIG_PELB", "CDS_LSS", "TER_T7"],
            "LLNZ9000",
            "Gıda sınıfı Lactococcus lactis, Usp45 ile salgılama: "
            "lizostafinin laktik asit bakterilerinde salgılandığı "
            "Turner 2007'de gösterilmiştir."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Mühendislik ürünü cilt "
                                           "probiyotiği",
             "Bizim tasarım": "Lactococcus lactis NZ9000",
             "Etiket": "Öneri"},
            {"Özellik": "Öldürücü", "Makale": "Lizostafin",
             "Bizim tasarım": "Lizostafin, olgun bölge (P10547 248-493)",
             "Etiket": "Literatür"},
            {"Özellik": "Algılayıcı", "Makale": "S. aureus AIP sinyali "
                                                "(agr sistemi)",
             "Bizim tasarım": "YOK - sürekli üretim",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Salgı sinyali", "Makale": "Konağın kendi sinyali",
             "Bizim tasarım": "Usp45 (L. lactis)", "Etiket": "Öneri"},
            {"Özellik": "Hedef darlığı", "Makale": "Yalnızca S. aureus",
             "Bizim tasarım": "Aynı", "Etiket": "Literatür"},
        ],
    ),
    # ------------------------------------------------------------ Diyabet
    Disease(
        key="diyabet",
        name="Tip 2 diyabet (glukoz toleransı bozukluğu)",
        short="Diyabet",
        logic="Bağırsakta GLP-1 salgılayarak insülin salımını destekler.",
        mechanism=(
            "GLP-1 bir inkretindir: yemekten sonra bağırsaktan salgılanır "
            "ve pankreastan insülin salımını artırır. Tip 2 diyabette "
            "inkretin yanıtı zayıflamıştır; GLP-1 benzeri ilaçlar bu "
            "yüzden işe yarar ama enjeksiyonla verilir.\n\n"
            "Duan ve ark. GLP-1(1-37)'yi bağırsakta salgılayan bir "
            "commensal bakteri kullandı ve bağırsak epitel hücrelerinin "
            "glukoza duyarlı şekilde insülin salgılamaya başladığını "
            "gösterdi: hedef pankreas değil, bağırsak epitelinin "
            "kendisidir."),
        papers=[
            "Duan FF ve ark. Engineered commensal bacteria reprogram "
            "intestinal cells into glucose-responsive insulin-secreting "
            "cells for the treatment of diabetes. Diabetes 64:1794-1803 "
            "(2015). doi:10.2337/db14-0635 · PMID 25626737",
            "Duan F ve ark. Secretion of insulinotropic proteins by "
            "commensal bacteria: rewiring the gut to treat diabetes. "
            "Applied and Environmental Microbiology 74:7437-7438 (2008). "
            "doi:10.1128/AEM.01019-08 · PMID 18836005",
        ],
        pathway=_gv(
            'a [label="Bakteri bağırsakta\\nGLP-1(1-37) salgılar"];',
            'b [label="Bağırsak epitel\\nhücresi uyarılır"];',
            'c [label="Hücre glukoza duyarlı\\ninsülin salgılar"];',
            'd [label="Kan glukozu\\ndüşer"];',
            "a -> b -> c -> d;"),
        pathway_caption=(
            "GLP-1(1-37) olgun GLP-1'in öncüsüdür; epitel hücresinin "
            "yeniden programlanmasını sağlayan sürüm budur."),
        cassette=GLP1_KASET,
        biomarker="Oral glukoz tolerans testi; epitelde insülin immünboyama",
        gap=("İki fark var. Birincisi konak: makalede *Lactobacillus "
             "gasseri* kullanılır, katalogda o suş yok, tasarım "
             "L. lactis üzerine kuruldu (`Öneri`). İkincisi ve daha "
             "önemlisi: GLP-1(1-37) 37 amino asitlik bir peptittir. "
             "Peptit terapötiklerinde darboğaz ifade değil "
             "**kararlılıktır** - bakteri içinde proteazlara, "
             "salgılandıktan sonra DPP-4'e açıktır. Uygulamanın protein "
             "hesapları bu kısalıkta bir peptidin yarı ömrü hakkında bir "
             "şey söylemez; bu, aracın bir sınırıdır."),
        strategies=_gram_poz_merdiven(
            GLP1_KASET,
            ["PRO_T7", "RBS_PET", "SIG_DSBA", "CDS_GLP1", "TER_T7"],
            "LLNZ9000",
            "Gıda sınıfı Lactococcus lactis, Usp45 ile salgılama."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Lactobacillus gasseri "
                                           "ATCC 33323",
             "Bizim tasarım": "Lactococcus lactis NZ9000 - katalogda "
                              "L. gasseri yok",
             "Etiket": "Öneri"},
            {"Özellik": "Terapötik", "Makale": "GLP-1(1-37)",
             "Bizim tasarım": "GLP-1(1-37), P01275'in 92-128 bölgesi",
             "Etiket": "Literatür"},
            {"Özellik": "Hedef hücre", "Makale": "Bağırsak epiteli",
             "Bizim tasarım": "Aynı", "Etiket": "Literatür"},
            {"Özellik": "Salgı sinyali", "Makale": "Konağa özgü",
             "Bizim tasarım": "Usp45 (L. lactis)", "Etiket": "Öneri"},
            {"Özellik": "Peptit kararlılığı", "Makale": "Çalışmada ölçülür",
             "Bizim tasarım": "Hesaplanamaz - aracın sınırı",
             "Etiket": "Öneri (eksik)"},
        ],
    ),
    # --------------------------------------------------------------- AI-2
    Disease(
        key="ai2",
        name="Patojen haberleşmesini bozma (AI-2 quorum söndürme)",
        short="AI-2",
        logic="Ortak quorum sinyali AI-2'yi ortamdan çekerek patojeni susturur.",
        mechanism=(
            "Quorum sensing sinyallerinin çoğu türe özgüdür, AI-2 ise "
            "**türler arası ortak bir dildir**: Salmonella, V. cholerae, "
            "EHEC ve daha pek çok bakteri virülans ve biyofilm "
            "kararlarını AI-2 yoğunluğuna göre verir.\n\n"
            "LsrK kinazı AI-2'yi fosforile eder; fosforile AI-2 hücreye "
            "çekilir ve ortamdan silinir. Fazla LsrK üreten bir bakteri, "
            "komşu patojenin okuyacağı sinyali tüketir. Xavier ve Bassler "
            "bu müdahalenin işe yaradığını gösterdi; Roy ve ark. enzimi "
            "doğrudan bir söndürme aracı olarak kullandı."),
        papers=[
            "Xavier KB, Bassler BL. Interference with AI-2-mediated "
            "bacterial cell-cell communication. Nature 437:750-753 "
            "(2005). doi:10.1038/nature03960 · PMID 16193054",
            "Roy V ve ark. Cross species quorum quenching using a native "
            "AI-2 processing enzyme. ACS Chemical Biology 5:223-232 "
            "(2010). doi:10.1021/cb9002738 · PMID 20025244",
        ],
        pathway=_gv(
            'a [label="Patojen AI-2\\nsalgılar"];',
            'b [label="Bizim bakteri\\nfazla LsrK üretir"];',
            'c [label="AI-2 fosforile edilip\\nhücreye çekilir"];',
            'd [label="Ortamda AI-2\\nkalmaz"];',
            'e [label="Patojen kalabalığı\\nalgılayamaz"];',
            'f [label="Virülans ve biyofilm\\naçılmaz"];',
            "a -> c; b -> c -> d -> e -> f;"),
        pathway_caption=(
            "Öldürme yok: yalnızca sinyal çekilir. Direnç gelişimi için "
            "seçilim baskısı yaratmaması bu yaklaşımın en güçlü yanıdır."),
        cassette=AI2_KASET,
        biomarker="V. harveyi BB170 biyoraportörüyle ortam AI-2 aktivitesi",
        gap=("İki sınır var. LsrK fosforilasyon için ATP harcar ve ürün "
             "hücre içinde birikir: doymuş bir sistemde söndürme durur, "
             "yani etki sınırsız değildir. İkincisi daha önemli: AI-2 "
             "ortak bir dildir, onu tüketmek **yararlı komşuların** "
             "haberleşmesini de bozar. Etki patojene özgü değildir ve "
             "mikrobiyota üzerindeki sonucu bu sayfada "
             "hesaplanamaz."),
        strategies=_ecoli_merdiven(
            AI2_KASET, "ORI_PBR322",
            "Sinyal yıkan tek enzimli tasarım: ürün ne protein ne "
            "molekül, ortamdan çekilen bir sinyaldir."),
        comparison=[
            {"Özellik": "Konak", "Makale": "E. coli (laboratuvar suşları)",
             "Bizim tasarım": "E. coli Nissle 1917 ΔdapA - terapötik "
                              "bağlam",
             "Etiket": "Öneri"},
            {"Özellik": "Enzim", "Makale": "LsrK (AI-2 kinaz)",
             "Bizim tasarım": "lsrK (E. coli, UniProt P77432)",
             "Etiket": "Literatür"},
            {"Özellik": "Mekanizma", "Makale": "AI-2'yi fosforile edip çek",
             "Bizim tasarım": "Aynı", "Etiket": "Literatür"},
            {"Özellik": "Hedef özgüllüğü", "Makale": "Türler arası (geniş)",
             "Bizim tasarım": "Geniş - yararlı bakterileri de etkiler",
             "Etiket": "Öneri (eksik)"},
            {"Özellik": "Promotör", "Makale": "İndüklenebilir ifade",
             "Bizim tasarım": "J23100 (sabit, güçlü)", "Etiket": "Öneri"},
        ],
    ),
    Disease(
        key="oksidatif",
        name="Oksidatif bağırsak hasarı (kolit)",
        short="Oksidatif hasar",
        logic="Bağırsakta serbest radikali söndüren bir laktik asit "
              "bakterisi.",
        mechanism=(
            "İltihaplı bağırsakta nötrofil ve makrofajlar bol miktarda "
            "superoksit (O2⁻) üretir. Superoksit kendi başına da "
            "zarar verir, ama asıl sorun zincirin devamıdır: hücrenin "
            "kendi enzimleri onu hidrojen perokside (H2O2) çevirir ve "
            "H2O2 demirle buluşunca hidroksil radikali doğar - "
            "hücrenin hiçbir enzimiyle söndüremediği, DNA'yı ve zarı "
            "doğrudan kesen bir tür.\n\n"
            "Fikir: radikali üretildiği yerde, bağırsak lümeninde "
            "söndüren bir bakteri. Ağızdan verilen antioksidanın "
            "ulaşamadığı yer tam olarak burasıdır."),
        papers=[
            "LeBlanc JG ve ark. Improvement of an experimental colitis "
            "in rats by lactic acid bacteria producing superoxide "
            "dismutase. Inflammatory Bowel Diseases 12:1069-1077 "
            "(2006). doi:10.1097/01.mib.0000235101.09231.9e "
            "· PMID 17075345",
        ],
        pathway=_gv(
            'a [label="İltihap bölgesinde\\nnötrofil ve makrofaj"];',
            'b [label="Superoksit (O2⁻)\\nsalınır"];',
            'c [label="SOD (sodA)\\nO2⁻ → H2O2"];',
            'd [label="Katalaz (katE)\\nH2O2 → H2O + O2"];',
            'e [label="Hidroksil radikali\\nhiç oluşmaz"];',
            "a -> b -> c -> d -> e;"),
        pathway_caption=(
            "İki enzimin SIRASI tasarımın kendisidir. Yalnızca SOD "
            "koymak superoksiti H2O2'de bırakır ve H2O2 de zararlıdır; "
            "yalnızca katalaz koymak superoksite hiç dokunmaz. Zincirin "
            "iki basamağı birden kapatılmazsa sorun yer değiştirir, "
            "çözülmez."),
        cassette=SOD_KASET,
        biomarker=("Dışkıda kalprotektin; doku malondialdehit (MDA) ve "
                   "miyeloperoksidaz (MPO) etkinliği"),
        gap=("Katalaz (katE) 753 amino asitlik büyük bir proteindir ve "
             "hem grubu gerektirir; Lactococcus lactis hem sentezini "
             "sınırlı yapar, bu yüzden katalaz etkinliği konakta düşük "
             "kalabilir. Makale yalnızca SOD üreten bir suş kullandı - "
             "iki enzimi birlikte koymak bizim önerimizdir ve "
             "deneysel olarak sınanmamıştır. Ayrıca bu iki enzim hücre "
             "içinde çalışır; lümendeki radikale ulaşmaları için "
             "bakterinin orada ve canlı olması gerekir."),
        strategies=_gram_poz_merdiven(
            SOD_KASET,
            ["PRO_LAC", "RBS_B0034", "CDS_SODA", "CDS_KATE", "TER_B0015"],
            "LLNZ9000",
            "Laktik asit bakterisi: iltihaplı bağırsakta hayatta kalır "
            "ve gıda sınıfı olduğu için düzenleyici açıdan uygundur."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Lactococcus lactis",
             "Bizim tasarım": "Lactococcus lactis NZ9000",
             "Etiket": "Literatür"},
            {"Özellik": "Enzim", "Makale": "Yalnızca SOD",
             "Bizim tasarım": "SOD (sodA, P0A4J2) + katalaz "
                              "(katE, P21179)",
             "Etiket": "Öneri"},
            {"Özellik": "Gerekçe",
             "Makale": "Superoksiti söndürmek",
             "Bizim tasarım": "H2O2'yi de söndürmek: SOD tek başına "
                              "radikali H2O2'de bırakır",
             "Etiket": "Öneri"},
            {"Özellik": "Konum", "Makale": "Hücre içi",
             "Bizim tasarım": "Hücre içi", "Etiket": "Literatür"},
            {"Özellik": "Model", "Makale": "Sıçanda deneysel kolit",
             "Bizim tasarım": "—", "Etiket": "Literatür"},
        ],
    ),

]

# IBD'ye ikinci bir terapötik protein: elafin. Aynı hastalık, bambaşka
# bir mantık - IL-10 bağışıklığı yatıştırır, elafin hasarı yapan enzimi
# durdurur. İkisini yan yana görmek "hangi terapötik?" sorusunu
# sorulabilir hale getirir.
_ibd = by_key("ibd")
if _ibd:
    _ibd.strategies += [
        Strategy("D", "IL-10 yerine elafin (proteaz inhibitörü yolu)",
                 ELAFIN_KASET + ["MRK_THYA_LL"], "LLNZ9000",
                 "Aynı konak, aynı salgı sinyali, farklı terapötik: "
                 "elafin nötrofil elastazını ve proteinaz-3'ü bağlar. "
                 "IL-10 bağışıklık yanıtını yatıştırır, elafin doku "
                 "yıkımını yapan enzimi durdurur - biri sinyali, öbürü "
                 "hasarı hedefler. Elafin 57 amino asitlik küçük ve "
                 "kararlı bir inhibitördür, bu yüzden salgılanması "
                 "IL-10'dan kolaydır.\n\n"
                 "Denetim burada da CORE-01 verir: Gram-pozitif "
                 "replikon açığı bu stratejide de duruyor."),
        Strategy("E", "Elafin, E. coli'de ön üretim",
                 ["PRO_T7", "RBS_PET", "SIG_DSBA", "CDS_ELAFIN", "TER_T7",
                  "ORI_PUC", "MRK_KAN"], "BL21DE3",
                 "Elafinin kendisi disülfür bağı taşıyan küçük bir "
                 "proteindir: DsbA sinyaliyle periplazmaya gönderilir, "
                 "çünkü disülfür bağları orada kurulur. Sitoplazmada "
                 "ifade edilse yanlış katlanır. Bu, IL-10 kasetinin "
                 "karşılaşmadığı bir sorundur ve sinyal peptidi "
                 "seçiminin neden proteine göre yapıldığını gösterir."),
    ]

# Amonyak ve oksalat tasarımlarına merdivenin eksik basamakları.
for _k, _kaset, _ori in (("amonyak", AMONYAK_KASET, "ORI_P15A"),
                         ("oksalat", OKSALAT_KASET, "ORI_P15A")):
    _d = by_key(_k)
    if _d and not any(s.key == "D" for s in _d.strategies):
        _d.strategies.append(Strategy(
            "D", "İdame yok (yalnızca konak bağımlılığı)",
            list(_kaset) + [_ori], "ECN_DAPA",
            "Karşılaştırma için: plazmiti tutan hiçbir seçilim yok. "
            "Antibiyotiksiz pasajlarda plazmit kaybolur. Konağın ΔdapA "
            "bağımlılığı bakterinin bağırsak dışında yaşamasını "
            "engeller ama plazmiti tutmaz - ikisi farklı şeydir.\n\n"
            "**Denetim burada CORE-02 hatası verir ve bu beklenen:** "
            "markersiz bir plazmitte transformasyondan sonra doğru "
            "koloniyi seçmenin yolu bile yoktur."))

# P. aeruginosa devresine antibiyotiksiz idame basamağı.
_psa = by_key("psa")
if _psa and not any(s.key == "D" for s in _psa.strategies):
    _psa.strategies.append(Strategy(
        "D", "Algılayıcı + öldürücü, ccdAB idamesi (antibiyotiksiz)",
        PSA_KASET + ["ORI_PBR322"] + _CCDAB, "ECN_DAPA",
        "A stratejisinin klinik hedefli sürümü: kanamisin markeri "
        "yerine toksin-antitoksin idamesi, konak da laboratuvar suşu "
        "yerine terapötik Nissle ΔdapA. Düşük kopyalı pBR322 orijini "
        "seçildi: iki kasetli bir devre yüksek kopyada hücreye ağır "
        "gelir ve algılayıcının arka plan sinyali yükselir."))


# ---------------------------------------------------------------------------
# Ürat, çölyak ve oral mukozit: üç yeni kayıt
#
# Üçü de aynı soruyu farklı yerden sorar: eksik olan ENZİM mi, enzimin
# substratına ERİŞİMİ mi, yoksa ürünün DOKUYA ULAŞMASI mı? Ürat
# kasetinde taşıyıcı olmadan enzim boşa çalışır, çölyakta enzimin
# peptitten önce bağırsağa varması gerekir, mukozitte protein
# tükürükle yıkandığı için bakterinin yerinde durması gerekir.
# ---------------------------------------------------------------------------

#: Hiperürisemi. Kaset kaynak makalenin son suşundaki GEN SIRASINI
#: izler (pucL-pucM-vhb-ygfU-katG) ve promotörü de makalenin
#: kullandığı Ptrc'dir. Kendi önerimiz (oksijensizlikte açılan PfnrS)
#: ayrı bir basamak olarak duruyor, kasete karıştırılmadı.
URAT_KASET = ["PRO_TRC", "RBS_B0034", "CDS_PUCL", "RBS_B0032", "CDS_PUCM",
              "RBS_B0032", "CDS_VHB", "RBS_B0032", "CDS_UACT",
              "RBS_B0032", "CDS_KATG", "TER_B0015"]

#: Çölyak. Salgılanan sürüm makalenin daha etkili bulduğu sürümdür;
#: hücre içinde tutulan sürüm D basamağında karşılaştırma olarak var.
COLYAK_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
                "CDS_PEP_MX", "TER_B0015"]

#: Oral mukozit. AG013'ün kaseti kromozoma yerleştirilmiştir, plazmite
#: değil; buradaki plazmit sürümü aracın yapabildiği şeydir ve fark
#: "bilinen açık" notunda açıkça yazılıdır.
MUKOZIT_KASET = ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "SIG_USP45",
                 "CDS_TFF1", "TER_B0015"]


DISEASES += [

    Disease(
        key="urat",
        name="Hiperürisemi ve gut",
        short="Hiperürisemi",
        logic="Bağırsakta ürik asidi allantoine parçalayan bir bakteri: "
              "insanın evrimde kaybettiği enzim basamağını geri koymak.",
        mechanism=(
            "Pürin yıkımı çoğu canlıda allantoine kadar gider; insanda "
            "gitmez. Ürat oksidaz geni hominoid soyunda sözde-gene "
            "dönüşmüştür, bu yüzden zincir ürik asitte durur ve son ürün "
            "olarak ürik asit atılır. Ürik asidin çözünürlüğü ise "
            "düşüktür: serum derişimi yaklaşık 6,8 mg/dL'yi geçtiğinde "
            "monosodyum ürat kristalleri eklemde çöker, gut krizini "
            "başlatan iltihap budur.\n\n"
            "Ürik asidin yaklaşık üçte ikisi böbrekten, kalanı "
            "bağırsaktan atılır. Hiperüriseminin en sık nedeni üretim "
            "fazlalığı değil, **böbrek atılımının azalmasıdır** - yani "
            "sorun kapasitededir. Fikir buradan çıkar: böbreği zorlamak "
            "yerine bağırsak yolunu güçlendirmek. Bağırsakta ürik asidi "
            "parçalayan bir bakteri, kandan bağırsağa geçen üratı "
            "tüketerek derişim eğimini atılım yönünde tutar.\n\n"
            "Damardan verilen ürikaz ilaçları (pegloticase gibi) etkilidir "
            "ama bağışıklık yanıtı ve infüzyon reaksiyonu sorunu vardır. "
            "Bağırsakta yerinde üretilen enzim kana karışmaz; aradaki "
            "fark budur."),
        papers=[
            "Zhao R, Li Z, Sun Y, Ge W, Wang M, Liu H, Xun L, Xia Y. "
            "Engineered Escherichia coli Nissle 1917 with urate oxidase "
            "and an oxygen-recycling system for hyperuricemia treatment. "
            "Gut Microbes 14:2070391 (2022). "
            "doi:10.1080/19490976.2022.2070391 · PMID 35491895",
        ],
        pathway=_gv(
            'a [label="Pürin yıkımı\\n(insanda ürat oksidaz yok)"];',
            'b [label="Ürik asit birikir\\nçözünürlüğü düşük"];',
            'c [label="UacT / YgfU\\nüratı hücreye alır"];',
            'd [label="PucL N-uç alanı\\nürat + O2 → HIU + H2O2"];',
            'e [label="PucM\\nHIU → OHCU"];',
            'f [label="PucL C-uç alanı\\nOHCU → (S)-allantoin"];',
            'g [label="KatG\\n2 H2O2 → O2 + 2 H2O", fillcolor="#E7F5E9"];',
            'h [label="Vhb\\naz oksijende çalışma", fillcolor="#E7F5E9"];',
            "a -> b -> c -> d -> e -> f;",
            'd -> g [label="H2O2"];',
            'g -> d [label="O2 geri", style=bold];',
            "h -> d [style=dashed];"),
        pathway_caption=(
            "Zincirin üç basamağı var ve ikisini aynı protein yapıyor: "
            "PucL'in N-uç alanı ürat oksidaz, C-uç alanı OHCU "
            "dekarboksilazdır; aradaki basamağı PucM kapatır. Kaynak "
            "makale ikinci ve üçüncü basamağın kendiliğinden de "
            "yürüyebildiğini, ama birincisinin vazgeçilmez olduğunu "
            "söyler.\n\n"
            "**Diyagramdaki kalın ok tasarımın can alıcı yeridir:** "
            "katalaz yalnızca zehir temizlemiyor, ürat oksidazın "
            "harcadığı oksijeni geri veriyor. Kalın bağırsak neredeyse "
            "anaerobiktir, yani orada sınırlayıcı olan enzim miktarı "
            "değil oksijendir - bu yüzden katalaz ve Vhb, kasete "
            "'ihtiyaten' eklenmiş yardımcılar değil, işin yürümesinin "
            "koşuludur.\n\n"
            "Soldaki iki kutu insanın kendisine ait: tasarımın yaptığı "
            "şey, insanda kopmuş olan zinciri bakteride tamamlamaktır."),
        cassette=URAT_KASET,
        biomarker=("Serum ürik asit (mg/dL); 24 saatlik idrar ürat "
                   "atılımı; dışkıda ürat/allantoin oranı (bakterinin "
                   "gerçekten çalıştığını gösteren doğrudan ölçüm)"),
        gap=("**Kasetin boyu:** beş gen, yalnızca kodlayan bölgeler "
             "5901 bp. Bu, bu araçtaki en büyük terapötik kasettir ve "
             "metabolik yükü yüksektir; makale de katalaz ile Vhb "
             "eklendiğinde büyüme hızının düştüğünü (son hücre "
             "yoğunluğunun düşmediğini) ölçmüştür. Yüksek kopyalı bir "
             "orijinle kurulmamalıdır.\n\n"
             "**Promotör sorunu:** Ptrc makalenin kullandığı promotördür "
             "ama bir lac operatörü taşır, yani lacI taşıyan bir konakta "
             "indükleyici gerekir. Bağırsakta IPTG yoktur. Denetim bunu "
             "kendiliğinden yakalar: Ptrc'li basamaklarda **REG-01** "
             "uyarısını (LacI baskılaması yetersiz kalabilir) görürsünüz. "
             "E basamağı bu yüzden PfnrS'li sürümü ayrıca sunuyor - bu "
             "bizim önerimizdir, makalede sınanmamıştır.\n\n"
             "**Orijin:** makale pBBR1MCS-2 (orta kopyalı, geniş konaklı) "
             "kullanır; katalogda pBBR1 replikonu yok, onun yerine düşük "
             "kopyalı p15A kullanıldı. Kopya sayısı farkı ifade düzeyini "
             "de değiştirir.\n\n"
             "**Dizi farkı:** makale kodon optimize edilmiş ve iki nokta "
             "mutasyonu taşıyan bir pucL kullanır; buradaki dizi yabanıl "
             "*B. subtilis* 168 pucL'idir (ayrıntı parçanın notunda).\n\n"
             "**En büyük bilinmeyen tasarımda değil, fizyolojide:** ürik "
             "asidin bağırsak lümenine geçişi ABCG2 taşıyıcısına bağlıdır. "
             "Bakteri ne kadar hızlı yıkarsa yıksın, lümene gelen ürat "
             "miktarı sınırlıysa serum etkisi de sınırlı kalır."),
        strategies=_ecoli_merdiven(
            URAT_KASET, "ORI_P15A",
            "Terapötik Nissle: makalenin konağı da E. coli Nissle "
            "1917'dir, buradaki suş onun ΔdapA içerme sürümüdür."),
        comparison=[
            {"Özellik": "Konak", "Makale": "E. coli Nissle 1917",
             "Bizim tasarım": "E. coli Nissle 1917 ΔdapA",
             "Etiket": "Literatür"},
            {"Özellik": "Genler",
             "Makale": "pucL + pucM (B. subtilis), ygfU + katG "
                       "(E. coli), vhb (Vitreoscilla)",
             "Bizim tasarım": "Aynı beş gen; uacT = ygfU (eş ad), "
                              "Vhb dizisi WP_019959060.1 ile birebir",
             "Etiket": "Literatür"},
            {"Özellik": "pucL sürümü",
             "Makale": "Kodon optimize + Asp44Val, Gln268Arg (pucL-M)",
             "Bizim tasarım": "Yabanıl B. subtilis 168 pucL (O32141)",
             "Etiket": "Fark"},
            {"Özellik": "Promotör", "Makale": "Ptrc",
             "Bizim tasarım": "Ptrc (A-D) · PfnrS (E basamağı)",
             "Etiket": "Öneri"},
            {"Özellik": "İskelet", "Makale": "pBBR1MCS-2 (orta kopya)",
             "Bizim tasarım": "p15A orijini (düşük kopya)",
             "Etiket": "Fark"},
            {"Özellik": "İdame",
             "Makale": "Plazmit markeri",
             "Bizim tasarım": "ccdAB toksin-antitoksin (A) ya da dapA "
                              "tamamlama (B)",
             "Etiket": "Öneri"},
            {"Özellik": "Model",
             "Makale": "Farede damardan ürik asit verilerek kurulan "
                       "yeni akut hiperürisemi modeli",
             "Bizim tasarım": "—", "Etiket": "Literatür"},
        ],
    ),

    Disease(
        key="colyak",
        name="Çölyak hastalığı (glutene bağlı enteropati)",
        short="Çölyak",
        logic="Glutenin insanın kesemediği bağını kesen bir enzimi, "
              "gıda sınıfı bir bakteriyle ince bağırsakta yerinde "
              "üretmek.",
        mechanism=(
            "Glutenin gliadin kısmı prolin ve glutamin bakımından "
            "olağandışı zengindir. İnsanın mide, pankreas ve fırça "
            "kenarı proteazları prolinden sonraki bağı verimli "
            "kesemez; bu yüzden sindirimden uzun, prolinden zengin "
            "peptitler **sağ çıkar**. En iyi çalışılmışı alfa-2 "
            "gliadinin 33 amino asitlik parçasıdır (33-mer).\n\n"
            "Bu peptit bağırsak epitelini geçtiğinde doku "
            "transglutaminazı (TG2) belirli glutaminleri glutamata "
            "çevirir (deamidasyon). Deamide peptit HLA-DQ2 ve DQ8'e "
            "yüksek ilgiyle bağlanır, CD4+ T hücreleri uyarılır ve "
            "ortaya çıkan iltihap villusları düzleştirir - emilim "
            "yüzeyi kaybolur.\n\n"
            "Tek tedavi ömür boyu glutensiz diyettir ve uyum gerçekte "
            "zordur: gluten, buğday içermeyen ürünlerde bile bulaşma "
            "yoluyla karşımıza çıkar. Fikir, eksik olan enzim "
            "etkinliğini yerine koymaktır. Ağızdan enzim vermenin "
            "sorunu mideyi geçmektir; enzimi bağırsakta ÜRETEN bir "
            "bakteri bu sorunu baştan atlar."),
        papers=[
            "Shan L, Molberg Ø, Parrot I, Hausch F, Filiz F, Gray GM, "
            "Sollid LM, Khosla C. Structural basis for gluten "
            "intolerance in celiac sprue. Science 297:2275-2279 (2002). "
            "doi:10.1126/science.1074129 · PMID 12351792",
            "Alvarez-Sieiro P, Martín MC, Redruello B, del Río B, "
            "Ladero V, Palanski BA, Khosla C, Fernández M, Álvarez MA. "
            "Generation of food-grade recombinant Lactobacillus casei "
            "delivering Myxococcus xanthus prolyl endopeptidase. "
            "Applied Microbiology and Biotechnology 98:6689-6700 "
            "(2014). doi:10.1007/s00253-014-5730-7 · PMID 24752841",
        ],
        pathway=_gv(
            'a [label="Gluten / gliadin\\nprolin + glutamin zengini"];',
            'b [label="İnsan proteazları\\nprolinden sonra kesemez"];',
            'c [label="33-mer sağ çıkar\\n(alfa-2 gliadin)"];',
            'd [label="TG2 deamidasyonu\\nGln → Glu"];',
            'e [label="HLA-DQ2/DQ8\\nsunum"];',
            'f [label="CD4+ T yanıtı\\nvillus atrofisi"];',
            'g [label="Prolil endopeptidaz\\n(M. xanthus PEP)", '
            'fillcolor="#E7F5E9"];',
            "a -> b -> c -> d -> e -> f;",
            'g -> c [label="33-mer\'i keser", style=bold];'),
        pathway_caption=(
            "Zincirin kırılma noktası **c kutusudur**: hastalık, "
            "peptitin sağ çıkmasıyla başlar. Enzim oraya müdahale eder, "
            "yani aşağı akıştaki hiçbir basamağı (TG2, HLA sunumu, T "
            "hücresi) hedeflemez - substratı ortadan kaldırır.\n\n"
            "Bunun pratik sonucu bir zamanlama kısıtıdır: enzimin "
            "peptiti, peptit epiteli geçmeden önce kesmesi gerekir. "
            "Yani etkinlik yalnızca enzimin katalitik hızına değil, "
            "bağırsak geçiş süresine ve enzimin o bölgede ne kadar "
            "bulunduğuna bağlıdır. Bir enzimin tüpte 33-mer'i kesmesi, "
            "hastada hasarı önlediğini göstermez.\n\n"
            "**b kutusu neden önemli:** burada bozuk olan insan "
            "proteazı değil; insanın hiç sahip olmadığı bir kesim "
            "yeteneği söz konusu. Bu yüzden tasarım bir eksiği "
            "tamamlamıyor, yeni bir etkinlik ekliyor."),
        cassette=COLYAK_KASET,
        biomarker=("Doku transglutaminaz IgA (anti-TG2) ve deamide "
                   "gliadin peptit antikorları; duodenum biyopsisinde "
                   "Marsh evresi. Enzimin kendisi için: kütle "
                   "spektrometresiyle artık 33-mer ölçümü"),
        gap=("**Konak farkı:** makale *Lactobacillus casei* kullanır "
             "(bugünkü adıyla *Lacticaseibacillus casei*); katalogda bu "
             "tür yok, onun yerine gıda sınıfı *Lactococcus lactis* "
             "NZ9000 kullanıldı. İkisi de laktik asit bakterisidir ama "
             "mide-bağırsak dayanıklılıkları aynı değildir ve makalenin "
             "benzetilmiş sindirim sistemi testi *L. casei* ile "
             "yapılmıştır.\n\n"
             "**Gram-pozitif replikon açığı (CORE-01):** A basamağında "
             "denetim hata verir, çünkü katalogda bu bakteride çoğalan "
             "bir orijin yok. Açık iskelettedir, kasette değil; "
             "Referanstan kes-yapıştır sayfasındaki pNZ8048 tarifi bu "
             "açığı kapatmak için duruyor.\n\n"
             "**Kayıt durumu:** PEP dizisi TrEMBL'dendir, yani elle "
             "incelenmemiştir. Çevirisi doğrulanmıştır ama küratörlüğü "
             "yapılmamıştır.\n\n"
             "**Kanıt düzeyi:** makale 33-mer'in benzetilmiş ince "
             "bağırsak ortamında yıkıldığını gösterir. Hastada klinik "
             "etkinlik gösterilmemiştir ve yazarlar da bunu açıkça "
             "yazar. Gluten sindirimi için geliştirilen enzimlerin "
             "klinik denemelerinde genel tablo şimdiye kadar "
             "karışıktır."),
        strategies=_gram_poz_merdiven(
            COLYAK_KASET,
            ["PRO_T7", "RBS_PET", "CDS_PEP_MX", "TER_T7"],
            "LLNZ9000",
            "Gıda sınıfı laktik asit bakterisi: enzimi yiyecekle "
            "birlikte alınıp ince bağırsakta üretmesi beklenir."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Lactobacillus casei BL23",
             "Bizim tasarım": "Lactococcus lactis NZ9000",
             "Etiket": "Fark"},
            {"Özellik": "Enzim",
             "Makale": "M. xanthus prolil endopeptidaz",
             "Bizim tasarım": "M. xanthus prolil endopeptidaz "
                              "(Q9X5N2, 689 aa)",
             "Etiket": "Literatür"},
            {"Özellik": "Konum",
             "Makale": "İki suş karşılaştırıldı: salgılanan ve hücre "
                       "içinde tutulan",
             "Bizim tasarım": "Salgılanan (A-C, Usp45 sinyali) ve "
                              "hücre içi (D) - aynı karşılaştırma",
             "Etiket": "Literatür"},
            {"Özellik": "Sonuç",
             "Makale": "Salgılanan suş 33-mer'i daha etkili yıktı ve "
                       "sindirim stresine dayandı",
             "Bizim tasarım": "Varsayılan basamak bu yüzden salgılanan "
                              "sürüm",
             "Etiket": "Literatür"},
            {"Özellik": "Gıda sınıfı seçilim",
             "Makale": "Gıda sınıfı sistem (antibiyotik markeri yok)",
             "Bizim tasarım": "thyA içermesi (timin bağımlılığı)",
             "Etiket": "Öneri"},
            {"Özellik": "Model",
             "Makale": "Benzetilmiş ince bağırsak ortamı (in vitro)",
             "Bizim tasarım": "—", "Etiket": "Literatür"},
        ],
    ),

    Disease(
        key="mukozit",
        name="Oral mukozit (ışın ve kemoterapiye bağlı)",
        short="Oral mukozit",
        logic="Ağız mukozasının üstünde durup iyileşmeyi hızlandıran "
              "peptidi sürekli salgılayan bir bakteri - gargara olarak.",
        mechanism=(
            "Baş-boyun ışın tedavisi ve yüksek doz kemoterapi, hızlı "
            "bölünen hücreleri ayırt etmeden vurur. Ağız epitelinin "
            "bazal hücreleri de bunlardandır: epitel incelir, sonra "
            "ülserleşir. Sonuç, hastanın yiyememesi, yüksek doz "
            "ağrı kesici gereksinimi, beslenme tüpü ve enfeksiyon "
            "riskidir.\n\n"
            "Buradaki asıl sorun konforun ötesindedir: ağrı yüzünden "
            "kanser tedavisine ara verilir ya da dozu düşürülür, bu da "
            "kanser sonucunu kötüleştirir. Yani mukozit, tedavinin "
            "kendisini sınırlayan bir yan etkidir.\n\n"
            "Trefoil peptitleri (TFF1/2/3) mukus jelinde duran, "
            "proteaza dirençli küçük peptitlerdir ve epitel onarımında "
            "**restitüsyonu** hızlandırır: hücreler yaranın üstüne göç "
            "eder, bariyer bölünmeyi beklemeden kapanır. Bölünmenin "
            "ışınla durdurulduğu bir dokuda işe yarayabilecek mekanizma "
            "tam olarak budur.\n\n"
            "Proteini doğrudan sürmenin sorunu basit: tükürük onu "
            "dakikalar içinde yıkar. Fikir, peptidi mukozanın üstünde "
            "duran canlı bir bakteriye sürekli ürettirmektir."),
        papers=[
            "Caluwaerts S, Vandenbroucke K, Steidler L, Neirynck S, "
            "Vanhoenacker P, Corveleyn S, Watkins B, Sonis S, Coulie B, "
            "Rottiers P. AG013, a mouth rinse formulation of "
            "Lactococcus lactis secreting human Trefoil Factor 1, "
            "provides a safe and efficacious therapeutic tool for "
            "treating oral mucositis. Oral Oncology 46:564-570 (2010). "
            "doi:10.1016/j.oraloncology.2010.04.008 · PMID 20542722",
            "Vandenbroucke K, Hans W, Van Huysse J, Neirynck S, "
            "Demetter P, Remaut E, Rottiers P, Steidler L. Active "
            "delivery of trefoil factors by genetically modified "
            "Lactococcus lactis prevents and heals acute colitis in "
            "mice. Gastroenterology 127:502-513 (2004). "
            "doi:10.1053/j.gastro.2004.05.020 · PMID 15300583",
        ],
        pathway=_gv(
            'a [label="Işın / kemoterapi"];',
            'b [label="Bazal keratinosit\\nbölünmesi durur"];',
            'c [label="Epitel incelir\\nve ülserleşir"];',
            'd [label="Ağrı, yiyememe\\nkanser tedavisine ara"];',
            'e [label="L. lactis mukozada\\nhTFF1 salgılar", '
            'fillcolor="#E7F5E9"];',
            'f [label="Restitüsyon\\nhücre göçüyle kapanma", '
            'fillcolor="#E7F5E9"];',
            "a -> b -> c -> d;",
            "e -> f;",
            'f -> c [label="bariyeri kapatır", style=bold];'),
        pathway_caption=(
            "Zincirin solu tasarımın değiştirebileceği bir şey değil: "
            "ışının bölünmeyi durdurması tedavinin kendisidir. "
            "Müdahale noktası **c kutusudur** ve seçilen mekanizma "
            "bunu belirler.\n\n"
            "Restitüsyonun burada seçilmesinin nedeni şudur: hücre "
            "bölünmesini uyaran bir büyüme faktörü, kanser tedavisi "
            "gören bir hastada istenmeyen bir şeydir. Trefoil "
            "peptitleri bölünmeyi değil GÖÇÜ hızlandırır, yani "
            "bariyeri yeni hücre yapmadan kapatır. Mekanizma seçimi, "
            "hastanın kanser hastası olmasının doğrudan sonucudur.\n\n"
            "Diyagramda bakteri zincirin içinde değil, yanında duruyor: "
            "bu tasarım bağırsağı kolonize etmeye çalışmaz, mukozanın "
            "üstünde geçici olarak kalır."),
        cassette=MUKOZIT_KASET,
        biomarker=("WHO ya da RTOG/NCI mukozit evresi; ülser yüzey "
                   "alanı; opioid kullanımı ve beslenme tüpü gereken "
                   "gün sayısı (hastanın gerçekten hissettiği sonuç)"),
        gap=("**En önemli fark plazmitin kendisidir.** AG013'ün suşu "
             "(sAGX0085) htff1 kasetini KROMOZOMUNA yerleştirmiştir, "
             "plazmite değil. Nedeni açık: klinik bir ürün için "
             "plazmit, seçilim baskısı ve kararlılık sorunu demektir. "
             "Bu araç plazmit tasarlar, dolayısıyla buradaki tasarım "
             "ürünün kendisi değil, ona giden yoldaki laboratuvar "
             "sürümüdür. Bir ürün yapılacaksa doğru cevap plazmit "
             "değildir - bunu gizlemek yerine yazıyoruz.\n\n"
             "**İçerme merdiveni burada yanlış soruyu yanıtlıyor:** "
             "thyA ve dapA bağımlılığı, bağırsakta kalıcı olmayı "
             "hedefleyen bakteriler için düşünülmüştür. AG013 ağızda "
             "topikal verilir ve kalıcı olması beklenmez; makalenin "
             "güvenlik çalışmaları suşun nötropenik koşullarda bile "
             "dolaşıma geçmediğini gösterir. Burada asıl soru "
             "kalıcılık değil, temas süresidir.\n\n"
             "**Gram-pozitif replikon açığı (CORE-01):** A basamağında "
             "denetim hata verir; katalogda L. lactis'te çoğalan bir "
             "orijin yok. pNZ8048 tarifi bu açığı kapatır.\n\n"
             "**Katlanma:** hTFF1 üç disülfür bağı taşır. "
             "L. lactis'te Usp45 sinyaliyle salgılandığında bağlar "
             "hücre dışında kurulur; E. coli'de ön üretim yapılırken "
             "sitoplazmada yanlış katlanır, bu yüzden B ve C "
             "basamaklarında DsbA sinyali gerekir."),
        strategies=_gram_poz_merdiven(
            MUKOZIT_KASET,
            ["PRO_T7", "RBS_PET", "SIG_DSBA", "CDS_TFF1", "TER_T7"],
            "LLNZ9000",
            "Gıda sınıfı L. lactis: ağız mukozasına gargara olarak "
            "verilir ve peptidi yerinde üretir."),
        comparison=[
            {"Özellik": "Konak", "Makale": "Lactococcus lactis sAGX0085",
             "Bizim tasarım": "Lactococcus lactis NZ9000",
             "Etiket": "Literatür"},
            {"Özellik": "Terapötik",
             "Makale": "İnsan trefoil faktör 1 (hTFF1)",
             "Bizim tasarım": "hTFF1 olgun zinciri (P04155, 61 aa)",
             "Etiket": "Literatür"},
            {"Özellik": "Kasetin yeri",
             "Makale": "Kromozoma kararlı yerleştirme",
             "Bizim tasarım": "Plazmit (araç plazmit tasarlar)",
             "Etiket": "Fark"},
            {"Özellik": "Veriliş yolu",
             "Makale": "Gargara (AG013), topikal",
             "Bizim tasarım": "Aynı kullanım hedeflenir",
             "Etiket": "Literatür"},
            {"Özellik": "İdame",
             "Makale": "Gerekmez (kromozomda)",
             "Bizim tasarım": "thyA içermesi ya da ön üretimde "
                              "antibiyotik markeri",
             "Etiket": "Fark"},
            {"Özellik": "Model",
             "Makale": "Hamsterda ışına bağlı oral mukozit; nötropenik "
                       "hamster ve sıçanda güvenlik",
             "Bizim tasarım": "—", "Etiket": "Literatür"},
        ],
    ),

]

# Çölyak kaydına dördüncü basamak: sinyal peptidini ÇIKARMAK. Makale
# iki suşu yan yana sınadı (salgılanan ve hücre içinde tutulan) ve
# salgılananı daha etkili buldu; o karşılaştırma buraya basamak olarak
# taşındı, çünkü "sinyal peptidi gerçekten gerekli mi?" sorusunun
# deneysel yanıtı olan ender durumlardan biri.
_colyak = by_key("colyak")
if _colyak and not any(s.key == "D" for s in _colyak.strategies):
    _colyak.strategies.append(Strategy(
        "D", "Enzim hücre içinde tutulur (sinyal peptidi yok)",
        ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034", "CDS_PEP_MX",
         "TER_B0015", "MRK_THYA_LL"], "LLNZ9000",
        "A basamağıyla aynı kaset, tek farkı Usp45 sinyalinin "
        "çıkarılmış olması: enzim sitoplazmada kalır. Makalenin "
        "sınadığı ikinci suş budur ve sonuç nettir - salgılanan suş "
        "33-mer'i daha etkili yıktı.\n\n"
        "Nedeni mekaniktir: 33-mer bağırsak lümeninde durur, bakterinin "
        "içine girmez. Enzim hücre içinde kalırsa substratına hiç "
        "ulaşamaz; tek şansı hücrenin parçalanıp enzimi salmasıdır, ki "
        "bu da denetlenemez. Ürat kasetindeki durumun tam tersi: orada "
        "substrat (ürik asit) bir taşıyıcıyla içeri alınabiliyordu, "
        "burada 33 amino asitlik bir peptit için öyle bir yol yok.\n\n"
        "Basamak merdivende duruyor çünkü sinyal peptidinin gerekli "
        "olup olmadığı genellikle tartışılır, burada ise deneysel "
        "yanıtı var.\n\n"
        "**Denetim burada da CORE-01 verir:** Gram-pozitif replikon "
        "açığı bu basamakta da duruyor."))


# Ürat kasetine beşinci basamak: promotörü değiştirmek. Ptrc bir lac
# operatörü taşır ve bağırsakta IPTG yoktur; PfnrS ise oksijensizlikte
# kendiliğinden açılır - yani kasetin açılma koşulu, kasetin çalışması
# gereken yerin kendisi olur. Bu bizim önerimizdir ve makalede
# sınanmamıştır, o yüzden varsayılan basamak değil ek basamak.
_urat = by_key("urat")
if _urat and not any(s.key == "E" for s in _urat.strategies):
    _urat.strategies.append(Strategy(
        "E", "PfnrS promotörü: bağırsakta kendiliğinden açılan sürüm",
        ["PRO_FNRS"] + URAT_KASET[1:] + ["ORI_P15A"] + _CCDAB, "ECN_DAPA",
        "A basamağıyla aynı kaset, tek farkı promotör. Ptrc güçlüdür "
        "ama bir lac operatörü taşır; lacI taşıyan bir konakta "
        "indükleyici gerekir ve bağırsakta IPTG yoktur. PfnrS, FNR "
        "düzenleyicisi aracılığıyla oksijen azaldığında açılır - "
        "yani kaset tam olarak çalışması gereken yerde, kalın "
        "bağırsağın anaerobik ortamında açılır ve laboratuvar "
        "koşullarında kapalı durur.\n\n"
        "Bedeli şudur: PfnrS'nin gücü Ptrc'den düşüktür ve beş genlik "
        "bir operonda bu, en uçtaki genin (katG) ifadesini "
        "sınırlayabilir. Oksijeni geri veren enzimi kısmak, oksijen "
        "kıtlığında çalışacak bir tasarımda ters etki yapar - "
        "basamağın bir karşılaştırma olarak durmasının nedeni budur.\n\n"
        "**Bu basamak bizim önerimizdir, makalede sınanmamıştır.**"))
