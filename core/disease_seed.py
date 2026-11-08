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
}


@functools.lru_cache(maxsize=1)
def disease_parts() -> tuple[dict[str, Any], ...]:
    out: list[dict[str, Any]] = []
    for rec in _raw().get("parts", []):
        incelenmis = rec.get("reviewed", True)
        notlar = [GENE_NOTES.get(rec["part_id"], "")]
        notlar.append(
            f"Dizi {rec['organism']} kaynaklıdır: UniProt "
            f"{rec['uniprot']} ({rec['uniprot_entry']}), nükleotid "
            f"{rec['nucleotide']}. {rec['length_aa']} amino asit / "
            f"{rec['length']} bp; çevirisinin UniProt'taki proteine eşit "
            f"olduğu doğrulandı."
            + ("" if incelenmis else
               " Kayıt TrEMBL durumundadır, yani otomatik üretilmiş ve "
               "elle incelenmemiştir."))
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
                     "İçerme korunur ama plazmit idamesi için yöntem yok."),
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
                     "direnci taşımaz - terapötik kullanımın şartı budur."),
            Strategy("B", "Kloramfenikol markerlı prototip",
                     IBD_KASET + ["REF_CAT_651"], "LLNZ9000",
                     "pNZ8048'den çıkarılmış cat markeri. Laboratuvarda "
                     "kurması kolay, hastada kullanılamaz."),
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
