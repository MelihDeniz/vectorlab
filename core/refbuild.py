"""Gerçek bir referans vektörü kesip yapıştırarak yeni vektör kurar.

Bu modül Biopython'a bağlı değildir; tanıma yerlerini ve kesim
konumlarını kendisi bulur, bu yüzden tek başına sınanabilir. Streamlit
sayfası (app_pages/refbuild.py) bunu kullanır ve sonucu uygulamanın
kendi denetleyicisi, haritası ve jel simülasyonuyla gösterir.

Kurulan her adım şunu yapar: iki farklı enzimin tek kesen yerleri
arasındaki bölge atılır, yerine insert yapıştırılır, iki kesim yeri
yeniden oluşturulur.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Kullandığımız enzimler: tanıma dizisi, kesimin sitenin kaçıncı
# nükleotidinden sonra olduğu (üst zincir) ve bıraktığı uç.
# Kaynak: REBASE / üretici katalogları (standart tanıma dizileri).
# Hepsi 5' cikinti birakir; cut_and_paste yapisi buna gore kurulmustur.
# tests/test_refbuild.py bu tablonun REBASE ile tutarli oldugunu sinar,
# yani elle yazilmis olmasina ragmen dogrulanmamis degildir.
ENZYMES: dict[str, dict[str, Any]] = {
    "HindIII": {"site": "AAGCTT", "cut": 1, "ovh": "5' AGCT"},
    "SalI":    {"site": "GTCGAC", "cut": 1, "ovh": "5' TCGA"},
    "BclI":    {"site": "TGATCA", "cut": 1, "ovh": "5' GATC",
                "dam": True},
    "EcoRI":   {"site": "GAATTC", "cut": 1, "ovh": "5' AATT"},
    "BamHI":   {"site": "GGATCC", "cut": 1, "ovh": "5' GATC"},
    "XbaI":    {"site": "TCTAGA", "cut": 1, "ovh": "5' CTAG",
                "dam": True},
    "SpeI":    {"site": "ACTAGT", "cut": 1, "ovh": "5' CTAG"},
    "XhoI":    {"site": "CTCGAG", "cut": 1, "ovh": "5' TCGA"},
    "NcoI":    {"site": "CCATGG", "cut": 1, "ovh": "5' CATG"},
}


@dataclass
class Step:
    """Tek bir kes-yapıştır adımının sonucu."""
    title: str
    left: str
    right: str
    sequence: str
    length: int
    left_site_pos: int        # 1 tabanlı, vektördeki tanıma yerinin başı
    right_site_pos: int
    removed: int              # atılan bölge (siteler dahil değil)
    insert_length: int
    insert_start: int         # ürün üzerinde insertin ilk nt'si (1 tabanlı)
    insert_end: int
    # (eski_baş, eski_son, yeni_baş): korunan bölgelerin koordinat eşlemesi
    segments: list[tuple[int, int, int]] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def find_sites(seq: str, enzyme: str) -> list[int]:
    """Dairesel dizide tanıma yerinin 1 tabanlı başlangıçları."""
    site = ENZYMES[enzyme]["site"]
    s = seq.upper()
    ext = s + s[: len(site) - 1]
    out = []
    i = ext.find(site)
    while i != -1:
        if i < len(s):
            out.append(i + 1)
        i = ext.find(site, i + 1)
    return out


def single_cutter(seq: str, enzyme: str) -> int:
    pos = find_sites(seq, enzyme)
    if len(pos) != 1:
        raise ValueError(f"{enzyme} vektörü {len(pos)} yerden kesiyor, "
                         f"tek kesen olması gerekir")
    return pos[0]


def cut_and_paste(seq: str, left: str, right: str, insert: str,
                  title: str = "") -> Step:
    """left ve right enzimlerinin tek kesen yerleri arasını insertle değiştir.

    Bölge left sitesinden right sitesine doğru (saat yönünde, gerekirse
    dizinin başını aşarak) alınır. İnsert, iki enzimin uçlarına göre yönlüdür:
    left tarafı insertin 5' ucudur.
    """
    s = seq.upper()
    L = len(s)
    ls = ENZYMES[left]["site"]
    rs = ENZYMES[right]["site"]
    a = single_cutter(s, left)          # 1 tabanlı
    b = single_cutter(s, right)
    a0 = a - 1                          # 0 tabanlı sitenin başı
    b0 = b - 1
    b_end = b0 + len(rs)                # right sitesinin sonrası (0 tabanlı)

    if b0 > a0:
        # Bölge dizinin başını aşmıyor.
        removed = b0 - (a0 + len(ls))
        new = s[:a0] + ls + insert.upper() + rs + s[b_end:]
        ins_start = a0 + len(ls) + 1
        segs = [(1, a0 + len(ls), 1),
                (b0 + 1, L, a0 + len(ls) + len(insert) + 1)]
        order = "dogrusal"
    else:
        # Bölge başlangıcı aşıyor: right sitesinden sonrası ... left sitesine
        # kadar olan kısım korunur; yeni dizi right sitesinden sonra başlar.
        kept = s[b_end:a0]
        removed = L - len(kept) - len(ls) - len(rs)
        new = kept + ls + insert.upper() + rs
        ins_start = len(kept) + len(ls) + 1
        segs = [(b_end + 1, a0 + len(ls), 1),
                (b0 + 1, b_end, len(new) - len(rs) + 1)]
        order = "sarar"

    st = Step(
        title=title or f"{left} + {right}",
        left=left, right=right, sequence=new, length=len(new),
        left_site_pos=a, right_site_pos=b, removed=removed,
        insert_length=len(insert),
        insert_start=ins_start, insert_end=ins_start + len(insert) - 1,
        segments=segs,
    )
    st.notes.append(order)
    # Kontroller: sonuçta iki enzim de yalnızca bir kez kesmeli.
    for e in (left, right):
        n = len(find_sites(new, e))
        if n != 1:
            st.notes.append(
                f"UYARI: {e} ürün içinde {n} yerden kesiyor; insert bu "
                f"enzimin yerini taşıyor olabilir.")
    return st


def map_position(step: Step, orig_pos: int) -> int | None:
    """Vektörün eski koordinatını (1 tabanlı) ürün koordinatına çevir.

    Atılan bölgenin içindeki konumlar için None döner.
    """
    for o_start, o_end, n_start in step.segments:
        if o_start <= orig_pos <= o_end:
            return n_start + (orig_pos - o_start)
    return None


def map_feature(step: Step, start: int, end: int,
                orig_length: int) -> tuple[int, int, str] | None:
    """Bir BÖLGENİN ürün koordinatındaki karşılığı (1 tabanlı, kapsayıcı).

    `start > end` bölgenin dizinin başını sardığı anlamına gelir; dönen
    değerde de aynı kural geçerlidir.

    Döner `(yeni_baş, yeni_son, durum)`:

    * `"tam"`   — bölge bütün ve kesintisiz korundu.
    * `"kesildi"` — bölgenin bir kısmı atıldı ya da ortasına insert
      girdi; yalnızca kalan en büyük parça döner.
    * `None`    — bölge tamamen atıldı.

    Bu ayrım süs değil. Bir markerin ya da lacZ'nin İÇİNE kaset takmak
    o geni bozar ve bozulması eleme ölçütüdür (insersiyonel
    inaktivasyon, mavi/beyaz tarama). `map_position` uçları ayrı ayrı
    taşıdığı için bozulan geni sapasağlam ve insertin üstüne binmiş
    gösterir; harita o zaman tasarımı yanlış anlatır.
    """
    araliklar = ([(start, end)] if start <= end
                 else [(start, orig_length), (1, end)])
    beklenen = sum(b - a + 1 for a, b in araliklar)

    parcalar: list[tuple[int, int]] = []
    for a0, b0 in araliklar:
        for o_start, o_end, n_start in step.segments:
            a = max(a0, o_start)
            b = min(b0, o_end)
            if a <= b:
                parcalar.append((n_start + (a - o_start),
                                 n_start + (b - o_start)))
    if not parcalar:
        return None

    # Ürün üzerinde uç uca gelen parçaları birleştir (saran bölge iki
    # parça olarak toplanır ama aslında tek parçadır).
    parcalar.sort()
    birlesik = [parcalar[0]]
    for a, b in parcalar[1:]:
        if a == birlesik[-1][1] + 1:
            birlesik[-1] = (birlesik[-1][0], b)
        else:
            birlesik.append((a, b))
    if (len(birlesik) > 1 and birlesik[0][0] == 1
            and birlesik[-1][1] == step.length):
        birlesik = [(birlesik[-1][0], birlesik[0][1])] + birlesik[1:-1]

    uzunluk = sum((b - a + 1) if a <= b else (step.length - a + 1 + b)
                  for a, b in birlesik)
    if len(birlesik) == 1 and uzunluk == beklenen:
        return birlesik[0][0], birlesik[0][1], "tam"
    a, b = max(birlesik,
               key=lambda p: (p[1] - p[0] + 1) if p[0] <= p[1]
               else (step.length - p[0] + 1 + p[1]))
    return a, b, "kesildi"


@dataclass
class Result:
    original_length: int
    steps: list[Step]
    final: str

    @property
    def final_length(self) -> int:
        return len(self.final)


def build(vector_seq: str, stage1_insert: str, stage2_insert: str) -> Result:
    """İki aşamalı kurulum: HindIII+SalI, sonra BclI+EcoRI."""
    s1 = cut_and_paste(vector_seq, "HindIII", "SalI", stage1_insert,
                       "1. aşama: PKU kasetini ekle")
    s2 = cut_and_paste(s1.sequence, "BclI", "EcoRI", stage2_insert,
                       "2. aşama: cat'i boz, ccdAB ekle")
    return Result(len(vector_seq), [s1, s2], s2.sequence)


def digest_sizes(seq: str, enzymes: list[str]) -> list[int]:
    """Dairesel diziyi verilen enzimlerle kes; fragman boyutları (büyükten küçüğe)."""
    L = len(seq)
    cuts = []
    for e in enzymes:
        for p in find_sites(seq, e):
            cuts.append((p - 1 + ENZYMES[e]["cut"]) % L)
    cuts = sorted(set(cuts))
    if not cuts:
        return [L]
    sizes = [(cuts[(i + 1) % len(cuts)] - c) % L or L
             for i, c in enumerate(cuts)]
    return sorted(sizes, reverse=True)


# ---------------------------------------------------------------------------
# Tarifler
# ---------------------------------------------------------------------------
# Her tarif gerçek bir referans vektörü, gerçek bir enzim çiftini ve
# gerçek bir terapötik kaseti birleştirir. Amaç yalnızca çeşitlilik
# değil: her tarif farklı bir şey öğretir ve sunumda ayrı bir soruya
# cevap verir.


@dataclass
class RecipeStep:
    left: str
    right: str
    insert_parts: list[str]
    title: str
    lesson: str = ""


@dataclass
class Recipe:
    key: str
    title: str
    vector_id: str
    disease: str              # disease_seed anahtarı ya da ""
    why_vector: str
    steps: list[RecipeStep]
    source: str = ""
    teaches: str = ""


def build_recipe(vector_seq: str,
                 steps: list[tuple[str, str, str, str]]) -> Result:
    """Sırayla kes-yapıştır uygula.

    `steps` her adım için (sol enzim, sağ enzim, insert dizisi, başlık)
    taşır. Her adım bir öncekinin ürünü üzerinde çalışır.
    """
    dizi = vector_seq
    sonuc: list[Step] = []
    for sol, sag, insert, baslik in steps:
        adim = cut_and_paste(dizi, sol, sag, insert, baslik)
        sonuc.append(adim)
        dizi = adim.sequence
    return Result(len(vector_seq), sonuc, dizi)


PKU_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_STLA", "RBS_B0032", "CDS_PHEP",
             "TER_B0015"]
CCD_MODUL = ["PRO_CCDAB", "MRK_CCDA", "MRK_CCDB", "TER_B0015"]

RECIPES: list[Recipe] = [
    Recipe(
        key="pacyc_pku",
        title="pACYC184 → PKU vektörü (iki aşama)",
        vector_id="PACYC184", disease="pku",
        why_vector=(
            "p15A replikonu düşük kopyalıdır ve ColE1 vektörleriyle bir "
            "arada durabilir; seçtiğimiz dört enzim vektörü tek yerden "
            "keser. İki antibiyotik direnci terapötik kullanım için zayıf "
            "noktadır ve revizyonu gerekçelendirir."),
        source="Chang & Cohen 1978, J Bacteriol 134:1141 (p15A vektörleri)",
        teaches=(
            "İki aşamalı kurulumda ikinci aşamanın enzimlerinin BİRİNCİ "
            "aşamanın ÜRÜNÜNDE seçilmesi. Bu tarif ilk kurulduğunda "
            "ikinci aşama BclI kullanıyordu; PKU kasetinin gerçek "
            "dizisi gelince BclI'nın üründe dört yerden kestiği "
            "görüldü (kasetin içinde üç tanıma yeri daha varmış) ve "
            "yerine NcoI geldi. Enzim vektörde tek kesiyor diye "
            "ürününde de tek keseceği varsayılamaz; uygulama bunu her "
            "aşamada yeniden hesaplar."),
        steps=[
            RecipeStep("HindIII", "SalI", PKU_KASET,
                       "1. aşama: PKU kasetini ekle",
                       "Tetrasiklin direncinin bir bölümü çıkar."),
            RecipeStep("NcoI", "EcoRI", CCD_MODUL,
                       "2. aşama: cat'i boz, ccdAB ekle",
                       "cat geninin 660 bp'sinden 219 bp'si yerinde kalır "
                       "ama kırpılmış CAT proteini işlev görmez: "
                       "kloramfenikol direnci biter, yerine "
                       "toksin-antitoksin idamesi gelir."),
        ],
    ),
    Recipe(
        key="puc19_amonyak",
        title="pUC19 → hiperammonemi vektörü (tek aşama)",
        vector_id="PUC19", disease="amonyak",
        why_vector=(
            "Yüksek kopyalı ve klonlaması en kolay vektör. Terapötik "
            "kullanım için uygun DEĞİLDİR (ampisilin direnci, çok yüksek "
            "kopya yükü); burada amaç kaseti hızlıca kurup doğrulamak."),
        source="Yanisch-Perron ve ark. 1985, Gene 33:103 (pUC19)",
        teaches=(
            "Mavi/beyaz tarama. Kaset lacZ-alfa'nın içindeki klonlama "
            "bölgesine girer, yani alfa komplementasyon bozulur: doğru "
            "klon BEYAZ koloni verir. Bu, jel çalıştırmadan önceki ilk "
            "elemedir ve yalnızca lacZ-dM15 taşıyan konakta çalışır."),
        steps=[
            RecipeStep("EcoRI", "BamHI",
                       ["PRO_FNRS", "RBS_B0034", "CDS_ARGA", "TER_B0015"],
                       "Hiperammonemi kasetini MCS'e tak",
                       "Klonlama bölgesinin bir kısmı ve lacZ-alfa "
                       "çerçevesi bozulur."),
        ],
    ),
    Recipe(
        key="pbr322_oksalat",
        title="pBR322 → hiperoksalüri vektörü (insersiyonel inaktivasyon)",
        vector_id="PBR322", disease="oksalat",
        why_vector=(
            "İki ayrı direnç markerı taşır ve BamHI bölgesi tetA geninin "
            "İÇİNDEDİR. Bu, mavi/beyaz taramanın olmadığı yerde ikinci "
            "bir eleme yöntemi verir."),
        source="Bolivar ve ark. 1977, Gene 2:95 (pBR322)",
        teaches=(
            "İnsersiyonel inaktivasyon. Kaset tetA'nın içine girdiği için "
            "doğru klon ampisiline DİRENÇLİ ama tetrasikline DUYARLI "
            "olur. Kolonileri iki plağa kopyalayıp bu farkı ararsınız - "
            "mavi/beyaz taramanın alternatifi budur."),
        steps=[
            RecipeStep("EcoRI", "BamHI",
                       ["PRO_FNRS", "RBS_B0034", "CDS_OXC", "RBS_B0032",
                        "CDS_FRC", "TER_B0015"],
                       "Oksalat yıkım operonunu tak",
                       "tetA bozulur: tetrasiklin duyarlılığı eleme "
                       "ölçütü olur."),
        ],
    ),
    Recipe(
        key="pseva_psa",
        title="pSEVA231 → P. aeruginosa algılayıcı devresi",
        vector_id="PSEVA231", disease="psa",
        why_vector=(
            "SEVA standardı: marker, orijin ve kargo üç ayrı modül "
            "hâlinde durur ve sınırları hep aynı enzimlerle kesilir. "
            "pBBR1 replikonu geniş konak aralıklıdır, yani devre E. coli "
            "dışındaki Gram-negatiflerde de denenebilir."),
        source="Silva-Rocha ve ark. 2013, Nucleic Acids Res 41:D666 (SEVA)",
        teaches=(
            "Modüler iskelet ve iki kasetli devre. Kargo bölgesinin iki "
            "yanındaki T0/T1 terminatörleri, devrenin iskeletle "
            "konuşmasını keser - algılayıcı devrelerde bu şarttır, yoksa "
            "iskeletten gelen okuma arka plan sinyali yaratır."),
        steps=[
            RecipeStep("EcoRI", "SpeI",
                       ["PRO_J23105", "RBS_B0032", "CDS_LASR", "TER_B0015",
                        "IGEM_R0062", "RBS_B0034", "CDS_PYOS5",
                        "IGEM_B0015"],
                       "Algılayıcı + öldürücü devreyi kargo bölgesine tak",
                       "İki kaset tek seferde girer."),
        ],
    ),
    Recipe(
        key="pnz_ibd",
        title="pNZ8048 → IBD vektörü (Gram-pozitif iskelet)",
        vector_id="PNZ8048", disease="ibd",
        why_vector=(
            "Lactococcus lactis replikonu taşır - Hastalık vektörü "
            "sayfasındaki IBD tasarımının eksik parçası tam olarak "
            "budur. E. coli plazmiti L. lactis'te çoğalmaz; iskeleti "
            "buradan almak o açığı kapatır."),
        source=("MK297322 (pNZ8048 türevi). DİKKAT: bu kayıt boş vektör "
                "değildir, içinde kargo vardır."),
        teaches=(
            "Konak organizmasının iskeleti belirlemesi. Aynı terapötik "
            "kaset, hangi bakteride çalışacaksa o bakterinin replikonuna "
            "takılmalıdır. Denetim, E. coli orijinli bir yapıyı "
            "L. lactis konağıyla eşleştirirseniz ORG-01 ile uyarır."),
        steps=[
            RecipeStep("EcoRI", "HindIII",
                       ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034",
                        "SIG_USP45", "CDS_IL10", "TER_B0015"],
                       "IL-10 salgılama kasetini tak",
                       "Usp45 sinyali proteini besiyerine çıkarır."),
        ],
    ),
]


def recipe_by_key(key: str) -> "Recipe | None":
    return next((r for r in RECIPES if r.key == key), None)
