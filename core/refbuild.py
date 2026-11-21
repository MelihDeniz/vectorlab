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
#
# Değerler elle yazılmadı: `Bio.Restriction` üzerinden REBASE'ten
# hesaplandı (`site`, `fst5`, `ovhg`) ve `tests/test_refbuild.py`
# içindeki `test_enzim_tablosu_rebase_ile_tutarli` her açılışta aynı
# hesabı tekrarlayıp tabloyla karşılaştırıyor.
#
# **Hepsi 5' çıkıntı bırakır.** `cut_and_paste` geometrisi buna göre
# kurulmuştur, o yüzden 3' çıkıntı bırakanlar (KpnI, SacI, PstI, SphI,
# ApaI, AatII, NsiI, SbfI) ve kör uç bırakanlar (SmaI, EcoRV, FspI)
# tabloda YOK. Eksiklik değil, sınır: bu enzimler Enzimler sayfasında
# bütün REBASE kütüphanesiyle birlikte zaten görünüyor, yalnızca
# kes-yapıştır tarifleri onları kullanamıyor.
#
# Tanıma dizileri belirsiz baz TAŞIMAZ. `find_sites` düz metin araması
# yapar: IUPAC kodu taşıyan bir site (SgrAI'nin CRCCGGYG'si gibi) hiç
# eşleşmez ve enzim sessizce "hiç kesmiyor" görünür. SgrAI bu yüzden
# tabloya alınmadı; bir test de bunu şart koşuyor.
ENZYMES: dict[str, dict[str, Any]] = {
    # --- 6 bazlık tanıma dizileri ---------------------------------------
    "AclI":    {"site": "AACGTT", "cut": 2, "ovh": "5' CG"},
    "AflII":   {"site": "CTTAAG", "cut": 1, "ovh": "5' TTAA"},
    "AgeI":    {"site": "ACCGGT", "cut": 1, "ovh": "5' CCGG"},
    "AvrII":   {"site": "CCTAGG", "cut": 1, "ovh": "5' CTAG"},
    "BamHI":   {"site": "GGATCC", "cut": 1, "ovh": "5' GATC"},
    "BclI":    {"site": "TGATCA", "cut": 1, "ovh": "5' GATC",
                "dam": True},
    "BglII":   {"site": "AGATCT", "cut": 1, "ovh": "5' GATC"},
    "BspEI":   {"site": "TCCGGA", "cut": 1, "ovh": "5' CCGG"},
    "BsrGI":   {"site": "TGTACA", "cut": 1, "ovh": "5' GTAC"},
    "BssHII":  {"site": "GCGCGC", "cut": 1, "ovh": "5' CGCG"},
    "ClaI":    {"site": "ATCGAT", "cut": 2, "ovh": "5' CG",
                "dam": True},
    "EagI":    {"site": "CGGCCG", "cut": 1, "ovh": "5' GGCC"},
    "EcoRI":   {"site": "GAATTC", "cut": 1, "ovh": "5' AATT"},
    "HindIII": {"site": "AAGCTT", "cut": 1, "ovh": "5' AGCT"},
    # KasI ve NarI aynı diziyi tanır ama FARKLI yerden keser: ikisi de
    # GGCGCC görür, KasI G^GCGCC, NarI GG^CGCC. Çıkıntıları da farklı
    # olduğu için ürünleri birbirine ligate OLMAZ - izoşizomer olmak
    # uyumlu olmak demek değildir.
    "KasI":    {"site": "GGCGCC", "cut": 1, "ovh": "5' GCGC"},
    "NarI":    {"site": "GGCGCC", "cut": 2, "ovh": "5' CG"},
    "MfeI":    {"site": "CAATTG", "cut": 1, "ovh": "5' AATT"},
    "MluI":    {"site": "ACGCGT", "cut": 1, "ovh": "5' CGCG",
                "dam": True},
    "NcoI":    {"site": "CCATGG", "cut": 1, "ovh": "5' CATG"},
    "NdeI":    {"site": "CATATG", "cut": 2, "ovh": "5' TA"},
    "NheI":    {"site": "GCTAGC", "cut": 1, "ovh": "5' CTAG"},
    "PciI":    {"site": "ACATGT", "cut": 1, "ovh": "5' CATG"},
    "SalI":    {"site": "GTCGAC", "cut": 1, "ovh": "5' TCGA"},
    "SpeI":    {"site": "ACTAGT", "cut": 1, "ovh": "5' CTAG"},
    "XbaI":    {"site": "TCTAGA", "cut": 1, "ovh": "5' CTAG",
                "dam": True},
    "XhoI":    {"site": "CTCGAG", "cut": 1, "ovh": "5' TCGA"},
    # --- 8 bazlık tanıma dizileri (nadir kesiciler) ----------------------
    # 8 bazlık site ortalama 65 kb'de bir geçer; bir plazmitte çoğu
    # zaman hiç kesmezler, bu da onları insert içine bölge eklemek
    # için elverişli kılar.
    "AscI":    {"site": "GGCGCGCC", "cut": 2, "ovh": "5' CGCG"},
    "NotI":    {"site": "GCGGCCGC", "cut": 2, "ovh": "5' GGCC"},
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
    # Uçlar aynıysa insert ters de girebilir. Farklı enzim kullanmak
    # yönü garanti etmeye YETMEZ: farklı tanıma dizileri aynı çıkıntıyı
    # bırakabilir (SalI/XhoI → 5' TCGA, XbaI/SpeI → 5' CTAG,
    # BamHI/BclI → 5' GATC). Bu, yönlü klonlamada en sinsi hatadır,
    # çünkü tasarım kağıt üzerinde doğru görünür.
    if ENZYMES[left]["ovh"] == ENZYMES[right]["ovh"]:
        st.notes.append(
            f"UYARI: {left} ve {right} AYNI çıkıntıyı bırakıyor "
            f"({ENZYMES[left]['ovh']}). Enzimler farklı olsa da uçlar "
            f"uyumlu olduğu için insert ters yönde de girebilir; "
            f"klonlardan yarısı yanlış yönde olur. Yönü garanti etmek "
            f"için çıkıntıları farklı olan bir çift seçin.")
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
            RecipeStep("BamHI", "SalI", PKU_KASET,
                       "1. aşama: PKU kasetini ekle",
                       "Tetrasiklin direncinin bir bölümü çıkar. "
                       "Sol enzim HindIII değil BamHI: PKU kasetinin "
                       "İÇİNDE bir HindIII yeri var, o enzimle "
                       "kesilse insert ortadan bölünürdü."),
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
            RecipeStep("HindIII", "SpeI",
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
            RecipeStep("XhoI", "BclI",
                       ["REF_NATIVE_LACTOBACILLUS_R_200", "RBS_B0034",
                        "SIG_USP45", "CDS_IL10", "TER_B0015"],
                       "IL-10 salgılama kasetini tak",
                       "Usp45 sinyali proteini besiyerine çıkarır. "
                       "Dokuz enzimimizden yalnızca BclI ve XhoI bu "
                       "kasetin dışında kalıyor: L. reuteri EF-Tu "
                       "promotörü hem EcoRI hem HindIII yeri taşıyor, "
                       "kalanlar da vektörde tek kesmiyor. Seçenek "
                       "tek çift kaldığında vektör seçimi bir "
                       "tercih değil zorunluluk olur."),
        ],
    ),
]


# Hastalık kasetleri. `core/disease_seed.py` ile aynı parçalar; burada
# ayrıca duruyorlar çünkü bu modül veritabanına ve kayıt defterine
# bağlı değildir - tek başına sınanabilmesi bilinçli bir karardır.
HCU_KASET = ["PRO_FNRS", "RBS_B0034", "CDS_MDEA", "TER_B0015"]
KOLERA_KASET = ["PRO_ARABAD", "RBS_B0034", "CDS_CQSA", "TER_B0015"]
LAKTOZ_KASET = ["RBS_B0034", "CDS_LACZ_LD", "TER_RRNB"]
MRSA_KASET = ["PRO_J23105", "RBS_B0034", "SIG_USP45", "CDS_LSS",
              "TER_B0015"]
AI2_KASET = ["RBS_B0034", "CDS_LSRK", "IGEM_B0015"]
CDIFF_KASET = ["PRO_J23100", "RBS_B0034", "SIG_AMYE", "CDS_BLAC_BL",
               "TER_B0015"]
GLP1_KASET = ["PRO_J23105", "RBS_B0034", "SIG_DSBA", "CDS_GLP1",
              "TER_B0015"]

RECIPES += [
    Recipe(
        key="pbsii_hcu",
        title="pBluescript II KS(+) → homosistinüri (uyumlu uç tuzağı)",
        vector_id="PBSIIKS", disease="homosistinuri",
        why_vector=(
            "Klonlama bölgesinde dokuz enzimimizden yedisi tek yerden "
            "kesiyor: seçenek bolluğu bu tarifi uç uyumluluğunu "
            "göstermek için en uygun yer yapıyor. f1 orijini de "
            "taşır, yani tek zincirli DNA üretebilir - dizileme ve "
            "mutagenez için tarihsel olarak önemliydi."),
        source="Short ve ark. 1988, Nucleic Acids Res 16:7583 (pBluescript)",
        teaches=(
            "**Uçların uyumluluğu, enzimlerin farklı olmasıyla aynı şey "
            "değildir.** Bu iskelette SalI ve XhoI'in ikisi de 5' TCGA "
            "bırakır; XbaI ile SpeI'in ikisi de 5' CTAG. Böyle bir "
            "çiftle keserseniz insert ters de girer ve klonlarınızın "
            "yarısı yanlış olur - tasarım kağıt üzerinde kusursuz "
            "göründüğü için bu hata en sinsi olanıdır. Tarif bu yüzden "
            "EcoRI (5' AATT) + XhoI (5' TCGA) kullanır: çıkıntılar "
            "farklı, yön garanti. Enzim tablosundaki 'Bıraktığı uç' "
            "sütununu karşılaştırarak kendiniz görebilirsiniz."),
        steps=[
            RecipeStep("EcoRI", "XhoI", HCU_KASET,
                       "MGL kasetini MCS'e tak",
                       "lacZ-alfa'nın klonlama bölgesi bozulur: "
                       "mavi/beyaz tarama burada da çalışır."),
        ],
    ),
    Recipe(
        key="ptrc_laktoz",
        title="pTrc99A → laktoz intoleransı (enzim seçimini insert belirler)",
        vector_id="PTRC99A", disease="laktoz",
        why_vector=(
            "Kendi lacI<sup>q</sup> genini taşır, yani konakta lacI "
            "olmasını beklemez: IPTG ile açılıp kapanan ifade her "
            "suşta çalışır. trc promotörü lac ve trp promotörlerinin "
            "melezidir ve ikisinden de güçlüdür."),
        source="Amann ve ark. 1988, Gene 69:301 (pTrc99A)",
        teaches=(
            "**Enzimi vektör değil insert seçer.** Bu tarifte NcoI "
            "kullanmak isterdik: tanıma dizisi CCATGG bir ATG İÇERİR, "
            "yani pTrc99A'nın NcoI bölgesi ribozom bağlanma bölgesine "
            "doğru mesafede duran hazır bir başlangıç kodonudur ve "
            "geni oraya takmak RBS aralığı hesabını gereksiz kılar.\n\n"
            "Ama **beta-galaktosidazın kendi içinde NcoI yeri var** "
            "(1007 amino asitlik bir gende bu şaşırtıcı değil): o "
            "enzimle keserseniz inserti ortadan ikiye bölersiniz, yani "
            "inserti hazırlayamazsınız. Aynı sorun HindIII, EcoRI ve "
            "SpeI'de de var. Geriye BamHI, SalI, BclI ve XbaI kalıyor; "
            "tarif BamHI + SalI kullanır. Uygulama bu elemeyi kendisi "
            "yapar: Enzimler sayfasında kasetin dizisini kesip hangi "
            "enzimin kaç yerden kestiğini görebilirsiniz."),
        steps=[
            RecipeStep("BamHI", "SalI", LAKTOZ_KASET,
                       "Beta-galaktosidazı trc promotörünün altına tak",
                       "Promotör ve RBS vektörden gelir; kaset "
                       "yalnızca geni ve terminatörü taşır."),
        ],
    ),
    Recipe(
        key="pbad24_kolera",
        title="pBAD24 → kolera (arabinozla ayarlanabilir ifade)",
        vector_id="PBAD24", disease="kolera",
        why_vector=(
            "araC düzenleyicisini ve araBAD promotörünü birlikte "
            "taşır. araBAD'ın üstünlüğü açık/kapalı olmaması: "
            "arabinoz derişimini değiştirerek ifade düzeyini "
            "kademeli ayarlayabilirsiniz.\n\n"
            "Aynı seride p15A orijinli pBAD33 de var ve düşük kopyalı "
            "olduğu için sızıntısı daha azdır; burada pBAD24 seçildi "
            "çünkü cqsA'nın içindeki tanıma yerleri pBAD33'te "
            "kullanılabilir tek bir enzim çifti bırakmıyordu. Vektör "
            "seçiminin inserte bağlı olması bu sayfanın tekrar eden "
            "dersidir."),
        source="Guzman ve ark. 1995, J Bacteriol 177:4121 (pBAD)",
        teaches=(
            "**Ayarlanabilir ifade neden önemli:** CAI-1 sinyali az "
            "olursa patojen kalabalığı algılamaz, çok olursa konak "
            "bakteri gereksiz metabolik yük taşır. Doğru dozu "
            "deneyle bulmak gerekir ve bunun için açık/kapalı bir "
            "promotör değil, kademeli bir promotör lazımdır. "
            "araBAD'ın bilinen sınırı da burada görülür: arabinoz "
            "bağırsakta da bulunur, yani 'kapalı' durum hastada "
            "garanti değildir."),
        steps=[
            RecipeStep("EcoRI", "SalI", KOLERA_KASET,
                       "cqsA kasetini araBAD'ın arkasına tak",
                       "Kaset kendi araBAD promotörünü de taşır: "
                       "vektörünkine güvenmek yerine kopyasını "
                       "getirir, böylece kaset başka iskelete de "
                       "taşınabilir."),
        ],
    ),
    Recipe(
        key="psb1c3_ai2",
        title="pSB1C3 → AI-2 söndürme (iGEM BioBrick standardı)",
        vector_id="PSB1C3", disease="ai2",
        why_vector=(
            "iGEM'in taşıyıcı vektörü. Klonlama bölgesi rastgele "
            "değildir: RFC10 standardının ön eki (EcoRI-XbaI) ve son "
            "eki (SpeI-PstI) her parçada aynı yerde durur, bu yüzden "
            "binlerce parça birbirine aynı iki enzimle takılabilir."),
        source=("Shetty ve ark. 2008, J Biol Eng 2:5 (BioBrick "
                "standardı, RFC10)"),
        teaches=(
            "**Standart, parça sayısıyla çarpılan bir kazançtır.** "
            "Her parça aynı ön ek ve son eki taşıdığı için bir "
            "parçayı ötekinin arkasına takmak hep aynı iki enzimle "
            "yapılır; katalogdaki iGEM parçalarının tamamı bu yüzden "
            "birbirine uyar. Bedeli **skar**dır: XbaI ve SpeI uçları "
            "birleşince araya sekiz bazlık bir iz (TACTAGAG) kalır ve "
            "bu iz protein füzyonlarında çerçeveyi bozar. RFC25 gibi "
            "sonraki standartlar tam olarak bu sorunu çözmek için "
            "çıktı."),
        steps=[
            RecipeStep("EcoRI", "SpeI", AI2_KASET,
                       "lsrK parçasını BioBrick ön/son eki arasına tak",
                       "Kesim yerleri standardın kendi sınırlarıdır: "
                       "EcoRI ön ekte, SpeI son ekte."),
        ],
    ),
    Recipe(
        key="pbbr1_mrsa",
        title="pBBR1MCS → MRSA (geniş konak aralıklı replikon)",
        vector_id="PBBR1MCS", disease="mrsa",
        why_vector=(
            "pBBR1 replikonu Gram-negatif bakterilerin çoğunda "
            "çoğalır: E. coli, Pseudomonas, Rhizobium, Bordetella. "
            "ColE1 ya da p15A gibi dar konaklı orijinlerle bunu "
            "yapamazsınız - tasarımı başka bir bakteride denemek "
            "istiyorsanız iskeleti baştan seçmeniz gerekir."),
        source="Antoine & Locht 1992, Mol Microbiol 6:1785 (pBBR1)",
        teaches=(
            "**Replikonun konak aralığı tasarımın sınırını çizer.** "
            "Aynı kaset, replikon dar konaklıysa yalnızca E. coli'de "
            "denenebilir; geniş konaklıysa onlarca bakteride. "
            "Lizostafin örneğinde bu doğrudan işe yarar: enzimi "
            "hangi cilt bakterisinin salgılayacağı henüz belli "
            "değilse, iskeleti şimdiden geniş konaklı seçmek "
            "tasarımı yeniden kurmaktan kurtarır. Sınırı da açıktır: "
            "pBBR1 Gram-negatiflerde çalışır, Gram-pozitiflerde "
            "çalışmaz - L. lactis için pNZ8048 tarifine bakın."),
        steps=[
            RecipeStep("HindIII", "SalI", MRSA_KASET,
                       "Lizostafin salgı kasetini tak",
                       "Usp45 sinyali kasette taşınır, iskeletten "
                       "gelmez."),
        ],
    ),
    Recipe(
        key="pub110_cdiff",
        title="pUB110 → C. difficile riski (Gram-pozitif, B. subtilis)",
        vector_id="PUB110", disease="cdiff",
        why_vector=(
            "Staphylococcus aureus plazmitidir ve Bacillus "
            "subtilis'te de çoğalır: Gram-pozitif salgılama "
            "konaklarının klasik iskeleti. Yuvarlanan halka (rolling "
            "circle) ile çoğalır, yani kopyalanma mekanizması "
            "E. coli plazmitlerinden tamamen farklıdır. Aynı aileden "
            "pC194 de katalogda var ama dokuz enzimimizden yalnızca "
            "üçü onda tek kesiyor ve hiçbir çift bu kaseti "
            "yerleştirmeye yetmiyor."),
        source="McKenzie ve ark. 1986, Plasmid 15:93 (pUB110 dizisi)",
        teaches=(
            "**Yuvarlanan halka replikasyonunun bedeli:** bu "
            "plazmitler tek zincirli bir ara ürün üzerinden çoğalır "
            "ve büyük insertlerle kararsızlaşır, parçalanır ya da "
            "çoklu formlar oluşturur. E. coli'de 10 kb'lık bir insert "
            "sorun değilken burada olabilir. Katalogdaki Gram-pozitif "
            "orijin açığının sebebi de budur: bu plazmitlerin "
            "çoğalma başlangıcı (dso) kayıtlarda ayrı bir öznitelik "
            "olarak işaretlenmez, rep geninin içinde durur - "
            "uygulamanın parça kataloğuna sınırı uydurulmadan "
            "eklenemez ve bu açık kapatılmadan yazılıdır."),
        steps=[
            RecipeStep("XbaI", "BamHI", CDIFF_KASET,
                       "Beta-laktamaz salgı kasetini tak",
                       "AmyE sinyali B. subtilis'in kendi salgı "
                       "yolunu kullanır."),
        ],
    ),
    Recipe(
        key="pkd46_diyabet",
        title="pKD46 → diyabet (sıcaklığa duyarlı replikon)",
        vector_id="PKD46", disease="diyabet",
        why_vector=(
            "repA101ts orijini 30 °C'de çoğalır, 37 °C'de çoğalmaz. "
            "Lambda Red rekombineering için tasarlanmıştır ve asıl "
            "amacı plazmitin işi bitince KOLAYCA ATILMASIDIR: "
            "kültürü 37 °C'ye alırsınız, plazmit kendiliğinden "
            "kaybolur."),
        source=("Datsenko & Wanner 2000, PNAS 97:6640 (pKD46 ve lambda "
                "Red)"),
        teaches=(
            "**Plazmiti kaybettirmek de bir tasarım aracıdır.** "
            "Sayfanın geri kalanı plazmiti TUTMAK üzerine kurulu "
            "(ccdAB, dapA, thyA); burası tersi: sıcaklığa duyarlı "
            "replikon, istediğiniz anda bırakmanızı sağlar. Terapötik "
            "bağlamda bu doğal bir güvenlik anahtarıdır - 37 °C insan "
            "vücut sıcaklığıdır, yani böyle bir plazmit hastada "
            "kendiliğinden seyrelir.\n\n"
            "Ama aynı özellik tedaviyi de bitirir: GLP-1 üretimi "
            "plazmitle birlikte kaybolur. İçerme ile etkinlik "
            "arasındaki gerilimin en çıplak hali budur ve tasarımcı "
            "bu ikisi arasında seçim yapmak zorundadır."),
        steps=[
            RecipeStep("SalI", "BclI", GLP1_KASET,
                       "GLP-1 salgı kasetini tak",
                       "Enzim SIRASI burada belirleyici: bölge sol "
                       "enzimden sağ enzime doğru saat yönünde "
                       "alınır, bu yüzden BclI+SalI yazmak "
                       "plazmitin neredeyse tamamını atar, "
                       "SalI+BclI ise yalnızca 10 bp'yi."),
        ],
    ),
]


def recipe_by_key(key: str) -> "Recipe | None":
    return next((r for r in RECIPES if r.key == key), None)
