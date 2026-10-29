"""mRNA ikincil yapısı: başlama bölgesini kapatan saç tokaları.

Neden önemli: ribozomun küçük alt birimi mRNA'ya Shine-Dalgarno
dizisinden tutunur ve başlama kodonuna kadar kayar. Bu bölge kendi
üzerine katlanıp çift zincirli bir saç tokası (hairpin) yaparsa ribozom
bağlanamaz. Promotör güçlü, kodonlar optimize, RBS doğru olduğu halde
protein alınamayan durumların önemli bir kısmı buradan çıkar - ve bu
plazmit haritasında hiç görünmez.

Hesap nasıl yapılıyor: sap (stem) bölgesinin serbest enerjisi
Biopython'un RNA en yakın komşu tablosundan (`RNA_NN3`, Xia ve ark.
1998) gelir, ilmek (loop) için ise standart biçimde yaklaşık bir ceza
eklenir.

ÖNEMLİ SINIR: bu tam bir katlanma hesabı DEĞİLDİR. Gerçek bir katlanma
programı (ViennaRNA, mfold) bütün olası yapıları birlikte değerlendirip
en düşük enerjili olanı bulur; burada yalnızca tek tek saç tokaları
taranır. Çıkan sayı mutlak bir ΔG değeri olarak değil, "bu bölgede
güçlü bir yapı var mı" sorusunun göreli yanıtı olarak okunmalıdır.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from Bio.SeqUtils import MeltingTemp as _mt

from core.seq_utils import clean_seq

T_KELVIN = 310.15          # 37 °C
RT = 0.0019872 * T_KELVIN  # kcal/mol

_NN = _mt.RNA_NN3

# Hangi bazlar eşleşir. G-U salınımı (wobble) RNA'da gerçektir ve
# yapıyı kararlı kılar, bu yüzden dahil edilir.
PAIRS = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"),
         ("G", "T"), ("T", "G")}

MIN_LOOP = 3               # 3 nükleotidden kısa ilmek fiziksel olarak olmaz

# Shine-Dalgarno çekirdeği; RBS'in ribozomla eşleşen kısmı.
SD_CORE = "AGGAGG"


@dataclass
class Hairpin:
    stem_start: int        # 1 tabanlı, sapın 5' ucu
    stem_end: int          # sapın 3' ucu (ilmeğin sonrası dahil)
    loop_start: int
    loop_end: int
    stem_len: int
    loop_len: int
    dg: float              # kcal/mol, negatif = kararlı
    sequence: str


@dataclass
class StructureReport:
    window_start: int
    window_end: int
    window_seq: str
    hairpins: list[Hairpin] = field(default_factory=list)
    strongest: Hairpin | None = None
    blocking_rbs: Hairpin | None = None
    blocking_start: Hairpin | None = None
    sd_position: int | None = None
    start_position: int | None = None
    notes: list[dict[str, str]] = field(default_factory=list)


def _pairs(a: str, b: str) -> bool:
    return (a, b) in PAIRS


def _stack_dg(top: str, bottom: str) -> float:
    """İki komşu baz çiftinin istifleme enerjisi (37 °C, kcal/mol).

    En yakın komşu tabloları her istifi tek bir kanonik biçimde tutar:
    `TC/AG` tabloda yoktur ama onunla aynı fiziksel istif olan `GA/CT`
    vardır (iki zincir ters çevrilince aynı şeye bakıyoruz). Bu yüzden
    anahtar bulunamazsa simetrik eşdeğeri denenir - bu yapılmazsa
    istiflerin yarısı sessizce kaybolur ve saç tokaları olduğundan çok
    daha kararsız görünür.
    """
    key = f"{top}/{bottom}"
    if key in _NN:
        dh, ds = _NN[key]
        return dh - T_KELVIN * ds / 1000.0
    flip = f"{bottom[::-1]}/{top[::-1]}"
    if flip in _NN:
        dh, ds = _NN[flip]
        return dh - T_KELVIN * ds / 1000.0
    # Geriye kalan yalnızca G-U içeren bazı istifler: tabloda yoklar,
    # yok saymak yerine küçük bir kararlılık katkısı veriyoruz.
    return -0.5


def _loop_penalty(n: int) -> float:
    """İlmek entropisi - yaklaşık.

    3-8 nükleotidlik ilmekler için sabit bir ceza, daha uzunlar için
    Jacobson-Stockmayer biçiminde logaritmik artış kullanılır. Gerçek
    katlanma programları ilmek dizisine göre de düzeltme yapar; burada
    o düzeltme yok.
    """
    if n < MIN_LOOP:
        return 99.0
    taban = 4.5
    if n <= 8:
        return taban
    return taban + 1.75 * RT * math.log(n / 8.0)


def find_hairpins(
    seq: str, min_stem: int = 4, max_dg: float = -3.0,
    max_loop: int = 12,
) -> list[Hairpin]:
    """Diziden saç tokalarını çıkar (en kararlıdan başlayarak)."""
    s = clean_seq(seq)
    n = len(s)
    out: list[Hairpin] = []
    gorulen: set[tuple[int, int]] = set()

    for i in range(n):
        for j in range(i + MIN_LOOP + 1, min(n, i + MIN_LOOP + 1 + max_loop)):
            if not _pairs(s[i], s[j]):
                continue
            # (i, j) ilmeği kapatan çift; dışa doğru uzat.
            uzunluk = 1
            while (i - uzunluk >= 0 and j + uzunluk < n
                   and _pairs(s[i - uzunluk], s[j + uzunluk])):
                uzunluk += 1
            if uzunluk < min_stem:
                continue
            bas, son = i - uzunluk + 1, j + uzunluk - 1
            if (bas, son) in gorulen:
                continue
            gorulen.add((bas, son))

            dg = 0.0
            for k in range(uzunluk - 1):
                p, q = i - k, j + k
                dg += _stack_dg(s[p - 1] + s[p], s[q + 1] + s[q])
            dg += _loop_penalty(j - i - 1)
            if dg > max_dg:
                continue
            out.append(Hairpin(
                stem_start=bas + 1, stem_end=son + 1,
                loop_start=i + 2, loop_end=j,
                stem_len=uzunluk, loop_len=j - i - 1,
                dg=round(dg, 2), sequence=s[bas:son + 1],
            ))
    return sorted(out, key=lambda h: h.dg)


def _overlaps(h: Hairpin, start: int, end: int) -> bool:
    return not (h.stem_end < start or h.stem_start > end)


def _sd_score(seq: str, start0: int) -> int:
    """Başlama kodonunun 4-20 nt öncesinde AGGAGG'a en iyi benzerlik."""
    bolge = seq[max(0, start0 - 20):max(0, start0 - 4)]
    en_iyi = 0
    for k in range(max(0, len(bolge) - 5)):
        en_iyi = max(en_iyi, sum(1 for a, b in zip(bolge[k:k + 6], SD_CORE)
                                 if a == b))
    return en_iyi


def guess_start(seq: str) -> int | None:
    """Gerçek başlama kodonunu tahmin et (1 tabanlı).

    Dizideki ilk ATG'yi almak yanlış sonuç verir: promotör dizilerinin
    içinde de ATG geçer (örneğin -10 kutusu TATAAT'ın hemen ardından
    gelen G ile). Bunun yerine her aday puanlanır - önünde
    Shine-Dalgarno benzeri bir dizi olması ve arkasında uzun bir açık
    okuma çerçevesi açması aranır.
    """
    s = clean_seq(seq)
    n = len(s)
    adaylar: list[tuple[float, int]] = []
    for i in range(n - 2):
        if s[i:i + 3] not in ("ATG", "GTG", "TTG"):
            continue
        # Bu adaydan başlayan çerçevede ilk durdurma kodonuna kadar uzunluk
        uzunluk = 0
        for j in range(i, n - 2, 3):
            if s[j:j + 3] in ("TAA", "TAG", "TGA"):
                break
            uzunluk += 1
        sd = _sd_score(s, i)
        # SD varlığı belirleyici, ORF uzunluğu ikincil ölçüt.
        puan = (3.0 if sd >= 5 else 1.0 if sd == 4 else 0.0) * 100 + uzunluk
        if s[i:i + 3] == "ATG":
            puan += 5          # ATG, GTG/TTG'ye göre daha olası
        adaylar.append((puan, i + 1))
    if not adaylar:
        return None
    return max(adaylar)[1]


def analyze_start_region(
    seq: str, start_codon_pos: int | None = None,
    upstream: int = 45, downstream: int = 35,
) -> StructureReport:
    """Başlama kodonunun çevresindeki yapıyı incele.

    `start_codon_pos` verilmezse ilk ATG/GTG/TTG aranır. Pencere
    varsayılan olarak başlama kodonunun 45 nt öncesinden 35 nt
    sonrasına kadardır: ribozomun kapladığı alan kabaca budur.
    """
    s = clean_seq(seq)
    if not s:
        return StructureReport(0, 0, "")

    if start_codon_pos is None:
        start_codon_pos = guess_start(s)
    if start_codon_pos is None:
        rep = StructureReport(1, len(s), s)
        rep.notes.append({
            "kind": "uyari", "title": "Başlama kodonu bulunamadı",
            "detail": "Dizide ATG, GTG ya da TTG yok. Doğru zinciri ve "
                      "doğru bölgeyi verdiğinizden emin olun."})
        return rep

    bas = max(1, start_codon_pos - upstream)
    son = min(len(s), start_codon_pos + 2 + downstream)
    pencere = s[bas - 1:son]

    rep = StructureReport(window_start=bas, window_end=son,
                          window_seq=pencere,
                          start_position=start_codon_pos)
    rep.hairpins = find_hairpins(pencere)
    if rep.hairpins:
        rep.strongest = rep.hairpins[0]

    # Shine-Dalgarno yalnızca başlama kodonunun 4-20 nükleotid öncesinde
    # aranır. Daha geniş aranırsa AGGAGG'a 4/6 benzeyen rastgele diziler
    # her yerde eşleşir ve yanlış konum bildirilir; bu aralık ise
    # ribozomun gerçekten tutunduğu aralıktır.
    start_rel0 = start_codon_pos - bas          # 0 tabanlı
    ara_bas = max(0, start_rel0 - 20)
    ara_son = max(0, start_rel0 - 4)
    bolge = pencere[ara_bas:ara_son]
    sd_rel = -1
    if bolge:
        tam = bolge.find(SD_CORE)
        if tam >= 0:
            sd_rel = ara_bas + tam
        else:
            en_iyi, en_iyi_skor = -1, 0
            for k in range(max(0, len(bolge) - 5)):
                skor = sum(1 for a, b in zip(bolge[k:k + 6], SD_CORE)
                           if a == b)
                if skor > en_iyi_skor:
                    en_iyi, en_iyi_skor = k, skor
            if en_iyi_skor >= 5:
                sd_rel = ara_bas + en_iyi
    if sd_rel >= 0:
        rep.sd_position = bas + sd_rel
        for h in rep.hairpins:
            if _overlaps(h, sd_rel + 1, sd_rel + 6):
                rep.blocking_rbs = h
                break

    start_rel = start_codon_pos - bas
    for h in rep.hairpins:
        if _overlaps(h, start_rel + 1, start_rel + 3):
            rep.blocking_start = h
            break

    rep.notes = _interpret(rep)
    return rep


def _interpret(rep: StructureReport) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []

    if not rep.hairpins:
        out.append({
            "kind": "bilgi", "title": "Başlama bölgesi açık görünüyor",
            "detail": "Bu pencerede kayda değer bir saç tokası bulunamadı. "
                      "Ribozomun bağlanmasını engelleyen bir yapı yok "
                      "demektir - çeviri başlangıcı tarafında sorun "
                      "beklenmiyor."})
        return out

    en = rep.strongest
    out.append({
        "kind": "bilgi", "title": "En kararlı yapı",
        "detail": (
            f"{en.stem_len} baz çiftlik sap, {en.loop_len} nükleotidlik "
            f"ilmek, yaklaşık {en.dg} kcal/mol. Pencerede "
            f"{rep.window_start + en.stem_start - 1}–"
            f"{rep.window_start + en.stem_end - 1}. nükleotidler arasında. "
            f"Bu sayı yaklaşıktır; gerçek bir katlanma programı farklı "
            f"bir değer verebilir.")})

    if rep.blocking_rbs:
        h = rep.blocking_rbs
        out.append({
            "kind": "hata" if h.dg < -8 else "uyari",
            "title": "Saç tokası ribozom bağlanma bölgesini kapatıyor",
            "detail": (
                f"Yaklaşık {h.dg} kcal/mol'lük bir yapı Shine-Dalgarno "
                f"dizisinin üstüne geliyor. Ribozom oraya tutunamazsa "
                f"çeviri hiç başlamaz; promotör ne kadar güçlü olursa "
                f"olsun protein alamazsınız. Çözüm, RBS ile başlama "
                f"kodonu arasındaki diziyi sessiz kodon değişiklikleriyle "
                f"ya da farklı bir RBS seçerek açmaktır.")})

    if rep.blocking_start:
        h = rep.blocking_start
        out.append({
            "kind": "hata" if h.dg < -8 else "uyari",
            "title": "Saç tokası başlama kodonunu kapatıyor",
            "detail": (
                f"Yaklaşık {h.dg} kcal/mol'lük bir yapı başlama kodonunun "
                f"üstüne geliyor. Genin ilk 10-15 kodonunu eşanlamlı "
                f"kodonlarla değiştirmek bu yapıyı çoğu zaman açar ve "
                f"proteini değiştirmez - kodon optimizasyonunun en çok "
                f"işe yaradığı yer burasıdır.")})

    if not rep.blocking_rbs and not rep.blocking_start and en.dg < -10:
        out.append({
            "kind": "bilgi", "title": "Güçlü yapı var ama kritik yerde değil",
            "detail": "Pencerede kararlı bir saç tokası var, ancak ne "
                      "Shine-Dalgarno dizisinin ne de başlama kodonunun "
                      "üstüne geliyor. Çeviri başlangıcını engellemesi "
                      "beklenmez."})

    if rep.sd_position is None:
        out.append({
            "kind": "uyari", "title": "Shine-Dalgarno dizisi bulunamadı",
            "detail": "Pencerede AGGAGG'a benzeyen bir dizi yok. Başlama "
                      "kodonundan 5-13 nükleotid önce bir ribozom "
                      "bağlanma bölgesi olmalı; yoksa çeviri verimi çok "
                      "düşük kalır."})
    return out


def dg_scale_note() -> str:
    return ("Negatif değer kararlı yapı demektir. Kabaca: -3 kcal/mol'ün "
            "üstü önemsiz, -3 ile -8 arası orta, -8'in altı güçlü. Bu "
            "eşikler yaklaşıktır ve yapının nerede olduğu, ne kadar "
            "kararlı olduğundan daha belirleyicidir.")
