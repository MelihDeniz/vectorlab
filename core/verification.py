"""Doğrulama: dizileme planı ve okunan diziyi tasarımla karşılaştırma.

Klonlamanın son adımı dizilemedir ve iki yerde tökezlenir.

Birincisi kapsama. Sanger okuması ~800-1000 bp'dir ve primerin hemen
ardındaki ilk 30-50 baz güvenilir okunmaz. Tek primerle uzun bir
inserti okutup "doğruladım" sanmak sık yapılan hatadır; aradaki bölge
hiç görülmemiştir. Burası, verilen uzunluğu baştan sona kapatmak için
kaç primer gerektiğini ve nereye konacaklarını hesaplar.

İkincisi karşılaştırma. Okunan dizi ile tasarım arasındaki tek bir baz
farkı haritada da jelde de görünmez ama proteini bozabilir. Hizalama
Biopython'un `PairwiseAligner` modülüyle yapılır ve farklar tek tek
listelenir; kodlayan bölgedeyse amino asit karşılığı da çıkarılır.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from Bio import Align
from Bio.Data import CodonTable

from core.seq_utils import (
    clean_seq, gc_percent, primer_dimer, primer_hairpin,
    revcomp, tm_nn,
)

_T11 = CodonTable.unambiguous_dna_by_id[11]

# Sanger okumasının gerçekçi sınırları.
READ_LENGTH = 900          # güvenilir okunan toplam uzunluk
DEAD_ZONE = 50             # primerin hemen ardındaki okunamayan bölge
OVERLAP = 150              # komşu okumaların güvenli örtüşmesi


@dataclass
class SeqPrimer:
    name: str
    sequence: str
    position: int          # kalıptaki başlangıç (1 tabanlı)
    strand: int
    tm: float | None
    gc: float
    covers_from: int
    covers_to: int
    warnings: list[str] = field(default_factory=list)


@dataclass
class CoveragePlan:
    length: int
    primers: list[SeqPrimer] = field(default_factory=list)
    gaps: list[tuple[int, int]] = field(default_factory=list)
    notes: list[dict[str, str]] = field(default_factory=list)

    @property
    def covered(self) -> int:
        return self.length - sum(b - a + 1 for a, b in self.gaps)


@dataclass
class Difference:
    position: int          # tasarımdaki konum (1 tabanlı)
    kind: str              # "degisim" | "eksik" | "fazla"
    expected: str
    observed: str
    codon_no: int | None = None
    aa_change: str = ""


@dataclass
class AlignReport:
    identity: float
    aligned_from: int
    aligned_to: int
    read_length: int
    differences: list[Difference] = field(default_factory=list)
    reversed_read: bool = False
    notes: list[dict[str, str]] = field(default_factory=list)


def _note(kind: str, title: str, detail: str) -> dict[str, str]:
    return {"kind": kind, "title": title, "detail": detail}


# ---------------------------------------------------------------------------
# Dizileme planı
# ---------------------------------------------------------------------------


def _pick_primer(seq: str, around: int, strand: int = 1,
                 target_tm: float = 58.0) -> tuple[str, int] | None:
    """`around` konumu civarında iyi bir primer bul (1 tabanlı döner)."""
    n = len(seq)
    en_iyi, en_iyi_puan = None, -1e9
    for kaydir in range(-40, 41, 2):
        bas = around + kaydir
        if bas < 1 or bas + 24 > n:
            continue
        for uzunluk in (18, 20, 22, 24):
            if bas + uzunluk - 1 > n:
                continue
            aday = seq[bas - 1:bas - 1 + uzunluk]
            if strand < 0:
                aday = revcomp(aday)
            tm = tm_nn(aday)
            if tm is None:
                continue
            gc = gc_percent(aday)
            puan = -abs(tm - target_tm) * 2 - abs(gc - 50) * 0.3
            if aday[-1] in "GC":
                puan += 2
            for b in "ACGT":
                if b * 5 in aday:
                    puan -= 10
            # Kalıpta birden çok yere bağlanıyorsa işe yaramaz
            if seq.count(aday[-15:]) > 1:
                puan -= 25
            if puan > en_iyi_puan:
                en_iyi, en_iyi_puan = (aday, bas), puan
    return en_iyi


def sequencing_plan(
    seq: str, read_length: int = READ_LENGTH,
    dead_zone: int = DEAD_ZONE, overlap: int = OVERLAP,
    both_strands: bool = False,
) -> CoveragePlan:
    """Diziyi baştan sona okutmak için primer yürüyüşü planla."""
    s = clean_seq(seq)
    plan = CoveragePlan(length=len(s))
    if len(s) < 60:
        plan.notes.append(_note("uyari", "Dizi çok kısa",
                                "Dizileme planı için en az 60 nt gerekir."))
        return plan

    yararli = read_length - dead_zone          # primer sonrası okunan
    adim = max(yararli - overlap, 200)
    konumlar: list[int] = []
    k = 1
    while k < len(s):
        konumlar.append(k)
        if k + yararli >= len(s):
            break
        k += adim

    for i, konum in enumerate(konumlar, start=1):
        bulunan = _pick_primer(s, konum)
        if not bulunan:
            continue
        dizi, bas = bulunan
        kapsam_bas = min(bas + len(dizi) + dead_zone, len(s))
        kapsam_son = min(bas + len(dizi) + read_length, len(s))
        tm = tm_nn(dizi)
        uyari = []
        if tm is not None and not (54 <= tm <= 66):
            uyari.append(f"Tm {tm:.1f} °C - 55-65 aralığı dışında")
        if s.count(dizi[-15:]) > 1:
            uyari.append("kalıpta birden çok yere bağlanıyor")
        toka = primer_hairpin(dizi)
        if toka["risk"] != "yok":
            uyari.append(f"saç tokası (puan {toka['score']})")
        kendi = primer_dimer(dizi)
        if kendi["risk"] != "yok":
            uyari.append(f"kendi dimeri (puan {kendi['score']})")
        plan.primers.append(SeqPrimer(
            name=f"seq_F{i}", sequence=dizi, position=bas, strand=1,
            tm=tm, gc=gc_percent(dizi),
            covers_from=kapsam_bas, covers_to=kapsam_son, warnings=uyari,
        ))

    if both_strands:
        for i, konum in enumerate(reversed(konumlar), start=1):
            bulunan = _pick_primer(s, min(konum + yararli, len(s) - 25),
                                   strand=-1)
            if not bulunan:
                continue
            dizi, bas = bulunan
            plan.primers.append(SeqPrimer(
                name=f"seq_R{i}", sequence=dizi, position=bas, strand=-1,
                tm=tm_nn(dizi), gc=gc_percent(dizi),
                covers_from=max(1, bas - read_length),
                covers_to=max(1, bas - dead_zone), warnings=[],
            ))

    # Kapsanmayan bölgeler
    kapali = [False] * (len(s) + 1)
    for p in plan.primers:
        for x in range(p.covers_from, min(p.covers_to, len(s)) + 1):
            kapali[x] = True
    bosluk_bas = None
    for x in range(1, len(s) + 1):
        if not kapali[x] and bosluk_bas is None:
            bosluk_bas = x
        elif kapali[x] and bosluk_bas is not None:
            plan.gaps.append((bosluk_bas, x - 1))
            bosluk_bas = None
    if bosluk_bas is not None:
        plan.gaps.append((bosluk_bas, len(s)))

    oran = round(100 * plan.covered / max(len(s), 1), 1)
    plan.notes.append(_note(
        "bilgi" if oran > 99 else "uyari", "Kapsama",
        f"{len(plan.primers)} primerle dizinin %{oran}'i okunuyor "
        f"({plan.covered}/{len(s)} nt). Okuma başına {read_length} nt "
        f"varsayıldı ve primerin hemen ardındaki {dead_zone} nt "
        f"güvenilmez sayıldı."))
    if plan.gaps:
        plan.notes.append(_note(
            "uyari", "Okunmayan bölgeler",
            ", ".join(f"{a}–{b}" for a, b in plan.gaps[:6])
            + ". Buralara ayrıca primer gerekir; yoksa o bölgedeki bir "
              "mutasyon fark edilmez."))
    plan.notes.append(_note(
        "bilgi", "İki zincir",
        "Kritik bölgeleri (insertin tamamı, birleşme noktaları) iki "
        "zincirden de okutun. Tek zincirli okumada kromatogramın "
        "belirsiz kaldığı noktalar, karşı zincirdeki okumayla "
        "kesinleşir."))
    return plan


# ---------------------------------------------------------------------------
# Okunan diziyi karşılaştırma
# ---------------------------------------------------------------------------


def _aligner() -> Align.PairwiseAligner:
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.match_score = 2
    a.mismatch_score = -1
    a.open_gap_score = -6
    a.extend_gap_score = -1
    return a


def compare_read(
    design: str, read: str, cds_start: int | None = None,
    trim: int = 0,
) -> AlignReport:
    """Dizileme okumasını tasarımla hizala ve farkları çıkar.

    `cds_start` verilirse (1 tabanlı), kodlayan bölgedeki farkların
    amino asit karşılığı da hesaplanır - sessiz mi, değiştiren mi,
    yoksa erken durdurma mı.
    """
    d = clean_seq(design)
    r = clean_seq(read)
    if trim:
        r = r[trim:]
    rep = AlignReport(identity=0.0, aligned_from=0, aligned_to=0,
                      read_length=len(r))
    if len(d) < 20 or len(r) < 20:
        rep.notes.append(_note("uyari", "Dizi çok kısa",
                               "Karşılaştırma için en az 20 nt gerekir."))
        return rep

    al = _aligner()
    duz = al.align(d, r)
    ters = al.align(d, revcomp(r))
    if len(ters) and len(duz) and ters[0].score > duz[0].score:
        r = revcomp(r)
        best = ters[0]
        rep.reversed_read = True
    elif len(duz):
        best = duz[0]
    else:
        rep.notes.append(_note("hata", "Hizalanamadı",
                               "Okuma tasarımla eşleşmiyor."))
        return rep

    hedef_bloklar, sorgu_bloklar = best.aligned
    if len(hedef_bloklar) == 0:
        rep.notes.append(_note("hata", "Hizalanamadı", ""))
        return rep

    rep.aligned_from = int(hedef_bloklar[0][0]) + 1
    rep.aligned_to = int(hedef_bloklar[-1][1])
    # Bloklarin kendileri de numpy; donguye girmeden sade int'e cevir.
    hedef_bloklar = [(int(a), int(b)) for a, b in hedef_bloklar]
    sorgu_bloklar = [(int(a), int(b)) for a, b in sorgu_bloklar]

    eslesme = 0
    toplam = 0
    farklar: list[Difference] = []
    onceki_hedef_son = None
    onceki_sorgu_son = None
    for (hb, hs), (sb, ss) in zip(hedef_bloklar, sorgu_bloklar):
        # Bloklar arası boşluk = indel
        if onceki_hedef_son is not None:
            if hb > onceki_hedef_son:
                farklar.append(Difference(
                    position=int(onceki_hedef_son) + 1, kind="eksik",
                    expected=d[onceki_hedef_son:hb], observed="-"))
                toplam += hb - onceki_hedef_son
            if sb > onceki_sorgu_son:
                farklar.append(Difference(
                    position=int(onceki_hedef_son) + 1, kind="fazla",
                    expected="-", observed=r[onceki_sorgu_son:sb]))
                toplam += sb - onceki_sorgu_son
        for k in range(hs - hb):
            toplam += 1
            if d[hb + k] == r[sb + k]:
                eslesme += 1
            else:
                farklar.append(Difference(
                    position=int(hb) + k + 1, kind="degisim",
                    expected=d[hb + k], observed=r[sb + k]))
        onceki_hedef_son, onceki_sorgu_son = int(hs), int(ss)

    rep.identity = round(100.0 * eslesme / max(toplam, 1), 2)

    # Kodlayan bölgedeki farkların amino asit karşılığı
    if cds_start:
        for f in farklar:
            if f.kind != "degisim" or f.position < cds_start:
                continue
            ofset = f.position - cds_start
            kodon_no = ofset // 3
            icinde = ofset % 3
            kb = cds_start - 1 + kodon_no * 3
            eski_kodon = d[kb:kb + 3]
            if len(eski_kodon) < 3:
                continue
            yeni_kodon = (eski_kodon[:icinde] + f.observed
                          + eski_kodon[icinde + 1:])
            eski_aa = ("*" if eski_kodon in _T11.stop_codons
                       else _T11.forward_table.get(eski_kodon, "?"))
            yeni_aa = ("*" if yeni_kodon in _T11.stop_codons
                       else _T11.forward_table.get(yeni_kodon, "?"))
            f.codon_no = kodon_no + 1
            f.aa_change = (f"{eski_aa}{kodon_no + 1}{yeni_aa}"
                           if eski_aa != yeni_aa
                           else f"sessiz ({eski_aa}{kodon_no + 1})")

    rep.differences = farklar
    rep.notes = _interpret_read(rep, d)
    return rep


def _interpret_read(rep: AlignReport, design: str) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    if rep.reversed_read:
        out.append(_note(
            "bilgi", "Okuma ters zincirden",
            "Okuma tasarımın ters tümleyeniyle eşleşti; karşılaştırma "
            "çevrilerek yapıldı. Geri primerle okuttuysanız bu beklenen "
            "durumdur."))

    kapsam = rep.aligned_to - rep.aligned_from + 1
    out.append(_note(
        "bilgi", "Hizalama",
        f"Okuma tasarımın {rep.aligned_from}–{rep.aligned_to}. "
        f"nükleotidlerine oturdu ({kapsam} nt), benzerlik %{rep.identity}. "
        f"Tasarımın geri kalanı bu okumayla doğrulanmadı."))

    if not rep.differences:
        out.append(_note(
            "bilgi", "Fark yok",
            "Okunan bölgede tasarımla fark bulunamadı."))
        return out

    degisim = [f for f in rep.differences if f.kind == "degisim"]
    indel = [f for f in rep.differences if f.kind != "degisim"]

    anlamli = [f for f in degisim if f.aa_change
               and not f.aa_change.startswith("sessiz")]
    durdurma = [f for f in anlamli if f.aa_change.endswith("*")]

    if durdurma:
        out.append(_note(
            "hata", "Erken durdurma kodonu",
            ", ".join(f.aa_change for f in durdurma)
            + ". Protein bu noktada kesilir; bu klon kullanılamaz."))
    if anlamli and not durdurma:
        out.append(_note(
            "uyari", "Amino asit değişimi",
            ", ".join(f.aa_change for f in anlamli[:8])
            + ". PCR sırasında girmiş olabilir. Başka bir koloni okutun; "
              "aynı değişiklik birden çok kolonide çıkıyorsa kalıpta "
              "vardır."))
    if indel:
        out.append(_note(
            "hata", "Eksik ya da fazla baz",
            f"{len(indel)} yerde uzunluk farkı var. Kodlayan bölgede üçün "
            f"katı olmayan bir ekleme/silme okuma çerçevesini kaydırır ve "
            f"proteini tamamen bozar."))
    if degisim and not anlamli and not indel:
        out.append(_note(
            "bilgi", "Yalnızca sessiz değişimler",
            f"{len(degisim)} baz farkı var ama hiçbiri amino asidi "
            f"değiştirmiyor. Protein tasarlandığı gibi."))
    return out


# ---------------------------------------------------------------------------
# Kromatogram (.ab1) okuma
# ---------------------------------------------------------------------------
# Dizileme firması size iki şey gönderir: metin dizi ve kromatogram.
# Metin dizi, makinenin her konumda verdiği en iyi tahmindir ve kalite
# bilgisini taşımaz. Oysa okumanın başı ve sonu neredeyse her zaman
# güvenilmezdir; oradaki bir "fark" gerçek mutasyon değil, okuma
# hatasıdır. Kalite değerleri olmadan ikisini ayırt edemezsiniz.


@dataclass
class Chromatogram:
    sequence: str
    qualities: list[int]
    length: int
    trim_start: int           # 1 tabanlı, güvenilir bölgenin başı
    trim_end: int
    mean_quality: float
    trimmed_sequence: str
    notes: list[dict[str, str]] = field(default_factory=list)


def _mott_trim(qual: list[int], cutoff: int = 20) -> tuple[int, int]:
    """Düşük kaliteli uçları kırp (Mott'un kayan toplam yöntemi).

    Her baz için `kalite - cutoff` puanı toplanır: eşiğin üstündeki
    bazlar toplamı büyütür, altındakiler küçültür. Toplam negatife
    düşerse sıfırlanır ve başlangıç ileri kaydırılır; en yüksek toplamı
    veren aralık güvenilir bölgedir. Sabit eşikle kesmekten daha
    iyidir, çünkü ortadaki tek bir kötü bazın okumayı ikiye bölmesini
    engeller.
    """
    if not qual:
        return (0, 0)
    en_iyi_toplam = toplam = 0.0
    en_iyi = (0, 0)
    bas = 0
    for i, q in enumerate(qual):
        toplam += (q - cutoff)
        if toplam < 0:
            toplam = 0.0
            bas = i + 1
        elif toplam > en_iyi_toplam:
            en_iyi_toplam = toplam
            en_iyi = (bas, i + 1)
    return en_iyi


def read_chromatogram(data: bytes, quality_cutoff: int = 20) -> Chromatogram | None:
    """ABI (.ab1) dosyasından diziyi ve kalite değerlerini oku."""
    from io import BytesIO

    from Bio import SeqIO

    try:
        rec = SeqIO.read(BytesIO(data), "abi")
    except Exception:
        return None

    dizi = str(rec.seq).upper()
    qual = list(rec.letter_annotations.get("phred_quality", []))
    if not qual:
        qual = [0] * len(dizi)

    bas, son = _mott_trim(qual, quality_cutoff)
    kirpik = dizi[bas:son]
    ortalama = round(sum(qual) / len(qual), 1) if qual else 0.0

    c = Chromatogram(
        sequence=dizi, qualities=qual, length=len(dizi),
        trim_start=bas + 1, trim_end=son, mean_quality=ortalama,
        trimmed_sequence=kirpik,
    )

    oran = round(100 * len(kirpik) / max(len(dizi), 1), 1)
    c.notes.append(_note(
        "bilgi" if oran > 60 else "uyari", "Güvenilir bölge",
        f"Okuma {len(dizi)} baz; kalite eşiği Q{quality_cutoff} ile "
        f"güvenilir bölge {bas + 1}–{son} ({len(kirpik)} baz, "
        f"%{oran}). Ortalama kalite Q{ortalama}.\n\n"
        "Karşılaştırmayı bu bölgeyle yapın: okumanın başındaki ve "
        "sonundaki farklar neredeyse her zaman okuma hatasıdır, "
        "gerçek mutasyon değil."))

    if ortalama < 20:
        c.notes.append(_note(
            "uyari", "Okuma kalitesi düşük",
            f"Ortalama Q{ortalama}. Q20 bir bazda %1 hata demektir; "
            f"bunun altı dizilemenin başarısız olduğunu gösterir. "
            f"Kalıp derişimi, primer kalitesi ya da plazmitteki saç "
            f"tokaları sebep olabilir - aynı bölgeyi karşı zincirden "
            f"okutun."))
    if len(kirpik) < 200:
        c.notes.append(_note(
            "uyari", "Güvenilir bölge çok kısa",
            "200 bazdan kısa güvenilir okuma, yapıyı doğrulamak için "
            "yeterli değildir. Dizilemeyi tekrarlayın."))
    return c
