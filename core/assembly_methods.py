"""Modern montaj yöntemleri: Golden Gate, Gibson ve bölge-yönelimli mutagenez.

Restriksiyon + ligasyon tek yöntem değil. Bu modül üç yaygın alternatifi
simüle eder ve her birinin kendi tuzaklarını kontrol eder:

Golden Gate  Tip IIS enzim tanıma dizisinin DIŞINDAN keser, yani çıkıntıyı
             siz tasarlarsınız. Birden çok parça tek tüpte, doğru sırayla
             ve iz bırakmadan birleşir. Tuzak: iki birleşme noktası aynı
             (ya da birbirinin tersi) çıkıntıyı alırsa parçalar yanlış
             sırada bağlanır ve bunu jelde fark edemezsiniz.

Gibson       Uçları örtüşen parçalar tek karışımda birleşir; enzim bölgesi
             hiç gerekmez. Tuzak: örtüşme bölgesi çok kısa ya da çok
             düşük Tm'liyse birleşmez, tekrar eden dizilere denk gelirse
             yanlış yere bağlanır.

Mutagenez    Tek baz ya da tek kodon değiştirmek için primer çifti.
             Tuzak: değişiklik primerin ucuna denk gelirse bağlanma
             zayıflar; ortada olmalı.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from core import enzymes as ez
from core.seq_utils import clean_seq, gc_percent, revcomp, tm_nn

# Golden Gate'te sık kullanılan Tip IIS enzimleri.
TYPE_IIS_COMMON = ["BsaI", "BsmBI", "BbsI", "SapI", "AarI", "PaqCI"]

# Enzimin verimli kesmesi için tanıma dizisinin 5' ucuna konan
# koruyucu bazlar. Tam ucundaki tanıma dizisi kesilmez.
SPACER = "TT"


@dataclass
class Junction:
    index: int                 # kaçıncı birleşme noktası (1 tabanlı)
    left_part: str
    right_part: str
    overhang: str
    problems: list[str] = field(default_factory=list)


@dataclass
class AssemblyPrimer:
    name: str
    sequence: str
    tm_anneal: float | None
    length: int
    tail: str
    anneal: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class AssemblyPlan:
    method: str
    parts: list[str]
    junctions: list[Junction] = field(default_factory=list)
    primers: list[AssemblyPrimer] = field(default_factory=list)
    product_length: int = 0
    notes: list[dict[str, str]] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(n["kind"] == "hata" for n in self.notes)


def _note(kind: str, title: str, detail: str) -> dict[str, str]:
    return {"kind": kind, "title": title, "detail": detail}


# ---------------------------------------------------------------------------
# Golden Gate
# ---------------------------------------------------------------------------


def _overhang_problems(oh: str) -> list[str]:
    """Tek bir çıkıntının kendi içindeki sorunları."""
    p: list[str] = []
    if oh == revcomp(oh):
        p.append("palindromik - kendi kendine bağlanabilir")
    if len(set(oh)) == 1:
        p.append("tek baz tekrarı - kayarak yanlış hizalanabilir")
    gc = gc_percent(oh)
    if gc == 0:
        p.append("hiç G/C yok - bağlanma zayıf kalır")
    elif gc == 100:
        p.append("tamamen G/C - yanlış eşleşmeye toleransı artar")
    return p


def golden_gate(
    parts: list[tuple[str, str]],
    enzyme: str = "BsaI",
    circular: bool = True,
) -> AssemblyPlan:
    """Tip IIS tabanlı tek tüp montajı planla.

    Çıkıntılar parçaların birleştiği yerdeki DOĞAL diziden alınır -
    böylece üründe iz (scar) kalmaz. `parts` sırası montaj sırasıdır.
    """
    temiz = [(ad, clean_seq(s)) for ad, s in parts if clean_seq(s)]
    plan = AssemblyPlan(method="Golden Gate", parts=[a for a, _ in temiz])
    if len(temiz) < 2:
        plan.notes.append(_note("hata", "En az iki parça gerekli",
                                "Golden Gate birden çok parçayı tek tüpte "
                                "birleştirmek için kullanılır."))
        return plan

    geo = ez.geometry(enzyme)
    if geo is None:
        plan.notes.append(_note("hata", f"{enzyme} tanınmadı",
                                "Tip IIS bir enzim seçin."))
        return plan
    if not geo.outside:
        plan.notes.append(_note(
            "hata", f"{enzyme} Tip IIS değil",
            f"{enzyme} tanıma dizisinin içinden kesiyor. Golden Gate için "
            f"dışından kesen bir enzim gerekir: "
            f"{', '.join(TYPE_IIS_COMMON)}."))
        return plan

    oh_len = abs(geo.overhang) or 4
    site = geo.site

    # --- Parçaların içinde enzim bölgesi var mı
    icte = []
    for ad, s in temiz:
        n = len(ez.sites_in(s, [enzyme], circular=False).get(enzyme, []))
        if n:
            icte.append(f"{ad} ({n})")
    if icte:
        plan.notes.append(_note(
            "hata", f"{enzyme} parçaların içinde kesiyor",
            f"{', '.join(icte)}. Golden Gate'te enzim reaksiyon boyunca "
            f"etkin kalır; parçanın içindeki bölge de kesilir ve montaj "
            f"dağılır. Bu bölgeleri sessiz kodon değişiklikleriyle "
            f"temizleyin (Dizi analizi → Kodon optimizasyonu) ya da başka "
            f"bir Tip IIS enzim seçin."))

    # --- Birleşme noktaları ve çıkıntılar
    n = len(temiz)
    sinir = n if circular else n - 1
    for i in range(sinir):
        sol_ad, sol = temiz[i]
        sag_ad, sag = temiz[(i + 1) % n]
        # Çıkıntı, sağ parçanın ilk oh_len bazından alınır: ürün dizisi
        # değişmez, yalnızca kesim orada olur.
        oh = sag[:oh_len].upper()
        j = Junction(index=i + 1, left_part=sol_ad, right_part=sag_ad,
                     overhang=oh)
        j.problems = _overhang_problems(oh)
        plan.junctions.append(j)

    # --- Çıkıntılar birbirinden ayırt edilebiliyor mu
    gorulen: dict[str, int] = {}
    for j in plan.junctions:
        for onceki_oh, onceki_i in gorulen.items():
            if j.overhang == onceki_oh:
                j.problems.append(
                    f"{onceki_i}. birleşme noktasıyla AYNI çıkıntı")
            elif j.overhang == revcomp(onceki_oh):
                j.problems.append(
                    f"{onceki_i}. birleşme noktasının TERSİ - parçalar ters "
                    f"yönde bağlanabilir")
        gorulen[j.overhang] = j.index

    cakisan = [j for j in plan.junctions if any("birleşme" in p
                                                for p in j.problems)]
    if cakisan:
        plan.notes.append(_note(
            "hata", "Çıkıntılar birbirinden ayırt edilemiyor",
            "Aynı ya da birbirinin tersi çıkıntı taşıyan birleşme "
            "noktaları var: "
            + ", ".join(f"{j.index}. ({j.overhang})" for j in cakisan)
            + ". Bu durumda parçalar yanlış sırada birleşir ve ürün "
              "jelde doğru boyutta göründüğü için hata dizilemeye kadar "
              "fark edilmez. Parça sınırlarını birkaç baz kaydırarak "
              "farklı çıkıntılar elde edin."))

    zayif = [j for j in plan.junctions if j.problems and j not in cakisan]
    if zayif:
        plan.notes.append(_note(
            "uyari", "Bazı çıkıntılar zayıf",
            "; ".join(f"{j.index}. ({j.overhang}): {', '.join(j.problems)}"
                      for j in zayif)))

    # --- Primerler
    for i, (ad, s) in enumerate(temiz):
        onceki = temiz[i - 1][1] if (circular or i > 0) else ""
        # İleri primer: bu parçanın başına kendi çıkıntısını getirir.
        f_anneal = _tune(s, 60.0)
        f_tail = SPACER + site + "A" + s[:oh_len]
        # Geri primer: bir sonraki parçanın çıkıntısını taşır.
        sonraki = temiz[(i + 1) % len(temiz)][1]
        r_anneal = _tune(revcomp(s), 60.0)
        r_tail = SPACER + site + "A" + revcomp(sonraki[:oh_len])
        for yon, tail, anneal in (("F", f_tail, f_anneal),
                                  ("R", r_tail, r_anneal)):
            tm = tm_nn(anneal)
            plan.primers.append(AssemblyPrimer(
                name=f"{ad}_{yon}", sequence=tail + anneal, tm_anneal=tm,
                length=len(tail + anneal), tail=tail, anneal=anneal,
                warnings=_primer_warn(anneal, tm),
            ))
        _ = onceki

    kisa = [ad for ad, s in temiz if len(s) < 40]
    if kisa:
        plan.notes.append(_note(
            "uyari", "Çok kısa parçalar PCR ile çoğaltılamaz",
            f"{', '.join(kisa)} 40 bp'den kısa. Bu uzunlukta primer "
            f"bağlanma bölgesi kalmaz - nitekim bu parçaların primerleri "
            f"düşük Tm uyarısı veriyor. Kısa parçaları (promotör, RBS, "
            f"kısa terminatör) PCR yerine sentetik oligo olarak sipariş "
            f"edip birbirine tavlayarak çift zincirli hale getirin; "
            f"çıkıntıları doğrudan oligonun içine yazarsınız."))

    plan.product_length = sum(len(s) for _, s in temiz)
    plan.notes.append(_note(
        "bilgi", "Reaksiyon",
        f"Bütün parçalar + {enzyme} + T4 DNA ligaz tek tüpte. Döngü: "
        f"37 °C 5 dk / 16 °C 5 dk, 25-30 tekrar; sonra 60 °C 5 dk "
        f"(enzimi kapatmak için). Kesilen ve doğru birleşen ürün artık "
        f"tanıma dizisi taşımadığı için reaksiyonun sonunda birikir - "
        f"yöntemin püf noktası budur."))
    return plan


# ---------------------------------------------------------------------------
# Gibson
# ---------------------------------------------------------------------------


def gibson(
    parts: list[tuple[str, str]],
    overlap: int = 20,
    circular: bool = True,
) -> AssemblyPlan:
    """Örtüşen uçlarla montaj planla ve örtüşmeleri denetle."""
    temiz = [(ad, clean_seq(s)) for ad, s in parts if clean_seq(s)]
    plan = AssemblyPlan(method="Gibson", parts=[a for a, _ in temiz])
    if len(temiz) < 2:
        plan.notes.append(_note("hata", "En az iki parça gerekli", ""))
        return plan

    n = len(temiz)
    sinir = n if circular else n - 1
    ortusmeler: list[str] = []
    for i in range(sinir):
        sol_ad, sol = temiz[i]
        sag_ad, sag = temiz[(i + 1) % n]
        if len(sol) < overlap or len(sag) < overlap:
            plan.notes.append(_note(
                "hata", "Parça örtüşmeden kısa",
                f"{sol_ad} ya da {sag_ad}, {overlap} bp örtüşme için "
                f"yeterince uzun değil."))
            continue
        ort = sol[-overlap:]
        j = Junction(index=i + 1, left_part=sol_ad, right_part=sag_ad,
                     overhang=ort)
        tm = tm_nn(ort)
        gc = gc_percent(ort)
        if tm is not None and tm < 48:
            j.problems.append(f"örtüşme Tm düşük ({tm:.1f} °C)")
        if gc < 30:
            j.problems.append(f"GC düşük (%{gc})")
        elif gc > 70:
            j.problems.append(f"GC yüksek (%{gc})")
        if ort in ortusmeler:
            j.problems.append("aynı örtüşme başka bir noktada da var")
        ortusmeler.append(ort)
        plan.junctions.append(j)

        # Birleşme noktası primerleri
        f_anneal = _tune(sag, 62.0)
        plan.primers.append(AssemblyPrimer(
            name=f"{sag_ad}_F", sequence=ort + f_anneal,
            tm_anneal=tm_nn(f_anneal), length=len(ort) + len(f_anneal),
            tail=ort, anneal=f_anneal,
            warnings=_primer_warn(f_anneal, tm_nn(f_anneal)),
        ))
        r_anneal = _tune(revcomp(sol), 62.0)
        plan.primers.append(AssemblyPrimer(
            name=f"{sol_ad}_R", sequence=revcomp(sag[:overlap]) + r_anneal,
            tm_anneal=tm_nn(r_anneal), length=overlap + len(r_anneal),
            tail=revcomp(sag[:overlap]), anneal=r_anneal,
            warnings=_primer_warn(r_anneal, tm_nn(r_anneal)),
        ))

    sorunlu = [j for j in plan.junctions if j.problems]
    if sorunlu:
        plan.notes.append(_note(
            "uyari", "Bazı örtüşmeler zayıf",
            "; ".join(f"{j.index}. ({j.left_part}→{j.right_part}): "
                      f"{', '.join(j.problems)}" for j in sorunlu)
            + ". Örtüşme uzunluğunu artırmak çoğu sorunu çözer; parça "
              "sınırını birkaç baz kaydırmak da GC'yi dengeler."))

    plan.product_length = sum(len(s) for _, s in temiz)
    plan.notes.append(_note(
        "bilgi", "Reaksiyon",
        "Parçaları eşit molde karıştırın, Gibson/HiFi karışımını ekleyin "
        "ve 50 °C'de 15-60 dakika tutun. Ekzonükleaz uçları yiyip tek "
        "zincirli kuyruk açar, örtüşen uçlar eşleşir, polimeraz boşluğu "
        "doldurur, ligaz kapatır. Enzim bölgesi gerekmez - bu yüzden "
        "restriksiyon haritası kalabalık yapılar için elverişlidir."))
    if overlap < 15:
        plan.notes.append(_note(
            "uyari", "Örtüşme kısa",
            f"{overlap} bp örtüşme çoğu protokol için kısadır; 20-40 bp "
            f"aralığı daha güvenlidir."))
    return plan


# ---------------------------------------------------------------------------
# Bölge-yönelimli mutagenez
# ---------------------------------------------------------------------------


def mutagenesis(
    template: str, position: int, replacement: str,
    flank: int = 15,
) -> AssemblyPlan:
    """Tek baz ya da kodon değişimi için primer çifti tasarla.

    `position` 1 tabanlıdır ve değiştirilecek ilk nükleotidi gösterir.
    `replacement` yerine konacak baz ya da kodondur.
    """
    s = clean_seq(template)
    rep = clean_seq(replacement)
    plan = AssemblyPlan(method="Bölge-yönelimli mutagenez", parts=["kalıp"])
    if not s or not rep:
        plan.notes.append(_note("hata", "Eksik girdi",
                                "Kalıp dizi ve yerine konacak baz gerekli."))
        return plan
    if position < 1 or position + len(rep) - 1 > len(s):
        plan.notes.append(_note(
            "hata", "Konum dizinin dışında",
            f"Dizi {len(s)} nt; {position}. konumdan itibaren {len(rep)} "
            f"baz değiştirilemez."))
        return plan

    i = position - 1
    eski = s[i:i + len(rep)]
    if eski == rep:
        plan.notes.append(_note(
            "uyari", "Değişiklik yok",
            f"{position}. konumda zaten {rep} var."))

    sol = s[max(0, i - flank):i]
    sag = s[i + len(rep):i + len(rep) + flank]
    if len(sol) < 8 or len(sag) < 8:
        plan.notes.append(_note(
            "uyari", "Değişiklik dizinin ucuna çok yakın",
            "Primerin iki yanında en az 10-15 eşleşen baz olmalı; "
            "burada bir taraf kısa kaldı. Kalıbı uzatın ya da uç "
            "bölgeler için başka bir yöntem kullanın."))

    fwd = sol + rep + sag
    rev = revcomp(fwd)
    tm_f = tm_nn(fwd)
    plan.primers = [
        AssemblyPrimer(name="mut_F", sequence=fwd, tm_anneal=tm_f,
                       length=len(fwd), tail="", anneal=fwd,
                       warnings=_primer_warn(fwd, tm_f)),
        AssemblyPrimer(name="mut_R", sequence=rev, tm_anneal=tm_nn(rev),
                       length=len(rev), tail="", anneal=rev,
                       warnings=_primer_warn(rev, tm_nn(rev))),
    ]
    yeni = s[:i] + rep + s[i + len(rep):]
    plan.product_length = len(yeni)
    plan.junctions = [Junction(index=1, left_part=f"{position}. konum",
                               right_part=f"{eski} → {rep}",
                               overhang=rep)]
    plan.notes.append(_note(
        "bilgi", "Değişiklik",
        f"{position}. konumdaki `{eski}` yerine `{rep}` gelecek. "
        f"Primerler değişikliği ORTADA taşıyor: iki yanda "
        f"{len(sol)} ve {len(sag)} eşleşen baz var. Değişiklik primerin "
        f"ucuna yakın olursa bağlanma zayıflar ve mutasyon girmez."))
    plan.notes.append(_note(
        "bilgi", "Protokol",
        "Kalıp plazmitle PCR yapın (iki primer de aynı noktadan ters "
        "yönde uzar), sonra ürünü DpnI ile kesin. DpnI yalnızca "
        "metillenmiş DNA'yı keser; kalıp plazmit E. coli'den geldiği "
        "için metillidir ve parçalanır, PCR'da yeni sentezlenen "
        "mutant DNA metilsizdir ve sağ kalır. Bu adımı atlarsanız "
        "kolonilerin çoğu değişmemiş kalıp olur."))
    return plan


# ---------------------------------------------------------------------------
# Yardımcılar
# ---------------------------------------------------------------------------


def _tune(seq: str, target_tm: float, lo: int = 18, hi: int = 32) -> str:
    """Hedef Tm'e en yakın uzunlukta bağlanma bölgesi seç."""
    s = clean_seq(seq)
    if len(s) <= lo:
        return s
    en_iyi, en_iyi_fark = s[:lo], 1e9
    for n in range(lo, min(hi, len(s)) + 1):
        tm = tm_nn(s[:n])
        if tm is None:
            continue
        fark = abs(tm - target_tm)
        if fark < en_iyi_fark:
            en_iyi, en_iyi_fark = s[:n], fark
    return en_iyi


def _primer_warn(core: str, tm: float | None) -> list[str]:
    w: list[str] = []
    if tm is not None and tm < 52:
        w.append(f"Tm düşük ({tm:.1f} °C)")
    if tm is not None and tm > 72:
        w.append(f"Tm yüksek ({tm:.1f} °C)")
    gc = gc_percent(core)
    if gc < 35:
        w.append(f"GC düşük (%{gc})")
    elif gc > 70:
        w.append(f"GC yüksek (%{gc})")
    if core[-1:] not in ("G", "C"):
        w.append("3' ucu G/C değil - uzamanın başlaması zayıflayabilir")
    for b in "ACGT":
        if b * 5 in core:
            w.append(f"{b*5} tekrarı var")
            break
    return w
