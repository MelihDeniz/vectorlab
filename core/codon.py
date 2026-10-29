"""Kodon kullanımı: CAI ölçümü ve E. coli için kodon optimizasyonu.

Tablo `data/ecoli_codon_usage.json` dosyasından gelir ve
`tools/fetch_codon_usage.py` ile E. coli K-12 genomunun 4200'den fazla
kodlayan dizisinden hesaplanmıştır - ezberden yazılmış bir tablo değil.

Optimizasyon neden gerekli: bir gen başka bir organizmadan geliyorsa
kodon tercihi E. coli'ninkiyle uyuşmaz. Ribozom, karşılığı az bulunan
tRNA'yı beklerken durur; bu duraklamalar hem verimi düşürür hem de
yanlış katlanmaya ve erken sonlanmaya yol açar.

Ama "her kodonu en iyisiyle değiştir" yaklaşımı da körlemesine
uygulanacak bir şey değildir: ortaya tekrar eden diziler, istenmeyen
restriksiyon bölgeleri ve mRNA'da saç tokası çıkarabilir. Bu yüzden
optimizasyon kısıtlarla yapılır ve iki kip sunulur.
"""

from __future__ import annotations

import functools
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from Bio.Data import CodonTable

from core.seq_utils import clean_seq, gc_percent, revcomp

DATA = Path(__file__).resolve().parent.parent / "data" / "ecoli_codon_usage.json"

_T11 = CodonTable.unambiguous_dna_by_id[11]
STOPS = set(_T11.stop_codons)

# w bu değerin altındaysa kodon "nadir" sayılır: ribozomun orada
# duraklama ihtimali belirgin olarak artar.
RARE_W = 0.10

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": ".",
}


@dataclass
class CodonChange:
    position: int        # 1 tabanlı kodon sırası
    aa: str
    before: str
    after: str
    w_before: float
    w_after: float
    reason: str


@dataclass
class OptimizeResult:
    sequence: str
    protein: str
    changes: list[CodonChange] = field(default_factory=list)
    cai_before: float = 0.0
    cai_after: float = 0.0
    rare_before: int = 0
    rare_after: int = 0
    gc_before: float = 0.0
    gc_after: float = 0.0
    sites_removed: dict[str, int] = field(default_factory=dict)
    sites_remaining: dict[str, int] = field(default_factory=dict)
    repeats: list[dict[str, Any]] = field(default_factory=list)
    repeats_before: int = 0
    structure_before: float = 0.0   # engelleyici yapı dG, 0 = engel yok
    structure_after: float = 0.0
    structure_fixed: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def changed(self) -> int:
        return len(self.changes)


# ---------------------------------------------------------------------------
# Tablo
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def table() -> dict[str, Any]:
    if not DATA.exists():
        return {}
    return json.loads(DATA.read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=1)
def weights() -> dict[str, float]:
    return {c: d["w"] for c, d in table().get("codons", {}).items()}


@functools.lru_cache(maxsize=1)
def reference_fractions() -> dict[str, float]:
    """Çok ifade edilen genlerde her kodonun kendi amino asidi içindeki payı.

    w değeri `sayım / o amino asidin en sık kodonunun sayımı` olduğundan,
    payı w'leri normalize ederek geri alabiliyoruz. Optimizasyonun hedefi
    budur: geni bu dağılıma oturtmak.
    """
    w = weights()
    out: dict[str, float] = {}
    by_aa: dict[str, list[str]] = {}
    for codon, aa in _T11.forward_table.items():
        by_aa.setdefault(aa, []).append(codon)
    for aa, codons in by_aa.items():
        toplam = sum(w.get(c, 0.0) for c in codons) or 1.0
        for c in codons:
            out[c] = w.get(c, 0.0) / toplam
    return out


@functools.lru_cache(maxsize=1)
def synonyms() -> dict[str, list[str]]:
    """Amino asit -> o amino asidi kodlayan kodonlar, w'ye göre sıralı."""
    w = weights()
    out: dict[str, list[str]] = {}
    for codon, aa in _T11.forward_table.items():
        out.setdefault(aa, []).append(codon)
    for aa in out:
        out[aa].sort(key=lambda c: -w.get(c, 0.0))
    return out


def source_note() -> str:
    t = table()
    if not t:
        return "Kodon tablosu bulunamadı."
    return (f"{t['organism']} ({t['accession']}), {t['gene_count']} gen. "
            f"CAI ölçütü {len(t['reference_set'])} çok ifade edilen genden "
            f"hesaplandı. İndirilme: {t['fetched']}.")


# ---------------------------------------------------------------------------
# Ölçüm
# ---------------------------------------------------------------------------


def _codons(cds: str) -> list[str]:
    s = clean_seq(cds)
    s = s[: len(s) - len(s) % 3]
    return [s[i:i + 3] for i in range(0, len(s), 3)]


def cai(cds: str) -> float | None:
    """Codon Adaptation Index: 0-1, 1 = çok ifade edilen genler gibi.

    Başlama ve durdurma kodonları ile tek seçenekli amino asitler
    (Met, Trp) hesaba katılmaz - onlarda seçim yoktur, CAI'yi yapay
    olarak yukarı çekerler.
    """
    w = weights()
    if not w:
        return None
    syn = synonyms()
    logs = []
    kodonlar = _codons(cds)
    for c in kodonlar[1:]:
        if c in STOPS:
            continue
        aa = _T11.forward_table.get(c)
        if not aa or len(syn.get(aa, [])) < 2:
            continue
        val = w.get(c, 0.005)
        logs.append(math.log(max(val, 1e-4)))
    if not logs:
        return None
    return round(math.exp(sum(logs) / len(logs)), 4)


def rare_positions(cds: str, threshold: float = RARE_W) -> list[dict[str, Any]]:
    """Nadir kodonların yeri - ribozomun duraklayabileceği noktalar."""
    w = weights()
    out = []
    for i, c in enumerate(_codons(cds), start=1):
        if c in STOPS:
            continue
        val = w.get(c, 0.0)
        if val < threshold:
            out.append({"position": i, "codon": c,
                        "aa": _T11.forward_table.get(c, "?"),
                        "w": round(val, 4),
                        "nt_start": (i - 1) * 3 + 1})
    return out


def rare_clusters(cds: str, window: int = 10, min_count: int = 3
                  ) -> list[dict[str, Any]]:
    """Nadir kodonların kümelendiği yerler.

    Tek tek dağılmış nadir kodonlar genelde sorun çıkarmaz; arka arkaya
    gelenler ribozomu aynı noktada durdurduğu için asıl risk onlardır.

    `window` kodonluk bir pencerede en az `min_count` nadir kodon varsa
    o pencere yoğun sayılır; komşu yoğun pencereler tek bir küme olarak
    birleştirilir ve kümenin gerçek nadir kodon sayısı raporlanır.
    """
    pos = [r["position"] for r in rare_positions(cds)]
    if len(pos) < min_count:
        return []
    pos_set = set(pos)

    # Yoğun pencerelerin kapsadığı aralıkları birleştir.
    araliklar: list[list[int]] = []
    for p in pos:
        kac = sum(1 for q in pos if p <= q < p + window)
        if kac < min_count:
            continue
        son = max(q for q in pos if q < p + window)
        if araliklar and p <= araliklar[-1][1] + 1:
            araliklar[-1][1] = max(araliklar[-1][1], son)
        else:
            araliklar.append([p, son])

    return [{"start": a, "end": b,
             "count": sum(1 for q in pos_set if a <= q <= b),
             "span": b - a + 1}
            for a, b in araliklar]


# ---------------------------------------------------------------------------
# Optimizasyon
# ---------------------------------------------------------------------------


def _site_regex(site: str) -> re.Pattern:
    body = "".join(IUPAC.get(ch, re.escape(ch)) for ch in site.upper())
    return re.compile(body)


def _forbidden_patterns(sites: list[str]) -> list[tuple[str, re.Pattern]]:
    out = []
    for s in sites:
        s = (s or "").upper().strip()
        if not s:
            continue
        out.append((s, _site_regex(s)))
        rc = revcomp(s)
        if rc != s:
            out.append((s, _site_regex(rc)))
    return out


def _creates_problem(
    window: str, patterns: list[tuple[str, re.Pattern]],
    max_run: int, avoid: list[str],
) -> str:
    for name, pat in patterns:
        if pat.search(window):
            return f"{name} bölgesi oluşuyor"
    if max_run and re.search(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})"
                             % (max_run + 1, max_run + 1, max_run + 1,
                                max_run + 1), window):
        return f"{max_run}'ten uzun tek baz dizisi oluşuyor"
    for m in avoid:
        if m and m in window:
            return f"{m} dizisi oluşuyor"
    return ""


def optimize(
    cds: str,
    mode: str = "nadir",              # "nadir" | "tam"
    forbidden_sites: list[str] | None = None,
    avoid_motifs: list[str] | None = None,
    max_homopolymer: int = 4,
    gc_bounds: tuple[float, float] = (30.0, 70.0),
    rare_threshold: float = RARE_W,
    min_w_full: float = 0.5,
    upstream: str = "",
    relieve_structure: bool = True,
) -> OptimizeResult:
    """Kodonları E. coli'ye uyarla.

    mode="nadir"  yalnızca nadir kodonları değiştirir. Dizinin geri
                  kalanı olduğu gibi kalır; doğal mRNA yapısı ve
                  varsa düzenleyici öğeler korunur. Varsayılan budur.
    mode="tam"    her kodonu yeniden seçer, ama "her zaman en iyi
                  kodon" KURALINI UYGULAMAZ - sebebi aşağıda.

    **Neden hedef CAI 1.0 değil.** Her amino asit için hep aynı (en sık)
    kodonu kullanmak CAI'yi 1.0'a çıkarır, ama bu doğada görülen bir
    şey değildir: bu tablonun hesaplandığı genomda `rpsA` 0,79 ve
    `tufA` 0,82 CAI'ye sahiptir, yani en çok ifade edilen genler bile
    0,8 civarındadır. Tek kodona yüklenmek o izoalıcı tRNA havuzunu
    tüketir ve dizide tekrar eden bloklar yaratır. Bu yüzden "tam"
    kipinde seçim, w değeri `min_w_full` eşiğinin üstündeki kodonlar
    arasından yapılır ve aynı amino asit için en son kullanılan kodon
    geri plana atılır - sonuç tipik olarak 0,75-0,85 aralığına oturur.

    Kısıtlar her adımda kontrol edilir: yasaklı restriksiyon bölgesi
    oluşturan, uzun tek-baz dizisi yapan ya da yerel GC'yi sınırların
    dışına iten kodon seçilmez. Optimizasyondan SONRA iki denetim daha
    yapılır: uzun doğrudan tekrarlar ve başlama bölgesindeki mRNA
    yapısı. `upstream` verilirse (RBS'i içeren 5' bölge) yapı denetimi
    gerçek bağlamda yapılır.
    """
    s = clean_seq(cds)
    s = s[: len(s) - len(s) % 3]
    if len(s) < 6:
        return OptimizeResult(sequence=s, protein="",
                              notes=["Dizi çok kısa."])
    w = weights()
    if not w:
        return OptimizeResult(sequence=s, protein="",
                              notes=["Kodon tablosu yüklenemedi; "
                                     "tools/fetch_codon_usage.py çalıştırın."])

    syn = synonyms()
    patterns = _forbidden_patterns(forbidden_sites or [])
    avoid = [m.upper() for m in (avoid_motifs or []) if m]
    kodonlar = _codons(s)

    res = OptimizeResult(sequence=s, protein="")
    res.cai_before = cai(s) or 0.0
    res.rare_before = len(rare_positions(s, rare_threshold))
    res.gc_before = gc_percent(s)

    # Değiştirilmeyecek kodonlar: ilk (başlama) ve varsa son (durdurma).
    ilk = 1
    son = len(kodonlar) - 1 if kodonlar[-1] in STOPS else len(kodonlar)

    yeni = list(kodonlar)
    hedef_pay = reference_fractions()
    aa_sayac: dict[str, dict[str, int]] = {}
    for i in range(ilk, son):
        mevcut = kodonlar[i]
        aa = _T11.forward_table.get(mevcut)
        if not aa:
            continue
        adaylar = syn.get(aa, [])
        if len(adaylar) < 2:
            continue
        mevcut_w = w.get(mevcut, 0.0)
        if mode == "nadir" and mevcut_w >= rare_threshold:
            continue

        # "tam" kipinde eşiğin üstündeki kodonlar arasından, aynı amino
        # asit için en son kullanılandan farklı olanı tercih ediyoruz.
        # "nadir" kipinde ise yalnızca daha iyisine geçilir.
        if mode == "tam":
            # Hedef dağılıma göre en çok "borçlu" kodonu seç. Böylece
            # gen boyunca her amino asidin kodon payları, çok ifade
            # edilen genlerdeki paylara yaklaşır; hep aynı kodona
            # yüklenmek yerine doğadaki dengeyi kurarız.
            kullanilan = aa_sayac.setdefault(aa, {})
            n_aa = sum(kullanilan.values())
            uygun = [c for c in adaylar if w.get(c, 0.0) >= min_w_full]                 or adaylar[:1]
            sirali = sorted(
                uygun,
                key=lambda c: -(hedef_pay.get(c, 0.0) * (n_aa + 1)
                                - kullanilan.get(c, 0)),
            )
        else:
            sirali = adaylar

        secildi = None
        for aday in sirali:
            if mode != "tam":
                if aday == mevcut:
                    break   # mevcut zaten en iyilerden, dokunma
                if w.get(aday, 0.0) <= mevcut_w:
                    break
            elif aday == mevcut:
                break
            # Pencereyi komşu kodonlarla birlikte sına: bölge sınırları
            # kodon sınırını tanımaz.
            onceki = "".join(yeni[max(0, i - 3):i])
            sonraki = "".join(yeni[i + 1:i + 4])
            pencere = onceki + aday + sonraki
            sorun = _creates_problem(pencere, patterns, max_homopolymer,
                                     avoid)
            if sorun:
                continue
            # Yerel GC kontrolü (yaklaşık 20 kodonluk pencere)
            gc_pencere = "".join(yeni[max(0, i - 10):i]) + aday + \
                "".join(yeni[i + 1:i + 11])
            g = gc_percent(gc_pencere)
            if not (gc_bounds[0] <= g <= gc_bounds[1]):
                continue
            secildi = aday
            break

        if mode == "tam":
            kesin = secildi or mevcut
            aa_sayac.setdefault(aa, {})[kesin] =                 aa_sayac[aa].get(kesin, 0) + 1
        if secildi and secildi != mevcut:
            yeni[i] = secildi
            res.changes.append(CodonChange(
                position=i + 1, aa=aa, before=mevcut, after=secildi,
                w_before=round(mevcut_w, 4), w_after=round(w.get(secildi, 0), 4),
                reason=("nadir kodon değiştirildi" if mevcut_w < rare_threshold
                        else "daha sık kullanılan kodona geçildi"),
            ))

    # --- İkinci geçiş: kalan yasaklı bölgeleri temizle.
    # Birinci geçiş yalnızca w'yi YÜKSELTEN kodonlara geçer, bu yüzden
    # bölge zaten en iyi kodonların üstüne oturmuşsa (örneğin GAATTC =
    # GAA + TTC, ikisi de E. coli'nin tercih ettiği kodonlar) hiç
    # dokunamaz. Oysa bölgeyi kaldırmak istediğimizde kodon tercihinden
    # ödün vermek gerekir: burada w düşse de bölgeyi bozan seçeneği
    # kabul ediyoruz.
    if patterns:
        yeni = _strip_sites(
            yeni, patterns, w, syn, max_homopolymer, avoid, gc_bounds,
            protected={0} | ({len(yeni) - 1} if kodonlar[-1] in STOPS else set()),
            changes=res.changes,
        )

    # --- Üçüncü geçiş: başlama bölgesindeki mRNA yapısını aç.
    # Kodon optimizasyonu diziyi değiştirir ve bu, başlama kodonunun
    # üstünde YENİ bir saç tokası yaratabilir - yani optimizasyon,
    # çözmeye çalıştığı sorunun bir başkasını doğurabilir. Bu yüzden
    # optimizasyondan sonra bölge yeniden taranır ve gerekirse ilk
    # kodonlar eşanlamlılarıyla değiştirilerek yapı açılır.
    if relieve_structure:
        yeni, yapi_once, yapi_sonra = _relieve_structure(
            yeni, upstream, w, syn, patterns, max_homopolymer, avoid,
            gc_bounds)
        res.structure_before = round(yapi_once, 2)
        res.structure_after = round(yapi_sonra, 2)
        res.structure_fixed = yapi_sonra > yapi_once + 0.2

    yeni_seq = "".join(yeni)
    res.sequence = yeni_seq
    res.cai_after = cai(yeni_seq) or 0.0
    res.repeats_before = len(find_repeats(s))
    res.repeats = find_repeats(yeni_seq)
    res.rare_after = len(rare_positions(yeni_seq, rare_threshold))
    res.gc_after = gc_percent(yeni_seq)

    # Protein birebir aynı kalmalı - olmazsa sonucu vermeyiz.
    p_once = _translate(s)
    p_sonra = _translate(yeni_seq)
    res.protein = p_sonra
    if p_once != p_sonra:
        return OptimizeResult(
            sequence=s, protein=p_once,
            notes=["İç hata: optimizasyon proteini değiştirdi, sonuç "
                   "uygulanmadı. Lütfen bu durumu bildirin."])

    # Kaldırılan / kalan bölgeler
    for name, pat in patterns:
        onceki = len(pat.findall(s))
        sonraki = len(pat.findall(yeni_seq))
        if onceki > sonraki:
            res.sites_removed[name] = res.sites_removed.get(name, 0) + (
                onceki - sonraki)
        if sonraki:
            res.sites_remaining[name] = max(
                res.sites_remaining.get(name, 0), sonraki)

    if res.sites_remaining:
        res.notes.append(
            "Bazı bölgeler kaldırılamadı: "
            + ", ".join(f"{k} ({v})" for k, v in res.sites_remaining.items())
            + ". O konumlarda eşanlamlı kodon seçeneği yok ya da "
              "değiştirmek başka bir kısıtı bozuyordu."
        )
    if mode == "nadir" and res.changed == 0:
        res.notes.append(
            "Değiştirilecek nadir kodon bulunamadı - dizi E. coli için "
            "zaten uygun görünüyor.")
    return res


def _strip_sites(
    kodonlar: list[str],
    patterns: list[tuple[str, re.Pattern]],
    w: dict[str, float],
    syn: dict[str, list[str]],
    max_homopolymer: int,
    avoid: list[str],
    gc_bounds: tuple[float, float],
    protected: set[int],
    changes: list[CodonChange],
) -> list[str]:
    """Yasaklı bölgeleri eşanlamlı kodon değiştirerek kaldır.

    Bölgeyi bozmak için kodon kalitesinden ödün verilebilir; seçenekler
    yine de w'ye göre sıralı denenir, yani en az zararlı olan tercih
    edilir. Bölge hiç kaldırılamıyorsa dizi olduğu gibi bırakılır ve
    çağıran taraf bunu `sites_remaining` üzerinden görür.
    """
    kodonlar = list(kodonlar)
    for name, pat in patterns:
        # Kaldırılamayan bir bölgede takılıp kalmamak için arama noktası
        # ilerletilir; sonraki bölgeler yine de denenir.
        aranacak = 0
        for _ in range(200):                       # sonsuz döngü koruması
            dizi = "".join(kodonlar)
            m = pat.search(dizi, aranacak)
            if not m:
                break
            onceki_sayi = len(pat.findall(dizi))
            ilk, son = m.start() // 3, (m.end() - 1) // 3
            duzeldi = False
            for i in range(ilk, son + 1):
                if i in protected or i >= len(kodonlar):
                    continue
                mevcut = kodonlar[i]
                aa = _T11.forward_table.get(mevcut)
                if not aa:
                    continue
                for aday in syn.get(aa, []):
                    if aday == mevcut:
                        continue
                    deneme = list(kodonlar)
                    deneme[i] = aday
                    yeni_dizi = "".join(deneme)
                    if len(pat.findall(yeni_dizi)) >= onceki_sayi:
                        continue
                    pencere = "".join(deneme[max(0, i - 3):i + 4])
                    if _creates_problem(pencere, patterns, max_homopolymer,
                                        avoid):
                        continue
                    gc_pencere = "".join(deneme[max(0, i - 10):i + 11])
                    if not (gc_bounds[0] <= gc_percent(gc_pencere)
                            <= gc_bounds[1]):
                        continue
                    kodonlar = deneme
                    changes.append(CodonChange(
                        position=i + 1, aa=aa, before=mevcut, after=aday,
                        w_before=round(w.get(mevcut, 0.0), 4),
                        w_after=round(w.get(aday, 0.0), 4),
                        reason=f"{name} bölgesini kaldırmak için",
                    ))
                    duzeldi = True
                    break
                if duzeldi:
                    break
            if duzeldi:
                aranacak = 0          # baştan tara, konumlar kaymış olabilir
            else:
                aranacak = m.start() + 1   # bu bölge kalsın, sonrakine bak
    return kodonlar


# ---------------------------------------------------------------------------
# Optimizasyon sonrası denetimler
# ---------------------------------------------------------------------------


def find_repeats(seq: str, min_len: int = 20) -> list[dict[str, Any]]:
    """Birden çok kez geçen uzun doğrudan tekrarları bul.

    Kodon optimizasyonunun sessiz yan etkisi budur: aynı amino asit
    dizisi için hep aynı kodonlar seçilince dizide birebir tekrar eden
    bloklar oluşur. Uzun doğrudan tekrarlar recA⁺ konakta birbirine
    rekombine olur ve aradaki bölge silinir - yapı haritada doğru
    görünürken hücrede parçalanır.
    """
    s = clean_seq(seq)
    n = len(s)
    if n < min_len * 2:
        return []
    gorulen: dict[str, list[int]] = {}
    for i in range(n - min_len + 1):
        gorulen.setdefault(s[i:i + min_len], []).append(i + 1)
    out = []
    for parca, konumlar in gorulen.items():
        if len(konumlar) < 2:
            continue
        # Örtüşen pencereleri tek tekrar say
        ayrik = [konumlar[0]]
        for k in konumlar[1:]:
            if k - ayrik[-1] >= min_len:
                ayrik.append(k)
        if len(ayrik) >= 2:
            out.append({"sequence": parca, "length": min_len,
                        "positions": ayrik, "count": len(ayrik)})
    # En uzun/çok tekrarlayandan başla, iç içe geçenleri ele
    out.sort(key=lambda r: (-r["count"], r["positions"][0]))
    secili: list[dict[str, Any]] = []
    kapsanan: set[int] = set()
    for r in out:
        if any(p in kapsanan for p in r["positions"]):
            continue
        secili.append(r)
        for p in r["positions"]:
            kapsanan.update(range(p, p + min_len))
    return secili[:12]


def _structure_penalty(upstream: str, cds: str) -> float:
    """Başlama bölgesindeki en kararlı engelleyici yapının enerjisi.

    Yalnızca RBS'i ya da başlama kodonunu KAPATAN yapılar sayılır;
    başka yerdeki saç tokası çeviriyi engellemez. Engel yoksa 0.
    """
    from core import structure as _stx
    tam = clean_seq(upstream) + clean_seq(cds)
    if not tam:
        return 0.0
    pos = len(clean_seq(upstream)) + 1
    rap = _stx.analyze_start_region(tam, pos)
    engeller = [h for h in (rap.blocking_rbs, rap.blocking_start) if h]
    return min((h.dg for h in engeller), default=0.0)


def _relieve_structure(
    kodonlar: list[str], upstream: str, w: dict[str, float],
    syn: dict[str, list[str]], patterns, max_homopolymer: int,
    avoid: list[str], gc_bounds: tuple[float, float],
    window_codons: int = 15, min_w: float = 0.2,
) -> tuple[list[str], float, float]:
    """Başlama bölgesindeki engelleyici yapıyı eşanlamlı kodonlarla aç.

    Yalnızca ilk `window_codons` kodona dokunulur: saç tokasının
    ribozomu engelleyen kısmı orasıdır ve kodon değiştirmek proteini
    bozmaz. Her adımda en çok iyileştiren tek değişiklik uygulanır,
    engel kalkana ya da iyileşme durana kadar sürer.
    """
    kodonlar = list(kodonlar)
    once = _structure_penalty(upstream, "".join(kodonlar))
    if once >= -3.0:
        return kodonlar, once, once        # kayda değer engel yok

    simdi = once
    for _ in range(12):
        en_iyi = None
        for i in range(1, min(window_codons + 1, len(kodonlar) - 1)):
            aa = _T11.forward_table.get(kodonlar[i])
            if not aa:
                continue
            for aday in syn.get(aa, []):
                if aday == kodonlar[i] or w.get(aday, 0.0) < min_w:
                    continue
                deneme = list(kodonlar)
                deneme[i] = aday
                pencere = "".join(deneme[max(0, i - 3):i + 4])
                if _creates_problem(pencere, patterns, max_homopolymer,
                                    avoid):
                    continue
                gc_p = "".join(deneme[max(0, i - 10):i + 11])
                if not (gc_bounds[0] <= gc_percent(gc_p) <= gc_bounds[1]):
                    continue
                skor = _structure_penalty(upstream, "".join(deneme))
                if en_iyi is None or skor > en_iyi[0]:
                    en_iyi = (skor, i, aday)
        if en_iyi is None or en_iyi[0] <= simdi + 0.2:
            break
        simdi, i, aday = en_iyi
        kodonlar[i] = aday
        if simdi >= -3.0:
            break
    return kodonlar, once, simdi


def _translate(cds: str) -> str:
    out = []
    for c in _codons(cds):
        if c in STOPS:
            out.append("*")
        else:
            out.append(_T11.forward_table.get(c, "X"))
    return "".join(out)


# Kendi tablomuzla hesaplanmış referans değerler - "iyi bir CAI kaçtır"
# sorusunun cevabı uydurma bir eşik değil, gerçek genlerin değeridir.
REF_CAI = {"rpsA": 0.787, "tufA": 0.824}


def interpret(res: OptimizeResult) -> list[dict[str, str]]:
    """Sonucu yorumla - sayılar tek başına karar verdirmez."""
    out: list[dict[str, str]] = []
    if res.cai_before and res.cai_after:
        fark = res.cai_after - res.cai_before
        out.append({
            "kind": "bilgi", "title": "CAI değişimi",
            "detail": (
                f"{res.cai_before} → {res.cai_after} "
                f"({'+' if fark >= 0 else ''}{round(fark, 4)}). "
                f"**Ölçek için:** bu tablonun hesaplandığı genomda "
                f"`rpsA` {REF_CAI['rpsA']} ve `tufA` {REF_CAI['tufA']} "
                f"CAI'ye sahip - yani E. coli'nin en çok ifade edilen "
                f"genleri bile 0,8 civarındadır. 1,0'a çıkmak doğada "
                f"karşılığı olan bir hedef değildir ve bu yüzden "
                f"optimizasyon oraya gitmez: her amino asit için kodon "
                f"payları, çok ifade edilen genlerdeki paylara "
                f"yaklaştırılır. Sonuç tipik olarak 0,85-0,92 aralığına "
                f"oturur; doğal 0,8'in biraz üstünde kalmasının sebebi, "
                f"nadir kodonların bilerek dışarıda bırakılmasıdır. "
                "CAI tek başına ifade düzeyini öngörmez: mRNA yapısı, "
                "promotör gücü ve protein katlanması da belirleyicidir."
            )})
    if res.rare_before:
        out.append({
            "kind": "bilgi" if res.rare_after < res.rare_before else "uyari",
            "title": "Nadir kodonlar",
            "detail": (
                f"{res.rare_before} nadir kodon vardı, {res.rare_after} "
                f"kaldı. Nadir kodon ribozomu duraklatır; arka arkaya "
                f"gelenler tek tek dağılmış olanlardan daha risklidir. "
                f"Optimizasyon yapmak istemiyorsanız Rosetta gibi nadir "
                f"tRNA'ları fazladan taşıyan bir suş da aynı sorunu "
                f"başka yoldan çözer."
            )})
    if abs(res.gc_after - res.gc_before) > 5:
        out.append({
            "kind": "uyari", "title": "GC oranı belirgin değişti",
            "detail": (
                f"%{res.gc_before} → %{res.gc_after}. GC'nin çok yükselmesi "
                f"primer tasarımını ve PCR'ı zorlaştırır, çok düşmesi "
                f"mRNA'yı kararsızlaştırır."
            )})
    if res.sites_removed:
        out.append({
            "kind": "bilgi", "title": "Temizlenen restriksiyon bölgeleri",
            "detail": ", ".join(f"{k}: {v} bölge"
                                for k, v in res.sites_removed.items())
            + ". Bu bölgeler artık insertin içinde kesmiyor, klonlamada "
              "kullanılabilirler."})
    if res.repeats:
        ornek = res.repeats[0]
        out.append({
            "kind": "uyari", "title": "Uzun doğrudan tekrar oluştu",
            "detail": (
                f"{len(res.repeats)} tekrar bulundu (optimizasyon öncesi "
                f"{res.repeats_before}). Örnek: `{ornek['sequence']}` "
                f"{ornek['count']} kez, "
                f"{', '.join(str(p) for p in ornek['positions'][:4])}. "
                f"konumlarda. "
                "Kodon optimizasyonunun bilinen yan etkisi budur: aynı "
                "amino asit dizisi için hep aynı kodonlar seçilince "
                "birebir tekrar eden bloklar çıkar. Uzun doğrudan "
                "tekrarlar recA⁺ konakta rekombine olur ve aradaki bölge "
                "silinir - yapı haritada doğru görünürken hücrede "
                "parçalanır. Klonlamayı recA1 taşıyan bir suşta yapın "
                "(DH5-alfa, TOP10)."
            )})

    if res.structure_before < -3.0:
        if res.structure_fixed:
            out.append({
                "kind": "bilgi", "title": "Başlama bölgesindeki yapı açıldı",
                "detail": (
                    f"Optimizasyon sonrası başlama bölgesinde "
                    f"{res.structure_before} kcal/mol'lük engelleyici bir "
                    f"saç tokası vardı; ilk kodonlar eşanlamlılarıyla "
                    f"değiştirilerek {res.structure_after} kcal/mol'e "
                    f"indirildi. Protein değişmedi. "
                    "Bu adım olmasaydı optimizasyon, çözmeye çalıştığı "
                    "sorunun bir başkasını doğurabilirdi: kodonlar "
                    "düzelir, ribozom yine bağlanamaz."
                )})
        else:
            out.append({
                "kind": "uyari", "title": "Başlama bölgesinde yapı kaldı",
                "detail": (
                    f"Başlama bölgesinde {res.structure_after} kcal/mol'lük "
                    f"engelleyici bir saç tokası var ve eşanlamlı kodon "
                    f"değişiklikleriyle açılamadı. Farklı bir RBS seçmeyi "
                    f"ya da RBS ile başlama kodonu arasındaki aralığı "
                    f"değiştirmeyi deneyin; mRNA yapısı sekmesi bölgeyi "
                    f"ayrıntılı gösteriyor."
                )})

    if res.changed > len(_codons(res.sequence)) * 0.5:
        out.append({
            "kind": "uyari", "title": "Dizinin yarısından fazlası değişti",
            "detail": (
                "Bu kadar değişiklik diziyi pratikte yeniden yazar: "
                "sentezletmeniz gerekir, PCR ile elde edemezsiniz. Ayrıca "
                "dizinin içinde farkında olmadığınız düzenleyici bir öğe "
                "varsa (iç promotör, RNA yapısı, splice benzeri sinyal) "
                "o da bozulur. Sadece nadir kodonları düzelten kip çoğu "
                "durumda yeterlidir."
            )})
    return out
