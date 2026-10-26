"""Restriksiyon enzimi kütüphanesi.

Veri kaynağı REBASE'dir (Roberts ve ark., rebase.neb.com) ve Biopython
ile birlikte gelir - akademik kullanıma açıktır. Biopython 1088 enzim
taşır, bunların 623'ü ticari olarak satılmaktadır.

Modülün asıl işi, "bu enzim nereden keser" sorusunu gösterilebilir hale
getirmek: tanıma dizisindeki kesim geometrisi ve gerçek bir dizideki
kesim bağlamı, iki zincir birden çizilerek verilir.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

from Bio.Restriction import AllEnzymes, CommOnly, RestrictionBatch
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Metilasyon duyarlılığı (REBASE'de bu alan yok, küratörlü liste)
# ---------------------------------------------------------------------------

#: E. coli'den saflaştırılan plazmit Dam (GATC) ve Dcm (CCWGG) metilasyonu
#: taşır. Aşağıdaki enzimler bu yüzden beklenen bölgeyi kesmeyebilir.
METHYLATION_BLOCKED: dict[str, str] = {
    "XbaI": "Dam - TCTAGA bölgesi 5'-GATC ile komşuysa bloke olur",
    "ClaI": "Dam - ATCGAT bölgesi bağlama göre bloke olabilir",
    "MboI": "Dam - GATC her zaman metillenir, E. coli DNA'sını kesmez",
    "BclI": "Dam - TGATCA her zaman bloke",
    "BspHI": "Dam - bağlama göre bloke olabilir",
    "MluI": "Dam - bağlama göre etkilenebilir",
    "NruI": "Dam - bağlama göre etkilenebilir",
    "TaqI": "Dam - bağlama göre etkilenebilir",
    "XhoII": "Dam - bağlama göre bloke",
    "StuI": "Dcm - bağlama göre bloke olabilir",
    "AvaII": "Dcm - etkilenebilir",
    "SexAI": "Dcm - etkilenebilir",
    "EcoRII": "Dcm - CCWGG her zaman bloke",
    "ApaI": "Dcm - bağlama göre etkilenebilir",
    "PpuMI": "Dcm - bağlama göre etkilenebilir",
}

#: Dam/Dcm eksik suşlar - metilasyon sorunu yaşandığında plazmit burada
#: çoğaltılır.
DAM_DCM_FREE_STRAINS = ["JM110", "SCS110", "INV110", "GM2163"]

#: Golden Gate ve benzeri yöntemlerde kullanılan Tip IIS enzimleri:
#: tanıma dizisinin DIŞINDAN keserler, bu yüzden istediğiniz çıkıntıyı
#: tasarlayabilirsiniz.
TYPE_IIS = [
    "BsaI", "BsmBI", "BbsI", "SapI", "AarI", "BtgZI", "BspQI", "Esp3I",
    "BsmAI", "BciVI", "BfuAI", "BspMI", "BtsI", "EarI", "HgaI", "MlyI",
    "PaqCI", "SfaNI",
]

#: Ders ve rutin klonlamada en sık kullanılanlar - varsayılan listeler
#: bunlarla başlar, ama artık tüm kütüphane erişilebilir.
WORKHORSE = [
    "EcoRI", "EcoRV", "BamHI", "HindIII", "XhoI", "NotI", "NcoI", "NdeI",
    "XbaI", "SpeI", "PstI", "SalI", "SacI", "KpnI", "SmaI", "XmaI", "ApaI",
    "BglII", "ClaI", "HpaI", "MluI", "NheI", "PvuI", "PvuII", "ScaI", "SphI",
    "StuI", "AatII", "AflII", "AgeI", "AscI", "AvrII", "DraI", "FseI", "PacI",
    "PmeI", "SgrAI", "SnaBI", "SwaI", "BsaI", "BsmBI", "BbsI", "SapI",
]


# ---------------------------------------------------------------------------
# Kütüphane erişimi
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def _enzyme_index() -> dict[str, Any]:
    return {str(e): e for e in AllEnzymes}


@lru_cache(maxsize=1)
def _commercial_names() -> frozenset[str]:
    return frozenset(str(e) for e in CommOnly)


def get(name: str):
    """Enzim nesnesini ada göre getir; bulunamazsa None."""
    return _enzyme_index().get(name)


def library_size() -> dict[str, int]:
    return {"toplam": len(_enzyme_index()), "ticari": len(_commercial_names())}


# ---------------------------------------------------------------------------
# Kesim geometrisi
# ---------------------------------------------------------------------------


@dataclass
class CutGeometry:
    """Enzimin tanıma dizisini nereden kestiğinin sayısal tarifi."""

    name: str
    site: str
    top_cut: int        # üst zincirde kesim, tanıma dizisinin başından itibaren
    bottom_cut: int     # alt zincirde kesim (aynı koordinat sisteminde)
    overhang: int       # <0: 5' çıkıntı, >0: 3' çıkıntı, 0: kör
    overhang_seq: str
    outside: bool       # Tip IIS gibi tanıma dizisinin dışından mı kesiyor


def geometry(name: str) -> CutGeometry | None:
    e = get(name)
    # Tip IIS enzimleri Biopython'da "Ambiguous" sınıfındadır: kesim
    # konumları bilinir, ama çıkıntı DİZİSİ tanıma bölgesinde değil
    # çevredeki DNA'da olduğu için NNNN görünür. Golden Gate'in çalışma
    # mantığı budur, dolayısıyla bunları dışlamak yanlış olur.
    # Sadece kesim konumu hiç bilinmeyen enzimler (is_unknown) elenir.
    if e is None or e.is_unknown():
        return None
    top = int(e.fst5)
    ovhg = int(e.ovhg)
    # Biopython'da fst5 üst zincirdeki kesim noktasıdır; alt zincirdeki
    # kesim ondan `ovhg` kadar ötededir. Çıkıntı negatifse (5' çıkıntı)
    # alt zincir daha ileriden kesilir.
    bottom = top - ovhg
    return CutGeometry(
        name=name,
        site=str(e.site),
        top_cut=top,
        bottom_cut=bottom,
        overhang=ovhg,
        overhang_seq=str(getattr(e, "ovhgseq", "") or ""),
        outside=bool(top < 0 or top > len(str(e.site))),
    )


def overhang_label(ovhg: int) -> str:
    if ovhg == 0:
        return "kör uç"
    if ovhg < 0:
        return f"5' çıkıntı ({abs(ovhg)} nt)"
    return f"3' çıkıntı ({ovhg} nt)"


def _complement(s: str) -> str:
    return str(Seq(s).complement())


def _spaced(seq: str) -> str:
    """Bazları birer boşlukla ayır: kesim okları aralara hizalanabilsin."""
    return " ".join(seq)


def _marker_line(length: int, cut: int, char: str, prefix: int = 3) -> str:
    """`cut` kadar bazdan SONRAKİ boşluğa ok koy.

    `_spaced` ile her baz 2 karakter kaplar; i. bazdan sonraki boşluk
    2*i - 1 karakterindedir. prefix, "5'-" önekinin genişliğidir.
    """
    if cut <= 0:
        col = prefix - 1
    elif cut >= length:
        col = prefix + 2 * length - 1
    else:
        col = prefix + 2 * cut - 1
    return " " * max(0, col) + char


def duplex_view(top: str, top_cut: int, bottom_cut: int) -> str:
    """Çift zinciri hizalı çiz, iki kesimi ok ile göster.

    Bazlar dikeyde hizalı kalır - komplementerlik görünür olmalı.
    ▼ üst zincirdeki, ▲ alt zincirdeki kesimi gösterir.
    """
    comp = _complement(top)
    n = len(top)
    return "\n".join([
        _marker_line(n, top_cut, "▼"),
        f"5'-{_spaced(top)}-3'",
        f"3'-{_spaced(comp)}-5'",
        _marker_line(n, bottom_cut, "▲"),
    ])


def product_view(top: str, top_cut: int, bottom_cut: int) -> str:
    """Kesimden SONRA oluşan iki parçayı yapışkan uçlarıyla göster."""
    comp = _complement(top)
    a, b = sorted((top_cut, bottom_cut))
    # Sol parça
    sol_ust, sol_alt = top[:top_cut], comp[:bottom_cut]
    sag_ust, sag_alt = top[top_cut:], comp[bottom_cut:]
    w = max(len(sol_ust), len(sol_alt))
    return "\n".join([
        f"5'-{sol_ust.ljust(w)}      {sag_ust}-3'",
        f"3'-{sol_alt.ljust(w)}      {sag_alt}-5'",
    ])


def site_diagram(name: str) -> str:
    """Tanıma dizisi üzerinde kesimi iki zincirde, hizalı göster."""
    g = geometry(name)
    if g is None:
        return ""
    top, bottom = g.top_cut, g.bottom_cut
    # Tanıma dizisinin dışından kesen enzimlerde (Tip IIS) araya N dolgusu
    pad = max(0, max(top, bottom) - len(g.site))
    shown = g.site + "N" * pad
    return duplex_view(shown, top, bottom)


def site_products(name: str) -> str:
    """Tanıma dizisi kesildiğinde oluşan uçlar."""
    g = geometry(name)
    if g is None:
        return ""
    pad = max(0, max(g.top_cut, g.bottom_cut) - len(g.site))
    shown = g.site + "N" * pad
    return product_view(shown, g.top_cut, g.bottom_cut)


# ---------------------------------------------------------------------------
# Gerçek bir dizideki kesim bağlamı
# ---------------------------------------------------------------------------


def cut_context(
    seq: str, name: str, position: int, flank: int = 12
) -> dict[str, str] | None:
    """Dizideki belirli bir kesimi iki zincirde, çevresiyle birlikte göster.

    `position`, Biopython'ın `search()` çıktısıdır: üst zincirde kesimden
    SONRAKİ ilk bazın 1 tabanlı konumu.
    """
    g = geometry(name)
    if g is None or not seq:
        return None

    L = len(seq)
    top_cut0 = position - 1              # 0 tabanlı kesim noktası
    bottom_cut0 = top_cut0 - g.overhang

    lo = max(0, min(top_cut0, bottom_cut0) - flank)
    hi = min(L, max(top_cut0, bottom_cut0) + flank)
    window = seq[lo:hi]

    t_rel = top_cut0 - lo
    b_rel = bottom_cut0 - lo

    # Çıkıntı dizisi: iki kesim arasında kalan bölge. Tip IIS'te bu
    # bölge çevredeki DNA'dan gelir, tanıma dizisinden değil.
    a, b = sorted((top_cut0, bottom_cut0))
    ovh = seq[a:b]

    return {
        "duplex": duplex_view(window, t_rel, b_rel),
        "products": product_view(window, t_rel, b_rel),
        "window_start": lo + 1,
        "window_end": hi,
        "truncated_left": lo > 0,
        "truncated_right": hi < L,
        "overhang_seq": ovh,
        "overhang_label": overhang_label(g.overhang),
        "top_cut_pos": top_cut0 + 1,
        "bottom_cut_pos": bottom_cut0 + 1,
    }


# ---------------------------------------------------------------------------
# Ayrıntı kartı
# ---------------------------------------------------------------------------


def details(name: str) -> dict[str, Any]:
    e = get(name)
    if e is None:
        return {}
    g = geometry(name)
    try:
        iso = sorted(str(x) for x in e.isoschizomers())
    except Exception:
        iso = []
    try:
        freq = float(e.frequency())
    except Exception:
        freq = None
    return {
        "name": name,
        "site": str(e.site),
        "size": int(e.size),
        "elucidate": e.elucidate() if not e.is_unknown() else "",
        "diagram": site_diagram(name),
        "overhang": overhang_label(g.overhang) if g else "",
        "overhang_seq": g.overhang_seq if g else "",
        "blunt": bool(e.is_blunt()),
        "type_iis": name in TYPE_IIS or (g.outside if g else False),
        "isoschizomers": iso,
        "commercial": name in _commercial_names(),
        "suppliers": len(getattr(e, "suppl", ()) or ()),
        "opt_temp": getattr(e, "opt_temp", None),
        "frequency": freq,
        "methylation": METHYLATION_BLOCKED.get(name, ""),
        "ambiguous": bool(e.is_ambiguous()),
    }


# ---------------------------------------------------------------------------
# Arama / süzme
# ---------------------------------------------------------------------------


def search(
    query: str = "",
    commercial_only: bool = True,
    site_length: int | None = None,
    overhang: str | None = None,   # "5", "3", "kor"
    type_iis_only: bool = False,
    limit: int | None = None,
) -> list[str]:
    """Kütüphanede enzim ara.

    `limit` varsayılan olarak YOKTUR. Daha önce 400'lük bir varsayılan
    vardı ve alfabetik sıralama yüzünden N'den sonraki enzimler (NdeI,
    NotI, PstI, XhoI…) listelerden sessizce düşüyordu.
    """
    names = sorted(_commercial_names() if commercial_only else _enzyme_index())
    q = query.strip().lower()
    out = []
    for n in names:
        e = _enzyme_index()[n]
        if q and q not in n.lower() and q not in str(e.site).lower():
            continue
        if site_length and int(e.size) != site_length:
            continue
        if type_iis_only and n not in TYPE_IIS:
            continue
        if overhang:
            if overhang == "kor" and not e.is_blunt():
                continue
            if overhang == "5" and not e.is_5overhang():
                continue
            if overhang == "3" and not e.is_3overhang():
                continue
        out.append(n)
        if limit is not None and len(out) >= limit:
            break
    return out


def compatible_with(name: str, commercial_only: bool = True) -> list[str]:
    """Aynı uca ligate edilebilen diğer enzimler.

    Yönlü klonlamada bunlardan kaçınmak, uyumlu uç arayan montajlarda ise
    bunları bulmak gerekir.
    """
    e = get(name)
    if e is None:
        return []
    out = []
    pool = _commercial_names() if commercial_only else set(_enzyme_index())
    for other in sorted(pool):
        if other == name:
            continue
        o = _enzyme_index()[other]
        try:
            if e.is_blunt() and o.is_blunt():
                out.append(other)
            elif (
                getattr(e, "ovhgseq", None)
                and getattr(o, "ovhgseq", None)
                and e.ovhgseq == o.ovhgseq
                and e.ovhg == o.ovhg
            ):
                out.append(other)
        except Exception:
            continue
    return out


def sites_in(seq: str, names: list[str], circular: bool = True) -> dict[str, list[int]]:
    """Verilen enzimlerin dizideki kesim konumları."""
    if not seq or not names:
        return {}
    try:
        batch = RestrictionBatch(names)
    except (ValueError, KeyError):
        return {}
    bio = Seq(seq)
    out: dict[str, list[int]] = {}
    for enz in batch:
        try:
            out[str(enz)] = sorted(enz.search(bio, linear=not circular))
        except (ValueError, TypeError):
            continue
    return out
