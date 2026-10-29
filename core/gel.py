"""Agaroz jel simülasyonu: sindirim sonucunu jelde nasıl göreceksin.

Fragman uzunluklarının listesi tek başına yanıltıcıdır. Jelde iki şey
olur ve ikisi de tabloda görünmez:

  - Birbirine yakın boyuttaki fragmanlar tek bant halinde çıkar. 2000 ve
    2100 bp'yi ayırt edemezsiniz; tabloda iki satır, jelde bir bant.
  - Bant parlaklığı molekül sayısıyla değil, o bandaki DNA KÜTLESİYLE
    orantılıdır. 100 bp'lik bir fragman, eşit molde olmasına rağmen
    3000 bp'lik fragmanın otuzda biri kadar boya bağlar ve çoğu zaman
    hiç görünmez.

Göç modeli: agarozda göç mesafesi kabaca log(boyut) ile doğrusaldır.
Jel yüzdesi, hangi boyut aralığının ayrıştığını belirler.
"""

from __future__ import annotations

import html
import math
from dataclasses import dataclass, field

# Jel yüzdesi -> (iyi ayrışan en büyük bp, en küçük bp).
# Bu aralığın dışındaki fragmanlar ya kuyunun dibinde sıkışır ya da
# boya cephesiyle birlikte jelden çıkar.
AGAROSE_RANGE: dict[float, tuple[int, int]] = {
    0.7: (12000, 800),
    1.0: (8000, 400),
    1.2: (6000, 300),
    1.5: (3000, 200),
    2.0: (1500, 100),
    3.0: (800, 50),
}

# Yaygın merdivenler. Üreticiler arasında küçük farklar olabilir;
# bunlar standart boyut setleridir.
LADDERS: dict[str, list[int]] = {
    "1 kb": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500],
    "1 kb plus": [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000,
                  700, 500, 400, 300, 200, 100],
    "100 bp": [1500, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
}

# Merdivende belirgin kalınlıkta basılan referans bantlar.
LADDER_BOLD: dict[str, set[int]] = {
    "1 kb": {3000, 1000},
    "1 kb plus": {1000, 500},
    "100 bp": {500},
}


@dataclass
class Band:
    size: int
    intensity: float            # 0-1 arası görünürlük
    label: str = ""
    merged_from: list[int] = field(default_factory=list)
    off_scale: str = ""         # "ust" | "alt" | ""


@dataclass
class Lane:
    name: str
    bands: list[Band]
    is_ladder: bool = False
    ladder_key: str = ""


# 50 bp'nin altındaki fragmanlar hiçbir agaroz jelde güvenilir biçimde
# görünmez; jel yüzdesini onlara göre seçmek, asıl bakılacak bantları
# ayrışmaz hale getirir.
VISIBLE_FLOOR = 50


def recommend_agarose(sizes: list[int]) -> float:
    """Fragman boyutlarına göre uygun jel yüzdesini seç.

    Ölçüt: en çok fragmanı ayırma aralığının içine alan yüzde. Eşitlik
    durumunda jelin aralığı, fragmanların boyut dağılımına en iyi oturan
    yüzde seçilir - yani jel veriye ortalanır.
    """
    real = [s for s in sizes if s and s >= VISIBLE_FLOOR]
    if not real:
        return 1.0

    hedef_log = sum(math.log10(s) for s in real) / len(real)

    def score(pct: float) -> tuple[int, float]:
        hi, lo = AGAROSE_RANGE[pct]
        icinde = sum(1 for s in real if lo <= s <= hi)
        merkez = (math.log10(hi) + math.log10(lo)) / 2
        return (icinde, -abs(hedef_log - merkez))

    return max(AGAROSE_RANGE, key=score)


def build_lane(
    name: str, sizes: list[int], agarose: float = 1.0,
    total_ng: float = 500.0, merge_tolerance: float = 0.06,
) -> Lane:
    """Fragman listesinden jel kulvarı kur: birleştir, parlaklık hesapla."""
    real = sorted((int(s) for s in sizes if s and s > 0), reverse=True)
    if not real:
        return Lane(name=name, bands=[])

    hi, lo = AGAROSE_RANGE.get(agarose, (8000, 400))
    toplam = sum(real)

    # Yakın boyutları birleştir - jelde zaten tek bant görünürler.
    gruplar: list[list[int]] = []
    for s in real:
        if gruplar and abs(gruplar[-1][0] - s) <= gruplar[-1][0] * merge_tolerance:
            gruplar[-1].append(s)
        else:
            gruplar.append([s])

    bands: list[Band] = []
    for g in gruplar:
        kutle_ng = total_ng * sum(g) / toplam
        # Görünürlük: ~5 ng altı çoğu jelde seçilmez, ~50 ng doyuma yakın.
        yogunluk = max(0.0, min(1.0, math.log10(max(kutle_ng, 0.5) / 2.0) / 1.6))
        temsil = g[0]
        off = ""
        if temsil > hi:
            off = "ust"
        elif temsil < lo:
            off = "alt"
        bands.append(Band(
            size=temsil,
            intensity=round(yogunluk, 3),
            label=f"{temsil} bp" if len(g) == 1 else f"{temsil} bp ×{len(g)}",
            merged_from=g if len(g) > 1 else [],
            off_scale=off,
        ))
    return Lane(name=name, bands=bands)


def ladder_lane(key: str = "1 kb") -> Lane:
    sizes = LADDERS.get(key, LADDERS["1 kb"])
    bold = LADDER_BOLD.get(key, set())
    return Lane(
        name=f"{key} merdiven",
        bands=[Band(size=s, intensity=1.0 if s in bold else 0.72,
                    label=f"{s}") for s in sizes],
        is_ladder=True, ladder_key=key,
    )


def lane_notes(lane: Lane, agarose: float) -> list[str]:
    """Kulvardaki sorunları düz Türkçe anlat."""
    out: list[str] = []
    hi, lo = AGAROSE_RANGE.get(agarose, (8000, 400))
    for b in lane.bands:
        if b.merged_from:
            out.append(
                f"{', '.join(str(x) + ' bp' for x in b.merged_from)} "
                f"fragmanları tek bant halinde çıkar - aralarındaki fark "
                f"%{agarose} jelde ayrışmayacak kadar küçük."
            )
        if b.off_scale == "ust":
            out.append(
                f"{b.size} bp, %{agarose} jelin ayırma sınırının ({hi} bp) "
                f"üstünde: kuyunun hemen altında sıkışır, boyutu "
                f"okunamaz. Daha düşük yüzdeli jel kullanın."
            )
        if b.off_scale == "alt":
            out.append(
                f"{b.size} bp, %{agarose} jelin alt sınırının ({lo} bp) "
                f"altında: boya cephesiyle birlikte jelden çıkabilir. "
                f"Daha yüksek yüzdeli jel kullanın."
            )
        if b.intensity < 0.2:
            out.append(
                f"{b.size} bp bandı çok sönük kalır: bu fragman toplam "
                f"DNA'nın küçük bir kısmı olduğu için az boya bağlar. "
                f"Kuyuya daha fazla DNA yükleyin ya da bu bandı "
                f"doğrulama ölçütü saymayın."
            )
    return out


# ---------------------------------------------------------------------------
# Çizim
# ---------------------------------------------------------------------------

_GEL_BG = "#17171C"
_BAND = "#E8F0D8"
_WELL = "#0C0C10"
_TEXT = "#C9C9D4"


def _y_for(size: int, agarose: float, top: float, bottom: float) -> float:
    hi, lo = AGAROSE_RANGE.get(agarose, (8000, 400))
    # Aralığın biraz dışına da yer bırak ki sınır bantları üst üste binmesin.
    hi_x, lo_x = hi * 1.8, lo / 1.8
    t = (math.log10(hi_x) - math.log10(max(size, 1))) / (
        math.log10(hi_x) - math.log10(lo_x))
    t = max(0.0, min(1.0, t))
    return top + t * (bottom - top)


def render_gel(
    lanes: list[Lane], agarose: float = 1.0, width: int = 760,
    title: str = "",
) -> str:
    """Jeli SVG olarak çiz."""
    lanes = [ln for ln in lanes if ln.bands]
    if not lanes:
        return "<svg xmlns='http://www.w3.org/2000/svg'></svg>"

    pad_l, pad_r, pad_t = 18, 18, 46
    gel_h = 330
    lane_w = max(52, (width - pad_l - pad_r) // max(len(lanes), 1))
    gel_w = lane_w * len(lanes)
    total_w = pad_l + gel_w + pad_r
    gel_top = pad_t
    gel_bottom = pad_t + gel_h
    well_y = gel_top + 10
    run_top = well_y + 20
    run_bottom = gel_bottom - 20

    p: list[str] = [
        f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {total_w} "
        f"{gel_bottom + 54}' width='100%' role='img' "
        f"aria-label='Agaroz jel simülasyonu'>",
        "<style>"
        ".gl{font:11px ui-sans-serif,system-ui,sans-serif;fill:" + _TEXT + "}"
        ".gn{font:600 12px ui-sans-serif,system-ui,sans-serif;fill:" + _TEXT + "}"
        ".gt{font:600 13px ui-sans-serif,system-ui,sans-serif;fill:currentColor}"
        "</style>",
    ]
    if title:
        p.append(f"<text x='{pad_l}' y='20' class='gt'>{html.escape(title)}"
                 f"</text>")
    p.append(f"<text x='{pad_l}' y='36' class='gl'>%{agarose} agaroz · "
             f"göç yukarıdan aşağıya</text>")
    p.append(f"<rect x='{pad_l}' y='{gel_top}' width='{gel_w}' height='{gel_h}'"
             f" rx='4' fill='{_GEL_BG}'/>")

    for i, ln in enumerate(lanes):
        lx = pad_l + i * lane_w
        cx = lx + lane_w / 2
        # Kuyu
        p.append(f"<rect x='{cx - lane_w * 0.31:.1f}' y='{well_y - 7}' "
                 f"width='{lane_w * 0.62:.1f}' height='9' rx='2' "
                 f"fill='{_WELL}'/>")
        p.append(f"<text x='{cx:.1f}' y='{gel_bottom + 18}' class='gn' "
                 f"text-anchor='middle'>{html.escape(ln.name[:16])}</text>")

        for b in ln.bands:
            y = _y_for(b.size, agarose, run_top, run_bottom)
            op = max(0.12, b.intensity)
            bw = lane_w * 0.6
            # Büyük fragmanlar jelde biraz daha yayılır.
            bh = 5.0 if b.size < 2000 else 6.0
            blur = 0.6 if b.intensity > 0.5 else 1.2
            p.append(
                f"<rect x='{cx - bw / 2:.1f}' y='{y - bh / 2:.1f}' "
                f"width='{bw:.1f}' height='{bh:.1f}' rx='{bh / 2:.1f}' "
                f"fill='{_BAND}' opacity='{op:.2f}' "
                f"style='filter:blur({blur}px)'>"
                f"<title>{html.escape(b.label)}</title></rect>"
            )
            if ln.is_ladder:
                p.append(f"<text x='{lx + 4}' y='{y + 3.5:.1f}' class='gl' "
                         f"font-size='9'>{html.escape(b.label)}</text>")
            elif b.intensity >= 0.18:
                p.append(f"<text x='{cx + bw / 2 + 4:.1f}' y='{y + 3.5:.1f}' "
                         f"class='gl' font-size='9'>"
                         f"{html.escape(b.label)}</text>")

    p.append(f"<text x='{pad_l}' y='{gel_bottom + 42}' class='gl'>"
             f"Bant parlaklığı o bantdaki DNA kütlesiyle orantılı çizildi; "
             f"kısa fragmanlar eşit molde bile sönük görünür.</text>")
    p.append("</svg>")
    return "".join(p)
