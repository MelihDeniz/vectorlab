"""Dairesel ve doğrusal plazmit haritası (SVG).

Metin ve iskelet `currentColor` kullanır, böylece Streamlit'in açık/koyu
temasına kendiliğinden uyar. Parça renkleri her iki temada da okunacak
şekilde seçilmiştir.
"""

from __future__ import annotations

import html
import math
from typing import Any

from core.parts_seed import CATEGORY_COLOR, CATEGORY_TR
from core.seq_utils import Feature


def render(svg: str) -> None:
    """Haritayı Streamlit'e bas.

    `st.html` bu sürümde <svg> içeriğini tamamen temizliyor (eleman kabı
    boş geliyor); `st.image` çalışır ama SVG'yi base64 görsele çevirdiği
    için `currentColor` ile gelen açık/koyu tema uyumu ve <title>
    ipuçları kaybolur. Satır içi markdown ikisini de korur.

    Haritadaki tüm metinler `html.escape` ile kaçırıldığı için parça
    adları buraya ham HTML taşıyamaz.
    """
    import streamlit as st

    st.markdown(svg, unsafe_allow_html=True)


def _pol(cx: float, cy: float, r: float, ang: float) -> tuple[float, float]:
    return cx + r * math.cos(ang), cy + r * math.sin(ang)


def _ang(pos: int, length: int) -> float:
    """1 tabanlı konumu açıya çevir; 12 yönünden başla, saat yönünde."""
    return 2 * math.pi * ((pos - 1) % length) / length - math.pi / 2


def _arc_band(
    cx: float, cy: float, r_in: float, r_out: float,
    a0: float, a1: float, arrow: int = 0,
) -> str:
    """İki yarıçap arasında kalan halka dilimi; arrow ile uç oku."""
    sweep = (a1 - a0) % (2 * math.pi)
    if sweep <= 0.0001:
        sweep = 0.0001
    large = 1 if sweep > math.pi else 0

    # Ok başı için dilimin ucundan bir pay ayır
    head = min(0.06, sweep * 0.35) if arrow else 0.0
    a_body = a0 + sweep - head if arrow > 0 else a0 + head
    if arrow > 0:
        b0, b1 = a0, a_body
    elif arrow < 0:
        b0, b1 = a_body, a0 + sweep
    else:
        b0, b1 = a0, a0 + sweep

    body_sweep = (b1 - b0) % (2 * math.pi)
    body_large = 1 if body_sweep > math.pi else 0

    x0o, y0o = _pol(cx, cy, r_out, b0)
    x1o, y1o = _pol(cx, cy, r_out, b1)
    x1i, y1i = _pol(cx, cy, r_in, b1)
    x0i, y0i = _pol(cx, cy, r_in, b0)

    d = (
        f"M {x0o:.2f} {y0o:.2f} "
        f"A {r_out:.2f} {r_out:.2f} 0 {body_large} 1 {x1o:.2f} {y1o:.2f} "
        f"L {x1i:.2f} {y1i:.2f} "
        f"A {r_in:.2f} {r_in:.2f} 0 {body_large} 0 {x0i:.2f} {y0i:.2f} Z"
    )

    if arrow:
        r_mid = (r_in + r_out) / 2
        if arrow > 0:
            tip_a, base_a = a0 + sweep, a_body
        else:
            tip_a, base_a = a0, a_body
        tx, ty = _pol(cx, cy, r_mid, tip_a)
        bx1, by1 = _pol(cx, cy, r_out, base_a)
        bx2, by2 = _pol(cx, cy, r_in, base_a)
        d += (
            f" M {bx1:.2f} {by1:.2f} L {tx:.2f} {ty:.2f} "
            f"L {bx2:.2f} {by2:.2f} Z"
        )
    return d


def short_label(name: str, limit: int = 24) -> str:
    """Harita etiketi için kısa ad.

    Parça adları açıklayıcı olduğu için uzun ("kanR / nptII (KanR,
    aminoglikozid fosfotransferaz)"). Haritada bunlar tuvalden taşar;
    parantez içini ve ikinci adı atıp kısaltıyoruz. Tam ad SVG `<title>`
    içinde kalır, fare üzerine gelince görünür.
    """
    s = name.split("(")[0].strip()
    s = s.split(" / ")[0].strip()
    s = s.rstrip(",").strip()
    if len(s) > limit:
        s = s[: limit - 1].rstrip() + "…"
    return s or name[:limit]


def _declutter(labels: list[dict], min_gap: float = 16.5) -> list[dict]:
    """Aynı tarafa düşen etiketleri dikeyde ayır."""
    for side in ("l", "r"):
        group = sorted([x for x in labels if x["side"] == side],
                       key=lambda x: x["y"])
        for i in range(1, len(group)):
            if group[i]["y"] - group[i - 1]["y"] < min_gap:
                group[i]["y"] = group[i - 1]["y"] + min_gap
        # Taşmayı yukarı doğru geri yay
        for i in range(len(group) - 2, -1, -1):
            if group[i + 1]["y"] - group[i]["y"] < min_gap:
                group[i]["y"] = group[i + 1]["y"] - min_gap
    return labels


def circular_map(
    name: str,
    length: int,
    features: list[Feature],
    sites: dict[str, list[int]] | None = None,
    size: int = 660,
    show_legend: bool = True,
) -> str:
    """Dairesel plazmit haritası döndür (SVG string)."""
    if length <= 0:
        return "<svg xmlns='http://www.w3.org/2000/svg'></svg>"

    sites = sites or {}
    # Tuval daireden geniş: etiketler iki yana açıldığı için yan boşluk gerekir.
    H = size
    W = int(size * 1.42)
    cx, cy = W / 2, H / 2
    R = size * 0.26          # iskelet yarıçapı
    band = size * 0.038      # özellik bandı kalınlığı
    label_r = R + band * 1.6
    label_edge = size * 0.40  # etiket metninin başladığı yatay mesafe

    parts: list[str] = []

    # --- İskelet
    parts.append(
        f'<circle cx="{cx}" cy="{cy}" r="{R}" fill="none" '
        f'stroke="currentColor" stroke-opacity="0.28" stroke-width="2.5"/>'
    )

    # --- Konum çizgileri (ölçek)
    step = _nice_step(length)
    for pos in range(0, length, step):
        a = _ang(pos + 1, length)
        x0, y0 = _pol(cx, cy, R - band * 0.55, a)
        x1, y1 = _pol(cx, cy, R - band * 0.95, a)
        parts.append(
            f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
            f'stroke="currentColor" stroke-opacity="0.3" stroke-width="1"/>'
        )
        tx, ty = _pol(cx, cy, R - band * 1.7, a)
        parts.append(
            f'<text x="{tx:.1f}" y="{ty:.1f}" font-size="9" fill="currentColor" '
            f'fill-opacity="0.45" text-anchor="middle" dominant-baseline="middle">'
            f'{pos if pos else 1}</text>'
        )

    # --- Özellikler
    labels: list[dict] = []
    for f in features:
        a0 = _ang(f.start, length)
        a1 = _ang(f.end + 1, length)
        if f.end < f.start:  # dairesel sarma
            a1 = _ang(f.end + 1 + length, length)
        color = CATEGORY_COLOR.get(f.category, "#71717A")
        inner = R - band / 2
        outer = R + band / 2
        d = _arc_band(cx, cy, inner, outer, a0, a1, arrow=f.strand)
        tip = html.escape(
            f"{f.name} | {CATEGORY_TR.get(f.category, f.category)} | "
            f"{f.start}..{f.end} ({'+' if f.strand > 0 else '-'})"
        )
        parts.append(
            f'<path d="{d}" fill="{color}" fill-opacity="0.85" '
            f'stroke="{color}" stroke-width="0.5"><title>{tip}</title></path>'
        )

        # Etiket konumu
        mid = a0 + ((a1 - a0) % (2 * math.pi)) / 2
        lx, ly = _pol(cx, cy, label_r, mid)
        labels.append({
            "x": lx, "y": ly, "angle": mid, "text": short_label(f.name),
            "full": f.name, "color": color,
            "side": "r" if math.cos(mid) >= 0 else "l",
        })

    # --- Etiketler (çakışma giderme + kılavuz çizgi)
    labels = _declutter(labels)
    for lb in labels:
        ax, ay = _pol(cx, cy, R + band * 0.6, lb["angle"])
        edge = cx + (label_edge if lb["side"] == "r" else -label_edge)
        anchor = "start" if lb["side"] == "r" else "end"
        text_x = edge + (6 if lb["side"] == "r" else -6)
        bend_x = cx + (R + band) * 1.12 * math.cos(lb["angle"])
        parts.append(
            f'<polyline points="{ax:.1f},{ay:.1f} {bend_x:.1f},{lb["y"]:.1f} '
            f'{edge:.1f},{lb["y"]:.1f}" fill="none" stroke="{lb["color"]}" '
            f'stroke-opacity="0.5" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{text_x:.1f}" y="{lb["y"]:.1f}" font-size="11.5" '
            f'fill="currentColor" text-anchor="{anchor}" '
            f'dominant-baseline="middle">{html.escape(lb["text"])}'
            f'<title>{html.escape(lb["full"])}</title></text>'
        )

    # --- Restriksiyon bölgeleri
    site_items = [(e, p) for e, pos in sites.items() for p in pos]
    for enz, pos in site_items:
        a = _ang(pos, length)
        x0, y0 = _pol(cx, cy, R + band * 0.55, a)
        x1, y1 = _pol(cx, cy, R + band * 1.15, a)
        parts.append(
            f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{x1:.1f}" y2="{y1:.1f}" '
            f'stroke="currentColor" stroke-opacity="0.55" stroke-width="1.2"/>'
        )

    # --- Merkez bilgi
    parts.append(
        f'<text x="{cx}" y="{cy - 10}" font-size="17" font-weight="600" '
        f'fill="currentColor" text-anchor="middle">{html.escape(name)}</text>'
    )
    parts.append(
        f'<text x="{cx}" y="{cy + 12}" font-size="12.5" fill="currentColor" '
        f'fill-opacity="0.6" text-anchor="middle">{length} bp</text>'
    )
    if site_items:
        parts.append(
            f'<text x="{cx}" y="{cy + 30}" font-size="10.5" '
            f'fill="currentColor" fill-opacity="0.45" text-anchor="middle">'
            f'{len(site_items)} kesim bölgesi</text>'
        )

    legend, legend_h = _legend_svg(features, W, H) if show_legend else ("", 0)
    total_h = H + legend_h

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {total_h}" '
        f'width="100%" role="img" style="max-width:{W}px;display:block;'
        f'margin:0 auto;font-family:system-ui,-apple-system,sans-serif">'
        + "".join(parts) + legend + "</svg>"
    )


def _legend_svg(features: list[Feature], width: int, top: int) -> tuple[str, int]:
    """Kategori açıklaması; (svg, kapladığı yükseklik) döndürür."""
    cats: list[str] = []
    for f in features:
        if f.category not in cats:
            cats.append(f.category)
    if not cats:
        return "", 0

    out: list[str] = []
    x, rows = 14.0, 1
    y = top + 14
    for c in cats:
        color = CATEGORY_COLOR.get(c, "#71717A")
        label = CATEGORY_TR.get(c, c)
        w_needed = 16 + len(label) * 6.2
        if x + w_needed > width - 14 and x > 14.0:
            x = 14.0
            y += 16
            rows += 1
        out.append(
            f'<rect x="{x:.0f}" y="{y - 7:.0f}" width="9" height="9" rx="2" '
            f'fill="{color}" fill-opacity="0.85"/>'
        )
        out.append(
            f'<text x="{x + 14:.0f}" y="{y:.0f}" font-size="10.5" '
            f'fill="currentColor" fill-opacity="0.7" '
            f'dominant-baseline="middle">{html.escape(label)}</text>'
        )
        x += w_needed
    return "".join(out), rows * 16 + 12


def _nice_step(length: int) -> int:
    for s in (100, 250, 500, 1000, 2000, 2500, 5000, 10000, 20000):
        if length / s <= 12:
            return s
    return 50000


def linear_map(
    name: str,
    length: int,
    features: list[Feature],
    sites: dict[str, list[int]] | None = None,
    width: int = 900,
) -> str:
    """Doğrusal harita - kaset düzenini okumak dairesel haritadan kolaydır."""
    if length <= 0:
        return "<svg xmlns='http://www.w3.org/2000/svg'></svg>"
    sites = sites or {}
    pad = 26
    track_y = 86
    h = 150
    usable = width - 2 * pad

    def px(pos: int) -> float:
        return pad + usable * (pos - 1) / max(length, 1)

    parts = [
        f'<line x1="{pad}" y1="{track_y}" x2="{width - pad}" y2="{track_y}" '
        f'stroke="currentColor" stroke-opacity="0.3" stroke-width="2"/>'
    ]

    step = _nice_step(length)
    for pos in range(0, length + 1, step):
        x = px(max(pos, 1))
        parts.append(
            f'<line x1="{x:.1f}" y1="{track_y + 16}" x2="{x:.1f}" '
            f'y2="{track_y + 22}" stroke="currentColor" stroke-opacity="0.3"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{track_y + 34}" font-size="9.5" '
            f'fill="currentColor" fill-opacity="0.5" text-anchor="middle">'
            f'{pos if pos else 1}</text>'
        )

    lanes: list[int] = []
    for f in sorted(features, key=lambda z: z.start):
        x0, x1 = px(f.start), px(max(f.end, f.start + 1))
        if f.end < f.start:
            x1 = px(length)
        w = max(x1 - x0, 2.5)
        lane = 0
        while lane < len(lanes) and lanes[lane] > x0 - 4:
            lane += 1
        if lane == len(lanes):
            lanes.append(0)
        lanes[lane] = x1
        y = track_y - 11 - lane * 24
        color = CATEGORY_COLOR.get(f.category, "#71717A")
        notch = min(7.0, w * 0.3)
        if f.strand > 0:
            d = (f"M {x0:.1f} {y:.1f} L {x1 - notch:.1f} {y:.1f} "
                 f"L {x1:.1f} {y + 7:.1f} L {x1 - notch:.1f} {y + 14:.1f} "
                 f"L {x0:.1f} {y + 14:.1f} Z")
        else:
            d = (f"M {x1:.1f} {y:.1f} L {x0 + notch:.1f} {y:.1f} "
                 f"L {x0:.1f} {y + 7:.1f} L {x0 + notch:.1f} {y + 14:.1f} "
                 f"L {x1:.1f} {y + 14:.1f} Z")
        tip = html.escape(f"{f.name} | {f.start}..{f.end}")
        parts.append(
            f'<path d="{d}" fill="{color}" fill-opacity="0.85">'
            f'<title>{tip}</title></path>'
        )
        if w > 46:
            fits = short_label(f.name, limit=max(3, int((w - 10) / 6.2)))
            parts.append(
                f'<text x="{(x0 + x1) / 2:.1f}" y="{y + 7:.1f}" font-size="10" '
                f'fill="#fff" text-anchor="middle" dominant-baseline="middle" '
                f'pointer-events="none">{html.escape(fits)}</text>'
            )
        h = max(h, abs(y) + 60)

    for enz, pos_list in sites.items():
        for pos in pos_list:
            x = px(pos)
            parts.append(
                f'<line x1="{x:.1f}" y1="{track_y - 6}" x2="{x:.1f}" '
                f'y2="{track_y + 6}" stroke="currentColor" '
                f'stroke-opacity="0.6" stroke-width="1.2"/>'
            )

    top = min([track_y - 11 - i * 24 for i in range(max(len(lanes), 1))] + [40])
    vb_top = top - 34
    vb_h = (track_y + 48) - vb_top
    parts.insert(
        0,
        f'<text x="{pad}" y="{vb_top + 16}" font-size="13" font-weight="600" '
        f'fill="currentColor">{html.escape(name)} - {length} bp</text>',
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'viewBox="0 {vb_top} {width} {vb_h}" width="100%" role="img" '
        f'style="display:block;font-family:system-ui,-apple-system,sans-serif">'
        + "".join(parts) + "</svg>"
    )
