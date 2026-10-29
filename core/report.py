"""Tasarımı tek sayfalık, yazdırılabilir bir rapora çevirir.

Sunumda ya da ödevde gereken şey ekran görüntüsü değil, kendi başına
duran bir belgedir: harita, parça listesi, denetim sonucu, protokol ve
dizinin nereden geldiği. Üretilen HTML tek dosyadır (harita SVG olarak
içine gömülür), tarayıcıda açılır ve oradan PDF'e basılabilir.
"""

from __future__ import annotations

import datetime as _dt
import html
from typing import Any

from core import protocol as pr
from core import provenance
from core.parts_seed import CATEGORY_TR

_CSS = """
:root{--ink:#16161a;--mut:#5b5b66;--line:#d8d8e0;--bg:#fff;--accent:#1f6feb;
--err:#b4232b;--warn:#8a5a00;--ok:#1a6b3c}
*{box-sizing:border-box}
body{margin:0;padding:32px 40px;background:var(--bg);color:var(--ink);
font:14px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",sans-serif;
max-width:900px;margin-inline:auto}
h1{font-size:25px;margin:0 0 4px;letter-spacing:-.01em}
h2{font-size:17px;margin:30px 0 10px;padding-bottom:5px;
border-bottom:1px solid var(--line)}
h3{font-size:14px;margin:18px 0 6px}
p{margin:8px 0}
.sub{color:var(--mut);margin:0 0 20px;font-size:13px}
.grid{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0}
.card{flex:1 1 130px;border:1px solid var(--line);border-radius:8px;
padding:9px 12px}
.card .k{font-size:11px;color:var(--mut);text-transform:uppercase;
letter-spacing:.04em}
.card .v{font-size:19px;font-weight:600;margin-top:2px}
table{border-collapse:collapse;width:100%;margin:10px 0;font-size:13px}
th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--line);
vertical-align:top}
th{font-weight:600;color:var(--mut);font-size:11px;text-transform:uppercase;
letter-spacing:.04em}
.f{border-left:3px solid var(--line);padding:7px 12px;margin:9px 0;
background:#fafafc;border-radius:0 6px 6px 0}
.f.hata{border-left-color:var(--err)}
.f.uyari{border-left-color:var(--warn)}
.f.bilgi{border-left-color:var(--accent)}
.f b{display:block;margin-bottom:2px}
.f span{color:var(--mut);font-size:12.5px}
.map{text-align:center;margin:16px 0}
.map svg{max-width:100%;height:auto}
.step{margin:14px 0;page-break-inside:avoid}
.step .t{font-weight:600}
.step .c{color:var(--mut);font-size:12.5px;margin-top:3px}
code{background:#f2f2f6;padding:1px 5px;border-radius:4px;font-size:12.5px}
footer{margin-top:34px;padding-top:12px;border-top:1px solid var(--line);
color:var(--mut);font-size:12px}
@media print{body{padding:0;max-width:none}h2{page-break-after:avoid}
.card{break-inside:avoid}}
"""


def _esc(x: Any) -> str:
    return html.escape(str(x if x is not None else ""))


def _cards(pairs: list[tuple[str, Any]]) -> str:
    out = ["<div class='grid'>"]
    for k, v in pairs:
        out.append(f"<div class='card'><div class='k'>{_esc(k)}</div>"
                   f"<div class='v'>{_esc(v)}</div></div>")
    out.append("</div>")
    return "".join(out)


def _table(headers: list[str], rows: list[list[Any]]) -> str:
    if not rows:
        return "<p class='sub'>Kayıt yok.</p>"
    out = ["<table><thead><tr>"]
    out += [f"<th>{_esc(h)}</th>" for h in headers]
    out.append("</tr></thead><tbody>")
    for r in rows:
        out.append("<tr>" + "".join(f"<td>{_esc(c)}</td>" for c in r) + "</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def build_html(
    *,
    name: str,
    length: int,
    gc: float,
    host_name: str = "",
    purpose: str = "",
    parts: list[dict[str, Any]] | None = None,
    map_svg: str = "",
    findings: list[Any] | None = None,
    protein_report: Any = None,
    proto: pr.Protocol | None = None,
    gel_svg: str = "",
    insert_name: str = "",
    insert_bp: int = 0,
    enzymes: list[str] | None = None,
    author: str = "",
) -> str:
    """Tek dosyalık HTML rapor üret."""
    parts = parts or []
    findings = findings or []
    enzymes = enzymes or []
    bugun = _dt.date.today().strftime("%d.%m.%Y")

    p: list[str] = [
        "<!doctype html><html lang='tr'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        f"<title>{_esc(name)} — vektör raporu</title>",
        f"<style>{_CSS}</style></head><body>",
        f"<h1>{_esc(name)}</h1>",
        f"<p class='sub'>Bakteriyel vektör tasarım raporu · {bugun}"
        + (f" · {_esc(author)}" if author else "")
        + f" · {provenance.describe()}</p>",
    ]

    kartlar = [("Uzunluk", f"{length} bp"), ("GC", f"%{gc}")]
    if host_name:
        kartlar.append(("Konak", host_name))
    if purpose:
        kartlar.append(("Amaç", purpose))
    if insert_bp:
        kartlar.append(("Insert", f"{insert_bp} bp"))
    if enzymes:
        kartlar.append(("Kesim", " + ".join(enzymes)))
    p.append(_cards(kartlar))

    if map_svg:
        p.append("<h2>Plazmit haritası</h2>")
        p.append(f"<div class='map'>{map_svg}</div>")

    if parts:
        p.append("<h2>Parçalar</h2>")
        p.append(_table(
            ["#", "Parça", "Tür", "Uzunluk", "Kaynak"],
            [[i + 1, x.get("name", ""),
              CATEGORY_TR.get(x.get("category", ""), x.get("category", "")),
              f"{x.get('length') or 0} bp", x.get("accession") or "—"]
             for i, x in enumerate(parts)],
        ))

    if findings:
        p.append("<h2>Tasarım denetimi</h2>")
        sayim = {"hata": 0, "uyari": 0, "bilgi": 0}
        for f in findings:
            sayim[getattr(f, "severity", "bilgi")] = sayim.get(
                getattr(f, "severity", "bilgi"), 0) + 1
        p.append(_cards([("Hata", sayim["hata"]), ("Uyarı", sayim["uyari"]),
                         ("Bilgi", sayim["bilgi"])]))
        for f in findings:
            sev = getattr(f, "severity", "bilgi")
            p.append(
                f"<div class='f {_esc(sev)}'><b>{_esc(getattr(f, 'title', ''))}"
                f" <code>{_esc(getattr(f, 'rule_id', ''))}</code></b>"
                f"<span>{_esc(getattr(f, 'detail', ''))}</span></div>"
            )

    if protein_report is not None:
        r = protein_report
        p.append("<h2>İfade edilecek protein</h2>")
        p.append(_cards([
            ("Uzunluk", f"{r.length_aa} aa"),
            ("Kütle", f"{r.mw_kda} kDa"),
            ("pI", r.pi),
            ("GRAVY", r.gravy),
            ("A280 (1 mg/mL)", r.a280_for_1mg_ml if r.a280_for_1mg_ml
             else "ölçülemez"),
        ]))
        for n in r.notes:
            p.append(
                f"<div class='f {_esc(n['kind'])}'><b>{_esc(n['title'])}</b>"
                f"<span>{_esc(n['detail'])}</span></div>"
            )

    if gel_svg:
        p.append("<h2>Beklenen jel görüntüsü</h2>")
        p.append(f"<div class='map'>{gel_svg}</div>")

    if proto is not None:
        p.append("<h2>Laboratuvar planı</h2>")
        if proto.warnings:
            for w in proto.warnings:
                p.append(f"<div class='f uyari'><b>Eksik</b>"
                         f"<span>{_esc(w)}</span></div>")
        for s in proto.steps:
            p.append(f"<div class='step'><div class='t'>{s.no}. "
                     f"{_esc(s.title)}</div><p>{_esc(s.detail)}</p>")
            if s.reagents:
                p.append(_table(["Bileşen", "Miktar", "Not"],
                                [[r.name, r.amount, r.note]
                                 for r in s.reagents]))
            if s.conditions:
                p.append(f"<div class='c'><b>Koşullar:</b> "
                         f"{_esc(s.conditions)}</div>")
            if s.warning:
                p.append(f"<div class='c'><b>Dikkat:</b> "
                         f"{_esc(s.warning)}</div>")
            p.append("</div>")

    p.append(
        "<footer>Diziler NCBI Nucleotide ve iGEM Registry'den (CC BY 4.0) "
        "alınmıştır; restriksiyon enzimi verisi REBASE'ten gelir. "
        "Protokoldeki derişimler standart E. coli çalışma değerleridir, "
        "laboratuvarınızın kendi protokolü önceliklidir. Bu rapor bir "
        "tasarım belgesidir; deney sonucu değildir.</footer>"
        "</body></html>"
    )
    return "".join(p)
