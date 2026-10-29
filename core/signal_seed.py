"""Sinyal peptitlerini kataloğa verir.

Veri `data/signal_peptides.json` içinde durur ve
`tools/fetch_signal_peptides.py` ile üretilir. Kesim noktası
UniProt'un `Signal` öznitelik aralığından alınır, nükleotid dizisi
NCBI'dan gelir ve çevirisinin UniProt'un bildirdiği peptide eşit
olduğu doğrulanmıştır.

Sinyal peptidi neden kendi kategorisinde: afinite etiketinden
işlevsel olarak farklıdır. Proteinin N ucunda OLMAK ZORUNDADIR,
taşıma sırasında kesilip atılır ve hangi yoldan taşınacağını belirler.
His etiketi gibi bir yere eklenip orada kalan bir parça değildir.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Any

DATA = Path(__file__).resolve().parent.parent / "data" / "signal_peptides.json"

#: Taşıma yoluna göre ne beklenmeli.
PATHWAY_NOTE = {
    "Sec (post-translasyonel)": (
        "Protein sitoplazmada KATLANMADAN, açık zincir hâlinde zara "
        "taşınır ve periplazmada katlanır. Hızlı katlanan ya da "
        "kofaktör bağlayan proteinler bu yolda tıkanır."),
    "SRP (ko-translasyonel)": (
        "Protein daha ribozomdan çıkarken zara yönlendirilir, yani "
        "sitoplazmada bekleme süresi neredeyse yoktur. Sitoplazmada "
        "toplanıp çökelen hedeflerde Sec yolundan daha iyi sonuç "
        "verebilir."),
    "Tat (twin-arginine)": (
        "Tek fark yaratan yol: KATLANMIŞ proteini taşır. Kofaktör "
        "bağlamış ya da sitoplazmada katlanması gereken proteinler "
        "için başka seçenek yoktur. Buna karşılık kapasitesi düşüktür, "
        "yüksek verim beklemeyin."),
    "Sec (Gram-pozitif)": (
        "Gram-pozitif bakterilerde dış zar yoktur: protein periplazmada "
        "sıkışmaz, doğrudan besiyerine salgılanır. Saflaştırma çok "
        "kolaylaşır - hücreyi parçalamanız bile gerekmez."),
}


@functools.lru_cache(maxsize=1)
def _raw() -> dict[str, Any]:
    if not DATA.exists():
        return {"peptides": []}
    return json.loads(DATA.read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=1)
def signal_parts() -> tuple[dict[str, Any], ...]:
    out: list[dict[str, Any]] = []
    for rec in _raw().get("peptides", []):
        yol_notu = PATHWAY_NOTE.get(rec["pathway"], "")
        gram_pozitif = "Gram-pozitif" in rec["pathway"]
        out.append({
            "part_id": rec["part_id"],
            "name": f"{rec['name']} ({rec['length_aa']} aa)",
            "category": "signal",
            "sequence": rec["sequence"],
            "seq_status": "full",
            "length": rec["length"],
            "strength": rec["pathway"],
            "induction": None,
            "requires": (["Gram-pozitif konak (B. subtilis vb.)"]
                         if gram_pozitif
                         else ["Gram-negatif konak (E. coli vb.)"]),
            "incompatible_with": [],
            "copy_number": None,
            "inc_group": None,
            "notes": (
                f"{rec['note']}\n\n**Taşıma yolu — {rec['pathway']}:** "
                f"{yol_notu}\n\n"
                f"Peptit dizisi `{rec['peptide']}`; taşıma sırasında "
                f"kesilir, olgun proteinde kalmaz. Genin başına, "
                f"başlama kodonunun yerine gelir - araya başka bir şey "
                f"girerse sinyal tanınmaz.\n\n"
                f"Kesim noktası UniProt {rec['uniprot']} "
                f"({rec['uniprot_entry']}) kaydındaki işaretli aralıktan "
                f"alındı; nükleotid dizisinin çevirisi bu peptide eşit "
                f"olduğu doğrulandı."
            ),
            "accession": (f"UniProt {rec['uniprot']} · "
                          f"{rec['nucleotide']}"),
            "pmid": "",
            "source": "signal",
            "pathway": rec["pathway"],
            "gram_positive": gram_pozitif,
        })
    return tuple(out)


def by_pathway() -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for p in signal_parts():
        out.setdefault(p["pathway"], []).append(p)
    return out


def source_note() -> str:
    return _raw().get("note", "")
