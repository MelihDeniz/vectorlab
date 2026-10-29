"""Referans vektörlerden çıkarılmış parçaları kataloğa verir.

Veri `data/extracted_parts.json` içinde durur ve
`tools/extract_parts_from_vectors.py` ile üretilir. Bu parçaların elle
yazılanlardan farkı, her birinin **nereden geldiğini koordinatıyla
birlikte** taşımasıdır: `JX560326:381..1196` gibi. Yani dizi tartışmaya
açık değil, kaydı açıp bakabilirsiniz.

Aynı dizi birden çok vektörde geçiyorsa tek parça olarak tutulur ve
notunda hangi vektörlerde bulunduğu yazar - örneğin SEVA'nın `neo`
markerı 35 iskelette aynıdır, bu da koleksiyonun modüler olduğunu
doğrudan gösterir.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Any

DATA = Path(__file__).resolve().parent.parent / "data" / "extracted_parts.json"


@functools.lru_cache(maxsize=1)
def _raw() -> dict[str, Any]:
    if not DATA.exists():
        return {"parts": []}
    return json.loads(DATA.read_text(encoding="utf-8"))


@functools.lru_cache(maxsize=1)
def extracted_parts() -> tuple[dict[str, Any], ...]:
    """Katalog kaydı biçiminde, vektörlerden çıkarılmış parçalar."""
    out: list[dict[str, Any]] = []
    for rec in _raw().get("parts", []):
        dizi = (rec.get("sequence") or "").upper()
        if not dizi:
            continue
        paylasim = rec.get("also_in") or []
        notlar = [
            f"{rec['source_vector']} vektörünün GenBank öznitelik "
            f"tablosundan çıkarıldı ({rec['feature_type']} tipi). "
            f"Koordinat ve aksesyon kaynak alanında yazılı, yani diziyi "
            f"kaydın kendisinden doğrulayabilirsiniz."
        ]
        if paylasim:
            n = len(paylasim) + 1
            ilk = ", ".join(paylasim[:4])
            notlar.append(
                f"Aynı dizi {n} vektörde birebir aynı: {rec['source_vector']}, "
                f"{ilk}" + (" ve diğerleri." if len(paylasim) > 4 else ".")
                + " Bu, parçanın gerçekten standart bir modül olduğunu "
                  "gösterir."
            )
        if rec.get("note"):
            notlar.append(f"Kayıttaki not: {rec['note']}")

        out.append({
            "part_id": rec["part_id"],
            "name": rec["name"],
            "category": rec["category"],
            "sequence": dizi,
            "seq_status": "full",
            "length": len(dizi),
            "strength": None,
            "induction": None,
            "requires": [],
            "incompatible_with": [],
            "copy_number": None,
            "inc_group": None,
            "notes": "\n\n".join(notlar),
            "accession": rec["accession"],
            "pmid": "",
            "source": "vector",
        })
    return tuple(out)


def source_note() -> str:
    d = _raw()
    return d.get("note", "")


def generated() -> str:
    return _raw().get("generated", "")
