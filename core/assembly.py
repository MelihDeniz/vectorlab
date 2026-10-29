"""Parça listesinden vektör montajı.

Stub parçalar (dizisi gömülü olmayan uzun CDS'ler) doğru uzunlukta N bloğu
olarak yerleştirilir ve açıkça işaretlenir. Böylece harita, uzunluk ve
kopya sayısı hesapları doğru kalır; sadece o bölgenin dizi düzeyinde
analizi (enzim kesimi vb.) yapılamaz - rapor bunu söyler.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core.seq_utils import Feature, clean_seq, gc_percent


@dataclass
class AssembledVector:
    name: str
    sequence: str
    features: list[Feature] = field(default_factory=list)
    circular: bool = True
    stub_regions: list[dict[str, Any]] = field(default_factory=list)
    part_order: list[dict[str, Any]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def gc(self) -> float:
        return gc_percent(self.sequence.replace("N", ""))

    @property
    def resolved_fraction(self) -> float:
        """Dizisi gerçekten bilinen bölgenin oranı (%)."""
        if not self.sequence:
            return 0.0
        known = sum(1 for c in self.sequence if c != "N")
        return round(100 * known / len(self.sequence), 1)

    def has_category(self, category: str) -> bool:
        return any(p["category"] == category for p in self.part_order)

    def parts_of(self, category: str) -> list[dict[str, Any]]:
        return [p for p in self.part_order if p["category"] == category]


#: Tipik bir bakteriyel ifade vektörünün mantıklı parça sırası.
#: Tasarımcı sayfası yeni yapı kurarken bu sırayı önerir.
CANONICAL_ORDER = [
    # Bileşik parça kendi promotörünü ve terminatörünü taşır, bu yüzden
    # kasetin başında durur ve tek başına çalışır.
    "device",
    "promoter",
    "operator",
    "rbs",
    # Sinyal peptidi proteinin N ucunda OLMAK ZORUNDA: taşıma
    # makinesi onu ilk 20-30 amino asitte arar. Etiketten de önce
    # gelir, çünkü etiket araya girerse sinyal tanınmaz.
    "signal",
    "tag",
    "protease_site",
    "mcs",
    "linker",
    "cds",
    "reporter",
    "terminator",
    "ori",
    "marker",
    "regulator",
]


def assemble(
    name: str,
    parts: list[dict[str, Any]],
    insert_seq: str = "",
    insert_name: str = "insert",
    insert_after: str | None = None,
    circular: bool = True,
) -> AssembledVector:
    """Parçaları verilen sırada birleştir.

    insert_after: hangi parçanın ardından insert yerleştirilecek
                  (None ise MCS varsa onun yerine, yoksa sona eklenir).
    """
    insert_seq = clean_seq(insert_seq)
    chunks: list[str] = []
    features: list[Feature] = []
    stubs: list[dict[str, Any]] = []
    order: list[dict[str, Any]] = []
    cursor = 0

    def place(seq: str, label: str, category: str, part_id: str, is_stub: bool):
        nonlocal cursor
        if not seq:
            return
        start = cursor + 1
        end = cursor + len(seq)
        chunks.append(seq)
        features.append(Feature(
            start=start, end=end, strand=1, name=label,
            category=category, part_id=part_id,
        ))
        if is_stub:
            stubs.append({
                "name": label, "part_id": part_id,
                "start": start, "end": end, "length": len(seq),
            })
        cursor = end

    # Insert nereye girecek?
    mcs_index = next(
        (i for i, p in enumerate(parts) if p["category"] == "mcs"), None
    )
    if insert_after:
        anchor = next(
            (i for i, p in enumerate(parts) if p["part_id"] == insert_after), None
        )
    elif mcs_index is not None:
        anchor = mcs_index
    else:
        anchor = len(parts) - 1 if parts else None

    for i, p in enumerate(parts):
        is_stub = p.get("seq_status") != "full"
        seq = (
            "N" * int(p.get("length") or 0)
            if is_stub
            else clean_seq(p.get("sequence") or "")
        )
        place(seq, p["name"], p["category"], p["part_id"], is_stub)
        order.append({
            "part_id": p["part_id"], "name": p["name"],
            "category": p["category"], "seq_status": p.get("seq_status"),
            "length": len(seq),
        })

        if insert_seq and anchor is not None and i == anchor:
            place(insert_seq, insert_name, "cds", "INSERT", False)
            order.append({
                "part_id": "INSERT", "name": insert_name, "category": "cds",
                "seq_status": "full", "length": len(insert_seq),
            })

    if insert_seq and anchor is None:
        place(insert_seq, insert_name, "cds", "INSERT", False)
        order.append({
            "part_id": "INSERT", "name": insert_name, "category": "cds",
            "seq_status": "full", "length": len(insert_seq),
        })

    return AssembledVector(
        name=name,
        sequence="".join(chunks),
        features=features,
        circular=circular,
        stub_regions=stubs,
        part_order=order,
    )


def sort_canonical(parts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Parçaları tipik vektör düzenine göre sırala."""
    rank = {c: i for i, c in enumerate(CANONICAL_ORDER)}
    return sorted(parts, key=lambda p: rank.get(p["category"], 99))


def to_genbank(vec: AssembledVector, organism: str = "Escherichia coli") -> str:
    """Yapıyı GenBank formatında yaz (SnapGene/Benchling bunu okur)."""
    from datetime import date

    topology = "circular" if vec.circular else "linear"
    today = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        # GenBank LOCUS satiri sutun duyarlidir: ad 13-28, uzunluk
        # 30-40 arasinda saga dayali, "bp" 42-43, molekul tipi 48,
        # topoloji 56. Kaydirirsaniz Biopython "malformed locus line"
        # uyarisi verir ve alanlari yanlis okur.
        f"LOCUS       {vec.name[:16]:<16} {vec.length:>11} bp    DNA     "
        f"{topology:<8} SYN {today}",
        f"DEFINITION  {vec.name} - VectorLab ile tasarlandi.",
        "SOURCE      synthetic DNA construct",
        f"  ORGANISM  {organism}",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{vec.length}",
        f'                     /organism="{organism}"',
        '                     /mol_type="other DNA"',
    ]

    kind = {
        "promoter": "promoter", "operator": "protein_bind", "rbs": "RBS",
        "cds": "CDS", "tag": "CDS", "protease_site": "misc_feature",
        "linker": "misc_feature", "terminator": "terminator",
        "ori": "rep_origin", "marker": "CDS", "regulator": "CDS",
        "mcs": "misc_feature", "reporter": "CDS",
    }

    for f in sorted(vec.features, key=lambda x: x.start):
        gb_type = kind.get(f.category, "misc_feature")
        loc = f"{f.start}..{f.end}"
        if f.strand < 0:
            loc = f"complement({loc})"
        lines.append(f"     {gb_type:<15} {loc}")
        lines.append(f'                     /label="{f.name}"')
        if f.part_id:
            lines.append(f'                     /note="part_id:{f.part_id}"')

    lines.append("ORIGIN")
    seq = vec.sequence.lower()
    for i in range(0, len(seq), 60):
        block = seq[i : i + 60]
        groups = " ".join(block[j : j + 10] for j in range(0, len(block), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines)


def to_fasta(vec: AssembledVector, width: int = 70) -> str:
    header = f">{vec.name} | {vec.length} bp | {'dairesel' if vec.circular else 'dogrusal'}"
    body = "\n".join(
        vec.sequence[i : i + width] for i in range(0, len(vec.sequence), width)
    )
    return f"{header}\n{body}\n"
