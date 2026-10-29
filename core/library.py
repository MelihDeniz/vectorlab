"""Kombinatoryal kütüphane: tek tasarım yerine bir seri.

Gerçek laboratuvarda ifade düzeyi tek bir tasarımla ayarlanmaz. Promotör
ve ribozom bağlanma bölgesi değiştirilerek bir seri kurulur, hepsi aynı
anda denenir ve işe yarayan seçilir. Bunun sebebi basit: ifade düzeyini
önceden hesaplayamazsınız, ama bir aralığı tarayabilirsiniz.

Bu modül o seriyi kurar. Sabit tuttuğunuz parçaların üstüne, değişken
bıraktığınız kategorilerde bütün birleşimleri üretir, her birini
denetimden geçirir ve beklenen ifade düzeyine göre sıralar.

**Beklenen ifade düzeyi bir ÖLÇÜM DEĞİLDİR.** Üç kaynaktan türetilmiş
sıralamadır ve her satırda hangisinden geldiği yazılıdır:

  1. Anderson promotörlerinde sigma70 konsensüsüne uzaklık - diziden
     hesaplanır, en nesnel olanı budur.
  2. Katalogdaki güç etiketi ("çok güçlü", "orta", "zayıf").
  3. iGEM Registry başlığındaki strong / medium / weak nitelemesi.

Bunlar ifade düzeyiyle aynı yönde gider ama kaç kat fark edeceğini
söylemez. Serinin amacı zaten budur: hesaplayamadığınız için
tarıyorsunuz.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from typing import Any

from core.assembly import assemble, sort_canonical
from core.validate import summarize, validate

# Ordinal basamaklar. Sayılar yalnızca sıralama içindir.
TIER_LABELS = [
    (0.85, "çok yüksek"), (0.65, "yüksek"), (0.40, "orta"),
    (0.20, "düşük"), (0.0, "çok düşük"),
]

_TR_STRENGTH = [
    (r"çok güçlü", 0.95), (r"güçlü", 0.78), (r"orta[- ]güçlü", 0.65),
    (r"orta", 0.45), (r"zayıf", 0.15), (r"ayarlanabilir", 0.60),
]
_EN_STRENGTH = [
    (r"\bstrong\b", 0.85), (r"\bmedium\b", 0.45),
    (r"\bweaker\b", 0.10), (r"\bweak\b", 0.20),
]


@dataclass
class Variant:
    name: str
    part_ids: list[str]
    length: int
    gc: float
    expression: float | None          # 0-1 ordinal
    tier: str
    basis: list[str] = field(default_factory=list)
    errors: int = 0
    warnings: int = 0
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.errors == 0


def tier_for(score: float | None) -> str:
    if score is None:
        return "bilinmiyor"
    for esik, ad in TIER_LABELS:
        if score >= esik:
            return ad
    return "çok düşük"


def strength_score(part: dict[str, Any]) -> tuple[float | None, str]:
    """Parçanın göreli gücü ve bu sayının nereden geldiği.

    Öncelik sırası bilinçli: diziden hesaplanan değer, elle yazılmış
    etiketten güvenilirdir.
    """
    pid = part.get("part_id", "")
    ad = part.get("name", "")

    # 1) Anderson ailesi: konsensüs uzaklığından hesapla
    if re.search(r"J231\d\d", pid) or re.search(r"J231\d\d", ad):
        from core.igem_seed import consensus_distance
        d = consensus_distance(part.get("sequence", ""))
        if d:
            # 12 olası fark; sıfır fark = konsensüsün kendisi
            skor = max(0.0, min(1.0, (12 - d["total"]) / 12))
            return round(skor, 3), f"konsensüs uzaklığı {d['total']}/12"

    # 2) Katalogdaki Türkçe güç etiketi
    etiket = (part.get("strength") or "")
    for desen, skor in _TR_STRENGTH:
        if re.search(desen, etiket, re.I):
            return skor, f"katalog etiketi: {etiket[:28]}"

    # 3) iGEM başlığındaki niteleme
    for desen, skor in _EN_STRENGTH:
        if re.search(desen, ad, re.I):
            m = re.search(desen, ad, re.I)
            return skor, f"Registry başlığı: {m.group(0)}"

    return None, "bilinmiyor"


def combined_expression(
    parts: list[dict[str, Any]]
) -> tuple[float | None, list[str]]:
    """Promotör ve RBS güçlerini birleştir.

    Çarpım kullanılıyor: transkripsiyon ve çeviri ardışık adımlardır,
    biri sıfıra yaklaşırsa diğeri telafi edemez. Bu, mekanizmaya uygun
    ama yine de kaba bir yaklaşımdır.
    """
    kaynaklar: list[str] = []
    carpim = 1.0
    bulundu = False
    for kategori in ("promoter", "rbs"):
        ilgili = [p for p in parts if p["category"] == kategori]
        if not ilgili:
            continue
        skorlar = []
        for p in ilgili:
            s, kaynak = strength_score(p)
            if s is not None:
                skorlar.append(s)
                kaynaklar.append(
                    f"{p['name'].split('(')[0].strip()[:22]}: {kaynak}")
        if skorlar:
            carpim *= max(skorlar)
            bulundu = True
    if not bulundu:
        return None, kaynaklar
    return round(carpim, 3), kaynaklar


def build(
    fixed_part_ids: list[str],
    variable: dict[str, list[str]],
    host: dict[str, Any] | None = None,
    purpose: str = "ifade",
    insert_seq: str = "",
    max_variants: int = 240,
    name_prefix: str = "var",
) -> tuple[list[Variant], list[str]]:
    """Bütün birleşimleri üret, denetle ve sırala.

    `variable` sözlüğü {kategori: [parça kimlikleri]} biçimindedir.
    Dönen ikinci değer, üretim sırasında çıkan uyarılardır.
    """
    from core.db import get_parts

    uyarilar: list[str] = []
    kategoriler = [k for k, v in variable.items() if v]
    if not kategoriler:
        return [], ["Değişken bırakılan kategori yok."]

    toplam = 1
    for k in kategoriler:
        toplam *= len(variable[k])
    if toplam > max_variants:
        uyarilar.append(
            f"{toplam} birleşim çıkıyor; ilk {max_variants} tanesi "
            f"üretildi. Değişken sayısını azaltın ya da her kategoriden "
            f"daha az parça seçin.")

    sabit = get_parts(fixed_part_ids)
    varyantlar: list[Variant] = []

    for i, secim in enumerate(itertools.product(
            *[variable[k] for k in kategoriler])):
        if i >= max_variants:
            break
        ids = list(fixed_part_ids) + list(secim)
        parcalar = sort_canonical(sabit + get_parts(list(secim)))
        vec = assemble(f"{name_prefix}-{i + 1:03d}", parcalar,
                       insert_seq=insert_seq)
        bulgular = validate(vec, host=host, purpose=purpose,
                            insert_seq=insert_seq)
        sayim = summarize(bulgular)
        skor, kaynaklar = combined_expression(parcalar)

        etiketler = []
        for k, pid in zip(kategoriler, secim):
            rec = next((p for p in parcalar if p["part_id"] == pid), None)
            if rec:
                etiketler.append(rec["name"].split("(")[0].strip()[:18])

        varyantlar.append(Variant(
            name=f"{name_prefix}-{i + 1:03d}",
            part_ids=ids, length=vec.length, gc=vec.gc,
            expression=skor, tier=tier_for(skor), basis=kaynaklar,
            errors=sayim["hata"], warnings=sayim["uyari"],
            notes=[b.title for b in bulgular if b.severity == "hata"],
        ))
        varyantlar[-1].name += " · " + " + ".join(etiketler)

    varyantlar.sort(key=lambda v: (-(v.expression or -1), v.errors))
    return varyantlar, uyarilar


def span(variants: list[Variant]) -> dict[str, Any]:
    """Serinin kapsadığı aralık - kütüphanenin asıl değeri budur."""
    skorlar = [v.expression for v in variants if v.expression is not None]
    if not skorlar:
        return {"min": None, "max": None, "fold": None, "tiers": []}
    en_az, en_cok = min(skorlar), max(skorlar)
    return {
        "min": en_az, "max": en_cok,
        "fold": round(en_cok / en_az, 1) if en_az > 0 else None,
        "tiers": sorted({v.tier for v in variants}),
    }


def to_csv(variants: list[Variant]) -> str:
    satirlar = ["ad,parcalar,uzunluk,gc,beklenen_ifade,basamak,hata,uyari"]
    for v in variants:
        satirlar.append(
            f'"{v.name}","{" ".join(v.part_ids)}",{v.length},{v.gc},'
            f'{v.expression if v.expression is not None else ""},'
            f'{v.tier},{v.errors},{v.warnings}')
    return "\n".join(satirlar) + "\n"
