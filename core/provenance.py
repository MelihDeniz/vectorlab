"""Sürüm ve veri kimliği.

İki soruya cevap verir:

  - Bu yapıyı hangi sürüm üretti?
  - O sırada hangi veri dosyaları yüklüydü?

İkincisi önemsiz görünür ama değildir: NCBI kayıtları revize edilebilir
ve iGEM Registry güncellenir. Altı ay sonra aynı tasarımı açtığınızda
farklı bir sonuç çıkarsa, sebebin kodda mı veride mi olduğunu bilmek
istersiniz. Her veri dosyasının SHA-256 özeti `data/CHECKSUMS.txt`
içinde durur; `tools/update_checksums.py` onu yeniler, testler de
tutarlılığı sınar.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

VERSION = "1.0.0"

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
CHECKSUM_FILE = DATA / "CHECKSUMS.txt"

# Özeti alınan dosyalar: uygulamanın davranışını belirleyen veriler.
# Veritabanı dahil değildir, o her açılışta bunlardan yeniden kurulur.
TRACKED = ["igem_parts.json", "ecoli_codon_usage.json"]
TRACKED_DIRS = ["reference_vectors"]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for blok in iter(lambda: f.read(65536), b""):
            h.update(blok)
    return h.hexdigest()


def tracked_files() -> list[Path]:
    out: list[Path] = []
    for ad in TRACKED:
        p = DATA / ad
        if p.exists():
            out.append(p)
    for d in TRACKED_DIRS:
        klasor = DATA / d
        if klasor.exists():
            out.extend(sorted(klasor.glob("*.json")))
    return out


def compute_checksums() -> dict[str, str]:
    """Her izlenen dosyanın SHA-256 özeti (yol: depo köküne göre)."""
    return {str(p.relative_to(ROOT)).replace("\\", "/"): _sha256(p)
            for p in tracked_files()}


def write_checksums() -> int:
    toplam = compute_checksums()
    satirlar = [
        "# VectorLab veri dosyalarinin SHA-256 ozetleri.",
        "# Yenilemek icin: python tools/update_checksums.py",
        f"# surum: {VERSION}",
        "",
    ]
    satirlar += [f"{ozet}  {yol}" for yol, ozet in sorted(toplam.items())]
    CHECKSUM_FILE.write_text("\n".join(satirlar) + "\n", encoding="utf-8")
    return len(toplam)


def read_checksums() -> dict[str, str]:
    if not CHECKSUM_FILE.exists():
        return {}
    out: dict[str, str] = {}
    for satir in CHECKSUM_FILE.read_text(encoding="utf-8").splitlines():
        satir = satir.strip()
        if not satir or satir.startswith("#"):
            continue
        parca = satir.split(None, 1)
        if len(parca) == 2:
            out[parca[1].strip()] = parca[0].strip()
    return out


def verify() -> dict[str, list[str]]:
    """Diskteki dosyalarla kayıtlı özetleri karşılaştır."""
    kayitli = read_checksums()
    simdiki = compute_checksums()
    return {
        "degisen": sorted(y for y in kayitli
                          if y in simdiki and kayitli[y] != simdiki[y]),
        "eksik": sorted(y for y in kayitli if y not in simdiki),
        "yeni": sorted(y for y in simdiki if y not in kayitli),
    }


def data_fingerprint() -> str:
    """Bütün veri dosyalarını tek bir kısa kimliğe indirger."""
    toplam = compute_checksums()
    h = hashlib.sha256()
    for yol in sorted(toplam):
        h.update(yol.encode("utf-8"))
        h.update(toplam[yol].encode("ascii"))
    return h.hexdigest()[:12]


def stamp() -> str:
    """Kaydedilen yapılara yazılan sürüm izi."""
    return json.dumps(
        {"version": VERSION, "data": data_fingerprint()},
        ensure_ascii=False, separators=(",", ":"))


def describe() -> str:
    durum = verify()
    n = len(compute_checksums())
    if durum["degisen"] or durum["eksik"]:
        return (f"VectorLab {VERSION} · veri kimliği {data_fingerprint()} · "
                f"DİKKAT: {len(durum['degisen'])} dosya değişmiş, "
                f"{len(durum['eksik'])} dosya eksik")
    return (f"VectorLab {VERSION} · {n} veri dosyası · "
            f"kimlik {data_fingerprint()}")
