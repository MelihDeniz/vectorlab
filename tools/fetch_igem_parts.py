"""iGEM Registry'den seçili parçaları çeker ve data/igem_parts.json'a yazar.

Kaynak:  https://api.registry.igem.org/v1   (registry.igem.org'un kendi API'si)
Lisans:  parça kayıtları Creative Commons Attribution 4.0 (cc-by-4.0) ile
         dağıtılır. Betik her parçanın lisansını API'den okur ve kayda
         yazar; cc-by olmayan bir kayıt gelirse atlanır ve ekrana basılır.

Çalıştırmak için ağ gerekir:
    python tools/fetch_igem_parts.py
Çıktı dosyası depoda saklanır, uygulama çalışırken ağa çıkmaz.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
import time
from pathlib import Path

import requests

API = "https://api.registry.igem.org/v1"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli parca katalogu)"}
OUT = Path(__file__).resolve().parent.parent / "data" / "igem_parts.json"
VECTORS = Path(__file__).resolve().parent.parent / "data" / "reference_vectors"

# Bu kayitlar parca degil tam plazmittir: referans vektor klasorune,
# NCBI'dan gelenlerle ayni bicimde yazilirlar.
BACKBONES = {"pSB1C3": "PSB1C3", "pSB1A3": "PSB1A3",
             "pSB1K3": "PSB1K3", "pSB4C5": "PSB4C5"}

# Çekilecek parçalar. Kapsam bilinçli olarak bakteriyel (E. coli) sistemler.
WANTED: list[str] = (
    # Anderson sabit promotör ailesi - güç serisi
    [f"BBa_J231{n:02d}" for n in range(0, 20)]
    # Ribozom bağlanma bölgeleri - güç serisi
    + ["BBa_B0029", "BBa_B0030", "BBa_B0031", "BBa_B0032", "BBa_B0033",
       "BBa_B0034", "BBa_B0035"]
    # Terminatörler
    + ["BBa_B0010", "BBa_B0012", "BBa_B0014", "BBa_B0015", "BBa_B1002",
       "BBa_B1006"]
    # Düzenlenebilir promotörler
    + ["BBa_R0010", "BBa_R0011", "BBa_R0040", "BBa_R0051", "BBa_R0062",
       "BBa_R0063", "BBa_R0065", "BBa_R0071", "BBa_R0080", "BBa_R0082"]
    # Düzenleyici proteinler
    + ["BBa_C0012", "BBa_C0040", "BBa_C0051", "BBa_C0062", "BBa_C0079",
       "BBa_C0080"]
    # Reporter'lar ve sinyal proteinleri
    + ["BBa_E0040", "BBa_E1010", "BBa_E0020", "BBa_E0030", "BBa_C0061",
       "BBa_C0161"]
    # Bozunma (degradation) etiketleri
    + ["BBa_M0050", "BBa_M0051", "BBa_M0052"]
    # Sık kullanılan bileşik parçalar
    + ["BBa_I13504", "BBa_J04450", "BBa_E0240", "BBa_F2620", "BBa_I719005"]
    # BioBrick iskeletleri
    + ["pSB1C3", "pSB1A3", "pSB1K3", "pSB4C5"]
)


# API hız sınırı uyguluyor (429). İstekler arasına sabit bir bekleme koyar ve
# 429 gelirse Retry-After'a uyarak artan aralıklarla tekrar dener - aksi halde
# listenin yarısı sessizce "bulunamadı" olarak düşer.
PAUSE = 1.2
MAX_TRY = 6


def get(path: str, **params):
    delay = 5.0
    for attempt in range(1, MAX_TRY + 1):
        r = requests.get(f"{API}{path}", params=params or None,
                         headers=HEADERS, timeout=45)
        if r.status_code != 429:
            r.raise_for_status()
            time.sleep(PAUSE)
            return r.json()
        wait = float(r.headers.get("Retry-After") or delay)
        print(f"    hız sınırı, {wait:.0f} sn bekleniyor "
              f"(deneme {attempt}/{MAX_TRY})")
        sys.stdout.flush()
        time.sleep(wait)
        delay = min(delay * 2, 60)
    raise RuntimeError(f"429 sürüyor: {path} {params}")


def licenses() -> dict[str, dict]:
    out, page = {}, 1
    while True:
        d = get("/licenses", page=page, pageSize=50)
        rows = d.get("data") or []
        for row in rows:
            out[row["uuid"]] = row
        if len(rows) < 50:
            break
        page += 1
    return out


def write_backbones(records: list[dict]) -> int:
    """BioBrick iskeletlerini referans vektor kaydi olarak yaz."""
    VECTORS.mkdir(parents=True, exist_ok=True)
    n = 0
    for rec in records:
        vid = BACKBONES.get(rec["igem_id"])
        if not vid:
            continue
        (VECTORS / f"{vid}.json").write_text(json.dumps({
            "vector_id": vid,
            "accession": rec["igem_id"],
            "version": rec["igem_id"],
            "definition": rec["title"],
            "organism": "synthetic construct",
            "genbank_topology": (rec.get("topology") or "").lower(),
            "length": rec["length"],
            "sequence": rec["sequence"],
            "features": [],
            "fetched": _dt.date.today().isoformat(),
            "source": f"iGEM Registry {rec['igem_id']} "
                      f"({rec['license'].upper()})",
        }, ensure_ascii=False, indent=1), encoding="utf-8")
        n += 1
    return n


def main() -> int:
    lic = licenses()
    kept, missing, skipped, failed = [], [], [], []

    # Daha önce çekilmiş kayıtlar varsa onları koru (429 sonrası devam etmek
    # için): aynı parçayı iki kez indirmeye gerek yok.
    have: dict[str, dict] = {}
    if OUT.exists():
        try:
            have = {p["igem_id"]: p
                    for p in json.loads(OUT.read_text(encoding="utf-8"))["parts"]}
        except Exception:
            have = {}

    for name in WANTED:
        if name in have:
            kept.append(have[name])
            continue
        try:
            hits = (get("/parts", name=name).get("data") or [])
        except Exception as exc:
            print(f"  {name}: istek hatası {exc}")
            failed.append(name)
            continue
        exact = [h for h in hits if h["name"].lower() == name.lower()]
        if not exact:
            missing.append(name)
            continue
        try:
            det = get(f"/parts/{exact[0]['uuid']}")
        except Exception as exc:
            print(f"  {name}: detay hatası {exc}")
            failed.append(name)
            continue
        spdx = (lic.get(det.get("licenseUUID")) or {}).get("spdxID", "")
        if not spdx.startswith("cc-by"):
            skipped.append((name, spdx or "lisans bilinmiyor"))
            continue
        seq = (det.get("sequence") or "").upper()
        kept.append({
            "igem_id": det["name"],
            "slug": det["slug"],
            "title": det.get("title") or "",
            "description": det.get("description") or "",
            "role": (det.get("role") or {}).get("label") or "",
            "role_so": (det.get("role") or {}).get("accession") or "",
            "topology": (det.get("topology") or {}).get("label") or "",
            "sequence": seq,
            "length": len(seq),
            "license": spdx,
            "url": f"https://registry.igem.org/parts/{det['slug']}",
        })
        print(f"  {det['name']:<12} {len(seq):>5} bp  {spdx:<12} "
              f"{(det.get('title') or '')[:50]}")
        sys.stdout.flush()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "source": "iGEM Registry of Standard Biological Parts",
        "api": API,
        "license_note": ("Parça kayıtları Creative Commons Attribution 4.0 "
                         "(CC BY 4.0) ile kullanılır; kaynak: iGEM Foundation, "
                         "registry.igem.org."),
        "parts": kept,
    }, ensure_ascii=False, indent=1), encoding="utf-8")

    print(f"\n{len(kept)} parça yazıldı -> {OUT}")
    if missing:
        print(f"kayıtta yok ({len(missing)}): {', '.join(missing)}")
    if failed:
        print(f"ÇEKİLEMEDİ ({len(failed)}), betiği tekrar çalıştır: "
              f"{', '.join(failed)}")
    if skipped:
        print(f"lisans uygun değil: {skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
