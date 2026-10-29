"""data/ altindaki veri dosyalarinin SHA-256 ozetlerini yeniler.

Veri cekme betiklerinden birini calistirdiktan sonra bunu da
calistirin; yoksa testler ozetlerin tutmadigini bildirir.

    python tools/update_checksums.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from core import provenance  # noqa: E402

if __name__ == "__main__":
    onceki = provenance.verify()
    n = provenance.write_checksums()
    print(f"{n} dosyanin ozeti yazildi -> {provenance.CHECKSUM_FILE}")
    if onceki["degisen"]:
        print("Degisen dosyalar:")
        for y in onceki["degisen"]:
            print(f"  {y}")
    if onceki["yeni"]:
        print(f"Yeni eklenen: {len(onceki['yeni'])} dosya")
    print(provenance.describe())
