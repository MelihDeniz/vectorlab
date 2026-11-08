"""PKU vektörünün dizisi eksik iki parçasını NCBI'dan çeker ve doğrular.

Çekilen parçalar:

  CDS_STLA  Photorhabdus luminescens PAL geni (stlA)
  PRO_FNRS  E. coli K-12 fnrS geninin yukarısındaki promotör bölgesi
  MRK_CCDA  F plazmiti (AP001918) ccdA antitoksin geni
  PRO_CCDAB ccdA'nın yukarısındaki doğal ccdAB promotör bölgesi

Dizi ezberden yazılmaz. Her parça adlı bir kayıttan gelir ve yazılmadan
önce sınanır; sınav geçmezse dosyaya hiçbir şey yazılmaz:

  stlA      çerçeve 3'ün katı, geçerli başlangıç kodonu, stop ile bitiş,
            iç stop yok, protein uzunluğu PAL aralığında (450-650 aa).
            PAL ailesinin MIO bölgesindeki Ala-Ser-Gly motifi aranır;
            bulunmazsa uyarı verilir (hata değil, motif yazımı ailelere
            göre değişir).
  PfnrS     bölgenin içinde FNR bağlanma kutusu (TTGAT-N4-ATCAA,
            en çok 2 uyuşmazlıkla) bulunmalıdır; promotör oksijensiz
            ortamda FNR ile açıldığı için bu kutunun olmaması bölgenin
            yanlış seçildiği anlamına gelir.

ccdA   AP001918 kaydında ccdB'nin (katalogdaki MRK_CCDB) hemen önündeki,
            aynı zincirdeki 219 bp'lik CDS aranır; çevirisi 72 aa olmalı.

PfnrS ve ccdAB promotörünün SINIRI bir karardır, kayıttan okunan bir gerçek değil:
varsayılan olarak fnrS geninin 5' ucunun yukarısındaki 100 bp alınır
(`--upstream` ile değişir); ccdAB için 150 bp (`--ccd-upstream`). Sonucu makalenin ek tablosundaki diziyle
karşılaştırın.

Çıktı data/pku_parts.json dosyasıdır; uygulama bir sonraki açılışta
dizisi eksik parçaları (stub) otomatik olarak tam parçayla değiştirir.

    python tools/fetch_pku_parts.py
    python tools/fetch_pku_parts.py --upstream 150
"""

from __future__ import annotations

import argparse
import datetime as _dt
import itertools
import json
import re
import sys
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HEADERS = {"User-Agent": "vectorlab/0.1 (+egitim amacli PKU vektoru)"}
OUT = ROOT / "data" / "pku_parts.json"

#: Photorhabdus luminescens subsp. laumondii TTO1 kromozomu. Betik indirdiği
#: kaydın başlığında "Photorhabdus" geçtiğini ayrıca doğrular.
PLUM_ACCESSION = "NC_005126.1"
#: E. coli K-12 MG1655 (NCBI Taxonomy)
ECOLI_TAXID = "511145"

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AAS[i]
         for i, (a, b, c) in enumerate(itertools.product(_BASES, repeat=3))}
STARTS = ("ATG", "GTG", "TTG")


# ---------------------------------------------------------------------------
# Saf yardımcılar (ağ gerektirmez, tests/test_pku.py bunları sınar)
# ---------------------------------------------------------------------------


def translate(seq: str) -> str:
    """Tablo 11 (bakteriyel). İlk kodon geçerli başlangıçsa Met okunur."""
    prot = "".join(CODON[seq[i:i + 3]] for i in range(0, len(seq) - len(seq) % 3, 3))
    if seq[:3] in STARTS:
        prot = "M" + prot[1:]
    return prot


def check_cds(seq: str, min_aa: int, max_aa: int) -> list[str]:
    """CDS'in sorunlarını döndür; boş liste = geçti."""
    sorun: list[str] = []
    if not seq or set(seq) - set("ACGT"):
        return ["dizi boş ya da ACGT dışı karakter içeriyor"]
    if len(seq) % 3:
        sorun.append(f"uzunluk 3'ün katı değil ({len(seq)})")
        return sorun
    if seq[:3] not in STARTS:
        sorun.append(f"geçerli başlangıç kodonu yok ({seq[:3]})")
    prot = translate(seq)
    if prot[-1] != "*":
        sorun.append("stop kodonuyla bitmiyor")
    if "*" in prot[:-1]:
        sorun.append(f"iç stop kodonu var (aa {prot.index('*') + 1})")
    n = len(prot) - 1
    if not (min_aa <= n <= max_aa):
        sorun.append(f"protein uzunluğu {n} aa, beklenen {min_aa}-{max_aa}")
    return sorun


def fnr_box(seq: str, max_mismatch: int = 2) -> tuple[int, str, int] | None:
    """TTGAT-N4-ATCAA kutusunu ara. (konum 0-bazlı, dizi, uyuşmazlık)"""
    model = "TTGAT" + "NNNN" + "ATCAA"
    en_iyi = None
    for i in range(len(seq) - len(model) + 1):
        pencere = seq[i:i + len(model)]
        mm = sum(1 for a, b in zip(model, pencere) if a != "N" and a != b)
        if mm <= max_mismatch and (en_iyi is None or mm < en_iyi[2]):
            en_iyi = (i, pencere, mm)
    return en_iyi


def upstream_window(chrstart: int, chrstop: int, n: int) -> tuple[int, int, int]:
    """Gen 5' ucunun yukarısındaki n bp için (seq_start, seq_stop, strand).

    NCBI Gene özetinde koordinatlar 0-bazlıdır ve eksi zincirdeki genlerde
    chrstart > chrstop olur (chrstart genin 5' ucudur). efetch ise 1-bazlı
    ve kapalı aralık ister; strand=2 eksi zincirdir.
    """
    if chrstart <= chrstop:                        # artı zincir
        return chrstart - n + 1, chrstart, 1
    return chrstart + 2, chrstart + 1 + n, 2       # eksi zincir


# ---------------------------------------------------------------------------
# Ağ
# ---------------------------------------------------------------------------


def _get(path: str, **params) -> requests.Response:
    params.setdefault("tool", "vectorlab")
    for deneme in range(4):
        r = requests.get(f"{EUTILS}/{path}", headers=HEADERS, params=params,
                         timeout=300)
        if r.status_code == 429:
            time.sleep(2 + 2 * deneme)
            continue
        r.raise_for_status()
        time.sleep(0.4)                            # NCBI hız sınırı
        return r
    raise RuntimeError("NCBI hız sınırına takıldı, biraz sonra yeniden deneyin")


def fetch_stla() -> dict:
    baslik = _get("esummary.fcgi", db="nuccore", id=PLUM_ACCESSION,
                  retmode="json").json()
    uid = baslik["result"]["uids"][0]
    title = baslik["result"][uid]["title"]
    if "photorhabdus" not in title.lower():
        raise RuntimeError(f"{PLUM_ACCESSION} beklenen organizma değil: {title}")
    print(f"kayıt: {PLUM_ACCESSION} - {title}")
    print("kodlayan diziler indiriliyor (birkaç MB)…")
    fasta = _get("efetch.fcgi", db="nuccore", id=PLUM_ACCESSION,
                 rettype="fasta_cds_na", retmode="text").text
    if not fasta.startswith(">"):
        raise RuntimeError(f"FASTA gelmedi: {fasta[:120]!r}")

    adaylar: list[tuple[str, str]] = []
    bas, parca = "", []
    for satir in fasta.splitlines() + [">"]:
        if satir.startswith(">"):
            if bas and re.search(
                    r"\[gene=stlA\]|phenylalanine ammonia[- ]lyase",
                    bas, re.I):
                adaylar.append((bas, "".join(parca).upper()))
            bas, parca = satir, []
        else:
            parca.append(satir.strip())

    if not adaylar:
        raise RuntimeError("stlA / PAL kaydı bulunamadı")
    # Yalnızca 'gene=stlA' olanları öne al; yoksa tek adayı kabul et.
    kesin = [a for a in adaylar if "[gene=stlA]" in a[0]]
    secilen = kesin if kesin else adaylar
    if len(secilen) != 1:
        for b, s in secilen:
            print("  aday:", b[:140], len(s), "bp")
        raise RuntimeError("birden fazla aday var, elle seçim gerekir")
    bas, seq = secilen[0]

    sorun = check_cds(seq, 450, 650)
    if sorun:
        raise RuntimeError("stlA doğrulaması geçmedi: " + "; ".join(sorun))
    prot = translate(seq)
    if "ASG" not in prot:
        print("UYARI: PAL ailesinin MIO motifi (Ala-Ser-Gly) bulunamadı; "
              "dizinin gerçekten PAL olduğunu elle doğrulayın.")
    loc = re.search(r"\[location=([^\]]+)\]", bas)
    tag = re.search(r"\[locus_tag=([^\]]+)\]", bas)
    return {
        "sequence": seq,
        "accession": f"{PLUM_ACCESSION}: {loc.group(1) if loc else '?'}"
                     + (f" ({tag.group(1)})" if tag else ""),
        "notes": (
            "Fenilalanini trans-sinnamata çevirir (Phe -> TCA). Makaledeki "
            "SYNB1618'de üç kopya stlA kromozoma yerleştirilmiştir.\n\n"
            f"Dizi: {PLUM_ACCESSION} kaydındaki kodlayan dizi, {len(seq)} bp, "
            f"{len(prot) - 1} aa; başlangıç kodonu, çerçeve ve iç stop "
            "yokluğu çevirerek sınandı. Başlık: " + bas[:160]),
    }


def fetch_pfnrs(upstream: int) -> dict:
    ara = _get("esearch.fcgi", db="gene", retmode="json",
               term=f"fnrS[Gene Name] AND {ECOLI_TAXID}[Taxonomy ID]").json()
    ids = ara["esearchresult"]["idlist"]
    if len(ids) != 1:
        raise RuntimeError(f"fnrS için {len(ids)} gen kaydı bulundu: {ids}")
    ozet = _get("esummary.fcgi", db="gene", id=ids[0], retmode="json").json()
    g = ozet["result"][ids[0]]
    gi = g["genomicinfo"][0]
    acc, bas, son = gi["chraccver"], int(gi["chrstart"]), int(gi["chrstop"])
    print(f"fnrS: {acc} {bas}..{son} (Gene ID {ids[0]})")

    s, e, strand = upstream_window(bas, son, upstream)
    fasta = _get("efetch.fcgi", db="nuccore", id=acc, rettype="fasta",
                 retmode="text", seq_start=s, seq_stop=e, strand=strand).text
    seq = "".join(fasta.splitlines()[1:]).upper()
    if len(seq) != upstream or set(seq) - set("ACGT"):
        raise RuntimeError(f"beklenmeyen bölge: {len(seq)} bp, {seq[:40]}")

    kutu = fnr_box(seq)
    if not kutu:
        raise RuntimeError(
            f"{upstream} bp içinde FNR kutusu (TTGAT-N4-ATCAA) bulunamadı; "
            "bölge yanlış seçilmiş olabilir. --upstream ile aralığı "
            f"değiştirin. Alınan dizi: {seq}")
    konum, dizi, mm = kutu
    print(f"FNR kutusu: {dizi} (konum {konum + 1}, {mm} uyuşmazlık)")
    yon = "artı" if strand == 1 else "eksi"
    return {
        "sequence": seq,
        "accession": f"{acc}: {s}..{e} ({yon} zincir, fnrS 5' ucunun "
                     f"yukarısı, Gene ID {ids[0]})",
        "notes": (
            "Oksijen yokken aktifleşen FNR proteini bu promotörü açar. "
            "Makalede SYNB1618'in tüm terapötik genleri bu promotörün "
            "kontrolündedir.\n\n"
            f"Dizi: fnrS geninin 5' ucunun yukarısındaki {upstream} bp. "
            "**Sınır bir karardır**, kayıttan okunan bir gerçek değildir: "
            "makalenin ek tablosundaki PfnrS dizisiyle karşılaştırın. "
            f"Bölgede FNR kutusu bulundu: {dizi} ({mm} uyuşmazlık)."),
    }


def parse_genbank(text: str) -> tuple[list[dict], str]:
    """Düz GenBank kaydından (öznitelikler, dizi). Yalnızca basit
    `a..b` ve `complement(a..b)` konumları okunur; join'ler atlanır."""
    feats: list[dict] = []
    seq: list[str] = []
    cur = None
    in_origin = False
    for line in text.splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if line.startswith("//"):
            break
        if in_origin:
            seq.append(re.sub(r"[^acgtn]", "", line.lower()))
            continue
        m = re.match(r"^     (\S+)\s+(complement\()?(\d+)\.\.(\d+)\)?\s*$", line)
        if m:
            cur = {"key": m.group(1), "strand": -1 if m.group(2) else 1,
                   "start": int(m.group(3)), "end": int(m.group(4)), "q": ""}
            feats.append(cur)
        elif re.match(r"^     \S", line):
            cur = None                       # okunamayan konumlu öznitelik
        elif cur is not None and line.startswith(" " * 21):
            cur["q"] += line.strip() + " "
    return feats, "".join(seq).upper()


def find_ccdab(feats: list[dict], seq: str, upstream: int) -> dict:
    """ccdB'yi (306 bp) bul, hemen önündeki ccdA'yı (219 bp) ve onun
    yukarısındaki promotör bölgesini döndür. Belirsizlikte hata verir."""
    cds = [f for f in feats if f["key"] == "CDS" and f["strand"] == 1]
    ccdb = [f for f in cds
            if f["end"] - f["start"] + 1 == 306
            and re.search(r"ccdB|letD|CcdB", f["q"], re.I)]
    if len(ccdb) != 1:
        raise RuntimeError(f"ccdB (306 bp) {len(ccdb)} kez bulundu")
    b = ccdb[0]
    ccda = [f for f in cds
            if f["end"] - f["start"] + 1 == 219
            and -10 <= b["start"] - f["end"] <= 100]
    if len(ccda) != 1:
        raise RuntimeError(
            f"ccdB'nin önünde 219 bp'lik {len(ccda)} CDS bulundu")
    a = ccda[0]
    a_seq = seq[a["start"] - 1:a["end"]]
    sorun = check_cds(a_seq, 72, 72)
    if sorun:
        raise RuntimeError("ccdA doğrulaması geçmedi: " + "; ".join(sorun))
    if a["start"] - 1 < upstream:
        raise RuntimeError("ccdA yukarısında yeterli dizi yok")
    pro = seq[a["start"] - 1 - upstream:a["start"] - 1]
    if len(pro) != upstream or set(pro) - set("ACGT"):
        raise RuntimeError("promotör bölgesi geçersiz")
    return {"ccda": a_seq, "ccda_loc": (a["start"], a["end"]),
            "ccdb_loc": (b["start"], b["end"]), "promoter": pro,
            "promoter_loc": (a["start"] - upstream, a["start"] - 1)}


def fetch_ccdab(upstream: int) -> dict:
    acc = "AP001918"
    print(f"{acc} (F plazmiti) indiriliyor…")
    gb = _get("efetch.fcgi", db="nuccore", id=acc, rettype="gb",
              retmode="text").text
    feats, seq = parse_genbank(gb)
    if not seq:
        raise RuntimeError("GenBank kaydı okunamadı")
    r = find_ccdab(feats, seq, upstream)
    a0, a1 = r["ccda_loc"]
    p0, p1 = r["promoter_loc"]
    print(f"ccdA: {a0}..{a1}; ccdB: {r['ccdb_loc'][0]}..{r['ccdb_loc'][1]}")
    return {
        "MRK_CCDA": {
            "sequence": r["ccda"],
            "accession": f"{acc} (F plazmiti) ccdA: {a0}..{a1}",
            "notes": (
                "CcdA, CcdB toksinine bağlanıp etkisiz kılan kararsız bir "
                "antitoksindir. ccdB ile birlikte plazmit idamesi sağlar: "
                "plazmiti kaybeden hücrede antitoksin yıkılır, toksin kalır "
                "ve hücre ölür. Katalogdaki MRK_CCDB (ccdB) ile çift olarak "
                "kullanılır.\n\n"
                f"Dizi: {acc} kaydında ccdB'nin ({r['ccdb_loc'][0]}.."
                f"{r['ccdb_loc'][1]}) hemen önündeki 219 bp'lik CDS; "
                "çevirisi 72 aa, başlangıç/stop/çerçeve sınandı."),
        },
        "PRO_CCDAB": {
            "sequence": r["promoter"],
            "accession": f"{acc} (F plazmiti): {p0}..{p1} (ccdA yukarısı)",
            "notes": (
                "ccdAB operonunun kendi promotörü. Toksin-antitoksin "
                "oranını oto-regülasyonla korur; yerine sentetik bir "
                "promotör koymak dengeyi bozabilir.\n\n"
                f"Dizi: ccdA başlangıcının yukarısındaki {upstream} bp "
                f"({p0}..{p1}). **Sınır bir karardır**, kayıttan okunan "
                "bir gerçek değildir."),
        },
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--upstream", type=int, default=100,
                    help="PfnrS için fnrS yukarısından alınacak bp (100)")
    ap.add_argument("--ccd-upstream", type=int, default=150,
                    help="ccdAB promotörü için ccdA yukarısından bp (150)")
    args = ap.parse_args()

    parts = {}
    parts["CDS_STLA"] = fetch_stla()
    parts["PRO_FNRS"] = fetch_pfnrs(args.upstream)
    parts.update(fetch_ccdab(args.ccd_upstream))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "fetched": _dt.date.today().isoformat(),
        "note": ("tools/fetch_pku_parts.py tarafından üretildi; elle "
                 "düzenlemeyin."),
        "parts": parts,
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"yazıldı: {OUT}")

    from core import provenance
    n = provenance.write_checksums()
    print(f"{n} dosyanın özeti güncellendi")
    print("Uygulamayı yeniden başlatın: parçalar otomatik olarak tam olur.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
