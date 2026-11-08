"""PKU (fenilketonüri) canlı bakteri vektörü için katalog parçaları ve konaklar.

Bu dosya iki şeyi katalogun geri kalanından AYRI tutar:

  1. Dizisi depoda zaten bulunan kaynaktan alınmış parçalar
     (`pheP`, `dapA`): dizi `data/_cache/NC_000913.3_cds.fna` içindeki
     E. coli K-12 MG1655 kayıtlarından kopyalandı, elle hiçbir baz
     yazılmadı. Her birinin başlangıç/stop kodonu ve çerçevesi
     çevrilerek sınandı (bkz. tests/test_pku.py).

  2. Dizisi henüz çekilmemiş parçalar (`stlA`/PAL, PfnrS promotörü,
     `ccdA` antitoksini ve ccdAB promotörü):
     bunlar `stub` olarak kayıtlıdır; yani haritada doğru yerde ve
     işaretli durur, ama dizi gömülü değildir. `tools/fetch_pku_parts.py`
     çalıştırılınca `data/pku_parts.json` oluşur ve bu parçalar
     otomatik olarak `full` olur. Ezberden dizi yazmak, bu aracın
     önlemeye çalıştığı sessiz hatanın ta kendisidir.

Literatür kaynağı: Isabella ve ark., Nature Biotechnology 36, 857-864
(2018), doi:10.1038/nbt.4222 (SYNB1618 suşu).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from core.parts_seed import _p

DATA = Path(__file__).resolve().parent.parent / "data" / "pku_parts.json"

PAPER = ("Isabella VM ve ark. Development of a synthetic live bacterial "
         "therapeutic for the human metabolic disease phenylketonuria. "
         "Nat Biotechnol 36, 857-864 (2018). doi:10.1038/nbt.4222")

_PHEP = (
        "GTGAAAAACGCGTCAACCGTATCGGAAGATACTGCGTCGAATCAAGAGCCGACGCTTCAT"
        "CGCGGATTACATAACCGTCATATTCAACTGATTGCGTTGGGTGGCGCAATTGGTACTGGT"
        "CTGTTTCTTGGCATTGGCCCGGCGATTCAGATGGCGGGTCCGGCTGTATTGCTGGGCTAC"
        "GGCGTCGCCGGGATCATCGCTTTCCTGATTATGCGCCAGCTTGGCGAAATGGTGGTTGAG"
        "GAGCCGGTATCCGGTTCATTTGCCCACTTTGCCTATAAATACTGGGGACCGTTTGCGGGC"
        "TTCCTCTCTGGCTGGAACTACTGGGTAATGTTCGTGCTGGTGGGAATGGCAGAGCTGACC"
        "GCTGCGGGCATCTATATGCAGTACTGGTTCCCGGATGTTCCAACGTGGATTTGGGCTGCC"
        "GCCTTCTTTATTATCATCAACGCCGTTAACCTGGTGAACGTGCGCTTATATGGCGAAACC"
        "GAGTTCTGGTTTGCGTTGATTAAAGTGCTGGCAATCATCGGTATGATCGGCTTTGGCCTG"
        "TGGCTGCTGTTTTCTGGTCACGGCGGCGAGAAAGCCAGTATCGACAACCTCTGGCGCTAC"
        "GGTGGTTTCTTCGCCACCGGCTGGAATGGGCTGATTTTGTCGCTGGCGGTAATTATGTTC"
        "TCCTTCGGCGGTCTGGAGCTGATTGGGATTACTGCCGCTGAAGCGCGCGATCCGGAAAAA"
        "AGCATTCCAAAAGCGGTAAATCAGGTGGTGTATCGCATCCTGCTGTTTTACATCGGTTCA"
        "CTGGTGGTTTTACTGGCGCTCTATCCGTGGGTGGAAGTGAAATCCAACAGTAGCCCGTTT"
        "GTGATGATTTTCCATAATCTCGACAGCAACGTGGTAGCTTCTGCGCTGAACTTCGTCATT"
        "CTGGTAGCATCGCTGTCAGTGTATAACAGCGGGGTTTACTCTAACAGCCGCATGCTGTTT"
        "GGCCTTTCTGTGCAGGGTAATGCGCCGAAGTTTTTGACTCGCGTCAGCCGTCGCGGTGTG"
        "CCGATTAACTCGCTGATGCTTTCCGGAGCGATCACTTCGCTGGTGGTGTTAATCAACTAT"
        "CTGCTGCCGCAAAAAGCGTTTGGTCTGCTGATGGCGCTGGTGGTAGCAACGCTGCTGTTG"
        "AACTGGATTATGATCTGTCTGGCGCATCTGCGTTTTCGTGCAGCGATGCGACGTCAGGGG"
        "CGTGAAACACAGTTTAAGGCGCTGCTCTATCCGTTCGGCAACTATCTCTGCATTGCCTTC"
        "CTCGGCATGATTTTGCTGCTGATGTGCACGATGGATGATATGCGCTTGTCAGCGATCCTG"
        "CTGCCGGTGTGGATTGTATTCCTGTTTATGGCATTTAAAACGCTGCGTCGGAAATAA"
    )

_DAPA = (
        "ATGTTCACGGGAAGTATTGTCGCGATTGTTACTCCGATGGATGAAAAAGGTAATGTCTGT"
        "CGGGCTAGCTTGAAAAAACTGATTGATTATCATGTCGCCAGCGGTACTTCGGCGATCGTT"
        "TCTGTTGGCACCACTGGCGAGTCCGCTACCTTAAATCATGACGAACATGCTGATGTGGTG"
        "ATGATGACGCTGGATCTGGCTGATGGGCGCATTCCGGTAATTGCCGGGACCGGCGCTAAC"
        "GCTACTGCGGAAGCCATTAGCCTGACGCAGCGCTTCAATGACAGTGGTATCGTCGGCTGC"
        "CTGACGGTAACCCCTTACTACAATCGTCCGTCGCAAGAAGGTTTGTATCAGCATTTCAAA"
        "GCCATCGCTGAGCATACTGACCTGCCGCAAATTCTGTATAATGTGCCGTCCCGTACTGGC"
        "TGCGATCTGCTCCCGGAAACGGTGGGCCGTCTGGCGAAAGTAAAAAATATTATCGGAATC"
        "AAAGAGGCAACAGGGAACTTAACGCGTGTAAACCAGATCAAAGAGCTGGTTTCAGATGAT"
        "TTTGTTCTGCTGAGCGGCGATGATGCGAGCGCGCTGGACTTCATGCAATTGGGCGGTCAT"
        "GGGGTTATTTCCGTTACGGCTAACGTCGCAGCGCGTGATATGGCCCAGATGTGCAAACTG"
        "GCAGCAGAAGGGCATTTTGCCGAGGCACGCGTTATTAATCAGCGTCTGATGCCATTACAC"
        "AACAAACTATTTGTCGAACCCAATCCAATCCCGGTGAAATGGGCATGTAAGGAACTGGGT"
        "CTTGTGGCGACCGATACGCTGCGCCTGCCAATGACACCAATCACCGACAGTGGTCGTGAG"
        "ACGGTCAGAGCGGCGCTTAAGCATGCCGGTTTGCTGTAA"
    )

#: Dizisi çekilene kadar stub olarak tutulan parçalar. `length` TAHMİNDİR;
#: dizi gelince gerçek uzunluğun yerine geçer.
_STUBS = {
    "CDS_STLA": {
        "name": "stlA / PAL (P. luminescens, fenilalanin amonyak liyaz)",
        "category": "cds",
        "length": 1572,
        "accession": "Photorhabdus luminescens subsp. laumondii TT01 "
                     "genomu BX470251 - stlA; tools/fetch_pku_parts.py ile çekilir",
        "pmid": "16079333",
        "notes": (
            "Fenilalanini trans-sinnamata çevirir (Phe -> TCA). Makaledeki "
            "SYNB1618 suşunda üç kopya stlA kromozoma yerleştirilmiştir; "
            "kaynak patent metninde Photorhabdus luminescens PAL'ı 'PAL3' "
            "diye geçer.\n\n"
            "Geni tanımlayan makale: Williams JS, Thomas M, Clarke DJ (2005) "
            "Microbiology 151:2543-2550, doi:10.1099/mic.0.28136-0 - "
            "Proteobacteria'da bildirilen ilk PAL aktivitesi.\n\n"
            "**Dizi henüz katalogda değil.** `python tools/fetch_pku_parts.py` "
            "NCBI'dan çeker ve çevirisini doğrular. Buradaki 1572 bp "
            "uzunluğu (523 aa varsayımı) bir TAHMİNDİR, doğrulanmadı; dizi "
            "gelince yerini gerçek uzunluk alır."),
    },
    "MRK_CCDA": {
        "name": "ccdA (antitoksin, F plazmiti)",
        "category": "marker",
        "length": 219,
        "strength": "toksin-antitoksin idamesi (antibiyotiksiz)",
        "requires": ["ccdB (toksin) aynı plazmitte, ccdA'dan sonra"],
        "accession": "AP001918 (F plazmiti) ccdA - tools/fetch_pku_parts.py ile çekilir",
        "notes": (
            "CcdA, CcdB toksinine bağlanıp etkisiz kılan kararsız bir "
            "antitoksindir. ccdB ile birlikte plazmit idamesi sağlar: "
            "plazmiti kaybeden hücrede antitoksin yıkılır, toksin kalır "
            "ve hücre ölür. Katalogdaki MRK_CCDB (ccdB) ile çift olarak "
            "kullanılır.\n\n"
            "**Dizi henüz katalogda değil.** `python tools/fetch_pku_parts.py` "
            "AP001918 kaydından, ccdB'nin hemen önündeki CDS'i bulup çeker "
            "ve çevirisini doğrular. 219 bp uzunluğu (72 aa) bilinen "
            "protein boyuna dayanır ama bu kayıtta doğrulanmadı."),
    },
    "PRO_CCDAB": {
        "name": "ccdAB promotörü (F plazmiti, doğal)",
        "category": "promoter",
        "length": 150,
        "strength": "oto-regüle (CcdA/CcdB kompleksi baskılar)",
        "accession": "AP001918 (F plazmiti) ccdA yukarısı - tools/fetch_pku_parts.py ile çekilir",
        "notes": (
            "ccdAB operonunun kendi promotörü. Toksin-antitoksin oranını "
            "oto-regülasyonla korur; yerine sentetik bir promotör "
            "koymak antitoksin/toksin dengesini bozabilir.\n\n"
            "**Dizi henüz katalogda değil.** Promotör sınırı bir KARARDIR: "
            "script ccdA başlangıcının yukarısındaki 150 bp'yi alır. "
            "Buradaki 150 bp uzunluğu TAHMİNDİR."),
    },
    "PRO_FNRS": {
        "name": "PfnrS (anaerobik FNR bağımlı promotör)",
        "category": "promoter",
        "length": 100,
        "strength": "oksijen yokken açık",
        "induction": "Anaerobiyoz (düşük oksijen) - FNR aktivasyonu",
        "accession": "E. coli K-12 MG1655 NC_000913.3 - fnrS geninin yukarısı; "
                     "tools/fetch_pku_parts.py ile çekilir",
        "pmid": "20070527",
        "notes": (
            "Oksijen yokken FNR proteini aktifleşir ve bu promotörü açar; "
            "bu yüzden terapötik genler yalnızca bağırsağın oksijensiz "
            "ortamında ifade edilir. Makalede SYNB1618'in tüm terapötik "
            "genleri bu promotörün kontrolündedir.\n\n"
            "PfnrS, E. coli'nin fnrS küçük RNA geninin kendi promotörüdür; "
            "sentetik bir dizi değildir. fnrS'i tanımlayan iki bağımsız "
            "makale: Durand S, Storz G (2010) Mol Microbiol 75:1215-1231, "
            "doi:10.1111/j.1365-2958.2010.07044.x ve Boysen A ve ark. (2010) "
            "J Biol Chem 285:10690-10702, doi:10.1074/jbc.M109.089755.\n\n"
            "**Dizi henüz katalogda değil.** Promotör sınırı bir KARARDIR "
            "(kayıttan okunan bir gerçek değil): script fnrS geninin "
            "yukarısındaki 100 bp'yi alır ve içinde FNR kutusu aradığını "
            "doğrular. Sınırı makalenin ek tablosundaki diziyle "
            "karşılaştırın. Buradaki 100 bp uzunluğu TAHMİNDİR."),
    },
}


def _build(overrides: dict[str, Any] | None = None) -> list[dict[str, Any]]:
    """Parça kayıtlarını kur; `overrides` çekilmiş dizileri taşır."""
    overrides = overrides or {}
    parts: list[dict[str, Any]] = [
        _p(
            part_id="CDS_PHEP",
            name="pheP (fenilalanin simporteri, E. coli K-12)",
            category="cds",
            sequence=_PHEP,
            strength="Phe alımı",
            requires=["Phe'yi dönüştüren bir enzim (PAL/stlA)"],
            notes=(
                "Yüksek afiniteli fenilalanin taşıyıcısı (PheP, "
                "H+/Phe simporteri). Makalede PAL ile birlikte ifade "
                "edilir: enzim ne kadar iyi olursa olsun, Phe hücreye "
                "yeterince girmezse hız sınırlayıcı adım taşımadır. "
                "SYNB1618'de iki kromozomal pheP kopyası vardır.\n\n"
                "Dizi: E. coli K-12 MG1655 (NC_000913.3) pheP CDS'i, "
                "1377 bp, 458 aa; depodaki CDS dosyasından kopyalandı. "
                "**GTG ile başlar**: bakterilerde başlangıçta Met okunur, "
                "denetim kuralları bunu alternatif başlama kodonu olarak "
                "tanır. Nissle 1917 ile K-12 arasında nokta farkı "
                "olabilir; protein düzeyinde ayrıca karşılaştırılmadı."),
            accession="NC_000913.3: 601959..603335 (b0576) · UniProt P24207",
            source="pku",
        ),
        _p(
            part_id="MRK_DAPA",
            name="dapA (DAP sentezi, tamamlama markeri)",
            category="marker",
            sequence=_DAPA,
            strength="besinsel seçilim (antibiyotiksiz)",
            requires=["dapA geni silinmiş (DAP bağımlı) konak"],
            notes=(
                "4-hidroksi-tetrahidrodipikolinat sentaz: bakteri hücre "
                "duvarı için gereken diaminopimelat (DAP) yolunun ilk "
                "adımı. dapA'sı silinmiş bir konak, besiyerine DAP "
                "eklenmedikçe büyüyemez; plazmit dapA'yı geri verirse "
                "plazmitsiz hücre ölür. Yani antibiyotik olmadan plazmit "
                "idamesi sağlar.\n\n"
                "**Bunun bedeli biyolojik içermedir.** Makaledeki "
                "SYNB1618'in dapA delesyonu, bakterinin DAP bulunmayan "
                "bağırsakta çoğalamamasını sağlayan güvenlik katmanıdır. "
                "dapA plazmite geri konursa bu katman kalkar: plazmiti "
                "taşıyan hücre DAP'sız da büyür. Denetim bunu BIO-01 ile "
                "bildirir.\n\n"
                "Dizi: E. coli K-12 MG1655 (NC_000913.3) dapA CDS'i, "
                "879 bp, 292 aa, kodlayan zincir; depodaki CDS "
                "dosyasından kopyalandı. Katalogdaki öteki marker CDS'leri "
                "gibi promotörsüzdür; gerçek yapıda kendi promotörü "
                "gerekir."),
            accession=("NC_000913.3: complement(2598882..2599760) (b2478) "
                       "· UniProt P0A6L2"),
            source="pku",
        ),
    ]

    for pid, spec in _STUBS.items():
        ov = overrides.get(pid)
        if ov and ov.get("sequence"):
            parts.append(_p(
                part_id=pid, name=spec["name"], category=spec["category"],
                sequence=ov["sequence"],
                strength=spec.get("strength"), induction=spec.get("induction"),
                requires=spec.get("requires", []),
                notes=ov.get("notes") or spec["notes"].split(
                    "**Dizi henüz katalogda değil.**")[0].strip(),
                accession=ov.get("accession") or spec["accession"],
                source="pku",
            ))
        else:
            parts.append(_p(
                part_id=pid, name=spec["name"], category=spec["category"],
                seq_status="stub", length=spec["length"],
                strength=spec.get("strength"), induction=spec.get("induction"),
                requires=spec.get("requires", []),
                notes=spec["notes"], accession=spec["accession"],
                source="pku",
            ))
    return parts


def pku_parts(path: Path | None = None) -> list[dict[str, Any]]:
    """Katalog parçaları; çekilmiş diziler varsa stub'ların yerini alır."""
    p = path or DATA
    overrides: dict[str, Any] = {}
    if p.exists():
        try:
            overrides = json.loads(p.read_text(encoding="utf-8")).get("parts", {})
        except (OSError, json.JSONDecodeError):
            overrides = {}
    return _build(overrides)


# ---------------------------------------------------------------------------
# Konaklar
# ---------------------------------------------------------------------------

PKU_HOSTS: list[dict[str, Any]] = [
    {
        "host_id": "ECN",
        "name": "E. coli Nissle 1917",
        "organism": "E. coli",
        "use": "terapötik",
        "t7_polymerase": False,
        "laciq": False,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "yabanıl tip (probiyotik)",
        "notes": (
            "Oral probiyotik olarak kullanılan, terapötik canlı bakteri "
            "çalışmalarının standart kasası. İnsan bağırsağını kalıcı "
            "kolonize etmez; bu yüzden tedavi kalıcı değildir, doz "
            "tekrarlanır. Klonlama suşu değildir: recA ve endA yabanıl "
            "olduğundan plazmit hazırlığını DH5-alfa gibi bir klonlama "
            "suşunda yapıp Nissle'a aktarmak daha güvenlidir. Yabanıl "
            "suşta kriptik plazmitler bildirilmiştir; orijin "
            "uyumluluğunu deneyle doğrulayın."),
    },
    {
        "host_id": "ECN_DAPA",
        "name": "E. coli Nissle 1917 ΔdapA",
        "organism": "E. coli",
        "use": "terapötik",
        "t7_polymerase": False,
        "laciq": False,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "E. coli Nissle 1917 ΔdapA (DAP bağımlı)",
        "notes": (
            "dapA geni silinmiş Nissle: besiyerine diaminopimelat (DAP) "
            "eklenmedikçe büyüyemez. Isabella ve ark. (2018) SYNB1618'de "
            "bu delesyonu biyolojik içerme için kullanır: DAP'ın "
            "bulunmadığı bağırsakta bakteri çoğalamaz. Üretimde "
            "besiyerine DAP eklenir (derişimi protokolünüzden alın)."),
    },
]
