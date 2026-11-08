"""PKU (terapötik canlı bakteri) vektörü: parçalar, konaklar, denetim kuralları.

Üç şey sınanır:

  1. Dizi bütünlüğü. pheP ve dapA elle yazılmadı; depodaki E. coli K-12
     CDS dosyasından kopyalandı. Test, katalogdaki dizinin o dosyadaki
     kayıtla BİREBİR aynı olduğunu ve çevirinin geçerli bir ORF verdiğini
     doğrular.
  2. Dizisi eksik parçaların (stub) dürüstçe işaretlenmesi ve
     `data/pku_parts.json` varsa tam parçaya dönüşmesi.
  3. Yeni denetim kuralları (THR, ANA, SEL, BIO). Özellikle BIO-01: dapA'yı
     plazmite koymak biyolojik içermeyi kaldırır ve uygulama bunu söylemelidir.
"""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import pytest
from Bio.Seq import Seq

from core import db, protocol
from core.assembly import assemble
from core.pku_seed import PKU_HOSTS, pku_parts
from core.validate import validate

ROOT = Path(__file__).resolve().parent.parent
CDS_CACHE = ROOT / "data" / "_cache" / "NC_000913.3_cds.fna"

KASET = ["PRO_FNRS", "RBS_B0034", "CDS_STLA", "RBS_B0034", "CDS_PHEP",
         "TER_B0015"]


@pytest.fixture(scope="module", autouse=True)
def _veritabani():
    db.init_db()


def _araç():
    """tools/fetch_pku_parts.py'yi modül olarak yükle (ağ çağrısı yapmaz)."""
    spec = importlib.util.spec_from_file_location(
        "fetch_pku_parts", ROOT / "tools" / "fetch_pku_parts.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _kur(part_ids):
    return assemble("test", db.get_parts(part_ids))


def _kurallar(part_ids, host_id, purpose="ifade"):
    vec = _kur(part_ids)
    return {f.rule_id: f for f in validate(
        vec, host=db.get_host(host_id), purpose=purpose)}


def _genom_kaydi(gen: str) -> str:
    """CDS önbelleğinden [gene=...] kaydının dizisi."""
    seq, aktif = [], False
    for satir in CDS_CACHE.read_text(encoding="utf-8").splitlines():
        if satir.startswith(">"):
            if aktif:
                break
            aktif = f"[gene={gen}]" in satir
        elif aktif:
            seq.append(satir.strip())
    return "".join(seq)


# ---------------------------------------------------------------------------
# 1. Dizi bütünlüğü
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("pid,gen,aa", [("CDS_PHEP", "pheP", 458),
                                        ("MRK_DAPA", "dapA", 292)])
def test_cds_gecerli_orf_veriyor(pid, gen, aa):
    p = db.get_part(pid)
    assert p["seq_status"] == "full"
    assert len(p["sequence"]) == p["length"]
    # cds=True: geçerli başlangıç, 3'ün katı, sonda stop, iç stop yok
    prot = str(Seq(p["sequence"]).translate(table=11, cds=True))
    assert len(prot) == aa


@pytest.mark.skipif(not CDS_CACHE.exists(), reason="CDS önbelleği yok")
@pytest.mark.parametrize("pid,gen", [("CDS_PHEP", "pheP"),
                                     ("MRK_DAPA", "dapA")])
def test_dizi_genom_kaydiyla_birebir_ayni(pid, gen):
    assert db.get_part(pid)["sequence"] == _genom_kaydi(gen)


def test_pheP_gtg_ile_baslar_ve_uygulama_bunu_taniyor():
    # Doğal pheP GTG ile başlar; ATG'ye çevirmedik. Başlangıç kodonu
    # olarak tanınmazsa ORF taraması geni kaçırır.
    from core.seq_utils import STARTS
    assert db.get_part("CDS_PHEP")["sequence"][:3] == "GTG"
    assert "GTG" in STARTS


# ---------------------------------------------------------------------------
# 2. Stub'lar ve çekilmiş dizi
# ---------------------------------------------------------------------------


def test_eksik_parcalar_isaretli_ve_uzunluklu():
    for pid in ("CDS_STLA", "PRO_FNRS", "MRK_CCDA", "PRO_CCDAB"):
        p = db.get_part(pid)
        assert p is not None, pid
        if p["seq_status"] != "full":
            assert p["length"] and p["length"] > 0
            assert "tools/fetch_pku_parts.py" in (p["notes"] + p["accession"])
            assert p["sequence"] == ""


def test_cekilmis_dizi_stub_yerine_gecer(tmp_path):
    f = tmp_path / "pku_parts.json"
    f.write_text(json.dumps({"parts": {
        "PRO_FNRS": {"sequence": "ACGT" * 25, "accession": "TEST:1..100",
                     "notes": "test"},
    }}), encoding="utf-8")
    parts = {p["part_id"]: p for p in pku_parts(f)}
    assert parts["PRO_FNRS"]["seq_status"] == "full"
    assert parts["PRO_FNRS"]["length"] == 100
    assert parts["CDS_STLA"]["seq_status"] == "stub"      # dokunulmadı


def test_bozuk_json_stub_olarak_kalir(tmp_path):
    f = tmp_path / "pku_parts.json"
    f.write_text("{bozuk", encoding="utf-8")
    assert {p["part_id"]: p["seq_status"] for p in pku_parts(f)}[
        "PRO_FNRS"] in ("stub",)


def test_init_db_tekrar_calisinca_kopya_olusturmaz():
    n1 = len([p for p in db.list_parts() if p["source"] == "pku"])
    db.init_db()
    db.init_db()
    n2 = len([p for p in db.list_parts() if p["source"] == "pku"])
    assert n1 == n2 == 6


def test_konaklar_veritabaninda_ve_genotip_dapA_delesyonunu_tasiyor():
    h = db.get_host("ECN_DAPA")
    assert h and h["organism"] == "E. coli" and h["use"] == "terapötik"
    assert "ΔdapA" in h["genotype"]
    assert db.get_host("ECN")
    assert {h["host_id"] for h in PKU_HOSTS} == {"ECN", "ECN_DAPA"}


# ---------------------------------------------------------------------------
# Araç betiğinin saf işlevleri (ağ yok)
# ---------------------------------------------------------------------------


def test_arac_cds_denetimi_bozuk_diziyi_yakalar():
    a = _araç()
    assert a.check_cds("ATG" + "GCT" * 200 + "TAA", 150, 250) == []
    assert a.check_cds("ATG" + "GCT" * 5 + "TAA" + "GCT" * 200 + "TAA",
                       150, 250)                       # iç stop
    assert a.check_cds("ATG" + "GCT" * 200, 150, 250)  # stop yok
    assert a.check_cds("CCC" + "GCT" * 200 + "TAA", 150, 250)  # başlangıç yok
    assert a.check_cds("ATG" + "GCT" * 20 + "TAA", 150, 250)   # çok kısa
    assert a.check_cds("ATGGC", 1, 9)                  # çerçeve


def test_arac_fnr_kutusu():
    a = _araç()
    assert a.fnr_box("GGGTTGATTTAAATCAAGGG")[2] == 0
    assert a.fnr_box("GGGTTGATTTAAATCCAGGG")[2] == 1    # 1 uyuşmazlık
    assert a.fnr_box("GCGCGCGCGCGCGCGCGCGC") is None


def test_arac_yukari_bolge_koordinatlari():
    a = _araç()
    # artı zincir: gen 0-bazlı 1000'de başlıyor -> 1-bazlı 901..1000
    assert a.upstream_window(1000, 1500, 100) == (901, 1000, 1)
    # eksi zincir: 5' uç 0-bazlı 1000 -> 1-bazlı 1002..1101, ters zincir
    assert a.upstream_window(1000, 500, 100) == (1002, 1101, 2)


# ---------------------------------------------------------------------------
# 3. Denetim kuralları
# ---------------------------------------------------------------------------


def test_revize_vektor_dapA_icerme_kaybini_bildirir():
    r = _kurallar(KASET + ["ORI_P15A", "MRK_DAPA"], "ECN_DAPA")
    assert "BIO-01" in r and r["BIO-01"].severity == "uyari"
    assert "SEL-01" not in r
    assert "ANA-01" in r
    # Tasarım yapısal olarak sağlam: hiçbir bulgu 'hata' olmamalı.
    assert [x for x in r.values() if x.severity == "hata"] == []


def test_dapA_markeri_normal_konakta_secilim_saglamaz():
    r = _kurallar(KASET + ["ORI_P15A", "MRK_DAPA"], "DH5A")
    assert r["SEL-01"].severity == "hata"
    assert "BIO-01" not in r


def test_dap_bagimli_konak_markersiz_tasarimda_uyarilir():
    r = _kurallar(KASET + ["ORI_P15A"], "ECN_DAPA")
    assert r["CORE-02"].severity == "hata"      # idame yöntemi yok
    assert "SEL-02" in r
    assert "BIO-01" not in r                    # içerme korunuyor


def test_terapotik_bakteride_antibiyotik_direnci_uyarilir():
    r = _kurallar(KASET + ["ORI_P15A", "MRK_KAN"], "ECN")
    assert r["BIO-02"].severity == "uyari"
    assert "BIO-03" in r
    # Aynı vektör normal bir ifade suşunda BIO-02 üretmemeli.
    assert "BIO-02" not in _kurallar(KASET + ["ORI_P15A", "MRK_KAN"], "DH5A")


def test_tasiyici_enzimsiz_ve_enzim_tasiyicisiz():
    taban = ["PRO_FNRS", "RBS_B0034", None, "TER_B0015", "ORI_P15A",
             "MRK_KAN"]
    sadece_phep = [x or "CDS_PHEP" for x in taban]
    sadece_stla = [x or "CDS_STLA" for x in taban]
    assert "THR-01" in _kurallar(sadece_phep, "ECN")
    assert "THR-02" in _kurallar(sadece_stla, "ECN")
    assert "THR-01" not in _kurallar(KASET + ["ORI_P15A", "MRK_KAN"], "ECN")


def test_anaerobik_promotor_olmadan_ana01_yok():
    r = _kurallar(["PRO_J23100", "RBS_B0034", "RPT_SFGFP", "TER_B0015",
                   "ORI_P15A", "MRK_CAT"], "DH5A")
    assert "ANA-01" not in r


# ---------------------------------------------------------------------------
# Protokol
# ---------------------------------------------------------------------------


def test_protokol_dapA_seciliminde_antibiyotik_yok_demez():
    ids = KASET + ["ORI_P15A", "MRK_DAPA"]
    p = protocol.build("pPKU", 5000, ids, [], host_name="Nissle ΔdapA")
    metin = " ".join(p.warnings)
    assert "dapA tamamlamasıyla" in metin
    assert "seçilim markerı yok" not in metin


# ---------------------------------------------------------------------------
# Toksin-antitoksin (ccdAB) idamesi
# ---------------------------------------------------------------------------

TA = KASET + ["ORI_P15A", "PRO_CCDAB", "MRK_CCDA", "MRK_CCDB"]


def test_ta_tasarimi_icermeyi_korur_ve_hatasizdir():
    r = _kurallar(TA, "ECN_DAPA")
    assert "TA-04" in r                          # nasıl çalıştığı bilgisi
    assert "BIO-01" not in r                     # içerme bozulmadı
    assert "BIO-02" not in r                     # antibiyotik yok
    assert "MRK-05" not in r                     # ccdB + ccdA birlikte
    assert "PRO-01" not in r                     # iki kaset terminatörle ayrık
    assert [x for x in r.values() if x.severity == "hata"] == []


def test_ccdB_tek_basina_hala_hata():
    r = _kurallar(KASET + ["ORI_P15A", "MRK_CCDB"], "ECN_DAPA")
    assert r["MRK-05"].severity == "hata"


def test_antitoksin_toksinsiz_uyarilir():
    r = _kurallar(KASET + ["ORI_P15A", "MRK_CCDA"], "ECN_DAPA")
    assert r["TA-01"].severity == "uyari"


def test_toksin_antitoksinden_once_gelirse_uyarilir():
    ters = KASET + ["ORI_P15A", "PRO_CCDAB", "MRK_CCDB", "MRK_CCDA"]
    assert "TA-02" in _kurallar(ters, "ECN_DAPA")
    assert "TA-02" not in _kurallar(TA, "ECN_DAPA")


def test_ccdB_dirençli_konakta_ta_calismaz():
    r = _kurallar(TA, "DB31")
    assert r["TA-03"].severity == "hata"


def test_iki_promotor_ayni_kasette_ise_pro01_sürer():
    r = _kurallar(["PRO_FNRS", "PRO_CCDAB", "RBS_B0034", "CDS_PHEP",
                   "TER_B0015", "ORI_P15A", "MRK_KAN"], "ECN")
    assert "PRO-01" in r


# ---------------------------------------------------------------------------
# Araç: GenBank okuma ve ccdAB bulma (sentetik kayıt, ağ yok)
# ---------------------------------------------------------------------------


def _sentetik_gb():
    import random
    rnd = random.Random(7)
    dolgu = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    ccda = "ATG" + "GCT" * 71 + "TAA"                      # 219 bp, 72 aa
    ccdb = "ATG" + "AAA" * 100 + "TAG"                       # 306 bp, 101 aa
    onc = dolgu(400)
    ara = dolgu(6)
    seq = onc + ccda + ara + ccdb + dolgu(200)
    a0 = len(onc) + 1
    a1 = a0 + 218
    b0 = a1 + len(ara) + 1
    b1 = b0 + 305
    gb = (
        "LOCUS       TEST\n"
        "FEATURES             Location/Qualifiers\n"
        f"     CDS             {a0}..{a1}\n"
        '                     /gene="ccdA"\n'
        f"     CDS             {b0}..{b1}\n"
        '                     /gene="ccdB"\n'
        "     CDS             join(1..5,9..12)\n"
        '                     /gene="atlanir"\n'
        "ORIGIN\n"
        + "".join(f"{i + 1:>9} {seq[i:i + 60].lower()}\n"
                  for i in range(0, len(seq), 60)) + "//\n"
    )
    return gb, seq, ccda, onc


def test_arac_genbank_okur_ve_ccdab_bulur():
    a = _araç()
    gb, seq, ccda, onc = _sentetik_gb()
    feats, okunan = a.parse_genbank(gb)
    assert okunan == seq
    assert len([f for f in feats if f["key"] == "CDS"]) == 2   # join atlandı
    r = a.find_ccdab(feats, okunan, 150)
    assert r["ccda"] == ccda
    assert r["promoter"] == onc[-150:]


def test_arac_ccdab_belirsizlikte_hata_verir():
    a = _araç()
    gb, seq, ccda, onc = _sentetik_gb()
    feats, okunan = a.parse_genbank(gb)
    with pytest.raises(RuntimeError):
        a.find_ccdab([f for f in feats if "ccdB" not in f["q"]], okunan, 150)
    with pytest.raises(RuntimeError):
        a.find_ccdab(feats, okunan, 5000)                  # yukarısı yetmiyor
