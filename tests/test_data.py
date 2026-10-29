"""Veri bütünlüğü ve denetim kuralları.

Bu testler uygulamanın *verisinin* doğruluğunu sınar. Diziler dışarıdan
çekildiği için bozulabilir ya da eksik inebilir; kayıt uzunlukları,
koordinatlar ve bilinen bölgeler burada kontrol edilir.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from core import db, igem_seed as ig, ref_vectors as rv
from core.assembly import assemble, sort_canonical, to_genbank
from core.parts_seed import CATEGORIES
from core.seq_utils import clean_seq, restriction_sites
from core.validate import validate

DATA = Path(__file__).resolve().parent.parent / "data"


@pytest.fixture(scope="module", autouse=True)
def _veritabani():
    db.init_db()


# ---------------------------------------------------------------------------
# Katalog bütünlüğü
# ---------------------------------------------------------------------------


def test_her_parcanin_uzunlugu_dizisiyle_tutar():
    for p in db.list_parts():
        if p["seq_status"] == "full":
            assert len(p["sequence"]) == p["length"], p["part_id"]


def test_parca_dizilerinde_gecersiz_karakter_yok():
    for p in db.list_parts():
        if p["seq_status"] == "full":
            assert set(p["sequence"]) <= set("ACGTN"), p["part_id"]


def test_kategoriler_taninir():
    for p in db.list_parts():
        assert p["category"] in CATEGORIES, p["part_id"]


def test_parca_kimlikleri_benzersiz():
    ids = [p["part_id"] for p in db.list_parts()]
    assert len(ids) == len(set(ids))


def test_igem_parcalari_cc_by_lisansli():
    for p in ig.igem_parts():
        assert "CC-BY" in p["accession"].upper(), p["part_id"]


def test_konaklarin_organizmasi_yazili():
    for h in db.list_hosts():
        assert h["organism"]


# ---------------------------------------------------------------------------
# Referans vektörler
# ---------------------------------------------------------------------------


def test_vektor_uzunluklari_dosyayla_tutar():
    for vid in rv.available():
        rec = rv.load(vid)
        assert len(rec["sequence"]) == rec["length"], vid


def test_vektor_dizilerinde_gecersiz_karakter_yok():
    for vid in rv.available():
        rec = rv.load(vid)
        assert set(rec["sequence"]) <= set("ACGTN"), vid


def test_her_vektorun_kaynagi_yazili():
    for vid in rv.available():
        rec = rv.load(vid)
        assert rec["accession"] and rec["source"], vid


def test_puc19_ecori_bolgesi_kanonik_konumda():
    """pUC19'un standart numaralandırmasında EcoRI 396. konumdadır.
    Bu, kaydın doğru sürüm ve doğru yönde olduğunu gösterir."""
    assert rv.region("PUC19", 396, 401) == "GAATTC"


def test_puc19_bilinen_tek_kesenler():
    seq = rv.load("PUC19")["sequence"]
    for enz in ("EcoRI", "BamHI", "SalI", "PstI", "HindIII", "KpnI", "SacI"):
        konumlar = restriction_sites(seq, [enz])[enz]
        assert len(konumlar) == 1, f"{enz}: {len(konumlar)} kesim"


def test_puc19_mcs_tanima_konumlari_yayinla_tutar():
    """Yayımlanmış pUC19 haritasındaki MCS koordinatları."""
    from core import enzymes as ez
    seq = rv.load("PUC19")["sequence"]
    beklenen = {"EcoRI": 396, "SacI": 402, "KpnI": 408, "SmaI": 412,
                "BamHI": 417, "XbaI": 423, "SalI": 429, "PstI": 435,
                "SphI": 441, "HindIII": 447}
    bulunan = ez.recognition_sites(seq, list(beklenen))
    for enz, konum in beklenen.items():
        assert bulunan[enz] == [konum], f"{enz}: {bulunan[enz]}"


def test_kesim_konumu_tanima_artı_ofset():
    """Kesim konumu = tanıma başı + üst zincir kesim ofseti.

    Bu sözleşme açıkça sınanıyor çünkü yayımlanmış haritalarla
    karşılaştırırken aradaki fark hata sanılabilir."""
    from core import enzymes as ez
    seq = rv.load("PUC19")["sequence"]
    for enz in ("EcoRI", "SacI", "KpnI", "SmaI", "BamHI", "HindIII"):
        g = ez.geometry(enz)
        tanima = ez.recognition_sites(seq, [enz])[enz][0]
        kesim = restriction_sites(seq, [enz])[enz][0]
        assert kesim == tanima + g.top_cut, enz


def test_puc19_notI_kesmez():
    seq = rv.load("PUC19")["sequence"]
    assert restriction_sites(seq, ["NotI"])["NotI"] == []


def test_dairesel_bolge_sarmasi():
    rec = rv.load("PUC19")
    L = rec["length"]
    sarmali = rv.region("PUC19", L - 2, 3)
    assert sarmali == rec["sequence"][-3:] + rec["sequence"][:3]


def test_seva_marker_hanesi_kayitla_tutarli():
    """pSEVA[marker][orijin] şeması kayıtların öznitelik tablosundan
    okundu; burada o eşleşmeyi doğruluyoruz."""
    beklenen = {"1": "bla", "2": "neo", "3": "cat",
                "4": "aadA", "5": "tetA", "6": "acc3"}
    for hane, gen in beklenen.items():
        vid = f"PSEVA{hane}11"
        if vid not in rv.available():
            continue
        etiketler = " ".join(f["label"] for f in rv.genbank_features(vid))
        assert gen in etiketler, f"{vid} içinde {gen} yok"


def test_seva_orijin_hanesi_kayitla_tutarli():
    for vid, iz in (("PSEVA211", "R6K"), ("PSEVA221", "trfA")):
        if vid not in rv.available():
            continue
        etiketler = " ".join(f["label"] for f in rv.genbank_features(vid))
        assert iz in etiketler, f"{vid} içinde {iz} yok"


def test_her_vektorun_organizmasi_var():
    for vid in rv.available():
        assert rv.load(vid).get("organizma")


# ---------------------------------------------------------------------------
# Kodon tablosu dosyası
# ---------------------------------------------------------------------------


def test_kodon_tablosu_dosyasi_tam():
    yol = DATA / "ecoli_codon_usage.json"
    assert yol.exists(), "tools/fetch_codon_usage.py çalıştırılmamış"
    d = json.loads(yol.read_text(encoding="utf-8"))
    assert len(d["codons"]) == 64
    assert d["gene_count"] > 3000
    assert len(d["reference_set"]) >= 40


# ---------------------------------------------------------------------------
# Montaj ve dışa aktarım
# ---------------------------------------------------------------------------


def test_montaj_uzunlugu_parcalarin_toplami():
    ids = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"]
    parts = sort_canonical(db.get_parts(ids))
    vec = assemble("t", parts)
    assert vec.length == sum(p["length"] for p in parts)


def test_genbank_ciktisi_bicimsel_olarak_gecerli():
    ids = ["PRO_T7LAC", "RBS_PET", "ORI_COLE1", "MRK_KAN"]
    vec = assemble("pTest", sort_canonical(db.get_parts(ids)))
    gb = to_genbank(vec)
    assert gb.startswith("LOCUS")
    assert "ORIGIN" in gb and gb.rstrip().endswith("//")
    from io import StringIO
    from Bio import SeqIO
    rec = SeqIO.read(StringIO(gb), "genbank")
    assert len(rec.seq) == vec.length


# ---------------------------------------------------------------------------
# Denetim kuralları
# ---------------------------------------------------------------------------


def _dogrula(ids, host_id, purpose="ifade"):
    vec = assemble("t", sort_canonical(db.get_parts(ids)))
    return validate(vec, host=db.get_host(host_id), purpose=purpose)


def test_orijin_yoksa_hata():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "MRK_KAN"], "BL21DE3")
    assert any(x.rule_id == "CORE-01" for x in f)


def test_marker_yoksa_hata():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "ORI_PUC"], "BL21DE3")
    assert any(x.rule_id == "CORE-02" for x in f)


def test_t7_promotoru_klonlama_susunda_uyarir():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "ORI_PUC", "MRK_KAN"], "DH5A")
    assert any(x.rule_id.startswith("PRO") or x.rule_id.startswith("HST")
               for x in f)


def test_r6k_pir_olmayan_konakta_hata():
    f = _dogrula(["PRO_J23100", "RBS_B0034", "ORI_R6K", "MRK_KAN"], "DH5A")
    assert any(x.rule_id == "ORI-05" and x.severity == "hata" for x in f)


def test_r6k_pir_konagda_sorun_yok():
    f = _dogrula(["PRO_J23100", "RBS_B0034", "ORI_R6K", "MRK_KAN"], "PIR1")
    assert not any(x.rule_id == "ORI-05" for x in f)


def test_e_coli_orijini_bacillusta_hata():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "ORI_PUC", "MRK_KAN"], "BSUB168")
    assert any(x.rule_id == "ORG-01" and x.severity == "hata" for x in f)


def test_e_coli_yapisi_e_coli_konakta_organizma_hatasi_vermez():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "ORI_PUC", "MRK_KAN"], "BL21DE3")
    assert not any(x.rule_id.startswith("ORG") for x in f)


def test_bilesik_parca_eksik_eleman_uyarisi_bastirir():
    f = _dogrula(["IGEM_J04450", "ORI_PUC", "MRK_BLA"], "BL21DE3")
    assert any(x.rule_id == "DEV-01" for x in f)
    assert not any(x.rule_id in ("CORE-03", "CORE-04") for x in f)


def test_her_bulgu_gerekli_alanlari_tasir():
    f = _dogrula(["PRO_T7LAC", "RBS_PET", "ORI_PUC", "MRK_KAN"], "DH5A")
    for x in f:
        assert x.severity in ("hata", "uyari", "bilgi")
        assert x.rule_id and x.title and x.detail


# ---------------------------------------------------------------------------
# Sürüm ve veri kimliği
# ---------------------------------------------------------------------------


def test_veri_ozetleri_tutuyor():
    """Veri dosyaları kayıtlı özetleriyle aynı mı.

    Tutmazsa ya bir dosya elle değiştirilmiştir ya da çekme betiği
    yeniden çalıştırılıp `tools/update_checksums.py` atlanmıştır.
    """
    from core import provenance
    durum = provenance.verify()
    assert not durum["degisen"], f"değişmiş: {durum['degisen']}"
    assert not durum["eksik"], f"eksik: {durum['eksik']}"
    assert not durum["yeni"], f"özeti yazılmamış: {durum['yeni']}"


def test_veri_kimligi_kararli():
    from core import provenance
    assert provenance.data_fingerprint() == provenance.data_fingerprint()
    assert len(provenance.data_fingerprint()) == 12


def test_kaydedilen_yapi_surum_izi_tasir():
    from core import provenance
    db.save_construct("test_surum_izi", "DH5A", ["ORI_PUC", "MRK_BLA"])
    kayitlar = {c["name"]: c for c in db.list_constructs()}
    assert "test_surum_izi" in kayitlar
    iz = kayitlar["test_surum_izi"].get("provenance", "")
    assert provenance.VERSION in iz
    db.delete_construct("test_surum_izi")
