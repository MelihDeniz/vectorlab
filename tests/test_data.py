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


#: Beş konak organizması için beklenen tablolar ve genomları.
#: Erişim numaraları NCBI'ya sorularak teyit edildi.
KODON_TABLOLARI = {
    "ecoli": ("NC_000913.3", 3000),
    "bsubtilis": ("NC_000964.3", 3000),
    "pputida": ("NC_002947.4", 4000),
    "llactis": ("NC_009004.1", 2000),
    "cglutamicum": ("NC_003450.3", 2500),
}


@pytest.mark.parametrize("anahtar", sorted(KODON_TABLOLARI))
def test_kodon_tablosu_dosyasi_tam(anahtar):
    accession, en_az_gen = KODON_TABLOLARI[anahtar]
    yol = DATA / "codon_usage" / f"{anahtar}.json"
    assert yol.exists(), "tools/fetch_codon_usage.py çalıştırılmamış"
    d = json.loads(yol.read_text(encoding="utf-8"))
    assert len(d["codons"]) == 64
    assert d["accession"] == accession
    assert d["gene_count"] >= en_az_gen
    assert len(d["reference_set"]) >= 40
    # Referans kümenin nasıl seçildiği kayıtta durmalı: tablo sonradan
    # denetlenecekse ölçütün ne olduğu bilinmek zorunda.
    assert d["reference_rule"]
    assert sum(d["reference_classes"].values()) == len(d["reference_set"]) or \
        sum(d["reference_classes"].values()) >= 40

    # "İyi bir CAI kaçtır" ölçeği ÖLÇÜLMÜŞ olmalı, gömülü sayı değil.
    olcek = d["reference_cai"]
    assert olcek["gene_count"] >= 40
    assert 0.50 < olcek["median"] < 0.90, (
        f"{anahtar} referans CAI ortancası {olcek['median']} - çok ifade "
        "edilen genlerin 0,7 civarında olması beklenir")
    assert olcek["p10"] < olcek["median"] < olcek["p90"]


@pytest.mark.parametrize("anahtar,beklenen", [
    ("ecoli", "TAA"), ("bsubtilis", "TAA"), ("llactis", "TAA"),
    ("pputida", "TGA"), ("cglutamicum", "TAA"),
])
def test_baskin_durdurma_kodonu_dogru(anahtar, beklenen):
    """Durdurma kodonu dağılımı gerçek olmalı.

    Bu testin varlık sebebi gerçek bir hatadır: durdurma kodonları
    `count_codons` tarafından hiç sayılmıyordu (ilk ve son kodon
    atlanır), ama tabloda sıfırdan farklı sayılar görünüyordu - onlar
    iç çerçevede durdurma kodonu taşıyan üç selenoprotein geninden
    geliyordu. Sonuç: E. coli tablosu "TGA, durdurma kodonlarının
    %100'ü" diyordu, oysa E. coli'de baskın durdurma kodonu TAA'dır.

    Beklenen değerler yayımlanmış örüntüdür: AT'den zengin genomlarda
    TAA baskındır, GC'den zengin P. putida'da TGA'ya kayar.
    """
    d = json.loads((DATA / "codon_usage" / f"{anahtar}.json")
                   .read_text(encoding="utf-8"))
    stoplar = {c: v["count"] for c, v in d["codons"].items()
               if v["aa"] == "*"}
    assert sum(stoplar.values()) > 1000, (
        "durdurma kodonları gen sayısı kadar sayılmamış")
    assert max(stoplar, key=lambda c: stoplar[c]) == beklenen, stoplar
    # Durdurma kodonu anlamlı kodon istatistiğine girmemeli.
    for c, v in d["codons"].items():
        if v["aa"] == "*":
            assert v["per_1000"] == 0.0 and v["w"] == 0.0


def test_kodon_tablolari_birbirinden_farkli():
    """Beş tablo gerçekten farklı olmalı.

    Aynı çıkarlarsa bütün bu iş boşunadır - tek tablo yeterdi. Fark
    organizmaların GC içeriğinden gelir ve w değerlerinde görünür.
    """
    tablolar = {}
    for anahtar in KODON_TABLOLARI:
        d = json.loads((DATA / "codon_usage" / f"{anahtar}.json")
                       .read_text(encoding="utf-8"))
        tablolar[anahtar] = tuple(sorted(
            (c, v["w"]) for c, v in d["codons"].items()))
    assert len(set(tablolar.values())) == len(tablolar), (
        "iki organizmanın tablosu birebir aynı çıktı")


def test_gc_zengini_genomlar_gc_zengini_kodon_seciyor():
    """Bağımsız biyolojik denetim.

    P. putida (%61 GC) ve C. glutamicum (%54 GC) üçüncü konumda G ya da
    C tercih etmeli; B. subtilis (%44) ve L. cremoris (%35) etmemeli.
    Tablolar bunu kendiliğinden göstermiyorsa hesap yanlıştır.
    """
    def gc3_agirligi(anahtar: str) -> float:
        d = json.loads((DATA / "codon_usage" / f"{anahtar}.json")
                       .read_text(encoding="utf-8"))
        # En iyi kodonların (w = 1,0) kaçı G/C ile bitiyor?
        en_iyi = [c for c, v in d["codons"].items()
                  if v["w"] >= 1.0 and v["aa"] != "*"]
        gc_biten = [c for c in en_iyi if c[2] in "GC"]
        return len(gc_biten) / len(en_iyi)

    gc_zengini = [gc3_agirligi("pputida"), gc3_agirligi("cglutamicum")]
    at_zengini = [gc3_agirligi("bsubtilis"), gc3_agirligi("llactis")]
    assert min(gc_zengini) > max(at_zengini), (
        f"GC zengini {gc_zengini} AT zengini {at_zengini} - "
        "tablolar genomun GC içeriğini izlemiyor")


# ---------------------------------------------------------------------------
# Montaj ve dışa aktarım
# ---------------------------------------------------------------------------


def test_montaj_uzunlugu_parcalarin_toplami():
    ids = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"]
    parts = sort_canonical(db.get_parts(ids))
    vec = assemble("t", parts)
    assert vec.length == sum(p["length"] for p in parts)


def test_genbank_satirlari_80_sutunu_gecmiyor():
    """GenBank başlığı 80 sütunda sarılır.

    Biopython uzun satırı da okur, ama SnapGene ve Benchling biçimde
    daha katıdır. Kullanım amacı beyanı COMMENT alanına eklendiğinde
    satır 150 karaktere çıkıyordu; `_comment_lines` bunun için var.
    """
    ids = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"]
    vec = assemble("pTest", sort_canonical(db.get_parts(ids)))
    baslik = to_genbank(vec).split("ORIGIN")[0]
    uzun = [ln for ln in baslik.splitlines() if len(ln) > 80]
    assert not uzun, f"80 sütunu geçen satır(lar): {uzun}"


def test_disa_aktarilan_dosyalar_kullanim_beyanini_tasiyor():
    """Beyan dosyayla birlikte gitmeli.

    Dosya uygulamadan çıkıp başka birine ulaştığında, hangi amaçla
    üretildiğini kendi içinde söylemesi gerekir.
    """
    from core import notice
    from core.assembly import to_fasta

    ids = ["PRO_T7LAC", "RBS_PET", "ORI_COLE1", "MRK_KAN"]
    vec = assemble("pTest", sort_canonical(db.get_parts(ids)))

    # GenBank COMMENT alanını Biopython'un kendisi okuyabilmeli.
    from io import StringIO

    from Bio import SeqIO
    rec = SeqIO.read(StringIO(to_genbank(vec)), "genbank")
    yorum = rec.annotations.get("comment", "").replace("\n", " ")
    assert "arastirma ve egitim" in yorum
    assert "tani ya da tedavi" in yorum

    assert notice.EXPORT in to_fasta(vec).splitlines()[0]


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


# ---------------------------------------------------------------------------
# Vektörlerden çıkarılmış parçalar
# ---------------------------------------------------------------------------


def test_cikarilan_parcalar_kaynak_vektorde_birebir_bulunur():
    """Her çıkarılan parça, kaydettiği aksesyon ve koordinatta
    gerçekten o dizide olmalı. Bu, parçanın kanıtını sınar."""
    from core import extracted_seed as xs
    from core.seq_utils import revcomp
    aksesyon_to_vid = {rv.load(v)["accession"]: v for v in rv.available()}
    denetlenen = 0
    for p in xs.extracted_parts():
        ak, _, aralik = p["accession"].partition(":")
        vid = aksesyon_to_vid.get(ak)
        if not vid:
            continue
        ters = "ters zincir" in aralik
        aralik = aralik.split(" ")[0]
        bas, _, son = aralik.partition("..")
        dizi = rv.region(vid, int(bas), int(son))
        if ters:
            dizi = revcomp(dizi)
        assert dizi == p["sequence"], p["part_id"]
        denetlenen += 1
    assert denetlenen > 30, f"yalnızca {denetlenen} parça denetlendi"


def test_cikarilan_parcalar_katalogda():
    parcalar = {p["part_id"] for p in db.list_parts()
                if p["source"] == "vector"}
    assert len(parcalar) > 40


def test_seva_koleksiyonu_tam():
    seva = [v for v in rv.available() if v.startswith("PSEVA")]
    assert len(seva) >= 95, f"yalnızca {len(seva)} SEVA vektörü"


def test_seva_hane_tablolari_kayitla_tutar():
    """Marker, orijin ve kargo haneleri kayıtların öznitelik
    tablolarından türetildi; örneklerle doğrula."""
    beklenen = [("PSEVA232", "lacZa_19"), ("PSEVA246", "luxC"),
                ("PSEVA237R", "mcherry"), ("PSEVA528", "xylS"),
                ("PSEVA414", "lacIq")]
    for vid, iz in beklenen:
        if vid not in rv.available():
            continue
        etiketler = " ".join(f["label"] for f in rv.genbank_features(vid))
        assert iz in etiketler, f"{vid} içinde {iz} yok"


def test_seva_suzgeci_calisir():
    kanamisin = rv.filter_vectors(seva_marker="2")
    assert kanamisin
    for vid in kanamisin:
        assert rv.load(vid)["seva_marker"] == "2"
    assert len(rv.filter_vectors(seva_cargo="7R")) == 3


# ---------------------------------------------------------------------------
# Sinyal peptitleri
# ---------------------------------------------------------------------------


def test_sinyal_peptidi_cevirisi_uniprot_ile_tutar():
    """Nükleotid dizisinin çevirisi, UniProt'un bildirdiği peptide
    eşit olmalı. Bir baz kayması burada yakalanır."""
    import json as _json
    from core.signal_seed import signal_parts
    yol = DATA / "signal_peptides.json"
    assert yol.exists(), "tools/fetch_signal_peptides.py çalıştırılmamış"
    ham = {p["part_id"]: p
           for p in _json.loads(yol.read_text(encoding="utf-8"))["peptides"]}
    from Bio.Data import CodonTable
    t11 = CodonTable.unambiguous_dna_by_id[11]
    for p in signal_parts():
        kayit = ham[p["part_id"]]
        dna = p["sequence"]
        assert len(dna) == len(kayit["peptide"]) * 3
        cevrilen = []
        for i in range(0, len(dna), 3):
            k = dna[i:i + 3]
            cevrilen.append("M" if i == 0 and k in ("ATG", "GTG", "TTG")
                            else t11.forward_table.get(k, "X"))
        assert "".join(cevrilen) == kayit["peptide"], p["part_id"]


def test_sinyal_peptitleri_katalogda():
    sinyaller = db.list_parts("signal")
    assert len(sinyaller) >= 6
    for s in sinyaller:
        assert s["accession"].startswith("UniProt")


def test_sinyal_montaj_sirasinda_etiketten_once():
    """Sinyal peptidi N ucunda olmak zorunda; etiket araya girerse
    taşıma makinesi sinyali tanımaz."""
    from core.assembly import CANONICAL_ORDER
    assert CANONICAL_ORDER.index("signal") < CANONICAL_ORDER.index("tag")
    assert CANONICAL_ORDER.index("rbs") < CANONICAL_ORDER.index("signal")


def _sig(ids, host_id, purpose="ifade"):
    from core.assembly import assemble, sort_canonical
    vec = assemble("t", sort_canonical(db.get_parts(ids)))
    return validate(vec, host=db.get_host(host_id), purpose=purpose)


_TEMEL = ["PRO_T7LAC", "RBS_PET", "TER_T7", "ORI_COLE1", "MRK_KAN"]


def test_iki_sinyal_hata_verir():
    f = _sig(_TEMEL + ["SIG_PELB", "SIG_OMPA"], "BL21DE3")
    assert any(x.rule_id == "SIG-01" and x.severity == "hata" for x in f)


def test_gram_uyusmazligi_uyarir():
    f = _sig(_TEMEL + ["SIG_SACB"], "BL21DE3")
    assert any(x.rule_id == "SIG-02" for x in f)
    f2 = _sig(_TEMEL + ["SIG_PELB"], "BSUB168")
    assert any(x.rule_id == "SIG-02" for x in f2)


def test_uyumlu_sinyal_uyari_vermez():
    f = _sig(_TEMEL + ["SIG_PELB", "TAG_HIS6"], "BL21DE3")
    assert not any(x.rule_id.startswith("SIG") for x in f)


def test_tat_yolu_bilgilendirir():
    f = _sig(_TEMEL + ["SIG_TORA"], "BL21DE3")
    assert any(x.rule_id == "SIG-03" for x in f)


def test_sinyal_ve_gfp_uyarir():
    f = _sig(_TEMEL + ["SIG_PELB", "RPT_SFGFP"], "BL21DE3")
    assert any(x.rule_id == "SIG-04" for x in f)
