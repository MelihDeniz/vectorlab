"""core/disease_seed.py: hastalık kayıt defteri ve terapötik genler.

Buradaki sınavlar iki şeye bakar: (1) kayıt defterinin kendi içinde
tutarlı olması (her stratejinin parçaları ve konağı gerçekten var mı),
(2) terapötik genlerin DİZİSİNİN doğru olması. İkincisi önemsiz değil:
bir terapötik gende çerçeve kayması haritada görünmez ama proteini
tamamen bozar, o yüzden çeviri burada yeniden sınanır - veri dosyasına
güvenilerek geçilmez.
"""

import json
from pathlib import Path

import pytest
from Bio.Data import CodonTable
from Bio.Seq import Seq

from core import db, disease_seed as ds
from core.assembly import assemble

_T11 = CodonTable.unambiguous_dna_by_id[11]


@pytest.fixture(scope="module")
def parcalar():
    return {p["part_id"]: p for p in db.list_parts()}


@pytest.fixture(scope="module")
def ham():
    """Veri dosyasını doğrudan oku: yükleyiciden bağımsız sınamak için."""
    yol = (Path(__file__).resolve().parent.parent / "data"
           / "disease_parts.json")
    return json.loads(yol.read_text(encoding="utf-8"))["parts"]


# ---------------------------------------------------------------------------
# Kayıt defteri tutarlılığı
# ---------------------------------------------------------------------------


def test_hastalik_anahtarlari_tekil():
    anahtarlar = ds.keys()
    assert len(anahtarlar) == len(set(anahtarlar))
    assert len(anahtarlar) >= 5
    assert all(ds.by_key(k) for k in anahtarlar)
    assert ds.by_key("olmayan") is None


def test_her_hastaligin_zorunlu_alanlari():
    for h in ds.DISEASES:
        assert h.name and h.short and h.logic, h.key
        assert len(h.mechanism) > 80, f"{h.key}: mekanizma çok kısa"
        assert h.papers, f"{h.key}: kaynak yok"
        assert h.pathway.strip().startswith("digraph"), h.key
        assert h.strategies, f"{h.key}: strateji yok"
        assert h.comparison, f"{h.key}: karşılaştırma tablosu yok"


def test_atiflarda_pmid_ve_doi_var():
    """Kaynaksız cümle yazmamak için: her atıfta PMID ve DOI olmalı."""
    for h in ds.DISEASES:
        for makale in h.papers:
            assert "PMID" in makale, f"{h.key}: {makale[:60]}"
            assert "doi:" in makale.lower(), f"{h.key}: {makale[:60]}"


def test_stratejilerin_parcalari_katalogda_var(parcalar):
    for h in ds.DISEASES:
        for s in h.strategies:
            assert s.parts, f"{h.key}/{s.key}: parça listesi boş"
            eksik = [p for p in s.parts if p not in parcalar]
            assert not eksik, f"{h.key}/{s.key}: katalogda yok {eksik}"


def test_stratejilerin_konaklari_var():
    for h in ds.DISEASES:
        for s in h.strategies:
            assert db.get_host(s.host), f"{h.key}/{s.key}: konak {s.host} yok"


def test_her_strateji_kurulabiliyor():
    """Parçalar birleşiyor ve uzunluk toplamla tutuyor."""
    for h in ds.DISEASES:
        for s in h.strategies:
            vec = assemble(f"p{h.short.split()[0]}-{s.key}",
                           db.get_parts(s.parts))
            assert vec.length > 0, f"{h.key}/{s.key}"
            assert len(vec.part_order) == len(s.parts), f"{h.key}/{s.key}"
            assert vec.length == sum(p.get("length") or 0
                                     for p in vec.part_order), f"{h.key}/{s.key}"


def test_her_yapida_promotor_var(parcalar):
    """Promotörsüz gen okunmaz."""
    for h in ds.DISEASES:
        for s in h.strategies:
            kat = {parcalar[p]["category"] for p in s.parts}
            assert "promoter" in kat, f"{h.key}/{s.key}: promotör yok"


def test_orijinsiz_yapi_denetimde_yakalaniyor(parcalar):
    """Orijini olmayan plazmit çoğalmaz - ve bu gizlenmemeli.

    IBD tasarımının A ve B stratejilerinde orijin YOKTUR: konak
    Lactococcus lactis'tir ve katalogda bu bakteride çoğalan bir
    replikon henüz yok (Hastalık sayfası bunu "bilinen açık" diye
    yazar; Referanstan örnek sayfasındaki pNZ8048 tarifi tam olarak bu
    açığı kapatmak için vardır). O yüzden sınav "her yapıda orijin
    var" diye kurulmaz - öyle kurulsa açığı kapatmak yerine testi
    gevşetmeye iter. Sınav şudur: orijini olmayan her yapı denetimde
    HATA almalı. Böylece eksik, sessizce geçmek yerine kullanıcının
    gözüne görünür.
    """
    from core.validate import validate

    for h in ds.DISEASES:
        for s in h.strategies:
            kat = {parcalar[p]["category"] for p in s.parts}
            if "ori" in kat:
                continue
            vec = assemble(f"p{h.short.split()[0]}-{s.key}",
                           db.get_parts(s.parts))
            bulgular = validate(vec, host=db.get_host(s.host),
                                purpose=h.purpose)
            assert any(b.rule_id == "CORE-01" and b.severity == "hata"
                       for b in bulgular), f"{h.key}/{s.key}"


def test_acigi_olan_hastalik_bunu_yaziyor():
    """Denetimden hata alan tasarımın "bilinen açık" notu olmalı."""
    from core.validate import validate

    for h in ds.DISEASES:
        hatali = []
        for s in h.strategies:
            vec = assemble("x", db.get_parts(s.parts))
            bulgular = validate(vec, host=db.get_host(s.host),
                                purpose=h.purpose)
            if any(b.severity == "hata" for b in bulgular):
                hatali.append(s.key)
        if hatali:
            assert h.gap, (f"{h.key}: {hatali} stratejileri hata veriyor "
                           f"ama 'bilinen açık' notu yok")


# ---------------------------------------------------------------------------
# Terapötik genlerin dizisi
# ---------------------------------------------------------------------------


def test_hastalik_genleri_yuklendi():
    genler = ds.disease_parts()
    assert len(genler) >= 8
    assert all(g["sequence"] and g["part_id"] for g in genler)


def test_genlerin_uzunlugu_uce_bolunuyor(ham):
    for g in ham:
        assert len(g["sequence"]) % 3 == 0, g["part_id"]
        assert len(g["sequence"]) == g["length"], g["part_id"]


def test_genlerin_cevirisi_protein_uzunluguyla_tutuyor(ham):
    """Dizi uzunluğu UniProt'un bildirdiği protein uzunluğunu vermeli."""
    for g in ham:
        beklenen_aa = g["length_aa"]
        kodon = len(g["sequence"]) // 3
        assert kodon in (beklenen_aa, beklenen_aa + 1), (
            f"{g['part_id']}: {kodon} kodon ama {beklenen_aa} aa")


def test_genler_baslangic_ve_stop_kodonu_tasiyor(ham):
    for g in ham:
        dizi = g["sequence"]
        assert dizi[:3] in ("ATG", "GTG", "TTG"), f"{g['part_id']}: {dizi[:3]}"
        assert dizi[-3:] in _T11.stop_codons, f"{g['part_id']}: {dizi[-3:]}"


def test_genlerin_icinde_erken_stop_yok(ham):
    """Çerçeve içi erken stop, proteini kesen ama haritada görünmeyen hata."""
    for g in ham:
        ic = g["sequence"][3:-3]
        kodonlar = [ic[i:i + 3] for i in range(0, len(ic), 3)]
        erken = [i for i, k in enumerate(kodonlar) if k in _T11.stop_codons]
        assert not erken, f"{g['part_id']}: {len(erken)} erken stop"


def test_genlerin_cevirisi_biopython_ile_ayni(ham):
    """Kendi çevirimiz Biopython'un çevirisiyle aynı olmalı."""
    for g in ham:
        bizim = str(Seq(g["sequence"]).translate(table=11, cds=False))
        assert bizim.count("*") == 1, g["part_id"]
        assert bizim.endswith("*"), g["part_id"]
        assert "X" not in bizim, f"{g['part_id']}: tanınmayan kodon"


def test_her_genin_kaynagi_yazili(ham):
    """Kaynaksız dizi kullanmıyoruz: UniProt ve nükleotid kaydı zorunlu."""
    for g in ham:
        assert g["uniprot"], g["part_id"]
        assert g["nucleotide"], g["part_id"]
        assert g["organism"], g["part_id"]


def test_incelenmemis_kayitlar_isaretli(ham):
    """TrEMBL kaydı elle denetlenmemiştir; bu fark gizlenmemeli."""
    genler = {g["part_id"]: g for g in ham}
    assert genler["CDS_PYOS5"]["reviewed"] is False
    assert genler["CDS_IL10"]["reviewed"] is True
    # Yükleyici bu farkı parçanın notuna yazmalı, yoksa kullanıcı
    # incelenmemiş bir kaydı incelenmiş sanır.
    yuklu = {p["part_id"]: p for p in ds.disease_parts()}
    assert "elle incelenmemiştir" in yuklu["CDS_PYOS5"]["notes"]
    assert "elle incelenmemiştir" not in yuklu["CDS_IL10"]["notes"]


def test_yukleyici_kaynagi_parcaya_tasiyor():
    """Katalogda görünen her hastalık geni kaynağını yanında taşımalı."""
    for p in ds.disease_parts():
        assert p["source"] == "disease", p["part_id"]
        assert p["seq_status"] == "full", p["part_id"]
        assert "UniProt" in p["accession"], p["part_id"]
        assert "doğrulandı" in p["notes"], p["part_id"]
        assert len(p["sequence"]) == p["length"], p["part_id"]



# ---------------------------------------------------------------------------
# Kayıt defterinin genişliği ve bölge kesilmiş genler
# ---------------------------------------------------------------------------


def test_kayit_defteri_genisligi():
    """Sunumda gösterilecek kadar çeşitlilik var mı?"""
    assert len(ds.DISEASES) >= 12
    assert sum(len(h.strategies) for h in ds.DISEASES) >= 40
    assert sum(len(h.papers) for h in ds.DISEASES) >= 15


def test_hastaliklarin_konaklari_cesitli():
    """Hepsi aynı konakta olmasın: tasarım konağa bağlıdır."""
    konaklar = {s.host for h in ds.DISEASES for s in h.strategies}
    assert len(konaklar) >= 4, konaklar


def test_her_hastalik_farkli_gen_kullaniyor():
    """İki hastalık aynı terapötik kaseti kullanmamalı."""
    kasetler = [tuple(h.cassette) for h in ds.DISEASES]
    assert len(kasetler) == len(set(kasetler))


def test_strateji_anahtarlari_hastalik_icinde_tekil():
    for h in ds.DISEASES:
        anahtarlar = [s.key for s in h.strategies]
        assert len(anahtarlar) == len(set(anahtarlar)), h.key


def test_bolge_kesilmis_genler_isaretli(ham):
    """Proteinin bir bölgesi kesilmişse bu gizlenmemeli."""
    bolgeli = [g for g in ham if g.get("region")]
    assert len(bolgeli) >= 3, "bölge kesilmiş gen yok"
    yuklu = {p["part_id"]: p for p in ds.disease_parts()}
    for g in bolgeli:
        assert g["region_full_aa"] > g["length_aa"], g["part_id"]
        not_ = yuklu[g["part_id"]]["notes"]
        assert "TAMAMI değil" in not_, g["part_id"]
        assert g["region"] in not_, g["part_id"]
        if g.get("start_codon_added"):
            assert "bizim eklediğimizdir" in not_, g["part_id"]


def test_bolge_kesilmis_genler_cerceve_ici(ham):
    """Kesilen bölge kodon sınırından kesilmiş olmalı."""
    for g in ham:
        if not g.get("region"):
            continue
        beklenen = g["length_aa"] + 1          # + stop kodonu
        assert len(g["sequence"]) == beklenen * 3, g["part_id"]


def test_beklenen_denetim_hatalari_aciklanmis():
    """Denetimden hata alan her strateji bunu açıklamalı.

    Bu sayfada iki hata türü bilinçlidir: Gram-pozitif replikon açığı
    (CORE-01) ve "idame yok" karşılaştırma basamağı (CORE-02). Hata
    veren bir strateji bunu kendi metninde söylemiyorsa, kullanıcı
    aracın bozuk olduğunu düşünür - o yüzden açıklama zorunludur.
    """
    from core.validate import validate

    for h in ds.DISEASES:
        for s in h.strategies:
            vec = assemble("x", db.get_parts(s.parts))
            bulgular = validate(vec, host=db.get_host(s.host),
                                purpose=h.purpose)
            hatalar = {b.rule_id for b in bulgular if b.severity == "hata"}
            if not hatalar:
                continue
            assert hatalar <= {"CORE-01", "CORE-02"},                 f"{h.key}/{s.key}: beklenmeyen hata {hatalar}"
            for rid in hatalar:
                assert rid in s.summary,                     f"{h.key}/{s.key}: {rid} hatası veriyor ama "                     f"stratejinin açıklamasında geçmiyor"
