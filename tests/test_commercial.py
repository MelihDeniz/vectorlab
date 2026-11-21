"""Sipariş listesi, konak taşıma, denetim izi ve komut satırı.

Bu dört parça tasarımdan SONRA gelen işi yapıyor ve hepsinin ortak bir
riski var: uydurmaya çok açık olmaları. Sipariş listesi fiyat
uydurabilir, konak taşıma "şöyle bir promotör olurdu" diyebilir,
denetim izi yalnızca başarıları yazabilir. Buradaki sınavlar tam
olarak bunları yasaklıyor:

* `test_fiyat_ve_katalog_numarasi_uydurulmuyor`
* `test_enzim_firmalari_rebase_ile_tutuyor`
* `test_tasinamayan_parca_gizlenmiyor`
* `test_hedef_parcalar_katalogda_gercekten_var`
* `test_denetim_izi_hatalari_de_yaziyor`

Geri kalanı sözleşme sınavları: çıkış kodları, CSV sütunları, tekil
oligo birleştirme gibi.
"""

from __future__ import annotations

import csv
import io
import json

import pytest

from core import audit, db, porting, procure
from core.assembly import assemble, sort_canonical
from core.seq_utils import Primer
from core.validate import RULE_GROUPS, RULE_IDS

KASET = ["PRO_T7LAC", "RBS_PET", "SIG_DSBA", "CDS_IL10", "TER_T7",
         "ORI_PUC", "MRK_KAN"]


@pytest.fixture(scope="module", autouse=True)
def _db():
    db.init_db()


def _primer(ad: str, dizi: str) -> Primer:
    return Primer(name=ad, sequence=dizi, tm=58.0, gc=50.0,
                  length=len(dizi))


# ===========================================================================
# Sipariş listesi
# ===========================================================================


def test_fiyat_ve_katalog_numarasi_uydurulmuyor():
    """Malzeme listesinde fiyat ya da katalog numarası OLMAMALI.

    İkisi de uydurmaya en açık veri türüdür: fiyat ülkeye ve yıla göre
    değişir, katalog numarası firmaya göre. Biri yanlış yazılırsa
    yanlış bütçe, öbürü yanlışsa yanlış ürün gelir. Modülün sözleşmesi
    "tüketimi hesapla, fiyatı kullanıcıdan al" - bu sınav o sözleşmeyi
    koruyor.
    """
    bom = procure.build("pTest", KASET, enzymes=["EcoRI"],
                        primers=[_primer("F", "ACGTACGTACGTACGTACGT")],
                        host_name="BL21(DE3)")
    for i in bom.items:
        assert not hasattr(i, "price"), i.name
        assert not hasattr(i, "catalog_no"), i.name
    # CSV'de fiyat sütunu BOŞ olmalı - doldurulmuş sanılmasın.
    satirlar = list(csv.DictReader(io.StringIO(procure.to_csv(bom))))
    assert satirlar, "CSV boş"
    assert "Birim fiyat" in satirlar[0]
    assert all(s["Birim fiyat"] == "" for s in satirlar)


def test_eksik_fiyat_toplama_katilmiyor():
    """Fiyatı girilmeyen kalem sıfır sayılmamalı.

    Eksik fiyatı sıfır saymak bütçeyi sessizce küçük gösterir; bütçe
    hatalarının en sık sebebi budur. Toplam bu durumda ALT SINIR'dır
    ve rapor bunu söylemek zorunda.
    """
    bom = procure.build("pTest", KASET, enzymes=["EcoRI"],
                        primers=[_primer("F", "ACGTACGTACGTACGTACGT")],
                        host_name="BL21(DE3)")
    h = procure.cost(bom, {"oligo": 1.0})          # yalnızca oligo fiyatı
    assert h["fiyati_girilmeyen"], "eksik fiyat listelenmedi"
    assert "ALT SINIR" in h["not"]
    # Toplam yalnızca oligolardan gelmeli.
    assert h["toplam"] == pytest.approx(bom.oligo_bases * 1.0)

    hepsi = procure.cost(bom, {"oligo": 1.0, "enzim": 1.0, "sus": 1.0,
                               "kit": 1.0, "kimyasal": 1.0})
    assert not hepsi["fiyati_girilmeyen"]
    assert "girildi" in hepsi["not"]


def test_enzim_firmalari_rebase_ile_tutuyor():
    """Satan firmalar REBASE'den gelmeli, ezberden değil."""
    from Bio import Restriction

    for ad in ("EcoRI", "BamHI", "BsaI"):
        bizim = procure.enzyme_suppliers(ad)
        rebase = list(getattr(Restriction, ad).supplier_list())
        assert bizim == rebase, ad
    # Olmayan enzim için boş liste, istisna değil.
    assert procure.enzyme_suppliers("OlmayanI") == []


def test_firma_sutunu_enzim_disinda_bos_kaliyor():
    """Kit ya da suş için "REBASE'de listelenmemiş" yazmak yanlış.

    REBASE bir restriksiyon enzimi veritabanıdır; miniprep kiti orada
    hiç aranmaz. Her satırda o notu basmak, aranıp bulunamadığı
    izlenimi verir.
    """
    bom = procure.build("pTest", KASET, enzymes=["EcoRI"],
                        host_name="BL21(DE3)")
    for i in bom.items:
        if i.kind == "enzim":
            assert i.supplier_text, i.name
        else:
            assert i.supplier_text == "", f"{i.kind}/{i.name}"


def test_ayni_oligo_bir_kez_siparis_ediliyor():
    """Aynı dizi iki adla geçerse tek satır olmalı.

    İki kez sipariş etmek çıplak gözle görülmeyen ama parayla ölçülen
    bir hatadır.
    """
    dizi = "ACGTACGTACGTACGTACGTAC"
    kalemler = procure.oligo_items([_primer("kPCR_F", dizi),
                                    _primer("seq_F", dizi),
                                    _primer("baska", "TTTTGGGGCCCCAAAATTTTGG")])
    assert len(kalemler) == 2, [k.name for k in kalemler]
    birlesik = [k for k in kalemler if k.sequence == dizi][0]
    assert "2 ayrı adla" in birlesik.note
    assert "kPCR_F" in birlesik.note and "seq_F" in birlesik.note

    # Sipariş sayfasında da tek satır olmalı.
    satirlar = list(csv.DictReader(io.StringIO(
        procure.oligo_order_csv([_primer("a", dizi), _primer("b", dizi)]))))
    assert len(satirlar) == 1


def test_oligo_olcegi_uzunlukla_degisiyor():
    """Uzun oligo daha büyük ölçek ve saflaştırma istemeli."""
    kisa = procure.oligo_scale(22)
    orta = procure.oligo_scale(75)
    uzun = procure.oligo_scale(120)
    cok_uzun = procure.oligo_scale(procure.SYNTHESIS_THRESHOLD + 1)
    assert kisa[0] == "25 nmol"
    assert orta[0] == "100 nmol"
    assert "HPLC" in uzun[1] or "PAGE" in uzun[1]
    # Oligo olarak sipariş edilmemeli, sentez önerilmeli.
    assert cok_uzun[0] == "sentez"
    assert "gen sentezi" in cok_uzun[2]


def test_primer_benzeri_her_girdi_kabul_ediliyor():
    """`catalog_anchors` sözlük döner, `colony_pcr` Primer döner.

    İkisini de aynı listede vermek mümkün olmalı, yoksa çağıran her
    yer kendi dönüşümünü yapar ve biri eksik yapar.
    """
    kalemler = procure.oligo_items([
        _primer("nesne", "ACGTACGTACGTACGTACGT"),
        {"name": "sozluk", "primer": "TTTTGGGGCCCCAAAATTTT"},
        "GGGGCCCCTTTTAAAAGGGG",
    ])
    assert len(kalemler) == 3
    assert all(k.sequence and k.quantity > 0 for k in kalemler)
    # Tm boş bırakılmamalı: hesaplanabiliyorsa hesaplanmalı.
    p = procure.as_primer({"name": "x", "primer": "ACGTACGTACGTACGTACGT"})
    assert p.tm is not None


def test_takvim_toplami_adimlarin_toplami():
    bom = procure.build("pTest", KASET, host_name="BL21(DE3)")
    assert bom.total_days == pytest.approx(sum(p.days for p in bom.phases))
    assert bom.total_days > 0
    # Dizileme ve ifade kapatılınca süre kısalmalı.
    kisa = procure.build("pTest", KASET, host_name="BL21(DE3)",
                         needs_sequencing=False, needs_expression=False,
                         needs_pcr=False)
    assert kisa.total_days < bom.total_days


def test_varsayimlar_yaziliyor():
    """Hesabın dayandığı varsayım gizlenmemeli."""
    bom = procure.build("pTest", KASET, enzymes=["EcoRI"],
                        host_name="BL21(DE3)")
    metin = " ".join(bom.assumptions)
    assert "Fiyat yok" in metin
    assert "katalog numarası yok" in metin
    assert "REBASE" in metin
    md = procure.to_markdown(bom)
    assert "Varsayımlar" in md and "Takvim" in md


# ===========================================================================
# Konak taşıma
# ===========================================================================


def test_hedef_parcalar_katalogda_gercekten_var():
    """Eşdeğerlik tablosundaki her parça katalogda bulunmalı.

    "Şöyle bir promotör olurdu" diye uydurulmuş bir kimlik, taşıma
    raporunu kullanılamaz kılar: kullanıcı o parçayı arar ve bulamaz.
    """
    for host_id, prof in porting.PROFILES.items():
        assert db.get_host(host_id), f"{host_id} konak değil"
        for alan in ("promoter", "signal", "marker", "ori", "rbs"):
            pid = getattr(prof, alan)
            if pid is None:
                continue
            p = db.get_part(pid)
            assert p, f"{host_id}.{alan}: {pid} katalogda yok"
            beklenen = {"promoter": "promoter", "signal": "signal",
                        "marker": "marker", "ori": "ori",
                        "rbs": "rbs"}[alan]
            assert p["category"] == beklenen, (
                f"{host_id}.{alan}: {pid} kategorisi {p['category']}")


def test_her_profilin_kategori_gerekcesi_var():
    """Değişen her parça kendi gerekçesini taşımalı.

    Konağın genel tanıtımını her satıra basmak tabloyu okunmaz
    kılıyordu; kullanıcının satırda görmek istediği şey "bu parça
    neden değişti" sorusunun yanıtı.
    """
    for host_id, prof in porting.PROFILES.items():
        for kat in ("promoter", "signal", "marker"):
            g = prof.reason_for(kat)
            assert g and g != prof.why, f"{host_id}.{kat}"
            assert len(g) > 40, f"{host_id}.{kat}: gerekçe çok kısa"


def test_tasinamayan_parca_gizlenmiyor():
    """Gram-pozitif konakta orijin taşınamaz ve bu söylenmeli.

    Sessizce en yakın orijini koymak, denetimden geçen ama
    laboratuvarda hiç çoğalmayan bir plazmit üretirdi.
    """
    r = porting.port(KASET, "LLNZ9000")
    assert not r.portable
    tikali = [s for s in r.blocked if s.category == "ori"]
    assert tikali, "orijin taşınamaz olarak işaretlenmedi"
    assert r.blockers, "ne yapılacağı yazılmadı"
    assert any("Referanstan" in b for b in r.blockers)
    # Taşınan kasette orijin HİÇ olmamalı - uydurulmuş bir orijin
    # konulmadığının kanıtı.
    kategoriler = {(db.get_part(p) or {}).get("category")
                   for p in r.ported_parts}
    assert "ori" not in kategoriler


def test_tikanan_tasimada_iskelet_adayi_veriliyor():
    """Taşınamıyorsa bari nereden alınacağı söylenmeli."""
    r = porting.port(KASET, "LLNZ9000")
    assert r.backbones, "iskelet adayı listelenmedi"
    assert "PNZ8048" in r.backbones


def test_t7_promotoru_de3_olmayan_konaga_tasinmiyor():
    """T7, polimerazı olmayan konakta çalışmaz - en sık ifade hatası."""
    r = porting.port(["PRO_T7LAC", "CDS_IL10"], "LLNZ9000")
    pro = [s for s in r.swaps if s.part_id == "PRO_T7LAC"][0]
    assert pro.action == "degisti"
    assert "T7 RNA polimeraz" in pro.reason
    # DE3 lizojeninde ise olduğu gibi kalmalı.
    r2 = porting.port(["PRO_T7LAC", "CDS_IL10"], "BL21DE3")
    pro2 = [s for s in r2.swaps if s.part_id == "PRO_T7LAC"][0]
    assert pro2.action == "kalir"


def test_ecoli_hedefinde_konaga_bagli_parca_degismiyor():
    """E. coli'den E. coli'ye taşımada parça değişmez.

    T7'siz bir kasetle sınanıyor: T7 promotörü E. coli'ye bile
    kendiliğinden taşınmaz, çünkü konağın DE3 lizojeni olması gerekir
    (bir sonraki sınav bunu ayrıca koruyor).
    """
    kaset = ["PRO_TRC", "RBS_B0034", "CDS_IL10", "TER_B0015",
             "ORI_PUC", "MRK_KAN"]
    r = porting.port(kaset, "DH5A")
    assert r.portable
    assert not r.changed
    assert r.ported_parts == kaset


def test_t7_de3_olmayan_ecoli_susuna_da_tasinmiyor():
    """DH5-alfa bir DE3 lizojeni değildir.

    Taşıma hedefi aynı tür olsa bile T7 promotörü SUŞA bağlıdır.
    Burada "aynı organizma, dolayısıyla sorun yok" diye geçmek, en
    sık yapılan ifade hatasını sessizce onaylamak olurdu - denetim
    motoru aynı hatayı PRO-03 ile yakalıyor.
    """
    r = porting.port(["PRO_T7LAC", "CDS_IL10"], "DH5A")
    assert not r.portable
    pro = [s for s in r.swaps if s.part_id == "PRO_T7LAC"][0]
    assert pro.action == "tasinamaz"
    assert "T7 RNA polimeraz" in pro.reason
    assert any("DE3" in b for b in r.blockers)


def test_rbs_dizisi_kaliyor_ama_gucu_tasinmiyor_deniyor():
    """RBS dizisi taşınır, gücü taşınmaz - bu ayrım yazılmalı."""
    r = porting.port(["RBS_B0034", "CDS_IL10"], "LLNZ9000")
    rbs = [s for s in r.swaps if s.category == "rbs"][0]
    assert rbs.action == "kalir"
    assert "gücü taşınmaz" in rbs.reason
    assert "DENEY" in rbs.reason.upper()


def test_kodon_tablosu_yoksa_cai_hesaplanmiyor():
    """Tablosu olmayan organizma için E. coli tablosuna DÜŞÜLMEMELİ.

    `codon.table_key_for_organism` bu kararı bilerek veriyor; taşıma
    modülü o kararı bozmamalı.
    """
    from core import codon

    # Katalogdaki her taşıma hedefinin tablosu var; yokmuş gibi
    # davranan bir organizmayla sınıyoruz.
    assert codon.table_key_for_organism("Saccharomyces cerevisiae") is None
    k = porting._codon_report(["CDS_IL10"], "Saccharomyces cerevisiae")
    assert k["tablo_yok"] is True
    assert "CAI hesaplanmadı" in k["not"]


def test_tasima_kodon_uyumunu_hedef_tabloyla_olcuyor():
    """İnsan geni L. lactis'e taşınınca CAI düşmeli - ve görülmeli."""
    r = porting.port(["CDS_IL10"], "LLNZ9000")
    k = r.codon_report
    assert not k["tablo_yok"]
    satir = [s for s in k["satirlar"] if s["part_id"] == "CDS_IL10"][0]
    # İnsan IL-10'u L. lactis'te E. coli'den daha kötü uyuyor.
    assert satir["cai"] < satir["cai_kaynak"]


def test_tasima_raporu_markdown_veriyor():
    md = porting.to_markdown(porting.port(KASET, "BSUB168"))
    assert "Konak taşıma" in md
    assert "Taşınamayanlar" in md
    assert "Kodon uyumu" in md


# ===========================================================================
# Denetim izi
# ===========================================================================


def test_denetim_izi_hatalari_de_yaziyor():
    """Yalnızca başarıları yazan bir rapor kanıt değildir.

    Hata veren kuralın dosyada görünmesi, uyum dosyasının bütün
    değeri. Bu yüzden hatası olan bir tasarımla sınanıyor.
    """
    # Orijini ve markeri olmayan kaset: CORE kuralları hata verir.
    r = audit.record("pEksik", ["PRO_T7LAC", "RBS_PET", "CDS_IL10"],
                     host="BL21DE3")
    assert r.errors, "hata beklenen tasarımda hata yok"
    assert "HATA" in r.verdict
    kimlikler = {f["rule_id"] for f in r.findings}
    assert kimlikler & {"CORE-01", "CORE-02"}
    # JSON ve Markdown ikisi de hatayı göstermeli.
    d = json.loads(audit.to_json(r))
    assert d["error_count"] == len(r.errors)
    assert any(f["severity"] == "hata" for f in d["findings"])
    assert "HATA" in audit.to_markdown(r)


def test_denetim_izi_uyariyi_yutmuyor():
    """Hata yoksa "geçti" denir, ama uyarı varsa o da yazılır."""
    r = audit.record("pTam", KASET, host="BL21DE3")
    if r.warnings and not r.errors:
        assert "uyarı" in r.verdict


def test_denetim_izi_veri_kimligini_tasiyor():
    """Hangi veriyle üretildiği dosyada olmalı."""
    from core import provenance

    r = audit.record("pTest", KASET, host="BL21DE3")
    assert r.data_fingerprint == provenance.data_fingerprint()
    assert r.app_version.startswith("VectorLab")
    assert len(r.sequence_sha256) == 64
    assert r.length > 0


def test_denetim_izi_parcanin_dizi_durumunu_yaziyor():
    """Parçanın dizisi tam mı kısmi mi - sonradan bilinmesi gerekir."""
    r = audit.record("pTest", KASET, host="BL21DE3")
    assert len(r.parts) == len(KASET)
    for p in r.parts:
        assert p["part_id"] and p["category"]
        assert "seq_status" in p and "accession" in p


def test_kural_katalogu_kaynakla_tutuyor():
    """`RULE_IDS`, `validate.py` içindeki kural kimlikleriyle aynı olmalı.

    Denetim izi "şu kurallar koştu" diyor; o iddianın doğru kalması
    için katalog kaynakla birlikte güncellenmek zorunda. Yeni bir
    kural eklenip katalog güncellenmezse bu sınav düşer.
    """
    import re
    from pathlib import Path

    kaynak = (Path(__file__).resolve().parent.parent / "core"
              / "validate.py").read_text(encoding="utf-8")
    # Katalogun kendi tanımını saymamak için RULE_IDS bloğunu çıkar.
    govde = kaynak.split("RULE_IDS: tuple[str, ...] = (", 1)
    assert len(govde) == 2, "RULE_IDS tanımı bulunamadı"
    kalan = govde[0] + govde[1].split(")\n", 1)[1]
    kaynaktaki = set(re.findall(r'"([A-Z]{3,4}-\d{2})"', kalan))
    assert kaynaktaki == set(RULE_IDS), {
        "kaynakta var katalogda yok": sorted(kaynaktaki - set(RULE_IDS)),
        "katalogda var kaynakta yok": sorted(set(RULE_IDS) - kaynaktaki),
    }
    # Her önekin okunur bir başlığı olmalı.
    onekler = {k.split("-")[0] for k in RULE_IDS}
    assert onekler == set(RULE_GROUPS)


def test_denetim_izi_uygunluk_belgesi_oldugunu_soylemiyor():
    """Bu dosya bir uygunluk belgesi değil ve öyle sunulmamalı."""
    r = audit.record("pTest", KASET, host="BL21DE3")
    assert r.notice_text, "kullanım amacı beyanı eklenmemiş"
    assert r.notice_text in audit.to_markdown(r)


# ===========================================================================
# Komut satırı
# ===========================================================================


def test_cli_cikis_kodu_sozlesmesi(tmp_path):
    """0 temiz, 1 hata var, 2 kullanım/okuma hatası.

    Boru hattında kullanılabilmesi bu sözleşmeye bağlı: bir CI adımı
    `denetle` çağırıp çıkış koduna bakar.
    """
    import cli
    from core import ref_vectors as rv

    temiz = tmp_path / "puc19.fasta"
    temiz.write_text(rv.to_fasta("PUC19"), encoding="utf-8")

    # Okunamayan dosya -> 2
    assert cli.main(["denetle", str(tmp_path / "yok.fasta")]) == 2
    # Boş dosya da okunamayan sayılır: boş diziyi denetlemek bütün
    # "yok" kurallarını tetikler ve sorunu tasarımdaymış gibi gösterir.
    bos = tmp_path / "bos.fasta"
    bos.write_text("", encoding="utf-8")
    assert cli.main(["denetle", str(bos)]) == 2
    # Gerçek plazmit, klonlama amacı -> hata yok
    assert cli.main(["denetle", str(temiz), "--konak", "DH5A"]) == 0
    # Yanlış konak -> kullanım hatası
    assert cli.main(["denetle", str(temiz), "--konak", "OLMAYAN"]) == 2


def test_cli_json_cikisi_gecerli(capsys, tmp_path):
    import cli
    from core import ref_vectors as rv

    yol = tmp_path / "puc19.fasta"
    yol.write_text(rv.to_fasta("PUC19"), encoding="utf-8")
    cli.main(["--json", "denetle", str(yol), "--konak", "DH5A"])
    d = json.loads(capsys.readouterr().out)
    assert d["komut"] == "denetle"
    assert d["sonuc"][0]["durum"] == "denetlendi"
    assert d["sonuc"][0]["uzunluk"] > 2000


def test_cli_tasima_tikaninca_hata_kodu_veriyor(capsys):
    """Taşıma tıkandıysa boru hattının fark etmesi gerekir."""
    import cli

    assert cli.main(["tasi", "LLNZ9000", "--parca", "ORI_PUC"]) == 1
    capsys.readouterr()
    assert cli.main(["tasi", "BL21DE3", "--parca", "PRO_T7LAC"]) == 0
    capsys.readouterr()


def test_cli_iz_hatali_tasarimda_hata_kodu_veriyor(capsys):
    import cli

    kod = cli.main(["--json", "iz", "--ad", "pEksik",
                    "--parca", "PRO_T7LAC", "--parca", "CDS_IL10"])
    d = json.loads(capsys.readouterr().out)
    assert kod == 1
    assert d["error_count"] > 0


def test_cli_oligo_csv_sutunlari(capsys, tmp_path):
    import cli
    from core import ref_vectors as rv

    yol = tmp_path / "puc19.fasta"
    yol.write_text(rv.to_fasta("PUC19"), encoding="utf-8")
    cli.main(["siparis", str(yol), "--bas", "400", "--son", "520",
              "--oligo-csv"])
    cikti = capsys.readouterr().out
    satirlar = list(csv.DictReader(io.StringIO(cikti)))
    assert satirlar, "oligo sayfası boş"
    for s in satirlar:
        assert s["Sequence"] and s["Sequence"] == s["Sequence"].upper()
        assert s["Scale"] and s["Purification"]
        assert int(s["Length"]) == len(s["Sequence"])


def test_cli_butun_komutlar_tanimli():
    """Belgelenen her alt komut gerçekten var olmalı."""
    import cli

    p = cli.build_parser()
    alt = [a for a in p._actions
           if hasattr(a, "choices") and isinstance(a.choices, dict)]
    assert alt, "alt komut bulunamadı"
    adlar = set(alt[0].choices)
    assert adlar == {"denetle", "primer", "siparis", "tasi", "iz",
                     "kimlik"}, adlar
    # Modül belgesi de aynı komutları saymalı.
    for ad in adlar:
        assert ad in cli.__doc__, f"{ad} belgede yok"
