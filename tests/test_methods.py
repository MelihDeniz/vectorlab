"""core/methods.py: numaralı tezgâh protokolleri.

Bu modülün amacı, bu işi hiç yapmamış birinin okuyup uygulayabilmesi.
Oradan iki risk doğuyor ve buradaki sınavlar o ikisini kovalıyor:

1. **Uydurma.** Protokol sayısı ve atıf, uydurmaya en açık veri
   türüdür ve yanlışı pahalıya patlar: yanlış bir inkübasyon süresi
   bir haftayı, uydurulmuş bir künye bir başvuruyu götürür. Atıfların
   PMID ve DOI taşıması şart koşuluyor.

2. **Eksik kontrol.** Yeni başlayan birinin düşeceği çukur genellikle
   protokolün kendisinde değil, atladığı kontroldedir. Kontrolsüz bir
   protokol, sonuç değil yorum üretir - bu yüzden her kaydın
   `kontrol` ve `tuzak` alanı zorunlu.
"""

from __future__ import annotations

import re

import pytest

from core import methods

PMID = re.compile(r"PMID \d{4,9}")
DOI = re.compile(r"doi:10\.\d{4,9}/\S+", re.I)


def test_kutuphane_bos_degil():
    assert len(methods.METHODS) >= 11


@pytest.mark.parametrize("key", sorted(methods.METHODS))
def test_her_yontemin_zorunlu_alanlari(key):
    m = methods.METHODS[key]
    assert m.key == key
    assert m.title and len(m.title) > 10
    assert m.amac and len(m.amac) > 40, "amaç yazılmamış"
    assert m.sure, "süre yazılmamış"


@pytest.mark.parametrize("key", sorted(methods.METHODS))
def test_her_yontemin_numarali_adimlari_var(key):
    """Kavramsal anlatım protokol değildir.

    "SDS-PAGE yap" bir adım değil, bir başlıktır. Uygulanabilir
    olması için adımların sayılabilir ve sırayla olması gerekir.

    Ölçüt adım BAŞINA uzunluk değil, ORTALAMA uzunluktur. İlk
    sürümde adım başına 40 karakter istemiştim ve sınav "570 nm
    civarında absorbansı oku." gibi KISA AMA TAM adımları da
    eledi - oysa onlarda eksik bir şey yok. Kovaladığımız şey tek
    tek kısa adımlar değil, protokolün bütününün el savuşturması.
    """
    m = methods.METHODS[key]
    assert len(m.steps) >= 6, f"{key}: {len(m.steps)} adım, çok az"
    for i, adim in enumerate(m.steps, start=1):
        assert len(adim) > 20, f"{key} adım {i}: parça cümle"
        assert adim.rstrip().endswith((".", ":", "*")), (
            f"{key} adım {i}: cümle bitmemiş")
    ortalama = sum(len(a) for a in m.steps) / len(m.steps)
    assert ortalama > 90, (
        f"{key}: adım başına ortalama {ortalama:.0f} karakter - "
        f"protokol değil başlık listesi")


@pytest.mark.parametrize("key", sorted(methods.METHODS))
def test_her_yontemin_malzemesi_yazili(key):
    """Neyin gerektiğini bilmeden protokole başlanmaz."""
    m = methods.METHODS[key]
    assert len(m.malzeme) >= 3, f"{key}: malzeme listesi eksik"


@pytest.mark.parametrize("key", sorted(methods.METHODS))
def test_her_yontemin_kontrolu_ve_tuzagi_yazili(key):
    """Kontrolsüz protokol sonuç değil yorum üretir.

    `tuzak` alanı da aynı sebeple zorunlu: yeni başlayan biri
    protokolü doğru uygulayıp yine de yanlış sonuca varabilir ve bu
    genellikle bilinen, tekrarlanan bir hatadır.
    """
    m = methods.METHODS[key]
    assert m.kontrol and len(m.kontrol) > 40, f"{key}: kontrol yazılmamış"
    assert m.tuzak and len(m.tuzak) > 40, f"{key}: tuzak yazılmamış"


@pytest.mark.parametrize("key", sorted(methods.METHODS))
def test_sayilarin_baglayici_olmadigi_yaziyor(key):
    """Protokoldeki sayı emir değil, başlangıç noktasıdır.

    Bunu söylememek, kullanıcıyı kaynağa bakmadan uygulamaya iter.
    Uyarı silinirse bu sınav düşer.
    """
    m = methods.METHODS[key]
    assert m.uyari, f"{key}: sayı uyarısı yok"
    assert "başlangıç" in m.uyari or "kaynağın" in m.uyari


def test_atiflar_pmid_ve_doi_tasiyor():
    """Kaynaksız protokol, tutmadığında bakacak yeri olmayan protokol.

    Bu deponun kuralı: künye ezberden yazılmaz. Atıfların hepsi
    Europe PMC'den PMID ve DOI'siyle alındı; sınav ikisinin de
    durduğunu şart koşuyor.
    """
    atifli = [m for m in methods.METHODS.values() if m.papers]
    assert len(atifli) >= 9, "çoğu yöntemin kaynağı yok"
    for m in atifli:
        for p in m.papers:
            assert PMID.search(p), f"{m.key}: PMID yok -> {p[:60]}"
            assert DOI.search(p), f"{m.key}: DOI yok -> {p[:60]}"


def test_atiflar_tekrarsiz_ve_tutarli():
    """Aynı makale iki yerde iki ayrı künyeyle yazılmamalı."""
    kunye: dict[str, str] = {}
    for m in methods.METHODS.values():
        for p in m.papers:
            pmid = PMID.search(p).group()
            if pmid in kunye:
                assert kunye[pmid] == p, (
                    f"{pmid} iki ayrı künyeyle yazılmış")
            kunye[pmid] = p
    assert len(kunye) >= 9


def test_protokolsuz_is_olabilir_ama_sessiz_degil():
    """Her işin numaralı protokolü olmak zorunda değil.

    Klonlama `core/protocol.py`'ye yönlendirilir, etkinlik tayini
    terapötik moleküle özgüdür. Önemli olan, protokolü olmayan işin
    bunu SÖYLEMESİ - boş bırakmak "unutulmuş" gibi durur.
    """
    from core import db, disease_seed as ds, roadmap as rm

    db.init_db()
    h = ds.by_key("ibd")
    r = rm.build("x", h.strategies[0].parts, host_id="LLNZ9000",
                 disease=h)
    gorevler = [t for p in r.phases for t in p.tasks]

    klonla = next(t for t in gorevler if t.key == "klonla")
    assert klonla.protokol is None
    assert "Laboratuvar sayfas" in klonla.nasil

    etkinlik = next(t for t in gorevler if t.key == "etkinlik")
    assert etkinlik.protokol is None
    assert "protokolü yok" in etkinlik.neden


def test_yol_haritasindaki_yontem_anahtarlari_gercek():
    """Olmayan bir protokole bağlanan görev, sessizce protokolsüz olur."""
    from core import db, disease_seed as ds, roadmap as rm

    db.init_db()
    for anahtar in ("ibd", "pku", "urat"):
        h = ds.by_key(anahtar)
        for s in h.strategies:
            r = rm.build("x", s.parts, host_id=s.host, disease=h)
            for p in r.phases:
                for t in p.tasks:
                    if t.method:
                        assert methods.get(t.method) is not None, (
                            f"{t.key} -> {t.method} diye bir protokol yok")


def test_protokoller_yol_haritasina_baglanmis():
    """Kütüphanede durup hiçbir göreve bağlanmayan protokol işe yaramaz."""
    from core import db, disease_seed as ds, roadmap as rm

    db.init_db()
    kullanilan: set[str] = set()
    for anahtar in ds.keys():
        h = ds.by_key(anahtar)
        for s in h.strategies:
            r = rm.build("x", s.parts, host_id=s.host, disease=h)
            for p in r.phases:
                for t in p.tasks:
                    if t.method:
                        kullanilan.add(t.method)
    bagsiz = set(methods.METHODS) - kullanilan
    assert not bagsiz, f"hiçbir göreve bağlanmayan protokol: {bagsiz}"


def test_markdown_ciktisi_tam():
    for m in methods.METHODS.values():
        md = methods.to_markdown(m)
        assert m.title in md
        for adim in m.steps:
            assert adim in md
        assert m.kontrol in md and m.tuzak in md
        assert m.uyari in md
        # Adımlar numaralı basılmalı.
        assert "\n1. " in md and f"\n{len(m.steps)}. " in md


def test_yol_haritasi_markdowninda_protokoller_var():
    """Yazdırıp tezgâha götüren kişi protokolü dosyada da bulmalı."""
    from core import db, disease_seed as ds, roadmap as rm

    db.init_db()
    h = ds.by_key("ibd")
    r = rm.build("x", h.strategies[0].parts, host_id="LLNZ9000",
                 disease=h)
    md = rm.to_markdown(r)
    for p in r.phases:
        for t in p.tasks:
            if t.protokol is None:
                continue
            for adim in t.protokol.steps:
                assert adim in md, f"{t.method}: adım dosyada yok"
            assert t.protokol.kontrol in md


def test_elektroporasyon_yalnizca_gram_pozitifte():
    """E. coli'ye glisinli elektroporasyon protokolü gereksiz."""
    from core import db, roadmap as rm

    db.init_db()
    gp = rm.build("x", ["PRO_T7LAC", "CDS_IL10", "MRK_THYA_LL"],
                  host_id="LLNZ9000")
    ec = rm.build("x", ["PRO_T7LAC", "CDS_IL10", "ORI_PUC", "MRK_KAN"],
                  host_id="BL21DE3")
    gp_yontem = {t.method for p in gp.phases for t in p.tasks}
    ec_yontem = {t.method for p in ec.phases for t in p.tasks}
    assert "elektroporasyon" in gp_yontem
    assert "elektroporasyon" not in ec_yontem
