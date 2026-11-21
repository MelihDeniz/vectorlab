"""Sipariş, maliyet ve denetim izi: tasarımdan sonra gelen iş.

Tasarım bitince iş bitmiyor. Sıradaki iki soru pratikte en çok para ve
en çok zaman kaybettiren yer:

1. **Ne satın alacağım?** Eksik sipariş edilen bir oligo bir haftaya,
   yanlış tampon bir deneye mal olur.
2. **Bu tasarımın nereden geldiğini altı ay sonra nasıl kanıtlarım?**
   Hangi veri sürümüyle, hangi denetim kurallarından geçerek üretildi?

Birincisine `core/procure.py`, ikincisine `core/audit.py` yanıt verir.

**Fiyat bu sayfada uydurulmaz.** Birim fiyat ülkeye, kuruma ve yıla
göre değişir; uydurulmuş bir fiyat yanlış bütçeye yol açar. Sayfa
TÜKETİMİ hesaplar (kaç baz oligo, kaç ünite enzim, kaç reaksiyon, kaç
iş günü) ve birim fiyatı siz girersiniz.
"""

import streamlit as st

from core import audit, db, primers as pm, procure
from core.assembly import assemble, sort_canonical
from core.i18n import t
from core.seq_utils import enzyme_info

st.caption(t(
    "Tasarımdan sonra iki soru kalır: ne satın alacağım ve bu tasarımın "
    "nereden geldiğini sonradan nasıl kanıtlarım. Bu sayfa ikisini de "
    "tek yerde yanıtlar. Fiyat uydurulmaz: tüketim hesaplanır, birim "
    "fiyatı siz girersiniz."
))

parts_by_id = {p["part_id"]: p for p in db.list_parts()}
kaset = [p for p in st.session_state.get("design_parts", [])
         if p in parts_by_id]
ad = st.session_state.get("design_name") or "pVL-01"
konak_id = st.session_state.get("design_host") or ""
konak = db.get_host(konak_id) if konak_id else None

if not kaset:
    st.info(t("Tasarımcı sayfasında henüz bir yapı yok. Bu sayfa oradaki "
              "yapıdan malzeme listesi ve denetim izi üretir."),
            icon=":material/science:")
    st.stop()

vec = assemble(ad, sort_canonical(db.get_parts(kaset)))

ozet = st.container(horizontal=True)
ozet.metric(t("Yapı"), ad[:18])
ozet.metric(t("Uzunluk"), f"{vec.length} bp")
ozet.metric(t("Konak"), (konak or {}).get("name", "—")[:18])

siparis, iz = st.tabs([t("Malzeme listesi"), t("Denetim izi")])

# ===========================================================================
# Malzeme listesi
# ===========================================================================

with siparis:
    ay = st.container(horizontal=True)
    enzimler = ay.multiselect(
        t("Kesim enzimleri"),
        sorted(enzyme_info()), default=[],
        help=t("Klonlamada kullanacağınız enzimler; ünite ihtiyacı ve "
               "satan firmalar buradan hesaplanır."))
    insert_bas = ay.number_input(t("Insert başı"), min_value=1,
                                 value=1, step=1)
    insert_son = ay.number_input(t("Insert sonu"), min_value=2,
                                 value=min(max(vec.length - 1, 2), 300),
                                 step=1)

    s2 = st.container(horizontal=True)
    dizileme = s2.checkbox(t("Dizileme yapılacak"), value=True)
    ifade = s2.checkbox(t("İfade denemesi yapılacak"), value=True)
    pcr = s2.checkbox(t("Insert PCR ile çoğaltılacak"), value=True)

    # Primerler tasarımın kendisinden türetilir: koloni PCR primerleri
    # ve katalog parçalarından çıkan dizileme çapaları. Elle primer
    # yazdırmak yerine aracın kendi hesabını kullanmak, sipariş
    # listesindeki dizinin tasarımla tutmasını garanti eder.
    primerler = []
    if vec.sequence and insert_son > insert_bas:
        try:
            k = pm.colony_pcr(vec.sequence, int(insert_bas), int(insert_son))
            primerler = [k.forward, k.reverse]
            primerler += pm.catalog_anchors(vec.sequence, db.list_parts())
        except Exception as exc:                        # pragma: no cover
            st.warning(t("Primer hesaplanamadı: ") + str(exc),
                       icon=":material/warning:")

    bom = procure.build(
        ad, part_ids=kaset, enzymes=enzimler, primers=primerler,
        host_name=(konak or {}).get("name", ""),
        needs_pcr=pcr, needs_sequencing=dizileme, needs_expression=ifade)

    m = st.container(horizontal=True)
    m.metric(t("Kalem"), len(bom.items))
    m.metric(t("Toplam oligo"), f"{bom.oligo_bases} baz")
    m.metric(t("Toplam süre"), f"{bom.total_days:g} " + t("iş günü"))

    TUR = {"oligo": t("Oligolar"), "enzim": t("Enzimler"),
           "sus": t("Suşlar"), "kit": t("Kitler"),
           "kimyasal": t("Kimyasallar")}
    for tur, baslik in TUR.items():
        kalemler = bom.of_kind(tur)
        if not kalemler:
            continue
        st.subheader(baslik)
        if tur == "oligo":
            st.dataframe(
                [{t("Ad"): i.name, t("Dizi"): i.sequence,
                  "nt": int(i.quantity), t("Ölçek"): i.scale,
                  t("Saflaştırma"): i.purification} for i in kalemler],
                hide_index=True, width="stretch")
        elif tur == "enzim":
            st.dataframe(
                [{t("Enzim"): i.name, t("Ünite"): i.quantity,
                  t("Satan firmalar (REBASE)"): i.supplier_text}
                 for i in kalemler],
                hide_index=True, width="stretch")
        else:
            st.dataframe(
                [{t("Kalem"): i.name,
                  t("Miktar"): f"{i.quantity:g} {i.unit}",
                  t("Not"): i.note} for i in kalemler],
                hide_index=True, width="stretch")

    st.subheader(t("Takvim"))
    st.dataframe(
        [{t("Aşama"): p.name, t("İş günü"): p.days, t("Not"): p.note}
         for p in bom.phases], hide_index=True, width="stretch")
    st.caption(t("Toplam adımların toplamıdır; paralelleştirme "
                 "varsayılmadı. Bir laboratuvarda iki adım aynı güne "
                 "sığabilir, ama planlarken üst sınırı bilmek işe yarar."))

    # -----------------------------------------------------------------
    # Fiyat - birim fiyat KULLANICIDAN
    # -----------------------------------------------------------------
    st.subheader(t("Maliyet"))
    st.caption(t(
        "Birim fiyatlar sizden gelir. Boş bırakılan kalem toplama "
        "KATILMAZ ve ayrı listelenir: eksik fiyatı sıfır saymak bütçeyi "
        "sessizce küçük gösterir."))
    f = st.container(horizontal=True)
    para = f.text_input(t("Para birimi"), value="TL", width=90)
    p_oligo = f.number_input(t("Oligo (baz başına)"), min_value=0.0,
                             value=0.0, step=0.1, format="%.2f")
    p_enzim = f.number_input(t("Enzim (ünite başına)"), min_value=0.0,
                             value=0.0, step=0.1, format="%.2f")
    p_sus = f.number_input(t("Yetkin hücre (dönüşüm)"), min_value=0.0,
                           value=0.0, step=1.0, format="%.2f")
    p_kit = f.number_input(t("Kit (kolon/okuma)"), min_value=0.0,
                           value=0.0, step=1.0, format="%.2f")
    p_kim = f.number_input(t("Kimyasal (şişe)"), min_value=0.0,
                           value=0.0, step=1.0, format="%.2f")

    fiyatlar = {k: v for k, v in (("oligo", p_oligo), ("enzim", p_enzim),
                                  ("sus", p_sus), ("kit", p_kit),
                                  ("kimyasal", p_kim)) if v > 0}
    if fiyatlar:
        hesap = procure.cost(bom, fiyatlar, currency=para or "TL")
        st.metric(t("Toplam"), f"{hesap['toplam']:.2f} {hesap['para_birimi']}")
        if hesap["fiyati_girilmeyen"]:
            st.warning(hesap["not"] + "\n\n"
                       + ", ".join(hesap["fiyati_girilmeyen"]),
                       icon=":material/warning:")
        else:
            st.success(hesap["not"], icon=":material/check_circle:")
        st.dataframe(hesap["satirlar"], hide_index=True, width="stretch")
    else:
        st.info(t("Birim fiyat girmediniz, bu yüzden toplam "
                  "hesaplanmadı. En az bir fiyat girin."),
                icon=":material/payments:")

    with st.expander(t("Varsayımlar")):
        for a in bom.assumptions:
            st.markdown(f"* {a}")

    d = st.container(horizontal=True)
    d.download_button(
        t("Malzeme listesi (CSV)"),
        data=procure.to_csv(bom).encode("utf-8-sig"),
        file_name=f"{ad}-malzeme.csv", mime="text/csv",
        icon=":material/download:")
    if primerler:
        d.download_button(
            t("Oligo sipariş sayfası (CSV)"),
            data=procure.oligo_order_csv(primerler).encode("utf-8-sig"),
            file_name=f"{ad}-oligo.csv", mime="text/csv",
            icon=":material/download:")
    d.download_button(
        t("Malzeme listesi (Markdown)"),
        data=procure.to_markdown(bom).encode("utf-8"),
        file_name=f"{ad}-malzeme.md", mime="text/markdown",
        icon=":material/download:")

# ===========================================================================
# Denetim izi
# ===========================================================================

with iz:
    st.caption(t(
        "Altı ay sonra ya da bir denetçi karşısında sorulacak soru şu: "
        "bu dizi hangi veriyle, hangi sürümle, hangi kurallardan "
        "geçerek üretildi ve hangilerinden GEÇEMEDİ? Bu dosya o soruya "
        "tek başına yanıt verir. Geçen kuralları değil, geçemeyenleri "
        "de yazar - yalnızca başarıları yazan bir rapor kanıt değildir."
    ))

    kayit = audit.record(ad, kaset, sequence=vec.sequence,
                         host=konak_id,
                         purpose=st.session_state.get("design_purpose",
                                                      "ifade"))
    if kayit.errors:
        st.error(kayit.verdict, icon=":material/report:")
    elif kayit.warnings:
        st.warning(kayit.verdict, icon=":material/warning:")
    else:
        st.success(kayit.verdict, icon=":material/check_circle:")

    if not kayit.data_intact:
        st.error(
            t("İzlenen veri dosyaları değişmiş: ")
            + ", ".join(kayit.data_changed)
            + t(". Bu kurulumla üretilen çıktı, depodaki veriyle aynı "
                "olmayabilir ve denetim izi bunu yazar."),
            icon=":material/warning:")

    st.dataframe(
        [{t("Alan"): t("Uygulama"), t("Değer"): kayit.app_version},
         {t("Alan"): t("Oluşturma"), t("Değer"): kayit.created},
         {t("Alan"): t("Veri parmak izi"),
          t("Değer"): kayit.data_fingerprint},
         {t("Alan"): t("Veri bütünlüğü"),
          t("Değer"): t("değişmemiş") if kayit.data_intact
          else t("DEĞİŞMİŞ")},
         {t("Alan"): t("Uzunluk"), t("Değer"): f"{kayit.length} bp"},
         {t("Alan"): t("Dizi SHA-256"), t("Değer"): kayit.sequence_sha256},
         {t("Alan"): t("Koşan kural sayısı"),
          t("Değer"): str(len(kayit.rules_run))}],
        hide_index=True, width="stretch")

    st.subheader(t("Parçaların kaynağı"))
    st.dataframe(
        [{t("Parça"): p["part_id"], t("Kategori"): p["category"],
          "bp": p["length"], t("Dizi durumu"): p["seq_status"],
          t("Kaynak kaydı"): p["accession"]} for p in kayit.parts],
        hide_index=True, width="stretch")

    st.subheader(t("Denetim bulguları"))
    if kayit.findings:
        st.dataframe(
            [{t("Önem"): f["severity"], t("Kural"): f["rule_id"],
              t("Bulgu"): f["title"]} for f in kayit.findings],
            hide_index=True, width="stretch")
    else:
        st.caption(t("Hiçbir kural bulgu üretmedi."))

    with st.expander(t("Koşan kurallar (%d)") % len(kayit.rules_run)):
        st.write(", ".join(kayit.rules_run))

    st.warning(t(
        "Bu dosya bir uygunluk belgesi DEĞİLDİR. Hiçbir düzenleyici "
        "kurum tarafından tanınmaz ve öyle sunulamaz; yaptığı iş "
        "tasarımın izlenebilir olmasıdır."
    ), icon=":material/gavel:")

    z = st.container(horizontal=True)
    z.download_button(
        t("Denetim izi (JSON)"),
        data=audit.to_json(kayit).encode("utf-8"),
        file_name=f"{ad}-denetim-izi.json", mime="application/json",
        icon=":material/download:")
    z.download_button(
        t("Denetim izi (Markdown)"),
        data=audit.to_markdown(kayit).encode("utf-8"),
        file_name=f"{ad}-denetim-izi.md", mime="text/markdown",
        icon=":material/download:")
