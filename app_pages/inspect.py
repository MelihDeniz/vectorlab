"""Plazmit denetle: var olan bir yapıyı yükleyip 42 kuraldan geçir."""

import streamlit as st

from core import db, enzymes as ez, inspect_plasmid as ip, ref_vectors as rv
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import circular_map, render as render_svg
from core.seq_utils import classify_cutters
from core.validate import summarize, validate
from core.i18n import t

SEV_BOX = {"hata": st.error, "uyari": st.warning, "bilgi": st.info}
SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

st.caption(
    t("Denetim motoru şimdiye kadar yalnızca parçalardan kurduğunuz "
    "yapılara uygulanabiliyordu. Burada elinizdeki bir plazmiti - "
    "hocanın verdiği, bir kayıttan indirdiğiniz ya da daha önce "
    "kurduğunuz - aynı 42 kuraldan geçirebilirsiniz.")
)

parts = db.list_parts()

kaynak = st.segmented_control(
    t("Plazmit nereden gelsin?"),
    ["Dosya yükle", "Dizi yapıştır", "Katalogdaki referans vektör"],
    default="Dosya yükle",
) or "Dosya yükle"

veri = None

if kaynak == "Dosya yükle":
    dosya = st.file_uploader(
        t("GenBank (.gb) ya da FASTA dosyası"),
        type=["gb", "gbk", "genbank", "fa", "fasta", "txt", "seq"],
        help=t("GenBank tercih edin: kaydın kendi öznitelik tablosu da "
             "kullanılır ve denetim çok daha doğru olur. SnapGene ve "
             "Benchling GenBank dışa aktarabilir."),
    )
    if dosya is not None:
        veri = ip.parse_upload(dosya.getvalue(), dosya.name)

elif kaynak == "Dizi yapıştır":
    c = st.container(horizontal=True, vertical_alignment="bottom")
    ad = c.text_input(t("Plazmit adı"), value="plazmitim", width=260)
    dairesel = c.toggle(t("Dairesel"), value=True)
    ham = st.text_area(t("Dizi (FASTA ya da düz)"), height=150,
                       placeholder=t("ATGC… ya da >adı ile başlayan FASTA"))
    if ham.strip():
        veri = ip.parse_upload(ham.encode("utf-8"), ad)
        veri["name"] = ad
        veri["circular"] = dairesel

else:
    vlist = rv.list_vectors()
    vid = st.selectbox(
        t("Referans vektör"), [v["vector_id"] for v in vlist],
        format_func=lambda i: next(
            f"{v['name']} — {v['length']} bp" for v in vlist
            if v["vector_id"] == i),
    )
    rec = rv.load(vid)
    veri = {"name": rec["name"], "sequence": rec["sequence"],
            "circular": True, "features": rec.get("features") or [],
            "source": "genbank"}

if veri is None:
    st.info(t("Bir plazmit yükleyin ya da seçin."), icon=":material/upload_file:")
    st.stop()
if veri.get("error"):
    st.error(veri["error"], icon=":material/error:")
    st.stop()
if not veri.get("sequence"):
    st.error(t("Dosyadan dizi okunamadı."), icon=":material/error:")
    st.stop()

# ---------------------------------------------------------------------------
# İnceleme
# ---------------------------------------------------------------------------

ayar = st.container(horizontal=True, vertical_alignment="bottom")
esik = ayar.slider(t("Katalog eşleşme eşiği (%)"), 85.0, 100.0, 92.0, 0.5,
                   width=260)
konak_id = ayar.selectbox(
    t("Konak suşu"), [h["host_id"] for h in db.list_hosts()],
    format_func=lambda h: (db.get_host(h) or {}).get("name", h), width=260)
amac = ayar.segmented_control(t("Amaç"), ["ifade", "klonlama"],
                              default="klonlama") or "klonlama"

ins = ip.inspect(veri, parts, min_identity=esik)
vec = ip.to_vector(ins)
bulgular = validate(vec, host=db.get_host(konak_id), purpose=amac,
                    annotation_complete=False)
sayim = summarize(bulgular)

m = st.container(horizontal=True)
m.metric(t("Plazmit"), ins.name[:18])
m.metric(t("Uzunluk"), f"{ins.length} bp")
m.metric(t("GC"), f"%{vec.gc}")
m.metric(t("Tanınan"), f"%{ins.coverage}")
m.metric(t("Hata"), sayim["hata"])
m.metric(t("Uyarı"), sayim["uyari"])

st.divider()

tab_den, tab_tani, tab_harita, tab_enz = st.tabs(
    [t("Denetim"), t("Tanınanlar"), t("Harita"), t("Enzim haritası")]
)

# ---------------------------------------------------------------------------

with tab_den:
    for n in ins.notes:
        SEV_BOX.get(n["kind"], st.info)(
            f"**{n['title']}** — {n['detail']}",
            icon=SEV_ICON.get(n["kind"], ":material/info:"))

    st.markdown(t("### Denetim bulguları"))
    if not bulgular:
        st.success(t("Bilinen bir sorun bulunamadı."),
                   icon=":material/check_circle:")
    for b in bulgular:
        with st.container(border=True):
            st.markdown(f"{SEV_ICON.get(b.severity, '')} **{b.title}**  "
                        f"`{b.rule_id}`")
            st.markdown(b.detail)
            if b.fix:
                st.markdown(f"**Ne yapmalı:** {b.fix}")
            if b.literature_query:
                st.caption(f"Literatür sorgusu: `{b.literature_query}`")

# ---------------------------------------------------------------------------

with tab_tani:
    st.markdown(t("**Denetime giren elemanlar**"))
    st.dataframe(
        [{"Tür": CATEGORY_TR.get(p["category"], p["category"]),
          "Eleman": p["name"],
          "Konum": f"{p['_span'][0]}–{p['_span'][1]}" if p.get("_span") else "—",
          "Benzerlik": (f"%{p['_identity']}" if p.get("_identity")
                        else "dosyadan"),
          "Uzunluk": p.get("length") or 0}
         for p in ins.derived_parts],
        width="stretch", hide_index=True,
        alt="Plazmitte tanınan elemanlar",
    )
    st.caption(
        t("\"dosyadan\" yazan satırlar, katalogda karşılığı bulunmayan ama "
        "kaydın kendi öznitelik tablosunda bildirilen elemanlardır. "
        "Bunlar olmasaydı denetim, var olan bir orijini yok sanabilirdi.")
    )

    if ins.unidentified:
        st.markdown(t("**Tanınmayan bölgeler**"))
        st.dataframe(
            [{"Başlangıç": g["start"], "Bitiş": g["end"],
              "Uzunluk": g["length"]} for g in ins.unidentified],
            width="stretch", hide_index=True,
            alt="Tanınmayan bölgelerin listesi",
        )
        en_buyuk = max(ins.unidentified, key=lambda g: g["length"])
        if st.button(t("En büyük bölgeyi dizi analizine gönder"),
                     icon=":material/search_insights:"):
            st.session_state.analyze_seq = ins.sequence[
                en_buyuk["start"] - 1:en_buyuk["end"]]
            st.switch_page("app_pages/analyze.py")

# ---------------------------------------------------------------------------

with tab_harita:
    tek = {k: v for k, v in
           classify_cutters(ins.sequence, ez.search("", commercial_only=True),
                            circular=ins.circular)["tek"]}
    goster = st.multiselect(t("Haritada gösterilecek tek kesen enzimler"),
                            sorted(tek), default=sorted(tek)[:8])
    render_svg(circular_map(ins.name, ins.length, ins.catalog_features,
                            {k: tek[k] for k in goster}))
    st.caption(
        t("Haritadaki etiketler katalog eşleşmelerinden gelir; dosyanın "
        "kendi öznitelikleri Tanınanlar sekmesinde listelenir.")
    )

# ---------------------------------------------------------------------------

with tab_enz:
    gruplar = classify_cutters(ins.sequence,
                               ez.search("", commercial_only=True),
                               circular=ins.circular)
    g = st.container(horizontal=True)
    g.metric(t("Tek kesen"), len(gruplar["tek"]))
    g.metric(t("İki kesen"), len(gruplar["cift"]))
    g.metric(t("Çok kesen"), len(gruplar["cok"]))
    g.metric(t("Kesmeyen"), len(gruplar["kesmez"]))
    tekler = {k: v for k, v in gruplar["tek"]}
    if tekler:
        tanimalar = ez.recognition_sites(ins.sequence, sorted(tekler),
                                         circular=ins.circular)
        st.dataframe(
            [{"Enzim": e, "Tanıma dizisi başı": tanimalar[e][0]
              if tanimalar[e] else "—", "Kesim": tekler[e][0]}
             for e in sorted(tekler)],
            width="stretch", hide_index=True,
            alt="Plazmiti bir kez kesen enzimler",
        )
        st.caption(
            t("Yayımlanmış haritalar genellikle tanıma dizisinin başını "
            "verir, uygulamalar ise kesim konumunu; ikisi de burada.")
        )
