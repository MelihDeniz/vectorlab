"""Hastalık vektörü: literatürden bir terapötik bakteri vektörünü al,
revize et, denetle.

Sayfa beş hastalığı aynı akışla ele alır: (1) hangi makale, (2) hastalık
mekanizması, (3) makaledeki vektör neydi ve biz neyi değiştirdik,
(4) hangi idame/içerme stratejisi, (5) tasarım denetimden geçiyor mu.
Her satır `Literatür` (makaleden) ya da `Öneri` (bizim kararımız) diye
etiketlenir; ikisi karışmaz.

Hastalıkların kendisi `core/disease_seed.py` içindeki kayıt
defterindedir; bu dosya yalnızca onu gösterir.
"""

import streamlit as st

from core import db, disease_seed as ds
from core.assembly import assemble, to_fasta, to_genbank
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import linear_map, render as render_svg
from core.seq_utils import restriction_sites
from core.validate import summarize, validate

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

parts_by_id = {p["part_id"]: p for p in db.list_parts()}

# ---------------------------------------------------------------------------
# Hastalık seçimi
# ---------------------------------------------------------------------------

st.caption(
    f"Katalogda {len(ds.DISEASES)} hastalık var ve her biri farklı bir "
    "tasarım desenini gösterir: substrat yıkan enzim, azot çeken "
    "biyosentez yolu, iki enzimli döngü, salgılanan protein terapötiği "
    "ve algılayıcı-öldürücü devre. Atıfların hepsi Europe PMC'den PMID "
    "ve DOI'siyle doğrulandı."
)

secim = st.selectbox(
    "Hastalık", ds.keys(),
    format_func=lambda k: f"{ds.by_key(k).name}  —  {ds.by_key(k).logic}",
)
hastalik = ds.by_key(secim)

st.subheader(hastalik.name)
st.markdown(hastalik.mechanism)

for makale in hastalik.papers:
    st.markdown(f"**Kaynak `Literatür`:** {makale}")

st.graphviz_chart(hastalik.pathway)
st.caption(hastalik.pathway_caption)

if hastalik.biomarker:
    st.caption(f"**Etkiyi ölçen biyobelirteç:** {hastalik.biomarker}")

st.divider()

# ---------------------------------------------------------------------------
# Karşılaştırma
# ---------------------------------------------------------------------------

st.subheader("Makaledeki vektör ve bizim revizyonumuz")
st.dataframe(hastalik.comparison, width="stretch", hide_index=True,
             alt=f"{hastalik.short} için makale ve revizyon karşılaştırması")
st.caption(
    "Makale sütunundaki her satır kaynaktan alınmıştır. Revizyonun işe "
    "yarayıp yaramayacağını hiçbir tasarım aracı hesaplayamaz; bu bir "
    "**tasarım önerisidir, deneysel olarak sınanmamıştır** ve sunumda "
    "böyle söylenmelidir."
)

if hastalik.gap:
    st.info(f"**Bilinen açık:** {hastalik.gap}",
            icon=":material/lightbulb:")

st.divider()

# ---------------------------------------------------------------------------
# Strateji ve denetim
# ---------------------------------------------------------------------------

st.subheader("İdame ve içerme stratejisi")

etiketler = {f"{s.key} · {s.title}": s for s in hastalik.strategies}
secili_ad = st.radio("Strateji", list(etiketler), key=f"str_{secim}")
strateji = etiketler[secili_ad]
st.markdown(strateji.summary)

konak = db.get_host(strateji.host)
parcalar = db.get_parts(strateji.parts)      # sıra elle kurulur
vec = assemble(f"p{hastalik.short.split()[0]}-{strateji.key}", parcalar)
bulgular = validate(vec, host=konak, purpose=hastalik.purpose)
sayim = summarize(bulgular)

m = st.container(horizontal=True)
m.metric("Uzunluk", f"{vec.length} bp")
m.metric("Dizisi bilinen", f"%{vec.resolved_fraction}")
m.metric("Konak", (konak or {}).get("name", strateji.host))
m.metric("Hata", sayim["hata"])
m.metric("Uyarı", sayim["uyari"])

st.caption(
    "Operon sırası elle kurulur: Tasarımcı parçaları kategoriye göre "
    "sıralar ve iki genli bir operonu bozar. Burada parçalar "
    "stratejide yazıldığı sırayla birleştirilir."
)

tab_harita, tab_denetim, tab_parca, tab_dizi = st.tabs(
    ["Harita", "Denetim", "Parçalar", "Dizi"]
)

with tab_harita:
    tek = {k: v for k, v in
           restriction_sites(vec.sequence, None, circular=True).items()
           if len(v) == 1}
    render_svg(linear_map(vec.name, vec.length, vec.features,
                          {k: tek[k] for k in sorted(tek)[:10]}))
    st.caption(
        "Doğrusal harita, kaset düzenini dairesel haritadan daha net "
        "gösterir. Oklar okuma yönüdür."
    )

with tab_denetim:
    if not bulgular:
        st.success("Bilinen bir sorun bulunamadı.",
                   icon=":material/check_circle:")
    for b in bulgular:
        with st.container(border=True):
            st.markdown(f"{SEV_ICON.get(b.severity, '')} **{b.title}**  "
                        f"`{b.rule_id}`")
            st.markdown(b.detail)
            if b.fix:
                st.markdown(f"**Ne yapmalı:** {b.fix}")

with tab_parca:
    st.dataframe(
        [{"Sıra": i + 1, "Parça": p["name"],
          "Tür": CATEGORY_TR.get(p["category"], p["category"]),
          "Uzunluk": p.get("length") or 0,
          "Dizi": "tam" if p.get("seq_status") == "full" else "eksik",
          "Kaynak": (parts_by_id.get(p["part_id"], {})
                     .get("accession", "") or "—")[:46]}
         for i, p in enumerate(vec.part_order)],
        width="stretch", hide_index=True,
        alt=f"{hastalik.short} yapısının parça listesi",
    )
    eksik = [p for p in vec.part_order
             if p.get("seq_status") != "full"]
    if eksik:
        st.warning(
            f"{len(eksik)} parçanın dizisi gömülü değil: "
            + ", ".join(p["name"].split("(")[0].strip() for p in eksik)
            + ". Haritada doğru uzunlukta N bloğu olarak durur; dizi "
              "düzeyinde analiz o bölgede yapılamaz.",
            icon=":material/pending:")

with tab_dizi:
    c = st.container(horizontal=True)
    c.download_button("GenBank indir", to_genbank(vec),
                      file_name=f"{vec.name}.gb",
                      mime="chemical/seq-na-genbank",
                      icon=":material/download:")
    c.download_button("FASTA indir", to_fasta(vec),
                      file_name=f"{vec.name}.fasta", mime="text/plain",
                      icon=":material/download:")
    if c.button("Tasarımcıya aktar", icon=":material/construction:"):
        from core.parts_seed import CATEGORIES
        for kat in CATEGORIES:
            st.session_state[f"ms_{kat}"] = [
                pid for pid in strateji.parts
                if parts_by_id.get(pid, {}).get("category") == kat]
        st.session_state.design_name = vec.name
        st.session_state.design_host = strateji.host
        st.session_state.design_purpose = hastalik.purpose
        st.switch_page("app_pages/designer.py")
    st.code(to_fasta(vec), language="text")

st.divider()
st.caption(
    "Terapötik genlerin dizileri UniProt ve NCBI'dan çekilmiş, "
    "çevirileri UniProt'taki proteine eşit olduğu doğrulanmıştır "
    "(`tools/fetch_disease_parts.py`). Hiçbir dizi ezberden yazılmadı."
)
