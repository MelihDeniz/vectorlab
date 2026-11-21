"""Konak taşıma sayfası: kaseti başka bir bakteriye uyarla.

Laboratuvarda sürekli sorulan ama hiçbir dizi aracının yanıtlamadığı
soru: *"E. coli'de kurduğum bu kaseti Lactococcus lactis'e taşıyabilir
miyim?"* Dizi araçları diziyi gösterir, konağı bilmez. Oysa bir
kasetin konağa bağlı parçaları var ve taşınırken değişmesi gerekir.

Sayfa üç şeyi ayrı ayrı gösterir ve üçüncüsünü gizlemez: neyin
değiştiği, neyin aynı kaldığı, neyin TAŞINAMADIĞI. Hesap
`core/porting.py` içindedir.
"""

import streamlit as st

from core import db, porting
from core.assembly import assemble, sort_canonical
from core.i18n import t
from core.validate import summarize, validate

st.caption(t(
    "Bir kaseti başka bir bakteriye taşımak diziyi kopyalamak değildir: "
    "promotör, sinyal peptidi, marker ve orijin konağa bağlıdır. Bu "
    "sayfa hangisinin değişmesi gerektiğini söyler ve değişemeyeni "
    "değişmiş gibi göstermez."
))

parts_by_id = {p["part_id"]: p for p in db.list_parts()}

# ---------------------------------------------------------------------------
# Kaynak kaset
# ---------------------------------------------------------------------------

kaynak = st.segmented_control(
    t("Kaset nereden gelsin?"),
    ["Tasarımcı'daki yapı", "Parça seç"],
    default="Tasarımcı'daki yapı",
    format_func=t) or "Tasarımcı'daki yapı"

if kaynak == "Tasarımcı'daki yapı":
    kaset = [p for p in st.session_state.get("design_parts", [])
             if p in parts_by_id]
    if not kaset:
        st.info(t("Tasarımcı sayfasında henüz bir yapı yok. Oradan parça "
                  "seçin ya da yukarıdan 'Parça seç'e geçin."),
                icon=":material/science:")
        st.stop()
    st.caption(t("Kaynak kaset: ") + " + ".join(kaset))
else:
    kaset = st.multiselect(
        t("Kasetin parçaları"), sorted(parts_by_id),
        format_func=lambda p: f"{p} — {parts_by_id[p]['name']}",
        default=["PRO_T7LAC", "RBS_PET", "SIG_DSBA", "TER_T7",
                 "ORI_PUC", "MRK_KAN"])
    if not kaset:
        st.info(t("Taşınacak kaset için en az bir parça seçin."),
                icon=":material/science:")
        st.stop()

# ---------------------------------------------------------------------------
# Hedef konak
# ---------------------------------------------------------------------------

konaklar = sorted(h["host_id"] for h in db.list_hosts())
# Taşıma profili tanımlı konaklar listenin başına alınır: E. coli'ye
# taşımak hiçbir şeyi değiştirmez ve kullanıcı bunu sayfanın asıl
# işi sanabilir.
profilli = [h for h in konaklar if h in porting.PROFILES]
digerleri = [h for h in konaklar if h not in porting.PROFILES]

hedef = st.selectbox(
    t("Hedef konak"), profilli + digerleri,
    format_func=lambda h: (
        f"{(db.get_host(h) or {}).get('name', h)}"
        + ("" if h in porting.PROFILES else t("  (taşıma profili yok)"))),
)

rapor = porting.port(kaset, hedef)

# ---------------------------------------------------------------------------
# Sonuç
# ---------------------------------------------------------------------------

o = st.container(horizontal=True)
o.metric(t("Hedef organizma"), rapor.target_organism or "—")
o.metric(t("Değişen parça"), len(rapor.changed))
o.metric(t("Taşınamayan"), len(rapor.blocked))

if rapor.portable:
    st.success(t("Kaset parça parça taşındı; taşınamayan parça yok."),
               icon=":material/check_circle:")
else:
    st.error(
        t("Taşıma tamamlanmadı: ") + str(len(rapor.blocked))
        + t(" parça taşınamadı. Bu bir arıza değil, kataloğun bilinen "
            "sınırı - aşağıda ne yapılacağı yazılı."),
        icon=":material/report:")

prof = porting.profile(hedef)
if prof:
    st.info(f"**{prof.organism}** (Gram {prof.gram}) — {prof.why}",
            icon=":material/biotech:")
else:
    st.info(t("Bu konak için taşıma profili yok: hedef E. coli olduğunda "
              "konağa bağlı parçalar değişmez, yalnızca T7 gereksinimi "
              "denetlenir."), icon=":material/info:")

st.subheader(t("Parça parça"))
st.dataframe(
    [{t("Parça"): s.part_id,
      t("Kategori"): s.category,
      "": s.symbol,
      t("Yerine"): s.replacement or "—",
      t("Gerekçe"): s.reason} for s in rapor.swaps],
    hide_index=True, width="stretch")

st.code(" + ".join(rapor.ported_parts), language=None)

if rapor.blockers:
    st.subheader(t("Taşınamayanlar ve ne yapmalı"))
    for b in rapor.blockers:
        st.warning(b, icon=":material/build:")

if rapor.backbones:
    with st.expander(
            t("%s içinde çoğalan referans vektörler (%d)")
            % (rapor.target_organism, len(rapor.backbones))):
        st.caption(t(
            "İskeleti parçadan kurmak yerine bütün olarak bu "
            "vektörlerden birinden alın; Referanstan kes-yapıştır "
            "sayfası bunu yapar."))
        st.write(", ".join(rapor.backbones))

# ---------------------------------------------------------------------------
# Kodon uyumu - dizi aynı kalır, uyum değişir
# ---------------------------------------------------------------------------

k = rapor.codon_report
if k:
    st.subheader(t("Kodon uyumu"))
    if k.get("tablo_yok"):
        st.warning(k["not"], icon=":material/warning:")
    elif not k.get("satirlar"):
        st.caption(t("Kasette kodon uyumu ölçülecek bir kodlayan dizi "
                     "(CDS) yok."))
    else:
        st.dataframe(
            [{t("Gen"): s["name"] or s["part_id"],
              t("CAI (hedef)"): round(s["cai"], 3),
              t("CAI (E. coli)"): round(s["cai_kaynak"], 3),
              t("Fark"): round(s["cai"] - s["cai_kaynak"], 3)}
             for s in k["satirlar"]],
            hide_index=True, width="stretch")
        st.caption(k["not"])
        dusen = [s for s in k["satirlar"]
                 if s["cai"] < s["cai_kaynak"] - 0.05]
        if dusen:
            st.warning(
                t("Şu genlerin hedef konaktaki kodon uyumu belirgin "
                  "biçimde düşük: ")
                + ", ".join(s["name"] or s["part_id"] for s in dusen)
                + t(". Dizi analizi sayfasındaki kodon optimizasyonunu "
                    "hedef konağı seçerek çalıştırın."),
                icon=":material/warning:")

# ---------------------------------------------------------------------------
# Taşınan kaseti denetle - taşıma işe yaradı mı, denetim ne diyor
# ---------------------------------------------------------------------------

st.subheader(t("Taşınan kasetin denetimi"))
st.caption(t(
    "Taşıma bittiğinde soru şu: yeni kaset denetimden geçiyor mu? "
    "Orijin taşınamadıysa CORE-01 hatası BEKLENİR ve bu, taşımanın "
    "değil kataloğun sınırıdır."))

vec = assemble("pTasima", sort_canonical(db.get_parts(rapor.ported_parts)))
bulgular = validate(vec, host=db.get_host(hedef),
                    purpose=st.session_state.get("design_purpose", "ifade"))
sayim = summarize(bulgular)
d = st.container(horizontal=True)
d.metric(t("Uzunluk"), f"{vec.length} bp")
d.metric(t("Hata"), sayim.get("hata", 0))
d.metric(t("Uyarı"), sayim.get("uyari", 0))

for b in bulgular:
    fn = {"hata": st.error, "uyari": st.warning}.get(b.severity, st.info)
    fn(f"**{b.rule_id} · {b.title}** — {b.detail}"
       + (f"\n\n*{b.fix}*" if b.fix else ""))

# ---------------------------------------------------------------------------
# Notlar ve indirme
# ---------------------------------------------------------------------------

with st.expander(t("Taşımanın sınırları")):
    for n in rapor.notes:
        st.markdown(f"* {n}")

st.download_button(
    t("Taşıma raporunu indir (Markdown)"),
    data=porting.to_markdown(rapor).encode("utf-8"),
    file_name=f"tasima-{hedef.lower()}.md",
    mime="text/markdown", icon=":material/download:")
