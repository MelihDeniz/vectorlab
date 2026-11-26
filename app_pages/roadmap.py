"""Gerçekleştirme yol haritası: tasarımdan tedavi adayına.

Tasarım bitince "artık tedavi var" olmuyor. Arada basamakları belli
bir yol var ve bu yolun en büyük sorunu basamakların **sırası**:
hayvan deneyine etik kurul onayından önce başlanmaz, içerme testi
yapılmadan canlı GDO hayvana verilmez, ifade gösterilmeden işlev
aranmaz. Sırayı bilmemek bir projeyi yanlış promotör seçmekten daha
kesin batırır.

Sayfa statik bir metin basmaz: yol haritası **tasarımın kendisinden**
türetilir (`core/roadmap.py`). Kasette sinyal peptidi varsa
salgılanma doğrulaması eklenir, konak Gram-pozitifse dönüşüm
optimizasyonu gelir, denetimden çıkan her bulgu proje riskine
dönüşür, hastalık kaydından gelen tasarımda in vivo okuma o
hastalığın biyobelirteci olur.
"""

import streamlit as st

from core import db, disease_seed as ds, procure, roadmap as rm
from core.assembly import assemble, sort_canonical
from core.i18n import t
from core.validate import validate

st.caption(t(
    "Tasarım bitince iş bitmiyor: sıradaki soru bu vektörü gerçekten "
    "nasıl yapacağınız ve hangi basamaklardan geçeceğiniz. Bu sayfa "
    "yol haritasını tasarımınızdan üretir - sinyal peptidi varsa "
    "salgılanma doğrulaması, Gram-pozitif konakta dönüşüm "
    "optimizasyonu, denetimden çıkan her bulgu da proje riski olarak "
    "eklenir."
))

parts_by_id = {p["part_id"]: p for p in db.list_parts()}

# ---------------------------------------------------------------------------
# Tasarım kaynağı
# ---------------------------------------------------------------------------

kaynak = st.segmented_control(
    t("Tasarım nereden gelsin?"),
    ["Tasarımcı'daki yapı", "Hastalık kayıt defterinden"],
    default="Tasarımcı'daki yapı",
    format_func=t) or "Tasarımcı'daki yapı"

hastalik = None
if kaynak == "Hastalık kayıt defterinden":
    h_key = st.selectbox(
        t("Hastalık"), ds.keys(),
        format_func=lambda k: ds.by_key(k).name)
    hastalik = ds.by_key(h_key)
    s_key = st.selectbox(
        t("Strateji"), [s.key for s in hastalik.strategies],
        format_func=lambda k: next(
            f"{s.key} · {s.title}" for s in hastalik.strategies
            if s.key == k))
    strateji = next(s for s in hastalik.strategies if s.key == s_key)
    kaset, konak_id = list(strateji.parts), strateji.host
    ad = f"p{hastalik.short.split()[0]}-{strateji.key}"
    amac = hastalik.purpose
else:
    kaset = [p for p in st.session_state.get("design_parts", [])
             if p in parts_by_id]
    konak_id = st.session_state.get("design_host") or ""
    ad = st.session_state.get("design_name") or "pVL-01"
    amac = st.session_state.get("design_purpose", "ifade")
    if not kaset:
        st.info(t("Tasarımcı sayfasında henüz bir yapı yok. Oradan "
                  "parça seçin ya da yukarıdan hastalık kayıt "
                  "defterinden bir tasarım alın."),
                icon=":material/science:")
        st.stop()

# ---------------------------------------------------------------------------
# Yol haritasını kur
# ---------------------------------------------------------------------------

vec = assemble(ad, sort_canonical(db.get_parts(kaset)))
bulgular = validate(vec, host=db.get_host(konak_id) if konak_id else None,
                    purpose=amac)

in_vivo = st.toggle(
    t("Hayvan modeli projenin kapsamında"), value=True,
    help=t("Kapalıysa in vivo aşaması ve ona bağlı izinler yol "
           "haritasından çıkarılır. Lisans projelerinin çoğu burada "
           "durur ve bunu başvuruda açıkça söylemek, "
           "söylememekten iyidir."))

yol = rm.build(ad, kaset, host_id=konak_id, disease=hastalik,
               findings=bulgular, in_vivo=in_vivo)

o = st.container(horizontal=True)
o.metric(t("Aşama"), len(yol.phases))
o.metric(t("Toplam süre"), f"{yol.ay:g} " + t("ay"))
o.metric(t("İzin"), len(yol.approvals))
o.metric(t("Risk"), len(yol.risks))

if yol.readout:
    st.success(t("Başarıyı ölçecek okuma: ") + yol.readout,
               icon=":material/target:")

if not in_vivo:
    st.info(t("Hayvan modeli kapsam dışı bırakıldı: in vivo aşaması "
              "ve ona bağlı izinler yol haritasında yok."),
            icon=":material/info:")

# ---------------------------------------------------------------------------
# Aşamalar
# ---------------------------------------------------------------------------

st.subheader(t("Basamak basamak"))
st.caption(t(
    "Her aşamanın bir **kapısı** var: o kapı açılmadan sonrakine "
    "geçilmez. Kapılar süslü değil, elemeli - ifade gösterilmeden "
    "işlev aranmaz, içerme ölçülmeden canlı GDO hayvana verilmez."
))

for p in yol.phases:
    with st.expander(
            t("Aşama {0} — {1}  ·  {2} iş günü").format(
                p.no, p.title, f"{p.gun:g}"),
            expanded=p.no <= 2):
        st.markdown(f"*{p.amac}*")
        for i, gorev in enumerate(p.tasks, start=1):
            st.markdown(t("**{0}.{1} · {2}**  `{3:g} gün`").format(p.no, i, gorev.title, gorev.gun))
            if gorev.kosul:
                st.caption(t("Bu adım şu yüzden var: ") + gorev.kosul)
            if gorev.approval:
                izin = rm.APPROVALS[gorev.approval]
                st.warning(t("Önce izin: ") + izin.title,
                           icon=":material/gavel:")
            st.markdown(t("**Nasıl:** ") + gorev.nasil)
            c = st.container(horizontal=True)
            c.markdown(t("**Gereken:** ") + gorev.gereken)
            c.markdown(t("**Çıktı:** ") + gorev.cikti)
            if gorev.neden:
                st.caption(t("Neden önemli: ") + gorev.neden)
            st.divider()
        if p.kapi:
            st.success(t("Sonraki aşamaya geçiş koşulu: ") + p.kapi,
                       icon=":material/door_front:")

# ---------------------------------------------------------------------------
# Takvim
# ---------------------------------------------------------------------------

st.subheader(t("Takvim"))
st.dataframe(
    [{t("Aşama"): f"{p.no}. {p.title}",
      t("İş günü"): p.gun,
      t("Ay"): round(p.gun / 21.0, 1),
      t("Görev"): len(p.tasks)} for p in yol.phases]
    + [{t("Aşama"): t("Toplam"), t("İş günü"): yol.gun,
        t("Ay"): yol.ay,
        t("Görev"): sum(len(p.tasks) for p in yol.phases)}],
    hide_index=True, width="stretch")
st.caption(t(
    "Süreler tahmindir ve kaynaktan gelmez: tek kişinin çalıştığı "
    "varsayılan, yerleşik akademik iş akışının alışıldık süreleridir. "
    "Toplam adımların toplamıdır; paralelleştirme varsayılmadı."
))

# ---------------------------------------------------------------------------
# İzinler
# ---------------------------------------------------------------------------

if yol.approvals:
    st.subheader(t("İzinler ve yasal dayanaklar"))
    for a in yol.approvals:
        with st.expander(f"{a.title} — {a.body}"):
            st.markdown(t("**Dayanak:** ") + a.basis)
            st.markdown(t("**Ne zaman:** ") + a.when)
            if a.note:
                st.markdown(a.note)
            if a.caution:
                st.warning(a.caution, icon=":material/warning:")

# ---------------------------------------------------------------------------
# Riskler
# ---------------------------------------------------------------------------

st.subheader(t("Risk yönetimi"))
if yol.risks:
    st.caption(t(
        "Denetimde görünen bir eksiğin başvuruda risk olarak "
        "yazılmaması, hakemin ilk yakalayacağı şeydir. Tersi de "
        "doğru: bilinen bir açığı yazıp önlemini söylemek başvuruyu "
        "zayıf değil GÜÇLÜ gösterir."
    ))
    for k in yol.risks:
        with st.expander(f"{k.title}"):
            st.markdown(t("**Etki:** ") + k.etki)
            st.markdown(t("**Önlem:** ") + k.onlem)
            if k.kaynak:
                st.caption(t("Kaynak: ") + k.kaynak)
else:
    st.info(t("Denetimden bulgu çıkmadı ve tasarımda bilinen bir "
              "risk deseni yok. Yine de başvuruya kendi "
              "risklerinizi yazın: risksiz proje olmaz, risksiz "
              "görünen başvuru inandırıcı olmaz."),
            icon=":material/info:")

# ---------------------------------------------------------------------------
# Notlar ve indirme
# ---------------------------------------------------------------------------

with st.expander(t("Yol haritasının sınırları")):
    for n in yol.notes:
        st.markdown(f"* {n}")

st.warning(t(
    "Bu sayfa bir düzenleyici danışmanlık değildir. Yasal dayanaklar "
    "künyesiyle yazılıdır ama mevzuat değişir - nitekim Biyogüvenlik "
    "Kurulu 2018'de kaldırılmıştır. Her izin için güncel durumu "
    "ilgili kurumdan kendiniz teyit edin."
), icon=":material/gavel:")

st.subheader(t("Proje başvurusu taslağı"))
st.caption(t(
    "Yol haritasından bir başvuru iskeleti üretilir: iş paketleri "
    "aşamalardan, takvim sürelerden, risk yönetimi denetim "
    "bulgularından, kaynakça hastalık kaydından gelir. Köşeli "
    "parantezli bölümler sizin yazmanız gerekenlerdir - araç ile "
    "sizin yazdığınız yerler karışmamalı."
))

bom = procure.build(ad, part_ids=kaset,
                    host_name=(db.get_host(konak_id) or {}).get("name", ""))

d = st.container(horizontal=True)
d.download_button(
    t("Yol haritası (Markdown)"),
    data=rm.to_markdown(yol).encode("utf-8"),
    file_name=f"{ad}-yol-haritasi.md", mime="text/markdown",
    icon=":material/download:")
d.download_button(
    t("Proje başvurusu taslağı (Markdown)"),
    data=rm.to_proposal(yol, disease=hastalik, bom=bom).encode("utf-8"),
    file_name=f"{ad}-basvuru-taslagi.md", mime="text/markdown",
    icon=":material/download:")

with st.expander(t("Başvuru taslağını burada oku")):
    st.markdown(rm.to_proposal(yol, disease=hastalik, bom=bom))
