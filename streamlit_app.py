"""VectorLab - bakteriyel vektör tasarım ve literatür aracı.

Çalıştırmak için:  streamlit run streamlit_app.py
"""

import streamlit as st

from core import db, i18n, notice, provenance
from core.env import SHARED_NOTICE, is_shared_deployment
from core.i18n import t

st.set_page_config(
    page_title="VectorLab",
    page_icon=":material/biotech:",
    layout="wide",
)

# Veritabanı ilk açılışta kurulur ve tohumlanır.
db.init_db()

# Dil, sayfalar çalışmadan ÖNCE ayarlanmalı: sayfa gövdesindeki her
# `t(...)` çağrısı etkin dile bakar. `core/i18n.py` Streamlit'e bağımlı
# değil, o yüzden bağlantı burada kuruluyor.
if "lang" not in st.session_state:
    st.session_state["lang"] = i18n.DEFAULT_LANG
i18n.set_language(st.session_state["lang"])

# --- Sayfalar arası paylaşılan durum
_DEFAULTS = {
    "design_name": "pVL-01",
    "design_parts": [],          # part_id listesi, montaj sırasıyla
    "design_host": "BL21DE3",
    "design_insert": "",
    "design_insert_name": "hedefGen",
    "design_purpose": "ifade",
    "design_enzymes": [],
    "analyze_seq": "",
    "lit_results": [],
}
for key, val in _DEFAULTS.items():
    if key not in st.session_state:
        st.session_state[key] = val

page = st.navigation(
    [
        st.Page(
            "app_pages/wizard.py",
            title=t("Sihirbaz"),
            icon=":material/auto_awesome:",
            default=True,
        ),
        st.Page(
            "app_pages/designer.py",
            title=t("Tasarımcı"),
            icon=":material/construction:",
        ),
        st.Page(
            "app_pages/disease.py",
            title=t("Hastalık vektörü"),
            icon=":material/medication:",
        ),
        st.Page(
            "app_pages/refbuild.py",
            title=t("Referanstan kes-yapıştır"),
            icon=":material/content_paste:",
        ),
        st.Page(
            "app_pages/vectors.py",
            title=t("Referans vektörler"),
            icon=":material/album:",
        ),
        st.Page(
            "app_pages/inspect.py",
            title=t("Plazmit denetle"),
            icon=":material/fact_check:",
        ),
        st.Page(
            "app_pages/analyze.py",
            title=t("Dizi analizi"),
            icon=":material/search_insights:",
        ),
        st.Page(
            "app_pages/primers.py",
            title=t("Primer tasarımı"),
            icon=":material/linear_scale:",
        ),
        st.Page(
            "app_pages/enzyme_lab.py",
            title=t("Enzimler"),
            icon=":material/content_cut:",
        ),
        st.Page(
            "app_pages/labwork.py",
            title=t("Laboratuvar"),
            icon=":material/science:",
        ),
        st.Page(
            "app_pages/parts.py",
            title=t("Parça kataloğu"),
            icon=":material/inventory_2:",
        ),
        st.Page(
            "app_pages/projects.py",
            title=t("Projeler"),
            icon=":material/folder_open:",
        ),
        st.Page(
            "app_pages/literature.py",
            title=t("Literatür"),
            icon=":material/menu_book:",
        ),
    ],
    position="top",
)

ust = st.container(horizontal=True, vertical_alignment="bottom")
with ust:
    st.title(page.title, icon=page.icon)
    # Dil seçicisi başlığın yanında. Değişince sayfa yeniden
    # çalıştırılır, çünkü `t(...)` çağrıları sayfa gövdesinde.
    secilen = st.selectbox(
        t("Dil"), list(i18n.LANGS),
        index=list(i18n.LANGS).index(st.session_state["lang"]),
        format_func=lambda k: i18n.LANGS[k],
        key="lang_secim", label_visibility="collapsed", width=130,
    )
if secilen != st.session_state["lang"]:
    st.session_state["lang"] = secilen
    st.rerun()

# Paylaşılan sunucuda çalışıyorsak ziyaretçi bunu en başta bilmeli.
if is_shared_deployment():
    st.info(SHARED_NOTICE, icon=":material/cloud:")

page.run()

# Sürüm ve veri kimliği: bir tasarımı sonradan açtığınızda aynı sonucu
# alıp almadığınızı ancak bu ikisi aynıysa bilebilirsiniz.
#
# Kullanım amacı beyanı buraya, her sayfanın altına konuyor: tek bir
# dosyada durursa kimse okumaz. Uzun hali DISCLAIMER.md içinde.
st.caption(f"{provenance.describe()}  \n{t(notice.FOOTER)}")

# İngilizce kipte kapsamı dürüstçe söyle. Arayüzün etiketleri ve
# açıklamaları çevrildi; parça kataloğundan, hastalık kayıt
# defterinden ve denetim kurallarından gelen BİLİMSEL METİN hâlâ
# Türkçe. "İngilizce destekliyor" deyip yarısını Türkçe bırakmak,
# kullanıcıyı yanıltmak olur.
if st.session_state["lang"] != i18n.DEFAULT_LANG:
    st.caption(
        ":material/translate: The interface is translated. Scientific "
        "content that comes from the data catalog - part notes, the "
        "disease registry, audit rule explanations and protocol text - "
        "is still in Turkish."
    )
