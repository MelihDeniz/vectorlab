"""VectorLab - bakteriyel vektör tasarım ve literatür aracı.

Çalıştırmak için:  streamlit run streamlit_app.py
"""

import streamlit as st

from core import db, provenance
from core.env import SHARED_NOTICE, is_shared_deployment

st.set_page_config(
    page_title="VectorLab",
    page_icon=":material/biotech:",
    layout="wide",
)

# Veritabanı ilk açılışta kurulur ve tohumlanır.
db.init_db()

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
            title="Sihirbaz",
            icon=":material/auto_awesome:",
            default=True,
        ),
        st.Page(
            "app_pages/designer.py",
            title="Tasarımcı",
            icon=":material/construction:",
        ),
        st.Page(
            "app_pages/vectors.py",
            title="Referans vektörler",
            icon=":material/album:",
        ),
        st.Page(
            "app_pages/inspect.py",
            title="Plazmit denetle",
            icon=":material/fact_check:",
        ),
        st.Page(
            "app_pages/analyze.py",
            title="Dizi analizi",
            icon=":material/search_insights:",
        ),
        st.Page(
            "app_pages/enzyme_lab.py",
            title="Enzimler",
            icon=":material/content_cut:",
        ),
        st.Page(
            "app_pages/labwork.py",
            title="Laboratuvar",
            icon=":material/science:",
        ),
        st.Page(
            "app_pages/parts.py",
            title="Parça kataloğu",
            icon=":material/inventory_2:",
        ),
        st.Page(
            "app_pages/literature.py",
            title="Literatür",
            icon=":material/menu_book:",
        ),
    ],
    position="top",
)

st.title(page.title, icon=page.icon)

# Paylaşılan sunucuda çalışıyorsak ziyaretçi bunu en başta bilmeli.
if is_shared_deployment():
    st.info(SHARED_NOTICE, icon=":material/cloud:")

page.run()

# Sürüm ve veri kimliği: bir tasarımı sonradan açtığınızda aynı sonucu
# alıp almadığınızı ancak bu ikisi aynıysa bilebilirsiniz.
st.caption(provenance.describe())
