"""Sihirbaz: birkaç soruya göre gerekçeli vektör önerisi."""

import streamlit as st

from core import db
from core.assembly import assemble, sort_canonical, to_genbank
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import circular_map, render as render_svg
from core.recommend import QUESTIONS, default_answers, recommend, visible_questions
from core.seq_utils import restriction_sites
from core.validate import summarize, validate

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

st.caption(
    "Birkaç soru soruyor, çalışan bir vektör öneriyor ve her parçayı neden "
    "seçtiğini yazıyor. Öneri bir başlangıç noktasıdır - Tasarımcı sayfasına "
    "aktarıp istediğini değiştirebilirsin."
)

# ---------------------------------------------------------------------------
# Sorular
# ---------------------------------------------------------------------------

answers = default_answers()
for q in QUESTIONS:
    key = f"wz_{q['key']}"
    if key in st.session_state:
        answers[q["key"]] = st.session_state[key]

left, right = st.columns([0.37, 0.63], gap="large")

with left:
    st.subheader("Sorular")
    shown = visible_questions(answers)
    for q in shown:
        opts = list(q["secenekler"].keys())
        key = f"wz_{q['key']}"
        if key not in st.session_state:
            st.session_state[key] = q["varsayilan"]
        st.radio(
            q["soru"],
            opts,
            format_func=lambda o, s=q["secenekler"]: s[o],
            key=key,
        )
        answers[q["key"]] = st.session_state[key]

    hidden = [q for q in QUESTIONS if q not in shown]
    if hidden:
        st.caption(
            f"{len(hidden)} soru bu seçimlerde anlamsız olduğu için "
            "gösterilmiyor."
        )

# ---------------------------------------------------------------------------
# Öneri
# ---------------------------------------------------------------------------

rec = recommend(answers)
parts = sort_canonical(db.get_parts(rec.part_ids))
host = db.get_host(rec.host_id)
vec = assemble("pVL-oneri", parts)
findings = validate(vec, host=host, purpose=rec.purpose)
counts = summarize(findings)

with right:
    m = st.container(horizontal=True)
    m.metric("Uzunluk", f"{vec.length} bp")
    m.metric("Parça", len(parts))
    m.metric("Hata", counts["hata"])
    m.metric("Uyarı", counts["uyari"])

    tab_oneri, tab_harita, tab_denetim, tab_protokol = st.tabs(
        ["Öneri ve gerekçeler", "Harita", "Denetim", "Protokol"]
    )

    with tab_oneri:
        with st.container(border=True):
            st.markdown(f"**Konak suşu — {host['name']}**")
            st.markdown(rec.host_reason)
            st.caption(f"Genotip: `{host['genotype']}`")

        for c in rec.choices:
            p = next((x for x in parts if x["part_id"] == c.part_id), None)
            if not p:
                continue
            with st.container(border=True):
                st.markdown(f"**{c.role} — {p['name']}**")
                st.markdown(c.reason)
                meta = []
                if p["strength"]:
                    meta.append(f"güç: {p['strength']}")
                if p["induction"]:
                    meta.append(f"indüksiyon: {p['induction']}")
                if p["copy_number"]:
                    meta.append(f"kopya sayısı: {p['copy_number']}")
                meta.append(f"{p['length']} bp")
                st.caption(" · ".join(meta))

        for n in rec.notes:
            st.info(n, icon=":material/lightbulb:")

    with tab_harita:
        single = {
            k: v for k, v in restriction_sites(vec.sequence).items() if len(v) == 1
        }
        render_svg(circular_map(vec.name, vec.length, vec.features, single))
        if single:
            st.caption(
                "Tek kesen enzimler: " + ", ".join(sorted(single)[:18])
                + (" …" if len(single) > 18 else "")
            )

    with tab_denetim:
        if not findings:
            st.success("Denetimden geçti - bilinen bir sorun yok.",
                       icon=":material/check_circle:")
        for f in findings:
            with st.container(border=True):
                st.markdown(
                    f"{SEV_ICON.get(f.severity, '')} **{f.title}**  "
                    f"`{f.rule_id}`"
                )
                st.markdown(f.detail)
                if f.fix:
                    st.markdown(f"**Ne yapmalı:** {f.fix}")

    with tab_protokol:
        st.markdown("**Bu yapıyla nasıl çalışırsın**")
        for i, adim in enumerate(rec.protocol, 1):
            st.markdown(f"{i}. {adim}")
        st.caption(
            "Bu adımlar genel başlangıç koşullarıdır; kendi proteinin için "
            "optimizasyon gerekebilir."
        )

    st.divider()
    actions = st.container(horizontal=True)
    if actions.button("Tasarımcıya aktar", type="primary",
                      icon=":material/construction:"):
        CATS = ["promoter", "operator", "rbs", "tag", "protease_site", "linker",
                "mcs", "reporter", "terminator", "ori", "marker", "regulator"]
        by_id = {p["part_id"]: p for p in parts}
        for cat in CATS:
            st.session_state[f"ms_{cat}"] = [
                pid for pid in rec.part_ids
                if by_id.get(pid, {}).get("category") == cat
            ]
        st.session_state.design_host = rec.host_id
        st.session_state.design_purpose = rec.purpose
        st.session_state.design_name = "pVL-oneri"
        st.switch_page("app_pages/designer.py")

    actions.download_button(
        "GenBank indir", to_genbank(vec),
        file_name=f"{vec.name}.gb", mime="chemical/seq-na-genbank",
        icon=":material/download:",
    )
