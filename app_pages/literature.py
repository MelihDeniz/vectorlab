"""Literatür: arama, kaynakça, metin düzeltme, katalog besleme."""

import streamlit as st

from core import db, literature as lit
from core.assembly import sort_canonical

tabs = st.tabs(
    ["Arama", "Tasarımı doğrula", "Kaynakça", "Metin düzelt", "Katalog besle"]
)


@st.cache_data(ttl=3600, max_entries=200, show_spinner=False)
def cached_search(query: str, limit: int, bacterial: bool, oa: bool, sort: str):
    return lit.search(query, limit=limit, bacterial_only=bacterial,
                      open_access_only=oa, sort=sort)


def render_results(results, key_prefix: str):
    for r in results:
        with st.container(border=True):
            url = lit.url_for(r)
            title = f"[{r['title']}]({url})" if url else r["title"]
            st.markdown(f"**{title}**")
            bits = [r["first_author"] + " ve ark." if r["first_author"] else "",
                    f"*{r['journal']}*" if r["journal"] else "",
                    r["year"]]
            line = " · ".join(b for b in bits if b)
            if r["cited_by"]:
                line += f" · {r['cited_by']} atıf"
            if r["is_oa"]:
                line += " · açık erişim"
            st.caption(line)
            if r["abstract"]:
                with st.expander("Özet"):
                    st.write(r["abstract"])
            cols = st.container(horizontal=True)
            if cols.button("Kaynakçaya ekle", key=f"{key_prefix}_save_{r['pmid']}",
                           icon=":material/bookmark_add:"):
                db.save_reference(r)
                st.toast(f"Kaynakçaya eklendi: {r['title'][:50]}…")


# ---------------------------------------------------------------------------
# 1. Arama
# ---------------------------------------------------------------------------

with tabs[0]:
    with st.form("lit_search"):
        q = st.text_input(
            "Sorgu",
            placeholder="örn. T7 promoter leaky expression pLysS",
        )
        c = st.container(horizontal=True, vertical_alignment="bottom")
        limit = c.number_input("Sonuç", 5, 50, 15, step=5)
        sort = c.segmented_control(
            "Sıralama", ["relevance", "cited", "recent"],
            default="relevance",
            format_func=lambda s: {"relevance": "Alaka", "cited": "Çok atıf",
                                   "recent": "Yeni"}[s],
        )
        bacterial = c.toggle("Bakteriyel konuyla sınırla", value=True)
        oa = c.toggle("Sadece açık erişim", value=False)
        submitted = st.form_submit_button("Ara", icon=":material/search:")

    if submitted and q.strip():
        with st.spinner("Europe PMC taranıyor…"):
            res = cached_search(q, int(limit), bacterial,
                                oa, sort or "relevance")
        if "error" in res:
            st.error(res["error"], icon=":material/cloud_off:")
        else:
            st.caption(
                f"{res['hit_count']} kayıt bulundu, ilk {len(res['results'])} "
                "gösteriliyor."
            )
            st.session_state.lit_results = res["results"]

    if st.session_state.lit_results:
        render_results(st.session_state.lit_results, "search")

# ---------------------------------------------------------------------------
# 2. Tasarımı doğrula
# ---------------------------------------------------------------------------

with tabs[1]:
    st.caption(
        "Tasarımcı sayfasındaki yapının her parçası için literatür taraması "
        "yapar. Denetim kurallarının çıkardığı sorgular da buraya düşer."
    )
    selected = sort_canonical(db.get_parts(st.session_state.design_parts))
    if not selected:
        st.info(
            "Önce Tasarımcı sayfasında bir yapı kurun.",
            icon=":material/construction:",
        )
    else:
        st.markdown(
            f"**{st.session_state.design_name}** — {len(selected)} parça"
        )
        queries = lit.queries_for_parts(selected)
        picked = st.multiselect(
            "Hangi parçalar için arama yapılsın?",
            [name for name, _ in queries],
            default=[name for name, _ in queries][:4],
        )
        if st.button("Literatürü tara", icon=":material/travel_explore:"):
            for name, query in queries:
                if name not in picked:
                    continue
                with st.container(border=True):
                    st.markdown(f"**{name}**")
                    st.caption(f"Sorgu: `{query}`")
                    res = cached_search(query, 3, True, False, "relevance")
                    if "error" in res:
                        st.error(res["error"])
                    elif not res["results"]:
                        st.caption("Sonuç bulunamadı.")
                    else:
                        for r in res["results"]:
                            url = lit.url_for(r)
                            st.markdown(
                                f"- [{r['title']}]({url}) — "
                                f"*{r['journal']}* {r['year']}"
                                + (f" · {r['cited_by']} atıf"
                                   if r["cited_by"] else "")
                            )
                            db.save_reference(r)
                            db.link_part_ref(
                                next(p["part_id"] for p in selected
                                     if p["name"] == name),
                                r["pmid"], note=query,
                            )
            st.success(
                "Bulunan kaynaklar kaydedildi ve ilgili parçalara bağlandı. "
                "Parça kataloğu sayfasında her parçanın altında görünürler.",
                icon=":material/check_circle:",
            )

# ---------------------------------------------------------------------------
# 3. Kaynakça
# ---------------------------------------------------------------------------

with tabs[2]:
    refs = db.list_references()
    if not refs:
        st.info(
            "Kaynakça boş. Arama sekmesinden kaynak ekleyin ya da "
            "'Tasarımı doğrula' sekmesinde taramayı çalıştırın.",
            icon=":material/bookmark:",
        )
    else:
        st.caption(f"{len(refs)} kaynak kayıtlı.")
        style = st.segmented_control(
            "Biçim", ["numarali", "bibtex", "markdown"],
            default="numarali",
            format_func=lambda s: {"numarali": "Numaralı liste",
                                   "bibtex": "BibTeX",
                                   "markdown": "Markdown"}[s],
        ) or "numarali"

        # Veritabanı satırlarını biçimleyicinin beklediği şekle getir
        norm = [
            {**r, "first_author": (r.get("authors") or "").split(",")[0].strip(),
             "cited_by": 0, "volume": "", "pages": ""}
            for r in refs
        ]
        text = {
            "numarali": lit.to_numbered,
            "bibtex": lit.to_bibtex,
            "markdown": lit.to_markdown_list,
        }[style](norm)

        st.code(text, language="bibtex" if style == "bibtex" else "text",
                wrap_lines=True)
        st.download_button(
            "İndir", text,
            file_name=f"kaynakca.{'bib' if style == 'bibtex' else 'txt'}",
            icon=":material/download:",
        )

        with st.expander("Kaynakları yönet"):
            for r in refs:
                row = st.container(horizontal=True,
                                   vertical_alignment="center")
                row.markdown(f"{r['title'][:70]}…")
                if row.button("Sil", key=f"delref_{r['pmid']}",
                              icon=":material/delete:"):
                    db.delete_reference(r["pmid"])
                    st.rerun()

# ---------------------------------------------------------------------------
# 4. Metin düzelt
# ---------------------------------------------------------------------------

with tabs[3]:
    st.caption(
        "Rapor ya da protokol metninizi moleküler biyoloji yazım kurallarına "
        "göre denetler: gen ve tür adlarının italiği, enzim ve suş yazımı, "
        "birim biçimleri."
    )
    text = st.text_area(
        "Metin", height=240,
        placeholder="Raporunuzun ilgili bölümünü buraya yapıştırın…",
    )
    if text.strip():
        issues = lit.check_nomenclature(text)
        if not issues:
            st.success("Bilinen bir yazım sorunu bulunamadı.",
                       icon=":material/check_circle:")
        else:
            auto = [i for i in issues if i.auto_fixable]
            manual = [i for i in issues if not i.auto_fixable]
            st.markdown(
                f"**{len(issues)} bulgu** — {len(auto)} tanesi otomatik "
                f"düzeltilebilir."
            )
            st.dataframe(
                [{"Tür": i.kind, "Bulunan": i.found, "Önerilen": i.suggestion,
                  "Otomatik": "evet" if i.auto_fixable else "hayır",
                  "Gerekçe": i.explanation}
                 for i in issues],
                width="stretch", hide_index=True,
                alt="Metinde bulunan yazım ve isimlendirme sorunları",
            )
            if manual:
                st.info(
                    "Türkçe ek alan terimler otomatik düzeltilmez: ünlü "
                    "uyumu ve kaynaştırma harfi yüzünden kökü değiştirip "
                    "eki olduğu gibi bırakmak bozuk sonuç verir. Bu "
                    "satırları siz uyarlayın.",
                    icon=":material/edit_note:",
                )
            fixed = lit.apply_fixes(text, issues)
            st.markdown("**Düzeltilmiş metin** (italikler Markdown `*` ile)")
            st.code(fixed, language="markdown", wrap_lines=True)
            st.download_button("Düzeltilmiş metni indir", fixed,
                               file_name="duzeltilmis_metin.md",
                               icon=":material/download:")

# ---------------------------------------------------------------------------
# 5. Katalog besle
# ---------------------------------------------------------------------------

with tabs[4]:
    st.caption(
        "Makale özetinden ya da yöntem bölümünden plazmit, suş, promotör, "
        "erişim numarası ve antibiyotik derişimi çıkarır. Çıkanlar "
        "ÖNERİDİR — siz onaylamadan kataloğa yazılmaz."
    )
    src = st.radio(
        "Kaynak", ["Metin yapıştır", "Son arama sonuçlarından"],
        horizontal=True,
    )

    corpus = []
    if src == "Metin yapıştır":
        t = st.text_area("Metin", height=200)
        if t.strip():
            corpus = [("elle girilen metin", "", t)]
    else:
        if not st.session_state.lit_results:
            st.info("Önce Arama sekmesinde bir tarama yapın.",
                    icon=":material/search:")
        else:
            for r in st.session_state.lit_results:
                if r["abstract"]:
                    corpus.append((r["title"], r["pmid"], r["abstract"]))

    if corpus:
        all_c = []
        for title, pmid, body in corpus:
            for c in lit.extract_candidates(body, pmid):
                all_c.append((title, c))

        if not all_c:
            st.caption("Tanınabilir bir öğe bulunamadı.")
        else:
            st.markdown(f"**{len(all_c)} aday bulundu**")
            kinds = sorted({c.kind for _, c in all_c})
            pick = st.pills("Türe göre süz", kinds, selection_mode="multi",
                            default=kinds)
            rows = [
                {"Tür": c.kind, "Değer": c.value,
                 "Ek bilgi": ", ".join(f"{k}: {v}" for k, v in c.extras.items()
                                       if v) or "—",
                 "Bağlam": c.context,
                 "PMID": c.pmid, "Makale": title[:50]}
                for title, c in all_c if c.kind in (pick or kinds)
            ]
            st.dataframe(
                rows, width="stretch", hide_index=True,
                alt="Literatürden çıkarılan parça adayları",
                column_config={
                    "Bağlam": st.column_config.TextColumn(width="large"),
                },
            )
            st.info(
                "Bu adaylar katalog kaydı DEĞİLDİR: dizi içermezler. "
                "Kullanmak istediğiniz bir öğeyi 'Parça kataloğu' "
                "sayfasındaki 'Yeni parça ekle' formuna, erişim "
                "numarasından çektiğiniz gerçek diziyle birlikte girin.",
                icon=":material/info:",
            )
