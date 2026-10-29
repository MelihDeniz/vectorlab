"""Referans vektörler: gerçek plazmitlerin tam dizileri üzerinde çalış."""

import streamlit as st

from core import cloning, db, enzymes as ez, ref_vectors as rv
from core.parts_seed import CATEGORIES, CATEGORY_TR
from core.plasmid_map import circular_map, linear_map, render as render_svg
from core.seq_utils import clean_seq, classify_cutters, gc_percent

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

all_parts = db.list_parts()
vectors = rv.list_vectors()

if not vectors:
    st.error(
        "Referans vektör verisi bulunamadı. `python "
        "tools/fetch_reference_vectors.py` komutunu bir kez çalıştırın.",
        icon=":material/cloud_off:",
    )
    st.stop()


@st.cache_data(ttl=3600, show_spinner=False)
def _features(vector_id: str, min_identity: float):
    return rv.catalog_features(vector_id, db.list_parts(), min_identity)


@st.cache_data(ttl=3600, show_spinner=False)
def _cutters(vector_id: str):
    rec = rv.load(vector_id)
    return classify_cutters(rec["sequence"], ez.search("", commercial_only=True))


@st.cache_data(ttl=3600, show_spinner=False)
def _commercial():
    return ez.search("", commercial_only=True)


# ---------------------------------------------------------------------------
# Vektör seçimi
# ---------------------------------------------------------------------------

labels = {v["vector_id"]: f"{v['name']} — {v['length']} bp · {v['kullanim']}"
          for v in vectors}

top = st.container(horizontal=True, vertical_alignment="bottom")
with top:
    organizmalar = ["tümü"] + rv.organisms()
    org = st.selectbox(
        "Hedef organizma", organizmalar,
        format_func=lambda o: "Tümü" if o == "tümü" else o, width=300,
        help="Vektörün hangi bakteride çoğalabildiği. E. coli plazmiti "
             "Bacillus'ta çoğalmaz; buradan süzerek doğru replikonu "
             "seçersiniz.",
    )
    ids = ([v["vector_id"] for v in vectors] if org == "tümü"
           else rv.by_organism(org))
    vid = st.selectbox("Vektör", ids, format_func=lambda i: labels[i],
                       width=330)
    benzerlik = st.slider(
        "Benzerlik eşiği (%)", 85.0, 100.0, 90.0, 0.5, width=240,
        help="Katalog parçaları bu orandan daha benzer bölgelerde "
             "etiketlenir. %100 yalnızca birebir eşleşme demektir.",
    )

rec = rv.load(vid)
feats = _features(vid, benzerlik)

m = st.container(horizontal=True)
m.metric("Uzunluk", f"{rec['length']} bp")
m.metric("GC", f"%{rec['gc']}")
m.metric("Kopya sayısı", rec.get("kopya", "—"))
m.metric("Uyumsuzluk grubu", rec.get("inc_group", "—"))
m.metric("Organizma", rec.get("organizma", "E. coli"))
m.metric("Etiketlenen bölge", len(feats))

with st.container(border=True):
    st.markdown(f"**{rec['name']}** — {rec.get('ozet', '')}")
    c = st.container(horizontal=True)
    c.markdown(f"**Marker:** {', '.join(rec.get('markerlar') or ['—'])}")
    c.markdown(f"**Promotör:** {rec.get('promotor', '—')}")
    st.markdown(f"**Klonlama bölgesi:** {rec.get('mcs', '—')}")
    st.markdown(rec.get("notlar", ""))
    st.caption(
        f"Dizi kaynağı: {rec['source']} · kayıt {rec['accession']} · "
        f"{rec['length']} bp · indirilme {rec['fetched']}"
    )

act = st.container(horizontal=True)
act.download_button(
    "FASTA indir", rv.to_fasta(vid), file_name=f"{rec['name']}.fasta",
    mime="text/plain", icon=":material/download:",
)
if act.button("Dizi analizine gönder", icon=":material/search_insights:"):
    st.session_state.analyze_seq = rec["sequence"]
    st.switch_page("app_pages/analyze.py")
if act.button("Etiketlenen parçaları tasarımcıya aktar",
              icon=":material/construction:"):
    picked = {c: [] for c in CATEGORIES}
    for f in feats:
        if f.part_id and f.part_id not in picked.get(f.category, []):
            picked.setdefault(f.category, []).append(f.part_id)
    for cat, pids in picked.items():
        st.session_state[f"ms_{cat}"] = pids
    st.session_state.design_name = f"{rec['name']}-türevi"
    st.switch_page("app_pages/designer.py")

st.divider()

tab_map, tab_feat, tab_enz, tab_clone, tab_part = st.tabs(
    ["Harita", "Öznitelikler", "Enzim haritası", "Klonlama planı",
     "Parça çıkar"]
)

# ---------------------------------------------------------------------------
# 1. Harita
# ---------------------------------------------------------------------------

with tab_map:
    single = {k: v for k, v in _cutters(vid)["tek"]}
    goster = st.multiselect(
        "Haritada gösterilecek tek kesen enzimler",
        sorted(single),
        default=sorted(single)[:10],
        help=f"{rec['name']} içinde tam olarak bir kez kesen "
             f"{len(single)} ticari enzim var.",
    )
    sites = {k: single[k] for k in goster}
    render_svg(circular_map(rec["name"], rec["length"], feats, sites))
    render_svg(linear_map(rec["name"], rec["length"], feats, sites))
    st.caption(
        "Etiketler katalog parçalarının bu dizide bulunduğu yerlerden "
        "gelir; GenBank'ın kendi öznitelikleri Öznitelikler sekmesinde."
    )

# ---------------------------------------------------------------------------
# 2. Öznitelikler
# ---------------------------------------------------------------------------

with tab_feat:
    st.markdown("**Katalog eşleşmeleri**")
    if not feats:
        st.info("Bu eşikte katalog parçası eşleşmedi.",
                icon=":material/search_off:")
    else:
        st.dataframe(
            [{
                "Başlangıç": f.start, "Bitiş": f.end,
                "Yön": "+" if f.strand > 0 else "−",
                "Uzunluk": f.length,
                "Benzerlik": f"%{f.identity}",
                "Tür": CATEGORY_TR.get(f.category, f.category),
                "Parça": f.name,
            } for f in feats],
            width="stretch", hide_index=True,
            alt=f"{rec['name']} üzerinde bulunan katalog parçaları tablosu",
        )
        kismi = [f for f in feats if f.identity < 100]
        if kismi:
            st.caption(
                "%100'ün altındaki satırlar nokta farkı taşıyan "
                "varyantlardır - aynı gen, farklı plazmit soyu. Örneğin "
                "`bla` pBR322'de birebir, pUC soyunda birkaç nükleotid "
                "farklıdır. Hizalama yapılmadığı için yalnızca aynı "
                "uzunluktaki farklar yakalanır."
            )

    gb = rv.genbank_features(vid)
    st.markdown(f"**Kaydın kendi öznitelikleri** ({len(gb)})")
    if not gb:
        st.caption(
            "Bu GenBank kaydı öznitelik tablosu taşımıyor - pUC19 ve pUC18 "
            "kayıtları yalnızca diziden oluşur. Haritadaki etiketler bu "
            "yüzden katalog taramasından geliyor."
        )
    else:
        st.dataframe(
            [{"Tür": g["type"], "Başlangıç": g["start"], "Bitiş": g["end"],
              "Yön": "+" if g["strand"] > 0 else "−",
              "Etiket": g["label"][:70]} for g in gb],
            width="stretch", hide_index=True,
            alt=f"{rec['name']} GenBank öznitelik tablosu",
        )

# ---------------------------------------------------------------------------
# 3. Enzim haritası
# ---------------------------------------------------------------------------

with tab_enz:
    gruplar = _cutters(vid)
    g = st.container(horizontal=True)
    g.metric("Tek kesen", len(gruplar["tek"]))
    g.metric("İki kesen", len(gruplar["cift"]))
    g.metric("Çok kesen", len(gruplar["cok"]))
    g.metric("Kesmeyen", len(gruplar["kesmez"]))
    st.caption(
        "Tek kesen enzimler klonlama için kullanılabilir; kesmeyenler "
        "inserti bozmadığı için insert içine eklenecek bölgeler olarak işe "
        "yarar."
    )

    tek = {k: v for k, v in gruplar["tek"]}
    if tek:
        tanimalar = ez.recognition_sites(rec["sequence"], sorted(tek))
        secili = st.selectbox(
            "Kesimi incele", sorted(tek),
            format_func=lambda n: (
                f"{n} — tanıma {tanimalar[n][0] if tanimalar[n] else '?'}, "
                f"kesim {tek[n][0]}"),
        )
        konum = tek[secili][0]
        if tanimalar[secili]:
            bas, son = ez.site_span(rec["sequence"], secili, konum)
            st.caption(
                f"Tanıma dizisi {bas}..{son}, üst zincirdeki kesim "
                f"{konum}. bazdan önce. Yayımlanmış plazmit haritaları "
                f"genellikle TANIMA başlangıcını verir, bu yüzden iki "
                f"sayı farklıdır - ikisi de doğrudur."
            )
        ctx = ez.cut_context(rec["sequence"], secili, konum, flank=12)
        if ctx:
            st.code(ctx["duplex"], language="text")
            st.caption(
                f"Üst zincir {ctx['top_cut_pos']}. bazdan önce, alt zincir "
                f"{ctx['bottom_cut_pos']}. bazdan önce kesilir · çıkıntı "
                f"`{ctx['overhang_seq'] or '—'}` ({ctx['overhang_label']})"
            )
            st.markdown("Kesimden sonra:")
            st.code(ctx["products"], language="text")

    with st.expander(f"Kesmeyen enzimler ({len(gruplar['kesmez'])})"):
        st.caption(
            "Bu enzimlerin tanıma dizisi vektörde hiç yok. Inserte primerle "
            "böyle bir bölge eklerseniz vektörü parçalamadan "
            "kullanabilirsiniz."
        )
        st.markdown(", ".join(f"`{n}`" for n, _ in gruplar["kesmez"][:150])
                    + (" …" if len(gruplar["kesmez"]) > 150 else ""))

# ---------------------------------------------------------------------------
# 4. Klonlama planı
# ---------------------------------------------------------------------------

with tab_clone:
    tek = {k: v for k, v in _cutters(vid)["tek"]}
    sec = st.container(horizontal=True)
    secenekler = _commercial()
    enz = sec.multiselect(
        "Kesim enzimleri (yönlü klonlama için iki farklı enzim)",
        secenekler,
        default=[e for e in sorted(tek) if e in ("EcoRI", "BamHI")][:2],
        max_selections=3,
    )
    kaynak = st.radio(
        "Insert dizisi", ["Tasarımcıdaki insert", "Buraya yapıştır"],
        horizontal=True,
    )
    if kaynak == "Tasarımcıdaki insert":
        ins = clean_seq(st.session_state.get("design_insert", ""))
        ad = st.session_state.get("design_insert_name") or "insert"
        st.caption(f"{ad}: {len(ins)} nt" if ins
                   else "Tasarımcı sayfasında insert girilmemiş.")
    else:
        ins = clean_seq(st.text_area("Insert", height=110,
                                     placeholder="ATGAGC…"))
        ad = "insert"

    p = cloning.plan(rec["sequence"], rec["name"], enz, ins, ad)

    if enz:
        k = st.container(horizontal=True)
        k.metric("İskelet", f"{p.backbone_length} bp")
        k.metric("Atılan parça", f"{p.removed_length} bp")
        k.metric("Beklenen yapı", f"{p.product_length} bp")
        if ins:
            k.metric("Insert GC", f"%{gc_percent(ins)}")

        for n in p.notes:
            box = {"hata": st.error, "uyari": st.warning,
                   "bilgi": st.info}[n.severity]
            box(f"**{n.title}** — {n.detail}", icon=SEV_ICON[n.severity])

        frag = cloning.fragments(rec["sequence"], enz)
        st.markdown("**Kesimden çıkan fragmanlar**")
        st.dataframe(
            [{"Başlangıç": f["start"], "Bitiş": f["end"],
              "Uzunluk (bp)": f["length"], "Not": f["note"]} for f in frag],
            width="stretch", hide_index=True,
            alt=f"{rec['name']} sindirim fragmanları tablosu",
        )
        st.caption(
            "Jelde bu boyutları görmeniz beklenir. Yakın boyutlu fragmanlar "
            "aynı bantta görünür; 50 bp'nin altındakiler çoğu jelde "
            "kaybolur."
        )

# ---------------------------------------------------------------------------
# 5. Parça çıkar
# ---------------------------------------------------------------------------

with tab_part:
    st.caption(
        "Vektörün bir bölgesini katalog parçası olarak kaydedin. Kayıt "
        "kaynağını (aksesyon ve koordinat) taşır, böylece dizinin nereden "
        "geldiği sonradan da belli olur."
    )

    secim = ["— elle koordinat gir —"] + [
        f"{f.name}  [{f.start}..{f.end} {'+' if f.strand > 0 else '−'}, "
        f"%{f.identity}]" for f in feats
    ]
    hangi = st.selectbox("Bölge", secim)

    if hangi == secim[0]:
        c = st.container(horizontal=True)
        bas = c.number_input("Başlangıç", 1, rec["length"], 1, width=140)
        bit = c.number_input("Bitiş", 1, rec["length"],
                             min(300, rec["length"]), width=140)
        yon = c.selectbox("Yön", [1, -1],
                          format_func=lambda v: "+" if v > 0 else "−",
                          width=100)
        on_ad = f"{rec['name']} {bas}..{bit}"
        on_kat = "cds"
    else:
        f = feats[secim.index(hangi) - 1]
        bas, bit, yon = f.start, f.end, f.strand
        on_ad = f"{f.name} ({rec['name']})"
        on_kat = f.category
        st.caption(f"Seçilen bölge: {bas}..{bit}, "
                   f"{'+' if yon > 0 else '−'} zincir, {f.length} bp")

    dizi = rv.region(vid, bas, bit, yon)
    if dizi:
        st.code(
            "\n".join(dizi[i:i + 60] for i in range(0, min(len(dizi), 300), 60))
            + ("\n…" if len(dizi) > 300 else ""),
            language="text",
        )
        st.caption(f"{len(dizi)} nt · GC %{gc_percent(dizi)}")

    form = st.container(horizontal=True, vertical_alignment="bottom")
    ad_yeni = form.text_input("Parça adı", value=on_ad, width=320)
    kat = form.selectbox(
        "Kategori", CATEGORIES,
        index=CATEGORIES.index(on_kat) if on_kat in CATEGORIES else 0,
        format_func=lambda c: CATEGORY_TR.get(c, c), width=220,
    )
    pid = form.text_input(
        "Parça kimliği", value=f"REF_{vid}_{bas}_{bit}", width=220
    )
    not_yeni = st.text_area(
        "Not", height=80,
        value=f"{rec['name']} ({rec['accession']}) üzerindeki "
              f"{bas}..{bit} bölgesinden çıkarıldı.",
    )

    if st.button("Kataloğa ekle", icon=":material/add:", type="primary"):
        if not dizi:
            st.error("Bölge boş.", icon=":material/error:")
        elif not pid.strip() or not ad_yeni.strip():
            st.error("Kimlik ve ad boş olamaz.", icon=":material/error:")
        else:
            db.upsert_part({
                "part_id": pid.strip(),
                "name": ad_yeni.strip(),
                "category": kat,
                "sequence": dizi,
                "seq_status": "full",
                "notes": not_yeni.strip(),
                "accession": f"{rec['accession']}:{bas}..{bit}"
                             f"{'' if yon > 0 else ' (ters zincir)'}",
                "source": "user",
            })
            st.success(
                f"`{pid.strip()}` kataloğa eklendi - Parça kataloğu ve "
                f"Tasarımcı sayfalarında görünür.",
                icon=":material/check_circle:",
            )
