"""Referans vektörler: gerçek plazmitlerin tam dizileri üzerinde çalış."""

import streamlit as st

from core import cloning, db, enzymes as ez, ref_vectors as rv
from core.parts_seed import CATEGORIES, CATEGORY_TR
from core.plasmid_map import circular_map, linear_map, render as render_svg
from core.seq_utils import clean_seq, classify_cutters, gc_percent
from core.i18n import t

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

all_parts = db.list_parts()
vectors = rv.list_vectors()

if not vectors:
    st.error(
        t("Referans vektör verisi bulunamadı. `python "
        "tools/fetch_reference_vectors.py` komutunu bir kez çalıştırın."),
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

facet = rv.seva_facets()

f1 = st.container(horizontal=True, vertical_alignment="bottom")
with f1:
    organizmalar = ["tümü"] + rv.organisms()
    org = f1.selectbox(
        t("Hedef organizma"), organizmalar,
        format_func=lambda o: "Tümü" if o == "tümü" else o, width=260,
        help=t("Vektörün hangi bakteride çoğalabildiği. E. coli plazmiti "
             "Bacillus'ta çoğalmaz; buradan süzerek doğru replikonu "
             "seçersiniz."),
    )
    kullanim = f1.selectbox(t("Amaç"), ["tümü", "klonlama", "ifade",
                                     "genom düzenleme"],
                            format_func=lambda u: u.capitalize(), width=180)
    arama = f1.text_input(t("Ara"), placeholder=t("pUC, KanR, p15A…"), width=200)

st.caption(
    t("Katalogda {0} referans vektör var; {1} farklı kargo ve altı marker ile SEVA koleksiyonunun tamamı dahil. Aşağıdaki hanelerle süzebilirsiniz.").format(len(vectors), len(facet['cargo']))
)

f2 = st.container(horizontal=True, vertical_alignment="bottom")
with f2:
    s_marker = f2.selectbox(
        t("SEVA markerı"), ["tümü"] + list(facet["marker"]),
        format_func=lambda k: ("Tümü" if k == "tümü"
                               else f"{k} — {facet['marker'][k]}"), width=240)
    s_ori = f2.selectbox(
        t("SEVA orijini"), ["tümü"] + list(facet["ori"]),
        format_func=lambda k: ("Tümü" if k == "tümü"
                               else f"{k} — {facet['ori'][k]}"), width=240)
    s_cargo = f2.selectbox(
        t("SEVA kargosu"), ["tümü"] + list(facet["cargo"]),
        format_func=lambda k: ("Tümü" if k == "tümü"
                               else f"{k} — {facet['cargo'][k]}"), width=300)

ids = rv.filter_vectors(
    organism=None if org == "tümü" else org,
    use=None if kullanim == "tümü" else kullanim,
    seva_marker=None if s_marker == "tümü" else s_marker,
    seva_ori=None if s_ori == "tümü" else s_ori,
    seva_cargo=None if s_cargo == "tümü" else s_cargo,
    search=arama,
)
if not ids:
    st.info(t("Bu süzgeçlerle eşleşen vektör yok. Filtreleri gevşetin."),
            icon=":material/search_off:")
    st.stop()

top = st.container(horizontal=True, vertical_alignment="bottom")
with top:
    vid = st.selectbox(t("Vektör ({0})").format(len(ids)), ids,
                       format_func=lambda i: labels[i], width=400)
    benzerlik = st.slider(
        t("Benzerlik eşiği (%)"), 85.0, 100.0, 90.0, 0.5, width=240,
        help=t("Katalog parçaları bu orandan daha benzer bölgelerde "
             "etiketlenir. %100 yalnızca birebir eşleşme demektir."),
    )

rec = rv.load(vid)
feats = _features(vid, benzerlik)

m = st.container(horizontal=True)
m.metric(t("Uzunluk"), f"{rec['length']} bp")
m.metric(t("GC"), f"%{rec['gc']}")
m.metric(t("Kopya sayısı"), rec.get("kopya", "—"))
m.metric(t("Uyumsuzluk grubu"), rec.get("inc_group", "—"))
m.metric(t("Organizma"), rec.get("organizma", "E. coli"))
m.metric(t("Etiketlenen bölge"), len(feats))

with st.container(border=True):
    st.markdown(f"**{rec['name']}** — {rec.get('ozet', '')}")
    c = st.container(horizontal=True)
    c.markdown(f"**Marker:** {', '.join(rec.get('markerlar') or ['—'])}")
    c.markdown(t("**Promotör:** {0}").format(rec.get('promotor', '—')))
    st.markdown(t("**Klonlama bölgesi:** {0}").format(rec.get('mcs', '—')))
    st.markdown(rec.get("notlar", ""))
    st.caption(
        t("Dizi kaynağı: {0} · kayıt {1} · {2} bp · indirilme {3}").format(rec['source'], rec['accession'], rec['length'], rec['fetched'])
    )

act = st.container(horizontal=True)
act.download_button(
    t("FASTA indir"), rv.to_fasta(vid), file_name=f"{rec['name']}.fasta",
    mime="text/plain", icon=":material/download:",
)
if act.button(t("Dizi analizine gönder"), icon=":material/search_insights:"):
    st.session_state.analyze_seq = rec["sequence"]
    st.switch_page("app_pages/analyze.py")
if act.button(t("Etiketlenen parçaları tasarımcıya aktar"),
              icon=":material/construction:"):
    picked = {c: [] for c in CATEGORIES}
    for f in feats:
        if f.part_id and f.part_id not in picked.get(f.category, []):
            picked.setdefault(f.category, []).append(f.part_id)
    for cat, pids in picked.items():
        st.session_state[f"ms_{cat}"] = pids
    st.session_state.design_name = t("{0}-türevi").format(rec['name'])
    st.switch_page("app_pages/designer.py")

st.divider()

tab_map, tab_feat, tab_enz, tab_clone, tab_part = st.tabs(
    [t("Harita"), t("Öznitelikler"), t("Enzim haritası"), t("Klonlama planı"),
     t("Parça çıkar")]
)

# ---------------------------------------------------------------------------
# 1. Harita
# ---------------------------------------------------------------------------

with tab_map:
    single = {k: v for k, v in _cutters(vid)["tek"]}
    goster = st.multiselect(
        t("Haritada gösterilecek tek kesen enzimler"),
        sorted(single),
        default=sorted(single)[:10],
        help=t("{0} içinde tam olarak bir kez kesen {1} ticari enzim var.").format(rec['name'], len(single)),
    )
    sites = {k: single[k] for k in goster}
    render_svg(circular_map(rec["name"], rec["length"], feats, sites))
    render_svg(linear_map(rec["name"], rec["length"], feats, sites))
    st.caption(
        t("Etiketler katalog parçalarının bu dizide bulunduğu yerlerden "
        "gelir; GenBank'ın kendi öznitelikleri Öznitelikler sekmesinde.")
    )

# ---------------------------------------------------------------------------
# 2. Öznitelikler
# ---------------------------------------------------------------------------

with tab_feat:
    st.markdown(t("**Katalog eşleşmeleri**"))
    if not feats:
        st.info(t("Bu eşikte katalog parçası eşleşmedi."),
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
            alt=t("{0} üzerinde bulunan katalog parçaları tablosu").format(rec['name']),
        )
        kismi = [f for f in feats if f.identity < 100]
        if kismi:
            st.caption(
                t("%100'ün altındaki satırlar nokta farkı taşıyan "
                "varyantlardır - aynı gen, farklı plazmit soyu. Örneğin "
                "`bla` pBR322'de birebir, pUC soyunda birkaç nükleotid "
                "farklıdır. Hizalama yapılmadığı için yalnızca aynı "
                "uzunluktaki farklar yakalanır.")
            )

    gb = rv.genbank_features(vid)
    st.markdown(t("**Kaydın kendi öznitelikleri** ({0})").format(len(gb)))
    if not gb:
        st.caption(
            t("Bu GenBank kaydı öznitelik tablosu taşımıyor - pUC19 ve pUC18 "
            "kayıtları yalnızca diziden oluşur. Haritadaki etiketler bu "
            "yüzden katalog taramasından geliyor.")
        )
    else:
        st.dataframe(
            [{"Tür": g["type"], "Başlangıç": g["start"], "Bitiş": g["end"],
              "Yön": "+" if g["strand"] > 0 else "−",
              "Etiket": g["label"][:70]} for g in gb],
            width="stretch", hide_index=True,
            alt=t("{0} GenBank öznitelik tablosu").format(rec['name']),
        )

# ---------------------------------------------------------------------------
# 3. Enzim haritası
# ---------------------------------------------------------------------------

with tab_enz:
    gruplar = _cutters(vid)
    g = st.container(horizontal=True)
    g.metric(t("Tek kesen"), len(gruplar["tek"]))
    g.metric(t("İki kesen"), len(gruplar["cift"]))
    g.metric(t("Çok kesen"), len(gruplar["cok"]))
    g.metric(t("Kesmeyen"), len(gruplar["kesmez"]))
    st.caption(
        t("Tek kesen enzimler klonlama için kullanılabilir; kesmeyenler "
        "inserti bozmadığı için insert içine eklenecek bölgeler olarak işe "
        "yarar.")
    )

    tek = {k: v for k, v in gruplar["tek"]}
    if tek:
        tanimalar = ez.recognition_sites(rec["sequence"], sorted(tek))
        secili = st.selectbox(
            t("Kesimi incele"), sorted(tek),
            format_func=lambda n: (
                t("{0} — tanıma {1}, kesim {2}").format(n, tanimalar[n][0] if tanimalar[n] else '?', tek[n][0])),
        )
        konum = tek[secili][0]
        if tanimalar[secili]:
            bas, son = ez.site_span(rec["sequence"], secili, konum)
            st.caption(
                t("Tanıma dizisi {0}..{1}, üst zincirdeki kesim {2}. bazdan önce. Yayımlanmış plazmit haritaları genellikle TANIMA başlangıcını verir, bu yüzden iki sayı farklıdır - ikisi de doğrudur.").format(bas, son, konum)
            )
        ctx = ez.cut_context(rec["sequence"], secili, konum, flank=12)
        if ctx:
            st.code(ctx["duplex"], language="text")
            st.caption(
                t("Üst zincir {0}. bazdan önce, alt zincir {1}. bazdan önce kesilir · çıkıntı `{2}` ({3})").format(ctx['top_cut_pos'], ctx['bottom_cut_pos'], ctx['overhang_seq'] or '—', ctx['overhang_label'])
            )
            st.markdown(t("Kesimden sonra:"))
            st.code(ctx["products"], language="text")

    with st.expander(t("Kesmeyen enzimler ({0})").format(len(gruplar['kesmez']))):
        st.caption(
            t("Bu enzimlerin tanıma dizisi vektörde hiç yok. Inserte primerle "
            "böyle bir bölge eklerseniz vektörü parçalamadan "
            "kullanabilirsiniz.")
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
        t("Kesim enzimleri (yönlü klonlama için iki farklı enzim)"),
        secenekler,
        default=[e for e in sorted(tek) if e in ("EcoRI", "BamHI")][:2],
        max_selections=3,
    )
    kaynak = st.radio(
        t("Insert dizisi"), ["Tasarımcıdaki insert", "Buraya yapıştır"],
        horizontal=True,
    format_func=t)
    if kaynak == "Tasarımcıdaki insert":
        ins = clean_seq(st.session_state.get("design_insert", ""))
        ad = st.session_state.get("design_insert_name") or "insert"
        st.caption(f"{ad}: {len(ins)} nt" if ins
                   else "Tasarımcı sayfasında insert girilmemiş.")
    else:
        ins = clean_seq(st.text_area(t("Insert"), height=110,
                                     placeholder=t("ATGAGC…")))
        ad = "insert"

    p = cloning.plan(rec["sequence"], rec["name"], enz, ins, ad)

    if enz:
        k = st.container(horizontal=True)
        k.metric(t("İskelet"), f"{p.backbone_length} bp")
        k.metric(t("Atılan parça"), f"{p.removed_length} bp")
        k.metric(t("Beklenen yapı"), f"{p.product_length} bp")
        if ins:
            k.metric(t("Insert GC"), f"%{gc_percent(ins)}")

        for n in p.notes:
            box = {"hata": st.error, "uyari": st.warning,
                   "bilgi": st.info}[n.severity]
            box(f"**{n.title}** — {n.detail}", icon=SEV_ICON[n.severity])

        frag = cloning.fragments(rec["sequence"], enz)
        st.markdown(t("**Kesimden çıkan fragmanlar**"))
        st.dataframe(
            [{"Başlangıç": f["start"], "Bitiş": f["end"],
              "Uzunluk (bp)": f["length"], "Not": f["note"]} for f in frag],
            width="stretch", hide_index=True,
            alt=t("{0} sindirim fragmanları tablosu").format(rec['name']),
        )
        st.caption(
            t("Jelde bu boyutları görmeniz beklenir. Yakın boyutlu fragmanlar "
            "aynı bantta görünür; 50 bp'nin altındakiler çoğu jelde "
            "kaybolur.")
        )

# ---------------------------------------------------------------------------
# 5. Parça çıkar
# ---------------------------------------------------------------------------

with tab_part:
    st.caption(
        t("Vektörün bir bölgesini katalog parçası olarak kaydedin. Kayıt "
        "kaynağını (aksesyon ve koordinat) taşır, böylece dizinin nereden "
        "geldiği sonradan da belli olur.")
    )

    secim = ["— elle koordinat gir —"] + [
        f"{f.name}  [{f.start}..{f.end} {'+' if f.strand > 0 else '−'}, "
        f"%{f.identity}]" for f in feats
    ]
    hangi = st.selectbox(t("Bölge"), secim)

    if hangi == secim[0]:
        c = st.container(horizontal=True)
        bas = c.number_input(t("Başlangıç"), 1, rec["length"], 1, width=140)
        bit = c.number_input(t("Bitiş"), 1, rec["length"],
                             min(300, rec["length"]), width=140)
        yon = c.selectbox(t("Yön"), [1, -1],
                          format_func=lambda v: "+" if v > 0 else "−",
                          width=100)
        on_ad = f"{rec['name']} {bas}..{bit}"
        on_kat = "cds"
    else:
        f = feats[secim.index(hangi) - 1]
        bas, bit, yon = f.start, f.end, f.strand
        on_ad = f"{f.name} ({rec['name']})"
        on_kat = f.category
        st.caption(t("Seçilen bölge: {0}..{1}, {2} zincir, {3} bp").format(bas, bit, '+' if yon > 0 else '−', f.length))

    dizi = rv.region(vid, bas, bit, yon)
    if dizi:
        st.code(
            "\n".join(dizi[i:i + 60] for i in range(0, min(len(dizi), 300), 60))
            + ("\n…" if len(dizi) > 300 else ""),
            language="text",
        )
        st.caption(f"{len(dizi)} nt · GC %{gc_percent(dizi)}")

    form = st.container(horizontal=True, vertical_alignment="bottom")
    ad_yeni = form.text_input(t("Parça adı"), value=on_ad, width=320)
    kat = form.selectbox(
        t("Kategori"), CATEGORIES,
        index=CATEGORIES.index(on_kat) if on_kat in CATEGORIES else 0,
        format_func=lambda c: CATEGORY_TR.get(c, c), width=220,
    )
    pid = form.text_input(
        t("Parça kimliği"), value=f"REF_{vid}_{bas}_{bit}", width=220
    )
    not_yeni = st.text_area(
        t("Not"), height=80,
        value=t("{0} ({1}) üzerindeki {2}..{3} bölgesinden çıkarıldı.").format(rec['name'], rec['accession'], bas, bit),
    )

    if st.button(t("Kataloğa ekle"), icon=":material/add:", type="primary"):
        if not dizi:
            st.error(t("Bölge boş."), icon=":material/error:")
        elif not pid.strip() or not ad_yeni.strip():
            st.error(t("Kimlik ve ad boş olamaz."), icon=":material/error:")
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
                t("`{0}` kataloğa eklendi - Parça kataloğu ve Tasarımcı sayfalarında görünür.").format(pid.strip()),
                icon=":material/check_circle:",
            )
