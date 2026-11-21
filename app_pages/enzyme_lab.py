"""Enzim kütüphanesi: nereden keser, ne uç bırakır, neyle uyumlu."""

import streamlit as st

from core import db, enzymes as ez, notice
from core.assembly import assemble, sort_canonical
from core.seq_utils import clean_seq
from core.i18n import t

lib = ez.library_size()
st.caption(
    t("REBASE'den gelen {0} enzim, {1} tanesi ticari olarak satılıyor; Biopython ile birlikte dağıtılır. REBASE'in istediği künye: {2}").format(lib['toplam'], lib['ticari'], notice.REBASE_CITATION)
)


@st.cache_data(ttl=3600, show_spinner=False)
def _search(q, ticari, boy, cikinti, iis):
    return ez.search(q, commercial_only=ticari, site_length=boy,
                     overhang=cikinti, type_iis_only=iis)


@st.cache_data(ttl=3600, show_spinner=False)
def _details(name):
    return ez.details(name)


@st.cache_data(ttl=3600, show_spinner=False)
def _compatible(name):
    return ez.compatible_with(name)


tab_lib, tab_cut, tab_comp = st.tabs(
    [t("Kütüphane"), t("Dizide kesim"), t("Uyumlu uçlar")]
)

# ---------------------------------------------------------------------------
# 1. Kütüphane
# ---------------------------------------------------------------------------

with tab_lib:
    f = st.container(horizontal=True, vertical_alignment="bottom")
    q = f.text_input(t("Ara"), placeholder=t("EcoRI ya da GAATTC…"), width=240)
    boy = f.selectbox(
        t("Tanıma dizisi"), [None, 4, 5, 6, 7, 8],
        format_func=lambda v: "hepsi" if v is None else f"{v} bp", width=130,
    )
    cikinti = f.selectbox(
        t("Uç tipi"), [None, "5", "3", "kor"],
        format_func=lambda v: {
            None: "hepsi", "5": "5' çıkıntı", "3": "3' çıkıntı", "kor": "kör uç"
        }[v], width=150,
    )
    ticari = f.toggle(t("Sadece satılanlar"), value=True)
    iis = f.toggle(t("Sadece Tip IIS"), value=False)

    sonuc = _search(q, ticari, boy, cikinti, iis)
    st.caption(f"{len(sonuc)} enzim")

    if not sonuc:
        st.info(t("Eşleşme yok. Filtreleri gevşetin."), icon=":material/search_off:")
    else:
        secili = st.selectbox(
            t("Enzim"), sonuc,
            format_func=lambda n: f"{n}  —  {_details(n)['site']}",
        )
        d = _details(secili)

        m = st.container(horizontal=True)
        m.metric(t("Tanıma dizisi"), d["site"])
        m.metric(t("Uç"), d["overhang"])
        m.metric(t("Optimum sıcaklık"),
                 f"{d['opt_temp']} °C" if d["opt_temp"] else "—")
        if d["frequency"]:
            m.metric(t("Beklenen sıklık"), f"~1/{int(d['frequency'])} bp")

        st.markdown(t("**Nereden keser**"))
        st.code(d["diagram"], language="text")
        st.caption(
            t("▼ üst zincirdeki, ▲ alt zincirdeki kesim. İkisi aynı sütundaysa "
            "kör uç, kaymışsa yapışkan uç oluşur.")
        )

        st.markdown(t("**Kesimden sonra oluşan uçlar**"))
        st.code(ez.site_products(secili), language="text")

        bilgi = []
        if d["type_iis"]:
            bilgi.append(
                "**Tip IIS:** tanıma dizisinin DIŞINDAN keser. Çıkıntı "
                "dizisi tanıma bölgesinde değil çevredeki DNA'da olduğu "
                "için diyagramda N görünür — Golden Gate montajının "
                "çalışma mantığı budur: çıkıntıyı siz tasarlarsınız, "
                "böylece birden fazla parça tek tüpte doğru sırayla "
                "birleşir ve bölge üründe kalmaz."
            )
        if d["blunt"]:
            bilgi.append(
                "**Kör uç:** yönlü klonlama yapılamaz, insert iki yönde de "
                "girer. Ligasyon verimi yapışkan uca göre düşüktür."
            )
        if d["methylation"]:
            bilgi.append(
                t("**Metilasyon:** {0}. E. coli'den saflaştırılan plazmitte bu bölge kesilmeyebilir; gerekirse {1} gibi dam⁻/dcm⁻ bir suşta çoğaltın.").format(d['methylation'], ', '.join(ez.DAM_DCM_FREE_STRAINS[:2]))
            )
        if not d["commercial"]:
            bilgi.append(
                "**Ticari olarak satılmıyor:** REBASE'de kayıtlı ama "
                "katalogdan sipariş edemezsiniz. İzoşizomerlerinden birini "
                "arayın."
            )
        for b in bilgi:
            st.warning(b, icon=":material/info:")

        if d["isoschizomers"]:
            st.markdown(
                "**İzoşizomerler** (aynı diziyi tanır): "
                + ", ".join(f"`{x}`" for x in d["isoschizomers"])
            )
            st.caption(
                t("Aynı bölgeyi tanırlar ama kesim noktaları ya da metilasyon "
                "duyarlılıkları farklı olabilir.")
            )

# ---------------------------------------------------------------------------
# 2. Dizide kesim
# ---------------------------------------------------------------------------

with tab_cut:
    kaynak = st.radio(
        t("Dizi kaynağı"), ["Tasarımcıdaki yapı", "Dizi yapıştır"],
        horizontal=True,
    format_func=t)

    seq = ""
    if kaynak == "Tasarımcıdaki yapı":
        ids = st.session_state.get("design_parts", [])
        if not ids:
            st.info(t("Tasarımcı sayfasında henüz bir yapı yok."),
                    icon=":material/construction:")
        else:
            vec = assemble(
                st.session_state.get("design_name", "yapı"),
                sort_canonical(db.get_parts(ids)),
                insert_seq=clean_seq(st.session_state.get("design_insert", "")),
            )
            seq = vec.sequence
            st.caption(f"{vec.name} — {len(seq)} bp")
    else:
        seq = clean_seq(
            st.text_area(t("Dizi"), height=120, placeholder=t("ATGC…"))
        )

    if seq:
        secim = st.multiselect(
            t("Enzimler"), ez.search("", commercial_only=True),
            default=[e for e in ["EcoRI", "BamHI", "HindIII", "XhoI"]
                     if e in ez.search("", commercial_only=True)][:4],
        )
        dairesel = st.toggle(t("Dairesel"), value=True, key="cut_circ")

        if secim:
            kesimler = ez.sites_in(seq, secim, circular=dairesel)
            toplam = sum(len(v) for v in kesimler.values())
            st.caption(t("{0} kesim bulundu.").format(toplam))

            for enz, konumlar in sorted(kesimler.items()):
                with st.container(border=True):
                    d = _details(enz)
                    if not konumlar:
                        st.markdown(
                            t("**{0}** (`{1}`) — bu dizide kesmiyor").format(enz, d['site'])
                        )
                        st.caption(
                            t("Insert içinde kullanmak için güvenli: diziyi "
                            "parçalamaz.")
                        )
                        continue
                    tanima = ez.recognition_sites(seq, [enz],
                                                  circular=dairesel)[enz]
                    st.markdown(
                        t("**{0}** (`{1}`, {2}) — {3} kesim").format(enz, d['site'], d['overhang'], len(konumlar))
                    )
                    st.caption(
                        "Tanıma dizisi başlangıcı: "
                        + ", ".join(str(k) for k in tanima[:12])
                        + (" … · " if len(tanima) > 12 else " · ")
                        + "üst zincirde kesim: "
                        + ", ".join(str(k) for k in konumlar[:12])
                        + (" …" if len(konumlar) > 12 else "")
                    )
                    hangi = konumlar[0]
                    if len(konumlar) > 1:
                        hangi = st.selectbox(
                            t("Hangi kesim gösterilsin?"), konumlar,
                            key=f"sel_{enz}",
                        )
                    ctx = ez.cut_context(seq, enz, hangi, flank=12)
                    if ctx:
                        st.code(ctx["duplex"], language="text")
                        st.caption(
                            t("Üst zincir {0}. bazdan önce, alt zincir {1}. bazdan önce kesilir · çıkıntı `{2}` ({3})").format(ctx['top_cut_pos'], ctx['bottom_cut_pos'], ctx['overhang_seq'] or '—', ctx['overhang_label'])
                        )
                        st.markdown(t("Kesimden sonra:"))
                        st.code(ctx["products"], language="text")

# ---------------------------------------------------------------------------
# 3. Uyumlu uçlar
# ---------------------------------------------------------------------------

with tab_comp:
    st.caption(
        t("Aynı çıkıntıyı bırakan enzimlerin ürünleri birbirine ligate "
        "edilebilir. Yönlü klonlamada bunlardan KAÇINMAK gerekir — iki uç "
        "uyumluysa insert ters de girer ve vektör kendi üzerine kapanır. "
        "Buna karşılık elinizde olmayan bir enzimin yerine geçecek bir uç "
        "ararken tam olarak bu liste lazımdır.")
    )
    hepsi = ez.search("", commercial_only=True)
    hedef = st.selectbox(t("Enzim"), hepsi, index=hepsi.index("XbaI")
                         if "XbaI" in hepsi else 0)
    uyumlular = _compatible(hedef)
    d = _details(hedef)
    st.markdown(
        t("**{0}** `{1}` → çıkıntı `{2}` ({3})").format(hedef, d['site'], d['overhang_seq'] or 'kör', d['overhang'])
    )
    if not uyumlular:
        st.info(t("Ticari enzimler arasında uyumlu uç bırakan başka enzim yok."),
                icon=":material/block:")
    else:
        st.markdown(
            t("**{0} uyumlu enzim:** ").format(len(uyumlular))
            + ", ".join(f"`{x}`" for x in uyumlular[:60])
            + (" …" if len(uyumlular) > 60 else "")
        )
        if d["blunt"]:
            st.warning(
                t("Bu bir kör uç enzimi: TÜM kör uç enzimlerinin ürünleriyle "
                "ligate olur. Kör uç ligasyonunda yön kontrolü yoktur."),
                icon=":material/warning:",
            )
