"""Enzim kütüphanesi: nereden keser, ne uç bırakır, neyle uyumlu."""

import streamlit as st

from core import db, enzymes as ez
from core.assembly import assemble, sort_canonical
from core.seq_utils import clean_seq

lib = ez.library_size()
st.caption(
    f"REBASE'den gelen {lib['toplam']} enzim, {lib['ticari']} tanesi ticari "
    "olarak satılıyor. Kaynak: Roberts ve ark., rebase.neb.com — Biopython "
    "ile birlikte dağıtılır."
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
    ["Kütüphane", "Dizide kesim", "Uyumlu uçlar"]
)

# ---------------------------------------------------------------------------
# 1. Kütüphane
# ---------------------------------------------------------------------------

with tab_lib:
    f = st.container(horizontal=True, vertical_alignment="bottom")
    q = f.text_input("Ara", placeholder="EcoRI ya da GAATTC…", width=240)
    boy = f.selectbox(
        "Tanıma dizisi", [None, 4, 5, 6, 7, 8],
        format_func=lambda v: "hepsi" if v is None else f"{v} bp", width=130,
    )
    cikinti = f.selectbox(
        "Uç tipi", [None, "5", "3", "kor"],
        format_func=lambda v: {
            None: "hepsi", "5": "5' çıkıntı", "3": "3' çıkıntı", "kor": "kör uç"
        }[v], width=150,
    )
    ticari = f.toggle("Sadece satılanlar", value=True)
    iis = f.toggle("Sadece Tip IIS", value=False)

    sonuc = _search(q, ticari, boy, cikinti, iis)
    st.caption(f"{len(sonuc)} enzim")

    if not sonuc:
        st.info("Eşleşme yok. Filtreleri gevşetin.", icon=":material/search_off:")
    else:
        secili = st.selectbox(
            "Enzim", sonuc,
            format_func=lambda n: f"{n}  —  {_details(n)['site']}",
        )
        d = _details(secili)

        m = st.container(horizontal=True)
        m.metric("Tanıma dizisi", d["site"])
        m.metric("Uç", d["overhang"])
        m.metric("Optimum sıcaklık",
                 f"{d['opt_temp']} °C" if d["opt_temp"] else "—")
        if d["frequency"]:
            m.metric("Beklenen sıklık", f"~1/{int(d['frequency'])} bp")

        st.markdown("**Nereden keser**")
        st.code(d["diagram"], language="text")
        st.caption(
            "▼ üst zincirdeki, ▲ alt zincirdeki kesim. İkisi aynı sütundaysa "
            "kör uç, kaymışsa yapışkan uç oluşur."
        )

        st.markdown("**Kesimden sonra oluşan uçlar**")
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
                f"**Metilasyon:** {d['methylation']}. E. coli'den "
                "saflaştırılan plazmitte bu bölge kesilmeyebilir; "
                f"gerekirse {', '.join(ez.DAM_DCM_FREE_STRAINS[:2])} gibi "
                "dam⁻/dcm⁻ bir suşta çoğaltın."
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
                "Aynı bölgeyi tanırlar ama kesim noktaları ya da metilasyon "
                "duyarlılıkları farklı olabilir."
            )

# ---------------------------------------------------------------------------
# 2. Dizide kesim
# ---------------------------------------------------------------------------

with tab_cut:
    kaynak = st.radio(
        "Dizi kaynağı", ["Tasarımcıdaki yapı", "Dizi yapıştır"],
        horizontal=True,
    )

    seq = ""
    if kaynak == "Tasarımcıdaki yapı":
        ids = st.session_state.get("design_parts", [])
        if not ids:
            st.info("Tasarımcı sayfasında henüz bir yapı yok.",
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
            st.text_area("Dizi", height=120, placeholder="ATGC…")
        )

    if seq:
        secim = st.multiselect(
            "Enzimler", ez.search("", commercial_only=True),
            default=[e for e in ["EcoRI", "BamHI", "HindIII", "XhoI"]
                     if e in ez.search("", commercial_only=True)][:4],
        )
        dairesel = st.toggle("Dairesel", value=True, key="cut_circ")

        if secim:
            kesimler = ez.sites_in(seq, secim, circular=dairesel)
            toplam = sum(len(v) for v in kesimler.values())
            st.caption(f"{toplam} kesim bulundu.")

            for enz, konumlar in sorted(kesimler.items()):
                with st.container(border=True):
                    d = _details(enz)
                    if not konumlar:
                        st.markdown(
                            f"**{enz}** (`{d['site']}`) — bu dizide kesmiyor"
                        )
                        st.caption(
                            "Insert içinde kullanmak için güvenli: diziyi "
                            "parçalamaz."
                        )
                        continue
                    tanima = ez.recognition_sites(seq, [enz],
                                                  circular=dairesel)[enz]
                    st.markdown(
                        f"**{enz}** (`{d['site']}`, {d['overhang']}) — "
                        f"{len(konumlar)} kesim"
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
                            "Hangi kesim gösterilsin?", konumlar,
                            key=f"sel_{enz}",
                        )
                    ctx = ez.cut_context(seq, enz, hangi, flank=12)
                    if ctx:
                        st.code(ctx["duplex"], language="text")
                        st.caption(
                            f"Üst zincir {ctx['top_cut_pos']}. bazdan önce, "
                            f"alt zincir {ctx['bottom_cut_pos']}. bazdan "
                            f"önce kesilir · çıkıntı `{ctx['overhang_seq'] or '—'}` "
                            f"({ctx['overhang_label']})"
                        )
                        st.markdown("Kesimden sonra:")
                        st.code(ctx["products"], language="text")

# ---------------------------------------------------------------------------
# 3. Uyumlu uçlar
# ---------------------------------------------------------------------------

with tab_comp:
    st.caption(
        "Aynı çıkıntıyı bırakan enzimlerin ürünleri birbirine ligate "
        "edilebilir. Yönlü klonlamada bunlardan KAÇINMAK gerekir — iki uç "
        "uyumluysa insert ters de girer ve vektör kendi üzerine kapanır. "
        "Buna karşılık elinizde olmayan bir enzimin yerine geçecek bir uç "
        "ararken tam olarak bu liste lazımdır."
    )
    hepsi = ez.search("", commercial_only=True)
    hedef = st.selectbox("Enzim", hepsi, index=hepsi.index("XbaI")
                         if "XbaI" in hepsi else 0)
    uyumlular = _compatible(hedef)
    d = _details(hedef)
    st.markdown(
        f"**{hedef}** `{d['site']}` → çıkıntı "
        f"`{d['overhang_seq'] or 'kör'}` ({d['overhang']})"
    )
    if not uyumlular:
        st.info("Ticari enzimler arasında uyumlu uç bırakan başka enzim yok.",
                icon=":material/block:")
    else:
        st.markdown(
            f"**{len(uyumlular)} uyumlu enzim:** "
            + ", ".join(f"`{x}`" for x in uyumlular[:60])
            + (" …" if len(uyumlular) > 60 else "")
        )
        if d["blunt"]:
            st.warning(
                "Bu bir kör uç enzimi: TÜM kör uç enzimlerinin ürünleriyle "
                "ligate olur. Kör uç ligasyonunda yön kontrolü yoktur.",
                icon=":material/warning:",
            )
