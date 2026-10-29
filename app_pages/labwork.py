"""Laboratuvar: protokol, jel simülasyonu, protein özellikleri, rapor."""

import base64

import streamlit as st

from core import (
    assembly_methods as am, db, enzymes as ez, gel, protein,
    protocol as pr, ref_vectors as rv, report,
)
from core.assembly import assemble, sort_canonical
from core.cloning import fragments, plan
from core.plasmid_map import circular_map, render as render_svg
from core.seq_utils import clean_seq, find_orfs, translate_cds
from core.validate import validate

SEV_BOX = {"hata": st.error, "uyari": st.warning, "bilgi": st.info}
SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

# Klonlamada ilk akla gelen enzimler; varsayılan seçim bunların arasından,
# yapıyı tam olarak bir kez kesenlerden yapılır.
YAYGIN = ["EcoRI", "BamHI", "HindIII", "XhoI", "NdeI", "SalI", "NcoI",
          "XbaI", "KpnI", "SacI", "NotI", "SpeI", "PstI", "BglII", "SphI"]


@st.cache_data(ttl=1800, show_spinner=False)
def _tek_kesenler(sequence: str) -> list[str]:
    """Yapıyı tam olarak bir kez kesen yaygın enzimler."""
    bulunan = ez.sites_in(sequence, YAYGIN, circular=True)
    return [e for e in YAYGIN if len(bulunan.get(e, [])) == 1]


@st.cache_data(ttl=1800, show_spinner=False)
def _ticari() -> list[str]:
    return ez.search("", commercial_only=True)

# ---------------------------------------------------------------------------
# Kaynak: tasarımcıdaki yapı ya da referans vektör
# ---------------------------------------------------------------------------

kaynak = st.segmented_control(
    "Hangi yapı üzerinde çalışıyoruz?",
    ["Tasarımcıdaki yapı", "Referans vektör"],
    default="Tasarımcıdaki yapı",
) or "Tasarımcıdaki yapı"

vec = None
seq = ""
ad = ""
uzunluk = 0
gc = 0.0
parcalar: list[dict] = []
bulgular: list = []
harita_svg = ""

if kaynak == "Tasarımcıdaki yapı":
    ids = st.session_state.get("design_parts", [])
    if not ids:
        st.info(
            "Tasarımcı sayfasında henüz yapı yok. Oradan parça seçin ya da "
            "Sihirbaz'ı çalıştırın; sonra buraya dönün.",
            icon=":material/construction:",
        )
    else:
        parcalar = sort_canonical(db.get_parts(ids))
        insert = clean_seq(st.session_state.get("design_insert", ""))
        vec = assemble(
            st.session_state.get("design_name") or "pVL-01", parcalar,
            insert_seq=insert,
            insert_name=st.session_state.get("design_insert_name") or "insert",
        )
        seq, ad, uzunluk, gc = vec.sequence, vec.name, vec.length, vec.gc
        host = db.get_host(st.session_state.get("design_host", ""))
        bulgular = validate(
            vec, host=host, insert_seq=insert,
            cloning_enzymes=st.session_state.get("design_enzymes", []),
            purpose=st.session_state.get("design_purpose", "ifade"),
        )
        harita_svg = circular_map(ad, uzunluk, vec.features, {})
else:
    vlist = rv.list_vectors()
    vid = st.selectbox(
        "Vektör", [v["vector_id"] for v in vlist],
        format_func=lambda i: next(
            f"{v['name']} — {v['length']} bp" for v in vlist
            if v["vector_id"] == i),
    )
    rec = rv.load(vid)
    seq, ad, uzunluk, gc = rec["sequence"], rec["name"], rec["length"], rec["gc"]
    feats = rv.catalog_features(vid, db.list_parts())
    harita_svg = circular_map(ad, uzunluk, feats, {})
    parcalar = [
        {"name": f.name, "category": f.category, "length": f.length,
         "accession": f"{rec['accession']}:{f.start}..{f.end}"}
        for f in feats
    ]

if not seq:
    st.stop()

insert_seq = clean_seq(st.session_state.get("design_insert", ""))
insert_ad = st.session_state.get("design_insert_name") or "insert"
enzimler = st.session_state.get("design_enzymes", [])

ust = st.container(horizontal=True)
ust.metric("Yapı", ad)
ust.metric("Uzunluk", f"{uzunluk} bp")
ust.metric("GC", f"%{gc}")
if insert_seq:
    ust.metric("Insert", f"{len(insert_seq)} bp")

st.divider()

tab_pro, tab_gel, tab_prot, tab_montaj, tab_rap = st.tabs(
    ["Protokol", "Jel simülasyonu", "Protein", "Montaj yöntemleri", "Rapor"]
)

# ---------------------------------------------------------------------------
# 1. Protokol
# ---------------------------------------------------------------------------

with tab_pro:
    c = st.container(horizontal=True, vertical_alignment="bottom")
    tekler = _tek_kesenler(seq)
    secilen = c.multiselect(
        "Kesim enzimleri", _ticari(),
        default=[e for e in enzimler if e in _ticari()][:2] or tekler[:2],
        help="Tasarımcıda seçtikleriniz buraya gelir. Seçim yoksa bu "
             "yapıyı tam olarak bir kez kesen enzimlerden ikisi "
             "varsayılan olarak gelir.",
    )
    if tekler:
        c.caption(f"Tek kesen: {', '.join(tekler[:8])}")
    oran = c.slider("Ligasyon molar oranı (insert:vektör)", 1.0, 10.0, 3.0,
                    0.5, width=260)
    vektor_ng = c.number_input("Ligasyona konan vektör (ng)", 10, 200, 50,
                               step=10, width=200)

    frag = [f["length"] for f in fragments(seq, secilen)] if secilen else []
    iskelet = max(frag) if frag else uzunluk

    p = pr.build(
        vector_name=ad, vector_bp=uzunluk,
        part_ids=[x.get("part_id", "") for x in parcalar],
        enzymes=secilen, insert_bp=len(insert_seq), insert_name=insert_ad,
        host_name=(db.get_host(st.session_state.get("design_host", "")) or
                   {}).get("name", ""),
        purpose=st.session_state.get("design_purpose", "ifade"),
        backbone_bp=iskelet, fragment_sizes=frag,
    )

    for w in p.warnings:
        st.warning(w, icon=":material/warning:")

    if insert_seq and secilen:
        lig = pr.ligation_amounts(iskelet, len(insert_seq), vektor_ng, oran)
        m = st.container(horizontal=True)
        m.metric("Vektör", f"{lig['vector_ng']} ng")
        m.metric("Insert", f"{lig['insert_ng']} ng")
        m.metric("Oran", f"{oran:g}:1")
        st.caption(f"Hesap: {lig['formula']}")

    for s in p.steps:
        with st.container(border=True):
            st.markdown(f"**{s.no}. {s.title}**")
            st.markdown(s.detail)
            if s.reagents:
                st.dataframe(
                    [{"Bileşen": r.name, "Miktar": r.amount, "Not": r.note}
                     for r in s.reagents],
                    width="stretch", hide_index=True,
                    alt=f"{s.title} adımının reaksiyon kurulumu",
                )
            if s.conditions:
                st.markdown(f"**Koşullar:** {s.conditions}")
            if s.warning:
                st.warning(s.warning, icon=":material/priority_high:")

    st.download_button(
        "Protokolü Markdown olarak indir", pr.to_markdown(p),
        file_name=f"{ad}-protokol.md", mime="text/markdown",
        icon=":material/download:",
    )

# ---------------------------------------------------------------------------
# 2. Jel
# ---------------------------------------------------------------------------

with tab_gel:
    st.caption(
        "Sindirimden ne çıkacağını fragman listesi tam anlatmaz: yakın "
        "boyutlar tek bant olur, kısa fragmanlar az boya bağladığı için "
        "sönük kalır. Burası onu gösterir."
    )
    g = st.container(horizontal=True, vertical_alignment="bottom")
    jel_enz = g.multiselect(
        "Sindirim enzimleri", _ticari(),
        default=(secilen or _tek_kesenler(seq))[:2],
        help="Tanı amaçlı sindirimde genelde iki tek kesen enzim "
             "kullanılır: çıkan iki bandın boyutu yapının doğru "
             "kurulduğunu gösterir.",
    )
    merdiven = g.selectbox("Merdiven", list(gel.LADDERS), width=170)
    yukleme = g.number_input("Kuyuya yüklenen DNA (ng)", 50, 2000, 500,
                             step=50, width=220)

    if not jel_enz:
        st.info("Enzim seçin.", icon=":material/science:")
    else:
        boyutlar = [f["length"] for f in fragments(seq, jel_enz)]
        onerilen = gel.recommend_agarose(boyutlar)
        yuzde = st.select_slider(
            "Agaroz yüzdesi", list(gel.AGAROSE_RANGE), value=onerilen,
            format_func=lambda v: f"%{v}",
            help=f"Bu fragmanlar için önerilen: %{onerilen}",
        )
        kulvar = gel.build_lane(
            " + ".join(jel_enz)[:16], boyutlar, yuzde, float(yukleme))
        kesilmemis = gel.build_lane("kesilmemiş", [uzunluk], yuzde,
                                    float(yukleme))
        render_svg(gel.render_gel(
            [gel.ladder_lane(merdiven), kesilmemis, kulvar], yuzde,
            title=f"{ad} · {' + '.join(jel_enz)}",
        ))
        st.dataframe(
            [{"Bant": b.label, "Boyut (bp)": b.size,
              "Görünürlük": b.intensity,
              "Durum": {"ust": "jelin üstünde sıkışır",
                        "alt": "boya cephesinde kaybolur"}.get(
                            b.off_scale, "ayrışır")}
             for b in kulvar.bands],
            width="stretch", hide_index=True,
            alt="Jeldeki bantların boyut ve görünürlük tablosu",
            column_config={"Görünürlük": st.column_config.ProgressColumn(
                "Görünürlük", min_value=0.0, max_value=1.0, format="percent")},
        )
        for n in gel.lane_notes(kulvar, yuzde):
            st.warning(n, icon=":material/info:")
        st.caption(
            "Kesilmemiş plazmit karşılaştırma için çizildi. Gerçek jelde "
            "kesilmemiş dairesel plazmit üç bant verir (süperkıvrımlı, "
            "açık dairesel, doğrusal) ve süperkıvrımlı olan kendi "
            "boyutundan daha hızlı göçer - buradaki tek bant bu yüzden "
            "basitleştirmedir."
        )

# ---------------------------------------------------------------------------
# 3. Protein
# ---------------------------------------------------------------------------

with tab_prot:
    st.caption(
        "İfade edeceğiniz proteinin kütlesi, pI'ı ve sönüm katsayısı; "
        "yani SDS-PAGE'de bandı nerede arayacağınız ve derişimi nasıl "
        "ölçeceğiniz."
    )
    kay = st.radio("Protein dizisi nereden gelsin?",
                   ["Insert dizisinden çevir", "Protein dizisi yapıştır"],
                   horizontal=True)

    aa = ""
    if kay == "Insert dizisinden çevir":
        if not insert_seq:
            st.info("Tasarımcı sayfasında insert girilmemiş.",
                    icon=":material/science:")
        else:
            orfs = find_orfs(insert_seq, min_aa=30, circular=False)
            if not orfs:
                st.warning(
                    "Insert içinde 30 aa'dan uzun açık okuma çerçevesi "
                    "bulunamadı. Dizi doğru yönde mi, başlama kodonu var mı?",
                    icon=":material/warning:",
                )
            else:
                sec = st.selectbox(
                    "Okuma çerçevesi", orfs,
                    format_func=lambda o: o.label(),
                )
                aa = sec.protein
                st.caption(
                    f"{sec.start}..{sec.end} ({'+' if sec.strand > 0 else '−'}), "
                    f"başlama kodonu {sec.start_codon}, "
                    f"{'durdurma kodonu var' if sec.has_stop else 'DURDURMA KODONU YOK'}"
                )
    else:
        aa = st.text_area("Protein dizisi (tek harf kodu)", height=110,
                          placeholder="MKWVTFISLLFLFSSAYS…")

    ph = st.slider("Saflaştırma tamponunun pH'ı", 5.0, 9.5, 8.0, 0.1,
                   width=280)
    rep = protein.analyze(aa, working_ph=ph) if aa else None

    if rep is None and aa:
        st.warning("Dizi çok kısa ya da standart amino asit içermiyor.",
                   icon=":material/warning:")
    elif rep:
        m = st.container(horizontal=True)
        m.metric("Uzunluk", f"{rep.length_aa} aa")
        m.metric("Kütle", f"{rep.mw_kda} kDa")
        m.metric("pI", rep.pi)
        m.metric("GRAVY", rep.gravy)
        m.metric("Kararsızlık", rep.instability)
        if rep.a280_for_1mg_ml:
            m.metric("A280 (1 mg/mL)", rep.a280_for_1mg_ml)
        st.caption(protein.sds_ladder_note(rep.mw_kda))
        for n in rep.notes:
            SEV_BOX.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}",
                icon=SEV_ICON.get(n["kind"], ":material/info:"))
        with st.expander("Amino asit bileşimi"):
            st.dataframe(
                [{"Amino asit": k, "Sayı": v,
                  "Oran": round(v / rep.length_aa, 4)}
                 for k, v in sorted(rep.aa_counts.items(),
                                    key=lambda x: -x[1])],
                width="stretch", hide_index=True,
                alt="Proteinin amino asit bileşimi",
            )

# ---------------------------------------------------------------------------
# 5. Montaj yöntemleri
# ---------------------------------------------------------------------------

with tab_montaj:
    st.caption(
        "Restriksiyon + ligasyon tek yol değil. Golden Gate birden çok "
        "parçayı tek tüpte ve iz bırakmadan birleştirir; Gibson hiç enzim "
        "bölgesi gerektirmez; mutagenez tek bazı değiştirir. Her birinin "
        "kendi tuzağı var ve burada onlar denetleniyor."
    )

    yontem = st.segmented_control(
        "Yöntem", ["Golden Gate", "Gibson", "Mutagenez"],
        default="Golden Gate", key="mont_yontem",
    ) or "Golden Gate"

    if yontem == "Mutagenez":
        kalip = st.text_area(
            "Kalıp dizi", height=110,
            value=insert_seq[:600] if insert_seq else "",
            placeholder="Değiştireceğiniz bölgeyi kapsayan diziyi yapıştırın…",
            key="mont_kalip",
        )
        k = st.container(horizontal=True, vertical_alignment="bottom")
        konum = k.number_input("Değişecek ilk nükleotid", 1,
                               max(len(clean_seq(kalip)), 1), 1,
                               key="mont_konum", width=220)
        yerine = k.text_input("Yerine gelecek baz ya da kodon",
                              value="GAT", key="mont_yerine", width=220)
        kanat = k.slider("Her iki yanda eşleşen baz", 10, 25, 15,
                         key="mont_kanat", width=240)
        plan_m = am.mutagenesis(kalip, int(konum), yerine, flank=int(kanat))
    else:
        st.markdown("**Birleştirilecek parçalar** — sıra montaj sırasıdır")
        katalog = [x for x in db.list_parts() if x["seq_status"] == "full"]
        secim = st.multiselect(
            "Katalogdan parça seç",
            [x["part_id"] for x in katalog],
            default=[p for p in st.session_state.get("design_parts", [])
                     if (db.get_part(p) or {}).get("seq_status") == "full"][:4],
            format_func=lambda i: (db.get_part(i) or {}).get("name", i)[:46],
            key="mont_parcalar",
        )
        parca_listesi = []
        for pid in secim:
            rec = db.get_part(pid)
            if rec and rec["sequence"]:
                parca_listesi.append((rec["name"].split("(")[0].strip()[:14],
                                      rec["sequence"]))
        if insert_seq and st.checkbox(
                f"Inserti de ekle ({insert_ad}, {len(insert_seq)} bp)",
                value=bool(insert_seq), key="mont_insert"):
            parca_listesi.append((insert_ad[:14], insert_seq))

        o = st.container(horizontal=True, vertical_alignment="bottom")
        dairesel = o.toggle("Dairesel ürün", value=True, key="mont_dairesel")
        if yontem == "Golden Gate":
            enz = o.selectbox("Tip IIS enzim", am.TYPE_IIS_COMMON,
                              key="mont_enz", width=200)
            plan_m = am.golden_gate(parca_listesi, enz, circular=dairesel)
        else:
            ort = o.slider("Örtüşme uzunluğu (bp)", 15, 45, 25,
                           key="mont_ort", width=260)
            plan_m = am.gibson(parca_listesi, overlap=int(ort),
                               circular=dairesel)

    if not plan_m.parts and yontem != "Mutagenez":
        st.info("Önce parça seçin.", icon=":material/category:")
    else:
        mm = st.container(horizontal=True)
        mm.metric("Yöntem", plan_m.method)
        mm.metric("Parça", len(plan_m.parts))
        mm.metric("Ürün", f"{plan_m.product_length} bp")
        mm.metric("Birleşme", len(plan_m.junctions))
        mm.metric("Durum", "geçti" if plan_m.ok else "sorunlu")

        for n in plan_m.notes:
            {"hata": st.error, "uyari": st.warning,
             "bilgi": st.info}.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}")

        if plan_m.junctions:
            st.markdown("**Birleşme noktaları**")
            st.dataframe(
                [{"#": j.index, "Sol": j.left_part, "Sağ": j.right_part,
                  ("Çıkıntı" if yontem == "Golden Gate" else "Örtüşme"):
                      j.overhang,
                  "Sorun": "; ".join(j.problems) if j.problems else "—"}
                 for j in plan_m.junctions],
                width="stretch", hide_index=True,
                alt="Montaj birleşme noktaları ve sorunları",
            )

        if plan_m.primers:
            st.markdown("**Sipariş edilecek primerler**")
            st.dataframe(
                [{"Ad": p2.name, "Dizi": p2.sequence, "Uzunluk": p2.length,
                  "Bağlanma Tm": p2.tm_anneal,
                  "Uyarı": "; ".join(p2.warnings) if p2.warnings else "—"}
                 for p2 in plan_m.primers],
                width="stretch", hide_index=True,
                alt="Montaj için tasarlanan primerler",
            )
            st.download_button(
                "Primerleri CSV olarak indir",
                "ad,dizi,uzunluk,tm\n" + "\n".join(
                    f"{p2.name},{p2.sequence},{p2.length},{p2.tm_anneal}"
                    for p2 in plan_m.primers),
                file_name=f"{ad}-primerler.csv", mime="text/csv",
                icon=":material/download:",
            )
            st.caption(
                "Primerlerin 5' ucundaki kuyruk, montajın kendisini "
                "sağlayan kısımdır; Tm sütunu yalnızca kalıba bağlanan "
                "bölümün değeridir, PCR koşulunu ona göre kurun."
            )


# ---------------------------------------------------------------------------
# 6. Rapor
# ---------------------------------------------------------------------------

with tab_rap:
    st.caption(
        "Haritayı, parça listesini, denetim sonucunu, protein özelliklerini "
        "ve protokolü tek bir HTML dosyasında toplar. Tarayıcıda açıp "
        "oradan PDF'e basabilirsiniz."
    )
    r = st.container(horizontal=True, vertical_alignment="bottom")
    yazar = r.text_input("Rapora yazılacak isim (isteğe bağlı)", width=300)
    icerik = r.multiselect(
        "Rapora girecek bölümler",
        ["Harita", "Parçalar", "Denetim", "Protein", "Jel", "Protokol"],
        default=["Harita", "Parçalar", "Denetim", "Protokol"],
    )

    html_rapor = report.build_html(
        name=ad, length=uzunluk, gc=gc,
        host_name=(db.get_host(st.session_state.get("design_host", "")) or
                   {}).get("name", ""),
        purpose=st.session_state.get("design_purpose", ""),
        parts=parcalar if "Parçalar" in icerik else [],
        map_svg=harita_svg if "Harita" in icerik else "",
        findings=bulgular if "Denetim" in icerik else [],
        protein_report=rep if ("Protein" in icerik and rep) else None,
        proto=p if "Protokol" in icerik else None,
        gel_svg=(gel.render_gel(
            [gel.ladder_lane("1 kb"),
             gel.build_lane(" + ".join(jel_enz)[:16],
                            [f["length"] for f in fragments(seq, jel_enz)],
                            gel.recommend_agarose(
                                [f["length"] for f in fragments(seq, jel_enz)]))],
            gel.recommend_agarose(
                [f["length"] for f in fragments(seq, jel_enz)]),
            title=ad)
            if ("Jel" in icerik and jel_enz) else ""),
        insert_name=insert_ad, insert_bp=len(insert_seq),
        enzymes=secilen, author=yazar,
    )

    st.download_button(
        "Raporu HTML olarak indir", html_rapor,
        file_name=f"{ad}-rapor.html", mime="text/html",
        icon=":material/download:", type="primary",
    )
    st.caption(f"Rapor {len(html_rapor) // 1024} KB · harita ve jel "
               "dosyanın içine gömülür, ayrı görsel gerekmez.")

    with st.expander("Raporu burada önizle"):
        # st.iframe kaynak olarak adres ister, ham HTML almaz; raporu
        # data: adresine gömüp öyle veriyoruz - böylece dosyayı diske
        # yazmadan ve kaldırılmış bir API'ye dayanmadan gösterilebiliyor.
        src = ("data:text/html;charset=utf-8;base64,"
               + base64.b64encode(html_rapor.encode("utf-8")).decode("ascii"))
        st.iframe(src, height=620, alt=f"{ad} raporunun önizlemesi")
