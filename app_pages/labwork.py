"""Laboratuvar: protokol, jel simülasyonu, protein özellikleri, rapor."""

import base64

import streamlit as st

from core import (
    assembly_methods as am, db, enzymes as ez, gel, protein,
    protocol as pr, ref_vectors as rv, report, troubleshoot as ts,
    verification as vf,
)
from core.assembly import assemble, sort_canonical
from core.cloning import fragments, plan
from core.plasmid_map import circular_map, render as render_svg
from core.seq_utils import clean_seq, find_orfs, translate_cds
from core.validate import validate
from core.i18n import t

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
    t("Hangi yapı üzerinde çalışıyoruz?"),
    ["Tasarımcıdaki yapı", "Referans vektör"],
    default="Tasarımcıdaki yapı",
    format_func=t) or "Tasarımcıdaki yapı"

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
            t("Tasarımcı sayfasında henüz yapı yok. Oradan parça seçin ya da "
            "Sihirbaz'ı çalıştırın; sonra buraya dönün."),
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
        t("Vektör"), [v["vector_id"] for v in vlist],
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
ust.metric(t("Yapı"), ad)
ust.metric(t("Uzunluk"), f"{uzunluk} bp")
ust.metric(t("GC"), f"%{gc}")
if insert_seq:
    ust.metric(t("Insert"), f"{len(insert_seq)} bp")

st.divider()

(tab_pro, tab_gel, tab_prot, tab_montaj, tab_dogrula, tab_sorun,
 tab_rap) = st.tabs(
    [t("Protokol"), t("Jel simülasyonu"), t("Protein"), t("Montaj yöntemleri"),
     t("Dizileme doğrulaması"), t("Sorun giderme"), t("Rapor")]
)

# ---------------------------------------------------------------------------
# 1. Protokol
# ---------------------------------------------------------------------------

with tab_pro:
    c = st.container(horizontal=True, vertical_alignment="bottom")
    tekler = _tek_kesenler(seq)
    secilen = c.multiselect(
        t("Kesim enzimleri"), _ticari(),
        default=[e for e in enzimler if e in _ticari()][:2] or tekler[:2],
        help=t("Tasarımcıda seçtikleriniz buraya gelir. Seçim yoksa bu "
             "yapıyı tam olarak bir kez kesen enzimlerden ikisi "
             "varsayılan olarak gelir."),
    )
    if tekler:
        c.caption(f"Tek kesen: {', '.join(tekler[:8])}")
    oran = c.slider(t("Ligasyon molar oranı (insert:vektör)"), 1.0, 10.0, 3.0,
                    0.5, width=260)
    vektor_ng = c.number_input(t("Ligasyona konan vektör (ng)"), 10, 200, 50,
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
        m.metric(t("Vektör"), f"{lig['vector_ng']} ng")
        m.metric(t("Insert"), f"{lig['insert_ng']} ng")
        m.metric(t("Oran"), f"{oran:g}:1")
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
                    alt=t("{0} adımının reaksiyon kurulumu").format(s.title),
                )
            if s.conditions:
                st.markdown(t("**Koşullar:** {0}").format(s.conditions))
            if s.warning:
                st.warning(s.warning, icon=":material/priority_high:")

    st.download_button(
        t("Protokolü Markdown olarak indir"), pr.to_markdown(p),
        file_name=f"{ad}-protokol.md", mime="text/markdown",
        icon=":material/download:",
    )

# ---------------------------------------------------------------------------
# 2. Jel
# ---------------------------------------------------------------------------

with tab_gel:
    st.caption(
        t("Sindirimden ne çıkacağını fragman listesi tam anlatmaz: yakın "
        "boyutlar tek bant olur, kısa fragmanlar az boya bağladığı için "
        "sönük kalır. Burası onu gösterir.")
    )
    g = st.container(horizontal=True, vertical_alignment="bottom")
    jel_enz = g.multiselect(
        t("Sindirim enzimleri"), _ticari(),
        default=(secilen or _tek_kesenler(seq))[:2],
        help=t("Tanı amaçlı sindirimde genelde iki tek kesen enzim "
             "kullanılır: çıkan iki bandın boyutu yapının doğru "
             "kurulduğunu gösterir."),
    )
    merdiven = g.selectbox(t("Merdiven"), list(gel.LADDERS), width=170)
    yukleme = g.number_input(t("Kuyuya yüklenen DNA (ng)"), 50, 2000, 500,
                             step=50, width=220)

    if not jel_enz:
        st.info(t("Enzim seçin."), icon=":material/science:")
    else:
        boyutlar = [f["length"] for f in fragments(seq, jel_enz)]
        onerilen = gel.recommend_agarose(boyutlar)
        yuzde = st.select_slider(
            t("Agaroz yüzdesi"), list(gel.AGAROSE_RANGE), value=onerilen,
            format_func=lambda v: f"%{v}",
            help=t("Bu fragmanlar için önerilen: %{0}").format(onerilen),
        )
        kulvar = gel.build_lane(
            " + ".join(jel_enz)[:16], boyutlar, yuzde, float(yukleme))
        kesilmemis = gel.uncut_lane("kesilmemiş", uzunluk, yuzde,
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
                t("Görünürlük"), min_value=0.0, max_value=1.0, format="percent")},
        )
        for n in gel.lane_notes(kulvar, yuzde):
            st.warning(n, icon=":material/info:")
        st.caption(
            t("Kesilmemiş plazmit karşılaştırma kulvarı üç bant olarak "
            "çizildi: süperkıvrımlı, doğrusal ve açık dairesel. "
            "Süperkıvrımlı biçim sıkı paketlendiği için kendi boyutundan "
            "hızlı göçer, bu yüzden o banttan merdivene bakarak boyut "
            "okuyamazsınız. Kesiminizin tam olup olmadığını anlamanın "
            "yolu, bu üçlünün kesilmiş kulvarda KAYBOLMUŞ olmasıdır.")
        )

# ---------------------------------------------------------------------------
# 3. Protein
# ---------------------------------------------------------------------------

with tab_prot:
    st.caption(
        t("İfade edeceğiniz proteinin kütlesi, pI'ı ve sönüm katsayısı; "
        "yani SDS-PAGE'de bandı nerede arayacağınız ve derişimi nasıl "
        "ölçeceğiniz.")
    )
    kay = st.radio(t("Protein dizisi nereden gelsin?"),
                   ["Insert dizisinden çevir", "Protein dizisi yapıştır"],
                   horizontal=True, format_func=t)

    aa = ""
    if kay == "Insert dizisinden çevir":
        if not insert_seq:
            st.info(t("Tasarımcı sayfasında insert girilmemiş."),
                    icon=":material/science:")
        else:
            orfs = find_orfs(insert_seq, min_aa=30, circular=False)
            if not orfs:
                st.warning(
                    t("Insert içinde 30 aa'dan uzun açık okuma çerçevesi "
                    "bulunamadı. Dizi doğru yönde mi, başlama kodonu var mı?"),
                    icon=":material/warning:",
                )
            else:
                sec = st.selectbox(
                    t("Okuma çerçevesi"), orfs,
                    format_func=lambda o: o.label(),
                )
                aa = sec.protein
                st.caption(
                    t("{0}..{1} ({2}), başlama kodonu {3}, {4}").format(sec.start, sec.end, '+' if sec.strand > 0 else '−', sec.start_codon, 'durdurma kodonu var' if sec.has_stop else 'DURDURMA KODONU YOK')
                )
    else:
        aa = st.text_area(t("Protein dizisi (tek harf kodu)"), height=110,
                          placeholder=t("MKWVTFISLLFLFSSAYS…"))

    ph = st.slider(t("Saflaştırma tamponunun pH'ı"), 5.0, 9.5, 8.0, 0.1,
                   width=280)
    rep = protein.analyze(aa, working_ph=ph) if aa else None

    if rep is None and aa:
        st.warning(t("Dizi çok kısa ya da standart amino asit içermiyor."),
                   icon=":material/warning:")
    elif rep:
        m = st.container(horizontal=True)
        m.metric(t("Uzunluk"), f"{rep.length_aa} aa")
        m.metric(t("Kütle"), f"{rep.mw_kda} kDa")
        m.metric(t("pI"), rep.pi)
        m.metric(t("GRAVY"), rep.gravy)
        m.metric(t("Kararsızlık"), rep.instability)
        if rep.a280_for_1mg_ml:
            m.metric(t("A280 (1 mg/mL)"), rep.a280_for_1mg_ml)
        st.caption(protein.sds_ladder_note(rep.mw_kda))
        for n in rep.notes:
            SEV_BOX.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}",
                icon=SEV_ICON.get(n["kind"], ":material/info:"))
        with st.expander(t("Amino asit bileşimi")):
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
        t("Restriksiyon + ligasyon tek yol değil. Golden Gate birden çok "
        "parçayı tek tüpte ve iz bırakmadan birleştirir; Gibson hiç enzim "
        "bölgesi gerektirmez; mutagenez tek bazı değiştirir. Her birinin "
        "kendi tuzağı var ve burada onlar denetleniyor.")
    )

    yontem = st.segmented_control(
        t("Yöntem"), ["Golden Gate", "Gibson", "Mutagenez"],
        default="Golden Gate", key="mont_yontem",
    ) or "Golden Gate"

    if yontem == "Mutagenez":
        kalip = st.text_area(
            t("Kalıp dizi"), height=110,
            value=insert_seq[:600] if insert_seq else "",
            placeholder=t("Değiştireceğiniz bölgeyi kapsayan diziyi yapıştırın…"),
            key="mont_kalip",
        )
        k = st.container(horizontal=True, vertical_alignment="bottom")
        konum = k.number_input(t("Değişecek ilk nükleotid"), 1,
                               max(len(clean_seq(kalip)), 1), 1,
                               key="mont_konum", width=220)
        yerine = k.text_input(t("Yerine gelecek baz ya da kodon"),
                              value="GAT", key="mont_yerine", width=220)
        kanat = k.slider(t("Her iki yanda eşleşen baz"), 10, 25, 15,
                         key="mont_kanat", width=240)
        plan_m = am.mutagenesis(kalip, int(konum), yerine, flank=int(kanat))
    else:
        st.markdown(t("**Birleştirilecek parçalar** — sıra montaj sırasıdır"))
        katalog = [x for x in db.list_parts() if x["seq_status"] == "full"]
        secim = st.multiselect(
            t("Katalogdan parça seç"),
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
                t("Inserti de ekle ({0}, {1} bp)").format(insert_ad, len(insert_seq)),
                value=bool(insert_seq), key="mont_insert"):
            parca_listesi.append((insert_ad[:14], insert_seq))

        o = st.container(horizontal=True, vertical_alignment="bottom")
        dairesel = o.toggle(t("Dairesel ürün"), value=True, key="mont_dairesel")
        if yontem == "Golden Gate":
            enz = o.selectbox(t("Tip IIS enzim"), am.TYPE_IIS_COMMON,
                              key="mont_enz", width=200)
            plan_m = am.golden_gate(parca_listesi, enz, circular=dairesel)
        else:
            ort = o.slider(t("Örtüşme uzunluğu (bp)"), 15, 45, 25,
                           key="mont_ort", width=260)
            plan_m = am.gibson(parca_listesi, overlap=int(ort),
                               circular=dairesel)

    if not plan_m.parts and yontem != "Mutagenez":
        st.info(t("Önce parça seçin."), icon=":material/category:")
    else:
        mm = st.container(horizontal=True)
        mm.metric(t("Yöntem"), plan_m.method)
        mm.metric(t("Parça"), len(plan_m.parts))
        mm.metric(t("Ürün"), f"{plan_m.product_length} bp")
        mm.metric(t("Birleşme"), len(plan_m.junctions))
        mm.metric(t("Durum"), "geçti" if plan_m.ok else "sorunlu")

        for n in plan_m.notes:
            {"hata": st.error, "uyari": st.warning,
             "bilgi": st.info}.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}")

        if plan_m.junctions:
            st.markdown(t("**Birleşme noktaları**"))
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
            st.markdown(t("**Sipariş edilecek primerler**"))
            st.dataframe(
                [{"Ad": p2.name, "Dizi": p2.sequence, "Uzunluk": p2.length,
                  "Bağlanma Tm": p2.tm_anneal,
                  "Uyarı": "; ".join(p2.warnings) if p2.warnings else "—"}
                 for p2 in plan_m.primers],
                width="stretch", hide_index=True,
                alt="Montaj için tasarlanan primerler",
            )
            st.download_button(
                t("Primerleri CSV olarak indir"),
                "ad,dizi,uzunluk,tm\n" + "\n".join(
                    f"{p2.name},{p2.sequence},{p2.length},{p2.tm_anneal}"
                    for p2 in plan_m.primers),
                file_name=f"{ad}-primerler.csv", mime="text/csv",
                icon=":material/download:",
            )
            st.caption(
                t("Primerlerin 5' ucundaki kuyruk, montajın kendisini "
                "sağlayan kısımdır; Tm sütunu yalnızca kalıba bağlanan "
                "bölümün değeridir, PCR koşulunu ona göre kurun.")
            )


# ---------------------------------------------------------------------------
# 6. Dizileme doğrulaması
# ---------------------------------------------------------------------------

with tab_dogrula:
    alt_plan, alt_oku = st.tabs([t("Dizileme planı"), t("Okumayı karşılaştır")])

    with alt_plan:
        st.caption(
            t("Sanger okuması ~900 nükleotiddir ve primerin hemen ardındaki "
            "ilk 50 baz güvenilir okunmaz. Tek primerle uzun bir yapıyı "
            "okutup 'doğruladım' demek, aradaki bölgeyi hiç görmemek "
            "demektir. Burası baştan sona kapatmak için kaç primer "
            "gerektiğini hesaplar.")
        )
        d = st.container(horizontal=True, vertical_alignment="bottom")
        okuma_uz = d.slider(t("Okuma uzunluğu (nt)"), 400, 1200,
                            vf.READ_LENGTH, 50, key="vf_read", width=240)
        olu = d.slider(t("Primer sonrası okunmayan (nt)"), 0, 120,
                       vf.DEAD_ZONE, 10, key="vf_dead", width=260)
        iki_zincir = d.toggle(t("İki zincirden de okut"), value=False,
                              key="vf_both")

        splan = vf.sequencing_plan(seq, read_length=int(okuma_uz),
                                   dead_zone=int(olu),
                                   both_strands=iki_zincir)
        sm = st.container(horizontal=True)
        sm.metric(t("Primer"), len(splan.primers))
        sm.metric(t("Kapsanan"), f"{splan.covered} / {splan.length} nt")
        sm.metric(t("Kapsama"), f"%{round(100 * splan.covered / max(splan.length, 1), 1)}")
        sm.metric(t("Boşluk"), len(splan.gaps))

        for n in splan.notes:
            {"hata": st.error, "uyari": st.warning,
             "bilgi": st.info}.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}")

        if splan.primers:
            st.dataframe(
                [{"Ad": p3.name, "Dizi": p3.sequence,
                  "Konum": p3.position,
                  "Yön": "+" if p3.strand > 0 else "−",
                  "Tm": round(p3.tm, 1) if p3.tm else None,
                  "GC": f"%{p3.gc}",
                  "Okur": f"{p3.covers_from}–{p3.covers_to}",
                  "Uyarı": "; ".join(p3.warnings) if p3.warnings else "—"}
                 for p3 in splan.primers],
                width="stretch", hide_index=True,
                alt="Dizileme primerleri ve kapsadıkları bölgeler",
            )
            st.download_button(
                t("Dizileme primerlerini indir"),
                "ad,dizi,konum,yon,tm\n" + "\n".join(
                    f"{p3.name},{p3.sequence},{p3.position},"
                    f"{'+' if p3.strand > 0 else '-'},{p3.tm}"
                    for p3 in splan.primers),
                file_name=f"{ad}-dizileme-primerleri.csv", mime="text/csv",
                icon=":material/download:",
            )

    with alt_oku:
        st.caption(
            t("Dizileme sonucunu yapıştırın: tasarımla hizalanır ve her "
            "fark tek tek çıkarılır. Kodlayan bölgedeki farkların amino "
            "asit karşılığı da hesaplanır - sessiz mi, değiştiren mi, "
            "yoksa erken durdurma mı.")
        )
        ab1 = st.file_uploader(
            t("Kromatogram dosyası (.ab1) — isteğe bağlı"),
            type=["ab1"], key="vf_ab1",
            help=t("Firmanın gönderdiği kromatogramı yüklerseniz kalite "
                 "değerleri okunur ve güvenilmez uçlar kırpılır. Metin "
                 "dizi kalite bilgisi taşımaz: uçlardaki bir 'fark' "
                 "gerçek mutasyon değil, okuma hatası olabilir."),
        )
        krom = None
        if ab1 is not None:
            krom = vf.read_chromatogram(ab1.getvalue())
            if krom is None:
                st.error(t("Dosya .ab1 olarak okunamadı."),
                         icon=":material/error:")
            else:
                q = st.container(horizontal=True)
                q.metric(t("Okuma"), f"{krom.length} baz")
                q.metric(t("Güvenilir bölge"),
                         f"{krom.trim_start}–{krom.trim_end}")
                q.metric(t("Ortalama kalite"), f"Q{krom.mean_quality}")
                for n in krom.notes:
                    {"hata": st.error, "uyari": st.warning,
                     "bilgi": st.info}.get(n["kind"], st.info)(
                        f"**{n['title']}** — {n['detail']}")
                st.line_chart(
                    {"kalite": krom.qualities}, height=160,
                    alt="Okuma boyunca baz çağırma kalitesi")

        okuma = st.text_area(
            t("Okunan dizi (FASTA ya da düz)"), height=130,
            value=(krom.trimmed_sequence if krom else ""),
            placeholder=t("Dizileme firmasından gelen diziyi yapıştırın "
                        "ya da yukarıdan .ab1 yükleyin…"),
            key="vf_okuma",
        )
        k2 = st.container(horizontal=True, vertical_alignment="bottom")
        kirp = k2.number_input(t("Baştan atılacak baz"), 0, 200, 0,
                               key="vf_trim", width=220,
                               help=t("Kromatogramın başındaki güvenilmez "
                                    "bölgeyi atmak için."))
        cds_var = k2.toggle(t("Kodlayan bölge var"), value=bool(insert_seq),
                            key="vf_cds")
        cds_bas = None
        if cds_var:
            cds_bas = k2.number_input(
                t("Başlama kodonunun konumu"), 1, max(len(seq), 1),
                min(1, len(seq)) or 1, key="vf_cdsbas", width=250)

        if clean_seq(okuma):
            rap2 = vf.compare_read(seq, okuma, cds_start=cds_bas,
                                   trim=int(kirp))
            om = st.container(horizontal=True)
            om.metric(t("Benzerlik"), f"%{rap2.identity}")
            om.metric(t("Hizalanan"), f"{rap2.aligned_from}–{rap2.aligned_to}")
            om.metric(t("Okuma"), f"{rap2.read_length} nt")
            om.metric(t("Fark"), len(rap2.differences))

            for n in rap2.notes:
                {"hata": st.error, "uyari": st.warning,
                 "bilgi": st.info}.get(n["kind"], st.info)(
                    f"**{n['title']}** — {n['detail']}")

            with st.expander(t("Hizalamayı göster"), expanded=bool(
                    rap2.differences)):
                st.code(
                    vf.alignment_view(seq, okuma, cds_start=cds_bas,
                                      trim=int(kirp)),
                    language="text")

            if rap2.differences:
                st.dataframe(
                    [{"Konum": f2.position,
                      "Tür": {"degisim": "değişim", "eksik": "eksik baz",
                              "fazla": "fazla baz"}[f2.kind],
                      "Tasarımda": f2.expected, "Okumada": f2.observed,
                      "Kodon": f2.codon_no or "—",
                      "Protein": f2.aa_change or "—"}
                     for f2 in rap2.differences],
                    width="stretch", hide_index=True,
                    alt="Okuma ile tasarım arasındaki farklar",
                )
        else:
            st.info(t("Karşılaştırmak için bir okuma yapıştırın."),
                    icon=":material/content_paste:")


# ---------------------------------------------------------------------------
# 7. Sorun giderme
# ---------------------------------------------------------------------------

with tab_sorun:
    alt_kesim, alt_belirti = st.tabs([t("Çift kesim uyumluluğu"),
                                      t("Belirtiden yola çık")])

    with alt_kesim:
        st.caption(
            t("İki enzimi aynı tüpte kullanmak her zaman mümkün değildir. "
            "Doğrulanabilir olanlar burada denetleniyor; tampon seçimi "
            "üreticiye özgü ölçüm verisi olduğu için uydurulmuyor.")
        )
        cift = st.multiselect(
            t("Birlikte kullanılacak enzimler"), _ticari(),
            default=(secilen or _tek_kesenler(seq))[:2],
            key="ts_cift",
        )
        for c in ts.double_digest(cift):
            {"hata": st.error, "uyari": st.warning,
             "bilgi": st.info}.get(c.kind, st.info)(
                f"**{c.title}** — {c.detail}")

    with alt_belirti:
        st.caption(
            t("Belirtiyi seçin: olası sebepler sıralanır. Uygulama "
            "tasarımınıza bakarak hangi sebebin sizde gerçekten var "
            "olduğunu işaretler - geri kalanı sizin kontrol etmeniz "
            "gereken şeylerdir.")
        )
        belirtiler = ts.symptoms()
        secili_b = st.selectbox(
            t("Ne oldu?"), [b.key for b in belirtiler],
            format_func=lambda k: next(b.label for b in belirtiler
                                       if b.key == k),
            key="ts_belirti",
        )
        sonuc = ts.evaluate(
            secili_b,
            part_ids=[x.get("part_id", "") for x in parcalar],
            host=db.get_host(st.session_state.get("design_host", "")),
            enzymes=secilen,
            purpose=st.session_state.get("design_purpose", "ifade"),
        )
        var = [c for c in sonuc.causes if c.applies is True]
        yok = [c for c in sonuc.causes if c.applies is False]
        belirsiz = [c for c in sonuc.causes if c.applies is None]

        if var:
            st.markdown(t("### Sizin tasarımınızda var"))
            for c in var:
                st.error(f"**{c.title}** — {c.detail}"
                         + (t("\n\n*Bakılacak:* {0}").format(c.check) if c.check else ""),
                         icon=":material/error:")
        if belirsiz:
            st.markdown(t("### Kontrol etmeniz gerekenler"))
            for c in belirsiz:
                with st.container(border=True):
                    st.markdown(f"**{c.title}**")
                    st.markdown(c.detail)
                    if c.check:
                        st.caption(t("Bakılacak: {0}").format(c.check))
        if yok:
            with st.expander(t("Elenen sebepler ({0})").format(len(yok))):
                for c in yok:
                    st.markdown(t("- **{0}** — tasarımınızda bu durum yok, bu sebebi eleyebilirsiniz.").format(c.title))
        if not sonuc.causes:
            st.info(t("Bu belirti için kayıt yok."), icon=":material/help:")


# ---------------------------------------------------------------------------
# 6. Rapor
# ---------------------------------------------------------------------------

with tab_rap:
    st.caption(
        t("Haritayı, parça listesini, denetim sonucunu, protein özelliklerini "
        "ve protokolü tek bir HTML dosyasında toplar. Tarayıcıda açıp "
        "oradan PDF'e basabilirsiniz.")
    )
    r = st.container(horizontal=True, vertical_alignment="bottom")
    yazar = r.text_input(t("Rapora yazılacak isim (isteğe bağlı)"), width=300)
    icerik = r.multiselect(
        t("Rapora girecek bölümler"),
        ["Harita", "Parçalar", "Denetim", "Protein", "Jel", "Protokol"],
        default=["Harita", "Parçalar", "Denetim", "Protokol"],
    format_func=t)

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
        t("Raporu HTML olarak indir"), html_rapor,
        file_name=f"{ad}-rapor.html", mime="text/html",
        icon=":material/download:", type="primary",
    )
    st.caption(t("Rapor {0} KB · harita ve jel dosyanın içine gömülür, ayrı görsel gerekmez.").format(len(html_rapor) // 1024))

    with st.expander(t("Raporu burada önizle")):
        # st.iframe kaynak olarak adres ister, ham HTML almaz; raporu
        # data: adresine gömüp öyle veriyoruz - böylece dosyayı diske
        # yazmadan ve kaldırılmış bir API'ye dayanmadan gösterilebiliyor.
        src = ("data:text/html;charset=utf-8;base64,"
               + base64.b64encode(html_rapor.encode("utf-8")).decode("ascii"))
        st.iframe(src, height=620, alt=t("{0} raporunun önizlemesi").format(ad))
