"""Dizi analizi: anotasyon, enzim haritası, ORF, primer, sindirim."""

import streamlit as st

from core import codon, db, enzymes as ez, structure as stx
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import circular_map, linear_map, render as render_svg
from core.seq_utils import (
    classify_cutters,
    clean_seq,
    codon_analysis,
    design_primers,
    digest,
    enzyme_info,
    find_orfs,
    gc_percent,
    rbs_spacing,
    restriction_sites,
    revcomp,
    seq_stats,
    translate_cds,
)
from core.i18n import t

src = st.container(horizontal=True, vertical_alignment="bottom")
with src:
    uploaded = st.file_uploader(
        t("FASTA ya da GenBank dosyası"), type=["fa", "fasta", "txt", "gb", "gbk", "seq"]
    )
    circular = st.toggle(t("Dairesel"), value=True)

if uploaded is not None:
    raw = uploaded.getvalue().decode("utf-8", errors="ignore")
    if "ORIGIN" in raw:  # kaba GenBank ayrıştırma: ORIGIN sonrası dizi
        raw = raw.split("ORIGIN", 1)[1]
    st.session_state.analyze_seq = raw

text = st.text_area(
    t("Dizi"),
    value=st.session_state.analyze_seq,
    height=160,
    placeholder=t("Diziyi yapıştırın ya da yukarıdan dosya yükleyin…"),
)
st.session_state.analyze_seq = text
seq = clean_seq(text)

if not seq:
    st.info(
        t("Analiz için bir dizi girin. Tasarımcı sayfasında kurduğunuz yapının "
        "FASTA çıktısını buraya yapıştırabilirsiniz."),
        icon=":material/science:",
    )
    st.stop()

stats = seq_stats(seq, circular)
m = st.container(horizontal=True)
m.metric(t("Uzunluk"), f"{stats['length']} bp")
m.metric(t("GC"), f"%{stats['gc']}")
m.metric(t("Kütle"), f"{stats['mw_kda']} kDa")
m.metric(t("Topoloji"), stats["topology"])
if stats["ambiguous"]:
    m.metric(t("Belirsiz baz"), stats["ambiguous"])

tabs = st.tabs(
    [t("Anotasyon ve harita"), t("Enzim haritası"), t("Sindirim"), t("ORF ve çeviri"),
     t("Primer tasarımı"), t("Kodon optimizasyonu"), t("mRNA yapısı")]
)

# ---------------------------------------------------------------------------
# Anotasyon
# ---------------------------------------------------------------------------

with tabs[0]:
    from core.seq_utils import annotate

    feats = annotate(seq, db.list_parts(), circular=circular)
    if not feats:
        st.warning(
            t("Katalogdaki parçalardan hiçbiri bu dizide bulunamadı. "
            "Anotasyon yalnızca dizisi tam olan parçalara karşı yapılır; "
            "uzun genler (direnç, reporter) katalogda dizi taşımadığı için "
            "eşleşmez."),
            icon=":material/search_off:",
        )
    else:
        st.caption(t("{0} eşleşme bulundu.").format(len(feats)))
        render_svg(
            circular_map("Analiz edilen dizi", len(seq), feats)
            if circular
            else linear_map("Analiz edilen dizi", len(seq), feats)
        )
        st.dataframe(
            [
                {
                    "Başlangıç": f.start, "Bitiş": f.end,
                    "Zincir": "+" if f.strand > 0 else "−",
                    "Uzunluk": f.length,
                    "Tür": CATEGORY_TR.get(f.category, f.category),
                    "Parça": f.name,
                }
                for f in feats
            ],
            width="stretch", hide_index=True,
            alt="Dizide bulunan katalog parçalarının listesi",
        )

# ---------------------------------------------------------------------------
# Enzim haritası
# ---------------------------------------------------------------------------

with tabs[1]:
    _tum = ez.search("", commercial_only=True)
    enz_set = st.multiselect(
        t("Taranacak enzimler"),
        _tum,
        default=[e for e in ez.WORKHORSE if e in _tum][:24],
        help=t("REBASE kütüphanesindeki {0} ticari enzimin hepsi aranabilir; varsayılan liste rutin klonlamada kullanılanlardır.").format(len(_tum)),
    )
    groups = classify_cutters(seq, enz_set, circular)

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(t("**Tek kesen** — klonlama için kullanılabilir"))
        if groups["tek"]:
            st.dataframe(
                [{"Enzim": n, "Konum": p[0],
                  "Uç": enzyme_info(n).get("overhang", ""),
                  "Tanıma": enzyme_info(n).get("site", "")}
                 for n, p in groups["tek"]],
                width="stretch", hide_index=True,
                alt="Diziyi bir kez kesen enzimler",
            )
        else:
            st.caption(t("Tek kesen enzim yok."))

        st.markdown(t("**Kesmeyen** — insertte güvenle kullanılabilir"))
        st.write(", ".join(n for n, _ in groups["kesmez"]) or "—")

    with c2:
        st.markdown(t("**İki kesen**"))
        st.write(
            ", ".join(f"{n} ({len(p)})" for n, p in groups["cift"]) or "—"
        )
        st.markdown(t("**Çok kesen**"))
        st.write(
            ", ".join(f"{n} ({len(p)})" for n, p in groups["cok"]) or "—"
        )

        with st.expander(t("Metilasyona duyarlı enzimler hakkında")):
            st.markdown(
                t("E. coli'den saflaştırılan plazmit Dam ve Dcm metilasyonu "
                "taşır. Aşağıdaki enzimler bu yüzden beklenen bölgeyi "
                "kesmeyebilir — sindirim boş çıkar ve nedeni kolay "
                "anlaşılmaz. Gerekirse plazmiti *dam⁻/dcm⁻* bir suşta "
                "(JM110, SCS110) çoğaltın.")
            )
            from core.seq_utils import METHYLATION_BLOCKED
            for k, v in METHYLATION_BLOCKED.items():
                st.markdown(f"- **{k}** — {v}")

# ---------------------------------------------------------------------------
# Sindirim
# ---------------------------------------------------------------------------

with tabs[2]:
    dig_enz = st.multiselect(
        t("Sindirim enzimleri"), ez.search("", commercial_only=True),
        default=[n for n, _ in classify_cutters(seq, None, circular)["tek"]][:2],
    )
    if not dig_enz:
        st.caption(t("Jel simülasyonu için en az bir enzim seçin."))
    else:
        frags = digest(seq, dig_enz, circular)
        st.caption(
            t("{0} fragman. Agaroz jelde büyükten küçüğe doğru sıralanır; 100 bp altındaki bantlar çoğu jelde görünmez.").format(len(frags))
        )
        st.dataframe(
            [{"Fragman": i + 1, "Uzunluk (bp)": f["length"],
              "Başlangıç": f["start"], "Bitiş": f["end"],
              "Not": f["note"] or ("jelde görünmeyebilir"
                                   if f["length"] < 100 else "")}
             for i, f in enumerate(frags)],
            width="stretch", hide_index=True,
            alt="Sindirim sonucu oluşan fragmanların uzunlukları",
        )
        st.bar_chart(
            {"Uzunluk (bp)": [f["length"] for f in frags]},
            height=220,
            alt="Fragman uzunluklarının çubuk grafiği",
        )

# ---------------------------------------------------------------------------
# ORF ve çeviri
# ---------------------------------------------------------------------------

with tabs[3]:
    o1, o2, o3 = st.columns([1, 1, 2])
    min_aa = o1.number_input(t("En kısa ORF (aa)"), 20, 500, 50, step=10)
    alt_start = o2.toggle(t("GTG/TTG başlangıçlarını da ara"), value=False)
    orfs = find_orfs(seq, min_aa=int(min_aa), circular=circular,
                     alt_starts=alt_start)
    o3.metric(t("Bulunan ORF"), len(orfs))

    if orfs:
        st.dataframe(
            [{"Başlangıç": o.start, "Bitiş": o.end,
              "Zincir": "+" if o.strand > 0 else "−",
              "Uzunluk (aa)": o.length_aa,
              "Başlama kodonu": o.start_codon,
              "Stop var": "evet" if o.has_stop else "hayır",
              "Protein (ilk 40 aa)": o.protein[:40]}
             for o in orfs],
            width="stretch", hide_index=True,
            alt="Dizide bulunan açık okuma çerçeveleri",
        )

    st.divider()
    st.markdown(t("**Seçili bölgeyi çevir**"))
    cds_text = st.text_area(
        t("CDS dizisi"), height=100,
        placeholder=t("ATG ile başlayan kodlayıcı diziyi yapıştırın…"),
    )
    if clean_seq(cds_text):
        res = translate_cds(cds_text)
        cc = st.container(horizontal=True)
        cc.metric(t("Aminoasit"), res["aa_len"])
        cc.metric(t("Kütle"), f"~{res['mw_kda']} kDa")
        for issue in res["issues"]:
            st.warning(issue, icon=":material/warning:")
        st.code(res["protein"], language="text")

        codons = codon_analysis(cds_text)
        if codons["total"]:
            st.markdown(
                t("**Nadir kodon:** {0}/{1} (%{2}) · yan yana nadir çift: {3}").format(len(codons['rare']), codons['total'], codons['rare_pct'], len(codons['tandem']))
            )
            if codons["rare_pct"] > 8:
                st.info(
                    t("Nadir kodon oranı yüksek. Rosetta(DE3) konağı bu "
                    "kodonların tRNA'larını sağlar; alternatif olarak "
                    "diziyi E. coli için kodon optimize edin."),
                    icon=":material/lightbulb:",
                )

    sd = rbs_spacing(seq)
    if sd:
        st.divider()
        st.markdown(t("**Shine-Dalgarno boşluk kontrolü**"))
        st.dataframe(
            [{"ATG konumu": s["atg_pos"], "SD motifi": s["motif"],
              "Boşluk (nt)": s["gap"], "Değerlendirme": s["verdict"]}
             for s in sd],
            width="stretch", hide_index=True,
            alt="Shine-Dalgarno motifi ile başlama kodonu arası boşluk",
        )

# ---------------------------------------------------------------------------
# Primer
# ---------------------------------------------------------------------------

with tabs[4]:
    st.caption(
        t("Hedef bölgeyi çoğaltacak primer çifti tasarlar. Uçlara restriksiyon "
        "bölgesi eklerseniz enzimin verimli kesmesi için koruyucu bazlar da "
        "eklenir.")
    )
    p1, p2, p3 = st.columns(3)
    region = p1.text_input(t("Hedef bölge (başlangıç-bitiş)"), value=f"1-{min(len(seq), 900)}")
    _primer_enz = ez.search("", commercial_only=True)
    f_enz = p2.selectbox(t("5' enzim"), ["— yok —"] + _primer_enz)
    r_enz = p3.selectbox(t("3' enzim"), ["— yok —"] + _primer_enz, index=0)
    tm_target = st.slider(t("Hedef Tm (°C)"), 52.0, 68.0, 60.0, 0.5)

    try:
        a, b = (int(x) for x in region.replace(" ", "").split("-"))
        target = seq[a - 1 : b]
    except ValueError:
        st.error(t("Bölgeyi `başlangıç-bitiş` biçiminde yazın, örn. `1-900`."))
        target = ""

    if target:
        res = design_primers(
            target,
            fwd_enzyme=None if f_enz.startswith("—") else f_enz,
            rev_enzyme=None if r_enz.startswith("—") else r_enz,
            target_tm=tm_target,
        )
        if "error" in res:
            st.error(res["error"])
        else:
            st.dataframe(
                [
                    {"Primer": p.name, "Dizi (5'→3')": p.sequence,
                     "Uzunluk": p.length, "Tm (°C)": p.tm, "GC %": p.gc,
                     "Kuyruk": p.tail or "—"}
                    for p in (res["forward"], res["reverse"])
                ],
                width="stretch", hide_index=True,
                alt="Tasarlanan ileri ve geri primerler",
            )
            cc = st.container(horizontal=True)
            cc.metric(t("Ürün uzunluğu"), f"{res['product_length']} bp")
            cc.metric(t("Önerilen bağlanma sıcaklığı"),
                      f"{res['annealing_temp']} °C"
                      if res["annealing_temp"] else "—")
            for p in (res["forward"], res["reverse"]):
                for w in p.warnings:
                    st.warning(f"{p.name}: {w}", icon=":material/warning:")
            for w in res["pair_warnings"]:
                st.error(w, icon=":material/error:")


# ---------------------------------------------------------------------------
# Kodon optimizasyonu
# ---------------------------------------------------------------------------

with tabs[5]:
    st.caption(
        t("Başka bir organizmadan gelen gen, konağın az bulundurduğu "
        "tRNA'lara denk gelen kodonlar taşıyorsa ribozom o noktalarda "
        "duraklar. Burası bunu ölçer ve istenirse düzeltir.")
    )

    # Hangi organizmaya göre optimize edildiği KULLANICININ SEÇİMİ olmak
    # zorunda. Tek bir E. coli tablosuyla çalışmak, L. cremoris'te
    # tasarlanan bir yapıyı sessizce yanlış yöne çekmek demekti.
    _org_etiket = {}
    for _a in codon.available():
        _t = codon.table(_a)
        _org_etiket[_a] = _t.get("organism", _a)

    # Varsayılan: tasarımın konağının organizması.
    _konak_id = st.session_state.get("design_host", "")
    _konak = next((h for h in db.list_hosts()
                   if h["host_id"] == _konak_id), None)
    _konak_org = (_konak or {}).get("organism", "")
    _varsayilan = codon.table_key_for_organism(_konak_org) or codon.DEFAULT_ORG

    _secenekler = sorted(_org_etiket, key=lambda a: _org_etiket[a])
    org = st.selectbox(
        t("Hangi organizmaya göre optimize edilsin?"),
        _secenekler,
        index=_secenekler.index(_varsayilan),
        format_func=lambda a: _org_etiket[a],
        key="kodon_org",
        help=t("Kodon tercihi organizmaya göre değişir. Tasarımınızın "
             "konağı neyse onu seçin - yanlış tabloyla optimize etmek, "
             "optimize etmemekten daha iyi değildir."),
    )
    if _konak_org and codon.table_key_for_organism(_konak_org) is None:
        st.warning(
            t("Tasarımın konağı **{0}** ve bu organizma için kodon tablomuz yok. Aşağıdaki seçim onun yerine kullanılacak; sonucu buna göre okuyun.").format(_konak_org),
            icon=":material/warning:",
        )
    elif _konak and codon.table_key_for_organism(_konak_org) != org:
        st.info(
            t("Tasarımın konağı **{0}** ({1}), ama optimizasyon **{2}** tablosuna göre yapılacak.").format(_konak['name'], _konak_org, _org_etiket[org]),
            icon=":material/info:",
        )
    st.caption(codon.source_note(org))

    if not seq:
        st.info(t("Önce bir dizi girin."), icon=":material/science:")
    else:
        kodon_kaynak = st.radio(
            t("Hangi bölge?"), ["Tüm dizi", "Seçtiğim açık okuma çerçevesi"],
            horizontal=True, key="kodon_kaynak",
        format_func=t)
        cds = seq
        if kodon_kaynak != "Tüm dizi":
            orf_listesi = find_orfs(seq, min_aa=40, circular=circular)
            if not orf_listesi:
                st.warning(t("40 aa'dan uzun açık okuma çerçevesi yok."),
                           icon=":material/warning:")
                cds = ""
            else:
                sec = st.selectbox(t("Çerçeve"), orf_listesi,
                                   format_func=lambda o: o.label(),
                                   key="kodon_orf")
                parca = seq[sec.start - 1:sec.end]
                cds = parca if sec.strand > 0 else revcomp(parca)

        if cds:
            mevcut_cai = codon.cai(cds, org=org)
            nadir = codon.rare_positions(cds, org=org)
            kume = codon.rare_clusters(cds, org=org)
            m = st.container(horizontal=True)
            m.metric(t("CAI"), mevcut_cai if mevcut_cai else "—")
            m.metric(t("Nadir kodon"), len(nadir))
            m.metric(t("Nadir küme"), len(kume))
            m.metric(t("GC"), f"%{gc_percent(cds)}")

            if kume:
                st.warning(
                    "Nadir kodonların kümelendiği yerler: "
                    + ", ".join(t("{0}–{1}. kodon ({2} adet)").format(k['start'], k['end'], k['count']) for k in kume[:6])
                    + ". Arka arkaya gelen nadir kodonlar, dağınık "
                      "olanlardan daha çok sorun çıkarır: ribozom aynı "
                      "noktada birikir.",
                    icon=":material/warning:",
                )

            st.divider()
            o = st.container(horizontal=True, vertical_alignment="bottom")
            kip = o.segmented_control(
                t("Kip"), ["nadir", "tam"],
                format_func=lambda k: ("Sadece nadir kodonlar"
                                       if k == "nadir" else "Tam optimizasyon"),
                default="nadir", key="kodon_kip",
            ) or "nadir"
            yasak = o.multiselect(
                t("İnsertten temizlenecek enzim bölgeleri"),
                ez.search("", commercial_only=True),
                default=[e for e in st.session_state.get("design_enzymes", [])
                         if e in ez.search("", commercial_only=True)][:2],
                help=t("Klonlamada kullanacağınız enzimleri seçin; "
                     "optimizasyon bu bölgeleri dizinin içinden kaldırır."),
                key="kodon_yasak",
            )

            siteler = []
            for e in yasak:
                d = ez.details(e)
                if d and d.get("site"):
                    siteler.append(d["site"])

            if st.button(t("Optimize et"), icon=":material/auto_fix_high:",
                         type="primary", key="kodon_calistir"):
                st.session_state["kodon_sonuc"] = codon.optimize(
                    cds, mode=kip, forbidden_sites=siteler, org=org)

            res = st.session_state.get("kodon_sonuc")
            if res and res.sequence:
                k = st.container(horizontal=True)
                k.metric(t("CAI"), res.cai_after,
                         delta=round(res.cai_after - res.cai_before, 4))
                k.metric(t("Nadir kodon"), res.rare_after,
                         delta=res.rare_after - res.rare_before,
                         delta_color="inverse")
                k.metric(t("GC"), f"%{res.gc_after}",
                         delta=round(res.gc_after - res.gc_before, 2))
                k.metric(t("Değişen kodon"), res.changed)

                for n in codon.interpret(res, org=org):
                    {"hata": st.error, "uyari": st.warning,
                     "bilgi": st.info}.get(n["kind"], st.info)(
                        f"**{n['title']}** — {n['detail']}")
                for n in res.notes:
                    st.warning(n, icon=":material/info:")

                st.code(res.sequence, language="text", wrap_lines=True)
                st.download_button(
                    t("Optimize diziyi indir"),
                    f">optimize_{org}_CAI{res.cai_after} "
                    f"[{_org_etiket[org]}]\n" + "\n".join(
                        res.sequence[i:i + 70]
                        for i in range(0, len(res.sequence), 70)) + "\n",
                    file_name="optimize.fasta", mime="text/plain",
                    icon=":material/download:",
                )
                if res.changes:
                    with st.expander(t("Değişen {0} kodon").format(len(res.changes))):
                        st.dataframe(
                            [{"Kodon no": c2.position, "Amino asit": c2.aa,
                              "Önce": c2.before, "Sonra": c2.after,
                              "w önce": c2.w_before, "w sonra": c2.w_after,
                              "Sebep": c2.reason} for c2 in res.changes],
                            width="stretch", hide_index=True,
                            alt="Optimizasyonda değişen kodonların listesi",
                        )


# ---------------------------------------------------------------------------
# mRNA ikincil yapısı
# ---------------------------------------------------------------------------

with tabs[6]:
    st.caption(
        t("Ribozom, Shine-Dalgarno dizisinden tutunup başlama kodonuna "
        "kayar. Bu bölge kendi üzerine katlanıp saç tokası yaparsa "
        "bağlanamaz: promotör güçlü ve kodonlar doğru olduğu halde "
        "protein alamamanın en sinsi sebebi budur, çünkü haritada "
        "görünmez.")
    )
    if not seq:
        st.info(t("Önce bir dizi girin."), icon=":material/science:")
    else:
        y = st.container(horizontal=True, vertical_alignment="bottom")
        tahmin = stx.guess_start(seq)
        elle = y.toggle(t("Başlama kodonunu elle gir"), value=False,
                        key="stx_elle")
        if elle:
            pos = y.number_input(t("Başlama kodonu konumu"), 1, max(len(seq), 1),
                                 tahmin or 1, key="stx_pos", width=220)
        else:
            pos = tahmin
            y.caption(t("Bulunan başlama kodonu: {0}. nükleotid").format(pos if pos else 'yok'))
        yukari = y.slider(t("Yukarı akış penceresi (nt)"), 20, 80, 45,
                          key="stx_up", width=240)

        rap = stx.analyze_start_region(seq, pos, upstream=yukari)
        if rap.start_position:
            s1 = st.container(horizontal=True)
            s1.metric(t("Başlama kodonu"), f"{rap.start_position}. nt")
            s1.metric(t("Shine-Dalgarno"),
                      f"{rap.sd_position}. nt" if rap.sd_position else "yok")
            if rap.sd_position:
                s1.metric(t("Aralık"),
                          f"{rap.start_position - rap.sd_position - 6} nt")
            s1.metric(t("Saç tokası"), len(rap.hairpins))
            if rap.strongest:
                s1.metric(t("En kararlı"), f"{rap.strongest.dg} kcal/mol")

        for n in rap.notes:
            {"hata": st.error, "uyari": st.warning,
             "bilgi": st.info}.get(n["kind"], st.info)(
                f"**{n['title']}** — {n['detail']}")

        if rap.hairpins:
            st.markdown(t("**Bulunan yapılar**"))
            st.dataframe(
                [{"Başlangıç": rap.window_start + h.stem_start - 1,
                  "Bitiş": rap.window_start + h.stem_end - 1,
                  "Sap (bç)": h.stem_len, "İlmek (nt)": h.loop_len,
                  "ΔG (kcal/mol)": h.dg, "Dizi": h.sequence}
                 for h in rap.hairpins[:12]],
                width="stretch", hide_index=True,
                alt="Başlama bölgesinde bulunan saç tokaları",
            )
            st.caption(stx.dg_scale_note())
        st.info(
            t("**Bu tam bir katlanma hesabı değildir.** Burada yalnızca tek "
            "tek saç tokaları taranır; ViennaRNA ya da mfold gibi "
            "programlar bütün olası yapıları birlikte değerlendirip en "
            "düşük enerjili katlanmayı bulur. Çıkan sayıyı mutlak bir "
            "ΔG olarak değil, 'bu bölgede güçlü bir yapı var mı' "
            "sorusunun göreli yanıtı olarak okuyun."),
            icon=":material/info:",
        )
