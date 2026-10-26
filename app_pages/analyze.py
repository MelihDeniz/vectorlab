"""Dizi analizi: anotasyon, enzim haritası, ORF, primer, sindirim."""

import streamlit as st

from core import db, enzymes as ez
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
    seq_stats,
    translate_cds,
)

src = st.container(horizontal=True, vertical_alignment="bottom")
with src:
    uploaded = st.file_uploader(
        "FASTA ya da GenBank dosyası", type=["fa", "fasta", "txt", "gb", "gbk", "seq"]
    )
    circular = st.toggle("Dairesel", value=True)

if uploaded is not None:
    raw = uploaded.getvalue().decode("utf-8", errors="ignore")
    if "ORIGIN" in raw:  # kaba GenBank ayrıştırma: ORIGIN sonrası dizi
        raw = raw.split("ORIGIN", 1)[1]
    st.session_state.analyze_seq = raw

text = st.text_area(
    "Dizi",
    value=st.session_state.analyze_seq,
    height=160,
    placeholder="Diziyi yapıştırın ya da yukarıdan dosya yükleyin…",
)
st.session_state.analyze_seq = text
seq = clean_seq(text)

if not seq:
    st.info(
        "Analiz için bir dizi girin. Tasarımcı sayfasında kurduğunuz yapının "
        "FASTA çıktısını buraya yapıştırabilirsiniz.",
        icon=":material/science:",
    )
    st.stop()

stats = seq_stats(seq, circular)
m = st.container(horizontal=True)
m.metric("Uzunluk", f"{stats['length']} bp")
m.metric("GC", f"%{stats['gc']}")
m.metric("Kütle", f"{stats['mw_kda']} kDa")
m.metric("Topoloji", stats["topology"])
if stats["ambiguous"]:
    m.metric("Belirsiz baz", stats["ambiguous"])

tabs = st.tabs(
    ["Anotasyon ve harita", "Enzim haritası", "Sindirim", "ORF ve çeviri",
     "Primer tasarımı"]
)

# ---------------------------------------------------------------------------
# Anotasyon
# ---------------------------------------------------------------------------

with tabs[0]:
    from core.seq_utils import annotate

    feats = annotate(seq, db.list_parts(), circular=circular)
    if not feats:
        st.warning(
            "Katalogdaki parçalardan hiçbiri bu dizide bulunamadı. "
            "Anotasyon yalnızca dizisi tam olan parçalara karşı yapılır; "
            "uzun genler (direnç, reporter) katalogda dizi taşımadığı için "
            "eşleşmez.",
            icon=":material/search_off:",
        )
    else:
        st.caption(f"{len(feats)} eşleşme bulundu.")
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
        "Taranacak enzimler",
        _tum,
        default=[e for e in ez.WORKHORSE if e in _tum][:24],
        help=f"REBASE kütüphanesindeki {len(_tum)} ticari enzimin hepsi "
             "aranabilir; varsayılan liste rutin klonlamada kullanılanlardır.",
    )
    groups = classify_cutters(seq, enz_set, circular)

    c1, c2 = st.columns(2)
    with c1:
        st.markdown("**Tek kesen** — klonlama için kullanılabilir")
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
            st.caption("Tek kesen enzim yok.")

        st.markdown("**Kesmeyen** — insertte güvenle kullanılabilir")
        st.write(", ".join(n for n, _ in groups["kesmez"]) or "—")

    with c2:
        st.markdown("**İki kesen**")
        st.write(
            ", ".join(f"{n} ({len(p)})" for n, p in groups["cift"]) or "—"
        )
        st.markdown("**Çok kesen**")
        st.write(
            ", ".join(f"{n} ({len(p)})" for n, p in groups["cok"]) or "—"
        )

        with st.expander("Metilasyona duyarlı enzimler hakkında"):
            st.markdown(
                "E. coli'den saflaştırılan plazmit Dam ve Dcm metilasyonu "
                "taşır. Aşağıdaki enzimler bu yüzden beklenen bölgeyi "
                "kesmeyebilir — sindirim boş çıkar ve nedeni kolay "
                "anlaşılmaz. Gerekirse plazmiti *dam⁻/dcm⁻* bir suşta "
                "(JM110, SCS110) çoğaltın."
            )
            from core.seq_utils import METHYLATION_BLOCKED
            for k, v in METHYLATION_BLOCKED.items():
                st.markdown(f"- **{k}** — {v}")

# ---------------------------------------------------------------------------
# Sindirim
# ---------------------------------------------------------------------------

with tabs[2]:
    dig_enz = st.multiselect(
        "Sindirim enzimleri", ez.search("", commercial_only=True),
        default=[n for n, _ in classify_cutters(seq, None, circular)["tek"]][:2],
    )
    if not dig_enz:
        st.caption("Jel simülasyonu için en az bir enzim seçin.")
    else:
        frags = digest(seq, dig_enz, circular)
        st.caption(
            f"{len(frags)} fragman. Agaroz jelde büyükten küçüğe doğru "
            "sıralanır; 100 bp altındaki bantlar çoğu jelde görünmez."
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
    min_aa = o1.number_input("En kısa ORF (aa)", 20, 500, 50, step=10)
    alt_start = o2.toggle("GTG/TTG başlangıçlarını da ara", value=False)
    orfs = find_orfs(seq, min_aa=int(min_aa), circular=circular,
                     alt_starts=alt_start)
    o3.metric("Bulunan ORF", len(orfs))

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
    st.markdown("**Seçili bölgeyi çevir**")
    cds_text = st.text_area(
        "CDS dizisi", height=100,
        placeholder="ATG ile başlayan kodlayıcı diziyi yapıştırın…",
    )
    if clean_seq(cds_text):
        res = translate_cds(cds_text)
        cc = st.container(horizontal=True)
        cc.metric("Aminoasit", res["aa_len"])
        cc.metric("Kütle", f"~{res['mw_kda']} kDa")
        for issue in res["issues"]:
            st.warning(issue, icon=":material/warning:")
        st.code(res["protein"], language="text")

        codons = codon_analysis(cds_text)
        if codons["total"]:
            st.markdown(
                f"**Nadir kodon:** {len(codons['rare'])}/{codons['total']} "
                f"(%{codons['rare_pct']}) · yan yana nadir çift: "
                f"{len(codons['tandem'])}"
            )
            if codons["rare_pct"] > 8:
                st.info(
                    "Nadir kodon oranı yüksek. Rosetta(DE3) konağı bu "
                    "kodonların tRNA'larını sağlar; alternatif olarak "
                    "diziyi E. coli için kodon optimize edin.",
                    icon=":material/lightbulb:",
                )

    sd = rbs_spacing(seq)
    if sd:
        st.divider()
        st.markdown("**Shine-Dalgarno boşluk kontrolü**")
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
        "Hedef bölgeyi çoğaltacak primer çifti tasarlar. Uçlara restriksiyon "
        "bölgesi eklerseniz enzimin verimli kesmesi için koruyucu bazlar da "
        "eklenir."
    )
    p1, p2, p3 = st.columns(3)
    region = p1.text_input("Hedef bölge (başlangıç-bitiş)", value=f"1-{min(len(seq), 900)}")
    _primer_enz = ez.search("", commercial_only=True)
    f_enz = p2.selectbox("5' enzim", ["— yok —"] + _primer_enz)
    r_enz = p3.selectbox("3' enzim", ["— yok —"] + _primer_enz, index=0)
    tm_target = st.slider("Hedef Tm (°C)", 52.0, 68.0, 60.0, 0.5)

    try:
        a, b = (int(x) for x in region.replace(" ", "").split("-"))
        target = seq[a - 1 : b]
    except ValueError:
        st.error("Bölgeyi `başlangıç-bitiş` biçiminde yazın, örn. `1-900`.")
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
            cc.metric("Ürün uzunluğu", f"{res['product_length']} bp")
            cc.metric("Önerilen bağlanma sıcaklığı",
                      f"{res['annealing_temp']} °C"
                      if res["annealing_temp"] else "—")
            for p in (res["forward"], res["reverse"]):
                for w in p.warnings:
                    st.warning(f"{p.name}: {w}", icon=":material/warning:")
            for w in res["pair_warnings"]:
                st.error(w, icon=":material/error:")
