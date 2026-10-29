"""Tasarımcı: parçalardan vektör kur, haritasını çiz, tasarımı denetle."""

import streamlit as st

from core import db, enzymes as ez, library as lib
from core.assembly import assemble, sort_canonical, to_fasta, to_genbank
from core.env import SHARED_NOTICE, is_shared_deployment
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import circular_map, linear_map, render as render_svg
from core.seq_utils import clean_seq, restriction_sites
from core.validate import summarize, validate

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

all_parts = db.list_parts()
parts_by_id = {p["part_id"]: p for p in all_parts}
hosts = db.list_hosts()

# ---------------------------------------------------------------------------
# Üst şerit: yapı kimliği
# ---------------------------------------------------------------------------

top = st.container(horizontal=True, vertical_alignment="bottom")
with top:
    st.session_state.design_name = st.text_input(
        "Yapı adı", value=st.session_state.design_name, width=220
    )
    host_ids = [h["host_id"] for h in hosts]
    host_labels = {h["host_id"]: f"{h['name']} — {h['use']}" for h in hosts}
    st.session_state.design_host = st.selectbox(
        "Konak suşu",
        host_ids,
        index=host_ids.index(st.session_state.design_host)
        if st.session_state.design_host in host_ids else 0,
        format_func=lambda h: host_labels[h],
        width=300,
    )
    st.session_state.design_purpose = st.segmented_control(
        "Amaç", ["ifade", "klonlama"],
        default=st.session_state.design_purpose,
    ) or st.session_state.design_purpose

host = db.get_host(st.session_state.design_host)
if host:
    st.caption(f"**Genotip:** `{host['genotype']}` · {host['notes']}")

st.divider()

left, right = st.columns([0.38, 0.62], gap="large")

# ---------------------------------------------------------------------------
# Sol: parça seçimi
# ---------------------------------------------------------------------------

with left:
    st.subheader("Parçalar")

    preset = st.selectbox(
        "Hazır şablondan başla",
        ["— seçim yok —", "pET tipi ifade vektörü", "pUC tipi klonlama vektörü",
         "Arabinozla indüklenen (P-BAD)", "Sabit ifade (sentetik devre)"],
    )
    presets = {
        "pET tipi ifade vektörü": (
            ["PRO_T7LAC", "RBS_PET", "MCS_PET28", "TER_T7", "ORI_COLE1",
             "MRK_KAN", "REG_LACI"], "BL21DE3", "ifade"),
        "pUC tipi klonlama vektörü": (
            ["PRO_LAC", "RPT_LACZA", "MCS_PUC19", "ORI_PUC", "MRK_BLA"],
            "DH5A", "klonlama"),
        "Arabinozla indüklenen (P-BAD)": (
            ["REG_ARAC", "PRO_ARABAD", "RBS_B0034", "TAG_HIS6", "PRS_TEV",
             "TER_B0015", "ORI_PBR322", "MRK_BLA"], "TOP10", "ifade"),
        "Sabit ifade (sentetik devre)": (
            ["PRO_J23105", "RBS_B0032", "RPT_SFGFP", "TER_B0015", "ORI_P15A",
             "MRK_CAT"], "DH5A", "ifade"),
    }
    # Kategori sırası; seçim durumu bu widget'ların kendisinde tutulur.
    CATS = ["device", "promoter", "operator", "rbs", "tag", "protease_site",
            "linker", "mcs", "reporter", "terminator", "ori", "marker",
            "regulator"]

    if preset in presets and st.button("Şablonu yükle", icon=":material/download:"):
        ids, h, purpose = presets[preset]
        # Widget key'leri doğrudan yazılmalı: key verilmiş bir multiselect
        # kendi durumunu korur ve `default` parametresini sonraki
        # çalışmalarda yok sayar.
        for cat in CATS:
            st.session_state[f"ms_{cat}"] = [
                pid for pid in ids
                if parts_by_id.get(pid, {}).get("category") == cat
            ]
        st.session_state.design_host = h
        st.session_state.design_purpose = purpose
        st.rerun()

    st.markdown("")

    for cat in CATS:
        cat_parts = [p for p in all_parts if p["category"] == cat]
        if not cat_parts:
            continue
        picked_here = st.session_state.get(f"ms_{cat}", [])
        label = CATEGORY_TR.get(cat, cat)
        with st.expander(
            f"{label}  ({len(picked_here)}/{len(cat_parts)})",
            expanded=bool(picked_here),
        ):
            st.multiselect(
                label,
                [p["part_id"] for p in cat_parts],
                format_func=lambda i: (
                    parts_by_id[i]["name"]
                    + ("" if parts_by_id[i]["seq_status"] == "full" else "  ·dizi eksik")
                ),
                key=f"ms_{cat}",
                label_visibility="collapsed",
            )

    # Seçimi tek yerde topla - başka hiçbir yerde yazılmaz.
    st.session_state.design_parts = [
        pid for cat in CATS for pid in st.session_state.get(f"ms_{cat}", [])
    ]

    st.markdown("")
    st.subheader("Insert")
    st.session_state.design_insert_name = st.text_input(
        "Insert adı", value=st.session_state.design_insert_name
    )
    raw_insert = st.text_area(
        "Insert dizisi (FASTA ya da düz dizi)",
        value=st.session_state.design_insert,
        height=130,
        placeholder="ATGAGC...",
    )
    st.session_state.design_insert = raw_insert
    insert = clean_seq(raw_insert)
    if insert:
        st.caption(f"{len(insert)} nt okundu.")

    st.session_state.design_enzymes = st.multiselect(
        "Klonlamada kullanılacak enzimler",
        ez.search("", commercial_only=True),
        default=st.session_state.design_enzymes,
        help="Seçerseniz insertin içinde kesip kesmediği ve uçların uyumu "
             "denetlenir.",
    )

# ---------------------------------------------------------------------------
# Montaj
# ---------------------------------------------------------------------------

selected = sort_canonical(db.get_parts(st.session_state.design_parts))
vec = assemble(
    st.session_state.design_name or "pVL-01",
    selected,
    insert_seq=insert,
    insert_name=st.session_state.design_insert_name or "insert",
)
findings = validate(
    vec, host=host, insert_seq=insert,
    cloning_enzymes=st.session_state.design_enzymes,
    purpose=st.session_state.design_purpose,
)
counts = summarize(findings)

# ---------------------------------------------------------------------------
# Sağ: harita + denetim
# ---------------------------------------------------------------------------

with right:
    if not selected:
        st.info(
            "Soldan parça seçin ya da hazır bir şablon yükleyin.",
            icon=":material/arrow_back:",
        )
    else:
        m = st.container(horizontal=True)
        m.metric("Uzunluk", f"{vec.length} bp")
        m.metric("GC", f"%{vec.gc}")
        m.metric("Dizisi bilinen", f"%{vec.resolved_fraction}")
        m.metric("Hata", counts["hata"], delta=None)
        m.metric("Uyarı", counts["uyari"], delta=None)

        tab_map, tab_lin, tab_check, tab_lib, tab_seq = st.tabs(
            ["Dairesel harita", "Doğrusal harita", "Denetim",
             "Kütüphane kur", "Dizi ve dışa aktarım"]
        )

        single = {
            k: v for k, v in restriction_sites(vec.sequence).items() if len(v) == 1
        }

        with tab_map:
            render_svg(circular_map(vec.name, vec.length, vec.features, single))
            if single:
                st.caption(
                    "Tek kesen enzimler: " + ", ".join(sorted(single)[:18])
                    + (" …" if len(single) > 18 else "")
                )

        with tab_lin:
            render_svg(linear_map(vec.name, vec.length, vec.features, single))
            st.caption(
                "Kaset düzenini okumak için doğrusal görünüm daha nettir; "
                "oklar yönü gösterir."
            )

        with tab_check:
            if not findings:
                st.success(
                    "Denetimden geçti - bilinen bir sorun bulunamadı.",
                    icon=":material/check_circle:",
                )
            for f in findings:
                with st.container(border=True):
                    st.markdown(
                        f"{SEV_ICON.get(f.severity, '')} **{f.title}**  "
                        f"`{f.rule_id}`"
                    )
                    st.markdown(f.detail)
                    if f.fix:
                        st.markdown(f"**Ne yapmalı:** {f.fix}")
                    if f.literature_query:
                        st.caption(
                            f"Literatür sorgusu: `{f.literature_query}`"
                        )


        with tab_lib:
            st.caption(
                "Gerçek laboratuvarda ifade düzeyi tek tasarımla "
                "ayarlanmaz: promotör ve RBS değiştirilerek bir seri "
                "kurulur, hepsi denenir, işe yarayan seçilir. Sebebi "
                "basit — ifade düzeyini önceden hesaplayamazsınız, ama "
                "bir aralığı tarayabilirsiniz."
            )

            degisken_kategoriler = ["promoter", "rbs", "ori", "marker",
                                    "terminator"]
            lk = st.container(horizontal=True, vertical_alignment="bottom")
            hangi = lk.multiselect(
                "Hangi kategoriler değişsin?", degisken_kategoriler,
                default=["promoter", "rbs"],
                format_func=lambda c: CATEGORY_TR.get(c, c),
            )
            lk.caption(
                "Seçmediğiniz kategoriler soldaki tasarımdan olduğu gibi "
                "alınır."
            )

            secenekler: dict[str, list[str]] = {}
            for kat in hangi:
                kat_parcalar = [p for p in all_parts
                                if p["category"] == kat
                                and p["seq_status"] == "full"]
                # Varsayilan: tasarimdaki secim + ayni kategoriden gucu
                # bilinen birkac kardes. Tek secenekle acilirsa sekme
                # bos gorunur ve kullanici ne yapacagini anlamaz.
                secili = [p["part_id"] for p in kat_parcalar
                          if p["part_id"] in
                          st.session_state.get(f"ms_{kat}", [])]
                gucu_bilinen = [
                    p["part_id"] for p in kat_parcalar
                    if p["part_id"] not in secili
                    and lib.strength_score(p)[0] is not None
                ]
                varsayilan = (secili + gucu_bilinen)[:4]
                secenekler[kat] = st.multiselect(
                    f"{CATEGORY_TR.get(kat, kat)} seçenekleri",
                    [p["part_id"] for p in kat_parcalar],
                    default=varsayilan,
                    format_func=lambda i: parts_by_id[i]["name"][:52],
                    key=f"lib_{kat}",
                )

            sabit = [pid for pid in st.session_state.design_parts
                     if parts_by_id.get(pid, {}).get("category") not in hangi]
            toplam = 1
            for v in secenekler.values():
                toplam *= max(len(v), 1)

            ozet = st.container(horizontal=True)
            ozet.metric("Sabit parça", len(sabit))
            ozet.metric("Birleşim", toplam if hangi else 0)

            if not hangi or toplam <= 1:
                st.info(
                    "En az bir kategoriyi değişken yapın ve o kategoriden "
                    "birden fazla parça seçin.",
                    icon=":material/category:")
            else:
                varyantlar, lib_uyari = lib.build(
                    fixed_part_ids=sabit, variable=secenekler,
                    host=host, purpose=st.session_state.design_purpose,
                    insert_seq=insert,
                    name_prefix=(st.session_state.design_name or "var"),
                )
                for u in lib_uyari:
                    st.warning(u, icon=":material/warning:")

                aralik = lib.span(varyantlar)
                a = st.container(horizontal=True)
                a.metric("Varyant", len(varyantlar))
                a.metric("Hatasız", sum(1 for v in varyantlar if v.ok))
                if aralik["fold"]:
                    a.metric("Taranan aralık", f"{aralik['fold']}×")
                a.metric("Basamak", len(aralik["tiers"]))

                st.dataframe(
                    [{"Varyant": v.name,
                      "Beklenen ifade": v.expression,
                      "Basamak": v.tier,
                      "Uzunluk": v.length,
                      "GC": v.gc,
                      "Hata": v.errors,
                      "Uyarı": v.warnings} for v in varyantlar],
                    width="stretch", hide_index=True,
                    alt="Kombinatoryal kütüphanedeki varyantlar",
                    column_config={
                        "Beklenen ifade": st.column_config.ProgressColumn(
                            "Beklenen ifade", min_value=0.0, max_value=1.0,
                            format="percent"),
                    },
                )

                st.info(
                    "**Beklenen ifade bir ölçüm değildir.** Üç kaynaktan "
                    "türetilmiş sıralamadır: Anderson promotörlerinde "
                    "sigma70 konsensüsüne uzaklık (diziden hesaplanır), "
                    "katalogdaki güç etiketi ve iGEM başlığındaki "
                    "strong/medium/weak nitelemesi. İfadeyle aynı yönde "
                    "gider ama kaç kat fark edeceğini söylemez - serinin "
                    "varlık sebebi de budur.",
                    icon=":material/functions:")

                if varyantlar and varyantlar[0].basis:
                    with st.expander("En yüksek varyantın dayanağı"):
                        for b in varyantlar[0].basis:
                            st.markdown(f"- {b}")

                hatali = [v for v in varyantlar if not v.ok]
                if hatali:
                    st.warning(
                        f"{len(hatali)} varyant denetimden geçemedi. "
                        "Örnek: " + "; ".join(
                            f"{v.name.split(' · ')[0]} — "
                            f"{', '.join(v.notes[:2])}"
                            for v in hatali[:3]),
                        icon=":material/warning:")

                d = st.container(horizontal=True)
                d.download_button(
                    "Kütüphaneyi CSV olarak indir", lib.to_csv(varyantlar),
                    file_name=f"{st.session_state.design_name}-kutuphane.csv",
                    mime="text/csv", icon=":material/download:",
                )
                secili_varyant = st.selectbox(
                    "Bir varyantı tasarımcıya yükle",
                    [v.name for v in varyantlar],
                )
                if st.button("Bu varyantı yükle",
                             icon=":material/download_for_offline:"):
                    hedef = next(v for v in varyantlar
                                 if v.name == secili_varyant)
                    for cat in CATS:
                        st.session_state[f"ms_{cat}"] = [
                            pid for pid in hedef.part_ids
                            if parts_by_id.get(pid, {}).get("category") == cat
                        ]
                    st.rerun()

        with tab_seq:
            c = st.container(horizontal=True)
            c.download_button(
                "GenBank indir", to_genbank(vec),
                file_name=f"{vec.name}.gb", mime="chemical/seq-na-genbank",
                icon=":material/download:",
            )
            c.download_button(
                "FASTA indir", to_fasta(vec),
                file_name=f"{vec.name}.fasta", mime="text/plain",
                icon=":material/download:",
            )
            if st.button("Yapıyı kaydet", icon=":material/save:"):
                db.save_construct(
                    vec.name, st.session_state.design_host,
                    st.session_state.design_parts, insert,
                    st.session_state.design_insert_name,
                )
                if is_shared_deployment():
                    st.warning(
                        f"`{vec.name}` kaydedildi, ancak " + SHARED_NOTICE,
                        icon=":material/cloud:",
                    )
                else:
                    st.success(f"`{vec.name}` kaydedildi.")

            st.dataframe(
                [
                    {
                        "Sıra": i + 1,
                        "Parça": p["name"],
                        "Tür": CATEGORY_TR.get(p["category"], p["category"]),
                        "Uzunluk": p["length"],
                        "Dizi": "tam" if p["seq_status"] == "full" else "eksik",
                    }
                    for i, p in enumerate(vec.part_order)
                ],
                width="stretch",
                hide_index=True,
                alt=f"{vec.name} yapısının parça sırası tablosu",
            )
            st.code(to_fasta(vec), language="text")
