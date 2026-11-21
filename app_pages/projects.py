"""Projeler: kayıtlı tasarımları aç, geçmişini gör, defter olarak ver.

Kaydetmek zaten vardı ama tek yönlüydü: bir tasarımı kaydedebiliyor,
sonra geri açamıyordunuz. Kayıt, geri alınamıyorsa kayıt değildir.

Bu sayfa üç şey yapar:

1. **Geri yükleme.** Kayıtlı bir tasarımı Tasarımcı'ya geri koyar.
   Tasarımcı parça seçimini `ms_<kategori>` widget'larında tuttuğu
   için geri yükleme o anahtarları doldurmak zorunda; `design_parts`
   yazmak yetmez, widget kendi durumunu korur ve onu yok sayar.

2. **Sürüm geçmişi.** Her kaydediş `construct_versions` tablosuna da
   yazılır. Asıl değerli olan sütun `provenance`: o sürümün hangi kod
   sürümü ve hangi veri dosyalarıyla üretildiğini söyler. "Altı ay önce
   aynı diziyi mi almıştım" sorusu ancak bu damga saklanırsa
   cevaplanabilir - ticari bir araçta aranan şey de budur.

3. **Laboratuvar defteri.** Tasarımı, parçalarını, denetim sonucunu ve
   damgasını tek bir HTML dosyasında verir. Yazdırılabilir, e-postayla
   gönderilebilir, deftere eklenebilir.
"""

import streamlit as st

from core import db, notice, report
from core.assembly import assemble, sort_canonical, to_fasta, to_genbank
from core.env import SHARED_NOTICE, is_shared_deployment
from core.parts_seed import CATEGORY_TR, DESIGNER_CATEGORIES
from core.plasmid_map import circular_map
from core.validate import summarize, validate
from core.i18n import t

SEV_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
            "bilgi": ":material/info:"}

parts_by_id = {p["part_id"]: p for p in db.list_parts()}
kayitlar = db.list_constructs()

if is_shared_deployment():
    st.info(SHARED_NOTICE, icon=":material/cloud:")

if not kayitlar:
    st.info(
        t("Henüz kayıtlı tasarım yok. **Tasarımcı** sayfasında bir yapı "
        "kurup \"Yapıyı kaydet\" düğmesine basın; buraya düşecek."),
        icon=":material/folder_open:",
    )
    st.stop()

st.caption(
    t("{0} kayıtlı tasarım. Her kaydediş bir sürüm bırakır ve her sürüm, hangi kod ve veri sürümüyle üretildiğinin damgasını taşır - aynı tasarımı sonradan açtığınızda aynı sonucu alıp almayacağınızı ancak o damga söyler.").format(len(kayitlar))
)

# ---------------------------------------------------------------------------
# Liste
# ---------------------------------------------------------------------------

st.dataframe(
    [{"Ad": c["name"],
      "Konak": c["host_id"] or "—",
      "Parça": len(c["part_ids"]),
      "Insert (nt)": len(c["insert_seq"]),
      "Sürüm": c["version_count"],
      "Güncelleme": c["updated_at"]}
     for c in kayitlar],
    width="stretch", hide_index=True,
    alt="Kayıtlı tasarımların listesi",
)

st.divider()

adlar = [c["name"] for c in kayitlar]
secili_ad = st.selectbox(t("Tasarım"), adlar, key="proje_secim")
kayit = db.get_construct(secili_ad)
if kayit is None:
    st.warning(t("Tasarım bulunamadı; liste yenilendi."),
               icon=":material/warning:")
    st.rerun()

# ---------------------------------------------------------------------------
# Tasarımı kur
# ---------------------------------------------------------------------------

parcalar = [parts_by_id[p] for p in kayit["part_ids"] if p in parts_by_id]
eksik = [p for p in kayit["part_ids"] if p not in parts_by_id]
vec = assemble(kayit["name"], sort_canonical(parcalar),
               insert_seq=kayit["insert_seq"],
               insert_name=kayit["insert_name"]) if parcalar else None

if eksik:
    st.warning(
        t("Bu tasarımdaki {0} parça katalogda yok: `{1}`. Katalog yenilendiyse parça kimliği değişmiş olabilir; yapı o parçalar olmadan kuruldu.").format(len(eksik), '`, `'.join(eksik)),
        icon=":material/warning:",
    )

ust = st.container(horizontal=True)
ust.metric(t("Uzunluk"), f"{vec.length} bp" if vec else "—")
ust.metric(t("GC"), f"%{vec.gc}" if vec else "—")
ust.metric(t("Parça"), len(parcalar))
ust.metric(t("Sürüm"), len(db.list_construct_versions(secili_ad)))

# ---------------------------------------------------------------------------
# Eylemler
# ---------------------------------------------------------------------------

eylem = st.container(horizontal=True, vertical_alignment="bottom")

if eylem.button(t("Tasarımcı'da aç"), icon=":material/construction:",
                type="primary"):
    # Tasarımcı parça seçimini widget anahtarlarında tutuyor; doğrudan
    # oraya yazmak zorundayız (bkz. modül açıklaması).
    for cat in DESIGNER_CATEGORIES:
        st.session_state[f"ms_{cat}"] = [
            pid for pid in kayit["part_ids"]
            if parts_by_id.get(pid, {}).get("category") == cat
        ]
    st.session_state.design_name = kayit["name"]
    st.session_state.design_host = kayit["host_id"] or "BL21DE3"
    st.session_state.design_parts = kayit["part_ids"]
    st.session_state.design_insert = kayit["insert_seq"]
    st.session_state.design_insert_name = kayit["insert_name"]
    st.switch_page("app_pages/designer.py")

if vec:
    eylem.download_button(
        t("GenBank"), to_genbank(vec), file_name=f"{vec.name}.gb",
        mime="text/plain", icon=":material/download:",
    )
    eylem.download_button(
        t("FASTA"), to_fasta(vec), file_name=f"{vec.name}.fasta",
        mime="text/plain", icon=":material/download:",
    )

# ---------------------------------------------------------------------------
# Not, ad, silme
# ---------------------------------------------------------------------------

with st.expander(t("Adı, notu, silme")):
    yeni_not = st.text_area(
        t("Not"), value=kayit["notes"], height=110, key="proje_not",
        help=t("Bu tasarımla ilgili kendi notunuz. Kaydedilince yeni bir "
             "sürüm olarak da saklanır."),
    )
    if st.button(t("Notu kaydet"), icon=":material/save:"):
        db.save_construct(kayit["name"], kayit["host_id"],
                          kayit["part_ids"], kayit["insert_seq"],
                          kayit["insert_name"], yeni_not)
        st.success(t("Not kaydedildi."))
        st.rerun()

    ad_satir = st.container(horizontal=True, vertical_alignment="bottom")
    yeni_ad = ad_satir.text_input(t("Yeni ad"), value=kayit["name"],
                                  key="proje_yeni_ad")
    if ad_satir.button(t("Yeniden adlandır"), icon=":material/edit:"):
        if db.rename_construct(kayit["name"], yeni_ad):
            st.session_state["proje_secim"] = yeni_ad.strip()
            st.success(f"`{kayit['name']}` → `{yeni_ad.strip()}`")
            st.rerun()
        else:
            st.error(
                t("Ad değiştirilemedi: boş, aynı ya da başka bir tasarım "
                "bu adı kullanıyor."), icon=":material/error:")

    st.markdown("")
    sil_onay = st.checkbox(
        t("`{0}` ve bütün sürüm geçmişini sil").format(kayit['name']),
        key="proje_sil_onay",
        help=t("Geri alınamaz. Silmeden önce defteri indirmek "
             "isteyebilirsiniz."),
    )
    if st.button(t("Sil"), icon=":material/delete:", disabled=not sil_onay):
        db.delete_construct(kayit["name"])
        st.session_state.pop("proje_secim", None)
        st.session_state.pop("proje_sil_onay", None)
        st.rerun()

# ---------------------------------------------------------------------------
# Parçalar ve denetim
# ---------------------------------------------------------------------------

alt_parca, alt_denetim, alt_surum, alt_defter = st.tabs(
    [t("Parçalar"), t("Denetim"), t("Sürüm geçmişi"), t("Laboratuvar defteri")])

with alt_parca:
    if parcalar:
        st.dataframe(
            [{"Kategori": CATEGORY_TR.get(p["category"], p["category"]),
              "Parça": p["name"], "Kimlik": p["part_id"],
              "Uzunluk": p["length"], "Kaynak": p.get("accession") or "—"}
             for p in sort_canonical(parcalar)],
            width="stretch", hide_index=True,
            alt="Tasarımdaki parçalar",
        )
    else:
        st.info(t("Bu tasarımda katalogdan parça yok."),
                icon=":material/info:")
    if kayit["insert_seq"]:
        st.caption(f"**Insert:** {kayit['insert_name'] or 'adsız'} — "
                   f"{len(kayit['insert_seq'])} nt")

with alt_denetim:
    if vec:
        konak = db.get_host(kayit["host_id"]) if kayit["host_id"] else None
        bulgular = validate(vec, host=konak,
                            insert_seq=kayit["insert_seq"])
        ozet = summarize(bulgular)
        o = st.container(horizontal=True)
        o.metric(t("Hata"), ozet["hata"])
        o.metric(t("Uyarı"), ozet["uyari"])
        o.metric(t("Bilgi"), ozet["bilgi"])
        for f in bulgular:
            with st.container(border=True):
                st.markdown(f"{SEV_ICON.get(f.severity, '')} "
                            f"**{f.title}**  `{f.rule_id}`")
                st.markdown(f.detail)
                if f.fix:
                    st.markdown(t("**Ne yapmalı:** {0}").format(f.fix))
        if not bulgular:
            st.success(t("Denetim bir şey bulmadı."),
                       icon=":material/check_circle:")
        st.caption(
            t("Hata vermemesi tasarımın çalışacağı anlamına gelmez: "
            "uygulama yalnızca kendisine öğretilen kuralları bilir.")
        )
    else:
        st.info(t("Denetlenecek bir yapı kurulamadı."),
                icon=":material/info:")

with alt_surum:
    surumler = db.list_construct_versions(secili_ad)
    st.caption(
        t("En yeni sürüm başta. `Damga` sütunu o sürümün hangi kod ve veri "
        "sürümüyle üretildiğini söyler; iki sürümün damgası farklıysa "
        "aradaki fark tasarımdan değil, altındaki veriden de gelebilir.")
    )
    st.dataframe(
        [{"Sürüm": s["version_id"], "Tarih": s["saved_at"],
          "Konak": s["host_id"] or "—", "Parça": len(s["part_ids"]),
          "Insert (nt)": len(s["insert_seq"]),
          "Not": s["notes"] or "—", "Damga": s["provenance"] or "—"}
         for s in surumler],
        width="stretch", hide_index=True,
        alt="Tasarımın sürüm geçmişi",
    )
    if len(surumler) >= 2:
        yeni, onceki = surumler[0], surumler[1]
        girdi = set(yeni["part_ids"]) - set(onceki["part_ids"])
        cikan = set(onceki["part_ids"]) - set(yeni["part_ids"])
        satir = []
        if girdi:
            satir.append("eklenen: `" + "`, `".join(sorted(girdi)) + "`")
        if cikan:
            satir.append("çıkarılan: `" + "`, `".join(sorted(cikan)) + "`")
        if len(yeni["insert_seq"]) != len(onceki["insert_seq"]):
            satir.append(f"insert {len(onceki['insert_seq'])} → "
                         f"{len(yeni['insert_seq'])} nt")
        if yeni["host_id"] != onceki["host_id"]:
            satir.append(f"konak {onceki['host_id']} → {yeni['host_id']}")
        st.markdown(
            t("**Son iki sürüm arasındaki fark:** ")
            + ("; ".join(satir) if satir
               else "parça ve insert aynı, yalnızca not ya da damga değişti")
        )

with alt_defter:
    st.caption(
        t("Tasarımı, parçalarını, denetim sonucunu ve damgasını tek bir "
        "HTML dosyasına yazar. Tarayıcıda açılır ve yazdırılabilir.")
    )
    if not vec:
        st.info(t("Defter için kurulabilir bir yapı gerekiyor."),
                icon=":material/info:")
    else:
        yazar = st.text_input(t("Hazırlayan (isteğe bağlı)"),
                              key="proje_yazar")
        konak = db.get_host(kayit["host_id"]) if kayit["host_id"] else None
        bulgular = validate(vec, host=konak, insert_seq=kayit["insert_seq"])
        harita = circular_map(vec.name, vec.length, vec.features)
        html = report.build_html(
            name=vec.name, length=vec.length, gc=vec.gc,
            host_name=(konak or {}).get("name", kayit["host_id"] or ""),
            parts=sort_canonical(parcalar),
            map_svg=harita,
            findings=bulgular,
            insert_name=kayit["insert_name"],
            insert_bp=len(kayit["insert_seq"]),
            author=yazar,
        )
        st.download_button(
            t("Defteri indir (HTML)"), html,
            file_name=f"{vec.name}-defter.html", mime="text/html",
            icon=":material/download:", type="primary",
        )
        st.caption(notice.FOOTER)
