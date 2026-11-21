"""Primer tasarımcısı: vektöre özel ve genel primer araçları.

Primer tasarımı uygulamada zaten vardı ama dağınıktı: bir sekmesi Dizi
analizinde, montaj primerleri Laboratuvar'da, dizileme primerleri yine
orada. Üstelik hepsi primeri yalnızca *hedef dizinin kendisine* göre
ölçüyordu - oysa primer laboratuvarda hedefe değil **bütün plazmite**
girer.

Bu sayfa ikisini birden yapar: primer işlerini tek yere toplar ve
hepsine vektörü gösterir. Hesap `core/primers.py` içindedir.
"""

import streamlit as st

from core import db, primers as pm, ref_vectors as rv
from core.assembly import assemble, sort_canonical
from core.i18n import t
from core.seq_utils import clean_seq, design_primers, enzyme_info, revcomp

KIND_FN = {"hata": st.error, "uyari": st.warning, "bilgi": st.info}
KIND_ICON = {"hata": ":material/error:", "uyari": ":material/warning:",
             "bilgi": ":material/info:"}


def notlari_yaz(notlar):
    for n in notlar:
        KIND_FN.get(n["kind"], st.info)(
            f"**{n['title']}** — {n['detail']}",
            icon=KIND_ICON.get(n["kind"]))


st.caption(t(
    "Primer tasarımı Tm ve GC'den ibaret değildir: primer plazmitte "
    "başka bir yere de bağlanıyorsa PCR ya yanlış ürün verir ya hiç "
    "vermez, ve bunu Tm söylemez. Buradaki araçların hepsi primeri "
    "seçtiğiniz vektörün TAMAMINA karşı denetler."
))

# ---------------------------------------------------------------------------
# Vektör kaynağı - bütün sekmeler bunu kullanır
# ---------------------------------------------------------------------------

kaynak = st.segmented_control(
    t("Vektör nereden gelsin?"),
    ["Tasarımcı'daki yapı", "Referans vektör", "Dizi yapıştır"],
    default="Tasarımcı'daki yapı",
    format_func=t) or "Tasarımcı'daki yapı"

vektor_dizi, vektor_ad, dairesel = "", "", True
insert_araligi = None
parts_by_id = {p["part_id"]: p for p in db.list_parts()}

if kaynak == "Tasarımcı'daki yapı":
    secili = st.session_state.get("design_parts", [])
    parcalar = [parts_by_id[p] for p in secili if p in parts_by_id]
    if not parcalar:
        st.info(
            t("Tasarımcı sayfasında henüz bir yapı yok. Oradan parça "
              "seçin, ya da yukarıdan başka bir kaynak kullanın."),
            icon=":material/science:")
    else:
        vec = assemble(
            st.session_state.get("design_name", "pVL-01"),
            sort_canonical(parcalar),
            insert_seq=clean_seq(st.session_state.get("design_insert", "")),
            insert_name=st.session_state.get("design_insert_name", "insert"),
        )
        vektor_dizi, vektor_ad = vec.sequence, vec.name
        ins = [f for f in vec.features if f.part_id == "INSERT"]
        if ins:
            insert_araligi = (ins[0].start, ins[0].end)

elif kaynak == "Referans vektör":
    vlist = rv.list_vectors()
    vid = st.selectbox(
        t("Referans vektör"), [v["vector_id"] for v in vlist],
        format_func=lambda i: next(
            f"{v['name']} — {v['length']} bp"
            for v in vlist if v["vector_id"] == i),
    )
    kayit = rv.load(vid)
    if kayit:
        vektor_dizi, vektor_ad = kayit["sequence"], kayit["name"]

else:
    c = st.container(horizontal=True, vertical_alignment="bottom")
    vektor_ad = c.text_input(t("Vektör adı"), value="plazmitim", width=240)
    dairesel = c.toggle(t("Dairesel"), value=True)
    ham = st.text_area(t("Dizi (FASTA ya da düz)"), height=130,
                       placeholder=t("ATGC… ya da >adı ile başlayan FASTA"))
    if ham.strip():
        govde = "".join(s for s in ham.splitlines()
                        if not s.startswith(">"))
        vektor_dizi = clean_seq(govde)

if vektor_dizi:
    u = st.container(horizontal=True)
    u.metric(t("Vektör"), vektor_ad)
    u.metric(t("Uzunluk"), f"{len(vektor_dizi)} bp")
    u.metric(t("Topoloji"), t("dairesel") if dairesel else t("doğrusal"))

st.divider()

sek_pcr, sek_koloni, sek_capa, sek_incele = st.tabs(
    [t("PCR primeri"), t("Koloni PCR"), t("Dizileme çapaları"),
     t("Primer incele")])

# ---------------------------------------------------------------------------
# 1. Bölgeyi çoğaltan PCR primeri
# ---------------------------------------------------------------------------

with sek_pcr:
    st.caption(t(
        "Vektörün bir bölgesini çoğaltacak primer çifti. Uçlara "
        "restriksiyon bölgesi eklerseniz enzimin verimli kesmesi için "
        "koruyucu bazlar da eklenir; ayrıca o enzimin bölgenin içinde "
        "kesip kesmediği denetlenir."
    ))
    if not vektor_dizi:
        st.info(t("Önce bir vektör seçin."), icon=":material/science:")
    else:
        c = st.container(horizontal=True, vertical_alignment="bottom")
        bas = c.number_input(t("Başlangıç"), min_value=1,
                            max_value=len(vektor_dizi), value=1, step=1)
        son = c.number_input(t("Bitiş"), min_value=1,
                            max_value=len(vektor_dizi),
                            value=min(600, len(vektor_dizi)), step=1)
        hedef_tm = c.number_input(t("Hedef Tm (°C)"), min_value=50.0,
                                 max_value=72.0, value=60.0, step=0.5)

        ez_listesi = [""] + sorted(
            e for e in ("EcoRI", "BamHI", "HindIII", "XhoI", "SalI",
                        "NdeI", "NcoI", "XbaI", "SpeI", "KpnI", "SacI",
                        "NotI", "PstI", "SphI", "AscI", "PacI"))
        c2 = st.container(horizontal=True, vertical_alignment="bottom")
        f_ez = c2.selectbox(t("5' enzim"), ez_listesi, index=0)
        r_ez = c2.selectbox(t("3' enzim"), ez_listesi, index=0)

        if son <= bas:
            st.warning(t("Bitiş konumu başlangıçtan büyük olmalı."),
                       icon=":material/warning:")
        else:
            hedef = vektor_dizi[int(bas) - 1:int(son)]
            res = design_primers(hedef, f_ez or None, r_ez or None,
                                 target_tm=float(hedef_tm),
                                 name_prefix="PCR")
            if "error" in res:
                st.error(res["error"], icon=":material/error:")
            else:
                st.dataframe(
                    [{t("Primer"): p.name, t("Dizi"): p.sequence,
                      t("Uzunluk"): p.length, "Tm": p.tm,
                      "GC": f"%{p.gc}",
                      t("Kuyruk"): p.tail or "—"}
                     for p in (res["forward"], res["reverse"])],
                    width="stretch", hide_index=True,
                    alt=t("Tasarlanan ileri ve geri primerler"))
                m = st.container(horizontal=True)
                m.metric(t("Ürün uzunluğu"), f"{res['product_length']} bp")
                m.metric(t("Önerilen bağlanma sıcaklığı"),
                         f"{res['annealing_temp']} °C"
                         if res["annealing_temp"] else "—")

                for p in (res["forward"], res["reverse"]):
                    for w in p.warnings:
                        st.warning(f"{p.name}: {w}",
                                   icon=":material/warning:")
                for w in res["pair_warnings"]:
                    st.error(w, icon=":material/error:")

                # Asıl katkı: primerleri bölgeye değil, BÜTÜN vektöre
                # karşı denetle.
                st.markdown(t("**Vektörün tamamına karşı özgüllük**"))
                for p in (res["forward"], res["reverse"]):
                    oz = pm.specificity(p.anneal, vektor_dizi,
                                        circular=dairesel)
                    for n in oz["notes"]:
                        KIND_FN.get(n["kind"], st.info)(
                            f"**{p.name}** — {n['detail']}",
                            icon=KIND_ICON.get(n["kind"]))

# ---------------------------------------------------------------------------
# 2. Koloni PCR
# ---------------------------------------------------------------------------

with sek_koloni:
    st.caption(t(
        "Ligasyondan çıkan kolonilerin bir kısmı kendi üzerine kapanmış "
        "BOŞ vektördür. Koloni PCR bunu ayırmak içindir ve primerler "
        "bilerek insertin değil, iskeletin dizisinden alınır: insertin "
        "içine bağlanan bir primer boş vektörde hiç ürün vermez, yani "
        "boş koloni ile tutmamış PCR birbirinden ayırt edilemez. "
        "İskeletten bakınca küçük bant boş, büyük bant dolu demektir."
    ))
    if not vektor_dizi:
        st.info(t("Önce bir vektör seçin."), icon=":material/science:")
    else:
        if insert_araligi:
            st.success(
                t("Insert Tasarımcı'daki yapıdan bulundu: "
                  f"{insert_araligi[0]}-{insert_araligi[1]}"),
                icon=":material/check_circle:")
        c = st.container(horizontal=True, vertical_alignment="bottom")
        i_bas = c.number_input(
            t("Insert başlangıcı"), min_value=1, max_value=len(vektor_dizi),
            value=insert_araligi[0] if insert_araligi else 1, step=1)
        i_son = c.number_input(
            t("Insert bitişi"), min_value=1, max_value=len(vektor_dizi),
            value=(insert_araligi[1] if insert_araligi
                   else min(500, len(vektor_dizi))), step=1)
        payanda = c.number_input(
            t("İskeletten pay (bp)"), min_value=20, max_value=300,
            value=60, step=10,
            help=t("Primer insertin bu kadar dışına yerleştirilir. "
                   "Büyütmek iki bandı da büyütür ama boş bandı jelde "
                   "daha görünür kılar."))

        if i_son <= i_bas:
            st.warning(t("Insert bitişi başlangıçtan büyük olmalı."),
                       icon=":material/warning:")
        else:
            r = pm.colony_pcr(vektor_dizi, int(i_bas), int(i_son),
                              circular=dairesel, offset=int(payanda))
            st.dataframe(
                [{t("Primer"): p.name, t("Dizi"): p.sequence,
                  t("Uzunluk"): p.length, "Tm": p.tm, "GC": f"%{p.gc}"}
                 for p in (r.forward, r.reverse)],
                width="stretch", hide_index=True,
                alt=t("Koloni PCR primerleri"))
            m = st.container(horizontal=True)
            m.metric(t("Insert giren koloni"), f"{r.with_insert} bp")
            m.metric(t("Boş vektör"), f"{r.without_insert} bp")
            m.metric(t("Bantlar arası fark"), f"{r.fark} bp")
            for p in (r.forward, r.reverse):
                for w in p.warnings:
                    st.warning(w, icon=":material/warning:")
            notlari_yaz(r.notes)

# ---------------------------------------------------------------------------
# 3. Dizileme çapaları
# ---------------------------------------------------------------------------

with sek_capa:
    st.caption(t(
        "Evrensel dizileme primerlerinin dizileri bu depoya ezberden "
        "yazılmadı. Bunun yerine primer vektörün KENDİ dizisinden "
        "kesilir: katalogda dizisi doğrulanmış bir promotör bu vektörde "
        "bulunuyorsa, oradan ileri okuyacak primer zaten vektörün "
        "içindedir. Böylece önerilen her primerin kaynağı bellidir."
    ))
    if not vektor_dizi:
        st.info(t("Önce bir vektör seçin."), icon=":material/science:")
    else:
        capalar = pm.catalog_anchors(vektor_dizi, db.list_parts(),
                                     circular=dairesel)
        if not capalar:
            st.info(
                t("Bu vektörde çapa olarak kullanılabilecek katalog "
                  "parçası bulunamadı. Çapa yalnızca dizisi tam olan "
                  "promotör, operatör ve çoklu klonlama bölgelerinden "
                  "çıkarılır."),
                icon=":material/info:")
        else:
            st.dataframe(
                [{t("Parça"): a["name"], t("Kimlik"): a["part_id"],
                  t("Konum"): f"{a['start']}-{a['end']}",
                  t("Primer"): a["primer"], "Tm": a["tm"],
                  "GC": f"%{a['gc']}"}
                 for a in capalar],
                width="stretch", hide_index=True,
                alt=t("Vektörde bulunan dizileme çapaları"))
            st.caption(t(
                "Sanger okuması primerin hemen ardındaki ilk ~50 bazı "
                "güvenilir okumaz; çapadan sonraki bölgeyi görmek "
                "istiyorsanız Laboratuvar sayfasındaki dizileme planı "
                "kaç primer gerektiğini hesaplar."))

# ---------------------------------------------------------------------------
# 4. Primer incele
# ---------------------------------------------------------------------------

with sek_incele:
    st.caption(t(
        "Elinizdeki bir primeri yapıştırın: Tm, GC, 3' ucu, kendi "
        "dimeri ve saç tokası ölçülür. Yukarıda bir vektör seçiliyse "
        "primerin o vektörde nereye bağlandığı da aranır."
    ))
    ham_primer = st.text_area(
        t("Primer dizisi"), height=90,
        placeholder=t("Bir ya da birden fazla primer, her satıra bir "
                      "tane…"))
    if ham_primer.strip():
        satirlar = [s.strip() for s in ham_primer.splitlines()
                    if s.strip() and not s.startswith(">")]
        for ham in satirlar[:10]:
            p = clean_seq(ham)
            if not p:
                continue
            with st.container(border=True):
                sonuc = pm.analyze(p, vektor_dizi, circular=dairesel)
                if "error" in sonuc:
                    st.error(sonuc["error"], icon=":material/error:")
                    continue
                o = sonuc["metrics"]
                st.code(o["sequence"], language="text", wrap_lines=True)
                m = st.container(horizontal=True)
                m.metric(t("Uzunluk"), f"{o['length']} nt")
                m.metric("Tm", f"{o['tm']} °C" if o["tm"] else "—")
                m.metric("GC", f"%{o['gc']}")
                m.metric(t("3' son 5 baz"), o["gc_clamp"])
                notlari_yaz(sonuc["notes"])
