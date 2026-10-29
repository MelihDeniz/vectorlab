"""Parça kataloğu: gözat, düzenle, ekle, dizi eksiklerini doldur."""

import streamlit as st

from core import db, igem_seed as ig
from core.parts_seed import CATEGORIES, CATEGORY_TR
from core.seq_utils import clean_seq, gc_percent

SOURCE_TR = {"seed": "VectorLab çekirdek katalog",
             "igem": "iGEM Registry (CC BY 4.0)",
             "vector": "referans vektörlerden çıkarıldı",
             "user": "sizin eklediğiniz"}

counts = db.category_counts()
total = sum(counts.values())

head = st.container(horizontal=True, vertical_alignment="bottom")
with head:
    search = st.text_input(
        "Ara", placeholder="ad, kod ya da not içinde ara…", width=320
    )
    cat = st.selectbox(
        "Kategori",
        ["tümü"] + CATEGORIES,
        format_func=lambda c: (
            "Tümü" if c == "tümü"
            else f"{CATEGORY_TR.get(c, c)} ({counts.get(c, 0)})"
        ),
        width=280,
    )
    kaynak = st.selectbox(
        "Kaynak", ["tümü", "seed", "igem", "user"],
        format_func=lambda k: "Tümü" if k == "tümü" else SOURCE_TR.get(k, k),
        width=280,
    )

parts = db.list_parts(None if cat == "tümü" else cat, search)
if kaynak != "tümü":
    parts = [p for p in parts if p["source"] == kaynak]
_stubs = sum(1 for p in db.list_parts() if p["seq_status"] != "full")
st.caption(
    f"{len(parts)} parça gösteriliyor · katalogda toplam {total} parça. "
    + (
        f"{_stubs} parçanın dizisi eksik; erişim numarasından gerçek diziyi "
        "çekip buradan doldurabilirsiniz."
        if _stubs
        else "Hepsinin gerçek dizisi gömülü; her parçanın kaynağı erişim "
             "numarası ve koordinatıyla yazılıdır."
    )
)

tab_list, tab_igem, tab_new = st.tabs(
    ["Katalog", "iGEM serileri", "Yeni parça ekle"]
)

with tab_list:
    for p in parts:
        full = p["seq_status"] == "full"
        badge = "" if full else "  ·  :orange[dizi eksik]"
        with st.expander(
            f"{p['name']}  —  {CATEGORY_TR.get(p['category'], p['category'])}"
            f"{badge}",
            icon=":material/check_circle:" if full else ":material/pending:",
        ):
            meta = st.container(horizontal=True)
            meta.metric("Uzunluk", f"{p['length'] or 0} bp")
            if full and p["sequence"]:
                meta.metric("GC", f"%{gc_percent(p['sequence'])}")
            if p["strength"]:
                meta.metric("Güç", p["strength"])
            if p["copy_number"]:
                meta.metric("Kopya sayısı", p["copy_number"])

            if p["induction"]:
                st.markdown(f"**İndüksiyon:** {p['induction']}")
            if p["requires"]:
                st.markdown("**Gerektirir:** " + ", ".join(p["requires"]))
            if p["notes"]:
                st.markdown(p["notes"])
            if p["accession"]:
                st.caption(f"Kaynak: {p['accession']}")

            if full:
                st.code(p["sequence"], language="text", wrap_lines=True)
            else:
                st.info(
                    "Bu parçanın dizisi gömülü değil. Uzunluk ve erişim "
                    "numarası kayıtlı; montajda doğru uzunlukta N bloğu "
                    "olarak yer alır. Gerçek diziyi aşağıya yapıştırırsanız "
                    "dizi düzeyinde analiz de yapılabilir hale gelir.",
                    icon=":material/info:",
                )
                with st.form(f"fill_{p['part_id']}"):
                    newseq = st.text_area(
                        "Gerçek dizi (FASTA ya da düz)", height=110,
                        key=f"ta_{p['part_id']}",
                    )
                    if st.form_submit_button("Diziyi kaydet",
                                             icon=":material/save:"):
                        s = clean_seq(newseq)
                        if not s:
                            st.error("Geçerli bir dizi okunamadı.")
                        else:
                            expected = p["length"] or 0
                            rec = dict(p)
                            rec["sequence"] = s
                            rec["seq_status"] = "full"
                            rec["source"] = "user"
                            db.upsert_part(rec)
                            if expected and abs(len(s) - expected) > 0:
                                st.warning(
                                    f"Kaydedildi, ancak uzunluk beklenenden "
                                    f"farklı: katalogda {expected} bp "
                                    f"yazıyordu, girilen {len(s)} bp. "
                                    "Doğru bölgeyi aldığınızı kontrol edin.",
                                    icon=":material/warning:",
                                )
                            else:
                                st.success("Dizi kaydedildi.")
                            st.rerun()

            refs = db.refs_for_part(p["part_id"])
            if refs:
                st.markdown("**Bağlı kaynaklar**")
                for r in refs:
                    st.markdown(
                        f"- [{r['title']}](https://pubmed.ncbi.nlm.nih.gov/"
                        f"{r['pmid']}/) — {r['journal']} {r['year']}"
                    )

            if p["source"] == "user":
                if st.button("Parçayı sil", key=f"del_{p['part_id']}",
                             icon=":material/delete:"):
                    db.delete_part(p["part_id"])
                    st.rerun()

with tab_igem:
    st.caption(
        f"Kaynak: {ig.source_name()}. {ig.license_note()} Her parçanın "
        "kaydında Registry kimliği ve lisansı yazılıdır."
    )

    aile = ig.anderson_family()
    st.subheader("Anderson sabit promotör ailesi")
    st.markdown(
        "Ailenin 20 üyesi 35 bp uzunluğunda ve yalnızca iki noktada "
        f"ayrışır: -35 kutusu ({ig.BOX_35[0]}..{ig.BOX_35[1]}. baz) ve -10 "
        f"kutusu ({ig.BOX_10[0]}..{ig.BOX_10[1]}. baz). Geri kalan 23 "
        "nükleotid hepsinde aynıdır, yani aradaki ifade farkı tamamen bu "
        "iki hexamerden gelir."
    )
    if aile:
        st.dataframe(
            [{
                "Parça": r["igem_id"],
                "-35 kutusu": r["box35"],
                "-35 farkı": r["d35"],
                "-10 kutusu": r["box10"],
                "-10 farkı": r["d10"],
                "Toplam fark": r["total"],
                "Konsensüs benzerliği": (12 - r["total"]) / 12,
            } for r in aile],
            width="stretch", hide_index=True,
            alt="Anderson promotör ailesinin konsensüs karşılaştırması",
            column_config={
                "Konsensüs benzerliği": st.column_config.ProgressColumn(
                    "Konsensüse benzerlik", min_value=0.0, max_value=1.0,
                    format="percent",
                ),
            },
        )
        st.info(
            f"**Bu sütun ölçüm değil, hesaplama.** sigma70 konsensüsü "
            f"(-35 `{ig.SIGMA70_35}`, -10 `{ig.SIGMA70_10}`) ile kaç "
            "nükleotid farkı olduğunu sayar. Konsensüse yakınlık güçle aynı "
            "yönde gider ama birebir değildir: spacer uzunluğu, yukarı akış "
            "elemanları ve ölçüm koşulları da etkiler. Gerçek ölçüm "
            "değerleri Registry'nin kendi parça sayfalarındadır.",
            icon=":material/functions:",
        )
        en_guclu = aile[0]
        st.caption(
            f"Üst sınır {en_guclu['igem_id']}: iki kutusu da konsensüsün "
            "tam kendisi. Devrenizde ifadeyi kademelendirmek istiyorsanız "
            "listenin iki ucundan ve ortasından birer promotör seçip aynı "
            "geni üçünden ifade edin."
        )

    st.subheader("Ribozom bağlanma bölgesi serisi")
    rbs = [p for p in ig.igem_parts() if p["category"] == "rbs"]
    st.dataframe(
        [{"Parça": p["igem_id"], "Uzunluk": p["length"],
          "Dizi": p["sequence"],
          "Registry tanımı": p["name"].split(" - ", 1)[-1]}
         for p in sorted(rbs, key=lambda x: x["igem_id"])],
        width="stretch", hide_index=True,
        alt="iGEM ribozom bağlanma bölgesi serisi",
    )
    st.caption(
        "Bu seride güçler birbirine göre tanımlıdır (Registry tanımındaki "
        "strong / medium / weak etiketleri). Promotörü sabit tutup yalnızca "
        "RBS'i değiştirmek, çeviri düzeyini ayarlamanın en temiz yoludur."
    )

    kalan = {}
    for p in ig.igem_parts():
        if p["category"] not in ("rbs",) and not p["igem_id"].startswith("BBa_J231"):
            kalan.setdefault(p["category"], []).append(p["igem_id"])
    st.subheader("Diğer iGEM parçaları")
    st.dataframe(
        [{"Kategori": CATEGORY_TR.get(c, c), "Sayı": len(v),
          "Parçalar": ", ".join(sorted(v))}
         for c, v in sorted(kalan.items())],
        width="stretch", hide_index=True,
        alt="Kategorilere göre diğer iGEM parçaları",
    )


with tab_new:
    st.caption(
        "Kendi parçanızı ekleyin. Dizi girmezseniz parça 'eksik' olarak "
        "kaydedilir ve sadece uzunluğuyla montaja girer."
    )
    with st.form("new_part"):
        c1, c2 = st.columns(2)
        pid = c1.text_input("Kod (benzersiz)", placeholder="PRO_KENDI_01")
        name = c2.text_input("Ad", placeholder="Kendi promotörüm")
        c3, c4 = st.columns(2)
        category = c3.selectbox(
            "Kategori", CATEGORIES,
            format_func=lambda c: CATEGORY_TR.get(c, c),
        )
        strength = c4.text_input("Güç / özellik", placeholder="güçlü")
        sequence = st.text_area("Dizi", height=120)
        c5, c6 = st.columns(2)
        accession = c5.text_input("Erişim numarası")
        pmid = c6.text_input("PMID")
        length_manual = st.number_input(
            "Dizi girmiyorsanız uzunluk (bp)", 0, 20000, 0, step=10
        )
        notes = st.text_area("Notlar", height=90)

        if st.form_submit_button("Kataloğa ekle", icon=":material/add:"):
            s = clean_seq(sequence)
            if not pid.strip() or not name.strip():
                st.error("Kod ve ad zorunlu.")
            elif db.get_part(pid.strip()):
                st.error(f"`{pid}` kodu zaten kullanılıyor.")
            elif not s and not length_manual:
                st.error("Ya bir dizi girin ya da uzunluk belirtin.")
            else:
                db.upsert_part({
                    "part_id": pid.strip(), "name": name.strip(),
                    "category": category, "sequence": s,
                    "seq_status": "full" if s else "stub",
                    "length": len(s) if s else int(length_manual),
                    "strength": strength or None, "notes": notes,
                    "accession": accession, "pmid": pmid, "source": "user",
                })
                st.success(f"`{name}` kataloğa eklendi.")
                st.rerun()

st.divider()
saved = db.list_constructs()
if saved:
    st.subheader("Kayıtlı yapılar")
    st.dataframe(
        [{"Ad": c["name"], "Konak": c["host_id"],
          "Parça sayısı": len(c["part_ids"]),
          "Insert (nt)": len(c["insert_seq"]),
          "Güncelleme": c["updated_at"]}
         for c in saved],
        width="stretch", hide_index=True,
        alt="Kaydedilmiş vektör yapıları",
    )
