"""Referanstan kes-yapıştır: gerçek bir vektörü restriksiyonla revize et.

Akış: literatürden gerçek bir vektör al, hangi enzimle nereden kestiğini
göster, çıkan ve eklenen parçaları söyle, orijinal ile yaptığımız
vektörü yan yana koy. Bütün şekiller uygulamanın kendi harita ve jel
çiziminden gelir; dizi hesabı `core/refbuild.py`'dedir.

Tarifler kayıt defterindedir (`refbuild.RECIPES`): her biri farklı bir
referans vektör, farklı bir enzim çifti ve farklı bir terapötik kaset
kullanır, ve her biri ayrı bir şey öğretir.
"""

import streamlit as st

from core import db, disease_seed as ds, gel, ref_vectors, refbuild
from core.assembly import assemble
from core.parts_seed import CATEGORY_TR
from core.plasmid_map import circular_map, render as render_svg
from core.seq_utils import Feature

st.caption(
    f"{len(refbuild.RECIPES)} tarif var. Sıfırdan kurmuyoruz: gerçek bir "
    "klonlama vektörü alınır, kesilip yeni parçalar yapıştırılır. "
    "Etiketler: `Literatür` makaleden, `Öneri` bizim kararımız."
)

tarif_key = st.selectbox(
    "Tarif", [r.key for r in refbuild.RECIPES],
    format_func=lambda k: refbuild.recipe_by_key(k).title,
)
tarif = refbuild.recipe_by_key(tarif_key)
rec = ref_vectors.load(tarif.vector_id)
if not rec:
    st.error(f"{tarif.vector_id} kaydı bulunamadı (data/reference_vectors).")
    st.stop()

hastalik = ds.by_key(tarif.disease) if tarif.disease else None
L0 = rec["length"]

# Kasetleri kur
kasetler = [assemble(f"kaset{i+1}", db.get_parts(s.insert_parts))
            for i, s in enumerate(tarif.steps)]
adimlar = [(s.left, s.right, k.sequence, s.title)
           for s, k in zip(tarif.steps, kasetler)]
try:
    res = refbuild.build_recipe(rec["sequence"], adimlar)
except ValueError as exc:
    st.error(f"Bu tarif kurulamadı: {exc}", icon=":material/error:")
    st.stop()

URUN_ADI = f"p{(hastalik.short if hastalik else 'VL').split()[0]}-" \
           f"{rec['name'].split()[0]}"

# ---------------------------------------------------------------------------
st.header("1. Referans vektör")
c1, c2 = st.columns([3, 2])
ORIJINAL = ref_vectors.catalog_features(tarif.vector_id, db.list_parts())
with c1:
    render_svg(circular_map(rec["name"], L0, ORIJINAL, {}, size=560))
with c2:
    st.markdown(
        f"**{rec['name']} — {L0} bp.** {rec.get('ozet', '')} "
        f"Dizi: NCBI `{rec.get('version', rec['accession'])}`.")
    st.markdown(f"**Neden bu vektör? `Öneri`** {tarif.why_vector}")
    if tarif.source:
        st.caption(f"Kaynak: {tarif.source}")
    if hastalik:
        st.info(f"**Bu tarif ne öğretiyor:** {tarif.teaches}",
                icon=":material/school:")

# ---------------------------------------------------------------------------
st.header("2. Hangi enzimler, neden")
satirlar = []
kaynak = rec["sequence"]
for i, (adim, basamak) in enumerate(zip(tarif.steps, res.steps), start=1):
    for yon, enz in (("sol", adim.left), ("sağ", adim.right)):
        bilgi = refbuild.ENZYMES[enz]
        pos = refbuild.find_sites(kaynak, enz)
        satirlar.append({
            "Enzim": enz, "Tanıma dizisi": bilgi["site"],
            "Bıraktığı uç": bilgi["ovh"],
            "Bu aşamada kesim sayısı": len(pos),
            "Konum (bp)": ", ".join(map(str, pos)) or "—",
            "Aşama": f"{i}. aşama {yon}",
            "Metilasyona duyarlı": "evet" if bilgi.get("dam") else "hayır",
        })
    kaynak = basamak.sequence
st.dataframe(satirlar, width="stretch", hide_index=True,
             alt="Kullanılan restriksiyon enzimleri tablosu")
st.caption(
    "Her aşamada iki FARKLI enzim kullanılır, ama **farklı olmaları "
    "yetmez**: ne bıraktıklarına bakın. SalI ile XhoI ikisi de "
    "5' TCGA, XbaI ile SpeI ikisi de 5' CTAG, BamHI ile BclI ikisi "
    "de 5' GATC bırakır - böyle bir çiftle keserseniz uçlar "
    "birbirine uyar ve insert ters de girer, klonlarınızın yarısı "
    "yanlış olur. 'Bıraktığı uç' sütunu bu yüzden tabloda duruyor; "
    "uygulama her aşamada kendisi de denetler ve uyumlu bir çift "
    "görürse uyarır.\n\n"
    "İkinci kontrol: kullanılan enzim **eklenecek kasetin içinde** "
    "kesmemeli. Kesiyorsa inserti o enzimlerle hazırlayamazsınız ve "
    "tarif kağıt üzerinde kalır. On iki tarifin hepsi bu iki "
    "sınavdan geçiyor (`tests/test_refbuild.py`)."
)
if any(refbuild.ENZYMES[e].get("dam")
       for a in tarif.steps for e in (a.left, a.right)):
    st.warning(
        "Tarifte Dam metilasyonuna duyarlı bir enzim var. E. coli'den "
        "saflaştırılan plazmitte o bölge kesilmeyebilir; plazmiti önce "
        "JM110 gibi dam⁻ bir suşta çoğaltın.",
        icon=":material/warning:")

# ---------------------------------------------------------------------------
st.header("3. Adım adım kes-yapıştır")
sekmeler = st.tabs([f"{i+1}. aşama: {s.title}"
                    for i, s in enumerate(tarif.steps)])
for sekme, adim, basamak, kaset in zip(sekmeler, tarif.steps, res.steps,
                                       kasetler):
    with sekme:
        m = st.container(horizontal=True)
        m.metric("Atılan bölge", f"{basamak.removed} bp")
        m.metric("Eklenen parça", f"{basamak.insert_length} bp")
        m.metric("Ürün", f"{basamak.length} bp")
        m.metric("Enzimler", f"{adim.left} + {adim.right}")
        if adim.lesson:
            st.markdown(f"**Ne oluyor:** {adim.lesson}")
        st.dataframe(
            [{"Sıra": i + 1, "Parça": p["name"][:52],
              "Uzunluk": p.get("length") or 0,
              "Dizi": "tam" if p.get("seq_status") == "full" else "eksik"}
             for i, p in enumerate(kaset.part_order)],
            width="stretch", hide_index=True,
            alt=f"{adim.title} adımında eklenen parçalar")
        for n in basamak.notes[1:]:
            st.warning(n)

# ---------------------------------------------------------------------------
st.header("4. Orijinal ve bizim vektör")


def _uzunluk(a: int, b: int, L: int) -> int:
    return (b - a + 1) if a <= b else (L - a + 1 + b)


def _tasi(start: int, end: int):
    """Bir bölgeyi bütün adımlardan geçir: (baş, son, durum, kalan)."""
    a, b, durum, L = start, end, "tam", L0
    for basamak in res.steps:
        m = refbuild.map_feature(basamak, a, b, L)
        if m is None:
            return None
        a, b, d = m
        L = basamak.length
        if d == "kesildi":
            durum = "kesildi"
    return a, b, durum, _uzunluk(a, b, L)


son_feats: list[Feature] = []
olanlar: list[dict] = []
for f in ORIJINAL:
    tam = _uzunluk(f.start, f.end, L0)
    tasinan = _tasi(f.start, f.end)
    if tasinan:
        a, b, durum, kalan = tasinan
        ad = f.name if durum == "tam" else f"{f.name.split('(')[0].strip()} (bozuldu)"
        son_feats.append(Feature(a, b, f.strand, ad, f.category))
    else:
        durum, kalan = "atıldı", 0
    olanlar.append({
        "Bölge": f.name[:44],
        "Tür": CATEGORY_TR.get(f.category, f.category),
        "Orijinal": tam,
        "Son vektörde": kalan,
        "Kalan": kalan / tam if tam else 0.0,
        "Durum": {"tam": "bütün korundu",
                  "kesildi": "kesildi — işlevi bitti",
                  "atıldı": "tamamen atıldı"}[durum],
    })

# Eklenen kasetlerin öznitelikleri: kendi adımından sonraki adımlardan geçir
for i, (kaset, basamak) in enumerate(zip(kasetler, res.steps)):
    for f in kaset.features:
        a = basamak.insert_start + f.start - 1
        b = basamak.insert_start + f.end - 1
        L = basamak.length
        atik = False
        for sonraki in res.steps[i + 1:]:
            m = refbuild.map_feature(sonraki, a, b, L)
            if m is None:
                atik = True
                break
            a, b, _d = m
            L = sonraki.length
        if not atik:
            son_feats.append(Feature(a, b, f.strand, f.name, f.category,
                                     f.part_id))
son_feats.sort(key=lambda f: f.start)

a, b = st.columns(2)
with a:
    st.subheader("Orijinal")
    render_svg(circular_map(rec["name"], L0, ORIJINAL, {}, size=520))
with b:
    st.subheader("Bizim vektör")
    render_svg(circular_map(URUN_ADI, res.final_length, son_feats, {},
                            size=520))

st.markdown("**Orijinal bölgelere ne oldu?**")
st.dataframe(
    olanlar, width="stretch", hide_index=True,
    alt="Referans vektörün bölgelerinin kurulan vektördeki durumu",
    column_config={
        "Orijinal": st.column_config.NumberColumn("Orijinal", format="%d bp"),
        "Son vektörde": st.column_config.NumberColumn("Son vektörde",
                                                      format="%d bp"),
        "Kalan": st.column_config.ProgressColumn(
            "Kalan oran", min_value=0.0, max_value=1.0, format="percent"),
    })
bozulan = [o["Bölge"] for o in olanlar if o["Durum"] != "bütün korundu"]
st.caption(
    "Sayılar diziden ölçülür, elle yazılmaz: her bölge adım adım "
    "koordinat eşlemesinden geçirilir ve ürün üzerinde ne kadarının "
    "kaldığına bakılır. **Kesilen bir gen bozulmuş gendir**: kalan "
    "parça çerçeve içinde bile olsa tam protein vermez. "
    + (f"Burada {len(bozulan)} bölge bozuldu: " + ", ".join(bozulan) + "."
       if bozulan else "Bu tarifte hiçbir bölge bozulmadı.")
)

st.dataframe(
    [{"": "Boyut", "Orijinal": f"{L0} bp",
      "Bizim vektör": f"{res.final_length} bp "
                      f"({res.final_length - L0:+d})"},
     {"": "Marker", "Orijinal": ", ".join(rec.get("markerlar") or ["—"]),
      "Bizim vektör": ", ".join(
          o["Bölge"].split("(")[0].strip() for o in olanlar
          if o["Tür"] == CATEGORY_TR.get("marker", "marker")
          and o["Durum"] == "bütün korundu") or "işlevli marker kalmadı"},
     {"": "Terapötik içerik", "Orijinal": "yok",
      "Bizim vektör": (hastalik.logic if hastalik else "kaset eklendi")},
     {"": "Eklenen parça sayısı", "Orijinal": "—",
      "Bizim vektör": f"{sum(len(k.part_order) for k in kasetler)} parça, "
                      f"{sum(k.length for k in kasetler)} bp"}],
    width="stretch", hide_index=True,
    alt="Orijinal ve kurulan vektörün karşılaştırması")

# ---------------------------------------------------------------------------
st.header("5. Jel simülasyonu")
ilk = tarif.steps[0]
son = tarif.steps[-1]
boyutlar = refbuild.digest_sizes(res.final, [son.left, son.right])
agaroz = gel.recommend_agarose(
    refbuild.digest_sizes(rec["sequence"], [ilk.left, ilk.right]) + boyutlar)
lanes = [gel.ladder_lane("1 kb"),
         gel.build_lane(f"{rec['name'][:10]}\n{ilk.left}+{ilk.right}",
                        refbuild.digest_sizes(rec["sequence"],
                                              [ilk.left, ilk.right]), agaroz)]
for i, (adim, basamak) in enumerate(zip(tarif.steps, res.steps), start=1):
    lanes.append(gel.build_lane(
        f"{i}. aşama\n{adim.left}+{adim.right}",
        refbuild.digest_sizes(basamak.sequence, [adim.left, adim.right]),
        agaroz))
render_svg(gel.render_gel(lanes, agaroz, title="Beklenen sindirim deseni"))
st.caption(
    "Teorik desen: gerçek jel kesim verimine ve yüklenen DNA miktarına "
    "bağlıdır. Doğru klonu ararken orijinal ile ürünün bant deseni "
    "arasındaki farka bakarsınız.")

# ---------------------------------------------------------------------------
st.header("6. Sonuç ve sınırlılıklar")
eksik = [p["name"] for k in kasetler for p in k.part_order
         if p.get("seq_status") != "full"]
notlar = [
    f"- Son vektör **{res.final_length} bp**; {rec['name']}'e göre "
    f"{res.final_length - L0:+d} bp.",
    "- Deneysel doğrulama yok: bu bir **tasarım önerisidir**.",
]
if hastalik and hastalik.gap:
    notlar.append(f"- Hastalık tarafındaki açık: {hastalik.gap}")
st.markdown("\n".join(notlar))
if eksik:
    st.warning(
        "Dizisi henüz çekilmemiş parçalar: " + "; ".join(sorted(set(eksik)))
        + ". Enzim kontrolü bu parçalar tamamlanınca yeniden yapılmalı "
          "(`python tools/fetch_pku_parts.py`).",
        icon=":material/download:")

c = st.container(horizontal=True)
c.download_button(
    "Son vektör FASTA",
    f">{URUN_ADI} {res.final_length} bp\n" + "\n".join(
        res.final[i:i + 70] for i in range(0, res.final_length, 70)) + "\n",
    file_name=f"{URUN_ADI}.fasta", mime="text/plain",
    icon=":material/download:")
if c.button("Plazmit denetle sayfasında aç", icon=":material/fact_check:"):
    st.session_state.analyze_seq = res.final
    st.switch_page("app_pages/inspect.py")
