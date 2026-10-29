"""Tasarımdan laboratuvar protokolü üretir - sayılarla.

Buradaki her değer tasarımın kendisinden hesaplanır: hangi antibiyotik
kaç µg/mL, ligasyona kaç ng insert, doğrulama jeli yüzde kaç, koloni
PCR ürünü kaç bp. Amaç, tasarımı "güzel bir harita" olmaktan çıkarıp
tezgâhta uygulanabilir bir plana çevirmek.

Derişimler standart E. coli çalışma derişimleridir; laboratuvarınızın
kendi protokolü farklıysa onunkini kullanın.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core import gel

# --- Seçilim antibiyotikleri: parça kimliği -> (ad, çalışma derişimi, not)
ANTIBIOTICS: dict[str, tuple[str, str, str]] = {
    "MRK_BLA": ("Ampisilin", "100 µg/mL",
                "Beta-laktamaz besiyerine salınır ve antibiyotiği çevrede "
                "bozar; uzun inkübasyonda plak üzerinde satellit koloniler "
                "çıkar. Karbenisilin (100 µg/mL) daha kararlıdır."),
    "MRK_KAN": ("Kanamisin", "50 µg/mL",
                "Satellit koloni yapmaz, seçilim temizdir."),
    "MRK_CAT": ("Kloramfenikol", "25 µg/mL",
                "Bakteriyostatiktir, öldürmez - koloniler ampisiline göre "
                "daha yavaş çıkar, bir gece yetmezse 24 saate kadar bekleyin."),
    "MRK_TET": ("Tetrasiklin", "12,5 µg/mL",
                "Işığa duyarlıdır; stok çözeltiyi karanlıkta saklayın ve "
                "plakları folyoya sarın."),
    "MRK_AADA": ("Spektinomisin", "50 µg/mL", ""),
}

# --- İndüksiyon: promotör parçası -> (indükleyici, derişim, nasıl)
INDUCTION: dict[str, tuple[str, str, str]] = {
    "PRO_T7LAC": ("IPTG", "0,4 mM",
                  "OD600 0,5-0,7'ye gelince ekleyin. Protein çökelti "
                  "(inclusion body) oluşturuyorsa 0,1 mM'e düşürüp "
                  "sıcaklığı 18 °C'ye çekin ve gece boyu indükleyin."),
    "PRO_T7": ("IPTG", "0,4 mM",
               "T7 promotöründe operatör yok: indüksiyon öncesi kaçak "
               "ifade belirgindir. Ürün toksikse pLysS taşıyan bir suş "
               "kullanın."),
    "PRO_TAC": ("IPTG", "0,5 mM", "Konağın kendi RNA polimerazı sürer; "
                "T7 kadar yüksek ifade vermez ama daha dengelidir."),
    "PRO_TRC": ("IPTG", "0,5 mM", "tac ile aynı mantık."),
    "PRO_LAC": ("IPTG", "0,1-1 mM",
                "Besiyerinde glukoz varsa katabolit baskılanması yüzünden "
                "indüksiyon zayıf kalır; glukozsuz besiyeri kullanın."),
    "PRO_ARABAD": ("L-arabinoz", "%0,02 (w/v)",
                   "Asıl üstünlüğü doz yanıtlı olmasıdır: %0,00002 ile "
                   "%0,2 arasında bir seri kurup ifade düzeyini "
                   "kademelendirebilirsiniz. Konak arabinozu tüketmemeli "
                   "(TOP10 gibi araBAD mutantı seçin)."),
    "PRO_RHABAD": ("L-ramnoz", "2 mM", "P-BAD'e benzer, doz yanıtlıdır."),
    "PRO_TET": ("Anhidrotetrasiklin (aTc)", "100 ng/mL",
                "Antibiyotik etkisi pratik olarak yoktur, bu yüzden "
                "devrelerde tetrasikline tercih edilir."),
    "PRO_PL": ("Sıcaklık", "30 °C → 42 °C",
               "cI857 duyarlı baskılayıcısı 42 °C'de bozulur. Kimyasal "
               "indükleyici gerekmez, ama 42 °C ısı şoku yanıtını da "
               "tetikler ve proteazları artırır."),
}

# Sabit promotörler indüksiyon istemez.
CONSTITUTIVE_PREFIX = ("PRO_J23", "IGEM_J231")


@dataclass
class Reagent:
    name: str
    amount: str
    note: str = ""


@dataclass
class Step:
    no: int
    title: str
    detail: str
    reagents: list[Reagent] = field(default_factory=list)
    conditions: str = ""
    warning: str = ""


@dataclass
class Protocol:
    title: str
    summary: str
    steps: list[Step] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Hesaplar
# ---------------------------------------------------------------------------


def ligation_amounts(
    vector_bp: int, insert_bp: int, vector_ng: float = 50.0,
    ratio: float = 3.0,
) -> dict[str, Any]:
    """Ligasyon için insert kütlesi.

    Molar oran kütleye şöyle çevrilir:
        ng insert = ng vektör × (insert bp / vektör bp) × oran

    Yani 3:1 molar oran, insert vektörden kısaysa kütlece ÇOK daha az
    insert demektir. Bu hesabı atlayıp eşit kütle koymak, ligasyonun en
    sık sessiz hatasıdır.
    """
    if vector_bp <= 0 or insert_bp <= 0:
        return {"error": "Vektör ve insert uzunluğu gerekli."}
    insert_ng = vector_ng * (insert_bp / vector_bp) * ratio
    return {
        "vector_ng": round(vector_ng, 1),
        "insert_ng": round(insert_ng, 1),
        "ratio": ratio,
        "vector_bp": vector_bp,
        "insert_bp": insert_bp,
        "formula": (f"{vector_ng:g} ng × ({insert_bp}/{vector_bp}) × "
                    f"{ratio:g} = {insert_ng:.1f} ng"),
    }


def digest_setup(
    dna_ng: float, enzyme_names: list[str], volume_ul: float = 50.0,
) -> list[Reagent]:
    """Sindirim reaksiyonunun kurulumu."""
    n = max(len(enzyme_names), 1)
    enzim_ul = 1.0
    tampon_ul = volume_ul / 10.0
    su = volume_ul - tampon_ul - n * enzim_ul - 5.0
    out = [
        Reagent("Plazmit DNA", f"{dna_ng:g} ng (~5 µL)",
                "Miniprep DNA'sı tipik olarak 100-300 ng/µL gelir."),
        Reagent("10× reaksiyon tamponu", f"{tampon_ul:g} µL",
                "Çift kesimde iki enzimin de iyi çalıştığı tamponu seçin; "
                "üreticinin uyumluluk tablosuna bakın."),
    ]
    for e in enzyme_names:
        out.append(Reagent(e, f"{enzim_ul:g} µL (~10 U)",
                           "Enzim hacmi toplam hacmin %10'unu geçmesin - "
                           "gliserol yıldız aktivitesi yapar."))
    out.append(Reagent("Nükleazsız su", f"{max(su, 0):g} µL'ye tamamla", ""))
    out.append(Reagent("Toplam", f"{volume_ul:g} µL", ""))
    return out


def colony_pcr_product(
    insert_bp: int, flank_bp: int = 200,
) -> dict[str, Any]:
    """Vektör primerleriyle koloni PCR'ında beklenen ürün boyutları."""
    bos = flank_bp
    dolu = insert_bp + flank_bp
    return {
        "empty_bp": bos,
        "insert_bp": dolu,
        "difference": dolu - bos,
        "resolvable": (dolu - bos) > 100,
    }


# ---------------------------------------------------------------------------
# Protokol kurma
# ---------------------------------------------------------------------------


def _antibiotic(part_ids: list[str]) -> tuple[str, str, str] | None:
    for pid in part_ids:
        if pid in ANTIBIOTICS:
            return ANTIBIOTICS[pid]
    return None


def _induction(part_ids: list[str]) -> tuple[str, str, str] | None:
    for pid in part_ids:
        if pid in INDUCTION:
            return INDUCTION[pid]
    for pid in part_ids:
        if any(pid.startswith(p) for p in CONSTITUTIVE_PREFIX):
            return ("yok", "-",
                    "Sabit (konstitütif) promotör: indükleyici "
                    "gerekmez, hücre promotörü sürekli okur. Buna karşılık "
                    "ifadeyi kapatamazsınız; ürün toksikse plazmit "
                    "kaybedilir ya da promotörde mutasyon seçilir.")
    return None


def build(
    vector_name: str,
    vector_bp: int,
    part_ids: list[str],
    enzymes: list[str],
    insert_bp: int = 0,
    insert_name: str = "insert",
    host_name: str = "",
    purpose: str = "ifade",
    backbone_bp: int = 0,
    fragment_sizes: list[int] | None = None,
) -> Protocol:
    """Tasarıma özgü klonlama + ifade protokolü üret."""
    ab = _antibiotic(part_ids)
    ind = _induction(part_ids)
    iskelet = backbone_bp or vector_bp
    p = Protocol(
        title=f"{vector_name} — klonlama ve ifade planı",
        summary=(
            f"{vector_bp} bp hedef yapı"
            + (f", {insert_bp} bp insert" if insert_bp else "")
            + (f", konak {host_name}" if host_name else "")
            + (f", kesim {' + '.join(enzymes)}" if enzymes else "")
        ),
    )
    no = 0

    def add(title, detail, reagents=None, conditions="", warning=""):
        nonlocal no
        no += 1
        p.steps.append(Step(no, title, detail, reagents or [], conditions,
                            warning))

    # --- 1. Vektör ve insert hazırlığı
    if enzymes:
        add(
            "Vektörü kes",
            f"{vector_name} iskeletini "
            f"{' ve '.join(enzymes)} ile kes. Kesim bittikten sonra "
            f"reaksiyonu ısıyla inaktive edin ya da kolondan geçirin; "
            f"kalan enzim ligasyonda yeni kurulan bölgeyi tekrar açar.",
            digest_setup(1000, enzymes),
            "37 °C, 1 saat",
            ("Tek enzimle kesiyorsanız kesilmiş vektörü alkalen fosfatazla "
             "(rSAP/CIP) muamele edin; yoksa vektör kendi üzerine kapanır "
             "ve plak boş kolonilerle dolar."
             if len(enzymes) < 2 else ""),
        )
        if insert_bp:
            add(
                f"{insert_name} fragmanını hazırla",
                f"Inserti aynı enzim çiftiyle kesin. PCR ürününden "
                f"çalışıyorsanız primerlerin 5' ucundaki koruyucu bazların "
                f"yerinde olduğundan emin olun: enzimler DNA'nın tam "
                f"ucundaki tanıma dizisini verimli kesmez.",
                digest_setup(1000, enzymes),
                "37 °C, 1 saat",
            )
        add(
            "Jelden kes ve saflaştır",
            "Kesilmiş vektörü ve inserti jelde yürütün, doğru bantları "
            "kesip kolonla saflaştırın. Jelden saflaştırma, kesilmemiş "
            "vektörü ve kesilip atılan küçük parçayı ortamdan çıkarır - "
            "boş koloni oranını en çok düşüren adım budur.",
            conditions=(
                f"%{gel.recommend_agarose(fragment_sizes or [iskelet, insert_bp])} "
                f"agaroz, 100 V, ~40 dakika"),
            warning="UV masasında uzun tutmayın; UV DNA'da kırık yapar ve "
                    "ligasyon verimi düşer. Mavi ışık masası varsa onu "
                    "tercih edin.",
        )

    # --- 2. Ligasyon
    if insert_bp:
        lig = ligation_amounts(iskelet, insert_bp)
        add(
            "Ligasyon kur",
            "T4 DNA ligazla vektör ve inserti birleştirin. Kütleleri "
            "molar orandan hesaplayın - eşit kütle koymak yanlıştır, "
            f"çünkü insert vektörden {'kısa' if insert_bp < iskelet else 'uzun'}.",
            [
                Reagent("Kesilmiş vektör", f"{lig['vector_ng']} ng",
                        f"{iskelet} bp iskelet"),
                Reagent(f"Kesilmiş {insert_name}", f"{lig['insert_ng']} ng",
                        f"3:1 molar oran → {lig['formula']}"),
                Reagent("10× T4 ligaz tamponu", "2 µL",
                        "ATP içerir ve donma-çözülmeye duyarlıdır; "
                        "küçük hacimlere bölerek saklayın."),
                Reagent("T4 DNA ligaz", "1 µL", ""),
                Reagent("Su", "20 µL'ye tamamla", ""),
            ],
            "16 °C'de gece boyu, ya da oda sıcaklığında 1 saat",
            "Yanına mutlaka bir kontrol kurun: insert YERİNE su koyun. O "
            "tüpten çıkan koloni sayısı, asıl tüpteki boş koloni "
            "arka planını gösterir.",
        )

    # --- 3. Transformasyon
    if ab:
        ad, derisim, abnot = ab
        add(
            "Transformasyon",
            f"Ligasyon karışımının 2-5 µL'sini kimyasal yetkin "
            f"{host_name or 'E. coli'} hücrelerine ekleyin. Buzda 30 "
            f"dakika bekletin, 42 °C'de 45 saniye ısı şoku verin, hemen "
            f"buza alın, 500 µL antibiyotiksiz besiyeri ekleyip 37 °C'de "
            f"1 saat çalkalayın (toparlanma), sonra plakaya yayın.",
            [Reagent(f"{ad} plağı", derisim, abnot)],
            "37 °C, gece boyu",
            "Toparlanma adımını atlamayın: direnç proteini daha "
            "üretilmemişken antibiyotiğe konan hücre ölür. Ampisilinde bu "
            "adım kısa tutulabilir ama kanamisinde zorunludur.",
        )

    # --- 4. Doğrulama
    if insert_bp:
        pcr = colony_pcr_product(insert_bp)
        add(
            "Kolonileri tara",
            f"Tek kolonileri kürdanla alıp koloni PCR kurun. Vektör "
            f"üzerindeki primerlerle insert taşıyan koloni "
            f"~{pcr['insert_bp']} bp, boş vektör ~{pcr['empty_bp']} bp "
            f"ürün verir."
            + ("" if pcr["resolvable"] else
               " İki ürün boyutu birbirine çok yakın - jelde ayırt etmek "
               "zor olabilir, doğrudan sindirimle doğrulamayı tercih edin."),
            conditions=f"%{gel.recommend_agarose([pcr['insert_bp'], pcr['empty_bp']])} agaroz",
            warning="Kürdanla aldığınız koloniyi önce yedek bir plağa "
                    "çizin; PCR'ı tutan koloniyi sonra bulamazsanız baştan "
                    "başlarsınız.",
        )
    if enzymes and fragment_sizes:
        pct = gel.recommend_agarose(fragment_sizes)
        add(
            "Sindirimle doğrula",
            "Pozitif kolonilerden miniprep yapın ve tanı amaçlı sindirim "
            "kurun. Beklenen bant boyutları: "
            + ", ".join(f"{s} bp" for s in sorted(fragment_sizes,
                                                  reverse=True)[:6])
            + ". Jel sekmesinde bu bantların nasıl görüneceğini "
              "çizdirebilirsiniz.",
            conditions=f"%{pct} agaroz",
        )
    add(
        "Diziyle doğrula",
        "Doğru görünen klonu Sanger dizilemeye gönderin. Insertin "
        "tamamını ve iki yandaki birleşme noktalarını okutun: PCR "
        "sırasında girmiş tek bir baz değişimi, haritada ve jelde "
        "görünmez ama proteini bozar. Okuma uzunluğu ~800-1000 bp "
        "olduğu için uzun insertlerde iç primer de gerekir.",
    )

    # --- 5. İfade
    if purpose == "ifade" and ind:
        indad, indder, indnot = ind
        if indad == "yok":
            add("İfade", indnot,
                conditions="37 °C, gece boyu")
        else:
            add(
                "İfadeyi indükle",
                f"Doğrulanmış klonu ifade suşuna aktarın, taze besiyerine "
                f"1:100 seyreltip 37 °C'de büyütün. OD600 0,5-0,7'ye "
                f"geldiğinde {indad} ekleyin.",
                [Reagent(indad, indder, indnot)],
                "İndüksiyon: 37 °C'de 3-4 saat, ya da 18 °C'de gece boyu",
                "İndüksiyondan hemen önce ve hasattan önce birer örnek "
                "alın. İkisini SDS-PAGE'de yan yana yüklemek, proteinin "
                "gerçekten üretildiğini gösteren en basit kanıttır.",
            )
            add(
                "Çözünürlüğü kontrol et",
                "Hücreleri parçaladıktan sonra santrifüjleyin ve "
                "ÜSTTEKİ sıvıyla DİPTEKİ çökeltiyi ayrı ayrı jele yükleyin. "
                "Protein yalnızca çökeltide görünüyorsa inclusion body "
                "oluşmuş demektir: indüksiyon sıcaklığını düşürün, "
                "indükleyiciyi azaltın ya da daha zayıf bir promotör/RBS "
                "seçin.",
            )

    # --- Genel uyarılar
    if not enzymes:
        p.warnings.append(
            "Kesim enzimi seçilmedi, bu yüzden protokolde kesim, ligasyon "
            "ve doğrulama sindirimi adımları yok. Tasarımcı ya da Referans "
            "vektörler sayfasından enzim seçerseniz protokol tamamlanır."
        )
    if not ab:
        p.warnings.append(
            "Yapıda tanınan bir seçilim markerı yok. Marker olmadan "
            "transformasyon sonrası plazmit alan hücreleri ayıramazsınız - "
            "protokoldeki transformasyon adımı bu yüzden eksik."
        )
    if purpose == "ifade" and not ind:
        p.warnings.append(
            "Promotör tanınmadığı için indüksiyon adımı üretilemedi. "
            "Katalogdan bir promotör seçerseniz indükleyici ve derişimi "
            "buraya eklenir."
        )
    return p


def to_markdown(p: Protocol) -> str:
    """Protokolü yazdırılabilir Markdown olarak ver."""
    out = [f"# {p.title}", "", p.summary, ""]
    if p.warnings:
        out.append("> **Eksikler**")
        for w in p.warnings:
            out.append(f"> - {w}")
        out.append("")
    for s in p.steps:
        out.append(f"## {s.no}. {s.title}")
        out.append("")
        out.append(s.detail)
        out.append("")
        if s.reagents:
            out.append("| Bileşen | Miktar | Not |")
            out.append("| --- | --- | --- |")
            for r in s.reagents:
                out.append(f"| {r.name} | {r.amount} | {r.note} |")
            out.append("")
        if s.conditions:
            out.append(f"**Koşullar:** {s.conditions}")
            out.append("")
        if s.warning:
            out.append(f"**Dikkat:** {s.warning}")
            out.append("")
    return "\n".join(out)
