"""Bakteriyel vektör tasarımı doğrulama motoru.

Her kural gerçek bir laboratuvar hatasına karşılık gelir. Çıktı üç
seviyeye ayrılır:

    hata   -> yapı çalışmaz ya da deney kurulamaz
    uyari  -> çalışır ama beklenenden kötü sonuç verir
    bilgi  -> dikkat edilmesi iyi olan nokta

`literature_query` alanı, literatür modülünün o bulgu için kaynak
aramasını sağlar.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core.assembly import AssembledVector
from core.parts_seed import CATEGORY_TR
from core.seq_utils import (
    codon_analysis,
    gc_percent,
    rbs_spacing,
    restriction_sites,
    METHYLATION_BLOCKED,
)

SEVERITY_ORDER = {"hata": 0, "uyari": 1, "bilgi": 2}


@dataclass
class Finding:
    severity: str
    rule_id: str
    title: str
    detail: str
    fix: str = ""
    literature_query: str = ""
    refs: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Promotör - konak uyumluluğu
# ---------------------------------------------------------------------------

#: Promotör part_id -> (gerekli konak özelliği, insan okur açıklama)
_PROMOTER_HOST_NEEDS = {
    "PRO_T7": ("t7_polymerase", "T7 RNA polimeraz"),
    "PRO_T7LAC": ("t7_polymerase", "T7 RNA polimeraz"),
}

#: Promotör part_id -> baskılayıcı gereksinimi
_PROMOTER_REPRESSOR = {
    "PRO_T7LAC": ("REG_LACI", "laciq", "LacI"),
    "PRO_LAC": ("REG_LACI", "laciq", "LacI"),
    "PRO_TAC": ("REG_LACI", "laciq", "LacI"),
    "PRO_TRC": ("REG_LACI", "laciq", "LacI"),
    "PRO_ARABAD": ("REG_ARAC", None, "AraC"),
    "PRO_TET": ("REG_TETR", None, "TetR"),
}

_VERY_STRONG = {"PRO_T7", "PRO_T7LAC", "PRO_PL"}

#: Konak suşlarının taşıdığı yardımcı plazmitler -> (marker, ori inc grubu)
_HOST_RESIDENT_PLASMIDS = {
    "BL21DE3PLYSS": ("kloramfenikol", "IncP15A", "pLysS"),
    "ROSETTA": ("kloramfenikol", "IncP15A", "pRARE"),
}

_MARKER_ANTIBIOTIC = {
    "MRK_BLA": "ampisilin",
    "MRK_KAN": "kanamisin",
    "MRK_CAT": "kloramfenikol",
    "MRK_TET": "tetrasiklin",
    "MRK_AADA": "spektinomisin",
}


def validate(
    vec: AssembledVector,
    host: dict[str, Any] | None = None,
    insert_seq: str = "",
    cloning_enzymes: list[str] | None = None,
    purpose: str = "ifade",
) -> list[Finding]:
    """Yapıyı denetle ve bulguları önem sırasına göre döndür.

    purpose: 'ifade' (protein üretimi) ya da 'klonlama' (sadece çoğaltma)
    """
    f: list[Finding] = []
    pid = {p["part_id"] for p in vec.part_order}
    cloning_enzymes = cloning_enzymes or []

    _check_host_organism(f, vec, host)
    _check_core_elements(f, vec, purpose)
    _check_promoter_host(f, vec, host, pid)
    _check_repressors(f, vec, host, pid)
    _check_origin(f, vec, host, pid)
    _check_markers(f, vec, host, pid)
    _check_host_choice(f, vec, host, pid, purpose)
    _check_reporter(f, vec, host, pid)
    _check_frame_and_tags(f, vec, insert_seq, pid)
    _check_insert_sequence(f, vec, insert_seq, host, cloning_enzymes)
    _check_physical(f, vec)
    _check_enzymes(f, vec, cloning_enzymes)

    return sorted(f, key=lambda x: (SEVERITY_ORDER.get(x.severity, 3), x.rule_id))


# ---------------------------------------------------------------------------
# 1. Temel elemanlar
# ---------------------------------------------------------------------------


def _check_core_elements(f, vec, purpose):
    if not vec.has_category("ori"):
        f.append(Finding(
            "hata", "CORE-01", "Replikasyon orijini yok",
            "Yapıda replikasyon orijini bulunmuyor. Plazmit bakteri "
            "içinde çoğalamaz; transformasyondan sonra hücre bölündükçe "
            "seyrelir ve kaybolur.",
            "Katalogdan bir orijin ekleyin: yüksek kopya için pUC, "
            "düşük kopya/kararlılık için pBR322 ya da p15A.",
            "plasmid origin of replication copy number Escherichia coli",
        ))
    if not vec.has_category("marker"):
        f.append(Finding(
            "hata", "CORE-02", "Seçilim markerı yok",
            "Antibiyotik direnç geni yok. Transformasyon sonrası plazmit "
            "alan hücreleri almayanlardan ayırmanın bir yolu olmaz.",
            "bla (AmpR), kanR (KanR) ya da cat (CmR) ekleyin.",
            "antibiotic selection marker plasmid maintenance E coli",
        ))
    # Bileşik parça (iGEM "device") kendi promotörünü, RBS'ini ve
    # terminatörünü içinde taşır. Bu durumda CORE-03/04/05 kurallarını
    # uygulamak yanlış alarm üretir: eleman eksik değil, tek parçanın
    # içinde geliyor.
    if purpose == "ifade" and vec.has_category("device"):
        adlar = ", ".join(p["name"].split(" - ")[0]
                          for p in vec.parts_of("device"))
        f.append(Finding(
            "bilgi", "DEV-01", "Bileşik parça kullanılıyor",
            f"Yapıda bileşik parça var ({adlar}). Bu parçalar promotör, "
            "ribozom bağlanma bölgesi, kodlayan dizi ve terminatörü "
            "birlikte taşır, bu yüzden promotör/RBS/terminatör eksikliği "
            "denetimi uygulanmadı.",
            "Ayrıca promotör ya da RBS eklediyseniz iki kaset üst üste "
            "biner - bileşik parçayı tek başına kullanın.",
            "BioBrick composite part device expression cassette",
        ))
        fazla = [c for c in ("promoter", "rbs", "terminator")
                 if vec.has_category(c)]
        if fazla:
            f.append(Finding(
                "uyari", "DEV-02", "Bileşik parçanın üstüne eleman eklenmiş",
                "Bileşik parça kendi " + ", ".join(
                    CATEGORY_TR.get(c, c).lower() for c in fazla)
                + " elemanını zaten içeriyor; ayrıca seçtikleriniz montajda "
                "onun önüne yerleşir ve ikinci bir yarım kaset oluşturur.",
                "Bileşik parçayı kullanacaksanız bu kategorilerdeki ayrı "
                "seçimleri kaldırın; kendi kasetinizi kuracaksanız bileşik "
                "parçayı çıkarın.",
                "composite part read-through tandem promoter interference",
            ))
    elif purpose == "ifade":
        if not vec.has_category("promoter"):
            f.append(Finding(
                "hata", "CORE-03", "Promotör yok",
                "İfade kaseti promotör içermiyor - transkripsiyon "
                "başlatılamaz, protein üretilmez.",
                "Konak T7 polimeraz taşıyorsa T7lac, taşımıyorsa tac/trc "
                "ya da sabit ifade için J23100 ekleyin.",
                "bacterial promoter strength protein expression",
            ))
        if not vec.has_category("rbs") and not vec.has_category("mcs"):
            f.append(Finding(
                "hata", "CORE-04", "Ribozom bağlanma bölgesi yok",
                "RBS olmadan mRNA üretilir ama ribozom bağlanamaz: "
                "transkripsiyon olur, çeviri olmaz.",
                "Promotör ile başlama kodonu arasına bir RBS ekleyin "
                "(B0034 güçlü, B0032 orta).",
                "Shine-Dalgarno ribosome binding site translation initiation",
            ))
        if not vec.has_category("terminator"):
            f.append(Finding(
                "uyari", "CORE-05", "Terminatör yok",
                "Transkripsiyon ifade kasetinin sonunda durmaz. "
                "Aşağı akıştaki bölgeler (orijin, marker) istenmeden "
                "okunur; bu hem plazmit kararlılığını hem ifade "
                "verimini düşürür.",
                "Kasetin sonuna T7 terminatörü (T7 promotörü için) ya da "
                "B0015 çift terminatörü (konak polimerazı için) ekleyin.",
                "transcription terminator plasmid stability read-through",
            ))


# ---------------------------------------------------------------------------
# 2. Promotör - konak
# ---------------------------------------------------------------------------


def _check_promoter_host(f, vec, host, pid):
    promoters = vec.parts_of("promoter")
    if len(promoters) > 1:
        names = ", ".join(p["name"] for p in promoters)
        f.append(Finding(
            "uyari", "PRO-01", "Birden fazla promotör",
            f"Yapıda {len(promoters)} promotör var ({names}). Aynı yönde "
            "sıralı promotörler birbirinin transkripsiyonunu engeller "
            "(promotör girişimi) ve ifade öngörülemez hale gelir.",
            "Her kaset için tek promotör bırakın; ikinci bir gen "
            "sürülecekse araya terminatör koyun ya da ayrı kaset kurun.",
            "promoter interference tandem promoters plasmid",
        ))

    if host is None:
        f.append(Finding(
            "bilgi", "PRO-02", "Konak suşu seçilmedi",
            "Konak seçilmeden promotör uyumluluğu, baskılayıcı varlığı ve "
            "seçilim çatışmaları denetlenemez.",
            "Tasarımcı sayfasından bir konak suşu seçin.",
        ))
        return

    for p in promoters:
        need = _PROMOTER_HOST_NEEDS.get(p["part_id"])
        if not need:
            continue
        attr, label = need
        if not host.get(attr):
            f.append(Finding(
                "hata", "PRO-03", f"{p['name']} bu konakta çalışmaz",
                f"{p['name']} transkripsiyon için {label} gerektirir, "
                f"ancak {host['name']} bunu taşımıyor. Yapı doğru kurulsa "
                "bile hiç protein üretilmez - bu en sık yapılan "
                "tasarım hatasıdır.",
                "Konağı BL21(DE3), Rosetta(DE3) ya da başka bir DE3 "
                "lizojenine çevirin. Alternatif olarak konağın kendi "
                "polimerazını kullanan tac/trc promotörüne geçin.",
                "T7 RNA polymerase DE3 lysogen expression host requirement",
                refs=["2199796"],
            ))

    # Yüksek kopya + çok güçlü promotör
    oris = vec.parts_of("ori")
    high_copy = any(o["part_id"] == "ORI_PUC" for o in oris)
    very_strong = [p for p in promoters if p["part_id"] in _VERY_STRONG]
    if high_copy and very_strong:
        f.append(Finding(
            "uyari", "PRO-04", "Yüksek kopya orijin + çok güçlü promotör",
            f"pUC orijini (500-700 kopya) ile {very_strong[0]['name']} "
            "birlikte kullanılıyor. Toplam ifade yükü hücrenin ribozom ve "
            "aminoasit kapasitesini aşar: büyüme durur, plazmit kaybı ve "
            "inklüzyon cismi oluşumu artar.",
            "Orijini pBR322 ya da p15A'ya düşürün; ya da indüksiyonu "
            "düşük IPTG (0.1 mM) ve düşük sıcaklıkta (18-25 C) yapın.",
            "plasmid copy number metabolic burden protein overexpression",
        ))


# ---------------------------------------------------------------------------
# 3. Baskılayıcılar
# ---------------------------------------------------------------------------


def _check_repressors(f, vec, host, pid):
    for p in vec.parts_of("promoter"):
        need = _PROMOTER_REPRESSOR.get(p["part_id"])
        if not need:
            continue
        reg_part, host_attr, label = need
        on_plasmid = reg_part in pid
        on_host = bool(host and host_attr and host.get(host_attr))
        if on_plasmid or on_host:
            continue

        if label == "LacI":
            # E. coli K-12 ve B türevlerinin kromozomunda yabanıl lacI vardır;
            # sorun baskılayıcının YOKLUĞU değil, YETERSİZLİĞİdir. Yüksek
            # kopyalı plazmitteki onlarca operatör kopyası kromozomun tek
            # kopyadan ürettiği LacI'yi tüketir (titrasyon).
            high_copy = any(o["part_id"] == "ORI_PUC" for o in vec.parts_of("ori"))
            detail = (
                f"{p['name']} baskılanmak için LacI'ye ihtiyaç duyar. "
                "Konağın kromozomunda yabanıl *lacI* bulunur, ancak tek "
                "kopyadan üretilen LacI miktarı sınırlıdır"
            )
            detail += (
                "; yüksek kopyalı orijin (pUC) kullandığınız için plazmit "
                "üzerindeki yüzlerce operatör kopyası baskılayıcıyı tüketir "
                "(titrasyon) ve promotör indüksiyon olmadan da açık kalır."
                if high_copy else
                " ve indüksiyon öncesi kaçak ifadeye yol açabilir."
            )
            f.append(Finding(
                "uyari", "REG-01", "LacI baskılaması yetersiz kalabilir",
                detail,
                "Plazmite *lacI* (tercihen lacI-q) ekleyin ya da lacI-q "
                "allelli bir konak seçin. Mavi/beyaz tarama yapıyorsanız "
                "kaçak ifade zaten istenen durumdur, bu uyarıyı "
                "yok sayabilirsiniz.",
                "lacI titration high copy plasmid leaky expression lacIq",
            ))
            continue

        f.append(Finding(
            "hata", "REG-01", f"{label} yok ({p['name']} için gerekli)",
            f"{p['name']} düzenlenebilmek için {label} proteinine ihtiyaç "
            f"duyar; ne plazmitte ne de konakta bulunuyor."
            + (" AraC hem baskılayıcı hem aktivatör olduğu için yokluğunda "
               "P-BAD indüksiyona hiç yanıt vermez."
               if label == "AraC" else
               " Baskılayıcı olmadan promotör sürekli açık kalır: toksik "
               "ürünler hücreyi öldürür ve indüksiyon deneyi anlamsızlaşır."),
            f"{label} genini plazmite ekleyin ya da bunu taşıyan bir "
            "konak seçin.",
            f"{label} repressor leaky expression inducible promoter E coli",
        ))


# ---------------------------------------------------------------------------
# 4. Orijin
# ---------------------------------------------------------------------------


def _check_origin(f, vec, host, pid):
    oris = vec.parts_of("ori")
    replicative = [o for o in oris if o["part_id"] != "ORI_F1"]

    if len(replicative) > 1:
        names = ", ".join(o["name"] for o in replicative)
        f.append(Finding(
            "hata", "ORI-01", "Aynı plazmitte birden fazla orijin",
            f"{names} birlikte bulunuyor. İki replikasyon orijini aynı "
            "molekülde çatışır: replikasyon kontrolü bozulur, plazmit "
            "yeniden düzenlenir ya da kaybolur.",
            "Tek bir orijin bırakın.",
            "plasmid two origins replication incompatibility",
        ))

    from core.db import get_part
    for o in replicative:
        rec = get_part(o["part_id"]) or {}
        if o["part_id"] == "ORI_R6K":
            if not host or "pir" not in (host.get("genotype", "") + host.get("notes", "")).lower():
                f.append(Finding(
                    "hata", "ORI-02", "R6K orijini pir+ konak gerektirir",
                    "R6K gamma orijini replikasyon için pi proteinine "
                    "(pir geni) ihtiyaç duyar. Seçilen konak bunu "
                    "taşımıyor - plazmit hiç çoğalmaz.",
                    "Çoğaltmayı EC100D pir+ gibi bir suşta yapın; ya da "
                    "amacınız suş inşası değilse normal bir orijin seçin.",
                    "R6K gamma origin pir protein conditional replication",
                ))
        if rec.get("requires"):
            for req in rec["requires"]:
                if "RepA" in req:
                    f.append(Finding(
                        "bilgi", "ORI-03", f"{o['name']} RepA gerektirir",
                        f"{o['name']} replikasyon için RepA proteinine "
                        "ihtiyaç duyar. Genelde orijin bölgesiyle birlikte "
                        "klonlanır; ayrı ayrı alındıysa eksik kalabilir.",
                        "Orijin kasetinin repA genini içerdiğini doğrulayın.",
                    ))

    # Konakta yerleşik plazmit ile orijin uyumsuzluğu
    if host and host["host_id"] in _HOST_RESIDENT_PLASMIDS:
        _, inc_group, plasmid_name = _HOST_RESIDENT_PLASMIDS[host["host_id"]]
        for o in replicative:
            rec = get_part(o["part_id"]) or {}
            if rec.get("inc_group") == inc_group:
                f.append(Finding(
                    "hata", "ORI-04",
                    f"Orijin, konağın {plasmid_name} plazmitiyle uyumsuz",
                    f"{host['name']} içinde {plasmid_name} yardımcı "
                    f"plazmiti bulunur ve o da {inc_group} grubundandır. "
                    "Aynı uyumsuzluk grubundan iki plazmit aynı hücrede "
                    "kararlı tutulamaz; biri kaybolur.",
                    "İfade plazmitinizin orijinini ColE1/pUC ailesine "
                    "çevirin ya da yardımcı plazmit taşımayan bir konak "
                    "(düz BL21(DE3)) seçin.",
                    "plasmid incompatibility group coexistence pACYC ColE1",
                ))


# ---------------------------------------------------------------------------
# 5. Markerlar
# ---------------------------------------------------------------------------


def _check_markers(f, vec, host, pid):
    markers = vec.parts_of("marker")
    antibiotics = [_MARKER_ANTIBIOTIC.get(m["part_id"]) for m in markers]
    antibiotics = [a for a in antibiotics if a]

    if len(set(antibiotics)) < len(antibiotics):
        f.append(Finding(
            "uyari", "MRK-01", "Aynı antibiyotik için iki marker",
            "Yapıda aynı antibiyotiğe direnç veren birden fazla gen var. "
            "Gereksiz yük; ayrıca tekrarlı diziler rekombinasyona açıktır.",
            "Fazla markerı çıkarın.",
        ))

    if "MRK_BLA" in pid:
        f.append(Finding(
            "bilgi", "MRK-02", "Ampisilin seçiliminin bilinen sınırı",
            "Beta-laktamaz periplazmaya salgılanır ve besiyerindeki "
            "ampisilini parçalar. Uzun inkübasyonda plazmitsiz "
            "satellit koloniler çıkar; yoğun kültürde (OD>1) seçilim "
            "baskısı kaybolur ve plazmit kaybı başlar.",
            "Plakları 16 saatten fazla bekletmeyin. Daha kararlı olan "
            "karbenisilin (50 ug/mL) kullanın. Büyük ölçekli ifade için "
            "kanamisin markerlı vektöre geçin.",
            "ampicillin satellite colonies beta-lactamase plasmid loss",
        ))

    if "MRK_TET" in pid and "PRO_TET" in pid:
        f.append(Finding(
            "hata", "MRK-03", "tetA markerı ile P-tet promotörü çatışıyor",
            "Aynı plazmitte hem tetrasiklin direnci (tetA) hem de "
            "tetrasiklin ile indüklenen promotör (P-tet) var. Seçilim "
            "için verilen tetrasiklin aynı zamanda promotörü indükler; "
            "indüksiyonu kontrol etmek imkansız hale gelir.",
            "Markerı kanamisine ya da kloramfenikole çevirin.",
            "tet regulatory system anhydrotetracycline induction selection",
        ))

    if host and host["host_id"] in _HOST_RESIDENT_PLASMIDS:
        resident_ab, _, plasmid_name = _HOST_RESIDENT_PLASMIDS[host["host_id"]]
        if resident_ab in antibiotics:
            f.append(Finding(
                "hata", "MRK-04",
                f"Marker, konağın {plasmid_name} plazmitiyle aynı",
                f"{host['name']} suşu {plasmid_name} plazmitini "
                f"{resident_ab} direnci ile taşır. İfade plazmitiniz de "
                "aynı direnci kullanıyorsa iki plazmiti birbirinden "
                "ayırt edemezsiniz ve yardımcı plazmitin kaybını "
                "fark etmezsiniz.",
                "İfade plazmitinizin markerını kanamisin ya da ampisiline "
                "çevirin.",
                "dual plasmid selection chloramphenicol pLysS Rosetta",
            ))

    if "MRK_CCDB" in pid:
        host_ok = bool(host and "ccdB" in host.get("notes", ""))
        if not host_ok:
            f.append(Finding(
                "hata", "MRK-05", "ccdB taşıyan vektör normal konakta çoğaltılamaz",
                "ccdB bir DNA giraz zehiridir. Standart suşlarda ifade "
                "edilince hücreyi öldürür, dolayısıyla vektörün kendisi "
                "çoğaltılamaz.",
                "Vektörü çoğaltmak için DB3.1 gibi gyrA mutantı dirençli "
                "suş kullanın. Insert klonladıktan sonra (ccdB yerine "
                "insert girdiğinde) normal suşa geçebilirsiniz.",
                "ccdB negative selection gyrase poison DB3.1",
            ))


# ---------------------------------------------------------------------------
# 6. Konak seçimi
# ---------------------------------------------------------------------------


def _check_host_choice(f, vec, host, pid, purpose):
    if not host:
        return
    is_expression_design = vec.has_category("promoter") and purpose == "ifade"

    if purpose == "klonlama" and host["use"] == "ifade" and host.get("recA"):
        f.append(Finding(
            "uyari", "HST-01", "Çoğaltma için recA+ ifade suşu seçildi",
            f"{host['name']} recA+ ve endA+ bir ifade suşudur. "
            "Rekombinasyon aktif olduğu için tekrarlı dizi içeren "
            "plazmitler zamanla yeniden düzenlenir; endonükleaz ise "
            "plazmit saflaştırma verimini düşürür.",
            "Klonlama ve plazmit üretimini DH5-alfa ya da TOP10'da yapın, "
            "sonra ifade suşuna aktarın.",
            "recA endA cloning strain plasmid stability",
        ))

    if is_expression_design and host["use"].startswith("klonlama"):
        needs_t7 = any(
            p["part_id"] in _PROMOTER_HOST_NEEDS for p in vec.parts_of("promoter")
        )
        if not needs_t7:
            f.append(Finding(
                "bilgi", "HST-02", "İfade için klonlama suşu seçildi",
                f"{host['name']} bir klonlama suşudur. Promotörünüz "
                "konağın kendi polimerazını kullandığı için bir miktar "
                "ifade olur, ancak proteaz eksik (ompT/lon) bir ifade "
                "suşunda verim belirgin yüksektir.",
                "Protein üretimi için BL21(DE3) ailesine geçin.",
            ))


# ---------------------------------------------------------------------------
# 7. Reporter
# ---------------------------------------------------------------------------


def _check_reporter(f, vec, host, pid):
    if "RPT_LACZA" in pid:
        if host and not host.get("lacz_dm15"):
            f.append(Finding(
                "hata", "RPT-01", "Mavi/beyaz tarama bu konakta çalışmaz",
                "lacZ-alfa parçası, aktif enzim oluşturmak için konaktaki "
                f"lacZ-delta-M15 allelini gerektirir. {host['name']} bunu "
                "taşımıyor: tüm koloniler beyaz kalır ve insert taşıyanları "
                "ayırt edemezsiniz.",
                "Taramayı DH5-alfa, JM109 ya da TOP10 gibi lacZ-dM15 "
                "suşlarında yapın. Besiyerine X-gal ve IPTG eklemeyi "
                "unutmayın.",
                "alpha complementation lacZ delta M15 blue white screening",
            ))
        else:
            f.append(Finding(
                "bilgi", "RPT-02", "Mavi/beyaz tarama için plak koşulları",
                "lacZ-alfa taraması için plaklara X-gal (40 ug/mL) ve "
                "IPTG (0.1 mM) eklenmelidir; aksi halde mavi renk "
                "oluşmaz ve tarama yapılamaz.",
                "",
            ))

    if "RPT_LACZ" in pid and "RPT_LACZA" in pid:
        f.append(Finding(
            "uyari", "RPT-03", "lacZ ve lacZ-alfa birlikte",
            "Tam lacZ varsa alfa-tamamlama taraması anlamsızdır: "
            "koloniler insert alsa da mavi kalır.",
            "İkisinden birini çıkarın.",
        ))


# ---------------------------------------------------------------------------
# 8. Okuma çerçevesi ve etiketler
# ---------------------------------------------------------------------------


def _check_frame_and_tags(f, vec, insert_seq, pid):
    tags = vec.parts_of("tag")
    order = vec.part_order
    idx = {p["part_id"]: i for i, p in enumerate(order)}

    if insert_seq and len(insert_seq) % 3:
        f.append(Finding(
            "hata", "FRM-01", "Insert uzunluğu 3'ün katı değil",
            f"Insert {len(insert_seq)} nt. Etiketle füzyon kurulduğunda "
            "okuma çerçevesi kayar: aşağı akıştaki her şey yanlış "
            "aminoasitlere çevrilir ve genelde erken stop kodonuna girer.",
            f"Insert uzunluğunu 3'ün katına getirin "
            f"({3 - len(insert_seq) % 3} nt ekleyin ya da "
            f"{len(insert_seq) % 3} nt çıkarın).",
            "reading frame fusion tag cloning frameshift",
        ))

    if insert_seq and tags:
        c_tags = [t for t in tags if idx.get(t["part_id"], 0) > idx.get("INSERT", 0)]
        stops = ("TAA", "TAG", "TGA")
        last_codon = insert_seq[-3:].upper() if len(insert_seq) >= 3 else ""
        if c_tags and last_codon in stops:
            f.append(Finding(
                "hata", "FRM-02", "Insert stop kodonu C-terminal etiketi kapatıyor",
                f"Insert {last_codon} stop kodonu ile bitiyor, ancak "
                f"ondan sonra C-terminal etiket ({c_tags[0]['name']}) "
                "geliyor. Çeviri stop kodonunda durur: etiket hiç "
                "üretilmez ve proteini Ni kolonunda yakalayamazsınız. "
                "Bu, pET-28a ile en sık yapılan hatadır.",
                "Insertin stop kodonunu çıkarın; çeviri vektördeki "
                "etiketten sonraki stop kodonunda sonlanacaktır.",
                "C-terminal His tag stop codon pET28 cloning error",
            ))
        if not c_tags and not any(
            idx.get(t["part_id"], 99) < idx.get("INSERT", 99) for t in tags
        ):
            pass
        n_tags = [t for t in tags if idx.get(t["part_id"], 99) < idx.get("INSERT", 99)]
        if n_tags and insert_seq[:3].upper() == "ATG":
            f.append(Finding(
                "bilgi", "FRM-03", "N-terminal etiketten sonra ikinci ATG",
                f"Insert kendi ATG'si ile başlıyor ve önünde N-terminal "
                f"etiket ({n_tags[0]['name']}) var. Çerçeve doğruysa bu "
                "sorun değil (fazladan bir metiyonin kalır), ancak "
                "etiket ile insert arasında kodon sayısının 3'ün katı "
                "olduğunu doğrulayın.",
                "Füzyon bölgesini çevirip çerçeveyi gözle kontrol edin.",
            ))

    # pET-28a MCS'i kendi C-terminal His etiketini İÇİNDE taşır; bu yüzden
    # ayrı bir 'tag' parçası olarak görünmez ve yukarıdaki kontrole yakalanmaz.
    if "MCS_PET28" in pid and insert_seq:
        stops = ("TAA", "TAG", "TGA")
        if insert_seq[-3:].upper() in stops:
            f.append(Finding(
                "uyari", "FRM-05",
                "pET-28a C-terminal His etiketi okunmayacak",
                f"Insert {insert_seq[-3:].upper()} stop kodonu ile bitiyor. "
                "pET-28a MCS'inin sonunda XhoI bölgesinden sonra bir "
                "C-terminal 6xHis etiketi bulunur; çeviri sizin stop "
                "kodonunuzda durduğu için bu etiket üretilmez. "
                "N-terminal His etiketini kullanıyorsanız sorun değil, "
                "ancak NcoI'dan girdiyseniz onu da kaybetmiş olursunuz - "
                "o durumda proteinde hiç etiket kalmaz.",
                "C-terminal etiketi istiyorsanız insertin stop kodonunu "
                "çıkarın. N-terminal etiketi kullanıyorsanız NdeI (NcoI "
                "değil) bölgesinden girdiğinizi doğrulayın.",
                "pET28a N-terminal C-terminal His tag frame NcoI NdeI XhoI",
            ))

    big_tags = [t for t in tags if t["part_id"] in ("TAG_GST", "TAG_MBP", "TAG_TRXA")]
    if big_tags and not vec.has_category("protease_site"):
        f.append(Finding(
            "uyari", "FRM-04", "Büyük füzyon etiketi için kesim bölgesi yok",
            f"{big_tags[0]['name']} büyük bir füzyon ortağıdır "
            "(26-42 kDa). Proteaz kesim bölgesi olmadan etiketi "
            "ayıramazsınız; GST ayrıca dimerleşerek proteinin oligomerik "
            "durumunu yanıltır ve kristalografi/aktivite ölçümlerini bozar.",
            "Etiket ile hedef arasına TEV ya da PreScission bölgesi "
            "ekleyin.",
            "fusion tag removal protease cleavage site GST dimerization",
        ))


# ---------------------------------------------------------------------------
# 9. Insert dizi analizi
# ---------------------------------------------------------------------------


def _check_insert_sequence(f, vec, insert_seq, host, cloning_enzymes):
    if not insert_seq:
        return

    codons = codon_analysis(insert_seq)
    if codons["total"] and codons["rare_pct"] > 8:
        rosetta = bool(host and host.get("rare_trna"))
        sev = "bilgi" if rosetta else "uyari"
        extra = (
            " Seçtiğiniz konak (Rosetta) bu kodonların tRNA'larını "
            "sağladığı için sorun beklenmez."
            if rosetta else
            " Ribozom bu kodonlarda durakladığı için çeviri yarıda "
            "kesilir, kısalmış ürünler ve düşük verim görülür."
        )
        f.append(Finding(
            sev, "INS-01", f"Nadir kodon oranı yüksek (%{codons['rare_pct']})",
            f"Insert {codons['total']} kodonun "
            f"{len(codons['rare'])} tanesi E. coli'de nadir. "
            f"{len(codons['tandem'])} yerde nadir kodonlar yan yana."
            f"{extra}",
            "Rosetta(DE3) konağına geçin ya da diziyi E. coli için kodon "
            "optimize edin.",
            "rare codon usage translational pausing Rosetta tRNA E coli",
        ))

    gc = gc_percent(insert_seq)
    if gc > 65:
        f.append(Finding(
            "uyari", "INS-02", f"Insert GC oranı yüksek (%{gc})",
            "Yüksek GC, PCR sırasında güçlü ikincil yapı oluşturur ve "
            "mRNA düzeyinde RBS erişimini engelleyebilir.",
            "PCR'de DMSO (%3-5) ya da betain kullanın; ifade zayıfsa "
            "5' bölgeyi kodon optimize edin.",
            "high GC content template PCR secondary structure expression",
        ))
    elif gc < 30 and len(insert_seq) > 100:
        f.append(Finding(
            "bilgi", "INS-03", f"Insert GC oranı düşük (%{gc})",
            "Düşük GC dizilerde primer bağlanma sıcaklığı düşer ve "
            "spesifik olmayan çoğaltma riski artar.",
            "Primerleri uzatarak Tm'i 58-62 C aralığına çekin.",
        ))

    spacing = rbs_spacing(vec.sequence)
    bad = [s for s in spacing if not s["ok"]]
    if bad and vec.has_category("rbs"):
        s = bad[0]
        f.append(Finding(
            "uyari", "INS-04", "RBS ile başlama kodonu arası boşluk uygun değil",
            f"{s['atg_pos']}. konumdaki ATG için Shine-Dalgarno motifi "
            f"({s['motif']}) ile arasındaki boşluk {s['gap']} nt - "
            f"{s['verdict']}. Optimum aralık 5-13 nt'dir.",
            "RBS ile ATG arasına uygun uzunlukta spacer ekleyin ya da "
            "boşluğu kısaltın.",
            "Shine-Dalgarno spacing translation initiation efficiency",
        ))

    for enz in cloning_enzymes:
        hits = restriction_sites(insert_seq, [enz], circular=False).get(enz, [])
        if hits:
            f.append(Finding(
                "hata", "INS-05", f"{enz} insertin içinde kesiyor",
                f"{enz} insert dizisinde {len(hits)} yerde kesim yapıyor "
                f"(konum: {', '.join(map(str, hits[:5]))}). Bu enzimle "
                "klonlama yaparsanız insert parçalanır.",
                "Başka bir enzim çifti seçin; ya da bölgeyi sessiz "
                "mutasyonla kaldırın (kodon değiştirip aminoasidi "
                "korumak); ya da Gibson/Golden Gate gibi dizi bağımsız "
                "bir montaj yöntemine geçin.",
                "internal restriction site cloning strategy silent mutation",
            ))


# ---------------------------------------------------------------------------
# 10. Fiziksel özellikler
# ---------------------------------------------------------------------------


def _check_physical(f, vec):
    if vec.stub_regions:
        total = sum(s["length"] for s in vec.stub_regions)
        names = ", ".join(s["name"] for s in vec.stub_regions[:4])
        f.append(Finding(
            "bilgi", "PHY-01", f"{len(vec.stub_regions)} parçanın dizisi eksik",
            f"{names} parçaları katalogda uzunluk ve erişim numarasıyla "
            f"kayıtlı, ancak dizi gömülü değil (toplam {total} bp, "
            f"yapının %{100 - vec.resolved_fraction}'i). Harita, uzunluk "
            "ve kopya sayısı hesapları doğru; ancak bu bölgelerde enzim "
            "kesimi ve ORF analizi yapılamaz.",
            "Gerçek dizileri NCBI'dan erişim numarasıyla çekip "
            "`Parça kataloğu` sayfasından doldurun; ya da elinizdeki "
            "vektörün GenBank dosyasını içe aktarın.",
        ))

    if vec.length > 15000:
        f.append(Finding(
            "uyari", "PHY-02", f"Plazmit büyük ({vec.length} bp)",
            "10-15 kb üzerindeki plazmitlerde kimyasal transformasyon "
            "verimi belirgin düşer ve plazmit kararlılığı azalır.",
            "Elektroporasyon kullanın; gereksiz bölgeleri çıkarın.",
            "plasmid size transformation efficiency electroporation",
        ))

    if vec.length and vec.length < 1500:
        f.append(Finding(
            "bilgi", "PHY-03", f"Plazmit çok küçük ({vec.length} bp)",
            "Tam bir plazmit genelde en az ~2 kb olur (orijin ~600 bp + "
            "marker ~800 bp + kaset). Eksik parça olabilir.",
            "Orijin ve marker parçalarının eklendiğini kontrol edin.",
        ))


# ---------------------------------------------------------------------------
# 11. Enzimler
# ---------------------------------------------------------------------------


def _check_enzymes(f, vec, cloning_enzymes):
    if len(cloning_enzymes) == 1:
        f.append(Finding(
            "uyari", "ENZ-01", "Tek enzimle klonlama yönsüzdür",
            f"Sadece {cloning_enzymes[0]} seçildi. Her iki uçta aynı "
            "uç olduğu için insert iki yönde de girer; kolonilerin "
            "yarısı ters yönde olur ve vektörün kendi üzerine kapanması "
            "(self-ligasyon) arka planı yükseltir.",
            "İki farklı enzim kullanarak yönlü klonlama yapın; "
            "mecbursanız vektörü alkalin fosfataz ile defosforile edin.",
            "directional cloning single enzyme self-ligation background",
        ))

    for enz in cloning_enzymes:
        if enz in METHYLATION_BLOCKED:
            f.append(Finding(
                "uyari", "ENZ-02", f"{enz} metilasyondan etkilenebilir",
                f"{METHYLATION_BLOCKED[enz]}. E. coli'den saflaştırılan "
                "plazmit Dam/Dcm metilasyonu taşıdığı için bu enzim "
                "beklenen bölgeyi kesmeyebilir - sindirim boş çıkar ve "
                "nedeni uzun süre anlaşılmaz.",
                "Plazmiti dam-/dcm- bir suşta (örn. JM110, SCS110) "
                "çoğaltın; ya da metilasyona duyarsız bir izoşizomer "
                "seçin.",
                "Dam Dcm methylation restriction enzyme blocked E coli",
            ))

    if len(cloning_enzymes) >= 2:
        from core.seq_utils import compatible_ends
        for i, a in enumerate(cloning_enzymes):
            for b in cloning_enzymes[i + 1 :]:
                if compatible_ends(a, b):
                    f.append(Finding(
                        "uyari", "ENZ-03",
                        f"{a} ve {b} uyumlu uç bırakıyor",
                        f"{a} ile {b} birbirine ligate edilebilen uçlar "
                        "üretir. Yönlü klonlama amaçlandıysa bu kaybolur: "
                        "insert ters de girebilir, vektör kendi üzerine "
                        "kapanabilir.",
                        "Uyumsuz uç veren bir enzim çifti seçin.",
                        "compatible cohesive ends restriction enzymes ligation",
                    ))


# ---------------------------------------------------------------------------
# Özet
# ---------------------------------------------------------------------------


def summarize(findings: list[Finding]) -> dict[str, int]:
    out = {"hata": 0, "uyari": 0, "bilgi": 0}
    for f in findings:
        if f.severity in out:
            out[f.severity] += 1
    return out


# ---------------------------------------------------------------------------
# Konak organizması uyumu
# ---------------------------------------------------------------------------
# Katalogdaki orijinlerin ve promotörlerin hangi organizmada çalıştığı.
# Bu, E. coli dışı konak seçildiğinde en sık yapılan hatayı yakalar:
# E. coli plazmiti Bacillus'ta çoğalmaz, ama harita gayet güzel görünür.

ORI_ORGANISMS: dict[str, set[str]] = {
    "ORI_COLE1": {"E. coli"},
    "ORI_PUC": {"E. coli"},
    "ORI_PBR322": {"E. coli"},
    "ORI_P15A": {"E. coli"},
    "ORI_PSC101": {"E. coli"},
    "ORI_R6K": {"E. coli"},          # ayrıca pir+ gerekir, bkz. ORI-05
    "ORI_F1": {"E. coli"},           # tek zincir üretimi, replikasyon değil
}

# E. coli'nin kendi makinesine bağlı promotörler. Gram-pozitif
# bakterilerde sigma faktörü ve ribozom bağlanma bölgesi farklıdır.
ECOLI_ONLY_PROMOTERS = {
    "PRO_T7", "PRO_T7LAC", "PRO_LAC", "PRO_TAC", "PRO_TRC",
    "PRO_ARABAD", "PRO_RHABAD", "PRO_PL",
}


def _check_host_organism(f, vec, host):
    """Konak organizması ile parçaların uyumu."""
    if not host:
        return
    organizma = host.get("organism", "E. coli")
    oris = vec.parts_of("ori")

    # --- R6K: pir+ konak şartı (E. coli içinde de geçerli)
    if any(o["part_id"] == "ORI_R6K" for o in oris):
        if not host.get("pir"):
            f.append(Finding(
                "hata", "ORI-05", "R6K orijini pir⁺ konak ister",
                f"Yapıda R6K orijini var ama {host['name']} pir geni "
                f"taşımıyor. R6K replikasyonu konağın ürettiği pi "
                f"proteinine bağlıdır; bu suşta plazmit hiç çoğalmaz ve "
                f"transformasyondan koloni çıkmaz.",
                "Plazmiti çoğaltmak için PIR1 gibi pir⁺ bir suş seçin. "
                "Eğer amacınız kromozoma entegrasyonsa bu zaten istenen "
                "durumdur: plazmit hedef suşta çoğalamaz, yalnızca "
                "entegre olanlar ayakta kalır.",
                "R6K gamma origin pir protein conditional replication",
            ))

    # --- Orijin organizmaya uyuyor mu
    if organizma != "E. coli":
        uyumsuz = [o for o in oris
                   if organizma not in ORI_ORGANISMS.get(
                       o["part_id"], {organizma})]
        if uyumsuz:
            adlar = ", ".join(o["name"].split("(")[0].strip()
                              for o in uyumsuz)
            f.append(Finding(
                "hata", "ORG-01", "Orijin bu organizmada çalışmaz",
                f"Konak {organizma}, ama yapıdaki orijin ({adlar}) "
                f"E. coli replikasyon makinesine göre tasarlanmıştır. "
                f"Plazmit bu konakta çoğalmaz - harita doğru görünse de "
                f"transformasyondan sonra kaybolur.",
                f"{organizma} için uygun bir replikon seçin. Referans "
                f"vektörler sayfasında organizmaya göre süzebilirsiniz: "
                f"Gram-pozitifler için pUB110 ya da pC194, Pseudomonas "
                f"için SEVA iskeletleri, Corynebacterium için pCG1. "
                f"Klonlamayı E. coli'de yapıp sonra aktaracaksanız iki "
                f"orijinli bir mekik vektör kurun.",
                f"shuttle vector replicon {organizma} plasmid host range",
            ))
        if not oris:
            pass   # CORE-01 zaten yakalıyor

        # --- Promotör organizmaya uyuyor mu
        yanlis_pro = [p for p in vec.parts_of("promoter")
                      if p["part_id"] in ECOLI_ONLY_PROMOTERS]
        if yanlis_pro:
            adlar = ", ".join(p["name"].split("(")[0].strip()
                              for p in yanlis_pro)
            f.append(Finding(
                "uyari", "ORG-02", "Promotör bu organizmaya göre değil",
                f"{adlar} E. coli'nin sigma faktörlerine ve düzenleyici "
                f"proteinlerine göre tasarlanmıştır. {organizma} bu "
                f"promotörü ya hiç okumaz ya da beklenmedik düzeyde "
                f"okur. Ribozom bağlanma bölgesi için de aynı sorun "
                f"geçerlidir: Shine-Dalgarno uzunluğu ve 16S rRNA "
                f"eşleşmesi türler arasında farklıdır.",
                f"{organizma} için tanımlanmış bir promotör kullanın. "
                f"Bu organizmada çalıştığı gösterilmiş bir vektörden "
                f"(Referans vektörler sayfası) promotör bölgesini "
                f"çıkarıp katalogunuza ekleyebilirsiniz.",
                f"{organizma} promoter ribosome binding site expression",
            ))

        # --- T7 sistemi E. coli dışında
        if any(p["part_id"] in ("PRO_T7", "PRO_T7LAC")
               for p in vec.parts_of("promoter")):
            f.append(Finding(
                "bilgi", "ORG-03", "T7 sistemi konakta kurulu olmalı",
                f"T7 promotörü yalnızca T7 RNA polimeraz üreten konakta "
                f"çalışır. {organizma} için bu, DE3 lizojeni gibi hazır "
                f"bir çözüm olmayabilir; polimerazı ayrıca bir plazmitle "
                f"ya da kromozoma yerleştirerek vermeniz gerekir.",
                "Konakta T7 polimeraz yoksa, organizmanın kendi "
                "polimerazıyla çalışan bir promotöre geçin.",
                f"T7 RNA polymerase expression system {organizma}",
            ))

    # --- E. coli konak + Gram-pozitif replikon
    else:
        gram_poz = [o for o in oris
                    if "E. coli" not in ORI_ORGANISMS.get(
                        o["part_id"], {"E. coli"})]
        if gram_poz:
            f.append(Finding(
                "uyari", "ORG-04", "Orijin E. coli'de çalışmayabilir",
                "Yapıdaki orijin E. coli dışı bir organizma için "
                "tanımlanmış. Klonlama adımlarını E. coli'de yapmak "
                "istiyorsanız yapıya ayrıca bir E. coli orijini "
                "ekleyin (mekik vektör).",
                "ColE1 ya da p15A orijinini ikinci replikon olarak "
                "ekleyin.",
                "shuttle vector dual origin E coli Gram-positive",
            ))
