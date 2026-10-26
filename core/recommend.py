"""Kural tabanlı vektör öneri motoru.

Kullanıcının birkaç sorusuna verdiği yanıtlardan çalışır bir bakteriyel
vektör önerir ve HER parça seçiminin gerekçesini yazar.

Tasarım ilkesi: motor asla "şunu kullan" demez, "şunu şu yüzden
kullan" der. Gerekçe olmadan öneri yapılmaz - ders bağlamında kararın
kendisi değil, kararın nedeni öğreticidir.

Sınır: bu bir kural motorudur. Sizin özel proteininizin katlanma
davranışını, gerçek toksisitesini ya da optimum indüksiyon koşulunu
bilemez; bilinen tuzaklara düşmeyen makul bir başlangıç noktası verir.
Son söz kullanıcıdadır.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

# ---------------------------------------------------------------------------
# Soru şeması
# ---------------------------------------------------------------------------

QUESTIONS: list[dict[str, Any]] = [
    {
        "key": "amac",
        "soru": "Bu vektörü ne için kuruyorsun?",
        "secenekler": {
            "ifade": "Protein üretmek (ifade)",
            "klonlama": "Sadece bir diziyi çoğaltmak / saklamak",
            "devre": "Sabit ifadeli genetik devre kurmak",
        },
        "varsayilan": "ifade",
    },
    {
        "key": "toksik",
        "soru": "Ürün konak için toksik olabilir mi?",
        "secenekler": {
            "hayir": "Hayır, sorun beklemiyorum",
            "bilmiyorum": "Bilmiyorum",
            "evet": "Evet ya da şüpheleniyorum",
        },
        "varsayilan": "bilmiyorum",
        "gorunur": lambda a: a.get("amac") != "klonlama",
    },
    {
        "key": "kaynak",
        "soru": "Gen hangi organizmadan geliyor?",
        "secenekler": {
            "ecoli": "E. coli ya da yakın bir bakteri",
            "bakteri": "Başka bir bakteri",
            "okaryot": "Ökaryot (insan, maya, bitki vb.)",
        },
        "varsayilan": "ecoli",
        "gorunur": lambda a: a.get("amac") != "klonlama",
    },
    {
        "key": "saflastirma",
        "soru": "Proteini saflaştıracak mısın?",
        "secenekler": {
            "yok": "Hayır, sadece ifade edilsin",
            "his": "Evet, standart afinite saflaştırması",
            "cozunurluk": "Evet, ama çözünürlük sorunu bekliyorum",
        },
        "varsayilan": "his",
        "gorunur": lambda a: a.get("amac") == "ifade",
    },
    {
        "key": "etiket_kesilsin",
        "soru": "Etiketin sonradan kesilmesi gerekiyor mu?",
        "secenekler": {
            "hayir": "Hayır, etiket kalabilir",
            "evet": "Evet, doğal protein lazım",
        },
        "varsayilan": "hayir",
        "gorunur": lambda a: (
            a.get("amac") == "ifade" and a.get("saflastirma") != "yok"
        ),
    },
    {
        "key": "disulfit",
        "soru": "Protein disülfit bağı içeriyor mu?",
        "secenekler": {
            "hayir": "Hayır ya da bilmiyorum",
            "evet": "Evet (antikor fragmanı, salgı proteini vb.)",
        },
        "varsayilan": "hayir",
        "gorunur": lambda a: a.get("amac") == "ifade",
    },
    {
        "key": "ikinci_plazmit",
        "soru": "Aynı hücrede ikinci bir plazmit bulunacak mı?",
        "secenekler": {
            "hayir": "Hayır, tek plazmit",
            "evet": "Evet, ko-ekspresyon ya da devre",
        },
        "varsayilan": "hayir",
    },
    {
        "key": "olcek",
        "soru": "Ne ölçekte çalışacaksın?",
        "secenekler": {
            "kucuk": "Tüp / küçük kültür (tipik ders deneyi)",
            "buyuk": "Büyük kültür, uzun inkübasyon",
        },
        "varsayilan": "kucuk",
    },
]


def default_answers() -> dict[str, str]:
    return {q["key"]: q["varsayilan"] for q in QUESTIONS}


def visible_questions(answers: dict[str, str]) -> list[dict[str, Any]]:
    """Önceki yanıtlara göre sorulması anlamlı olan soruları döndür."""
    out = []
    for q in QUESTIONS:
        gor = q.get("gorunur")
        if gor is None or gor(answers):
            out.append(q)
    return out


# ---------------------------------------------------------------------------
# Öneri
# ---------------------------------------------------------------------------


@dataclass
class Choice:
    part_id: str
    role: str
    reason: str


@dataclass
class Recommendation:
    choices: list[Choice] = field(default_factory=list)
    host_id: str = "BL21DE3"
    host_reason: str = ""
    purpose: Literal["ifade", "klonlama"] = "ifade"
    notes: list[str] = field(default_factory=list)
    protocol: list[str] = field(default_factory=list)

    @property
    def part_ids(self) -> list[str]:
        return [c.part_id for c in self.choices]


def recommend(a: dict[str, str]) -> Recommendation:
    """Yanıtlardan vektör önerisi üret."""
    amac = a.get("amac", "ifade")
    toksik = a.get("toksik", "bilmiyorum")
    kaynak = a.get("kaynak", "ecoli")
    saf = a.get("saflastirma", "his")
    kes = a.get("etiket_kesilsin", "hayir") == "evet"
    disulfit = a.get("disulfit", "hayir") == "evet"
    ikinci = a.get("ikinci_plazmit", "hayir") == "evet"
    buyuk = a.get("olcek", "kucuk") == "buyuk"

    r = Recommendation()
    add = lambda pid, role, reason: r.choices.append(Choice(pid, role, reason))

    # -----------------------------------------------------------------
    # 1) Yalnızca klonlama
    # -----------------------------------------------------------------
    if amac == "klonlama":
        r.purpose = "klonlama"
        add("PRO_LAC", "Promotör",
            "lac promotörü, lacZ-alfa parçasını sürerek mavi/beyaz tarama "
            "yapmanı sağlar. Burada amaç protein üretmek değil, insert "
            "alan kolonileri gözle ayırt edebilmek.")
        add("RPT_LACZA", "Tarama",
            "lacZ-alfa, konaktaki eksik lacZ-delta-M15 ile birleşip aktif "
            "enzim oluşturur. Insert MCS'e girince çerçeve bozulur ve "
            "koloni beyaz kalır - doğru klonu plakta görürsün.")
        add("MCS_PUC19", "Klonlama bölgesi",
            "pUC19 MCS'i 13 enzim için 10 bölge sunar ve lacZ-alfa "
            "çerçevesinin tam içindedir; mavi/beyaz tarama ancak böyle "
            "çalışır.")
        add("ORI_PUC", "Orijin",
            "Yüksek kopya (500-700). Amaç çoğaltma olduğu için bol plazmit "
            "verimi istiyorsun; ifade yükü olmadığından metabolik baskı "
            "sorun değil.")
        add("MRK_BLA" if not buyuk else "MRK_KAN", "Marker",
            "Ampisilin, kısa klonlama deneylerinin standardıdır ve ucuzdur."
            if not buyuk else
            "Kanamisin seçtim: uzun inkübasyonda beta-laktamaz besiyerindeki "
            "ampisilini parçalar ve plazmitsiz satellit koloniler çıkar.")
        r.host_id = "DH5A"
        r.host_reason = (
            "DH5-alfa klonlamanın standardı: recA1 eklenmiş insertin "
            "rekombinasyonla bozulmasını, endA1 ise plazmit saflaştırmada "
            "nükleaz bulaşmasını engeller. lacZ-delta-M15 taşıdığı için "
            "mavi/beyaz tarama burada çalışır."
        )
        r.protocol = [
            "Plaklara X-gal (40 µg/mL) ve IPTG (0.1 mM) ekle - yoksa mavi "
            "renk oluşmaz ve tarama yapılamaz.",
            "Beyaz kolonileri seç, koloni PCR ya da sindirimle doğrula.",
        ]
        r.notes.append(
            "Denetim bu yapıda 'LacI baskılaması yetersiz kalabilir' "
            "uyarısı verecek: yüksek kopyalı pUC orijini, kromozomdaki tek "
            "lacI kopyasının ürettiği baskılayıcıyı tüketir. Burada bu bir "
            "sorun DEĞİL, istenen durumdur - lacZ-alfa'nın bir miktar "
            "ifade edilmesi zaten mavi rengin oluşması için gerekir."
        )
        if ikinci:
            r.notes.append(
                "İkinci plazmit kullanacaksan bu yapının pUC orijini "
                "ColE1 ailesindendir; diğer plazmit p15A ya da pSC101 "
                "orijinli olmalı, yoksa ikisi aynı hücrede kararlı "
                "tutulamaz."
            )
        return r

    # -----------------------------------------------------------------
    # 2) Sabit ifadeli devre
    # -----------------------------------------------------------------
    if amac == "devre":
        r.purpose = "ifade"
        add("PRO_J23105", "Promotör",
            "Sabit (konstitütif) orta güçte sentetik promotör. Devrede "
            "indükleyici kimyasala bağlı kalmak istemezsin; J23105 "
            "referans J23100'ün yaklaşık dörtte biri gücünde, bu da "
            "devre elemanları için genelde daha kararlı davranış verir.")
        add("RBS_B0032", "RBS",
            "Orta güçte RBS. Devrelerde maksimum ifade değil, öngörülebilir "
            "ve dengeli ifade istenir; çok güçlü RBS konağı yorar ve "
            "devrenin dinamiğini bozar.")
        add("TER_B0015", "Terminatör",
            "B0015 çift terminatörü konağın kendi RNA polimerazını verimli "
            "durdurur. Sabit promotör kullandığın için T7 terminatörü "
            "burada yetersiz kalırdı.")
        add("ORI_P15A", "Orijin",
            "p15A, ColE1 ailesiyle UYUMLUDUR. Devrelerde ikinci bir "
            "plazmit neredeyse her zaman gerekir; bu orijin onunla "
            "yan yana durabilir." if ikinci else
            "p15A orta-düşük kopya (10-12). Sabit promotör sürekli açık "
            "olduğu için yüksek kopya gereksiz metabolik yük yaratırdı.")
        add("MRK_CAT", "Marker",
            "Kloramfenikol. Devre kurulumlarında ampisilin ve kanamisin "
            "genelde diğer modüllere ayrılır; üçüncü bağımsız seçilim "
            "katmanı işini görür.")
        r.host_id = "DH5A"
        r.host_reason = (
            "DH5-alfa: devre sabit promotörle sürüldüğü için T7 polimeraza "
            "ihtiyaç yok, ve klonlama suşu olması plazmiti kararlı tutar."
        )
        r.notes.append(
            "Devreye bir reporter eklemek istersen Tasarımcı sayfasından "
            "sfGFP ya da mCherry seçebilirsin; promotör gücünü ölçmenin "
            "en pratik yolu budur."
        )
        return r

    # -----------------------------------------------------------------
    # 3) Protein ifadesi - asıl karar ağacı
    # -----------------------------------------------------------------
    r.purpose = "ifade"

    # --- Promotör + düzenleyici + konak ailesi
    if toksik == "evet":
        add("REG_ARAC", "Düzenleyici",
            "AraC, P-BAD'ın hem baskılayıcısı hem aktivatörüdür. "
            "Olmadan promotör indüksiyona hiç yanıt vermez, bu yüzden "
            "zorunludur.")
        add("PRO_ARABAD", "Promotör",
            "Ürün toksik olabileceği için T7 yerine P-BAD seçtim. "
            "Arabinoz dozunu 0.00002-0.2 % arasında değiştirerek ifadeyi "
            "kademeli ayarlayabilirsin; T7 ise ya kapalı ya tam açıktır ve "
            "toksik üründe hücreyi öldürür.")
        add("RBS_B0032", "RBS",
            "Orta güçte RBS. Toksik üründe çeviriyi de kısmak, promotörü "
            "kısmak kadar işe yarar; güçlü RBS ile düşük indüksiyonun "
            "kazancını kaybedersin.")
        add("TER_B0015", "Terminatör",
            "P-BAD konağın kendi polimerazıyla okunur, bu yüzden T7 "
            "terminatörü değil B0015 çift terminatörü gerekir.")
        host = "TOP10"
        host_reason = (
            "TOP10: araD139 ve ara operonu delesyonu taşır, yani arabinozu "
            "karbon kaynağı olarak tüketmez. İndüksiyon seviyesi kültür "
            "boyunca sabit kalır - normal suşta arabinoz tükendikçe ifade "
            "düşer ve deney tekrarlanamaz hale gelir."
        )
    else:
        add("REG_LACI", "Düzenleyici",
            "lacI'yi plazmite koydum. Konakta lacI-q olsa bile plazmit "
            "üzerindeki kopya, indüksiyon öncesi kaçak ifadeyi belirgin "
            "düşürür.")
        add("PRO_T7LAC", "Promotör",
            "T7lac: bakteride elde edebileceğin en yüksek ifade. T7 "
            "polimeraz konağın kendi polimerazından çok daha hızlıdır, "
            "altına eklenen lac operatörü de indüksiyon öncesi kaçağı "
            "bastırır.")
        add("RBS_PET", "RBS",
            "T7 fag gen 10 bölgesinden türetilmiş çok güçlü RBS. "
            "pET-28a MCS'i ile birlikte tasarlandığı için RBS-ATG "
            "mesafesi zaten doğru; ayrıca ayarlaman gerekmez.")
        add("TER_T7", "Terminatör",
            "T7 promotörünü T7 terminatörü kapatır; ikisi aynı polimeraza "
            "göre tasarlanmıştır.")
        host = "BL21DE3"
        host_reason = (
            "BL21(DE3): T7 promotörünün çalışması için konağın T7 RNA "
            "polimeraz taşıması ZORUNLUDUR, bunu DE3 lizojeni sağlar. "
            "Ayrıca ompT ve lon proteazları eksik olduğu için rekombinant "
            "protein daha az parçalanır."
        )
        if toksik == "bilmiyorum":
            host = "BL21DE3PLYSS"
            host_reason = (
                "BL21(DE3) pLysS: toksisiteden emin olmadığın için "
                "güvenli tarafı seçtim. pLysS plazmiti T7 lizozimi "
                "üretir, bu da T7 polimerazı doğal olarak inhibe edip "
                "indüksiyon öncesi kaçak ifadeyi düşürür. Ürün zararsız "
                "çıkarsa düz BL21(DE3)'e geçip verimi artırabilirsin."
            )

    # --- Konağı özel durumlar geçersiz kılabilir
    if disulfit:
        host = "SHUFFLE"
        host_reason = (
            "SHuffle T7: disülfit bağı gerektiğini söyledin. Bu suşun "
            "sitoplazması oksitleyici hale getirilmiştir (trxB gor "
            "mutasyonu + sitoplazmik DsbC), yani disülfit bağları "
            "periplazmaya hedefleme olmadan doğru kurulabilir. Normal "
            "BL21'in sitoplazması indirgeyicidir ve bu bağlar oluşmaz."
        )
    elif kaynak == "okaryot":
        host = "ROSETTA"
        host_reason = (
            "Rosetta(DE3): ökaryot kaynaklı gen dediğin için seçtim. "
            "pRARE plazmiti E. coli'de nadir olan kodonların (AGA/AGG, "
            "CGG, AUA, CUA, CCC, GGA) tRNA'larını sağlar. Bu tRNA'lar "
            "olmadan ribozom o kodonlarda duraklar, kısalmış ürünler ve "
            "düşük verim görülür."
        )

    r.host_id = host
    r.host_reason = host_reason

    # --- Etiket ve kesim bölgesi
    if saf == "his":
        add("TAG_HIS6", "Etiket",
            "6xHis: sadece 6 aminoasit olduğu için proteinin katlanmasını "
            "ve işlevini nadiren etkiler, bu yüzden genelde kesmek bile "
            "gerekmez. Ni-NTA ile tek adımda saflaştırırsın.")
    elif saf == "cozunurluk":
        add("TAG_SUMO", "Etiket",
            "His-SUMO: çözünürlük sorunu beklediğin için seçtim. Hem "
            "çözünürlüğü artırır hem de SUMO proteazı (Ulp1) füzyonu "
            "hedefe HİÇ artık bırakmadan keser - MBP ve GST'de kesimden "
            "sonra birkaç aminoasit kalır. Kendi His etiketini de taşıdığı "
            "için ayrıca etiket eklemeye gerek yok.")

    if kes and saf == "his":
        add("PRS_TEV", "Kesim bölgesi",
            "TEV proteaz bölgesi: doğal protein istediğin için ekledim. "
            "Çok özgüldür ve kesimden sonra hedefin N-terminalinde sadece "
            "tek bir glisin kalır. Trombin bölgesine göre yan kesim riski "
            "çok daha düşüktür.")

    # --- MCS
    if toksik != "evet":
        add("MCS_PET28", "Klonlama bölgesi",
            "pET-28a MCS'i: NcoI/NdeI'den XhoI'ye kadar sekiz enzim "
            "sunar ve T7lac + pET RBS ile aynı bağlamda tasarlanmıştır.")
    else:
        add("MCS_PUC19", "Klonlama bölgesi",
            "pUC19 MCS'i: P-BAD kasetine pET MCS'ini takmanın bir anlamı "
            "yok (onun çerçevesi pET bağlamına göre kurulu). Bu MCS 10 "
            "bağımsız bölge sunar.")

    # --- Orijin
    # Konak zaten p15A orijinli bir yardımcı plazmit taşıyor olabilir
    # (pRARE/pLysS). O durumda p15A önermek üçüncü bir çakışma yaratır.
    konakta_p15a = host in ("ROSETTA", "BL21DE3PLYSS")

    if ikinci and konakta_p15a:
        yardimci = "pRARE" if host == "ROSETTA" else "pLysS"
        add("ORI_COLE1", "Orijin",
            f"ColE1: burada üç plazmit var. Seçtiğim konak {yardimci} "
            "yardımcı plazmitini zaten p15A orijiniyle taşıyor, bu yüzden "
            "senin vektörüne p15A veremem - ikisi aynı uyumsuzluk "
            f"grubundan olur ve biri kaybolur. ColE1, {yardimci}'nin "
            "p15A'sıyla uyumludur.")
        r.notes.append(
            f"Üç plazmitli bir kurulum yapıyorsun: konaktaki {yardimci} "
            "(p15A), bu vektör (ColE1) ve senin ikinci plazmitin. "
            "Üçüncüsü için pSC101 orijinini seç - ColE1 ve p15A'nın "
            "ikisiyle de uyumlu olan tek seçenek odur. Üç farklı "
            "antibiyotik de gerekir: kloramfenikol yardımcı plazmitin, "
            "kanamisin bu vektörün, üçüncüsüne spektinomisin ya da "
            "ampisilin kalır."
        )
    elif ikinci:
        add("ORI_P15A", "Orijin",
            "p15A: ikinci bir plazmit kullanacağını söyledin. Bu orijin "
            "ColE1 ailesiyle UYUMLUDUR, yani iki plazmit aynı hücrede "
            "kararlı durabilir. pUC ya da ColE1 seçseydim biri kaybolurdu."
            " İkinci plazmitini ColE1/pUC ailesinden seçebilirsin.")
    elif toksik == "evet":
        add("ORI_PBR322", "Orijin",
            "pBR322 (15-20 kopya): toksik üründe düşük kopya sayısı "
            "kaçak ifadenin mutlak miktarını da düşürür. Yüksek kopya "
            "orijin, promotörü ne kadar sıkı kontrol edersen et, "
            "toplam kaçağı yukarı çeker.")
    else:
        add("ORI_COLE1", "Orijin",
            "ColE1 (15-20 kopya): pET ailesinin orijini. İfade "
            "vektörlerinde yüksek kopya (pUC) tercih edilmez - güçlü "
            "promotörle birleşince hücre ribozom kapasitesini aşar, "
            "büyüme durur ve inklüzyon cismi oluşumu artar.")

    # --- Marker (konakta yerleşik plazmitle çakışmamalı)
    konakta_cm = host in ("ROSETTA", "BL21DE3PLYSS")
    if konakta_cm:
        marker, marker_reason = "MRK_KAN", (
            f"Kanamisin: seçtiğim konak ({'pRARE' if host == 'ROSETTA' else 'pLysS'}) "
            "yardımcı plazmitini kloramfenikol direnciyle taşıyor. "
            "Markerı kloramfenikol seçseydim iki plazmiti birbirinden "
            "ayırt edemezdin ve yardımcı plazmitin kaybını fark etmezdin."
        )
    elif buyuk:
        marker, marker_reason = "MRK_KAN", (
            "Kanamisin: büyük ölçekte çalışacağını söyledin. "
            "Beta-laktamaz periplazmaya salgılanıp besiyerindeki "
            "ampisilini parçalar; yoğun kültürde seçilim baskısı kaybolur "
            "ve plazmitsiz hücreler kültürü ele geçirir. Kanamisin "
            "fosfotransferazı hücre içinde kalır, bu sorun olmaz."
        )
    elif ikinci:
        marker, marker_reason = "MRK_KAN", (
            "Kanamisin: ikinci plazmit için ampisilin ya da kloramfenikolü "
            "boşta bıraktım; iki plazmitin markerı farklı olmak zorundadır, "
            "yoksa ikisinin de korunduğunu doğrulayamazsın."
        )
    else:
        marker, marker_reason = "MRK_KAN", (
            "Kanamisin: modern ifade vektörlerinin standardı. Ampisiline "
            "göre seçilim baskısı kültür boyunca daha kararlı kalır ve "
            "satellit koloni sorunu yoktur."
        )
    add(marker, "Marker", marker_reason)

    # -----------------------------------------------------------------
    # Protokol ve notlar
    # -----------------------------------------------------------------
    if toksik == "evet":
        r.protocol = [
            "Önce arabinoz doz taraması yap: 0.00002, 0.0002, 0.002, "
            "0.02, 0.2 % - toksik proteinde en yüksek doz en iyi verimi "
            "vermez.",
            "İndüksiyonu OD600 ~0.5'te başlat ve 18-25 °C'de gece boyu "
            "sürdür; düşük sıcaklık katlanmayı iyileştirir.",
            "Klonlamayı DH5-alfa'da yap, doğrulanmış plazmiti sonra "
            "ifade suşuna aktar.",
        ]
    else:
        r.protocol = [
            "İndüksiyon: OD600 0.4-0.6 arasında 0.1-1 mM IPTG.",
            "Çözünürlük sorunu görürsen sıcaklığı 37 °C'den 18-25 °C'ye "
            "düşür ve IPTG'yi 0.1 mM'ye indir - en sık işe yarayan iki "
            "müdahale bunlar.",
            "Klonlamayı DH5-alfa'da yap: ifade suşları recA+ olduğu için "
            "plazmiti uzun süre kararlı tutmaz.",
        ]

    if kaynak == "okaryot" and host != "ROSETTA":
        r.notes.append(
            "Gen ökaryot kaynaklı ama konağı başka bir gerekçeyle seçtim. "
            "Verim düşük çıkarsa Dizi analizi sayfasında nadir kodon "
            "oranına bak; %8'in üzerindeyse diziyi kodon optimize ettir."
        )
    if saf == "cozunurluk" and kes:
        r.notes.append(
            "SUMO zaten artık bırakmadan kesildiği için ayrıca TEV "
            "bölgesi eklemedim - iki kesim bölgesi gereksiz karmaşa olurdu."
        )
    if saf == "yok":
        r.notes.append(
            "Etiket eklemedim. Proteini jelde göreceksen sorun yok, ancak "
            "saflaştırma ya da Western blot gerekirse sonradan etiket "
            "eklemek yeni bir klonlama turu demektir - şimdiden "
            "düşünmekte fayda var."
        )
    r.notes.append(
        "Bu öneri bir kural motorundan çıktı: senin proteininin katlanma "
        "davranışını ve gerçek toksisitesini bilemez. Tasarımcı sayfasına "
        "aktarıp istediğin parçayı değiştirebilirsin; denetim her "
        "değişiklikten sonra yeniden çalışır."
    )
    return r
