"""Bakteriyel vektör parça kataloğu (başlangıç verisi).

Kapsam bilinçli olarak SADECE bakteriyel sistemlerdir (ağırlıklı E. coli).
Memeli / bitki / böcek / maya elemanları buraya girmez.

Dizi durumu (`seq_status`):
    full  -> gerçek, doğrulanabilir dizi gömülü. Kısa düzenleyici elemanlar,
             etiketler, proteaz bölgeleri, MCS'ler. Klonlama simülasyonunda
             birebir kullanılır.
    stub  -> dizi gömülü DEĞİL, sadece uzunluk + erişim numarası var.
             Uzun CDS'ler (direnç genleri, reporter'lar, ori bölgeleri) için.
             Montajda doğru uzunlukta N bloğu olarak yer alır ve açıkça
             işaretlenir; kullanıcı `Parça kataloğu` sayfasından ya da bir
             GenBank dosyası içe aktararak gerçek diziyi doldurabilir.

Bu ayrım kasıtlı: ezberden uzun gen dizisi yazmak sessiz hatalar üretir.
Her kaydın `accession` alanı dizinin nereden çekileceğini söyler, `pmid`
alanı ise literatür modülünün kaynakça üretmesini sağlar.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Kategoriler
# ---------------------------------------------------------------------------

CATEGORIES = [
    "promoter",
    "operator",
    "rbs",
    "cds",
    "tag",
    "protease_site",
    "linker",
    "terminator",
    "ori",
    "marker",
    "regulator",
    "mcs",
    "reporter",
    "device",
]

CATEGORY_TR = {
    "promoter": "Promotör",
    "operator": "Operatör",
    "rbs": "Ribozom bağlanma bölgesi",
    "cds": "Kodlayan dizi",
    "tag": "Afinite/füzyon etiketi",
    "protease_site": "Proteaz kesim bölgesi",
    "linker": "Bağlayıcı",
    "terminator": "Terminatör",
    "ori": "Replikasyon orijini",
    "marker": "Seçilim markerı",
    "regulator": "Düzenleyici protein",
    "mcs": "Çoklu klonlama bölgesi",
    "reporter": "Reporter",
    "device": "Bileşik parça (devre)",
}

# Harita çiziminde kategori renkleri (hem açık hem koyu temada okunur tonlar)
CATEGORY_COLOR = {
    "promoter": "#22C55E",
    "operator": "#16A34A",
    "rbs": "#EAB308",
    "cds": "#3B82F6",
    "tag": "#8B5CF6",
    "protease_site": "#A78BFA",
    "linker": "#A1A1AA",
    "terminator": "#EF4444",
    "ori": "#F97316",
    "marker": "#E11D48",
    "regulator": "#0891B2",
    "mcs": "#71717A",
    "reporter": "#2563EB",
    "device": "#7C3AED",
}


def _p(**kw):
    """Parça kaydı kur, eksik alanları varsayılanla doldur."""
    rec = {
        "part_id": None,
        "name": None,
        "category": None,
        "sequence": "",
        "seq_status": "full",
        "length": None,
        "strength": None,
        "induction": None,
        "requires": [],
        "incompatible_with": [],
        "copy_number": None,
        "inc_group": None,
        "notes": "",
        "accession": "",
        "pmid": "",
        "source": "seed",
    }
    rec.update(kw)
    if rec["seq_status"] == "full":
        rec["sequence"] = rec["sequence"].upper().replace(" ", "").replace("\n", "")
        rec["length"] = len(rec["sequence"])
    return rec


# ---------------------------------------------------------------------------
# PROMOTÖRLER
# ---------------------------------------------------------------------------

PROMOTERS = [
    _p(
        part_id="PRO_T7",
        name="T7 promotör",
        category="promoter",
        sequence="TAATACGACTCACTATAGGG",
        strength="çok güçlü",
        induction="IPTG (konaktaki lacUV5-T7 RNA pol kaseti üzerinden)",
        requires=["T7 RNA polimeraz (konak genomunda, örn. DE3 lizojeni)"],
        notes=(
            "Konak T7 RNA polimeraz taşımıyorsa HİÇ transkripsiyon olmaz. "
            "Klonlama suşlarında (DH5-alfa, TOP10) sessiz kalır - bu aslında "
            "toksik ürünleri klonlarken avantajdır."
        ),
        accession="T7 fag genomu: NC_001604",
        pmid="6283312",
    ),
    _p(
        part_id="PRO_T7LAC",
        name="T7lac promotör (T7 + lac operatör)",
        category="promoter",
        sequence="TAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCC",
        strength="çok güçlü, baskılanabilir",
        induction="IPTG",
        requires=[
            "T7 RNA polimeraz (konak genomunda)",
            "LacI (konak lacI-q ya da plazmit üzerinde)",
        ],
        notes=(
            "pET serisinin standart promotörü. T7 promotörünün hemen altına "
            "yerleştirilmiş lac operatörü, indüksiyon öncesi kaçak ifadeyi "
            "düşürür. Sıkı kontrol gerekiyorsa pLysS/pLysE konağı eklenir."
        ),
        accession="pET-28a: EU822321 benzeri",
        pmid="2199796",
    ),
    _p(
        part_id="PRO_LAC",
        name="lac promotör (pUC tipi, lacUV5)",
        category="promoter",
        sequence=(
            "TTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACAGGA"
        ),
        strength="orta",
        induction="IPTG veya laktoz (glukoz varlığında katabolik baskılanma)",
        requires=["LacI"],
        notes=(
            "Promotör + operatör + zayıf RBS bağlamını birlikte taşır. "
            "pUC19/pBluescript ailesinin lacZ-alfa ifadesini süren bölge. "
            "Glukoz içeren besiyerinde cAMP-CAP aktivasyonu düştüğü için "
            "beklenenden zayıf çalışır."
        ),
        accession="pUC19: L09137",
        pmid="2985470",
    ),
    _p(
        part_id="PRO_TAC",
        name="tac promotör (trp -35 / lac -10 hibrit)",
        category="promoter",
        sequence="TTGACAATTAATCATCGGCTCGTATAATGTGTGGA",
        strength="güçlü",
        induction="IPTG",
        requires=["LacI"],
        notes=(
            "trp promotörünün -35 bölgesi ile lac promotörünün -10 bölgesinin "
            "birleşimi; her ikisinden de güçlü. Konağın kendi RNA "
            "polimerazını kullandığı için T7 sistemine gerek yoktur - "
            "yani DH5-alfa gibi suşlarda da çalışır."
        ),
        accession="pDR540 / pGEX serisi",
        pmid="6269053",
    ),
    _p(
        part_id="PRO_TRC",
        name="trc promotör",
        category="promoter",
        sequence="TTGACAATTAATCATCCGGCTCGTATAATGTGTGGA",
        strength="güçlü",
        induction="IPTG",
        requires=["LacI"],
        notes=(
            "tac ile neredeyse aynı; -35 ve -10 arası boşluk 17 bp "
            "(tac'ta 16 bp). pTrc99A ve birçok ifade vektöründe kullanılır."
        ),
        accession="pTrc99A: U13872",
        pmid="2991595",
    ),
    _p(
        part_id="PRO_ARABAD",
        name="araBAD promotör (P-BAD)",
        category="promoter",
        sequence=(
            "ACTTTTCATACTCCCACCATTCAGAGAAGAAACCAATTGTCCATATTGCATCAGACATTG"
            "CCGTCACTGCGTCTTTTACTGGCTCTTCTCGCTAACCCAACCGGTAACCCCGCTTATTAA"
            "AAGCATTCTGTAACAAAGCGGGACCAAAGCCATGACAAAAACGCGTAACAAAAGTGTCTA"
            "TAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCGTCACACTTTGCTATGCCAT"
            "AGCATTTTTATCCATAAGATTAGCGGATCCTACCTGACGCTTTTTATCGCAACTCTCTAC"
            "TGTTTCTCCATACCCGTTTTTTTGGATGGAGTGAAACG"
        ),
        strength="ayarlanabilir (çok geniş aralık)",
        induction="L-arabinoz (0.00002 - 0.2 %)",
        requires=["AraC düzenleyicisi", "araE/araFGH taşıyıcı (homojen yanıt için)"],
        notes=(
            "Doz-yanıt eğrisi çok geniş: toksik ya da katlanması zor "
            "proteinlerde T7'ye göre çok daha kontrollü. araI1/araI2/araO2 "
            "operatörlerini ve CAP bölgesini içerir, bu yüzden tek bir "
            "'konsensüs' dizi olarak kısaltılamaz - gerçek diziyi erişim "
            "numarasından çekin. Hücreler arası tekdüzelik için "
            "arabinoz taşıyıcısı sabit ifade edilmelidir."
            "\n\n"
            "Dizi: araB-araC genler arası bölgesinin tamamı; araI1/araI2, CAP "
            "bölgesi, -35/-10 ve araB'nin doğal RBS'i dahil, 3' ucu araB'nin "
            "ATG'sinden hemen önce bitiyor. pBAD vektörlerinin bilinen "
            "'ACTGTTTCTCCATACCCGTTTTTTTGG' ucu dizide birebir doğrulandı. "
            "Katalogdaki 130 bp tahmini, gerçek bölge 338 bp."
        ),
        accession="NC_000913.3 (E. coli K-12 MG1655): complement(70049..70386)",
        pmid="7592312",
    ),
    _p(
        part_id="PRO_RHABAD",
        name="rhaBAD promotör (P-rhaBAD)",
        category="promoter",
        sequence=(
            "ACTGGCCTCCTGATGTCGTCAACACGGCGAAATAGTAATCACGAGGTCAGGTTCTTACCT"
            "TAAATTTTCGACGGAAAACCACGTAAAAAACGTCGATTTTTCAAGATACAGCGTGAATTT"
            "TCAGGAAATGCGGTGAGCATCACATCACCACAATTCAGCAAATTGTGAACATCATCACGT"
            "TCATCTTTCCCTGGTTGCCAATGGCCCATTTTCCTGTCAGTAACGAGAAGGTCGCGAATT"
            "CAGGCGCTTTTTAGACTGGTCGTAATGAAATTCAGCAGGATCACATT"
        ),
        strength="ayarlanabilir",
        induction="L-ramnoz",
        requires=["RhaR ve RhaS düzenleyicileri"],
        notes=(
            "araBAD'a alternatif; ramnoz E. coli tarafından karbon kaynağı "
            "olarak daha yavaş tüketildiği için indüksiyon daha uzun süre "
            "kararlı kalır. Toksik protein ifadesinde iyi seçenek."
            "\n\n"
            "Dizi: rhaB-rhaS genler arası bölgesinin tamamı (RhaS ve CRP "
            "bölgeleri, -35/-10, rhaB RBS'i); 3' ucu rhaB'nin ATG'sinden hemen "
            "önce bitiyor. Katalogdaki 160 bp tahmini, gerçek bölge 287 bp."
        ),
        accession="NC_000913.3 (E. coli K-12 MG1655): complement(4097449..4097735)",
        pmid="9463464",
    ),
    _p(
        part_id="PRO_TET",
        name="tetA promotör / tet operatör (P-tet)",
        category="promoter",
        sequence="TCCCTATCAGTGATAGAGATTGACATCCCTATCAGTGATAGAGATAATGAGCAC",
        strength="orta-güçlü",
        induction="anhidrotetrasiklin (aTc), ~100 ng/mL",
        requires=["TetR baskılayıcısı"],
        notes=(
            "İki tetO2 operatör kopyası içerir. Katabolik baskılanmadan "
            "etkilenmez, bu yüzden glukozlu besiyerinde lac tabanlı "
            "sistemlerden daha güvenilirdir. aTc antibiyotik etkisi "
            "göstermeyecek kadar düşük dozda kullanılır."
        ),
        accession="pASK75 / Tn10 tetA bölgesi",
        pmid="9092630",
    ),
    _p(
        part_id="PRO_J23100",
        name="J23100 (sabit, güçlü sentetik promotör)",
        category="promoter",
        sequence="TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC",
        strength="güçlü (Anderson serisi referansı, 1.00)",
        induction="yok - sabit (konstitütif)",
        notes=(
            "Anderson promotör kütüphanesinin referans üyesi. İndükleyici "
            "gerekmez; sabit düşük/orta ifade isteyen devrelerde ve "
            "düzenleyici proteinleri (LacI, TetR, AraC) sürmek için ideal. "
            "Seride J23119 en güçlü, J23113 en zayıf üyedir."
        ),
        accession="iGEM Registry: BBa_J23100",
        pmid="",
    ),
    _p(
        part_id="PRO_J23105",
        name="J23105 (sabit, orta sentetik promotör)",
        category="promoter",
        sequence="TTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGC",
        strength="orta (~0.24 x J23100)",
        induction="yok - sabit",
        notes=(
            "J23100 ile aynı iskelet, -35 ve -10 bölgelerinde tek baz "
            "farklarıyla zayıflatılmış. Düzenleyici protein ifadesi için "
            "J23100'den daha uygun olabilir."
        ),
        accession="iGEM Registry: BBa_J23105",
        pmid="",
    ),
    _p(
        part_id="PRO_PL",
        name="lambda P-L promotör (sıcaklıkla indüklenen)",
        category="promoter",
        sequence=(
            "GATAACCATCTGCGGTGATAAATTATCTCTGGCGGTGTTGACATAAATACCACTGGCGGT"
            "GATACTGAGCACATCAGCAG"
        ),
        strength="çok güçlü",
        induction="sıcaklık kaydırması 30 C -> 42 C (cI857 ile)",
        requires=["cI857 sıcaklığa duyarlı baskılayıcı"],
        notes=(
            "Kimyasal indükleyici gerektirmez - büyük ölçekli "
            "fermentasyonda maliyet avantajı. 42 C'ye çıkmak aynı zamanda "
            "ısı şoku yanıtını ve proteazları tetikler; bu yüzden hassas "
            "proteinlerde dezavantaj olabilir."
            "\n\n"
            "Dizi: oL3/oL2/oL1 operatörleri kaydın kendi anotasyonundan; pL "
            "yönünde okununca -35 = TTGACA ve 17 bp aralıkla -10 = GATACT "
            "çıkıyor, +1 lambda 35582'ye denk geliyor."
        ),
        accession="NC_001416 (lambda fajı): complement(35575..35654)",
        pmid="6292855",
    ),
]

# ---------------------------------------------------------------------------
# OPERATÖRLER
# ---------------------------------------------------------------------------

OPERATORS = [
    _p(
        part_id="OP_LACO1",
        name="lac operatör O1",
        category="operator",
        sequence="AATTGTGAGCGGATAACAATT",
        requires=["LacI"],
        notes=(
            "Simetrik olmayan doğal operatör. LacI tetramerinin bağlandığı "
            "ana bölge. İdeal simetrik varyant (lacOsym) daha sıkı bağlanır "
            "ama indüksiyon da zorlaşır."
        ),
        accession="E. coli lac operonu",
        pmid="4587270",
    ),
    _p(
        part_id="OP_TETO2",
        name="tet operatör O2",
        category="operator",
        sequence="TCCCTATCAGTGATAGAGA",
        requires=["TetR"],
        notes="19 bp palindromik operatör. Genelde 2 kopya tandem kullanılır.",
        accession="Tn10",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# RİBOZOM BAĞLANMA BÖLGELERİ
# ---------------------------------------------------------------------------

RBS_PARTS = [
    _p(
        part_id="RBS_B0034",
        name="B0034 (güçlü, sentetik)",
        category="rbs",
        sequence="AAAGAGGAGAAA",
        strength="güçlü (referans)",
        notes=(
            "iGEM'in fiili standardı. Shine-Dalgarno çekirdeği AGGAGG'e "
            "yakın. Başlama kodonuna kadar 6-8 nt boşluk bırakılmalı."
        ),
        accession="iGEM Registry: BBa_B0034",
        pmid="",
    ),
    _p(
        part_id="RBS_B0032",
        name="B0032 (orta, sentetik)",
        category="rbs",
        sequence="TCACACAGGAAAG",
        strength="orta (~0.3 x B0034)",
        notes="İfadeyi kasıtlı düşürmek, çözünürlüğü artırmak için kullanılır.",
        accession="iGEM Registry: BBa_B0032",
        pmid="",
    ),
    _p(
        part_id="RBS_B0031",
        name="B0031 (zayıf, sentetik)",
        category="rbs",
        sequence="TCACACAGGAAACC",
        strength="zayıf (~0.07 x B0034)",
        notes=(
            "Toksik ya da inklüzyon cismi yapan proteinlerde ifadeyi "
            "kısmak için. Düzenleyici proteinler için de uygun."
        ),
        accession="iGEM Registry: BBa_B0031",
        pmid="",
    ),
    _p(
        part_id="RBS_PET",
        name="pET tipi RBS (gen 10 kaynaklı)",
        category="rbs",
        sequence="AAGGAGATATACAT",
        strength="çok güçlü",
        notes=(
            "T7 fag gen 10 bölgesinden türetilmiş. Dizinin sonundaki ATG "
            "doğrudan NcoI bölgesinin (CCATGG) parçasıdır - yani RBS ile "
            "başlama kodonu arası boşluk sabittir, ayrı ayarlanmaz."
        ),
        accession="pET-28a bölgesi",
        pmid="2199796",
    ),
]

# ---------------------------------------------------------------------------
# ETİKETLER
# ---------------------------------------------------------------------------

TAGS = [
    _p(
        part_id="TAG_HIS6",
        name="6xHis etiketi",
        category="tag",
        sequence="CATCATCATCATCATCAC",
        notes=(
            "18 bp / 6 aa. Ni-NTA ya da Co-TALON ile saflaştırma. Hem N- hem "
            "C-terminalde çalışır. Küçük olduğu için genelde kesilmesi "
            "gerekmez. Denatüre koşullarda da (8 M üre) bağlanır - "
            "inklüzyon cisimlerinden saflaştırmada tek avantajlı etiket."
        ),
        accession="",
        pmid="3125348",
    ),
    _p(
        part_id="TAG_HIS10",
        name="10xHis etiketi",
        category="tag",
        sequence="CATCATCATCATCATCATCATCATCATCAC",
        notes=(
            "Daha güçlü bağlanma, daha az kaçak. Düşük ifadeli ya da "
            "kısmen gömülü etiketlerde 6xHis yetmediğinde kullanılır."
        ),
        accession="",
        pmid="",
    ),
    _p(
        part_id="TAG_FLAG",
        name="FLAG etiketi (DYKDDDDK)",
        category="tag",
        sequence="GACTACAAAGACGATGACGACAAG",
        notes=(
            "Antikor tabanlı saflaştırma ve Western blot tespiti. "
            "Hidrofilik olduğu için terminal bölgede yüzeyde kalır. "
            "Enterokinaz bölgesini (DDDDK) kendi içinde taşır."
        ),
        accession="",
        pmid="",
    ),
    _p(
        part_id="TAG_STREPII",
        name="Strep-tag II (WSHPQFEK)",
        category="tag",
        sequence="TGGAGCCACCCGCAGTTCGAAAAA",
        notes=(
            "Strep-Tactin ile saflaştırma; desthiobiyotin ile çok yumuşak "
            "koşullarda elüsyon. Metal iyonu kullanmadığı için metalloprotein "
            "ve metal-duyarlı enzimlerde His etiketine tercih edilir."
        ),
        accession="",
        pmid="8771702",
    ),
    _p(
        part_id="TAG_T7",
        name="T7 etiketi (MASMTGGQQMG)",
        category="tag",
        sequence="ATGGCTAGCATGACTGGTGGACAGCAAATGGGT",
        notes=(
            "N-terminal tespit etiketi; T7 fag gen 10 lider peptidi. "
            "Saflaştırma için değil, tespit ve çeviri başlangıcını "
            "güçlendirmek için kullanılır."
        ),
        accession="pET-3 serisi",
        pmid="",
    ),
    _p(
        part_id="TAG_GST",
        name="GST füzyon etiketi",
        category="tag",
        sequence=(
            "ATGTCCCCTATACTAGGTTATTGGAAAATTAAGGGCCTTGTGCAACCCACTCGACTTCTT"
            "TTGGAATATCTTGAAGAAAAATATGAAGAGCATTTGTATGAGCGCGATGAAGGTGATAAA"
            "TGGCGAAACAAAAAGTTTGAATTGGGTTTGGAGTTTCCCAATCTTCCTTATTATATTGAT"
            "GGTGATGTTAAATTAACACAGTCTATGGCCATCATACGTTATATAGCTGACAAGCACAAC"
            "ATGTTGGGTGGTTGTCCAAAAGAGCGTGCAGAGATTTCAATGCTTGAAGGAGCGGTTTTG"
            "GATATTAGATACGGTGTTTCGAGAATTGCATATAGTAAAGACTTTGAAACTCTCAAAGTT"
            "GATTTTCTTAGCAAGCTACCTGAAATGCTGAAAATGTTCGAAGATCGTTTATGTCATAAA"
            "ACATATTTAAATGGTGATCATGTAACCCATCCTGACTTCATGTTGTATGACGCTCTTGAT"
            "GTTGTTTTATACATGGACCCAATGTGCCTGGATGCGTTCCCAAAATTAGTTTGTTTTAAA"
            "AAACGTATTGAAGCTATCCCACAAATTGATAAGTACTTGAAATCCAGCAAGTATATAGCA"
            "TGGCCTTTGCAGGGCTGGCAAGCCACGTTTGGTGGTGGCGACCATCCTCCAAAA"
        ),
        notes=(
            "223 aa / ~26 kDa. Glutatyon-sefaroz ile saflaştırma. Büyük bir "
            "füzyon ortağı olduğu için çözünürlüğü belirgin artırır, ancak "
            "dimer oluşturur - oligomerik durum çalışmalarında yanıltır ve "
            "mutlaka kesilmesi gerekir. S. japonicum kaynaklı."
            "\n\n"
            "Dizi: Sj26 CDS'i, stop kodonu çıkarılmış (218 kodon). Katalogdaki "
            "223 aa / 669 bp değeri yanlıştı: gerçek Sj26 GST 218 aa."
        ),
        accession="M14654 (S. japonicum Sj26): 17..670",
        pmid="3047011",
    ),
    _p(
        part_id="TAG_MBP",
        name="MBP (maltoz bağlayan protein) füzyonu",
        category="tag",
        sequence=(
            "ATGAAAATCGAAGAAGGTAAACTGGTAATCTGGATTAACGGCGATAAAGGCTATAACGGT"
            "CTCGCTGAAGTCGGTAAGAAATTCGAGAAAGATACCGGAATTAAAGTCACCGTTGAGCAT"
            "CCGGATAAACTGGAAGAGAAATTCCCACAGGTTGCGGCAACTGGCGATGGCCCTGACATT"
            "ATCTTCTGGGCACACGACCGCTTTGGTGGCTACGCTCAATCTGGCCTGTTGGCTGAAATC"
            "ACCCCGGACAAAGCGTTCCAGGACAAGCTGTATCCGTTTACCTGGGATGCCGTACGTTAC"
            "AACGGCAAGCTGATTGCTTACCCGATCGCTGTTGAAGCGTTATCGCTGATTTATAACAAA"
            "GATCTGCTGCCGAACCCGCCAAAAACCTGGGAAGAGATCCCGGCGCTGGATAAAGAACTG"
            "AAAGCGAAAGGTAAGAGCGCGCTGATGTTCAACCTGCAAGAACCGTACTTCACCTGGCCG"
            "CTGATTGCTGCTGACGGGGGTTATGCGTTCAAGTATGAAAACGGCAAGTACGACATTAAA"
            "GACGTGGGCGTGGATAACGCTGGCGCGAAAGCGGGTCTGACCTTCCTGGTTGACCTGATT"
            "AAAAACAAACACATGAATGCAGACACCGATTACTCCATCGCAGAAGCTGCCTTTAATAAA"
            "GGCGAAACAGCGATGACCATCAACGGCCCGTGGGCATGGTCCAACATCGACACCAGCAAA"
            "GTGAATTATGGTGTAACGGTACTGCCGACCTTCAAGGGTCAACCATCCAAACCGTTCGTT"
            "GGCGTGCTGAGCGCAGGTATTAACGCCGCCAGTCCGAACAAAGAGCTGGCGAAAGAGTTC"
            "CTCGAAAACTATCTGCTGACTGATGAAGGTCTGGAAGCGGTTAATAAAGACAAACCGCTG"
            "GGTGCCGTAGCGCTGAAGTCTTACGAGGAAGAGTTGGCGAAAGATCCACGTATTGCCGCC"
            "ACCATGGAAAACGCCCAGAAAGGTGAAATCATGCCGAACATCCCGCAGATGTCCGCTTTC"
            "TGGTATGCCGTGCGTACTGCGGTGATCAACGCCGCCAGCGGTCGTCAGACTGTCGATGAA"
            "GCCCTGAAAGACGCGCAGACTCGTATCACCAAG"
        ),
        notes=(
            "~42 kDa. E. coli malE kaynaklı. Çözünürlük artırıcı olarak "
            "en etkili etiketlerden biri; amiloz reçinesi ile saflaştırılır. "
            "Periplazmik hedefleme için kendi sinyal peptidiyle de kullanılır."
            "\n\n"
            "Dizi: malE CDS'inden sinyal peptidi (kodon 1-26) çıkarılıp başa ATG "
            "eklendi - pMAL vektörlerindeki olgun MBP ile aynı: Met + "
            "MBP(27-396), 371 kodon."
        ),
        accession="NC_000913.3 malE: ATG + complement(4245232..4246341)",
        pmid="3049613",
    ),
    _p(
        part_id="TAG_SUMO",
        name="SUMO füzyonu (His-SUMO)",
        category="tag",
        sequence=(
            "ATGTCGGACTCAGAAGTCAATCAAGAAGCTAAGCCAGAGGTCAAGCCAGAAGTCAAGCCT"
            "GAGACTCACATCAATTTAAAGGTGTCCGATGGATCTTCAGAGATCTTCTTCAAGATCAAA"
            "AAGACCACTCCTTTAAGAAGGCTGATGGAAGCGTTCGCTAAAAGACAGGGTAAGGAAATG"
            "GACTCCTTAAGATTCTTGTACGACGGTATTAGAATTCAAGCTGATCAGACCCCTGAAGAT"
            "TTGGACATGGAGGATAACGATATTATTGAGGCTCACAGAGAACAGATTGGTGGT"
        ),
        notes=(
            "Çözünürlük artırır ve SUMO proteazı (Ulp1) füzyon bölgesini "
            "hedefe HİÇ artık bırakmadan keser - N-terminali doğal olması "
            "gereken proteinlerde en temiz seçenek."
            "\n\n"
            "Dizi: Smt3 CDS'inin 1-98. kodonları. Ulp1 Gly97-Gly98'den sonra "
            "kestiği için 99-101 (ATY) kuyruğu kasten dışarıda bırakıldı; böylece "
            "kesimden sonra hedef proteinde artık kalmaz."
        ),
        accession="NM_001180818 (S. cerevisiae SMT3): 1..294",
        pmid="15299233",
    ),
    _p(
        part_id="TAG_TRXA",
        name="Tiyoredoksin (TrxA) füzyonu",
        category="tag",
        sequence=(
            "ATGAGCGATAAAATTATTCACCTGACTGACGACAGTTTTGACACGGATGTACTCAAAGCG"
            "GACGGGGCGATCCTCGTCGATTTCTGGGCAGAGTGGTGCGGTCCGTGCAAAATGATCGCC"
            "CCGATTCTGGATGAAATCGCTGACGAATATCAGGGCAAACTGACCGTTGCAAAACTGAAC"
            "ATCGATCAAAACCCTGGCACTGCGCCGAAATATGGCATCCGTGGTATCCCGACTCTGCTG"
            "CTGTTCAAAAACGGTGAAGTGGCGGCAACCAAAGTGGGTGCACTGTCTAAAGGTCAGTTG"
            "AAAGAGTTCCTCGACGCTAACCTGGCG"
        ),
        notes=(
            "109 aa. Çözünürlüğü artırır ve sitoplazmada disülfit "
            "oluşumuna yardım eder. Isıya dayanıklı olduğu için ön "
            "saflaştırmada ısı adımı kullanılabilir."
            "\n\n"
            "Dizi: trxA CDS'i, stop kodonu çıkarılmış (109 kodon)."
        ),
        accession="NC_000913.3 trxA: 3965761..3966087",
        pmid="8365758",
    ),
]

# ---------------------------------------------------------------------------
# PROTEAZ KESİM BÖLGELERİ
# ---------------------------------------------------------------------------

PROTEASE_SITES = [
    _p(
        part_id="PRS_TEV",
        name="TEV proteaz bölgesi (ENLYFQ/G)",
        category="protease_site",
        sequence="GAAAACCTGTATTTTCAGGGC",
        notes=(
            "Çok özgül; kesimden sonra hedefin N-terminalinde tek bir "
            "glisin kalır. 4 C'de de çalışır. Kendisi His etiketli olarak "
            "üretilip Ni kolonunda ayrılabilir."
        ),
        accession="",
        pmid="2site",
    ),
    _p(
        part_id="PRS_THROMBIN",
        name="Trombin bölgesi (LVPR/GS)",
        category="protease_site",
        sequence="CTGGTGCCGCGCGGCAGC",
        notes=(
            "pET-28a'nın standart kesim bölgesi. Kesimden sonra GS artığı "
            "kalır. Özgüllüğü TEV'den düşük - hedef proteinde içsel "
            "benzer diziler varsa istenmeyen kesim olabilir."
        ),
        accession="pET-28a bölgesi",
        pmid="",
    ),
    _p(
        part_id="PRS_ENTEROKINASE",
        name="Enterokinaz bölgesi (DDDDK/)",
        category="protease_site",
        sequence="GACGACGACGACAAG",
        notes=(
            "Artık bırakmadan keser. FLAG etiketinin sonunda doğal olarak "
            "bulunur. Özgüllüğü sınırlı, yan kesim bildirilmiştir."
        ),
        accession="",
        pmid="",
    ),
    _p(
        part_id="PRS_PRESCISSION",
        name="PreScission / HRV 3C bölgesi (LEVLFQ/GP)",
        category="protease_site",
        sequence="CTGGAAGTTCTGTTCCAGGGGCCC",
        notes=(
            "4 C'de verimli çalışır - sıcaklığa hassas hedeflerde avantaj. "
            "GST etiketli proteaz formu kolonda kesim yapıp kendisi "
            "reçinede kalır."
        ),
        accession="",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# BAĞLAYICILAR
# ---------------------------------------------------------------------------

LINKERS = [
    _p(
        part_id="LNK_GS_SHORT",
        name="Esnek bağlayıcı (GGGGS)",
        category="linker",
        sequence="GGTGGCGGTGGCTCT",
        notes="Kısa esnek bağlayıcı. Füzyon ortaklarının bağımsız katlanması için.",
    ),
    _p(
        part_id="LNK_GS_LONG",
        name="Esnek bağlayıcı (GGGGS x3)",
        category="linker",
        sequence="GGTGGCGGTGGCTCTGGTGGCGGTGGCTCTGGTGGCGGTGGCTCT",
        notes=(
            "Uzun esnek bağlayıcı. Tekrarlı dizi olduğu için PCR ve "
            "dizileme sırasında kayma riski taşır; tekrarları farklı "
            "kodonlarla yazmak daha güvenli."
        ),
    ),
    _p(
        part_id="LNK_RIGID",
        name="Rijit bağlayıcı (EAAAK x2)",
        category="linker",
        sequence="GAAGCTGCTGCTAAAGAAGCTGCTGCTAAA",
        notes=(
            "Alfa-heliks yapıda, iki domaini sabit mesafede tutar. "
            "Domainler arası istenmeyen etkileşimi azaltır."
        ),
    ),
]

# ---------------------------------------------------------------------------
# TERMİNATÖRLER
# ---------------------------------------------------------------------------

TERMINATORS = [
    _p(
        part_id="TER_T7",
        name="T7 terminatör",
        category="terminator",
        sequence="CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        notes=(
            "T7 RNA polimeraza özgü. pET vektörlerinde ifade kasetini "
            "kapatır. Konak RNA polimerazını durdurmakta zayıftır - "
            "sabit promotörlü kasetlerde tek başına yetmez."
        ),
        accession="T7 fag Tphi terminatörü",
        pmid="",
    ),
    _p(
        part_id="TER_B0015",
        name="B0015 çift terminatör",
        category="terminator",
        sequence=(
            "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGT"
            "CGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATA"
        ),
        notes=(
            "B0010 (rrnB T1) + B0012 (TE) birleşimi. Konak RNA polimerazını "
            "verimli durdurur. iGEM standardı; sentetik devrelerde "
            "kasetler arası okuma-geçişini engellemek için kullanılır."
        ),
        accession="iGEM Registry: BBa_B0015",
        pmid="",
    ),
    _p(
        part_id="TER_RRNB",
        name="rrnB T1T2 terminatörü",
        category="terminator",
        sequence=(
            "CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCG"
            "GTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAA"
            "CGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAG"
            "AAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACA"
        ),
        notes=(
            "E. coli rrnB operonunun çift terminatörü. Çok güçlü; pET ve "
            "pTrc vektörlerinde kaset sınırı olarak kullanılır. İki ayrı "
            "terminatör ve aralarındaki bölgeyi içerdiği için uzun."
            "\n\n"
            "Dizi: rrfB (5S rRNA) geninin bitiminden hemen sonrası; T1 ve T2 saç "
            "tokaları ile U-izleri dizide birebir bulundu, bitiş T2'den 10 bp "
            "sonra. Katalogdaki 360 bp tahmini, gerçek T1+ara+T2 bölgesi 218 bp."
        ),
        accession="NC_000913.3 (rrnB, rrfB sonrası): 4171757..4171974",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# REPLİKASYON ORİJİNLERİ
# ---------------------------------------------------------------------------

ORIGINS = [
    _p(
        part_id="ORI_PUC",
        name="pUC ori (pMB1 türevi, yüksek kopya)",
        category="ori",
        sequence=(
            "CTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGC"
            "TACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTG"
            "GCTTCAGCAGAGCGCAGATACCAAATACTGTTCTTCTAGTGTAGCCGTAGTTAGGCCACC"
            "ACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGG"
            "CTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGG"
            "ATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAA"
            "CGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCG"
            "AAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGA"
            "GGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCT"
            "GACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCA"
            "GCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTT"
        ),
        copy_number="500-700",
        inc_group="ColE1 / IncB",
        incompatible_with=["ORI_PBR322", "ORI_COLE1"],
        notes=(
            "pMB1 orijininin RNA II bölgesinde tek nokta mutasyonu taşır; "
            "bu mutasyon kopya sayısı kontrolünü bozarak 37 C'de kopya "
            "sayısını çok yükseltir. Yüksek plazmit verimi (mini-prep) "
            "için ideal, ancak yüksek kopya + güçlü promotör birleşimi "
            "metabolik yük yaratır. 30 C'de kopya sayısı düşer."
            "\n\n"
            "Dizi: ColE1 (J01566) orijin bölgesi pBR322'ye %93.8 kimlikle "
            "hizalanıp ColE1'in ori noktası pBR322'nin kendi rep_origin "
            "anotasyonuna 3 bp içinde oturdu; aynı bölge pUC19'da "
            "complement(805..1459). Sınırlar RNA II +1'den 40 bp önce, ori "
            "noktasından 60 bp sonra. ORI_PBR322 ile aradaki TEK fark 152. "
            "bazdaki C->T: pUC'nin yüksek kopya mutasyonu, notlarda anlatılan "
            "mutasyonun ta kendisi."
        ),
        accession="L09137 (pUC19): complement(805..1459)",
        pmid="2985470",
    ),
    _p(
        part_id="ORI_PBR322",
        name="pBR322 ori (pMB1, düşük kopya)",
        category="ori",
        sequence=(
            "CTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGC"
            "TACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTG"
            "GCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACC"
            "ACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGG"
            "CTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGG"
            "ATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAA"
            "CGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCG"
            "AAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGA"
            "GGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCT"
            "GACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCA"
            "GCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTT"
        ),
        copy_number="15-20",
        inc_group="ColE1 / IncB",
        incompatible_with=["ORI_PUC", "ORI_COLE1"],
        notes=(
            "Rop/Rom proteini ile kopya sayısı baskılanır. Toksik ürünler "
            "ve büyük insertler için pUC orijininden daha kararlıdır."
            "\n\n"
            "Dizi: ORI_PUC ile aynı geometri (RNA II +1 -40 .. ori +60). İki dizi "
            "655 bp boyunca yalnızca 152. bazda ayrılıyor (burada C, pUC19'da T)."
        ),
        accession="J01749 (pBR322): complement(2472..3126)",
        pmid="",
    ),
    _p(
        part_id="ORI_COLE1",
        name="ColE1 ori",
        category="ori",
        sequence=(
            "CTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGC"
            "TACCAACGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTG"
            "GCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTCGGGCCACT"
            "ACTTCAAGAACTCTGTAGCACCGTTTGTGCCATCATCGCTCTGCTAATCCGGTTACCAGT"
            "GGCTGCTGCCAGTGGCGTTAAGGCGTGCCTTACCGGGTTGGACTCAAGACGATAGTTACC"
            "GGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCG"
            "AACGACCTACACCGAACTGAGATACCAACAGCGTGAGCTATGAGAAAGCGCCACGCTTCC"
            "CGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCAC"
            "GAGGGAGCTTCCAGGGGGAAACGCCTGGTAGCTTTATAGTCCTGTCGGGTTTCGCCACCT"
            "CTGACTTGAGCGTCTATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGC"
            "CTGCTACGTGGCCTTCTTCCTGTTCCTGGTCTTTTGCTCACATGTTCTTTCCGGCCT"
        ),
        copy_number="15-20",
        inc_group="ColE1",
        incompatible_with=["ORI_PUC", "ORI_PBR322"],
        notes=(
            "Klasik gevşek kontrollü orijin. pET vektörleri bu aileyi "
            "kullanır - bu yüzden pET + pUC tabanlı bir plazmit aynı "
            "hücrede kararlı tutulamaz."
            "\n\n"
            "Dizi: kaydın kendi anotasyonları (RNA II başlangıcı 641, RNA I "
            "645..751 ters zincir, replikasyon orijini 1197); bölge RNA II +1'den "
            "40 bp önce başlar, ori'den 60 bp sonra biter."
        ),
        accession="J01566 (ColE1): 601..1257",
        pmid="",
    ),
    _p(
        part_id="ORI_P15A",
        name="p15A ori",
        category="ori",
        sequence=(
            "GCGCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGAAGTG"
            "CTTCATGTGGCAGGAGAAAAAAGGCTGCACCGGTGCGTCAGCAGAATATGTGATACAGGA"
            "TATATTCCGCTTCCTCGCTCACTGACTCGCTACGCTCGGTCGTTCGACTGCGGCGAGCGG"
            "AAATGGCTTACGAACGGGGCGGAGATTTCCTGGAAGATGCCAGGAAGATACTTAACAGGG"
            "AAGTGAGAGGGCCGCGGCAAAGCCGTTTTTCCATAGGCTCCGCCCCCCTGACAAGCATCA"
            "CGAAATCTGACGCTCAAATCAGTGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGC"
            "GTTTCCCCCTGGCGGCTCCCTCGTGCGCTCTCCTGTTCCTGCCTTTCGGTTTACCGGTGT"
            "CATTCCGCTGTTATGGCCGCGTTTGTCTCATTCCACGCCTGACACTCAGTTCCGGGTAGG"
            "CAGTTCGCTCCAAGCTGGACTGTATGCACGAACCCCCCGTTCAGTCCGACCGCTGCGCCT"
            "TATCCGGTAACTATCGTCTTGAGTCCAACCCGGAAAGACATGCAAAAGCACCACTGGCAG"
            "CAGCCACTGGTAATTGATTTAGAGGAGTTAGTCTTGAAGTCATGCGCCGGTTAAGGCTAA"
            "ACTGAAAGGACAAGTTTTGGTGACTGCGCTCCTCCAAGCCAGTTACCTCGGTTCAAAGAG"
            "TTGGTAGCTCAGAGAACCTTCGAAAAACCGCCCTGCAAGGCGGTTTTTTCGTTTTCAGAG"
            "CAAGAGATTACGCGCAGACCAAAACGATCTCAAGAAGATCATCTTATTAATCAGATAAAA"
            "TATTTCTAGATTTCAGTGCAATTTATCTCTTCAAATGTAGCACCTGAAGTCAGCCCCATA"
            "CGATATAAGTTGT"
        ),
        copy_number="10-12",
        inc_group="IncP15A",
        notes=(
            "ColE1 ailesiyle UYUMLU - iki plazmitin aynı hücrede "
            "tutulması gereken ko-ekspresyon ve devre kurulumlarında "
            "standart seçimdir. pACYC184/pACYCDuet bu orijini taşır."
            "\n\n"
            "Dizi: pACYC184 kaydının kendi 'p15A origin of replication' "
            "rep_origin anotasyonu. Katalogdaki 570 bp tahmini, anotasyon 913 bp."
        ),
        accession="X06403 (pACYC184): 581..1493",
        pmid="",
    ),
    _p(
        part_id="ORI_PSC101",
        name="pSC101 ori",
        category="ori",
        sequence=(
            "AGCTTGCGAGGGTGCTACTTAAGCCTTTAGGGTTTTAAGGTCTGTTTTGTAGAGGAGCAA"
            "ACAGCGTTTGCGACATCCTTTTGTAATACTGCGGAACTGACTAAAGTAGTGAGTTATACA"
            "CAGGGCTGGGATCTATTCTTTTTATCTTTTTTTATTCTTTCTTTATTCTATAAATTATAA"
            "CCACTTGAATATAAACAAAAAAAACACACAAAGGTCTAGCGGAATTTACAGAGGGTCTAG"
            "CAGAATTTACAAGTTTTCCAGCAAAGGTCTAGCAGAATTTACAGATACCCACAACTCAAA"
            "GGAAAAGGACTAGTAATTATCATTGACTAGCCCATCTCAATTGGTATAGTGATTAAAATC"
            "ACCTAGACCAATTGAGATGTATGTCTGAATTAGTTGTTTTCAAAGCAAATGAACTAGCGA"
            "TTAGTCGCTATGACTTAACGGAGCATGAAACCAAGCTAATTTTATGCTGTGTGGCACTAC"
            "TCAACCCCACGATTGAAAACCCTACAAGGAAAGAACGGACGGTATCGTTCACTTATAACC"
            "AATACGCTCAGATGATGAACATCAGTAGGGAAAATGCTTATGGTGTATTAGCTAAAGCAA"
            "CCAGAGAGCTGATGACGAGAACTGTGGAAATCAGGAATCCTTTGGTTAAAGGCTTTGAGA"
            "TTTTCCAGTGGACAAACTATGCCAAGTTCTCAAGCGAAAAATTAGAATTAGTTTTTAGTG"
            "AAGAGATATTGCCTTATCTTTTCCAGTTAAAAAAATTCATAAAATATAATCTGGAACATG"
            "TTAAGTCTTTTGAAAACAAATACTCTATGAGGATTTATGAGTGGTTATTAAAAGAACTAA"
            "CACAAAAGAAAACTCACAAGGCAAATATAGAGATTAGCCTTGATGAATTTAAGTTCATGT"
            "TAATGCTTGAAAATAACTACCATGAGTTTAAAAGGCTTAACCAATGGGTTTTGAAACCAA"
            "TAAGTAAAGATTTAAACACTTACAGCAATATGAAATTGGTGGTTGATAAGCGAGGCCGCC"
            "CGACTGATACGTTGATTTTCCAAGTTGAACTAGATAGACAAATGGATCTCGTAACCGAAC"
            "TTGAGAACAACCAGATAAAAATGAATGGTGACAAAATACCAACAACCATTACATCAGATT"
            "CCTACCTACATAACGGACTAAGAAAAACACTACACGATGCTTTAACTGCAAAAATTCAGC"
            "TCACCAGTTTTGAGGCAAAATTTTTGAGTGACATGCAAAGTAAGTATGATCTCAATGGTT"
            "CGTTCTCATGGCTCACGCAAAAACAACGAACCACACTAGAGAACATACTGGCTAAATACG"
            "GAAGGATCTGA"
        ),
        copy_number="~5",
        inc_group="IncFII",
        requires=["RepA proteini (genelde orijin bölgesiyle birlikte gelir)"],
        notes=(
            "Çok düşük kopya, sıkı kontrollü. Toksik gen taşımak ve "
            "ColE1 + p15A ile üçüncü bir plazmit katmanı eklemek için "
            "kullanılır. Düşük kopya olduğu için plazmit verimi düşüktür."
            "\n\n"
            "Dizi: K00042 kaydının 'Rep proteini varlığında replikasyon için "
            "zorunlu bölge' (370..826), üç iteron, uyumsuzluk belirleyicisi ve "
            "repA geni (750..1700) birlikte. pSC101 orijini RepA olmadan "
            "çalışmadığı için parça ori + repA'yı BİRLİKTE içerir - bu yüzden "
            "1331 bp."
        ),
        accession="K00042 (pSC101 ori + repA): 370..1700",
        pmid="",
    ),
    _p(
        part_id="ORI_R6K",
        name="R6K gamma ori (koşullu)",
        category="ori",
        sequence=(
            "GATCGCTAGTTTGTTTTGACTCCATCCATTAGGGCTTCTAAAACGCCTTCTAAGGCCATG"
            "TCAGCCGTTAAGTGTTCCTGTGTCACTGAAAATTGCTTTGAGAGGCTCTAAGGGCTTCTC"
            "AGTGCGTTACATCCCTGGCTTGTTGTCCACAACCGTTAAACCTTAAAAGCTTTAAAAGCC"
            "TTATATATTCTTTTTTTTCTTATAAAACTTAAAACCTTAGAGGCTATTTAAGTTGCTGAT"
            "TTATATTAATTTTATTGTTCAAACATGAGAGCTTAGTACGTGAAACATGAGAGCTTAGTA"
            "CGTTAGCCATGAGAGCTTAGTACGTTAGCCATGAGGGTTTAGTTCGTTAAACATGAGAGC"
            "TTAGTACGTTAAACATGAGAGCTTAGTACGTGAAACATGAGAGCTTAGTACGTACTATCA"
            "ACAGGTTGAACTGCTGATCTTC"
        ),
        copy_number="15 (sadece pir+ konakta)",
        requires=["pi proteini (pir geni) - konak genomunda, örn. EC100D pir+"],
        notes=(
            "Normal suşlarda ÇOĞALAMAZ. Suş inşasında 'suicide vector' "
            "olarak kullanılır: pir+ suşta çoğaltılır, hedef suşa "
            "aktarıldığında replike olamaz ve sadece genoma entegre "
            "olanlar seçilebilir."
            "\n\n"
            "Dizi: M65025 kaydının kendi 'replication origin gamma region' "
            "anotasyonu. pi (p1) proteini aynı kayıtta 623..1540'ta ama kasten "
            "dışarıda: gamma ori'nin koşulluluğu pir'in konaktan gelmesine "
            "dayanır."
        ),
        accession="M65025 (R6K): 1..442",
        pmid="",
    ),
    _p(
        part_id="ORI_F1",
        name="f1 / M13 ori (tek zincir)",
        category="ori",
        sequence=(
            "CGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAG"
            "CGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTT"
            "TCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCA"
            "CCTCGACCCCAAAAAACTTGATTTGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATA"
            "GACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCA"
            "AACTGGAACAACACTCAACCCTATCTCGGGCTATTCTTTTGATTTATAAGGGATTTTGCC"
            "GATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAA"
            "CAAAATATTAACGTTTACAATTTAAATATTTGCTTATACAATCTTCCTGTTTTTGGGGCT"
            "TTTCTGATTATCAACCGGGGTACAT"
        ),
        notes=(
            "Çift zincir replikasyonu sağlamaz; yardımcı fag "
            "enfeksiyonunda tek zincirli DNA üretimi için. Eski Sanger "
            "dizileme ve bölge-yönelimli mutagenez protokolleri ile "
            "fag gösterimi (phage display) kütüphanelerinde kullanılır."
            "\n\n"
            "Dizi: M13 gen IV ile gen II arasındaki genler arası bölgenin tamamı; "
            "(+) ve (-) zincir orijinlerinin ikisini de içerir."
        ),
        accession="V00604 (faj M13): 5501..6005",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# SEÇİLİM MARKERLARI
# ---------------------------------------------------------------------------

MARKERS = [
    _p(
        part_id="MRK_BLA",
        name="bla (AmpR, beta-laktamaz)",
        category="marker",
        sequence=(
            "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCT"
            "GTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCA"
            "CGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCC"
            "GAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCC"
            "CGTGTTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTG"
            "GTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTA"
            "TGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATC"
            "GGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTT"
            "GATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATG"
            "CCTGCAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCT"
            "TCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGC"
            "TCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCT"
            "CGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTAC"
            "ACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCC"
            "TCACTGATTAAGCATTGGTAA"
        ),
        notes=(
            "Ampisilin 100 ug/mL ya da karbenisilin 50 ug/mL. "
            "DİKKAT: enzim periplazmaya salgılanır ve besiyerindeki "
            "ampisilini parçalar; uzun inkübasyonda plazmitsiz "
            "'satellit koloniler' çıkar ve yoğun kültürlerde seçilim "
            "baskısı kaybolur. Büyük ölçekli ifadede kanamisin tercih "
            "edilir. Karbenisilin daha kararlıdır."
            "\n\n"
            "Dizi: pBR322 kaydının kendi bla CDS anotasyonu; çeviri TEM-1 "
            "beta-laktamazı (286 aa)."
        ),
        accession="J01749 (pBR322) bla: complement(3293..4153)",
        pmid="",
    ),
    _p(
        part_id="MRK_KAN",
        name="kanR / nptII (KanR, aminoglikozid fosfotransferaz)",
        category="marker",
        sequence=(
            "ATGATTGAACAAGATGGATTGCACGCAGGTTCTCCGGCCGCTTGGGTGGAGAGGCTATTC"
            "GGCTATGACTGGGCACAACAGACAATCGGCTGCTCTGATGCCGCCGTGTTCCGGCTGTCA"
            "GCGCAGGGGCGCCCGGTTCTTTTTGTCAAGACCGACCTGTCCGGTGCCCTGAATGAACTG"
            "CAGGACGAGGCAGCGCGGCTATCGTGGCTGGCCACGACGGGCGTTCCTTGCGCAGCTGTG"
            "CTCGACGTTGTCACTGAAGCGGGAAGGGACTGGCTGCTATTGGGCGAAGTGCCGGGGCAG"
            "GATCTCCTGTCATCTCACCTTGCTCCTGCCGAGAAAGTATCCATCATGGCTGATGCAATG"
            "CGGCGGCTGCATACGCTTGATCCGGCTACCTGCCCATTCGACCACCAAGCGAAACATCGC"
            "ATCGAGCGAGCACGTACTCGGATGGAAGCCGGTCTTGTCGATCAGGATGATCTGGACGAA"
            "GAGCATCAGGGGCTCGCGCCAGCCGAACTGTTCGCCAGGCTCAAGGCGCGCATGCCCGAC"
            "GGCGAGGATCTCGTCGTGACCCATGGCGATGCCTGCTTGCCGAATATCATGGTGGAAAAT"
            "GGCCGCTTTTCTGGATTCATCGACTGTGGCCGGCTGGGTGTGGCGGACCGCTATCAGGAC"
            "ATAGCGTTGGCTACCCGTGATATTGCTGAAGAGCTTGGCGGCGAATGGGCTGACCGCTTC"
            "CTCGTGCTTTACGGTATCGCCGCTCCCGATTCGCAGCGCATCGCCTTCTATCGCCTTCTT"
            "GACGAGTTCTTCTGA"
        ),
        notes=(
            "Kanamisin 50 ug/mL. Enzim hücre içinde kalır, antibiyotik "
            "besiyerinde parçalanmaz - seçilim baskısı kültür boyunca "
            "korunur. Satellit koloni yapmaz. pET-28a ve çoğu modern "
            "ifade vektörünün markerı."
            "\n\n"
            "Dizi: Tn5 kaydının kendi nptII/neo CDS anotasyonu (264 aa + stop = "
            "795 bp). Katalogdaki 816 bp, vektörlerde kullanılan daha uzun "
            "alternatif başlangıca ait; burada kaydın anotasyonu esas alındı."
        ),
        accession="U00004 (Tn5) nptII: 1551..2345",
        pmid="",
    ),
    _p(
        part_id="MRK_CAT",
        name="cat (CmR, kloramfenikol asetiltransferaz)",
        category="marker",
        sequence=(
            "ATGGAGAAAAAAATCACTGGATATACCACCGTTGATATATCCCAATGGCATCGTAAAGAA"
            "CATTTTGAGGCATTTCAGTCAGTTGCTCAATGTACCTATAACCAGACCGTTCAGCTGGAT"
            "ATTACGGCCTTTTTAAAGACCGTAAAGAAAAATAAGCACAAGTTTTATCCGGCCTTTATT"
            "CACATTCTTGCCCGCCTGATGAATGCTCATCCGGAATTCCGTATGGCAATGAAAGACGGT"
            "GAGCTGGTGATATGGGATAGTGTTCACCCTTGTTACACCGTTTTCCATGAGCAAACTGAA"
            "ACGTTTTCATCGCTCTGGAGTGAATACCACGACGATTTCCGGCAGTTTCTACACATATAT"
            "TCGCAAGATGTGGCGTGTTACGGTGAAAACCTGGCCTATTTCCCTAAAGGGTTTATTGAG"
            "AATATGTTTTTCGTCTCAGCCAATCCCTGGGTGAGTTTCACCAGTTTTGATTTAAACGTG"
            "GCCAATATGGACAACTTCTTCGCCCCCGTTTTCACCATGGGCAAATATTATACGCAAGGC"
            "GACAAGGTGCTGATGCCGCTGGCGATTCAGGTTCATCATGCCGTCTGTGATGGCTTCCAT"
            "GTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGCGGGGCGTAA"
        ),
        notes=(
            "Kloramfenikol 25-34 ug/mL. Bakteriyostatik - koloniler yavaş "
            "büyür. p15A orijinli ko-ekspresyon plazmitlerinin (pACYCDuet, "
            "pLysS, Rosetta tRNA plazmiti) standart markerı."
            "\n\n"
            "Dizi: pACYC184'te cat CDS'i kaydın numaralandırma başlangıcını "
            "aştığı için kayıtta CDS olarak işaretli değil. İki zincirde 660 "
            "nt'lik tek açık okuma çerçevesi arandı; çeviri UniProt P62577 (Tn9 "
            "CAT, 219 aa) ile birebir aynı."
        ),
        accession="X06403 (pACYC184) cat: complement(join(3805..4245,1..219))",
        pmid="",
    ),
    _p(
        part_id="MRK_TET",
        name="tetA (TetR, tetrasiklin efluks pompası)",
        category="marker",
        sequence=(
            "ATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGC"
            "ATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGC"
            "ATCGCCAGTCACTATGGCGTGCTGCTAGCGCTATATGCGTTGATGCAATTTCTATGCGCA"
            "CCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTA"
            "CTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTAC"
            "GCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATC"
            "GCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTC"
            "GGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCTTGCAT"
            "GCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTA"
            "ATGCAGGAGTCGCATAAGGGAGAGCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTC"
            "AGCTCCTTCCGGTGGGCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTT"
            "ATCATGCAACTCGTAGGACAGGTGCCGGCAGCGCTCTGGGTCATTTTCGGCGAGGACCGC"
            "TTTCGCTGGAGCGCGACGATGATCGGCCTGTCGCTTGCGGTATTCGGAATCTTGCACGCC"
            "CTCGCTCAAGCCTTCGTCACTGGTCCCGCCACCAAACGTTTCGGCGAGAAGCAGGCCATT"
            "ATCGCCGGCATGGCGGCCGACGCGCTGGGCTACGTCTTGCTGGCGTTCGCGACGCGAGGC"
            "TGGATGGCCTTCCCCATTATGATTCTTCTCGCTTCCGGCGGCATCGGGATGCCCGCGTTG"
            "CAGGCCATGCTGTCCAGGCAGGTAGATGACGACCATCAGGGACAGCTTCAAGGATCGCTC"
            "GCGGCTCTTACCAGCCTAACTTCGATCACTGGACCGCTGATCGTCACGGCGATTTATGCC"
            "GCCTCGGCGAGCACATGGAACGGGTTGGCATGGATTGTAGGCGCCGCCCTATACCTTGTC"
            "TGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCGGGCCACCTCGACCTGA"
        ),
        notes=(
            "Tetrasiklin 10-15 ug/mL. Bir membran pompasıdır, enzim değil - "
            "aşırı ifadesi membran stresi yaratır. TetR tabanlı indükleyici "
            "sistemle (P-tet) aynı plazmitte kullanılmamalı: aTc hem "
            "seçilimi hem devreyi etkiler."
            "\n\n"
            "Dizi: pBR322 kaydının kendi tet (tetA) CDS anotasyonu, 396 aa."
        ),
        accession="J01749 (pBR322) tetA: 86..1276",
        pmid="",
    ),
    _p(
        part_id="MRK_AADA",
        name="aadA (SpecR/StrepR, aminoglikozid adeniltransferaz)",
        category="marker",
        sequence=(
            "ATGAGGGAAGCGGTGATCGCCGAAGTATCGACTCAACTATCAGAGGTAGTTGGCGTCATC"
            "GAGCGCCATCTCGAACCGACGTTGCTGGCCGTACATTTGTACGGCTCCGCAGTGGATGGC"
            "GGCCTGAAGCCACACAGTGATATTGATTTGCTGGTTACGGTGACCGTAAGGCTTGATGAA"
            "ACAACGCGGCGAGCTTTGATCAACGACCTTTTGGAAACTTCGGCTTCCCCTGGAGAGAGC"
            "GAGATTCTCCGCGCTGTAGAAGTCACCATTGTTGTGCACGACGACATCATTCCGTGGCGT"
            "TATCCAGCTAAGCGCGAACTGCAATTTGGAGAATGGCAGCGCAATGACATTCTTGCAGGT"
            "ATCTTCGAGCCAGCCACGATCGACATTGATCTGGCTATCTTGCTGACAAAAGCAAGAGAA"
            "CATAGCGTTGCCTTGGTAGGTCCAGCGGCGGAGGAACTCTTTGATCCGGTTCCTGAACAG"
            "GATCTATTTGAGGCGCTAAATGAAACCTTAACGCTATGGAACTCGCCGCCCGACTGGGCT"
            "GGCGATGAGCGAAATGTAGTGCTTACGTTGTCCCGCATTTGGTACAGCGCAGTAACCGGC"
            "AAAATCGCGCCGAAGGATGTCGCTGCCGACTGGGCAATGGAGCGCCTGCCGGCCCAGTAT"
            "CAGCCCGTCATACTTGAAGCTAGACAGGCTTATCTTGGACAAGAAGAAGATCGCTTGGCC"
            "TCGCGCGCAGATCAGTTGGAAGAATTTGTCCACTACGTGAAAGGCGAGATCACCAAGGTA"
            "GTCGGCAAATAA"
        ),
        notes=(
            "Spektinomisin 50-100 ug/mL. Dördüncü bağımsız seçilim "
            "katmanı gerektiğinde (Amp + Kan + Cm zaten kullanılıyorsa) "
            "devreye girer."
            "\n\n"
            "Dizi: sınıf 1 integron tipi aadA geni, 263 aa + stop."
        ),
        accession="U57972 aadA: 1..792",
        pmid="",
    ),
    _p(
        part_id="MRK_CCDB",
        name="ccdB (negatif seçilim)",
        category="marker",
        sequence=(
            "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTA"
            "CAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGT"
            "CTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGC"
            "TGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTG"
            "GCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGA"
            "ATATAA"
        ),
        requires=["Normal konak (ccdB'ye duyarlı); çoğaltma için DB3.1 gibi gyrA-r suş"],
        notes=(
            "DNA giraz zehiri - ifade edilirse hücreyi öldürür. Klonlama "
            "bölgesine yerleştirilir: insert girince ccdB kaybolur, "
            "böylece sadece doğru rekombinantlar büyür. Boş vektör "
            "arka planını neredeyse sıfıra indirir. Vektörün KENDİSİNİ "
            "çoğaltmak için dirençli özel suş gerekir."
            "\n\n"
            "Dizi: E. coli K-12 F plazmiti kaydının kendi ccdB(letD) CDS "
            "anotasyonu, 101 aa."
        ),
        accession="AP001918 (F plazmiti) ccdB/letD: 46560..46865",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# DÜZENLEYİCİLER
# ---------------------------------------------------------------------------

REGULATORS = [
    _p(
        part_id="REG_LACI",
        name="lacI (LacI baskılayıcısı)",
        category="regulator",
        sequence=(
            "GTGAAACCAGTAACGTTATACGATGTCGCAGAGTATGCCGGTGTCTCTTATCAGACCGTT"
            "TCCCGCGTGGTGAACCAGGCCAGCCACGTTTCTGCGAAAACGCGGGAAAAAGTGGAAGCG"
            "GCGATGGCGGAGCTGAATTACATTCCCAACCGCGTGGCACAACAACTGGCGGGCAAACAG"
            "TCGTTGCTGATTGGCGTTGCCACCTCCAGTCTGGCCCTGCACGCGCCGTCGCAAATTGTC"
            "GCGGCGATTAAATCTCGCGCCGATCAACTGGGTGCCAGCGTGGTGGTGTCGATGGTAGAA"
            "CGAAGCGGCGTCGAAGCCTGTAAAGCGGCGGTGCACAATCTTCTCGCGCAACGCGTCAGT"
            "GGGCTGATCATTAACTATCCGCTGGATGACCAGGATGCCATTGCTGTGGAAGCTGCCTGC"
            "ACTAATGTTCCGGCGTTATTTCTTGATGTCTCTGACCAGACACCCATCAACAGTATTATT"
            "TTCTCCCATGAAGACGGTACGCGACTGGGCGTGGAGCATCTGGTCGCATTGGGTCACCAG"
            "CAAATCGCGCTGTTAGCGGGCCCATTAAGTTCTGTCTCGGCGCGTCTGCGTCTGGCTGGC"
            "TGGCATAAATATCTCACTCGCAATCAAATTCAGCCGATAGCGGAACGGGAAGGCGACTGG"
            "AGTGCCATGTCCGGTTTTCAACAAACCATGCAAATGCTGAATGAGGGCATCGTTCCCACT"
            "GCGATGCTGGTTGCCAACGATCAGATGGCGCTGGGCGCAATGCGCGCCATTACCGAGTCC"
            "GGGCTGCGCGTTGGTGCGGATATCTCGGTAGTGGGATACGACGATACCGAAGACAGCTCA"
            "TGTTATATCCCGCCGTTAACCACCATCAAACAGGATTTTCGCCTGCTGGGGCAAACCAGC"
            "GTGGACCGCTTGCTGCAACTCTCTCAGGGCCAGGCGGTGAAGGGCAATCAGCTGTTGCCC"
            "GTCTCACTGGTGAAAAGAAAAACCACCCTGGCGCCCAATACGCAAACCGCCTCTCCCCGC"
            "GCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAG"
            "TGA"
        ),
        notes=(
            "lac/tac/trc/T7lac promotörlerini baskılar; IPTG ile serbest "
            "kalır. Konak lacI-q allelini taşımıyorsa plazmit üzerinde "
            "bulunması gerekir, yoksa indüksiyon öncesi ciddi kaçak "
            "ifade olur. lacI-q promotör mutasyonu normalden ~10 kat "
            "fazla baskılayıcı üretir."
            "\n\n"
            "Dizi: NP_414879'un coded_by alanı. Başlama kodonu GTG'dir (LacI "
            "proteini Val ile başlar) - bu bir hata değil."
        ),
        accession="NC_000913.3 lacI: complement(366428..367510)",
        pmid="",
    ),
    _p(
        part_id="REG_ARAC",
        name="araC (AraC düzenleyicisi)",
        category="regulator",
        sequence=(
            "ATGGCTGAAGCGCAAAATGATCCCCTGCTGCCGGGATACTCGTTTAACGCCCATCTGGTG"
            "GCGGGTTTAACGCCGATTGAGGCCAACGGTTATCTCGATTTTTTTATCGACCGACCGCTG"
            "GGAATGAAAGGTTATATTCTCAATCTCACCATTCGCGGTCAGGGGGTGGTGAAAAATCAG"
            "GGACGAGAATTTGTCTGCCGACCGGGTGATATTTTGCTGTTCCCGCCAGGAGAGATTCAT"
            "CACTACGGTCGTCATCCGGAGGCTCGCGAATGGTATCACCAGTGGGTTTACTTTCGTCCG"
            "CGCGCCTACTGGCATGAATGGCTTAACTGGCCGTCAATATTTGCCAATACGGGTTTCTTT"
            "CGCCCGGATGAAGCGCACCAGCCGCATTTCAGCGACCTGTTTGGGCAAATCATTAACGCC"
            "GGGCAAGGGGAAGGGCGCTATTCGGAGCTGCTGGCGATAAATCTGCTTGAGCAATTGTTA"
            "CTGCGGCGCATGGAAGCGATTAACGAGTCGCTCCATCCACCGATGGATAATCGGGTACGC"
            "GAGGCTTGTCAGTACATCAGCGATCACCTGGCAGACAGCAATTTTGATATCGCCAGCGTC"
            "GCACAGCATGTTTGCTTGTCGCCGTCGCGTCTGTCACATCTTTTCCGCCAGCAGTTAGGG"
            "ATTAGCGTCTTAAGCTGGCGCGAGGACCAACGCATTAGTCAGGCGAAGCTGCTTTTGAGC"
            "ACTACCCGGATGCCTATCGCCACCGTCGGTCGCAATGTTGGTTTTGACGATCAACTCTAT"
            "TTCTCGCGAGTATTTAAAAAATGCACCGGGGCCAGCCCGAGCGAGTTTCGTGCCGGTTGT"
            "GAAGAAAAAGTGAATGATGTAGCCGTCAAGTTGTCATAA"
        ),
        notes=(
            "Arabinoz yokken P-BAD'ı baskılar, arabinoz varken aktive "
            "eder - hem negatif hem pozitif düzenleyici. P-BAD kullanan "
            "her yapıda bulunmalıdır."
            "\n\n"
            "Dizi: NP_414606'nın coded_by alanı, 292 aa."
        ),
        accession="NC_000913.3 araC: 70387..71265",
        pmid="",
    ),
    _p(
        part_id="REG_TETR",
        name="tetR (TetR baskılayıcısı)",
        category="regulator",
        sequence=(
            "ATGTCTAGATTAGATAAAAGTAAAGTGATTAACAGCGCATTAGAGCTGCTTAATGAGGTC"
            "GGAATCGAAGGTTTAACAACCCGTAAACTCGCCCAGAAGCTAGGTGTAGAGCAGCCTACA"
            "TTGTATTGGCATGTAAAAAATAAGCGGGCTTTGCTCGACGCCTTAGCCATTGAGATGTTA"
            "GATAGGCACCATACTCACTTTTGCCCTTTAGAAGGGGAAAGCTGGCAAGATTTTTTACGT"
            "AATAACGCTAAAAGTTTTAGATGTGCTTTACTAAGTCATCGCGATGGAGCAAAAGTACAT"
            "TTAGGTACACGGCCTACAGAAAAACAGTATGAAACTCTCGAAAATCAATTAGCCTTTTTA"
            "TGCCAACAAGGTTTTTCACTAGAGAATGCATTATATGCACTCAGCGCTGTGGGGCATTTT"
            "ACTTTAGGTTGCGTATTGGAAGATCAAGAGCATCAAGTCGCTAAAGAAGAAAGGGAAACA"
            "CCTACTACTGATAGTATGCCGCCATTATTACGACAAGCTATCGAATTATTTGATCACCAA"
            "GGTGCAGAGCCAGCCTTCTTATTCGGCCTTGAATTGATCATATGCGGATTAGAAAAACAA"
            "CTTAAATGTGAAAGTGGGTCTTAA"
        ),
        notes=(
            "tetO'ya bağlanır, aTc ile serbest kalır. P-tet kullanan her "
            "yapıda gereklidir."
            "\n\n"
            "Dizi: Tn10 kaydının kendi tetR CDS anotasyonu (sınıf B TetR, 207 "
            "aa)."
        ),
        accession="J01830 (Tn10) tetR: complement(903..1526)",
        pmid="",
    ),
    _p(
        part_id="REG_ROP",
        name="rop / rom (kopya sayısı baskılayıcısı)",
        category="regulator",
        sequence=(
            "GTGACCAAACAGGAAAAAACCGCCCTTAACATGGCCCGCTTTATCAGAAGCCAGACATTA"
            "ACGCTTCTGGAGAAACTCAACGAGCTGGACGCGGATGAACAGGCAGACATCTGTGAATCG"
            "CTTCACGACCACGCTGATGAGCTTTACCGCAGCTGCCTCGCGCGTTTCGGTGATGACGGT"
            "GAAAACCTCTGA"
        ),
        notes=(
            "RNA I - RNA II eşleşmesini kuvvetlendirerek ColE1 ailesi "
            "orijinlerin kopya sayısını düşürür. pBR322'de vardır, "
            "pUC19'da silinmiştir - pUC'nin yüksek kopya sayısının "
            "ikinci nedeni budur."
            "\n\n"
            "Dizi: pBR322 kaydının kendi ROP CDS anotasyonu; başlama kodonu GTG."
        ),
        accession="J01749 (pBR322) rop: 1915..2106",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# REPORTER'LAR
# ---------------------------------------------------------------------------

REPORTERS = [
    _p(
        part_id="RPT_SFGFP",
        name="sfGFP (süper katlanan GFP)",
        category="reporter",
        sequence=(
            "ATGAGTAAAGGAGAAGAATTATTCACTGGAGTTGTCCCAATTTTAGTTGAATTAGATGGT"
            "GATGTTAATGGGCATAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCAACAAATGGA"
            "AAATTAACCTTAAAATTTATTTGTACTACTGGAAAACTACCTGTTCCATGGCCAACATTA"
            "GTCACTACTTTAACTTATGGTGTTCAATGTTTTTCAAGATACCCAGATCATATGAAACAG"
            "CATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGTACTATATCTTTC"
            "AAAGATGACGGGACTTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCTTAGTT"
            "AATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTTTAGGACATAAA"
            "TTGGAATACAACTTTAACTCACATAATGTATACATCACAGCAGACAAACAAAAGAATGGA"
            "ATCAAAGCTAACTTCAAAATTCGTCATAACGTTGAAGATGGAAGTGTTCAACTAGCAGAC"
            "CATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCTTATTACCAGACAACCATTAC"
            "CTGTCTACACAATCTGTTTTATCTAAAGATCCCAACGAAAAGCGTGACCATATGGTCCTT"
            "TTAGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA"
        ),
        notes=(
            "238 aa. Klasik GFP'den çok daha hızlı ve verimli katlanır; "
            "füzyon ortağı olarak kullanıldığında bile floresansı korur. "
            "E. coli'de kodon optimize varyantları yaygın. Promotör "
            "gücü ölçümü için standart reporter."
            "\n\n"
            "Dizi: avGFP (P42212) ile karşılaştırıldığında tam olarak 11 kanonik "
            "sfGFP mutasyonu çıkıyor (S30R, Y39N, F64L, S65T, F99S, N105T, Y145F, "
            "M153T, V163A, I171V, A206V). 'sfGFP' etiketli bir başka kayıt "
            "(JQ341914) bu imzayı taşımadığı için elendi - adı doğru olan dizinin "
            "doğru olduğunu göstermiyor."
        ),
        accession="MW132720 (sentetik sfGFP): 1..717",
        pmid="16369541",
    ),
    _p(
        part_id="RPT_MCHERRY",
        name="mCherry",
        category="reporter",
        sequence=(
            "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAG"
            "GTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGC"
            "CGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCC"
            "TTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCAC"
            "CCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGC"
            "GTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGAC"
            "GGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTA"
            "ATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGC"
            "GCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCT"
            "GAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTC"
            "AACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAA"
            "CGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"
        ),
        notes=(
            "Kırmızı floresan protein. sfGFP ile iki renkli deneylerde "
            "spektral çakışma az. Olgunlaşması GFP'den yavaştır - "
            "hızlı dinamik ölçümlerde bu bir sınırlamadır."
            "\n\n"
            "Dizi: Shaner ve ark. 2004'ün özgün kaydı; kaydın referansı bu "
            "parçanın pmid alanıyla aynı (15558047)."
        ),
        accession="AY678264 (mCherry, Shaner 2004): 1..711",
        pmid="15558047",
    ),
    _p(
        part_id="RPT_LACZ",
        name="lacZ (beta-galaktosidaz)",
        category="reporter",
        sequence=(
            "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCT"
            "GGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGC"
            "GAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGC"
            "TTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCT"
            "GAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATC"
            "TACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCG"
            "ACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACG"
            "CGAATTATTTTTGATGGCGTTAACTCGGCGTTTCATCTGTGGTGCAACGGGCGCTGGGTC"
            "GGTTACGGCCAGGACAGTCGTTTGCCGTCTGAATTTGACCTGAGCGCATTTTTACGCGCC"
            "GGAGAAAACCGCCTCGCGGTGATGGTGCTGCGCTGGAGTGACGGCAGTTATCTGGAAGAT"
            "CAGGATATGTGGCGGATGAGCGGCATTTTCCGTGACGTCTCGTTGCTGCATAAACCGACT"
            "ACACAAATCAGCGATTTCCATGTTGCCACTCGCTTTAATGATGATTTCAGCCGCGCTGTA"
            "CTGGAGGCTGAAGTTCAGATGTGCGGCGAGTTGCGTGACTACCTACGGGTAACAGTTTCT"
            "TTATGGCAGGGTGAAACGCAGGTCGCCAGCGGCACCGCGCCTTTCGGCGGTGAAATTATC"
            "GATGAGCGTGGTGGTTATGCCGATCGCGTCACACTACGTCTGAACGTCGAAAACCCGAAA"
            "CTGTGGAGCGCCGAAATCCCGAATCTCTATCGTGCGGTGGTTGAACTGCACACCGCCGAC"
            "GGCACGCTGATTGAAGCAGAAGCCTGCGATGTCGGTTTCCGCGAGGTGCGGATTGAAAAT"
            "GGTCTGCTGCTGCTGAACGGCAAGCCGTTGCTGATTCGAGGCGTTAACCGTCACGAGCAT"
            "CATCCTCTGCATGGTCAGGTCATGGATGAGCAGACGATGGTGCAGGATATCCTGCTGATG"
            "AAGCAGAACAACTTTAACGCCGTGCGCTGTTCGCATTATCCGAACCATCCGCTGTGGTAC"
            "ACGCTGTGCGACCGCTACGGCCTGTATGTGGTGGATGAAGCCAATATTGAAACCCACGGC"
            "ATGGTGCCAATGAATCGTCTGACCGATGATCCGCGCTGGCTACCGGCGATGAGCGAACGC"
            "GTAACGCGAATGGTGCAGCGCGATCGTAATCACCCGAGTGTGATCATCTGGTCGCTGGGG"
            "AATGAATCAGGCCACGGCGCTAATCACGACGCGCTGTATCGCTGGATCAAATCTGTCGAT"
            "CCTTCCCGCCCGGTGCAGTATGAAGGCGGCGGAGCCGACACCACGGCCACCGATATTATT"
            "TGCCCGATGTACGCGCGCGTGGATGAAGACCAGCCCTTCCCGGCTGTGCCGAAATGGTCC"
            "ATCAAAAAATGGCTTTCGCTACCTGGAGAGACGCGCCCGCTGATCCTTTGCGAATACGCC"
            "CACGCGATGGGTAACAGTCTTGGCGGTTTCGCTAAATACTGGCAGGCGTTTCGTCAGTAT"
            "CCCCGTTTACAGGGCGGCTTCGTCTGGGACTGGGTGGATCAGTCGCTGATTAAATATGAT"
            "GAAAACGGCAACCCGTGGTCGGCTTACGGCGGTGATTTTGGCGATACGCCGAACGATCGC"
            "CAGTTCTGTATGAACGGTCTGGTCTTTGCCGACCGCACGCCGCATCCAGCGCTGACGGAA"
            "GCAAAACACCAGCAGCAGTTTTTCCAGTTCCGTTTATCCGGGCAAACCATCGAAGTGACC"
            "AGCGAATACCTGTTCCGTCATAGCGATAACGAGCTCCTGCACTGGATGGTGGCGCTGGAT"
            "GGTAAGCCGCTGGCAAGCGGTGAAGTGCCTCTGGATGTCGCTCCACAAGGTAAACAGTTG"
            "ATTGAACTGCCTGAACTACCGCAGCCGGAGAGCGCCGGGCAACTCTGGCTCACAGTACGC"
            "GTAGTGCAACCGAACGCGACCGCATGGTCAGAAGCCGGGCACATCAGCGCCTGGCAGCAG"
            "TGGCGTCTGGCGGAAAACCTCAGTGTGACGCTCCCCGCCGCGTCCCACGCCATCCCGCAT"
            "CTGACCACCAGCGAAATGGATTTTTGCATCGAGCTGGGTAATAAGCGTTGGCAATTTAAC"
            "CGCCAGTCAGGCTTTCTTTCACAGATGTGGATTGGCGATAAAAAACAACTGCTGACGCCG"
            "CTGCGCGATCAGTTCACCCGTGCACCGCTGGATAACGACATTGGCGTAAGTGAAGCGACC"
            "CGCATTGACCCTAACGCCTGGGTCGAACGCTGGAAGGCGGCGGGCCATTACCAGGCCGAA"
            "GCAGCGTTGTTGCAGTGCACGGCAGATACACTTGCTGATGCGGTGCTGATTACGACCGCT"
            "CACGCGTGGCAGCATCAGGGGAAAACCTTATTTATCAGCCGGAAAACCTACCGGATTGAT"
            "GGTAGTGGTCAAATGGCGATTACCGTTGATGTTGAAGTGGCGAGCGATACACCGCATCCG"
            "GCGCGGATTGGCCTGAACTGCCAGCTGGCGCAGGTAGCAGAGCGGGTAAACTGGCTCGGA"
            "TTAGGGCCGCAAGAAAACTATCCCGACCGCCTTACTGCCGCCTGTTTTGACCGCTGGGAT"
            "CTGCCATTGTCAGACATGTATACCCCGTACGTCTTCCCGAGCGAAAACGGTCTGCGCTGC"
            "GGGACGCGCGAATTGAATTATGGCCCACACCAGTGGCGCGGCGACTTCCAGTTCAACATC"
            "AGCCGCTACAGTCAACAGCAACTGATGGAAACCAGCCATCGCCATCTGCTGCACGCGGAA"
            "GAAGGCACATGGCTGAATATCGACGGTTTCCATATGGGGATTGGTGGCGACGACTCCTGG"
            "AGCCCGTCAGTATCGGCGGAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTC"
            "TGGTGTCAAAAATAA"
        ),
        notes=(
            "1024 aa - oldukça uzun. X-gal ile mavi/beyaz koloni "
            "seçilimi ve Miller ünitesi ölçümü. Ekipman gerektirmediği "
            "için ders laboratuvarlarının klasiğidir."
            "\n\n"
            "Dizi: NP_414878'in coded_by alanı, 1024 aa."
        ),
        accession="NC_000913.3 lacZ: complement(363231..366305)",
        pmid="",
    ),
    _p(
        part_id="RPT_LACZA",
        name="lacZ-alfa (mavi/beyaz tarama için)",
        category="reporter",
        sequence=(
            "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTA"
            "CCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTT"
            "ACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAG"
            "GCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATG"
            "CGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGT"
            "ACAATCTGCTCTGATGCCGCATAG"
        ),
        requires=["lacZ-delta-M15 konak (DH5-alfa, JM109, TOP10)"],
        notes=(
            "lacZ'nin sadece N-terminal parçası. Konaktaki eksik "
            "lacZ-delta-M15 ile birleşerek (alfa-tamamlama) aktif enzim "
            "oluşturur. MCS bu parçanın İÇİNDE yer alır: insert girerse "
            "çerçeve bozulur, koloni beyaz kalır. "
            "pUC19'un mavi/beyaz seçilimi bu mekanizmadır. "
            "Konak lacZ-delta-M15 taşımıyorsa tarama ÇALIŞMAZ."
            "\n\n"
            "Dizi: pUC19 kaydında anotasyon yok; lacZ'nin 5' ucu pUC19'a eksi "
            "zincirde hizalanıp çerçeveye uyan ATG 469'da, ilk iç stop 146'da "
            "bulundu. Çeviri MTMIT (lacZ kodon 1-5) + polylinker + lacZ "
            "çerçevesinin devamı; MCS gerçekten parçanın içinde. Katalogdaki 450 "
            "bp tahmini, gerçek ORF 324 bp."
        ),
        accession="L09137 (pUC19) lacZ-alfa: complement(146..469)",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# ÇOKLU KLONLAMA BÖLGELERİ
# ---------------------------------------------------------------------------

MCS_PARTS = [
    _p(
        part_id="MCS_PUC19",
        name="pUC19 MCS",
        category="mcs",
        sequence="GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT",
        notes=(
            "Sıra: EcoRI - SacI - KpnI - SmaI/XmaI - BamHI - XbaI - SalI - "
            "PstI - SphI - HindIII. 13 enzim için 10 bölge. "
            "lacZ-alfa çerçevesinin içindedir, bu yüzden mavi/beyaz "
            "tarama ile birlikte çalışır. SmaI kör uç verir - yönlü "
            "klonlamada kullanılmaz."
        ),
        accession="pUC19: L09137",
        pmid="2985470",
    ),
    _p(
        part_id="MCS_PET28",
        name="pET-28a MCS (His-trombin-MCS-His)",
        category="mcs",
        sequence=(
            "CCATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCATATG"
            "GCTAGCATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGACAAGCTT"
            "GCGGCCGCACTCGAGCACCACCACCACCACCAC"
        ),
        notes=(
            "pET-28a'nın gerçek düzeni: NcoI - N-terminal 6xHis - trombin "
            "bölgesi - NdeI - T7 etiketi - BamHI/EcoRI/SacI/SalI/HindIII/"
            "NotI - XhoI - C-terminal 6xHis. "
            "ÖNEMLİ ÇERÇEVE KURALI: NcoI'dan girerseniz N-terminal His "
            "etiketini kaybedersiniz. XhoI'dan önce stop kodon koyarsanız "
            "C-terminal His etiketi okunmaz - bu en sık yapılan hatadır."
        ),
        accession="pET-28a(+)",
        pmid="2199796",
    ),
    _p(
        part_id="MCS_BIOBRICK",
        name="BioBrick ön/arka ek (RFC10)",
        category="mcs",
        sequence="GAATTCGCGGCCGCTTCTAGAG",
        notes=(
            "Ön ek: EcoRI - NotI - XbaI. Arka ek: SpeI - NotI - PstI. "
            "XbaI ve SpeI uyumlu uç verdiği için standart montaj "
            "sonrasında 'TACTAG' izi oluşur - bu dizi çerçeveyi kaydırır, "
            "bu yüzden RFC10 protein füzyonlarında kullanılmaz "
            "(onun için RFC25 vardır)."
        ),
        accession="iGEM RFC10",
        pmid="",
    ),
]

# ---------------------------------------------------------------------------
# Toplu liste
# ---------------------------------------------------------------------------

ALL_PARTS = (
    PROMOTERS
    + OPERATORS
    + RBS_PARTS
    + TAGS
    + PROTEASE_SITES
    + LINKERS
    + TERMINATORS
    + ORIGINS
    + MARKERS
    + REGULATORS
    + REPORTERS
    + MCS_PARTS
)


# ---------------------------------------------------------------------------
# Konak suşları
# ---------------------------------------------------------------------------

HOSTS = [
    {
        "host_id": "DH5A",
        "name": "DH5-alfa",
        "use": "klonlama",
        "t7_polymerase": False,
        "laciq": False,
        "lacz_dm15": True,
        "recA": False,
        "endA": False,
        "rare_trna": False,
        "genotype": "F- endA1 recA1 gyrA96 thi-1 relA1 glnV44 deoR phi80d lacZ-dM15",
        "notes": (
            "Klonlamanın standardı. recA1 eklenmiş insertin rekombinasyonla "
            "bozulmasını, endA1 ise plazmit saflaştırmada nükleaz "
            "bulaşmasını engeller. lacZ-dM15 taşıdığı için mavi/beyaz "
            "tarama çalışır. T7 polimeraz YOK - T7 promotörlü yapılar "
            "burada ifade edilmez."
        ),
    },
    {
        "host_id": "TOP10",
        "name": "TOP10",
        "use": "klonlama",
        "t7_polymerase": False,
        "laciq": False,
        "lacz_dm15": True,
        "recA": False,
        "endA": False,
        "rare_trna": False,
        "genotype": "F- mcrA D(mrr-hsdRMS-mcrBC) lacZ-dM15 araD139 D(araleu) recA1 endA1",
        "notes": (
            "DH5-alfa'ya benzer, ayrıca araD139 ve ara operonu delesyonu "
            "taşır - arabinoz tüketilmediği için P-BAD indüksiyonu daha "
            "kararlıdır. Metilasyona duyarlı bölgeleri (mcrA, mrr) "
            "silinmiş, bu yüzden dışarıdan gelen DNA'yı daha iyi kabul eder."
        ),
    },
    {
        "host_id": "BL21DE3",
        "name": "BL21(DE3)",
        "use": "ifade",
        "t7_polymerase": True,
        "laciq": True,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "F- ompT hsdSB(rB- mB-) gal dcm (DE3)",
        "notes": (
            "Protein ifadesinin standardı. DE3 lizojeni lacUV5 promotörü "
            "altında T7 RNA polimeraz taşır - T7/T7lac yapıları burada "
            "çalışır. ompT ve lon proteazları eksik olduğu için "
            "rekombinant protein daha az parçalanır. "
            "recA+ olduğundan tekrarlı dizi içeren plazmitlerin uzun "
            "süreli çoğaltılması için uygun DEĞİLDİR - klonlamayı "
            "DH5-alfa'da yapıp buraya aktarın."
        ),
    },
    {
        "host_id": "BL21DE3PLYSS",
        "name": "BL21(DE3) pLysS",
        "use": "ifade",
        "t7_polymerase": True,
        "laciq": True,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "BL21(DE3) + pLysS (CmR, p15A ori)",
        "notes": (
            "pLysS plazmiti T7 lizozimi üretir; bu protein T7 RNA "
            "polimerazı doğal olarak inhibe eder ve indüksiyon öncesi "
            "kaçak ifadeyi belirgin düşürür. Toksik proteinler için "
            "tercih edilir. pLysS kloramfenikol direnci ve p15A orijini "
            "taşır - ifade plazmitiniz CmR kullanıyorsa ya da p15A "
            "orijinliyse ÇATIŞIR."
        ),
    },
    {
        "host_id": "ROSETTA",
        "name": "Rosetta(DE3)",
        "use": "ifade",
        "t7_polymerase": True,
        "laciq": True,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": True,
        "genotype": "BL21(DE3) + pRARE (CmR, p15A ori)",
        "notes": (
            "pRARE plazmiti E. coli'de nadir olan kodonların (AGA/AGG, "
            "CGG, AUA, CUA, CCC, GGA) tRNA'larını sağlar. Ökaryot "
            "kaynaklı genleri kodon optimizasyonu yapmadan ifade "
            "etmek için. pRARE de CmR + p15A taşır - aynı çatışma riski."
        ),
    },
    {
        "host_id": "C41DE3",
        "name": "C41(DE3)",
        "use": "ifade",
        "t7_polymerase": True,
        "laciq": True,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "BL21(DE3) türevi, DE3 bölgesinde mutasyon",
        "notes": (
            "T7 polimeraz seviyesi düşük bir mutant. Membran proteinleri "
            "ve toksik hedefler için: yavaş ifade, daha iyi katlanma, "
            "daha az inklüzyon cismi."
        ),
    },
    {
        "host_id": "SHUFFLE",
        "name": "SHuffle T7",
        "use": "ifade",
        "t7_polymerase": True,
        "laciq": True,
        "lacz_dm15": False,
        "recA": True,
        "endA": True,
        "rare_trna": False,
        "genotype": "trxB gor suppressor mutant + sitoplazmik DsbC, (DE3)",
        "notes": (
            "Sitoplazma oksitleyici hale getirilmiş; disülfit bağı "
            "gerektiren proteinler (antikor fragmanları, bazı enzimler) "
            "periplazmik hedefleme olmadan doğru katlanabilir. "
            "Büyümesi normal BL21'den yavaştır."
        ),
    },
    {
        "host_id": "DB31",
        "name": "DB3.1",
        "use": "klonlama (ccdB dirençli)",
        "t7_polymerase": False,
        "laciq": False,
        "lacz_dm15": False,
        "recA": False,
        "endA": False,
        "rare_trna": False,
        "genotype": "gyrA462 (ccdB'ye dirençli)",
        "notes": (
            "Sadece ccdB taşıyan vektörleri çoğaltmak için. Giraz "
            "mutasyonu sayesinde ccdB zehirinden etkilenmez. "
            "Başka hiçbir amaçla kullanılmaz."
        ),
    },
]
