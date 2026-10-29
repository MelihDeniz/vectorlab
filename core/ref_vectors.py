"""Referans bakteriyel vektörler: gerçek, tam dizili plazmitler.

Diziler ve GenBank öznitelikleri `data/reference_vectors/*.json` içinde
durur; o dosyalar `tools/fetch_reference_vectors.py` ile NCBI'dan çekilir.
Uygulama çalışırken ağa çıkmaz.

Bu modüldeki `CURATED` sözlüğü editöryel katmandır: hangi plazmit ne işe
yarar, kopya sayısı, uyumsuzluk grubu, dikkat edilecekler. Dizi verisi
NCBI'dan, yorum buradan gelir - ikisi bilinçli olarak ayrı tutulur.
"""

from __future__ import annotations

import functools
import json
from pathlib import Path
from typing import Any

from core.seq_utils import (
    Feature,
    annotate,
    annotate_similar,
    gc_percent,
    revcomp,
)

DATA = Path(__file__).resolve().parent.parent / "data" / "reference_vectors"

# ---------------------------------------------------------------------------
# Editöryel katman
# ---------------------------------------------------------------------------
# `lit_query`: Literatür sayfasının Europe PMC'de arayacağı sorgu. Bilerek
# PMID yazmıyoruz - ezberden atıf numarası uydurmak sessiz hata üretir.

CURATED: dict[str, dict[str, Any]] = {
    "PUC19": {
        "name": "pUC19",
        "ozet": "Yüksek kopyalı standart klonlama vektörü; mavi/beyaz tarama.",
        "kullanim": "klonlama",
        "kopya": "yüksek (~500-700)",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac (lacUV5)",
        "mcs": "HindIII-SphI-PstI-SalI-XbaI-BamHI-SmaI-KpnI-SacI-EcoRI",
        "notlar": (
            "pBR322 iskeleti + M13mp19'dan gelen lac bölgesi. MCS, lacZ-alfa "
            "parçasının içine yerleştirilmiştir: insert girdiğinde alfa "
            "komplementasyon bozulur ve koloni X-gal üzerinde beyaz kalır. "
            "Bu taramanın çalışması için konağın lacZ-delta-M15 taşıması "
            "gerekir (DH5-alfa, TOP10). Kopya sayısının yüksekliği "
            "orijindeki nokta mutasyonundan gelir; rop/rom geni yoktur."
        ),
        "lit_query": "pUC19 cloning vector sequence Yanisch-Perron",
    },
    "PUC18": {
        "name": "pUC18",
        "ozet": "pUC19'un MCS'i ters yönde olan ikizi.",
        "kullanim": "klonlama",
        "kopya": "yüksek (~500-700)",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac (lacUV5)",
        "mcs": "EcoRI-SacI-KpnI-SmaI-BamHI-XbaI-SalI-PstI-SphI-HindIII",
        "notlar": (
            "Dizisi pUC19 ile neredeyse aynıdır, fark MCS'in yönüdür. "
            "Insertin lac promotörüne göre hangi yönde gireceğini seçmek "
            "istediğinizde ikisi arasından seçim yapılır."
        ),
        "lit_query": "pUC18 pUC19 vector Yanisch-Perron Vieira Messing",
    },
    "PBR322": {
        "name": "pBR322",
        "ozet": "Klonlamanın atası; iki ayrı direnç markerı taşır.",
        "kullanim": "klonlama",
        "kopya": "düşük (~15-20)",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)", "TetR (tetA)"],
        "promotor": "yok (kendi genlerinin promotörleri var)",
        "mcs": "yok - tek kesen bölgeler gen içlerine dağılmış",
        "notlar": (
            "Iki markerı olması insersiyonel inaktivasyona izin verir: "
            "BamHI ya da SalI ile tetA içine klonlarsanız koloni ampisiline "
            "dirençli ama tetrasikline duyarlı kalır. rop (rom) geni kopya "
            "sayısını baskılar - pUC19'un yüksek kopyalı olmasının bir "
            "nedeni de bu genin oradan çıkarılmış olmasıdır."
        ),
        "lit_query": "pBR322 construction characterization Bolivar Rodriguez",
    },
    "PACYC184": {
        "name": "pACYC184",
        "ozet": "p15A orijinli; ColE1 vektörlerle aynı hücrede durur.",
        "kullanim": "klonlama",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["CmR (cat)", "TetR (tetA)"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "Iki plazmitli kurulumların klasik ikinci üyesi: orijini ColE1 "
            "ile aynı uyumsuzluk grubunda olmadığı için pUC/pET ile birlikte "
            "tutulabilir. Rosetta ve pLysS suşlarının taşıdığı yardımcı "
            "plazmitler de p15A orijinlidir - üç plazmitli kurulumda bu "
            "çakışmaya dikkat edin."
        ),
        "lit_query": "pACYC184 p15A replicon Chang Cohen plasmid",
    },
    "PACYC177": {
        "name": "pACYC177",
        "ozet": "p15A orijinli, ampisilin + kanamisin dirençli.",
        "kullanim": "klonlama",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["AmpR (bla)", "KanR"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "pACYC184'ün kardeşi; marker çifti farklıdır. GenBank kaydı "
            "doğrusal olarak yatırılmıştır ama molekül daireseldir - "
            "uygulama bu kaydı dairesel olarak işler."
        ),
        "lit_query": "pACYC177 p15A cloning vector Chang Cohen",
    },
    "PBAD24": {
        "name": "pBAD24",
        "ozet": "Arabinozla indüklenen, sıkı kontrollü ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "araBAD (P-BAD), araC düzenleyicisiyle",
        "mcs": "NheI-KpnI-SmaI-XbaI-SalI-PstI-SphI-HindIII",
        "notlar": (
            "T7 sistemine göre asıl avantajı doz yanıtlı olmasıdır: "
            "arabinoz derişimini değiştirerek ifade düzeyini kademeli "
            "ayarlayabilirsiniz, IPTG'de bu kadar düzgün bir kademe yoktur. "
            "Konak arabinozu yıkmamalı (araBAD mutantı, örn. TOP10) ve "
            "glukoz katabolit baskılanması yaptığı için besiyerinde "
            "bulunmamalıdır."
        ),
        "lit_query": "pBAD24 tight regulation arabinose promoter Guzman",
    },
    "PBAD18": {
        "name": "pBAD18",
        "ozet": "pBAD24'ün MCS'li, çok yönlü kardeşi.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "araBAD (P-BAD)",
        "mcs": "pUC tipi genişletilmiş MCS",
        "notlar": (
            "pBAD24 ile aynı düzenleme mantığı; farkı klonlama bölgesidir. "
            "Başlangıç kodonunun promotöre uzaklığı pBAD24'te daha kısadır, "
            "pBAD18'de araya MCS girer."
        ),
        "lit_query": "pBAD18 arabinose expression vector Guzman Beckwith",
    },
    "PBAD33": {
        "name": "pBAD33",
        "ozet": "Arabinozla indüklenen, p15A orijinli ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "düşük (~10-12)",
        "inc_group": "p15A",
        "markerlar": ["CmR (cat)"],
        "promotor": "araBAD (P-BAD)",
        "mcs": "pBAD24 ile aynı dizilim",
        "notlar": (
            "Iki proteini ayrı plazmitlerden ifade edeceğiniz kurulumların "
            "arabinoz kolu: ColE1 orijinli bir pET ya da pUC ile birlikte "
            "durabilir ve markerı kloramfenikoldür, yani ampisilin ya da "
            "kanamisin ile çakışmaz."
        ),
        "lit_query": "pBAD33 arabinose vector p15A chloramphenicol Guzman",
    },
    "PTRC99A": {
        "name": "pTrc99A",
        "ozet": "IPTG ile indüklenen trc promotörlü ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "trc (trp -35 + lacUV5 -10)",
        "mcs": "NcoI-EcoRI-SmaI-BamHI-XbaI-SalI-PstI-HindIII",
        "notlar": (
            "T7 gerektirmez: trc promotörü konağın kendi RNA polimerazıyla "
            "çalışır, yani DE3 lizojeni olmayan suşlarda da ifade alırsınız. "
            "Plazmit kendi lacI genini taşır, bu yüzden yüksek kopyada bile "
            "baskılama ayakta kalır. Aşağı akışta rrnB T1T2 terminatörü var."
        ),
        "lit_query": "pTrc99A expression vector Amann Brosius trc promoter",
    },
    "PBSIIKS": {
        "name": "pBluescript II KS(+)",
        "ozet": "MCS'i T3 ve T7 promotörleri arasında; in vitro transkripsiyon.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "lac + MCS'i kuşatan T3 / T7",
        "mcs": "pUC tipi MCS, iki yanında T7 ve T3 promotörleri",
        "notlar": (
            "Inserti iki yönden de RNA'ya çevirebilirsiniz: bir yandan T7, "
            "diğer yandan T3 promotörü. f1 (M13) orijini sayesinde tek "
            "zincirli DNA üretilebilir - dizileme ve mutageneze yönelik "
            "tasarlanmıştır. lacZ-alfa taraması pUC19'da olduğu gibi çalışır."
        ),
        "lit_query": "pBluescript II phagemid vector Short Alting-Mees",
    },
    "PSC101": {
        "name": "pSC101",
        "ozet": "Doğal düşük kopyalı plazmit; üçüncü uyumsuzluk grubu.",
        "kullanim": "klonlama",
        "kopya": "çok düşük (~5)",
        "inc_group": "pSC101",
        "markerlar": ["TetR (tet)"],
        "promotor": "yok",
        "mcs": "yok",
        "notlar": (
            "Klonlama vektörü değil, doğal bir plazmit: moleküler biyolojide "
            "klonlanan ilk plazmit olduğu için ve orijini ColE1 ile p15A'nın "
            "ikisinden de ayrı olduğu için burada. Üç plazmitli kurulumlarda "
            "üçüncü orijin olarak bu replikon kullanılır. Kopya sayısı çok "
            "düşük olduğundan toksik ürünler için elverişlidir; buna karşılık "
            "plazmit izolasyonu az DNA verir."
        ),
        "lit_query": "pSC101 replicon copy number incompatibility Cohen Chang",
    },
    # --- iGEM BioBrick iskeletleri (dizi: iGEM Registry, CC BY 4.0) --------
    # Aşağıdaki içerik tahmin değil: her iskelet katalog parçalarıyla
    # taranarak doğrulandı. pSB1C3/1A3/1K3'te pUC orijini 212..866
    # konumunda, pSB4C5'te pSC101 orijini 689..2019 konumunda çıkıyor -
    # yani addaki "1" yüksek kopyayı, "4" pSC101'i gösteriyor.
    "PSB1C3": {
        "name": "pSB1C3",
        "ozet": "iGEM'in standart yüksek kopyalı montaj iskeleti (CmR).",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["CmR (cat)"],
        "promotor": "yok - parçanın kendi promotörü gelir",
        "mcs": "BioBrick ön eki (EcoRI-XbaI) ve arka eki (SpeI-PstI)",
        "notlar": (
            "iGEM dağıtımının varsayılan iskeleti: parçalar buraya ön ek ve "
            "arka ek arasına takılır. Insert bölgesi numaralandırmanın "
            "başlangıcını kuşatır - EcoRI ve XbaI dizinin sonunda, SpeI ve "
            "PstI başında görünür. RFC10 montajında ön ek EcoRI-XbaI, arka "
            "ek SpeI-PstI ile kesilir; XbaI ve SpeI uçları ligasyonda "
            "birleşerek üründe 6 bp'lik bir iz (scar) bırakır, bu yüzden "
            "protein füzyonları için uygun değildir."
        ),
        "lit_query": "BioBrick standard assembly RFC10 Knight registry",
    },
    "PSB1A3": {
        "name": "pSB1A3",
        "ozet": "pSB1C3'ün ampisilin dirençli eşi.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "pSB1C3 ile aynı iskelet, markerı ampisilindir. İki BioBrick "
            "yapısını aynı hücrede tutmak istediğinizde marker farkı işe "
            "yarar - ama orijinleri aynı olduğu için uyumsuzdur, uzun "
            "sürede biri kaybolur."
        ),
        "lit_query": "BioBrick standard assembly plasmid backbone registry",
    },
    "PSB1K3": {
        "name": "pSB1K3",
        "ozet": "pSB1C3'ün kanamisin dirençli eşi.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "ColE1 (pMB1 türevi)",
        "markerlar": ["KanR"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "Kanamisin seçilimi ampisiline göre daha temizdir: beta-laktamaz "
            "besiyerine salınıp antibiyotiği bozduğu için ampisilin "
            "plakalarında satellit koloniler çıkar, kanamisinde çıkmaz. "
            "Haritada direnç geni otomatik etiketlenmez: bu iskeletin "
            "kanamisin geni katalogdaki nptII varyantından yeterince "
            "farklı, benzerlik eşiğinin altında kalıyor."
        ),
        "lit_query": "BioBrick kanamycin backbone pSB1K3 registry",
    },
    "PSB4C5": {
        "name": "pSB4C5",
        "ozet": "Düşük kopyalı BioBrick iskeleti (pSC101 orijini, CmR).",
        "kullanim": "klonlama",
        "kopya": "düşük (pSC101)",
        "inc_group": "pSC101",
        "markerlar": ["CmR (cat)"],
        "promotor": "yok",
        "mcs": "BioBrick ön eki ve arka eki",
        "notlar": (
            "Yüksek kopyalı iskelet taşıyıcıya yük bindirdiğinde ya da "
            "ürün toksikse buna geçilir. Orijini pSC101 olduğu için "
            "pSB1 ailesiyle (ColE1) aynı hücrede durabilir - iki ayrı "
            "devreyi iki plazmitte taşımanın yolu bu."
        ),
        "lit_query": "pSB4C5 low copy BioBrick pSC101 backbone",
    },
}


# ---------------------------------------------------------------------------
# SEVA koleksiyonu
# ---------------------------------------------------------------------------
# Adlandırma sistematiktir: pSEVA[marker][orijin][kargo]. Aşağıdaki iki
# tablo TAHMİN DEĞİL - her kaydın kendi GenBank öznitelik tablosundan
# okundu. pSEVA111'de `bla` ve `R6K`, pSEVA211'de `neo` ve `R6K`,
# pSEVA221'de `neo` ve `trfA` görünüyor; haneler bu şekilde eşleşiyor.

SEVA_MARKER = {
    "1": ("AmpR (bla)", "ampisilin"),
    "2": ("KanR (neo)", "kanamisin"),
    "3": ("CmR (cat)", "kloramfenikol"),
    "4": ("SpecR/StrepR (aadA)", "spektinomisin / streptomisin"),
    "5": ("TetR (tetA)", "tetrasiklin"),
    "6": ("GmR (acc3)", "gentamisin"),
}

SEVA_ORI = {
    "1": ("R6K", "koşullu - yalnızca pir⁺ konakta çoğalır",
          "Replikasyonu için konağın *pi* proteinini (pir geni) "
          "üretmesi gerekir. Normal bir suşta hiç çoğalmaz; bu yüzden "
          "entegrasyon (suicide) vektörü olarak kullanılır: plazmit "
          "kromozoma girmezse kaybolur."),
    "2": ("RK2 (trfA)", "düşük kopya, çok geniş konak aralığı",
          "RK2/RP4 replikonu neredeyse bütün Gram-negatif bakterilerde "
          "çalışır. Kopya sayısı düşüktür, bu da kararlılık ve düşük "
          "metabolik yük demektir."),
    "3": ("pBBR1", "orta kopya, geniş konak aralığı",
          "Pseudomonas, Rhizobium ve birçok Gram-negatifte çalışır; "
          "ColE1 ile uyumsuz değildir, aynı hücrede durabilir."),
    "4": ("pRO1600/ColE1", "orta-yüksek kopya",
          "İki orijinli melez: E. coli'de ColE1, Pseudomonas'ta "
          "pRO1600 kolu çalışır. İki organizma arasında gidip gelen "
          "çalışmalar için elverişlidir."),
    "5": ("RSF1010", "düşük kopya, en geniş konak aralığı",
          "IncQ grubundan; bilinen en geniş konak aralıklı "
          "replikonlardan biridir. Kendi mobilizasyon genlerini de "
          "taşıdığı için konjugasyonla aktarılabilir."),
}


# Kargo hanesi. Bu tablo da TAHMİN DEĞİL: 99 kaydın öznitelik tabloları
# kargo hanesine göre gruplanıp okundu. "2" olanların hepsinde lacZa_19,
# "6" olanların hepsinde luxCDABE, "8" olanların hepsinde xylS ve Pm
# görünüyor.
SEVA_CARGO: dict[str, tuple[str, str]] = {
    "1": ("boş çoklu klonlama bölgesi",
          "Kargo yok, yalnızca standart klonlama bölgesi. Kendi "
          "kasetinizi kuracaksanız başlangıç noktası budur."),
    "2": ("lacZ-alfa (pUC19 dizilimi)",
          "Mavi/beyaz tarama yapılabilir. Klonlama bölgesi pUC19'un "
          "polilinker dizilimindedir; konağın lacZ-delta-M15 taşıması "
          "gerekir."),
    "2S": ("lacZ-alfa + I-SceI bölgesi",
           "Mavi/beyaz taramaya ek olarak I-SceI tanıma bölgesi var. "
           "I-SceI 18 bazlık bir diziyi tanır ve bakteri genomunda "
           "pratik olarak hiç bulunmaz; bu yüzden plazmiti hücre içinde "
           "tek bir yerden kesmek ve kromozoma entegrasyonu seçmek için "
           "kullanılır."),
    "3": ("lacZ-alfa (pUC18 dizilimi)",
          "2 numaralı kargoyla aynı, farkı klonlama bölgesinin ters "
          "yönde olmasıdır. Insertin promotöre göre hangi yönde "
          "gireceğini seçmek istediğinizde ikisi arasından seçersiniz."),
    "4": ("lacI-q + Ptrc promotörü",
          "IPTG ile indüklenen ifade kaseti. Plazmit kendi lacI-q "
          "genini taşıdığı için baskılama konaktan bağımsızdır - "
          "E. coli dışı bakterilerde çalışmanın yolu budur, çünkü o "
          "konaklarda kromozomal lacI yoktur."),
    "5": ("lacZ (tam gen)",
          "Alfa parçası değil, beta-galaktosidazın tamamı. Promotör "
          "gücü ölçmek için klasik muhabir: ONPG ile renk ölçümü "
          "yaparak ifadeyi sayısallaştırırsınız."),
    "5T": ("lacZ (terminatörlü sürüm)",
           "5 numaralı kargonun varyantı."),
    "6": ("luxCDABE ışıma operonu",
          "Hücre kendi ışığını üretir - substrat eklemeniz gerekmez, "
          "çünkü luxC, luxD ve luxE substratı da sentezler. Canlı "
          "hücrede zamana bağlı ifade izlemek için en pratik muhabirdir; "
          "floresan proteinlerin aksine olgunlaşma süresi yoktur."),
    "7C": ("CFP (camgöbeği floresan protein)",
           "Floresan muhabir. Aynı hücrede birden çok sinyali ayrı ayrı "
           "izlemek için 7C/7D/7R/7Y dörtlüsü farklı dalga boylarında "
           "çalışır."),
    "7D": ("DsRed2 (kırmızı floresan protein)",
           "Floresan muhabir; tetramerik yapısı füzyon proteinlerinde "
           "sorun çıkarabilir, serbest muhabir olarak kullanın."),
    "7R": ("mCherry (kırmızı floresan protein)",
           "Monomerik kırmızı muhabir. DsRed2'ye göre füzyonlarda daha "
           "güvenlidir."),
    "7Y": ("YFP (sarı floresan protein)",
           "Floresan muhabir; CFP ile birlikte FRET ölçümlerinde "
           "kullanılan çifttir."),
    "8": ("xylS + Pm promotörü",
          "m-toluat ile indüklenen ifade kaseti. Pseudomonas'ta çok "
          "kullanılır ve IPTG'ye göre üstünlüğü, indükleyicinin ucuz "
          "olması ve konağın onu tüketmemesidir. lacI/Ptrc sisteminin "
          "çalışmadığı konaklarda alternatiftir."),
}


def _seva_entry(code: str) -> dict[str, Any]:
    """pSEVA kodundan editöryel kayıt üret (örn. '241', '227C', '112S').

    Kod üç haneli: [marker][orijin][kargo]. Kargo hanesi harfle
    uzatılmış olabilir (2S, 5T, 7C...). Üç tablo da kayıtların kendi
    öznitelik tablolarından okunarak kuruldu.
    """
    m, o, kargo_kodu = code[0], code[1], code[2:]
    marker, antibiyotik = SEVA_MARKER[m]
    ori_ad, ori_ozet, ori_not = SEVA_ORI[o]
    kargo_ad, kargo_not = SEVA_CARGO.get(
        kargo_kodu, ("tanımlanmamış kargo", ""))

    kullanim = "ifade" if kargo_kodu in ("4", "8") else "klonlama"
    return {
        "name": f"pSEVA{code}",
        "ozet": f"SEVA: {marker.split()[0]} + {ori_ad} + {kargo_ad}.",
        "kullanim": kullanim,
        "kopya": ori_ozet,
        "inc_group": ori_ad,
        "markerlar": [marker],
        "promotor": ("Ptrc (IPTG)" if kargo_kodu == "4"
                     else "Pm (m-toluat)" if kargo_kodu == "8"
                     else "yok - kargo kendi promotörünü getirir"),
        "mcs": "SEVA standart çoklu klonlama bölgesi (PacI ve SpeI arasında)",
        "organizma": "geniş (Gram-negatif)",
        "seva_marker": m,
        "seva_ori": o,
        "seva_cargo": kargo_kodu,
        "notlar": (
            f"SEVA koleksiyonu, parçaları sabit tutup tek tek "
            f"değiştirebilmek için tasarlandı. Her iskelette aynı mimari "
            f"var: marker, orijin ve kargo üç ayrı modül halinde duruyor "
            f"ve aralarındaki sınırlar hep aynı enzimlerle kesiliyor. "
            f"Orijini değiştirmek istediğinizde yapının geri kalanını "
            f"yeniden kurmanız gerekmiyor - adın üç hanesi tam olarak bu "
            f"üç modülü gösterir.\n\n"
            f"**Marker ({m}):** {marker}, seçilim {antibiyotik} ile.\n\n"
            f"**Orijin ({o}):** {ori_ad}. {ori_not}\n\n"
            f"**Kargo ({kargo_kodu}):** {kargo_ad}. {kargo_not}\n\n"
            f"Kargo bölgesinin iki yanında T0 ve T1 terminatörleri var; "
            f"bu, kargo ile iskelet arasındaki okuma geçişini keser ve "
            f"ifadeyi ölçülebilir kılar. Ayrıca oriT taşıdığı için "
            f"konjugasyonla aktarılabilir - dönüştürülmesi zor "
            f"bakterilerde elektroporasyon yerine eşleşme kullanılır."
        ),
        "lit_query": "SEVA standard European vector architecture Silva-Rocha",
    }


# Diskteki her SEVA kaydı için editöryel kayıt üret - koleksiyon
# büyürse elle liste güncellemeye gerek kalmasın.
for _yol in sorted(DATA.glob("PSEVA*.json")):
    _kod = _yol.stem[5:]
    if len(_kod) >= 3 and _kod[0] in SEVA_MARKER and _kod[1] in SEVA_ORI:
        CURATED[f"PSEVA{_kod}"] = _seva_entry(_kod)


# ---------------------------------------------------------------------------
# Genom düzenleme ve özel amaçlı vektörler
# ---------------------------------------------------------------------------

CURATED.update({
    "PKD46": {
        "name": "pKD46",
        "ozet": "Lambda-Red rekombinaz plazmiti - kromozomda gen silme.",
        "kullanim": "genom düzenleme",
        "kopya": "düşük",
        "inc_group": "repA101ts (sıcaklığa duyarlı)",
        "markerlar": ["AmpR (bla)"],
        "promotor": "araBAD (P-BAD), araC ile",
        "mcs": "yok",
        "organizma": "E. coli",
        "notlar": (
            "Kromozomda gen silmenin standart yolu. Plazmit üç Red genini "
            "taşır (gam, bet, exo) ve arabinozla açılır: gam konağın "
            "doğrusal DNA'yı parçalayan RecBCD nükleazını bloke eder, "
            "bet ve exo ise kısa homoloji kollarıyla rekombinasyonu "
            "yapar. Böylece 40-50 bazlık homoloji taşıyan bir PCR ürünü "
            "kromozomdaki hedefin yerine geçebilir.\n\n"
            "Orijini sıcaklığa duyarlıdır: 30 °C'de çoğalır, 37 °C'nin "
            "üstünde kaybolur. İşiniz bittiğinde hücreyi 37 °C'de "
            "büyüterek plazmitten kurtulursunuz - ayrıca bir şey yapmaya "
            "gerek yoktur. Bu yüzden bütün çalışma 30 °C'de yürütülür."
        ),
        "lit_query": "one-step inactivation chromosomal genes PCR products "
                     "Datsenko Wanner",
    },
    "PKD4": {
        "name": "pKD4",
        "ozet": "FRT ile kuşatılmış kanamisin kaseti - kalıp plazmit.",
        "kullanim": "genom düzenleme",
        "kopya": "orta",
        "inc_group": "R6K (koşullu)",
        "markerlar": ["KanR"],
        "promotor": "kasetin kendi promotörü",
        "mcs": "yok",
        "organizma": "E. coli",
        "notlar": (
            "pKD46 ile birlikte kullanılır. Bu plazmitten, iki yanına "
            "hedef gene özgü 40-50 bazlık homoloji kolları eklenmiş "
            "primerlerle kanamisin kaseti çoğaltılır; ürün Red "
            "rekombinazı taşıyan hücreye verilince kromozomdaki hedefin "
            "yerine geçer.\n\n"
            "Kasetin iki yanındaki FRT bölgeleri işin ikinci yarısıdır: "
            "silme doğrulandıktan sonra FLP rekombinaz taşıyan bir "
            "plazmit (pCP20) verilir ve direnç kaseti kromozomdan "
            "çıkarılır. Geriye tek bir FRT izi kalır. Böylece aynı "
            "markerı bir sonraki silme için tekrar kullanabilirsiniz."
        ),
        "lit_query": "FRT flanked resistance cassette FLP recombinase "
                     "chromosomal deletion E coli",
    },
    "PKD3": {
        "name": "pKD3",
        "ozet": "FRT ile kuşatılmış kloramfenikol kaseti - kalıp plazmit.",
        "kullanim": "genom düzenleme",
        "kopya": "orta",
        "inc_group": "R6K (koşullu)",
        "markerlar": ["CmR (cat)"],
        "promotor": "kasetin kendi promotörü",
        "mcs": "yok",
        "organizma": "E. coli",
        "notlar": (
            "pKD4 ile aynı mantık, markerı kloramfenikoldür. Arka arkaya "
            "iki gen silecekseniz ikisini sırayla kullanmak işi "
            "kolaylaştırır: ilk silmede kanamisin, ikincisinde "
            "kloramfenikol seçilimi yaparsınız.\n\n"
            "Her iki kalıp plazmitin de orijini R6K'dır, yani ancak "
            "pir⁺ bir suşta çoğalır. Bu kasıtlıdır: PCR kalıbı olarak "
            "kullandığınız plazmitin hedef suşta çoğalmasını istemezsiniz, "
            "yoksa plazmitten gelen direnç, kromozoma giren kasetle "
            "karışır."
        ),
        "lit_query": "pKD3 chloramphenicol FRT cassette lambda Red",
    },
    "PBBR1MCS": {
        "name": "pBBR1MCS",
        "ozet": "Geniş konak aralıklı, orta kopyalı klonlama vektörü.",
        "kullanim": "klonlama",
        "kopya": "orta",
        "inc_group": "pBBR1",
        "markerlar": ["CmR (cat)"],
        "promotor": "lac (mavi/beyaz tarama ile)",
        "mcs": "pBluescript tipi MCS, lacZ-alfa içinde",
        "organizma": "geniş (Gram-negatif)",
        "notlar": (
            "E. coli dışındaki Gram-negatiflerde çalışmak istediğinizde "
            "ilk akla gelen iskeletlerden biri. Replikonu IncP, IncQ ve "
            "IncW gruplarının hiçbirine ait değildir, bu yüzden o "
            "gruplardaki plazmitlerle aynı hücrede durabilir.\n\n"
            "lacZ-alfa taşıdığı için E. coli'de mavi/beyaz tarama "
            "çalışır; klonlamayı E. coli'de yapıp doğruladıktan sonra "
            "asıl konağa aktarmak yaygın bir yoldur."
        ),
        "lit_query": "pBBR1MCS broad host range cloning vector Kovach",
    },
    "RSF1010": {
        "name": "RSF1010",
        "ozet": "Doğal IncQ plazmiti; en geniş konak aralıklı replikon.",
        "kullanim": "klonlama",
        "kopya": "düşük-orta",
        "inc_group": "IncQ (RSF1010)",
        "markerlar": ["SmR", "SulR"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "çok geniş (Gram-negatif, bazı Gram-pozitif)",
        "notlar": (
            "Klonlama vektörü değil, doğal bir plazmit - ama birçok "
            "geniş konaklı vektörün replikonu buradan gelir (SEVA'nın "
            "5 numaralı orijini dahil). Kendi replikasyon proteinlerini "
            "(repA, repB, repC) taşıdığı için konağın kendi replikasyon "
            "makinesine bağımlılığı azdır; konak aralığının bu kadar "
            "geniş olmasının sebebi budur.\n\n"
            "Mobilizasyon genleri de üzerinde olduğundan, yardımcı bir "
            "konjugasyon plazmiti varlığında başka bakterilere "
            "aktarılabilir."
        ),
        "lit_query": "RSF1010 IncQ broad host range replicon mobilization",
    },
    "PBHR1": {
        "name": "pBHR1",
        "ozet": "pBBR1 türevi, iki markerlı geniş konaklı vektör.",
        "kullanim": "klonlama",
        "kopya": "orta",
        "inc_group": "pBBR1",
        "markerlar": ["KanR", "CmR"],
        "promotor": "lac",
        "mcs": "MCS, lacZ-alfa içinde",
        "organizma": "geniş (Gram-negatif)",
        "notlar": (
            "pBBR1MCS ailesinden; farkı iki direnç markerı taşımasıdır. "
            "Biriyle seçilim yapıp diğerini insersiyonel inaktivasyon "
            "için kullanabilirsiniz, ya da konağa göre uygun olanı "
            "seçersiniz - bazı bakteriler kloramfenikole doğal dirençli "
            "olduğu için orada kanamisin gerekir."
        ),
        "lit_query": "pBHR1 broad host range vector pBBR1 derivative",
    },
})


# ---------------------------------------------------------------------------
# Gram-pozitif bakteriler
# ---------------------------------------------------------------------------

CURATED.update({
    "PUB110": {
        "name": "pUB110",
        "ozet": "B. subtilis'in klasik replikonu; kanamisin dirençli.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "pUB110 (rolling circle)",
        "markerlar": ["KanR/NeoR (aadD)", "BleoR (bleO)"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "Gram-pozitif (B. subtilis, S. aureus)",
        "notlar": (
            "*Bacillus subtilis* vektörlerinin çoğu bu replikondan "
            "türer. Çoğalma biçimi E. coli plazmitlerinden farklıdır: "
            "yuvarlanan halka (rolling circle) replikasyonu yapar ve "
            "ara ürün olarak tek zincirli DNA üretir. Bunun pratik "
            "sonucu şudur - uzun insertler (>5-6 kb) bu replikonda "
            "kararsızdır, silinmeye ve yeniden düzenlenmeye eğilimlidir. "
            "Büyük parçalar için theta replikasyonu yapan bir iskelet "
            "seçin.\n\n"
            "E. coli'de çoğalmaz; klonlamayı E. coli'de yapmak "
            "istiyorsanız iki orijinli bir mekik (shuttle) vektör "
            "gerekir."
        ),
        "lit_query": "pUB110 Bacillus subtilis rolling circle replication "
                     "plasmid stability",
    },
    "PC194": {
        "name": "pC194",
        "ozet": "Gram-pozitif kloramfenikol replikonu.",
        "kullanim": "klonlama",
        "kopya": "yüksek",
        "inc_group": "pC194 (rolling circle)",
        "markerlar": ["CmR (catA)"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "Gram-pozitif (B. subtilis, S. aureus)",
        "notlar": (
            "Gram-pozitif mekik vektörlerinin kloramfenikol kolu "
            "genellikle buradan alınır. pUB110 gibi yuvarlanan halka "
            "replikasyonu yapar, aynı kararlılık sınırı burada da "
            "geçerlidir.\n\n"
            "E. coli'deki `cat` geninden farklı bir alleldir; "
            "Gram-pozitif ribozom bağlanma bölgesiyle gelir, bu yüzden "
            "E. coli'de ifade edilmesi beklenmez."
        ),
        "lit_query": "pC194 chloramphenicol resistance Gram-positive "
                     "rolling circle plasmid",
    },
    "PE194": {
        "name": "pE194",
        "ozet": "Sıcaklığa duyarlı Gram-pozitif replikon; eritromisin.",
        "kullanim": "genom düzenleme",
        "kopya": "orta",
        "inc_group": "pE194 (rolling circle, ısıya duyarlı)",
        "markerlar": ["ErmR (ermC)"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "Gram-pozitif (B. subtilis, S. aureus)",
        "notlar": (
            "Asıl değeri sıcaklığa duyarlı olmasıdır: 30 °C'de çoğalır, "
            "42-45 °C'de çoğalamaz. Bu, Gram-pozitif bakterilerde "
            "kromozoma entegrasyonla gen değiştirmenin standart "
            "yöntemini mümkün kılar - plazmit izin verici sıcaklıkta "
            "çoğaltılır, sonra sıcaklık yükseltilip seçilim sürdürülür; "
            "ayakta kalabilen tek hücreler plazmiti kromozoma entegre "
            "etmiş olanlardır.\n\n"
            "pKD46'nın E. coli'de yaptığı işin Gram-pozitif karşılığı "
            "sayılabilir, yöntem farklı olsa da mantık aynıdır: işi "
            "bitince kendiliğinden gitsin."
        ),
        "lit_query": "pE194 temperature sensitive replicon Bacillus "
                     "allelic exchange",
    },
    "PNZ8048": {
        "name": "pNZ8048 türevi (kargo içerir)",
        "ozet": "Lactococcus nisinle indüklenen ifade vektörü.",
        "kullanim": "ifade",
        "kopya": "orta",
        "inc_group": "pSH71 / pWV01 türevi",
        "markerlar": ["CmR (cat)", "ErmR (ermC)"],
        "promotor": "nisinle indüklenen (NICE sistemi)",
        "mcs": "NcoI başlama kodonu üzerinde",
        "organizma": "Lactococcus lactis",
        "notlar": (
            "*Lactococcus lactis* için NICE (nisin controlled expression) "
            "sisteminin ifade plazmiti. İndükleyici nisindir - bu bir "
            "bakteriyosin, yani gıda sınıfı bir peptit; sistemin gıda "
            "biyoteknolojisinde tercih edilmesinin sebebi budur. "
            "Konağın nisK/nisR iki bileşenli algılayıcısını taşıması "
            "gerekir (NZ9000 gibi).\n\n"
            "DİKKAT: buradaki GenBank kaydı (MK297322) boş vektör "
            "değildir; içinde bir kargo (murine leptin) ve ek elemanlar "
            "var. Replikon, markerlar ve promotör mimarisi incelemek "
            "için doğrudur, ama bunu boş iskelet sanıp uzunluk hesabı "
            "yapmayın."
        ),
        "lit_query": "NICE nisin controlled expression Lactococcus lactis "
                     "pNZ8048",
    },
    "PWV01": {
        "name": "pWV01",
        "ozet": "Lactococcus doğal plazmiti; geniş Gram-pozitif replikon.",
        "kullanim": "klonlama",
        "kopya": "orta",
        "inc_group": "pWV01 (rolling circle)",
        "markerlar": ["yok (doğal plazmit)"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "Gram-pozitif (Lactococcus ve ötesi)",
        "notlar": (
            "Laktik asit bakterilerinde kullanılan vektörlerin büyük "
            "kısmı bu replikondan türer (pSH71, pNZ ailesi). Konak "
            "aralığı Gram-pozitifler arasında oldukça geniştir ve "
            "bazı koşullarda E. coli'de de çoğalabilir - bu yüzden "
            "mekik vektör kurmak için elverişlidir.\n\n"
            "Doğal plazmit olduğu için üzerinde direnç markerı yoktur; "
            "vektör kurarken marker dışarıdan eklenir."
        ),
        "lit_query": "pWV01 Lactococcus lactis replicon lactic acid "
                     "bacteria vector",
    },
    "PCG1": {
        "name": "pCG1",
        "ozet": "Corynebacterium glutamicum replikonu.",
        "kullanim": "klonlama",
        "kopya": "orta",
        "inc_group": "pCG1",
        "markerlar": ["yok (yalnızca rep bölgesi)"],
        "promotor": "yok",
        "mcs": "yok",
        "organizma": "Corynebacterium glutamicum",
        "notlar": (
            "*C. glutamicum* endüstriyel amino asit üretiminin "
            "(glutamat, lizin) ana organizmasıdır ve bu plazmit o "
            "organizmada kullanılan vektörlerin replikonudur - pEKEx ve "
            "pXMJ ailesi buradan türer.\n\n"
            "Kayıt yalnızca replikasyon bölgesini içerir; kullanılabilir "
            "bir vektör kurmak için marker ve klonlama bölgesi "
            "eklenmesi gerekir. E. coli orijiniyle birleştirilip mekik "
            "vektör yapılması standart yoldur, çünkü klonlama adımları "
            "E. coli'de çok daha hızlıdır."
        ),
        "lit_query": "pCG1 Corynebacterium glutamicum plasmid replicon "
                     "shuttle vector",
    },
})


# Varsayilan olarak E. coli kabul edilir; digerleri kendi kaydinda yazar.
for _vid, _rec in CURATED.items():
    _rec.setdefault("organizma", "E. coli")

ORDER = [
    # E. coli klonlama ve ifade
    "PUC19", "PUC18", "PBSIIKS", "PBR322", "PACYC184", "PACYC177",
    "PTRC99A", "PBAD24", "PBAD18", "PBAD33", "PSC101",
    # BioBrick iskeletleri
    "PSB1C3", "PSB1A3", "PSB1K3", "PSB4C5",
    # SEVA koleksiyonu ORDER'a asagida toplu eklenir.
    # Genis konak arali
    "PBBR1MCS", "PBHR1", "RSF1010",
    # Genom duzenleme
    "PKD46", "PKD4", "PKD3", "PE194",
    # Gram-pozitif ve diger bakteriler
    "PUB110", "PC194", "PNZ8048", "PWV01", "PCG1",
]


# SEVA kayitlari marker/orijin/kargo sirasiyla eklenir.
ORDER = ORDER + sorted(
    (p.stem for p in DATA.glob("PSEVA*.json")),
    key=lambda v: (v[5], v[6], v[7:]),
)


def organisms() -> list[str]:
    """Katalogdaki vektörlerin hedef organizmaları."""
    gorulen: list[str] = []
    for vid in available():
        rec = load(vid)
        o = (rec or {}).get("organizma", "E. coli")
        if o not in gorulen:
            gorulen.append(o)
    return gorulen


def seva_facets() -> dict[str, dict[str, str]]:
    """SEVA süzgeçleri için hane -> açıklama tabloları."""
    return {
        "marker": {k: v[0] for k, v in SEVA_MARKER.items()},
        "ori": {k: v[0] for k, v in SEVA_ORI.items()},
        "cargo": {k: v[0] for k, v in SEVA_CARGO.items()},
    }


def filter_vectors(
    organism: str | None = None,
    use: str | None = None,
    seva_marker: str | None = None,
    seva_ori: str | None = None,
    seva_cargo: str | None = None,
    search: str = "",
) -> list[str]:
    """Vektörleri süz. 126 kaydı tek listede sunmak kullanışsızdır."""
    out = []
    ara = (search or "").strip().lower()
    for vid in available():
        rec = load(vid)
        if not rec:
            continue
        if organism and rec.get("organizma") != organism:
            continue
        if use and rec.get("kullanim") != use:
            continue
        if seva_marker and rec.get("seva_marker") != seva_marker:
            continue
        if seva_ori and rec.get("seva_ori") != seva_ori:
            continue
        if seva_cargo and rec.get("seva_cargo") != seva_cargo:
            continue
        if ara:
            govde = " ".join(str(rec.get(k, "")) for k in
                             ("name", "ozet", "inc_group", "accession"))
            govde += " " + " ".join(rec.get("markerlar") or [])
            if ara not in govde.lower():
                continue
        out.append(vid)
    return out


def by_organism(organism: str) -> list[str]:
    return [v for v in available()
            if (load(v) or {}).get("organizma") == organism]


# ---------------------------------------------------------------------------
# Yükleme
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=32)
def load(vector_id: str) -> dict[str, Any] | None:
    """Tam kaydı oku (dizi dahil). Bulunamazsa None."""
    path = DATA / f"{vector_id}.json"
    if not path.exists():
        return None
    rec = json.loads(path.read_text(encoding="utf-8"))
    rec.update(CURATED.get(vector_id, {}))
    # Dosya adındaki kimlik her zaman kazanır.
    rec["vector_id"] = vector_id
    rec.setdefault("name", vector_id)
    rec["gc"] = gc_percent(rec["sequence"])
    return rec


def available() -> list[str]:
    """Diskte verisi olan vektör kimlikleri, editöryel sırayla."""
    have = {p.stem for p in DATA.glob("*.json")}
    return [v for v in ORDER if v in have] + sorted(have - set(ORDER))


def list_vectors() -> list[dict[str, Any]]:
    """Listeleme için özet kayıtlar (dizi taşımaz)."""
    out = []
    for vid in available():
        rec = load(vid)
        if not rec:
            continue
        out.append({k: v for k, v in rec.items()
                    if k not in ("sequence", "features")})
    return out


# ---------------------------------------------------------------------------
# Öznitelik katmanları
# ---------------------------------------------------------------------------


def _spans(f: Feature, length: int) -> list[tuple[int, int]]:
    if f.start <= f.end:
        return [(f.start, f.end)]
    return [(f.start, length), (1, f.end)]


def _overlaps(a: Feature, b: Feature, length: int) -> bool:
    """Iki öznitelik büyük ölçüde üst üste mi biniyor (dairesel)."""
    shared = 0
    for s1, e1 in _spans(a, length):
        for s2, e2 in _spans(b, length):
            shared += max(0, min(e1, e2) - max(s1, s2) + 1)
    return shared >= 0.6 * min(a.length, b.length)


def catalog_features(
    vector_id: str, parts: list[dict[str, Any]], min_identity: float = 90.0
) -> list[Feature]:
    """Katalog parçalarını vektörde bul; çakışanlardan en iyisini tut.

    `pUC ori` ile `pBR322 ori` arasında tek nükleotid fark vardır, bu
    yüzden ikisi de aynı bölgeye düşer. Aynı kategoriden çakışan
    eşleşmelerde önce tam eşleşme, sonra en yüksek benzerlik kazanır -
    yoksa harita aynı bölgeyi iki kez etiketler.
    """
    rec = load(vector_id)
    if not rec:
        return []
    L = rec["length"]
    exact = annotate(rec["sequence"], parts, circular=True)
    fuzzy = annotate_similar(rec["sequence"], parts, circular=True,
                             min_identity=min_identity)

    keep: list[Feature] = list(exact)
    for f in sorted(fuzzy, key=lambda x: (-x.identity, -x.length)):
        if any(k.part_id == f.part_id and k.start == f.start for k in keep):
            continue
        if any(k.category == f.category and _overlaps(k, f, L) for k in keep):
            continue
        keep.append(f)
    return sorted(keep, key=lambda f: (f.start, -f.length))


def genbank_features(vector_id: str) -> list[dict[str, Any]]:
    """Kaydın kendi öznitelik tablosu (NCBI'dan geldiği gibi)."""
    rec = load(vector_id)
    return list(rec.get("features") or []) if rec else []


# ---------------------------------------------------------------------------
# Bölge alma ve dışa aktarım
# ---------------------------------------------------------------------------


def region(vector_id: str, start: int, end: int, strand: int = 1) -> str:
    """1 tabanlı, kapsayıcı bölgeyi döndür; dairesel sarmayı destekler."""
    rec = load(vector_id)
    if not rec:
        return ""
    seq, L = rec["sequence"], rec["length"]
    start = max(1, min(L, int(start)))
    end = max(1, min(L, int(end)))
    sub = seq[start - 1:end] if start <= end else seq[start - 1:] + seq[:end]
    return revcomp(sub) if strand == -1 else sub


def to_fasta(vector_id: str) -> str:
    rec = load(vector_id)
    if not rec:
        return ""
    body = "\n".join(rec["sequence"][i:i + 70]
                     for i in range(0, rec["length"], 70))
    return (f">{rec['name']} | {rec['accession']} | {rec['length']} bp | "
            f"dairesel\n{body}\n")
