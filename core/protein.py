"""İfade edilecek proteinin özellikleri ve bunların pratik karşılığı.

Sayıların kendisi Biopython'un `Bio.SeqUtils.ProtParam` modülünden gelir
(ExPASy ProtParam ile aynı yöntemler). Bu modülün kattığı şey yorumdur:
hangi sayı laboratuvarda neye karşılık geliyor, hangi değer sorun
habercisi.

Dört soru için yazıldı:
  - SDS-PAGE'de bandı nerede arayacağım?
  - Saflaştırdıktan sonra derişimi nasıl ölçeceğim?
  - Çözünür ifade şansım ne, inclusion body riski var mı?
  - Tamponun pH'ı ne olmalı?
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ProtParam yalnızca 20 standart amino asidi tanır; seçenek/durdurma
# kodonları önce ayıklanmalı yoksa KeyError atar.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Yaygın afinite etiketlerinin protein dizisindeki karşılığı. Etiket
# bulunduğunda kütle katkısı ve saflaştırma yöntemi ayrıca raporlanır.
TAG_MOTIFS: dict[str, tuple[str, str]] = {
    "HHHHHH": ("His6", "IMAC (Ni-NTA ya da Co-TALON) ile saflaştırılır"),
    "HHHHHHHH": ("His8", "IMAC; His6'ya göre daha sıkı bağlanır"),
    "DYKDDDDK": ("FLAG", "anti-FLAG afinite reçinesi; pahalı ama çok seçici"),
    "EQKLISEEDL": ("Myc", "çoğunlukla saflaştırma değil, tespit etiketi"),
    "WSHPQFEK": ("Strep-tag II", "Strep-Tactin; desthiobiyotinle nazik elüsyon"),
}

# Proteaz kesim bölgeleri - etiketi sonradan koparmak için.
PROTEASE_MOTIFS: dict[str, tuple[str, str]] = {
    "ENLYFQ": ("TEV", "TEV proteazı ENLYFQ/G arasından keser"),
    "LVPRGS": ("Trombin", "LVPR/GS arasından keser"),
    "IEGR": ("Faktör Xa", "IEGR/ sonrasından keser"),
    "DDDDK": ("Enterokinaz", "DDDDK/ sonrasından keser"),
}


@dataclass
class ProteinReport:
    length_aa: int
    mw_da: float
    mw_kda: float
    pi: float
    extinction_reduced: int        # M-1 cm-1, sistinler indirgenmiş
    extinction_cystine: int        # M-1 cm-1, sistin köprüleri varsayımıyla
    a280_for_1mg_ml: float | None  # 1 mg/mL çözeltinin A280 değeri
    instability: float
    gravy: float
    aromaticity: float
    charge_ph7: float
    helix: float
    turn: float
    sheet: float
    aa_counts: dict[str, int]
    tags: list[tuple[str, str]] = field(default_factory=list)
    proteases: list[tuple[str, str]] = field(default_factory=list)
    notes: list[dict[str, str]] = field(default_factory=list)


def _clean_protein(seq: str) -> str:
    """Durdurma kodonu ve standart dışı harfleri ayıkla."""
    s = (seq or "").upper().replace(" ", "").replace("\n", "")
    if "*" in s:
        s = s.split("*", 1)[0]
    return "".join(c for c in s if c in STANDARD_AA)


def _note(kind: str, title: str, detail: str) -> dict[str, str]:
    return {"kind": kind, "title": title, "detail": detail}


def analyze(protein_seq: str, working_ph: float = 8.0) -> ProteinReport | None:
    """Protein dizisinden özellikleri hesapla ve yorumla.

    `working_ph` saflaştırma tamponunun pH'ıdır; IMAC için tipik olarak
    7.5-8.0. pI bu değere çok yakınsa protein çöker.
    """
    s = _clean_protein(protein_seq)
    if len(s) < 10:
        return None

    pa = ProteinAnalysis(s)
    eps_red, eps_cys = pa.molar_extinction_coefficient()
    mw = pa.molecular_weight()
    pi = pa.isoelectric_point()
    helix, turn, sheet = pa.secondary_structure_fraction()
    gravy = pa.gravy()
    instab = pa.instability_index()

    rep = ProteinReport(
        length_aa=len(s),
        mw_da=round(mw, 1),
        mw_kda=round(mw / 1000.0, 2),
        pi=round(pi, 2),
        extinction_reduced=int(eps_red),
        extinction_cystine=int(eps_cys),
        a280_for_1mg_ml=round(eps_red / mw, 3) if eps_red else None,
        instability=round(instab, 1),
        gravy=round(gravy, 3),
        aromaticity=round(pa.aromaticity(), 3),
        charge_ph7=round(pa.charge_at_pH(7.0), 2),
        helix=round(helix, 3),
        turn=round(turn, 3),
        sheet=round(sheet, 3),
        aa_counts=dict(pa.count_amino_acids()),
    )

    for motif, (name, how) in TAG_MOTIFS.items():
        if motif in s:
            rep.tags.append((name, how))
    # His8 bulunduysa His6 da eşleşir; uzun olanı tut.
    if any(t[0] == "His8" for t in rep.tags):
        rep.tags = [t for t in rep.tags if t[0] != "His6"]
    for motif, (name, how) in PROTEASE_MOTIFS.items():
        if motif in s:
            rep.proteases.append((name, how))

    rep.notes = _interpret(rep, working_ph)
    return rep


def _interpret(r: ProteinReport, working_ph: float) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []

    # --- SDS-PAGE
    out.append(_note(
        "bilgi", "SDS-PAGE'de bandı nerede ara",
        f"Beklenen kütle {r.mw_kda} kDa. Jelde bandı bu boyutun çevresinde "
        f"arayın. Gerçek göç her zaman birebir değildir: asidik ve prolince "
        f"zengin proteinler olduğundan ağır, çok bazik olanlar hafif "
        f"görünebilir. İndüksiyon sonrası örnekle indüksiyon öncesi örneği "
        f"yan yana yüklerseniz fark zaten gözle seçilir."
    ))

    # --- Derişim ölçümü
    if not r.extinction_reduced:
        out.append(_note(
            "uyari", "A280 ile derişim ölçülemez",
            "Dizide triptofan ve tirozin yok, bu yüzden 280 nm'de ışığı "
            "ölçülebilir biçimde soğurmaz. Nanodrop sıfıra yakın değer "
            "verir ve bu protein olmadığı anlamına GELMEZ. Derişimi "
            "Bradford ya da BCA gibi kolorimetrik bir yöntemle ölçün."
        ))
    else:
        out.append(_note(
            "bilgi", "A280 ile derişim ölçümü",
            f"Sönüm katsayısı {r.extinction_reduced} M⁻¹cm⁻¹. 1 mg/mL "
            f"çözelti A280 = {r.a280_for_1mg_ml} verir, yani "
            f"**derişim (mg/mL) ≈ A280 / {r.a280_for_1mg_ml}**. "
            + ("Dizide sistein var: disülfit köprüleri oluşursa katsayı "
               f"{r.extinction_cystine} M⁻¹cm⁻¹'e çıkar. Sitoplazmada "
               "ortam indirgeyici olduğu için genelde ilk değeri kullanın."
               if r.extinction_cystine != r.extinction_reduced else "")
        ))

    # --- Çözünürlük
    if r.gravy > 0.0:
        out.append(_note(
            "uyari", "Hidrofobik: inclusion body riski yüksek",
            f"GRAVY {r.gravy} (pozitif = ortalama hidrofobik). Bu tür "
            f"proteinler E. coli sitoplazmasında sıklıkla çökelti "
            f"(inclusion body) olarak birikir. Karşı önlemler: indüksiyonu "
            f"düşük sıcaklıkta yapın (16-20 °C, gece boyu), indükleyici "
            f"derişimini düşürün, daha zayıf bir promotör ya da RBS seçin. "
            f"Dördüncü bir yol daha var: proteini periplazmaya "
            f"göndermek. Orada ortam oksitleyicidir, şaperon dağılımı "
            f"farklıdır ve birçok hedef sitoplazmada çökerken orada "
            f"çözünür kalır - katalogdaki sinyal peptitlerinden birini "
            f"(pelB, OmpA, DsbA) kasetin başına ekleyin. Hiçbiri "
            f"işe yaramazsa çökeltiyi üre ile çözüp yeniden katlamak "
            f"son seçenektir, ama verim düşer."
        ))
    elif r.gravy < -0.4:
        out.append(_note(
            "bilgi", "Hidrofilik: çözünür ifade şansı iyi",
            f"GRAVY {r.gravy} (belirgin negatif). Protein ortalama olarak "
            f"suyu seviyor; standart 37 °C indüksiyonla çözünür kısımda "
            f"bulma şansınız yüksek."
        ))

    if r.instability > 40:
        out.append(_note(
            "uyari", "Kararsızlık indeksi yüksek",
            f"İndeks {r.instability} (40'ın üstü 'kararsız' sayılır). Bu "
            f"hesaplanmış bir tahmindir, ölçüm değil; ama hücre içinde "
            f"hızlı bozunma ihtimaline işaret eder. Hasat sonrası "
            f"örnekleri buzda tutun, liziz tamponuna proteaz inhibitörü "
            f"ekleyin ve saflaştırmayı aynı gün bitirin."
        ))

    # --- pH / pI
    fark = abs(r.pi - working_ph)
    if fark < 1.0:
        out.append(_note(
            "uyari", "Tampon pH'ı pI'ya çok yakın",
            f"pI {r.pi}, çalışma pH'ı {working_ph}. Bir proteinin net yükü "
            f"pI'da sıfırdır; molekiller birbirini itmediği için topaklanıp "
            f"çöker. Tamponun pH'ını pI'dan en az 1 birim uzağa alın "
            f"(örneğin pH {round(r.pi + 1.5, 1)} ya da "
            f"pH {round(r.pi - 1.5, 1)}), ya da 100-300 mM NaCl ekleyin."
        ))
    else:
        yon = "negatif" if working_ph > r.pi else "pozitif"
        kolon = "anyon değiştirici (Q)" if working_ph > r.pi else "katyon değiştirici (SP)"
        out.append(_note(
            "bilgi", "Tampon pH'ı ve iyon değişim kolonu",
            f"pI {r.pi}; pH {working_ph} tamponunda protein net {yon} "
            f"yüklüdür. İyon değişimiyle saflaştıracaksanız {kolon} "
            f"kolon uygundur."
        ))

    # --- Sistein
    cys = r.aa_counts.get("C", 0)
    if cys >= 2:
        out.append(_note(
            "bilgi", f"{cys} sistein var",
            "E. coli sitoplazması indirgeyicidir, disülfit köprüleri "
            "orada düzgün kurulmaz. Protein disülfit gerektiriyorsa "
            "SHuffle gibi oksitleyici sitoplazmalı bir suş kullanın ya da "
            "periplazmik hedefleme sinyali ekleyin - katalogda pelB, "
            "OmpA, PhoA, MalE ve DsbA sinyalleri var, kesim noktaları "
            "UniProt'tan doğrulanmış. Periplazma oksitleyici bir ortamdır "
            "ve DsbA/DsbC sistemi disülfitleri orada kurar. Liziz tamponuna DTT "
            "ya da TCEP eklerseniz köprüler açık kalır - IMAC reçinesiyle "
            "yüksek DTT'nin uyuşmadığını unutmayın (2 mM TCEP güvenlidir)."
        ))

    # --- Etiketler
    if r.tags:
        adlar = ", ".join(t[0] for t in r.tags)
        nasil = "; ".join(f"{t[0]}: {t[1]}" for t in r.tags)
        out.append(_note(
            "bilgi", f"Afinite etiketi bulundu ({adlar})",
            f"{nasil}. Etiketin kendi kütlesi yukarıdaki toplam kütleye "
            f"dahildir."
        ))
        if not r.proteases:
            out.append(_note(
                "bilgi", "Etiketi koparma yolu yok",
                "Yapıda proteaz kesim bölgesi görünmüyor; etiket son "
                "üründe kalır. Çoğu uygulamada sorun değildir, ama "
                "kristalografi ya da aktivite ölçümü yapacaksanız araya "
                "bir TEV bölgesi (ENLYFQ/G) eklemeyi düşünün."
            ))
    if r.proteases:
        out.append(_note(
            "bilgi", "Proteaz kesim bölgesi var",
            "; ".join(f"{p[0]}: {p[1]}" for p in r.proteases)
            + ". Kesimden sonra proteazı ve kopan etiketi ikinci bir "
            "IMAC geçişiyle uzaklaştırabilirsiniz (ikisi de kolona "
            "bağlanır, hedef protein akıp geçer)."
        ))

    return out


def sds_ladder_note(mw_kda: float) -> str:
    """Hangi jel yüzdesi bu boyut için uygun."""
    if mw_kda < 15:
        return ("%15-18 akrilamid jel kullanın; düşük yüzdeli jelde bu "
                "boyuttaki bant boya cephesiyle birlikte çıkar.")
    if mw_kda < 60:
        return "%12 akrilamid jel bu aralık için standarttır."
    if mw_kda < 150:
        return "%10 akrilamid jel kullanın."
    return ("%8 ya da daha düşük yüzdeli jel kullanın; yüksek yüzdeli "
            "jelde bu boyuttaki bant kuyuya yakın sıkışır.")
