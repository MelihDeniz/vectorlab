"""Konak taşıma: bir kaseti başka bir bakteriye uyarlama.

Bu, aracın en ayırt edici parçasıdır ve sorduğu soru pratikte sürekli
sorulur: *"E. coli'de kurduğum bu kaseti Lactococcus lactis'e
taşıyabilir miyim?"* Mevcut dizi araçları bu soruya yanıt vermez;
diziyi gösterir, konağı hiç bilmez. Oysa bir kasetin konağa bağlı
parçaları vardır ve taşınırken DEĞİŞMESİ gerekir:

* **Promotör** — E. coli RNA polimerazının tanıdığı -35/-10 dizisi
  Gram-pozitif bir bakteride okunmayabilir. T7 promotörü ise hiçbir
  yere taşınmaz: T7 polimerazı yalnızca DE3 lizojenleri üretir.
* **Sinyal peptidi** — salgı yolu organizmaya özgüdür. DsbA, E. coli
  periplazmasına gönderir; Gram-pozitifte periplazma yoktur.
* **Marker** — her antibiyotik direnci her bakteride çalışmaz.
* **Orijin (replikon)** — en sert kısıt. ColE1 Gram-pozitifte
  çoğalmaz; taşımanın gerçekten tıkandığı yer çoğu zaman burasıdır.
* **RBS** — dizi taşınır ama **gücü taşınmaz.** Sentetik iGEM
  RBS'lerinin ölçüleri E. coli'de yapılmıştır.
* **Kodon kullanımı** — dizi aynı kalır, ama CAI değişir.

Modülün tasarım ilkesi: **taşınamayanı taşınmış gibi göstermemek.**
Katalogda hedef konakta çoğalan bir orijin yoksa, bu bir `tasinamaz`
kaydı olarak döner ve kullanıcı Referanstan kes-yapıştır sayfasındaki
gerçek iskeletlere yönlendirilir. Sessizce en yakın orijini koymak,
denetimden geçen ama laboratuvarda hiç çoğalmayan bir plazmit
üretirdi - bu aracın önlemeye çalıştığı hatanın ta kendisi.

Eşdeğerlik tablosundaki her hedef parça katalogda gerçekten bulunan,
dizisi kaynağından doğrulanmış bir parçadır; "şöyle bir promotör
olurdu" diye uydurulmuş kayıt yoktur.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core import codon, db, ref_vectors as rv

#: Konağa bağlı kategoriler. Bunun dışındaki parça (CDS, terminatör,
#: MCS, operatör...) taşınırken olduğu gibi kalır - terminatörün
#: çalışması da organizmaya bağlıdır ama katalogda hedefe özgü
#: terminatör olmadığı için uydurma yapmak yerine not düşülür.
HOST_DEPENDENT = ("promoter", "rbs", "signal", "marker", "ori")


@dataclass
class HostProfile:
    """Bir konağa taşırken hangi parçanın yerine ne konacağı.

    `ori` None ise: katalogda bu organizmada çoğalan bir orijin
    PARÇASI yok. Bu bir eksiklik değil, kasıtlı bir boşluk - Gram
    pozitif plazmitler yuvarlanan halka ile çoğalır ve çoğalma
    başlangıcı (dso) GenBank kayıtlarında ayrı bir öznitelik olarak
    işaretlenmez, rep geninin içinde durur. Sınırı uydurmadan parça
    yapılamaz; iskelet bütün olarak referans vektörden alınır.
    """

    host_id: str
    organism: str
    gram: str
    promoter: str
    signal: str
    marker: str
    ori: str | None
    rbs: str | None
    why: str
    #: Kategori -> o DEĞİŞİKLİĞE özgü gerekçe. Konak düzeyindeki
    #: `why` tabloda her satıra tekrar basılırsa tablo okunmaz hale
    #: gelir; kullanıcının satırda görmek istediği şey "promotör
    #: neden değişti" sorusunun yanıtıdır, konağın genel tanıtımı
    #: değil. `why` bir kez, başlıkta gösterilir.
    reasons: dict[str, str] = field(default_factory=dict)

    def reason_for(self, category: str) -> str:
        return self.reasons.get(category, self.why)


#: Katalogda gerçekten bulunan parçalarla kurulmuş eşdeğerlik tablosu.
#: Promotör, sinyal ve marker seçimleri hastalık kayıt defterinin
#: kendi seçimleriyle aynıdır (bkz. core/disease_seed.py
#: `_gram_poz_merdiven` ve `_GP_MARKER`): iki yerde ayrı karar
#: verilse biri değiştiğinde öbürü sessizce tutarsız olurdu.
PROFILES: dict[str, HostProfile] = {
    "LLNZ9000": HostProfile(
        host_id="LLNZ9000", organism="Lactococcus lactis",
        gram="pozitif",
        promoter="REF_NATIVE_LACTOBACILLUS_R_200",
        signal="SIG_USP45",
        marker="MRK_THYA_LL",
        ori=None,
        rbs=None,
        why=("Gıda sınıfı laktik asit bakterisi. Promotör pNZ8048 "
             "türevinden çıkarılmış doğal *L. reuteri* EF-Tu "
             "promotörüdür; salgı sinyali Usp45, laktik asit "
             "bakterilerinde en çok kullanılan sinyaldir. Marker "
             "yerine thyA tamamlaması önerilir: konağın thyA'sı "
             "silinince plazmitsiz hücre timinsiz ortamda ölür ve "
             "antibiyotik gerekmez."),
        reasons={
            "promoter": ("Doğal *L. reuteri* EF-Tu promotörü, pNZ8048 "
                         "türevinin GenBank öznitelik tablosundan "
                         "çıkarıldı. Gram-pozitif bir bakteride "
                         "okunduğu bilinen bir promotör gerekiyor; "
                         "E. coli -35/-10 dizisi burada "
                         "okunmayabilir."),
            "signal": ("Usp45, laktik asit bakterilerinde en çok "
                       "kullanılan salgı sinyalidir ve Steidler'in "
                       "IL-10 salgılayan *L. lactis* suşu da onu "
                       "kullanır. Gram-pozitifte periplazma yoktur: "
                       "protein doğrudan besiyerine çıkar, bu yüzden "
                       "E. coli periplazma sinyali (DsbA, PelB) "
                       "burada anlamsızdır."),
            "marker": ("thyA tamamlaması: konağın thyA'sı silinince "
                       "plazmitsiz hücre timinsiz ortamda ölür ve "
                       "antibiyotik gerekmez. Gıda sınıfı bir "
                       "organizmada antibiyotik direnç geni "
                       "düzenleyici açıdan sorundur. **Konağın ΔthyA "
                       "olması gerekir** - yabanıl NZ9000'de bu "
                       "marker seçilim sağlamaz."),
        },
    ),
    "BSUB168": HostProfile(
        host_id="BSUB168", organism="Bacillus subtilis",
        gram="pozitif",
        promoter="REF_PHELP_PROMOTER_211",
        signal="SIG_AMYE",
        marker="REF_ERM_C_735",
        ori=None,
        rbs=None,
        why=("Endüstriyel salgı konağı. pHelp promotörü ve AmyE "
             "sinyali ikisi de Gram-pozitif kaynaklıdır; AmyE, "
             "B. subtilis alfa-amilazının sinyalidir ve endüstriyel "
             "enzim üretiminin standart yoludur. Marker erm(C): "
             "katalogda thyA'sı silinmiş bir B. subtilis suşu yok, "
             "o yüzden burada antibiyotiksiz idame kurulamıyor - "
             "bu bir eksik, tercih değil."),
        reasons={
            "promoter": ("pHelp promotörü, pNZ8048 türevinin "
                         "öznitelik tablosundan çıkarıldı: "
                         "Gram-pozitif kaynaklı, güçlü ve "
                         "kurucu (constitutive) bir promotör."),
            "signal": ("AmyE, *B. subtilis* alfa-amilazının "
                       "sinyalidir ve endüstriyel enzim üretiminin "
                       "standart yoludur. Gram-pozitifte protein "
                       "periplazmada sıkışmaz, doğrudan besiyerine "
                       "salgılanır - saflaştırma bu yüzden kolaydır."),
            "marker": ("pE194'ten çıkarılmış erm(C): Gram-pozitif "
                       "bakterilerde çalışan gerçek bir seçilim "
                       "markeri. Burada thyA içermesi "
                       "kullanılamıyor çünkü katalogda thyA'sı "
                       "silinmiş bir *B. subtilis* suşu yok - "
                       "antibiyotiksiz idame bu konakta henüz "
                       "kurulamıyor ve bu bir eksik, tercih değil."),
        },
    ),
    "PPUTIDA": HostProfile(
        host_id="PPUTIDA", organism="Pseudomonas putida",
        gram="negatif",
        promoter="REF_PM_53",
        signal="SIG_PELB",
        marker="MRK_KAN",
        ori=None,
        rbs=None,
        why=("Gram-negatif olduğu için salgı yolu E. coli'ye "
             "benzer ve PelB sinyali periplazmaya gönderir. Promotör "
             "Pm (pSEVA128'den): SEVA iskeletlerinin geniş konaklı "
             "promotörü. Orijin için ColE1 YETMEZ - Pseudomonas'ta "
             "çoğalmaz; pBBR1, RK2, pRO1600 ya da RSF1010 gerekir ve "
             "bunlar katalogda parça değil, SEVA referans "
             "vektörlerinin içinde bulunur."),
        reasons={
            "promoter": ("Pm promotörü (pSEVA128'den): SEVA "
                         "iskeletlerinin geniş konaklı promotörü, "
                         "XylS düzenleyicisiyle çalışır. T7 "
                         "promotörü taşınamaz çünkü *P. putida* "
                         "KT2440 bir DE3 lizojeni değildir."),
            "signal": ("PelB: *P. putida* Gram-negatiftir, yani dış "
                       "zarı ve periplazması vardır. Salgı yolu "
                       "E. coli'ye benzer, bu yüzden Gram-pozitif "
                       "sinyali (Usp45, AmyE) yerine periplazmaya "
                       "gönderen bir sinyal doğru olanıdır."),
            "marker": ("Kanamisin direnci *P. putida*'da çalışır; "
                       "ampisilin çalışmaz, çünkü Pseudomonas "
                       "türlerinin doğal beta-laktam direnci "
                       "seçilimi anlamsız kılar."),
        },
    ),
    "CGLUT": HostProfile(
        host_id="CGLUT", organism="Corynebacterium glutamicum",
        gram="pozitif",
        promoter="REF_PHELP_PROMOTER_211",
        signal="SIG_AMYE",
        marker="MRK_KAN",
        ori=None,
        rbs=None,
        why=("Endüstriyel amino asit üreticisi; endotoksin "
             "üretmediği ve neredeyse proteaz salgılamadığı için "
             "ürün saflığında E. coli'ye üstündür. Hücre duvarı "
             "mikolik asit içerdiği için dönüşüm zordur. Vektörler "
             "pCG1 replikonu üzerine kurulur ve katalogda pCG1 "
             "parçası yok - iskelet referans vektörden alınmalı."),
        reasons={
            "promoter": ("Gram-pozitif kaynaklı güçlü bir kurucu "
                         "promotör. *C. glutamicum* promotör "
                         "tercihleri E. coli'den belirgin biçimde "
                         "farklıdır ve katalogda bu türde ölçülmüş "
                         "bir promotör yok: seçim en yakın "
                         "Gram-pozitif seçenektir, ölçülmüş bir "
                         "eşdeğer değil."),
            "signal": ("AmyE, Gram-pozitif Sec yolunu kullanır. "
                       "*C. glutamicum*'un hücre duvarı mikolik "
                       "asit içerir ve salgı verimi bu yüzden "
                       "*B. subtilis*'tekiyle aynı olmayabilir."),
            "marker": ("Kanamisin direnci *C. glutamicum*'da "
                       "yerleşik olarak kullanılır."),
        },
    ),
}

#: E. coli konakları: taşıma hedefi E. coli ise konağa bağlı parça
#: değişmez, yalnızca T7 gereksinimi denetlenir.
ECOLI_HOSTS = {"DH5A", "TOP10", "BL21DE3", "BL21DE3PLYSS", "ROSETTA",
               "C41DE3", "SHUFFLE", "DB31", "BL21STAR", "TUNER",
               "LEMO21", "NICO21", "ARCTIC", "JM109", "XL1BLUE",
               "JM110", "PIR1", "MG1655", "ECN", "ECN_DAPA"}

#: T7 polimeraz gerektiren promotörler - hiçbir konağa kendiliğinden
#: taşınmaz, hedef konağın DE3 lizojeni olması gerekir.
T7_PROMOTERS = {"PRO_T7", "PRO_T7LAC", "TER_T7"}


@dataclass
class Swap:
    """Tek bir parçaya ne olduğu."""

    part_id: str
    category: str
    action: str                 # kalir | degisti | tasinamaz
    replacement: str = ""
    reason: str = ""

    @property
    def symbol(self) -> str:
        return {"kalir": "=", "degisti": "→", "tasinamaz": "✗"}.get(
            self.action, "?")


@dataclass
class PortReport:
    target_host: str
    target_organism: str
    swaps: list[Swap] = field(default_factory=list)
    ported_parts: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    backbones: list[str] = field(default_factory=list)
    codon_report: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    @property
    def changed(self) -> list[Swap]:
        return [s for s in self.swaps if s.action == "degisti"]

    @property
    def blocked(self) -> list[Swap]:
        return [s for s in self.swaps if s.action == "tasinamaz"]

    @property
    def portable(self) -> bool:
        """Taşıma kaseti bozmadan tamamlanabildi mi?

        Orijin taşınamıyorsa yanıt HAYIR'dır, ama kaset yine de
        kullanılabilir: iskeleti referans vektörden almak gerekir.
        `blockers` ne yapılacağını söyler.
        """
        return not self.blocked


def profile(host_id: str) -> HostProfile | None:
    """Hedef konağın taşıma profili; E. coli için None (değişiklik yok)."""
    return PROFILES.get(host_id)


def backbones_for(organism: str) -> list[str]:
    """Hedef organizmada çoğalan referans vektörler.

    Liste `core/ref_vectors.py` kataloğundan gelir, yani her biri
    gerçek, tam dizili, aksesyonu yazılı bir plazmittir. Organizma
    etiketleri kataloğun kendi alanlarıdır; "geniş" ve "çok geniş"
    konak aralığı etiketleri de hedefi kapsar sayılır çünkü SEVA
    replikonlarının tanımı budur.
    """
    hedef = [organism]
    if "Lactococcus" in organism or "Bacillus" in organism:
        hedef += ["Gram-pozitif (B. subtilis, S. aureus)",
                  "Gram-pozitif (Lactococcus ve ötesi)",
                  "çok geniş (Gram-negatif, bazı Gram-pozitif)"]
    if "Pseudomonas" in organism:
        hedef += ["geniş (Gram-negatif)",
                  "çok geniş (Gram-negatif, bazı Gram-pozitif)",
                  "E. coli (mobilize edilerek Pseudomonas)"]
    if "Corynebacterium" in organism:
        hedef += ["E. coli (mobilize edilerek Corynebacterium, "
                  "Pseudomonas ve diğerleri)"]
    out: list[str] = []
    for org in hedef:
        for vid in rv.by_organism(org):
            if vid not in out:
                out.append(vid)
    return out


def _codon_report(part_ids: list[str], organism: str) -> dict[str, Any]:
    """Hedef konakta kodon uyumu.

    Hedef organizma için ölçülmüş kodon tablosu yoksa `tablo_yok`
    döner ve CAI HESAPLANMAZ. E. coli tablosuna düşmek, yanlış
    tabloyla ölçmek olur - `codon.table_key_for_organism` bu yüzden
    None döndürür ve biz de o kararı burada bozmayız.
    """
    anahtar = codon.table_key_for_organism(organism)
    if anahtar is None:
        return {"tablo_yok": True, "organizma": organism,
                "not": (f"{organism} için ölçülmüş kodon tablomuz yok, "
                        f"bu yüzden CAI hesaplanmadı. E. coli tablosuyla "
                        f"ölçmek yanlış sayı üretirdi.")}
    satirlar = []
    for pid in part_ids:
        p = db.get_part(pid)
        if not p or p.get("category") != "cds":
            continue
        dizi = (p.get("sequence") or "").upper()
        if len(dizi) < 60 or len(dizi) % 3:
            continue
        satirlar.append({
            "part_id": pid, "name": p.get("name", ""),
            "cai": codon.cai(dizi, org=anahtar),
            "cai_kaynak": codon.cai(dizi, org="ecoli"),
        })
    olcek = codon.reference_scale(anahtar)
    return {"tablo_yok": False, "organizma": organism, "tablo": anahtar,
            "satirlar": satirlar, "referans_olcek": olcek,
            "not": ("CAI, hedef konağın kendi ölçülmüş tablosuyla "
                    "hesaplandı. Karşılaştırma sütunu aynı dizinin "
                    "E. coli tablosundaki değeridir: iki sayının "
                    "farkı, optimizasyonun bu taşımada ne kadar "
                    "gerekli olduğunu gösterir.")}


def port(part_ids: list[str], target_host: str) -> PortReport:
    """Parça listesini hedef konağa uyarla.

    Dönen rapor üç şeyi ayrı ayrı söyler: neyin DEĞİŞTİĞİ, neyin
    OLDUĞU GİBİ KALDIĞI ve neyin TAŞINAMADIĞI. Üçüncüsü boş değilse
    `portable` False olur ve `blockers` ne yapılacağını yazar.
    """
    konak = db.get_host(target_host)
    organizma = (konak or {}).get("organism", "")
    prof = profile(target_host)
    rapor = PortReport(target_host=target_host, target_organism=organizma)

    ecoli_hedef = target_host in ECOLI_HOSTS or prof is None

    for pid in part_ids:
        p = db.get_part(pid) or {}
        kat = p.get("category", "")

        # T7: konak DE3 lizojeni değilse hiçbir yere taşınmaz.
        if pid in T7_PROMOTERS and kat == "promoter":
            if konak and konak.get("t7_polymerase"):
                rapor.swaps.append(Swap(pid, kat, "kalir", reason=(
                    "Hedef konak T7 polimeraz üretiyor (DE3 lizojeni), "
                    "promotör olduğu gibi çalışır.")))
                rapor.ported_parts.append(pid)
                continue
            yerine = prof.promoter if prof else ""
            rapor.swaps.append(Swap(
                pid, kat, "degisti" if yerine else "tasinamaz",
                replacement=yerine, reason=(
                    "T7 promotörü transkripsiyon için T7 RNA "
                    "polimerazı gerektirir ve bu enzim yalnızca DE3 "
                    "lizojenlerinde bulunur. Hedef konak DE3 değil, "
                    "bu yüzden promotör değiştirilmek zorunda.")))
            if yerine:
                rapor.ported_parts.append(yerine)
            else:
                rapor.blockers.append(
                    f"{pid}: hedef konakta T7 polimeraz yok ve yerine "
                    f"konacak promotör tanımlı değil. DE3 lizojeni bir "
                    f"konak seçin ya da konağa uygun bir promotör "
                    f"elle koyun.")
            continue

        if ecoli_hedef or kat not in HOST_DEPENDENT:
            sebep = ("Konağa bağlı olmayan kategori; olduğu gibi taşınır."
                     if kat not in HOST_DEPENDENT else
                     "Hedef konak da E. coli, konağa bağlı parça "
                     "değişmiyor.")
            rapor.swaps.append(Swap(pid, kat, "kalir", reason=sebep))
            rapor.ported_parts.append(pid)
            continue

        # Buradan sonra: hedef E. coli değil ve parça konağa bağlı.
        assert prof is not None
        hedef_pid = {"promoter": prof.promoter, "signal": prof.signal,
                     "marker": prof.marker, "ori": prof.ori,
                     "rbs": prof.rbs}.get(kat)

        if kat == "ori" and hedef_pid is None:
            rapor.swaps.append(Swap(pid, kat, "tasinamaz", reason=(
                f"{organizma} içinde çoğalan bir orijin PARÇASI "
                f"katalogda yok. Gram-pozitif plazmitler yuvarlanan "
                f"halka ile çoğalır ve çoğalma başlangıcı (dso) "
                f"GenBank kayıtlarında ayrı bir öznitelik olarak "
                f"işaretlenmez, rep geninin içinde durur - sınırı "
                f"uydurmadan parça yapılamaz.")))
            rapor.blockers.append(
                f"Orijin ({pid}) taşınamıyor. İskeleti parçadan kurmak "
                f"yerine BÜTÜN OLARAK bir referans vektörden alın: "
                f"Referanstan kes-yapıştır sayfası bunu yapar. "
                f"{organizma} için uygun iskeletler aşağıda listelendi.")
            continue

        if kat == "rbs" and hedef_pid is None:
            # RBS dizisi taşınır, çünkü katalogda hedefe özgü ölçülmüş
            # bir RBS yok ve rastgele başka bir RBS koymak hiçbir şeyi
            # düzeltmez - yalnızca değişiklik yapılmış gibi gösterir.
            # Taşınmayan şey güçtür ve bu söylenmek zorunda.
            gerekce = ("Dizi taşınır ama **gücü taşınmaz**: sentetik "
                       "iGEM RBS'lerinin ölçüleri E. coli'de "
                       "yapılmıştır. Katalogda bu organizmada ölçülmüş "
                       "bir RBS yok, bu yüzden parça olduğu gibi "
                       "kalıyor; rastgele başka bir RBS koymak "
                       "yalnızca değişiklik yapılmış gibi gösterirdi. "
                       "Beklediğiniz ifade düzeyini DENEYLE ölçmeniz "
                       "gerekir.")
            if prof.gram == "pozitif":
                gerekce += (" Gram-pozitiflerde başlangıç bölgesinin "
                            "ribozoma sunumu farklıdır; fark burada "
                            "daha büyük olabilir.")
            if pid == "RBS_PET":
                gerekce += (" Ayrıca bu RBS pET'in T7 ifade "
                            "sisteminden gelir ve promotör "
                            "değiştirildi: ikisi birlikte ölçülmüş bir "
                            "çift olmaktan çıktı.")
            rapor.swaps.append(Swap(pid, kat, "kalir", reason=gerekce))
            rapor.ported_parts.append(pid)
            continue

        if hedef_pid is None or hedef_pid == pid:
            rapor.swaps.append(Swap(pid, kat, "kalir", reason=(
                "Hedef konak için ayrı bir eşdeğer tanımlı değil; "
                "parça olduğu gibi taşınıyor.")))
            rapor.ported_parts.append(pid)
            continue

        yeni = db.get_part(hedef_pid) or {}
        rapor.swaps.append(Swap(pid, kat, "degisti", replacement=hedef_pid,
                                reason=prof.reason_for(kat)))
        rapor.ported_parts.append(hedef_pid)
        if not yeni:
            rapor.blockers.append(
                f"{hedef_pid} katalogda bulunamadı - eşdeğerlik tablosu "
                f"ile katalog arasında tutarsızlık var.")

    if not ecoli_hedef:
        rapor.backbones = backbones_for(organizma)
    rapor.codon_report = _codon_report(part_ids, organizma)

    rapor.notes = [
        "Taşıma bir ÖNERİDİR, deneysel olarak sınanmamıştır. "
        "Değişen her parça katalogda dizisi doğrulanmış gerçek bir "
        "parçadır, ama kasetin hedef konakta beklendiği gibi "
        "çalışacağını hiçbir tasarım aracı hesaplayamaz.",
        "Terminatör de organizmaya bağlıdır. Katalogda hedefe özgü "
        "terminatör olmadığı için olduğu gibi taşınıyor; rho'dan "
        "bağımsız bir saç tokası çoğu bakteride iş görür ama verimi "
        "değişir.",
        "Kodlayan diziler (CDS) olduğu gibi taşınır - değişen şey "
        "dizi değil, o dizinin hedef konaktaki kodon uyumudur.",
    ]
    if rapor.blocked:
        rapor.notes.insert(0, (
            f"**Taşıma tamamlanmadı:** {len(rapor.blocked)} parça "
            f"taşınamadı. Bu bir arıza değil, kataloğun bilinen sınırı - "
            f"ne yapılacağı aşağıda yazılı."))
    return rapor


def to_markdown(r: PortReport) -> str:
    s = [f"# Konak taşıma: {r.target_organism} ({r.target_host})", ""]
    s.append(f"**Sonuç:** "
             + ("kaset parça parça taşındı." if r.portable else
                f"{len(r.blocked)} parça taşınamadı, iskelet referans "
                f"vektörden alınmalı."))
    prof = profile(r.target_host)
    if prof:
        s += ["", f"**Hedef konak:** {prof.why}"]
    s += ["", "| Parça | Kategori | | Yerine | Gerekçe |",
          "|---|---|---|---|---|"]
    for sw in r.swaps:
        s.append(f"| {sw.part_id} | {sw.category} | {sw.symbol} | "
                 f"{sw.replacement or '—'} | "
                 f"{sw.reason.replace(chr(10), ' ')} |")
    s += ["", f"**Taşınan kaset:** `{' + '.join(r.ported_parts)}`", ""]

    if r.blockers:
        s += ["## Taşınamayanlar ve ne yapmalı", ""]
        for b in r.blockers:
            s.append(f"* {b}")
        s.append("")
    if r.backbones:
        s += [f"## {r.target_organism} içinde çoğalan referans vektörler",
              "", ", ".join(r.backbones), ""]
    k = r.codon_report
    if k:
        s += ["## Kodon uyumu", ""]
        if k.get("tablo_yok"):
            s.append(k["not"])
        else:
            s += [f"Tablo: `{k['tablo']}`", "",
                  "| Gen | CAI (hedef) | CAI (E. coli) |", "|---|---|---|"]
            for sat in k["satirlar"]:
                s.append(f"| {sat['name'] or sat['part_id']} | "
                         f"{sat['cai']:.3f} | {sat['cai_kaynak']:.3f} |")
            s += ["", k["not"]]
        s.append("")
    s += ["## Notlar", ""]
    for n in r.notes:
        s.append(f"* {n}")
    return "\n".join(s) + "\n"
