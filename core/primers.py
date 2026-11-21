"""Vektöre özel primer tasarımı ve primer inceleme.

`core/seq_utils.py` içindeki `design_primers` bir primer çiftini
*hedef dizinin kendisine* göre tasarlar: Tm, GC, dimer, saç tokası.
Bunlar gerekli ama yeterli değil, çünkü primer laboratuvarda hedef
diziye değil **bütün plazmite** girer.

Bu modülün eklediği şey budur:

- **Özgüllük.** Primer plazmitte başka bir yere de bağlanıyorsa PCR ya
  yanlış ürün verir ya hiç vermez. Tm ve GC bunu söylemez; yalnızca
  diziyi taramak söyler.
- **Koloni PCR.** İnsertin iki yanındaki iskelet bölgesinden primer
  tasarlar, böylece tek bir PCR boş vektörü dolu vektörden ayırır.
  Vektöre özel olmasının sebebi açık: primerler insertin değil,
  vektörün dizisinden gelir.
- **Dizileme çapaları.** Vektörde zaten bulunan katalog parçalarını
  (T7, lac, araBAD promotörü gibi) arar ve oradan okumaya başlayacak
  primer önerir.
- **Primer inceleme.** Elinizdeki primeri yapıştırıp ölçtüren araç.

**Modelin sınırı.** Bağlanma taraması bir *hizalama* değildir: boşluk
(indel) aranmaz ve termodinamik hesap yapılmaz. Eşleşmeyen baz sayılır
ve 3' ucunun tam eşleşmesi şart koşulur. Gerekçe biyolojiktir -
polimeraz uzamayı 3' ucundan başlatır, o uç oturmuyorsa primer
bağlansa bile ürün çıkmaz. Bu yüzden 5' ucundaki uyumsuzluk bağışlanır,
3' ucundaki bağışlanmaz. Gerçek bir özgüllük denetimi BLAST'tır; bu
tarama onun yerine geçmez, yalnızca **aynı plazmit içindeki** ikinci
bağlanma yerini yakalar.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core.seq_utils import (
    Primer,
    clean_seq,
    gc_percent,
    primer_dimer,
    primer_hairpin,
    revcomp,
    tm_nn,
)

#: 3' ucunda tam eşleşmesi gereken baz sayısı. Polimeraz uzamayı
#: buradan başlatır; son üç baz oturmuyorsa uzama başlamaz.
THREE_PRIME_EXACT = 3

#: Bağlanma sayılırken izin verilen toplam uyumsuz baz.
DEFAULT_MAX_MISMATCH = 2


@dataclass
class BindingSite:
    position: int          # 1 tabanlı, primerin 5' ucunun oturduğu yer
    strand: int            # +1 üst zincir, -1 alt zincir
    mismatches: int
    end_position: int      # primerin 3' ucunun oturduğu yer


@dataclass
class ColonyPCR:
    forward: Primer
    reverse: Primer
    with_insert: int       # insert varken beklenen ürün (bp)
    without_insert: int    # boş vektörde beklenen ürün (bp)
    notes: list[dict[str, str]] = field(default_factory=list)

    @property
    def fark(self) -> int:
        return self.with_insert - self.without_insert


def _note(kind: str, title: str, detail: str) -> dict[str, str]:
    return {"kind": kind, "title": title, "detail": detail}


# ---------------------------------------------------------------------------
# Özgüllük
# ---------------------------------------------------------------------------


def binding_sites(primer: str, template: str, circular: bool = True,
                  max_mismatch: int = DEFAULT_MAX_MISMATCH
                  ) -> list[BindingSite]:
    """Primerin kalıpta bağlanabileceği yerler.

    Her konum için uyumsuz baz sayılır; `max_mismatch`'i aşan konumlar
    elenir ve 3' ucundaki son `THREE_PRIME_EXACT` bazın tam eşleşmesi
    şart koşulur (bkz. modül açıklamasındaki gerekçe).
    """
    p = clean_seq(primer)
    t = clean_seq(template)
    if not p or not t or len(p) > len(t):
        return []

    # Dairesel kalıpta primer başlangıç noktasının üstünden geçebilir.
    tarama = t + t[: len(p) - 1] if circular else t

    out: list[BindingSite] = []
    for strand in (1, -1):
        # Alt zincirde primer, kalıbın ters tümleyenine oturur. Aramayı
        # primerin ters tümleyeniyle üst zincirde yapmak aynı kapıya
        # çıkar ve koordinatlar kalıbın kendi numaralandırmasında kalır.
        ara = p if strand == 1 else revcomp(p)
        for i in range(len(tarama) - len(ara) + 1):
            pencere = tarama[i:i + len(ara)]
            # 3' uç: üst zincirde dizinin sonu, alt zincirde başı.
            uc = (pencere[-THREE_PRIME_EXACT:] if strand == 1
                  else pencere[:THREE_PRIME_EXACT])
            uc_primer = (ara[-THREE_PRIME_EXACT:] if strand == 1
                         else ara[:THREE_PRIME_EXACT])
            if uc != uc_primer:
                continue
            yanlis = sum(1 for a, b in zip(pencere, ara) if a != b)
            if yanlis > max_mismatch:
                continue
            bas = i % len(t) + 1
            son = (i + len(ara) - 1) % len(t) + 1
            out.append(BindingSite(
                position=bas if strand == 1 else son,
                strand=strand,
                mismatches=yanlis,
                end_position=son if strand == 1 else bas,
            ))
    return out


def specificity(primer: str, template: str, circular: bool = True,
                max_mismatch: int = DEFAULT_MAX_MISMATCH) -> dict[str, Any]:
    """Primer bu plazmitte tek bir yere mi bağlanıyor?"""
    yerler = binding_sites(primer, template, circular, max_mismatch)
    tam = [y for y in yerler if y.mismatches == 0]
    notlar: list[dict[str, str]] = []

    if not yerler:
        notlar.append(_note(
            "hata", "Primer bu dizide hiç bağlanmıyor",
            "Ne tam ne de yaklaşık bir bağlanma yeri bulundu. Primer "
            "başka bir kalıp için tasarlanmış ya da yönü ters olabilir."))
    elif len(yerler) == 1:
        notlar.append(_note(
            "bilgi", "Tek bağlanma yeri",
            f"Primer plazmitte yalnızca {yerler[0].position}. konumda "
            f"({'üst' if yerler[0].strand > 0 else 'alt'} zincir) "
            f"bağlanıyor. Özgüllük bakımından sorun yok."))
    else:
        ikincil = len(yerler) - len(tam) if tam else len(yerler) - 1
        notlar.append(_note(
            "uyari", f"{len(yerler)} bağlanma yeri",
            f"Primer plazmitte birden fazla yere oturuyor "
            f"({len(tam)} tam eşleşme, {max(ikincil, 0)} yaklaşık). "
            f"Konumlar: "
            + ", ".join(f"{y.position}"
                        f"{'+' if y.strand > 0 else '-'}"
                        f"{'' if y.mismatches == 0 else f'({y.mismatches} uyumsuz)'}"
                        for y in yerler[:8])
            + ". Birden fazla bağlanma yan ürün verir ya da asıl ürünü "
              "bastırır; primeri birkaç baz kaydırmak genelde yeter."))
    return {"sites": yerler, "unique": len(yerler) == 1, "notes": notlar}


# ---------------------------------------------------------------------------
# Koloni PCR
# ---------------------------------------------------------------------------


def _en_iyi_uzunluk(seq: str, target_tm: float,
                    lo: int = 18, hi: int = 30) -> str:
    """Hedef Tm'e en yakın uzunluk. `seq` 5'->3' yönünde verilir."""
    en_iyi, fark = seq[:lo], 1e9
    for n in range(lo, min(hi, len(seq)) + 1):
        t = tm_nn(seq[:n])
        if t is None:
            continue
        if abs(t - target_tm) < fark:
            en_iyi, fark = seq[:n], abs(t - target_tm)
    return en_iyi


def colony_pcr(vector: str, insert_start: int, insert_end: int,
               circular: bool = True, target_tm: float = 58.0,
               offset: int = 60, name_prefix: str = "kPCR") -> ColonyPCR:
    """İnsertin iki yanındaki İSKELETTEN koloni tarama primerleri.

    `insert_start` ve `insert_end` 1 tabanlı ve kapsayıcıdır.

    Primerler bilerek insertin dışına, vektörün kendi dizisine
    konur. Sebebi koloni PCR'ın ne işe yaradığıdır: ligasyondan çıkan
    kolonilerin bir kısmı kendi üzerine kapanmış boş vektördür ve
    insertin içine bağlanan bir primer boş vektörde hiç ürün vermez -
    yani boş koloni ile başarısız PCR birbirinden ayırt edilemez.
    İskelete bağlanan primer ikisinde de ürün verir: küçük bant boş,
    büyük bant dolu demektir. Bant yoksa PCR'ın kendisi tutmamıştır.
    """
    v = clean_seq(vector)
    n = len(v)
    notlar: list[dict[str, str]] = []
    if n == 0 or insert_start < 1 or insert_end > n or insert_start > insert_end:
        return ColonyPCR(
            forward=Primer("", "", None, 0.0, 0),
            reverse=Primer("", "", None, 0.0, 0),
            with_insert=0, without_insert=0,
            notes=[_note("hata", "Insert koordinatları geçersiz",
                         f"1-{n} aralığında bir bölge bekleniyor, "
                         f"{insert_start}-{insert_end} verildi.")])

    insert_uz = insert_end - insert_start + 1

    # İleri primer insertin SOLUNDAKİ iskelette, geri primer SAĞINDAKİ
    # iskelette. Dairesel plazmitte iki uç da sarılabilir.
    def dilim(bas: int, uzunluk: int) -> str:
        """1 tabanlı `bas` konumundan `uzunluk` baz; dairesel sarar."""
        bas0 = (bas - 1) % n
        if bas0 + uzunluk <= n:
            return v[bas0:bas0 + uzunluk]
        return v[bas0:] + v[: bas0 + uzunluk - n]

    ileri_bas = insert_start - offset
    if not circular and ileri_bas < 1:
        ileri_bas = 1
        notlar.append(_note(
            "uyari", "İleri primer insertin yakınına sıkıştı",
            "Doğrusal dizide insertin solunda yeterli iskelet yok; "
            "primer olabildiğince sola yerleştirildi."))
    ileri_bolge = dilim(ileri_bas, min(offset, 40))
    ileri_cekirdek = _en_iyi_uzunluk(ileri_bolge, target_tm)

    geri_bas = insert_end + 1
    geri_bolge = dilim(geri_bas, min(offset, 40))
    if not circular and insert_end + min(offset, 40) > n:
        notlar.append(_note(
            "uyari", "Geri primer insertin yakınına sıkıştı",
            "Doğrusal dizide insertin sağında yeterli iskelet yok."))
    geri_cekirdek = _en_iyi_uzunluk(revcomp(geri_bolge), target_tm)

    ileri = Primer(
        name=f"{name_prefix}_F", sequence=ileri_cekirdek,
        tm=tm_nn(ileri_cekirdek), gc=gc_percent(ileri_cekirdek),
        length=len(ileri_cekirdek))
    geri = Primer(
        name=f"{name_prefix}_R", sequence=geri_cekirdek,
        tm=tm_nn(geri_cekirdek), gc=gc_percent(geri_cekirdek),
        length=len(geri_cekirdek))

    # Ürün: ileri primerin 5' ucundan geri primerin 5' ucuna.
    sol_payanda = insert_start - ileri_bas
    sag_payanda = len(geri_cekirdek)
    dolu = sol_payanda + insert_uz + sag_payanda
    bos = sol_payanda + sag_payanda

    # Dimer ve saç tokası, tasarımcıdaki ile aynı ölçütle.
    for pr, etiket in ((ileri, "İleri"), (geri, "Geri")):
        kendi = primer_dimer(pr.sequence)
        if kendi["risk"] != "yok":
            pr.warnings.append(f"{etiket} primer {kendi['note']}")
        toka = primer_hairpin(pr.sequence)
        if toka["risk"] != "yok":
            pr.warnings.append(f"{etiket} primer: {toka['note']}")

    capraz = primer_dimer(ileri.sequence, geri.sequence)
    if capraz["risk"] != "yok":
        notlar.append(_note("uyari", "Primerler birbirine bağlanıyor",
                            capraz["note"]))

    # Özgüllük: ikisi de vektörde tek yere bağlanmalı.
    for pr, etiket in ((ileri, "İleri"), (geri, "Geri")):
        oz = specificity(pr.sequence, v, circular=circular)
        if not oz["unique"]:
            notlar.append(_note(
                "uyari", f"{etiket} primer vektörde tek yere bağlanmıyor",
                oz["notes"][0]["detail"]))

    # Jelde ayırt edilebilirlik: koloni PCR'ın bütün anlamı iki bandı
    # birbirinden ayırabilmektir.
    if dolu - bos < 100:
        notlar.append(_note(
            "uyari", "İki bant jelde ayrılmayabilir",
            f"Dolu {dolu} bp, boş {bos} bp - aradaki fark "
            f"{dolu - bos} bp. Agaroz jelde bu kadar yakın iki bant "
            f"aynı yerde görünür. Insert kısa olduğu için primerleri "
            f"daha içeriye almak işe yaramaz; bu durumda tarama "
            f"kesimle ya da dizilemeyle yapılır."))
    else:
        notlar.append(_note(
            "bilgi", "Beklenen bantlar",
            f"Insert giren koloni **{dolu} bp**, boş vektör "
            f"**{bos} bp** verir. Hiç bant yoksa PCR tutmamıştır - "
            f"boş koloniyle karıştırmayın."))

    return ColonyPCR(forward=ileri, reverse=geri,
                     with_insert=dolu, without_insert=bos, notes=notlar)


# ---------------------------------------------------------------------------
# Dizileme çapaları
# ---------------------------------------------------------------------------


def catalog_anchors(vector: str, parts: list[dict[str, Any]],
                    circular: bool = True, read_from: int = 20,
                    target_tm: float = 55.0) -> list[dict[str, Any]]:
    """Vektörde bulunan katalog parçalarından dizileme primeri öner.

    Evrensel dizileme primerlerinin (M13, T7, SP6) dizileri bu depoya
    EZBERDEN yazılmadı. Bunun yerine primer, vektörün kendi dizisinden
    kesilir: katalogda dizisi doğrulanmış bir promotör bu vektörde
    bulunuyorsa, onun hemen ardından okumaya başlayacak bir primer
    zaten vektörün içinde demektir. Böylece önerilen her primerin
    kaynağı bellidir ve uydurulmuş dizi olmaz.
    """
    v = clean_seq(vector)
    if not v:
        return []
    out: list[dict[str, Any]] = []
    for p in parts:
        dizi = clean_seq(p.get("sequence") or "")
        # Dizisi olmayan ya da çok kısa parçadan çapa çıkmaz.
        if len(dizi) < 18 or "N" in dizi:
            continue
        if p.get("category") not in ("promoter", "mcs", "operator"):
            continue
        arama = v + v[: len(dizi) - 1] if circular else v
        i = arama.find(dizi)
        if i < 0:
            continue
        bas = i % len(v) + 1
        son = (i + len(dizi) - 1) % len(v) + 1
        # Primer parçanın SON kısmından alınır ve ileriye okur.
        cekirdek = _en_iyi_uzunluk(dizi[-min(len(dizi), 30):], target_tm)
        out.append({
            "part_id": p.get("part_id", ""),
            "name": p.get("name", ""),
            "start": bas,
            "end": son,
            "primer": cekirdek,
            "tm": tm_nn(cekirdek),
            "gc": gc_percent(cekirdek),
            "reads_from": (son + read_from) % len(v) + 1,
        })
    return out


# ---------------------------------------------------------------------------
# Primer inceleme
# ---------------------------------------------------------------------------


def analyze(primer: str, template: str = "", circular: bool = True
            ) -> dict[str, Any]:
    """Elinizdeki bir primeri ölç ve yorumla.

    `template` verilirse özgüllük de denetlenir.
    """
    s = clean_seq(primer)
    notlar: list[dict[str, str]] = []
    if len(s) < 10:
        return {"error": "Primer çok kısa (<10 nt) ya da geçerli baz yok."}

    tm = tm_nn(s)
    gc = gc_percent(s)
    olcum = {
        "sequence": s, "length": len(s), "tm": tm, "gc": gc,
        "gc_clamp": s[-5:], "gc_clamp_count": sum(1 for c in s[-5:]
                                                  if c in "GC"),
    }

    if tm is None:
        notlar.append(_note("uyari", "Tm hesaplanamadı",
                            "Dizide standart olmayan baz olabilir."))
    elif not 55 <= tm <= 65:
        notlar.append(_note(
            "uyari", f"Tm {tm} °C",
            "55-65 °C dışında. Çok düşükse özgül olmayan bağlanma, çok "
            "yüksekse bağlanma sıcaklığını ayarlamak zorlaşır."))
    if not 40 <= gc <= 60:
        notlar.append(_note(
            "uyari", f"GC %{gc}",
            "%40-60 dışında. GC düşükse bağlanma zayıf, yüksekse "
            "özgül olmayan bağlanma ve saç tokası riski artar."))
    if s[-1] not in "GC":
        notlar.append(_note(
            "uyari", "3' ucu G/C değil",
            "Son baz A/T olduğunda uzamanın başlangıcı gevşektir. "
            "Mümkünse primeri bir baz kaydırın."))
    if olcum["gc_clamp_count"] > 3:
        notlar.append(_note(
            "uyari", "3' ucunda GC yoğunluğu yüksek",
            f"Son 5 bazın {olcum['gc_clamp_count']}'i G/C. Bu uç "
            "yanlış yerlere de sıkıca oturur ve yan ürün verir."))

    kendi = primer_dimer(s)
    if kendi["risk"] != "yok":
        notlar.append(_note("uyari", "Kendi dimeri", kendi["note"]))
    toka = primer_hairpin(s)
    if toka["risk"] != "yok":
        notlar.append(_note("uyari", "Saç tokası", toka["note"]))

    oz = None
    if template:
        oz = specificity(s, template, circular=circular)
        notlar.extend(oz["notes"])

    if not notlar:
        notlar.append(_note("bilgi", "Belirgin bir sorun yok",
                            "Ölçülen değerlerin hepsi olağan aralıkta."))
    return {"metrics": olcum, "notes": notlar,
            "specificity": oz, "dimer": kendi, "hairpin": toka}
