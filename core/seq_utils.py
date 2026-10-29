"""Dizi motoru: temizleme, enzim haritası, ORF, primer, anotasyon.

Ağır işin tamamını Biopython yapar; buradaki kod onu bakteriyel klonlama
iş akışına uygun hale getirir.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from Bio.Data import CodonTable
from Bio.Restriction import CommOnly, RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

# ---------------------------------------------------------------------------
# Temel temizleme / ölçüm
# ---------------------------------------------------------------------------

_IUPAC_DNA = set("ACGTNRYSWKMBDHV")
_AMBIG = {
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def clean_seq(raw: str) -> str:
    """FASTA başlıklarını, rakamları, boşlukları at; büyük harfe çevir."""
    if not raw:
        return ""
    lines = [ln for ln in raw.splitlines() if not ln.strip().startswith(">")]
    text = "".join(lines).upper()
    text = re.sub(r"[^A-Z]", "", text)
    return "".join(ch for ch in text if ch in _IUPAC_DNA)


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return round(gc_fraction(seq, ambiguous="ignore") * 100, 2)


def tm_nn(seq: str) -> float | None:
    """En yakın komşu modeliyle erime sıcaklığı (primer için)."""
    if len(seq) < 8:
        return None
    try:
        return round(float(mt.Tm_NN(Seq(seq), nn_table=mt.DNA_NN4)), 1)
    except (ValueError, IndexError):
        return None


def molecular_weight_ds(seq: str) -> float:
    """Çift zincirli dsDNA kütlesi (kDa), yaklaşık."""
    return round(len(seq) * 650 / 1000, 1)


def seq_stats(seq: str, circular: bool = True) -> dict[str, Any]:
    return {
        "length": len(seq),
        "gc": gc_percent(seq),
        "topology": "dairesel" if circular else "doğrusal",
        "mw_kda": molecular_weight_ds(seq),
        "ambiguous": sum(1 for c in seq if c in _AMBIG),
    }


# ---------------------------------------------------------------------------
# Restriksiyon analizi
# ---------------------------------------------------------------------------

#: Ders ve rutin klonlamada en çok kullanılan enzimler.
COMMON_ENZYMES = [
    "EcoRI", "EcoRV", "BamHI", "HindIII", "XhoI", "NotI", "NcoI", "NdeI",
    "XbaI", "SpeI", "PstI", "SalI", "SacI", "KpnI", "SmaI", "XmaI", "ApaI",
    "BglII", "ClaI", "HpaI", "MluI", "NheI", "PvuI", "PvuII", "ScaI", "SphI",
    "StuI", "AatII", "AflII", "AgeI", "AscI", "AvrII", "DraI", "FseI", "PacI",
    "PmeI", "SgrAI", "SnaBI", "SwaI",
]

#: Golden Gate / Tip IIS enzimleri (tanıma bölgesi dışında keser).
TYPE_IIS_ENZYMES = ["BsaI", "BsmBI", "BbsI", "SapI", "AarI", "BtgZI"]

#: Dam/Dcm metilasyonundan etkilenen, E. coli'de sorun çıkaran enzimler.
METHYLATION_BLOCKED = {
    "XbaI": "Dam - TCTAGA bölgesi 5'-GATC ile komşuysa bloke olur",
    "ClaI": "Dam - ATCGAT bölgesi bağlama göre bloke olabilir",
    "MboI": "Dam - GATC her zaman metillenir, E. coli DNA'sını kesmez",
    "BclI": "Dam - TGATCA her zaman bloke",
    "StuI": "Dcm - bağlama göre bloke olabilir",
    "AvaII": "Dcm - etkilenebilir",
    "SexAI": "Dcm - etkilenebilir",
}


def _batch(names: list[str] | None) -> RestrictionBatch:
    if not names:
        return RestrictionBatch(COMMON_ENZYMES)
    return RestrictionBatch(names)


def restriction_sites(
    seq: str, enzymes: list[str] | None = None, circular: bool = True
) -> dict[str, list[int]]:
    """Enzim -> kesim pozisyonları (1 tabanlı)."""
    if not seq:
        return {}
    bio = Seq(seq)
    out: dict[str, list[int]] = {}
    for enz in _batch(enzymes):
        try:
            hits = enz.search(bio, linear=not circular)
        except (ValueError, TypeError):
            continue
        out[str(enz)] = sorted(hits)
    return out


def classify_cutters(
    seq: str, enzymes: list[str] | None = None, circular: bool = True
) -> dict[str, list[tuple[str, list[int]]]]:
    """Enzimleri kesim sayısına göre grupla - klonlama bölgesi seçimi için."""
    sites = restriction_sites(seq, enzymes, circular)
    groups: dict[str, list[tuple[str, list[int]]]] = {
        "tek": [], "cift": [], "cok": [], "kesmez": [],
    }
    for name, pos in sorted(sites.items()):
        n = len(pos)
        key = "kesmez" if n == 0 else "tek" if n == 1 else "cift" if n == 2 else "cok"
        groups[key].append((name, pos))
    return groups


def digest(
    seq: str, enzyme_names: list[str], circular: bool = True
) -> list[dict[str, Any]]:
    """Sindirim sonucu fragmanlar - jel simülasyonu için."""
    cuts: list[int] = []
    for name in enzyme_names:
        for pos in restriction_sites(seq, [name], circular).get(name, []):
            cuts.append(pos)
    cuts = sorted(set(cuts))
    if not cuts:
        kind = "kesilmedi (dairesel)" if circular else "kesilmedi (doğrusal)"
        return [{"start": 1, "end": len(seq), "length": len(seq), "note": kind}]

    frags = []
    if circular:
        for i, c in enumerate(cuts):
            nxt = cuts[(i + 1) % len(cuts)]
            length = (nxt - c) % len(seq) or len(seq)
            frags.append({"start": c, "end": nxt, "length": length, "note": ""})
    else:
        bounds = [1, *cuts, len(seq)]
        for a, b in zip(bounds, bounds[1:]):
            if b - a > 0:
                frags.append({"start": a, "end": b, "length": b - a, "note": ""})
    return sorted(frags, key=lambda f: -f["length"])


def enzyme_info(name: str) -> dict[str, Any]:
    """Tek bir enzimin özellikleri."""
    try:
        enz = RestrictionBatch([name]).get(name)
    except (ValueError, KeyError):
        return {}
    blunt = "kör uç" if enz.is_blunt() else "çıkıntılı uç"
    end = "5' çıkıntı" if enz.is_5overhang() else (
        "3' çıkıntı" if enz.is_3overhang() else "kör"
    )
    return {
        "name": name,
        "site": str(enz.site),
        "elucidate": enz.elucidate(),
        "cut_type": blunt,
        "overhang": end,
        "size": enz.size,
        "methylation": METHYLATION_BLOCKED.get(name, ""),
    }


def compatible_ends(enz_a: str, enz_b: str) -> bool:
    """İki enzim birbirine ligate edilebilir uç bırakıyor mu?"""
    try:
        a = RestrictionBatch([enz_a]).get(enz_a)
        b = RestrictionBatch([enz_b]).get(enz_b)
    except (ValueError, KeyError):
        return False
    if a.is_blunt() and b.is_blunt():
        return True
    try:
        return a.elucidate() != b.elucidate() and a.ovhgseq == b.ovhgseq
    except AttributeError:
        return False


# ---------------------------------------------------------------------------
# ORF / çeviri
# ---------------------------------------------------------------------------

_TABLE11 = CodonTable.unambiguous_dna_by_id[11]  # Bakteriyel kodon tablosu
STARTS = ("ATG", "GTG", "TTG")  # Bakterilerde alternatif başlama kodonları


@dataclass
class Orf:
    start: int          # 1 tabanlı
    end: int
    strand: int         # +1 / -1
    frame: int
    length_nt: int
    length_aa: int
    start_codon: str
    protein: str
    has_stop: bool

    def label(self) -> str:
        yon = "+" if self.strand > 0 else "-"
        return f"ORF {self.start}..{self.end} ({yon}, {self.length_aa} aa)"


def find_orfs(
    seq: str, min_aa: int = 50, circular: bool = True, alt_starts: bool = False
) -> list[Orf]:
    """Her iki zincirde ORF tara. Dairesel ise sınırı aşanları da yakala."""
    if not seq:
        return []
    starts = STARTS if alt_starts else ("ATG",)
    L = len(seq)
    search_space = seq + seq[: min(L, 3000)] if circular else seq
    stops = {"TAA", "TAG", "TGA"}
    out: list[Orf] = []

    for strand in (1, -1):
        target = search_space if strand > 0 else revcomp(search_space)
        for i in range(len(target) - 2):
            codon = target[i : i + 3]
            if codon not in starts:
                continue
            j = i
            prot_codons = []
            hit_stop = False
            while j + 3 <= len(target):
                c = target[j : j + 3]
                if c in stops:
                    hit_stop = True
                    break
                prot_codons.append(c)
                j += 3
            if len(prot_codons) < min_aa:
                continue
            nt_len = len(prot_codons) * 3
            if nt_len > L:
                continue
            # Çeviri
            dna = "".join(prot_codons)
            try:
                protein = str(Seq(dna).translate(table=11))
            except Exception:
                continue
            protein = "M" + protein[1:] if protein else protein

            if strand > 0:
                s0 = i % L
                e0 = (i + nt_len - 1) % L
            else:
                # Ters zincirdeki i, revcomp uzayındaki konum
                s0 = (len(target) - 1 - i) % L
                e0 = (len(target) - 1 - (i + nt_len - 1)) % L

            out.append(
                Orf(
                    start=s0 + 1,
                    end=e0 + 1,
                    strand=strand,
                    frame=(i % 3) + 1,
                    length_nt=nt_len,
                    length_aa=len(prot_codons),
                    start_codon=codon,
                    protein=protein,
                    has_stop=hit_stop,
                )
            )

    # Aynı bitiş noktasına sahip iç içe ORF'lerden en uzununu tut
    best: dict[tuple[int, int], Orf] = {}
    for o in out:
        key = (o.end, o.strand)
        if key not in best or o.length_aa > best[key].length_aa:
            best[key] = o
    return sorted(best.values(), key=lambda o: -o.length_aa)


def translate_cds(seq: str) -> dict[str, Any]:
    """Bir CDS'i çevir ve yaygın sorunları raporla."""
    s = clean_seq(seq)
    issues = []
    if not s:
        return {"protein": "", "issues": ["Dizi boş."], "aa_len": 0}
    if len(s) % 3:
        issues.append(
            f"Uzunluk 3'ün katı değil ({len(s)} nt, kalan {len(s) % 3}) - "
            "çerçeve kayması var."
        )
    if s[:3] not in STARTS:
        issues.append(f"Başlama kodonu yok (ilk kodon: {s[:3]}).")
    prot = str(Seq(s[: len(s) - len(s) % 3]).translate(table=11))
    internal_stop = prot[:-1].find("*")
    if internal_stop >= 0:
        issues.append(
            f"İç stop kodonu: {internal_stop + 1}. aminoasit konumunda. "
            "Kodon optimizasyonu ya da çerçeve hatası olabilir."
        )
    if not prot.endswith("*"):
        issues.append("Stop kodonu ile bitmiyor - okuma çerçevesi devam ediyor.")
    return {
        "protein": prot,
        "issues": issues,
        "aa_len": len(prot.rstrip("*")),
        "mw_kda": round(len(prot.rstrip("*")) * 0.110, 1),
    }


# ---------------------------------------------------------------------------
# Kodon kullanımı
# ---------------------------------------------------------------------------

#: E. coli'de nadir kodonlar. Rosetta suşunun pRARE plazmiti bunların
#: tRNA'larını sağlar.
RARE_CODONS_ECOLI = {
    "AGG": ("Arg", "çok nadir - tandem AGG AGG neredeyse her zaman sorun"),
    "AGA": ("Arg", "çok nadir"),
    "CGA": ("Arg", "nadir"),
    "CGG": ("Arg", "nadir"),
    "ATA": ("Ile", "nadir"),
    "CTA": ("Leu", "nadir"),
    "CCC": ("Pro", "nadir"),
    "GGA": ("Gly", "nadir"),
    "TCA": ("Ser", "orta nadir"),
}


def codon_analysis(cds: str) -> dict[str, Any]:
    """Nadir kodon yükünü ölç - Rosetta gerekli mi kararı için."""
    s = clean_seq(cds)
    s = s[: len(s) - len(s) % 3]
    codons = [s[i : i + 3] for i in range(0, len(s), 3)]
    total = len(codons)
    if not total:
        return {"total": 0, "rare": [], "rare_pct": 0.0, "tandem": [], "counts": {}}

    counts: dict[str, int] = {}
    rare_positions = []
    for idx, c in enumerate(codons):
        counts[c] = counts.get(c, 0) + 1
        if c in RARE_CODONS_ECOLI:
            rare_positions.append({"pos_aa": idx + 1, "codon": c,
                                   "aa": RARE_CODONS_ECOLI[c][0],
                                   "note": RARE_CODONS_ECOLI[c][1]})

    tandem = []
    for i in range(len(codons) - 1):
        if codons[i] in RARE_CODONS_ECOLI and codons[i + 1] in RARE_CODONS_ECOLI:
            tandem.append({"pos_aa": i + 1,
                           "codons": f"{codons[i]} {codons[i + 1]}"})

    return {
        "total": total,
        "rare": rare_positions,
        "rare_pct": round(100 * len(rare_positions) / total, 1),
        "tandem": tandem,
        "counts": counts,
    }


# ---------------------------------------------------------------------------
# RBS boşluğu
# ---------------------------------------------------------------------------

_SD_MOTIFS = ["AGGAGG", "AAGGAG", "AGGAGA", "GAGGAG", "AGGAG", "AAGGA", "GGAGG"]


def rbs_spacing(seq: str, window: int = 60) -> list[dict[str, Any]]:
    """Shine-Dalgarno motifi ile ATG arası boşluğu ölç.

    Optimum aralık 5-13 nt'dir; bunun dışı çeviri verimini düşürür.
    """
    s = clean_seq(seq)
    results = []
    for m in re.finditer("ATG", s):
        atg = m.start()
        upstream = s[max(0, atg - window) : atg]
        best = None
        for motif in _SD_MOTIFS:
            idx = upstream.rfind(motif)
            if idx < 0:
                continue
            gap = len(upstream) - (idx + len(motif))
            if best is None or gap < best["gap"]:
                best = {"motif": motif, "gap": gap}
        if best is None:
            continue
        gap = best["gap"]
        if 5 <= gap <= 13:
            verdict, ok = "optimum aralıkta", True
        elif gap < 5:
            verdict, ok = "çok kısa - ribozom ATG'ye ulaşamaz", False
        else:
            verdict, ok = "çok uzun - çeviri verimi düşer", False
        results.append({
            "atg_pos": atg + 1, "motif": best["motif"],
            "gap": gap, "verdict": verdict, "ok": ok,
        })
    return results


# ---------------------------------------------------------------------------
# Primer tasarımı
# ---------------------------------------------------------------------------


@dataclass
class Primer:
    name: str
    sequence: str
    tm: float | None
    gc: float
    length: int
    tail: str = ""
    warnings: list[str] = field(default_factory=list)

    @property
    def anneal(self) -> str:
        return self.sequence[len(self.tail):]


def _primer_warnings(core: str, tm: float | None) -> list[str]:
    w = []
    if tm is not None and not 55 <= tm <= 65:
        w.append(f"Tm {tm} C - 55-65 C aralığının dışında")
    gc = gc_percent(core)
    if not 40 <= gc <= 60:
        w.append(f"GC %{gc} - %40-60 aralığının dışında")
    if core[-1] not in "GC":
        w.append("3' ucu G/C değil - uzamada kayma riski")
    if re.search(r"(A{5,}|T{5,}|G{4,}|C{5,})", core):
        w.append("Tek nükleotid tekrarı var - yanlış eşleşme riski")
    if gc_percent(core[-5:]) > 80:
        w.append("3' ucunda GC yoğunluğu çok yüksek - spesifik olmayan bağlanma")
    return w


def _tune_length(seq: str, target_tm: float, lo: int = 18, hi: int = 32) -> str:
    """Hedef Tm'e en yakın uzunluğu seç."""
    best, best_d = seq[:lo], 1e9
    for n in range(lo, min(hi, len(seq)) + 1):
        cand = seq[:n]
        t = tm_nn(cand)
        if t is None:
            continue
        d = abs(t - target_tm)
        if d < best_d:
            best, best_d = cand, d
    return best


def design_primers(
    target: str,
    fwd_enzyme: str | None = None,
    rev_enzyme: str | None = None,
    target_tm: float = 60.0,
    add_spacer: bool = True,
    name_prefix: str = "P",
) -> dict[str, Any]:
    """Hedef bölgeyi çoğaltacak primer çifti tasarla.

    Restriksiyon bölgesi eklenirse enzimin verimli kesmesi için 5' ucuna
    kısa bir koruyucu bazlar bloğu (spacer) konur.
    """
    s = clean_seq(target)
    if len(s) < 40:
        return {"error": "Hedef dizi primer tasarımı için çok kısa (<40 nt)."}

    # Dizisi bilinmeyen (stub) bölgelerde primer tasarlanamaz; kullanıcı
    # bunu "GC %0" uyarısından değil, doğrudan görmelidir.
    n_count = s.count("N")
    if n_count > len(s) * 0.5:
        return {
            "error": (
                f"Hedef bölgenin %{round(100 * n_count / len(s))}'i "
                "belirsiz baz (N) - bu bölge dizisi bilinmeyen bir parçaya "
                "denk geliyor. Primer tasarlamak için önce o parçanın "
                "gerçek dizisini katalogdan doldurun ya da dizinin "
                "bilinen bir bölgesini seçin."
            )
        }

    fwd_core = _tune_length(s, target_tm)
    rev_core = _tune_length(revcomp(s), target_tm)

    def tail_for(enz: str | None) -> str:
        if not enz:
            return ""
        info = enzyme_info(enz)
        site = info.get("site", "")
        if not site:
            return ""
        spacer = "AAA" if add_spacer else ""
        return spacer + site

    f_tail, r_tail = tail_for(fwd_enzyme), tail_for(rev_enzyme)
    fwd = Primer(
        name=f"{name_prefix}_F",
        sequence=f_tail + fwd_core,
        tm=tm_nn(fwd_core),
        gc=gc_percent(fwd_core),
        length=len(f_tail + fwd_core),
        tail=f_tail,
        warnings=_primer_warnings(fwd_core, tm_nn(fwd_core)),
    )
    rev = Primer(
        name=f"{name_prefix}_R",
        sequence=r_tail + rev_core,
        tm=tm_nn(rev_core),
        gc=gc_percent(rev_core),
        length=len(r_tail + rev_core),
        tail=r_tail,
        warnings=_primer_warnings(rev_core, tm_nn(rev_core)),
    )

    pair_warnings = []
    if fwd.tm and rev.tm and abs(fwd.tm - rev.tm) > 5:
        pair_warnings.append(
            f"Primer Tm farkı {abs(fwd.tm - rev.tm):.1f} C - 5 C'yi aşıyor, "
            "bağlanma sıcaklığını ikisine birden ayarlamak zor olur."
        )
    for enz in (fwd_enzyme, rev_enzyme):
        if enz and restriction_sites(s, [enz], circular=False).get(enz):
            pair_warnings.append(
                f"{enz} hedef dizinin İÇİNDE de kesiyor - bu enzimle "
                "klonlama insertinizi parçalar."
            )
    if fwd_enzyme and rev_enzyme and fwd_enzyme == rev_enzyme:
        pair_warnings.append(
            "İki uçta aynı enzim - yönlü klonlama yapılamaz, insert iki "
            "yönde de girer."
        )

    # Dimer ve sac tokasi: Tm/GC kadar belirleyici, cunku primerler
    # birbirini kalip olarak kullanirsa urun hic cikmaz.
    for primer, etiket in ((fwd, "İleri"), (rev, "Geri")):
        kendi = primer_dimer(primer.sequence)
        if kendi["risk"] != "yok":
            primer.warnings.append(f"{etiket} primer {kendi['note']}")
        toka = primer_hairpin(primer.sequence)
        if toka["risk"] != "yok":
            primer.warnings.append(f"{etiket} primer: {toka['note']}")
    capraz = primer_dimer(fwd.sequence, rev.sequence)
    if capraz["risk"] != "yok":
        pair_warnings.append(capraz["note"])

    anneal_tm = min(t for t in (fwd.tm, rev.tm) if t) - 5 if (fwd.tm and rev.tm) else None
    return {
        "forward": fwd,
        "reverse": rev,
        "product_length": len(s),
        "annealing_temp": round(anneal_tm, 1) if anneal_tm else None,
        "pair_warnings": pair_warnings,
        "dimer": capraz,
    }


# ---------------------------------------------------------------------------
# Anotasyon (parça kataloğuna karşı eşleştirme)
# ---------------------------------------------------------------------------


@dataclass
class Feature:
    start: int       # 1 tabanlı
    end: int
    strand: int
    name: str
    category: str
    part_id: str = ""
    identity: float = 100.0

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def annotate(
    seq: str, parts: list[dict[str, Any]], circular: bool = True, min_len: int = 10
) -> list[Feature]:
    """Dizide katalogdaki parçaları ara (iki zincir, dairesel sarma dahil).

    Yalnızca `seq_status == 'full'` olan parçalar aranabilir; stub parçalar
    dizi taşımadığı için eşleşme yapılamaz.
    """
    s = clean_seq(seq)
    if not s:
        return []
    L = len(s)
    searchable = [
        p for p in parts
        if p.get("seq_status") == "full" and len(p.get("sequence") or "") >= min_len
    ]
    if not searchable:
        return []

    max_part = max(len(p["sequence"]) for p in searchable)
    haystack = s + s[: min(max_part - 1, L)] if circular else s
    rc_haystack = revcomp(haystack)

    found: list[Feature] = []
    for p in searchable:
        needle = p["sequence"]
        n = len(needle)
        # İleri zincir
        for m in re.finditer(f"(?={re.escape(needle)})", haystack):
            st = m.start() % L
            found.append(Feature(
                start=st + 1, end=((st + n - 1) % L) + 1, strand=1,
                name=p["name"], category=p["category"], part_id=p["part_id"],
            ))
        # Ters zincir
        for m in re.finditer(f"(?={re.escape(needle)})", rc_haystack):
            # rc_haystack konumunu orijinal diziye çevir
            end_fwd = len(haystack) - 1 - m.start()
            st = (end_fwd - n + 1) % L
            found.append(Feature(
                start=st + 1, end=((st + n - 1) % L) + 1, strand=-1,
                name=p["name"], category=p["category"], part_id=p["part_id"],
            ))

    # Tekrarları temizle
    uniq: dict[tuple, Feature] = {}
    for f in found:
        uniq[(f.start, f.end, f.strand, f.part_id)] = f
    return sorted(uniq.values(), key=lambda f: (f.start, -f.length))


def annotate_similar(
    seq: str,
    parts: list[dict[str, Any]],
    circular: bool = True,
    min_identity: float = 90.0,
    min_len: int = 40,
) -> list[Feature]:
    """Katalog parçalarını *yaklaşık* eşleşmeyle ara.

    Gerçek plazmitlerde tam eşleşme çoğu zaman yetmez: aynı gen farklı
    vektör soylarında birkaç nükleotid farkla taşınır. Örneğin `bla`
    pBR322 ve pACYC177'de birebir eşleşirken pUC19, pBAD24, pTrc99A ve
    pBluescript'te nokta farklarıyla bulunur - bunlar pUC soyunun bilinen
    varyantlarıdır. Tam eşleşme arayan `annotate` bu bölgeleri hiç
    göstermez ve harita boş görünür.

    Yöntem: parçanın başından, ortasından ve sonundan alınan tohum
    k-mer'lerle aday başlangıç konumları bulunur, sonra parçanın tamamı
    boyunca uyuşmayan baz sayılır. Hizalama yapılmaz, yani yalnızca aynı
    uzunluktaki nokta farkları yakalanır; insersiyon/delesyon taşıyan
    varyantlar bulunamaz.
    """
    s = clean_seq(seq)
    if not s:
        return []
    L = len(s)
    searchable = [
        p for p in parts
        if p.get("seq_status") == "full" and len(p.get("sequence") or "") >= min_len
    ]
    if not searchable:
        return []

    max_part = max(len(p["sequence"]) for p in searchable)
    haystack = s + s[: min(max_part - 1, L)] if circular else s
    strands = ((1, haystack), (-1, revcomp(haystack)))
    seed = 20

    best: dict[tuple[int, int, str], Feature] = {}
    for p in searchable:
        needle = p["sequence"]
        n = len(needle)
        max_mm = int(n * (1.0 - min_identity / 100.0))
        offsets = sorted({0, max(0, n // 2 - seed // 2), n - seed})
        for strand, hay in strands:
            H = len(hay)
            starts: set[int] = set()
            for off in offsets:
                probe = needle[off:off + seed]
                if len(probe) < seed:
                    continue
                at = hay.find(probe)
                while at != -1:
                    st = at - off
                    if 0 <= st <= H - n:
                        starts.add(st)
                    at = hay.find(probe, at + 1)
            for st in starts:
                window = hay[st:st + n]
                mm = 0
                for a, b in zip(window, needle):
                    if a != b:
                        mm += 1
                        if mm > max_mm:
                            break
                if mm > max_mm:
                    continue
                identity = round(100.0 * (n - mm) / n, 1)
                # Ters zincir konumunu özgün diziye çevir
                if strand == 1:
                    real = st % L
                else:
                    real = (H - st - n) % L
                key = (real, strand, p["part_id"])
                prev = best.get(key)
                if prev is None or identity > prev.identity:
                    best[key] = Feature(
                        start=real + 1, end=((real + n - 1) % L) + 1,
                        strand=strand, name=p["name"], category=p["category"],
                        part_id=p["part_id"], identity=identity,
                    )
    return sorted(best.values(), key=lambda f: (f.start, -f.length))


# ---------------------------------------------------------------------------
# Primer dimer ve saç tokası denetimi
# ---------------------------------------------------------------------------
# Tm, GC ve 3' ucu bakmak yetmez: iki primer birbirine ya da kendine
# yapışırsa polimeraz kalıp yerine onları uzatır. Ürün yerine jelde
# 40-80 bp'lik bulanık bir bant görürsünüz. 3' ucundaki eşleşme
# özellikle tehlikelidir, çünkü uzama oradan başlar.

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def _pair_score(a: str, b: str) -> int:
    """Karşılıklı dizilen iki bölgenin bağ puanı (G/C 2, A/T 1)."""
    skor = 0
    for x, y in zip(a, b):
        if _COMPLEMENT.get(x) == y:
            skor += 2 if x in "GC" else 1
    return skor


def _best_duplex(a: str, b: str) -> tuple[int, int, str]:
    """İki primer arasındaki en güçlü hizalanmayı bul.

    `b` ters çevrilerek `a` üzerinde kaydırılır; her kaymada eşleşen
    baz çiftlerinin puanı toplanır. En yüksek puanlı kayma döner.
    Dönen üçlü: (puan, 3' ucundaki ardışık eşleşme, görsel).
    """
    b_rev = b[::-1]
    en_iyi = (0, 0, "")
    for kaydir in range(-(len(b) - 1), len(a)):
        ust_parca, alt_parca = [], []
        for i in range(len(a)):
            j = i - kaydir
            if 0 <= j < len(b_rev):
                ust_parca.append(a[i])
                alt_parca.append(b_rev[j])
        if len(ust_parca) < 3:
            continue
        ust, alt = "".join(ust_parca), "".join(alt_parca)
        skor = _pair_score(ust, alt)
        if skor < 8:
            continue
        # 3' ucundaki ardışık eşleşme: a'nın son bazından geriye doğru
        uc = 0
        if kaydir + len(b) - 1 >= len(a) - 1:
            for x, y in zip(reversed(ust), reversed(alt)):
                if _COMPLEMENT.get(x) == y:
                    uc += 1
                else:
                    break
        if skor > en_iyi[0]:
            cizim = "\n".join([
                f"5'-{ust}-3'",
                "   " + "".join("|" if _COMPLEMENT.get(x) == y else " "
                                for x, y in zip(ust, alt)),
                f"3'-{alt}-5'",
            ])
            en_iyi = (skor, uc, cizim)
    return en_iyi


def primer_dimer(
    primer_a: str, primer_b: str | None = None, threshold: int = 14,
    end_threshold: int = 3,
) -> dict[str, Any]:
    """Primer dimerini (ya da kendi dimerini) ölç ve yorumla.

    `threshold` toplam bağ puanı, `end_threshold` ise 3' ucundaki
    ardışık eşleşme sayısı için sınırdır. 3' ucu eşleşmesi daha azıyla
    bile sorun çıkarır, çünkü polimeraz oradan uzatmaya başlar.
    """
    a = clean_seq(primer_a)
    b = clean_seq(primer_b) if primer_b else a
    kendi = primer_b is None
    if len(a) < 6 or len(b) < 6:
        return {"score": 0, "end_match": 0, "risk": "yok", "diagram": "",
                "note": ""}

    skor, uc, cizim = _best_duplex(a, b)
    tur = "kendi dimeri" if kendi else "çapraz dimer"

    if uc >= end_threshold:
        risk = "yüksek"
        note = (
            f"{tur}: 3' ucunda {uc} ardışık baz eşleşiyor (bağ puanı "
            f"{skor}). Polimeraz uzamayı 3' ucundan başlattığı için bu "
            f"en kötü durumdur - primerler birbirini kalıp olarak "
            f"kullanır ve jelde ürün yerine 40-80 bp'lik bulanık bir "
            f"bant çıkar. Primerin 3' ucundaki birkaç bazı değiştirin "
            f"ya da bağlanma bölgesini birkaç baz kaydırın."
        )
    elif skor >= threshold:
        risk = "orta"
        note = (
            f"{tur}: bağ puanı {skor}, ama eşleşme 3' ucunda değil. "
            f"İç bölgedeki eşleşme uzamayı doğrudan başlatmaz; yine de "
            f"primerin bir kısmını bağlar ve etkin derişimi düşürür. "
            f"PCR zayıf çalışırsa bakılacak yerlerden biridir."
        )
    else:
        risk = "yok"
        note = ""
    return {"score": skor, "end_match": uc, "risk": risk,
            "diagram": cizim if risk != "yok" else "", "note": note,
            "kind": tur}


def primer_hairpin(primer: str, threshold: int = 12) -> dict[str, Any]:
    """Primerin kendi üzerine katlanıp katlanmadığına bak.

    Saç tokası yapan primer kalıba bağlanamaz; özellikle 3' ucu
    yapının içinde kalıyorsa uzama hiç başlamaz.
    """
    s = clean_seq(primer)
    n = len(s)
    en_iyi = (0, "", 0)
    for i in range(n):
        for j in range(i + 4, n):
            uzunluk = 0
            while (i - uzunluk >= 0 and j + uzunluk < n
                   and _COMPLEMENT.get(s[i - uzunluk]) == s[j + uzunluk]):
                uzunluk += 1
            if uzunluk < 3:
                continue
            sap_a = s[i - uzunluk + 1:i + 1]
            sap_b = s[j:j + uzunluk][::-1]
            skor = _pair_score(sap_a, sap_b)
            if skor > en_iyi[0]:
                uc_icinde = (j + uzunluk - 1) >= n - 3
                en_iyi = (skor, f"{sap_a} / {sap_b}", uc_icinde)
    skor, cizim, uc_icinde = en_iyi
    if skor < threshold:
        return {"score": skor, "risk": "yok", "note": ""}
    return {
        "score": skor,
        "risk": "yüksek" if uc_icinde else "orta",
        "note": (
            f"Primer kendi üzerine katlanıyor (bağ puanı {skor}, sap "
            f"{cizim})."
            + (" 3' ucu yapının içinde kalıyor: uzama hiç başlamayabilir."
               if uc_icinde else
               " 3' ucu serbest, yine de bağlanma verimi düşer.")
        ),
    }
