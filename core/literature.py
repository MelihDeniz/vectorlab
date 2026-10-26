"""Literatür modülü.

Dört işlev:
  1. ara()                  - Europe PMC üzerinden literatür taraması
  2. bibliography()         - kaynakça / atıf üretimi (BibTeX, numaralı liste)
  3. check_nomenclature()   - metin ve isimlendirme düzeltmesi
  4. extract_candidates()   - makale metninden parça adayı çıkarma

Europe PMC REST API anahtar gerektirmez ve PubMed kayıtlarını da kapsar.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

import requests

EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
TIMEOUT = 20

#: Bakteriyel vektör bağlamına daraltan ek sorgu terimleri.
BACTERIAL_FILTER = (
    '(("Escherichia coli" OR plasmid OR vector OR "expression system" '
    'OR cloning OR bacterial) NOT (mammalian OR "cell line" OR mice OR '
    'plant OR yeast OR baculovirus OR insect))'
)


# ---------------------------------------------------------------------------
# 1. Arama
# ---------------------------------------------------------------------------


def search(
    query: str,
    limit: int = 20,
    bacterial_only: bool = True,
    open_access_only: bool = False,
    sort: str = "relevance",
) -> dict[str, Any]:
    """Europe PMC'de ara. Hata durumunda {'error': ...} döndürür.

    sort: 'relevance' (varsayılan), 'cited' (çok atıf alan önce),
          'recent' (yeni önce).
    """
    if not query.strip():
        return {"results": [], "hit_count": 0}

    q = query.strip()
    if bacterial_only:
        q = f"({q}) AND {BACTERIAL_FILTER}"
    if open_access_only:
        q = f"({q}) AND OPEN_ACCESS:Y"

    params: dict[str, Any] = {
        "query": q,
        "format": "json",
        "resultType": "core",
        "pageSize": min(limit, 100),
    }
    # Europe PMC'de sort parametresini hiç göndermemek alaka sırası demektir;
    # göndermek sıralamayı alakadan koparır, bu yüzden sadece istendiğinde ekle.
    if sort == "cited":
        params["sort"] = "CITED desc"
    elif sort == "recent":
        params["sort"] = "P_PDATE_D desc"

    try:
        r = requests.get(
            EPMC_SEARCH,
            params=params,
            timeout=TIMEOUT,
            headers={"User-Agent": "VectorLab/1.0 (egitim amacli ders projesi)"},
        )
        r.raise_for_status()
        data = r.json()
    except requests.exceptions.Timeout:
        return {"error": "Europe PMC yanıt vermedi (zaman aşımı). Tekrar deneyin."}
    except requests.exceptions.RequestException as e:
        return {"error": f"Bağlantı hatası: {e}"}
    except ValueError:
        return {"error": "Europe PMC beklenmeyen bir yanıt döndürdü."}

    body = data.get("resultList", {}).get("result", [])
    return {
        "results": [_normalize(x) for x in body],
        "hit_count": data.get("hitCount", 0),
        "query_used": q,
    }


def _normalize(x: dict[str, Any]) -> dict[str, Any]:
    authors = x.get("authorString", "") or ""
    return {
        "pmid": x.get("pmid") or x.get("id") or "",
        "pmcid": x.get("pmcid", ""),
        "doi": x.get("doi", ""),
        "title": (x.get("title") or "").rstrip("."),
        "authors": authors,
        "first_author": authors.split(",")[0].strip() if authors else "",
        "journal": (
            x.get("journalInfo", {}).get("journal", {}).get("title", "")
            or x.get("journalTitle", "")
        ),
        "year": str(x.get("pubYear", "")),
        "volume": x.get("journalInfo", {}).get("volume", ""),
        "pages": x.get("pageInfo", ""),
        "abstract": x.get("abstractText", "") or "",
        "cited_by": x.get("citedByCount", 0),
        "is_oa": (x.get("isOpenAccess") == "Y"),
        "source_db": "europepmc",
    }


def url_for(ref: dict[str, Any]) -> str:
    if ref.get("doi"):
        return f"https://doi.org/{ref['doi']}"
    if ref.get("pmid"):
        return f"https://pubmed.ncbi.nlm.nih.gov/{ref['pmid']}/"
    return ""


# ---------------------------------------------------------------------------
# 2. Kaynakça
# ---------------------------------------------------------------------------


def _bibtex_key(ref: dict[str, Any]) -> str:
    first = re.sub(r"[^A-Za-z]", "", (ref.get("first_author") or "anon").split()[0])
    word = ""
    for w in re.findall(r"[A-Za-z]{5,}", ref.get("title", "")):
        word = w.lower()
        break
    return f"{first.lower()}{ref.get('year', '')}{word}" or "ref"


def to_bibtex(refs: list[dict[str, Any]]) -> str:
    out = []
    for r in refs:
        fields = [
            f"  title   = {{{r.get('title', '')}}}",
            f"  author  = {{{(r.get('authors') or '').replace(', ', ' and ')}}}",
            f"  journal = {{{r.get('journal', '')}}}",
            f"  year    = {{{r.get('year', '')}}}",
        ]
        if r.get("doi"):
            fields.append(f"  doi     = {{{r['doi']}}}")
        if r.get("pmid"):
            fields.append(f"  pmid    = {{{r['pmid']}}}")
        out.append("@article{" + _bibtex_key(r) + ",\n" + ",\n".join(fields) + "\n}")
    return "\n\n".join(out)


def to_numbered(refs: list[dict[str, Any]]) -> str:
    """Numaralı kaynakça (Vancouver benzeri)."""
    lines = []
    for i, r in enumerate(refs, 1):
        authors = r.get("authors", "")
        parts = [a.strip() for a in authors.split(",") if a.strip()]
        if len(parts) > 6:
            authors = ", ".join(parts[:6]) + ", et al."
        bits = [authors, r.get("title", "")]
        tail = r.get("journal", "")
        if r.get("year"):
            tail += f". {r['year']}"
        if r.get("volume"):
            tail += f";{r['volume']}"
        if r.get("pages"):
            tail += f":{r['pages']}"
        bits.append(tail)
        if r.get("doi"):
            bits.append(f"doi:{r['doi']}")
        lines.append(f"{i}. " + ". ".join(b for b in bits if b).replace("..", ".") + ".")
    return "\n".join(lines)


def to_markdown_list(refs: list[dict[str, Any]]) -> str:
    out = []
    for r in refs:
        link = url_for(r)
        title = r.get("title", "")
        head = f"[{title}]({link})" if link else title
        out.append(
            f"- {head}  \n  {r.get('first_author', '')} ve ark., "
            f"*{r.get('journal', '')}* {r.get('year', '')}"
            + (f" · {r['cited_by']} atıf" if r.get("cited_by") else "")
        )
    return "\n".join(out)


# ---------------------------------------------------------------------------
# 3. İsimlendirme / metin düzeltme
# ---------------------------------------------------------------------------


@dataclass
class TextIssue:
    kind: str
    found: str
    suggestion: str
    explanation: str
    position: int = 0
    auto_fixable: bool = True


#: Moleküler biyolojide italik yazılan bakteriyel gen adları.
GENE_NAMES = [
    "lacZ", "lacI", "lacY", "lacA", "bla", "recA", "recBCD", "endA", "araC",
    "araBAD", "araE", "rhaBAD", "rhaR", "rhaS", "tetR", "tetA", "cat", "aadA",
    "nptII", "gyrA", "gyrB", "ompT", "lon", "trxA", "malE", "dam", "dcm",
    "pir", "repA", "rop", "rom", "ccdB", "hsdR", "mcrA", "mcrBC", "mrr",
    "gal", "thi", "relA", "deoR", "glnV", "sulA", "dnaA", "oriC", "rrnB",
]

#: Bakteri tür adları - ilk kullanımda açık, sonra kısaltılmış, hep italik.
SPECIES = [
    "Escherichia coli", "Bacillus subtilis", "Pseudomonas putida",
    "Salmonella typhimurium", "Vibrio natriegens", "Corynebacterium glutamicum",
    "Lactococcus lactis", "Streptomyces coelicolor",
]

#: Restriksiyon enzimleri: güncel kurala göre TAMAMEN düz yazılır (italik yok).
ENZYME_NAMES = [
    "EcoRI", "EcoRV", "BamHI", "HindIII", "XhoI", "NotI", "NcoI", "NdeI",
    "XbaI", "SpeI", "PstI", "SalI", "SacI", "KpnI", "SmaI", "XmaI", "BglII",
    "ClaI", "NheI", "PvuII", "SphI", "BsaI", "BsmBI", "BbsI", "DpnI",
]


def check_nomenclature(text: str) -> list[TextIssue]:
    """Metni moleküler biyoloji yazım kurallarına göre denetle."""
    issues: list[TextIssue] = []
    if not text.strip():
        return issues

    def add(kind, found, sugg, expl, pos, auto=True):
        issues.append(TextIssue(kind, found, sugg, expl, pos, auto))

    # --- Tür adları italik değil
    for sp in SPECIES:
        for m in re.finditer(rf"(?<![*\w]){re.escape(sp)}(?![*\w])", text):
            add("tür adı", sp, f"*{sp}*",
                "Tür adları (cins + tür) italik yazılır.", m.start())
        abbr = sp.split()[0][0] + ". " + sp.split()[1]
        for m in re.finditer(rf"(?<![*\w]){re.escape(abbr)}(?![*\w])", text):
            add("tür adı", abbr, f"*{abbr}*",
                "Kısaltılmış tür adı da italik yazılır.", m.start())

    # --- Gen adları italik değil
    for g in GENE_NAMES:
        for m in re.finditer(rf"(?<![*\w\-]){re.escape(g)}(?![*\w\-])", text):
            add("gen adı", g, f"*{g}*",
                f"Gen adları italiktir; aynı adın protein karşılığı "
                f"({g[0].upper() + g[1:]}) düz ve büyük harfle başlar.",
                m.start())

    # --- Enzim adları yanlışlıkla italik
    for e in ENZYME_NAMES:
        for m in re.finditer(rf"\*{re.escape(e)}\*", text):
            add("enzim adı", f"*{e}*", e,
                "Güncel kurala göre restriksiyon enzimi adları tamamen "
                "düz yazılır; eski kaynaklardaki italik kullanım terk "
                "edilmiştir.", m.start())

    # --- Sıcaklık biçimi
    for m in re.finditer(r"(\d+)\s*(?:°\s*C|C\b(?!\w))", text):
        raw = m.group(0)
        fixed = f"{m.group(1)} °C"
        if raw != fixed:
            add("birim", raw, fixed,
                "SI yazımında sayı ile birim arasında boşluk bırakılır: 37 °C.",
                m.start())

    # --- Derişim biçimi
    for m in re.finditer(r"(\d+(?:[.,]\d+)?)\s*(ug|µg|mg|ng)\s*/\s*(ml|mL|l|L)",
                         text):
        amount, unit, vol = m.group(1), m.group(2), m.group(3)
        unit_fixed = "µg" if unit in ("ug", "µg") else unit
        vol_fixed = "mL" if vol.lower() == "ml" else "L"
        fixed = f"{amount} {unit_fixed}/{vol_fixed}"
        if m.group(0) != fixed:
            add("birim", m.group(0), fixed,
                "Hacim birimi büyük harfle (mL, L), mikro işareti µ ile "
                "yazılır ve sayıdan boşlukla ayrılır.", m.start())

    # --- bp / kb boşluğu
    for m in re.finditer(r"(\d+)\s*(bp|kb|Kb|KB)\b", text):
        fixed = f"{m.group(1)} {'bp' if m.group(2).lower() == 'bp' else 'kb'}"
        if m.group(0) != fixed:
            add("birim", m.group(0), fixed,
                "Sayı ile birim arasında boşluk olur; kilobaz 'kb' şeklinde "
                "yazılır.", m.start())

    # --- Suş adları
    # Desenler sondaki boşluğu YUTMAMALI; yoksa düzeltme kelimeleri birleştirir
    # ("BL21(DE3)susunda").
    strain_fixes = {
        r"\bDH5\s*-?\s*(?:a|alfa|alpha)\b": "DH5α",
        r"\bBL21\s*\(?\s*DE3\s*\)?": "BL21(DE3)",
        r"\bTop\s*10\b": "TOP10",
        r"\bJM\s*109\b": "JM109",
    }
    for pat, fixed in strain_fixes.items():
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            found = m.group(0).rstrip()
            if found != fixed:
                add("suş adı", found, fixed,
                    "Suş adları standart yazımla verilir; DE3 lizojeni "
                    "parantez içinde gösterilir.", m.start())

    # --- Plazmit adları
    for m in re.finditer(r"\b(P)(UC|ET|BR|ACYC|GEX|BAD)([0-9A-Za-z\-]*)", text):
        fixed = "p" + m.group(2) + m.group(3)
        add("plazmit adı", m.group(0), fixed,
            "Plazmit adları küçük 'p' ile başlar (pUC19, pET-28a).",
            m.start())

    # --- Sık yapılan yazım hataları
    typos = {
        r"\bEscherichia\s+Coli\b": "Escherichia coli",
        r"\bE\.?\s*Coli\b": "E. coli",
        r"\bshine[\s-]*dalgarno\b": "Shine-Dalgarno",
        r"\bplazmid\b": "plazmit",
        r"\bnukleotid\b": "nükleotid",
        r"\bkodlayici\b": "kodlayıcı",
        r"\bX\s*gal\b": "X-gal",
        r"\bkanamisin\s*direnci\s*geni\b": "kanamisin direnç geni",
    }
    for pat, fixed in typos.items():
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            if m.group(0) != fixed:
                add("yazım", m.group(0), fixed,
                    "Standart yazım.", m.start())

    # --- Türkçe karşılığı tercih edilen terimler
    # Bunlar OTOMATİK DÜZELTİLMEZ. Türkçe eklerde ünlü uyumu ve kaynaştırma
    # harfi vardır; kökü değiştirip eki olduğu gibi bırakmak bozuk sonuç
    # üretir ("ekspresyonu" -> "ifadeu"). Doğru çekimi yazar seçer.
    stem_fixes = {
        r"\bdizayn[a-zçğıöşüi]*\b": "tasarım",
        r"\bkonstrakt[a-zçğıöşüi]*\b": "yapı",
        r"\bekspresyon[a-zçğıöşüi]*\b": "ifade",
        r"\bamplifiye\b": "çoğaltılmış",
    }
    for pat, stem in stem_fixes.items():
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            add("terim", m.group(0), stem,
                f"Türkçe karşılığı tercih edilir ('{stem}'). Eki cümleye "
                "göre siz uyarlayın - otomatik uygulanmaz.",
                m.start(), auto=False)

    # --- İlk kullanımda kısaltma tanımı
    if re.search(r"\bE\.\s*coli\b", text) and not re.search(
        r"\bEscherichia\s+coli\b", text
    ):
        add("kısaltma", "E. coli", "Escherichia coli",
            "Kısaltılmış tür adı ilk kullanımda açık yazılır; sonraki "
            "kullanımlarda kısaltılabilir.", 0, auto=False)

    # Aynı bulgu tekrarlarını tek kayda indir
    seen, uniq = set(), []
    for i in issues:
        key = (i.kind, i.found, i.suggestion)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(i)
    return sorted(uniq, key=lambda i: i.position)


def apply_fixes(text: str, issues: list[TextIssue]) -> str:
    """Otomatik düzeltilebilir bulguları metne uygula."""
    out = text
    for i in sorted(issues, key=lambda x: -len(x.found)):
        if not i.auto_fixable:
            continue
        if i.kind in ("gen adı", "tür adı"):
            # Zaten italik olanlara dokunma
            out = re.sub(
                rf"(?<![*\w\-]){re.escape(i.found)}(?![*\w\-])",
                i.suggestion.replace("\\", r"\\"),
                out,
            )
        else:
            out = out.replace(i.found, i.suggestion)
    return out


# ---------------------------------------------------------------------------
# 4. Makale metninden parça adayı çıkarma
# ---------------------------------------------------------------------------

_ACCESSION = re.compile(r"\b(?:[A-Z]{1,2}\d{5,6}|NC_\d{6}|NP_\d{6,9}|[A-Z]{2}\d{6})\b")
_PLASMID = re.compile(r"\bp[A-Z][A-Za-z0-9]{1,12}(?:-[A-Za-z0-9]{1,6})?\b")
#: Promotör adlandırmasının gerçekten kullanılan biçimleri.
#: Önceki hali "p" ile başlayan her kelimeyi (phage, purified) yakalıyordu.
_PROMOTER = re.compile(
    r"(?:"
    r"\bP(?:BAD|bad|lac|tac|trc|tet|rha|L|R|araBAD|rhaBAD|T7)\b"          # Plac, PBAD
    r"|\bP[_\-](?:BAD|lac|tac|trc|tet|rha|T7|L|R)\b"                      # P_BAD, P-lac
    r"|\b(?:T7lac|T7)\s+promoter\b"                                        # T7 promoter
    r"|\b(?:lac|tac|trc|tet|ara|araBAD|rhaBAD|rrnB|lacUV5|J231\d\d)\s+promoter\b"
    r"|\bpromoter\s+(?:P[A-Za-z0-9]{1,8})\b"
    r")"
)
_ANTIBIOTIC = re.compile(
    r"\b(ampicillin|kanamycin|chloramphenicol|tetracycline|spectinomycin|"
    r"carbenicillin)\b(?:[^.]{0,30}?(\d+)\s*(?:µ|u)?g\s*/\s*m[lL])?",
    re.IGNORECASE,
)
# Sondaki \b kaldırıldı: ")" ile biten eşleşmelerde sınır sağlanamadığı için
# regex geri izleyip kapanış parantezini düşürüyordu ("BL21(DE3").
_STRAIN = re.compile(
    r"\b(?:BL21\s*\(?DE3\)?(?:\s*pLysS)?|DH5[αa]\b|TOP10\b|JM10[59]\b|"
    r"Rosetta\s*\(?DE3\)?|SHuffle(?:\s*T7)?\b|C4[13]\s*\(?DE3\)?)"
)


@dataclass
class Candidate:
    kind: str
    value: str
    context: str
    pmid: str = ""
    extras: dict[str, Any] = field(default_factory=dict)


def extract_candidates(text: str, pmid: str = "") -> list[Candidate]:
    """Abstract/tam metinden katalog adayı bilgileri çıkar.

    Çıkan sonuçlar ÖNERİDİR - kullanıcı onaylamadan kataloğa yazılmaz.
    """
    out: list[Candidate] = []
    if not text:
        return out

    def ctx(m) -> str:
        a, b = max(0, m.start() - 70), min(len(text), m.end() + 70)
        return "…" + text[a:b].replace("\n", " ").strip() + "…"

    seen: set[tuple[str, str]] = set()

    def push(kind, value, m, **extras):
        key = (kind, value.lower())
        if key in seen:
            return
        seen.add(key)
        out.append(Candidate(kind, value, ctx(m), pmid, extras))

    for m in _PLASMID.finditer(text):
        push("plazmit", m.group(0), m)
    for m in _STRAIN.finditer(text):
        push("suş", m.group(0), m)
    for m in _ACCESSION.finditer(text):
        push("erişim no", m.group(0), m)
    for m in _PROMOTER.finditer(text):
        v = m.group(0).strip()
        if len(v) > 2:
            push("promotör", v, m)
    for m in _ANTIBIOTIC.finditer(text):
        push("antibiyotik", m.group(1).lower(), m,
             konsantrasyon=(m.group(2) + " µg/mL") if m.group(2) else "")

    return out


# ---------------------------------------------------------------------------
# Tasarım doğrulama için literatür sorgusu üretimi
# ---------------------------------------------------------------------------


def queries_for_parts(parts: list[dict[str, Any]]) -> list[tuple[str, str]]:
    """Her parça için anlamlı bir literatür sorgusu üret.

    (parça adı, sorgu) çiftleri döndürür.
    """
    out = []
    for p in parts:
        base = p["name"].split("(")[0].strip()
        cat = p["category"]
        hint = {
            "promoter": "promoter induction expression Escherichia coli",
            "ori": "origin of replication copy number plasmid",
            "marker": "antibiotic resistance selection plasmid",
            "tag": "affinity tag purification fusion protein",
            "rbs": "ribosome binding site translation initiation",
            "terminator": "transcription terminator",
            "reporter": "reporter gene assay bacteria",
            "regulator": "repressor regulation inducible",
        }.get(cat, "plasmid vector cloning")
        out.append((p["name"], f"{base} {hint}"))
    return out
