"""Denetim izi: bir tasarımın nereden geldiğini kanıtlayan dosya.

Bir tasarımın kendisi yeterli değildir. Altı ay sonra -ya da bir
denetçi karşısında- sorulacak soru şudur: *bu dizi hangi veriyle, hangi
sürümle, hangi denetim kurallarından geçerek üretildi ve hangi
kurallardan GEÇEMEDİ?* Bu modül o soruya tek dosyayla yanıt verir.

Üç şeyi bir arada tutar:

1. **Kimlik** — uygulama sürümü, veri parmak izi (izlenen bütün veri
   dosyalarının SHA-256 özetinden türetilir) ve dosyaların o anda
   değişmemiş olduğunun doğrulaması. Veri dosyası değiştiyse
   parmak izi değişir, yani eski bir raporun hangi veriyle
   üretildiği belirsiz kalmaz.
2. **Tasarım** — parça listesi, her parçanın dizi durumu, uzunluğu
   ve aksesyonu; yapının kendi SHA-256'sı.
3. **Denetim sonucu** — hangi kural kimlikleri koştu, ne buldu,
   hangileri HATA verdi. Hata veren kural **gizlenmez**: uyum
   dosyasının değeri, geçenleri değil geçemeyenleri de yazmasından
   gelir. Yalnızca başarıları yazan bir rapor kanıt değildir.

Dosya iki biçimde verilir: `to_json` makine okuyacak sürüm (boru
hattında ve komut satırında kullanılır), `to_markdown` insan okuyacak
sürüm. İkisi de aynı sözlükten üretilir, yani ayrışamazlar.

**Bu dosya bir uygunluk belgesi değildir.** Hiçbir düzenleyici kurum
tarafından tanınmaz ve öyle sunulamaz; yaptığı iş, tasarımın
izlenebilir olmasıdır. `core/notice.py` içindeki kullanım amacı
beyanı her çıktıya eklenir ki bu ayrım çıktının kendisinde dursun.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
from dataclasses import asdict, dataclass, field
from typing import Any

from core import db, notice, provenance


@dataclass
class AuditRecord:
    """Tek bir tasarımın denetim izi."""

    design_name: str
    created: str
    app_version: str
    data_fingerprint: str
    data_intact: bool
    data_changed: list[str] = field(default_factory=list)
    host: str = ""
    purpose: str = ""
    length: int = 0
    sequence_sha256: str = ""
    parts: list[dict[str, Any]] = field(default_factory=list)
    rules_run: list[str] = field(default_factory=list)
    findings: list[dict[str, str]] = field(default_factory=list)
    notice_text: str = ""

    @property
    def errors(self) -> list[dict[str, str]]:
        return [f for f in self.findings if f["severity"] == "hata"]

    @property
    def warnings(self) -> list[dict[str, str]]:
        return [f for f in self.findings if f["severity"] == "uyari"]

    @property
    def verdict(self) -> str:
        """Tek satırlık sonuç - hatayı gizlemeden.

        "Geçti" demek için hata olmaması YETERLİ sayılır, ama uyarı
        varsa bu da yazılır: uyarıyı yutan bir sonuç satırı, raporu
        okuyan kişiyi yanlış güvene sürükler.
        """
        if self.errors:
            return (f"HATA: {len(self.errors)} kural hata verdi"
                    + (f", {len(self.warnings)} uyarı"
                       if self.warnings else ""))
        if self.warnings:
            return f"Hata yok, {len(self.warnings)} uyarı var"
        return "Hata ve uyarı yok"


def _part_row(p: dict[str, Any]) -> dict[str, Any]:
    """Denetim izine girecek parça satırı.

    Dizinin TAMAMI yazılmaz (dosya şişer), ama dizi durumu, uzunluğu
    ve aksesyonu yazılır: bir parçanın dizisi "full" mü yoksa kısmi mi
    olduğunu sonradan bilmek, o parçayla kurulmuş yapıya ne kadar
    güvenileceğini belirler.
    """
    return {
        "part_id": p.get("part_id", ""),
        "name": p.get("name", ""),
        "category": p.get("category", ""),
        "length": p.get("length"),
        "seq_status": p.get("seq_status", ""),
        "accession": p.get("accession", ""),
        "source": p.get("source", ""),
        "pmid": p.get("pmid", ""),
    }


def record(design_name: str, part_ids: list[str], sequence: str = "",
           host: str = "", purpose: str = "ifade",
           findings: list[Any] | None = None) -> AuditRecord:
    """Tasarımdan denetim izi kaydı üret.

    `findings` verilmezse denetim BURADA koşturulur; verilirse olduğu
    gibi alınır. İkinci yol, arayüzün kullanıcıya gösterdiği
    bulguların aynısının dosyaya yazılmasını sağlar - ekranda bir
    şey, dosyada başka bir şey görünmesi en kötü ihtimaldir.
    """
    from core.assembly import assemble
    from core.validate import RULE_IDS, validate

    parcalar = db.get_parts(part_ids)
    vec = assemble(design_name, parcalar)
    dizi = (sequence or vec.sequence or "").upper()

    if findings is None:
        findings = validate(vec, host=db.get_host(host) if host else None,
                            purpose=purpose)

    dogrulama = provenance.verify()
    # Katalogun kendisi, `validate()` içinde tanımlı bütün kural
    # kimlikleri. Bir testle kaynağa bağlıdır, yani "şu kurallar
    # koştu" satırı zamanla sessizce yanlış olmaz.
    kurallar = list(RULE_IDS)

    return AuditRecord(
        design_name=design_name,
        created=_dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        # `stamp()` makine için JSON döner; burada insan okuyacak
        # sürümü yazıyoruz, çünkü bu satır raporun başında duruyor.
        # Makine okuyacak hali zaten `data_fingerprint` alanında var.
        app_version=f"VectorLab {provenance.VERSION}",
        data_fingerprint=provenance.data_fingerprint(),
        data_intact=not (dogrulama.get("degisen") or dogrulama.get("eksik")),
        data_changed=list(dogrulama.get("degisen") or [])
        + [f"{x} (eksik)" for x in dogrulama.get("eksik") or []],
        host=host,
        purpose=purpose,
        length=len(dizi),
        sequence_sha256=hashlib.sha256(dizi.encode()).hexdigest(),
        parts=[_part_row(p) for p in parcalar],
        rules_run=kurallar,
        findings=[{"severity": f.severity, "rule_id": f.rule_id,
                   "title": f.title, "detail": f.detail,
                   "fix": getattr(f, "fix", "")} for f in findings],
        notice_text=notice.EXPORT if hasattr(notice, "EXPORT")
        else notice.FOOTER,
    )


def to_json(r: AuditRecord, indent: int = 1) -> str:
    g = asdict(r)
    # Türetilmiş alanlar da dosyaya girer: dosyayı okuyan bir araç
    # sonucu yeniden hesaplamak zorunda kalmasın.
    g["verdict"] = r.verdict
    g["error_count"] = len(r.errors)
    g["warning_count"] = len(r.warnings)
    return json.dumps(g, ensure_ascii=False, indent=indent)


def to_markdown(r: AuditRecord) -> str:
    s = [f"# Denetim izi — {r.design_name}", "",
         f"**Sonuç:** {r.verdict}", ""]
    s += ["| Alan | Değer |", "|---|---|",
          f"| Oluşturma | {r.created} |",
          f"| Uygulama | {r.app_version} |",
          f"| Veri parmak izi | `{r.data_fingerprint}` |",
          f"| Veri bütünlüğü | "
          + ("değişmemiş" if r.data_intact
             else f"**DEĞİŞMİŞ**: {', '.join(r.data_changed)}") + " |",
          f"| Konak | {r.host or '—'} |",
          f"| Amaç | {r.purpose} |",
          f"| Uzunluk | {r.length} bp |",
          f"| Dizi SHA-256 | `{r.sequence_sha256}` |",
          f"| Koşan kural sayısı | {len(r.rules_run)} |", ""]

    s += ["## Parçalar", "",
          "| Parça | Kategori | bp | Dizi durumu | Kaynak kaydı |",
          "|---|---|---|---|---|"]
    for p in r.parts:
        s.append(f"| {p['part_id']} | {p['category']} | "
                 f"{p['length'] if p['length'] is not None else '—'} | "
                 f"{p['seq_status'] or '—'} | {p['accession'] or '—'} |")
    s.append("")

    s += ["## Denetim bulguları", ""]
    if not r.findings:
        s.append("Hiçbir kural bulgu üretmedi.")
    else:
        s += ["| Önem | Kural | Bulgu |", "|---|---|---|"]
        for f in r.findings:
            s.append(f"| {f['severity']} | {f['rule_id']} | "
                     f"{f['title']} |")
    s.append("")
    s += ["## Koşan kurallar", "", ", ".join(r.rules_run) or "—", ""]
    if r.notice_text:
        s += ["---", "", r.notice_text, ""]
    return "\n".join(s) + "\n"
