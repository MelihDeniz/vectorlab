"""VectorLab komut satırı: toplu iş ve boru hattı kullanımı.

Arayüz tek tasarım için iyi; bir şirketin ya da bir laboratuvarın
asıl sorusu şu: *"elimdeki 400 plazmiti denetim motorundan nasıl
geçiririm?"* Tek tek tıklayarak geçirilmez. Bu modül aynı motorları
komut satırından ve JSON çıktısıyla erişilebilir kılar, böylece
VectorLab bir uygulama olmaktan çıkıp bir boru hattı parçası olur.

Alt komutlar:

    vectorlab-cli denetle  DOSYA...   GenBank/FASTA dosyalarını denetle
    vectorlab-cli primer   DOSYA      koloni PCR + dizileme çapaları
    vectorlab-cli siparis  DOSYA      malzeme listesi ve takvim
    vectorlab-cli tasi     DOSYA KONAK  kaseti başka konağa uyarla
    vectorlab-cli iz       DOSYA      denetim izi (uyum dosyası)
    vectorlab-cli kimlik              sürüm ve veri parmak izi

Her komut `--json` ile makine okunur çıktı verir; `--json` olmadan
insan okunur metin yazar. Çıkış kodu sözleşmesi boru hattında
kullanılmak için nettir:

    0  iş bitti, denetimde HATA yok
    1  iş bitti, denetimde en az bir HATA var
    2  dosya okunamadı ya da komut yanlış kullanıldı

Yani `vectorlab-cli denetle *.gb` bir CI adımında doğrudan kullanılır:
hata varsa adım kırmızıya döner. Uyarılar çıkış kodunu değiştirmez,
çünkü uyarı tasarımı geçersiz kılmaz - bu ayrımı kodun kendisine
yazmak, kullanıcıyı `grep` yapmaktan kurtarır.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

EXIT_OK = 0
EXIT_FINDINGS = 1
EXIT_USAGE = 2


# ---------------------------------------------------------------------------
# Girdi okuma
# ---------------------------------------------------------------------------


def read_plasmid(path: Path) -> dict[str, Any]:
    """Dosyayı `parse_upload` ile oku.

    GenBank, FASTA ve düz dizi, `core/inspect_plasmid.py` içindeki
    çözümleyiciyle okunur - arayüzün kullandığı çözümleyicinin
    aynısı. İki ayrı okuyucu olsa biri düzeltildiğinde öbürü sessizce
    eski kalırdı.
    """
    from core import inspect_plasmid as ip

    return ip.parse_upload(path.read_bytes(), path.name)


def read_sequence(path: Path) -> tuple[str, str]:
    """Dosyadan (ad, dizi) oku."""
    veri = read_plasmid(path)
    return (veri.get("name") or path.stem), (veri.get("sequence") or "")


# ---------------------------------------------------------------------------
# denetle
# ---------------------------------------------------------------------------


def cmd_denetle(args: argparse.Namespace) -> int:
    from core import db, inspect_plasmid as ip
    from core.validate import validate

    db.init_db()
    parts = db.list_parts()
    konak = db.get_host(args.konak) if args.konak else None
    if args.konak and konak is None:
        print(f"konak bulunamadi: {args.konak}", file=sys.stderr)
        return EXIT_USAGE

    cikti: list[dict[str, Any]] = []
    hata_var = False
    okunamayan = False

    for ham in args.dosya:
        yol = Path(ham)
        if not yol.exists():
            cikti.append({"dosya": ham, "durum": "okunamadi",
                          "sebep": "dosya yok"})
            okunamayan = True
            continue
        try:
            veri = read_plasmid(yol)
            # Dizi boşsa denetim KOŞTURULMAZ. Boş bir dizi bütün
            # "yok" kurallarını tetikler ve "orijin yok, marker yok,
            # promotör yok" diye beş uyarı basar - oysa sorun
            # tasarımda değil, dosyanın okunamamış olmasındadır.
            # İkisini karıştırmak, toplu işte en yanıltıcı çıktı olur.
            if not (veri.get("sequence") or "").strip():
                cikti.append({"dosya": ham, "durum": "okunamadi",
                              "sebep": "dosyada dizi bulunamadı "
                                       "(boş ya da tanınmayan biçim)"})
                okunamayan = True
                continue
            ins = ip.inspect(veri, parts, min_identity=args.esik)
            vec = ip.to_vector(ins)
            # `annotation_complete=False` ARAYÜZDEKİYLE AYNI ve
            # kasıtlı: dışarıdan gelen bir plazmitin bütün parçaları
            # katalogda tanınmayabilir, o yüzden "yok" diyen kurallar
            # (CORE-01 gibi) hataya değil uyarıya düşürülür. Aksi
            # halde müşterinin her plazmiti sahte hata üretirdi.
            bulgular = validate(vec, host=konak, purpose=args.amac,
                                annotation_complete=False)
        except Exception as exc:                        # pragma: no cover
            cikti.append({"dosya": ham, "durum": "okunamadi",
                          "sebep": str(exc)})
            okunamayan = True
            continue

        hatalar = [b for b in bulgular if _sev(b) == "hata"]
        uyarilar = [b for b in bulgular if _sev(b) == "uyari"]
        if hatalar:
            hata_var = True
        cikti.append({
            "dosya": ham,
            "durum": "denetlendi",
            "ad": ins.name,
            "uzunluk": ins.length,
            "dairesel": ins.circular,
            "tanimlanan_parca": len(ins.catalog_features),
            "kapsam_yuzde": ins.coverage,
            "hata": len(hatalar), "uyari": len(uyarilar),
            "bulgular": [_bulgu(b) for b in bulgular],
        })

    if args.json:
        print(json.dumps({"komut": "denetle", "sonuc": cikti},
                         ensure_ascii=False, indent=1))
    else:
        for r in cikti:
            if r["durum"] == "okunamadi":
                print(f"[OKUNAMADI] {r['dosya']}: {r['sebep']}")
                continue
            isaret = "HATA" if r["hata"] else ("uyari" if r["uyari"]
                                               else "temiz")
            print(f"[{isaret:6}] {r['dosya']}  {r['uzunluk']} bp, "
                  f"{r['tanimlanan_parca']} parça tanındı "
                  f"(%{r['kapsam_yuzde']:.0f} kapsam), "
                  f"{r['hata']} hata / {r['uyari']} uyarı")
            for b in r["bulgular"]:
                if b["severity"] in ("hata", "uyari"):
                    print(f"           {b['severity']:5} {b['rule_id']}  "
                          f"{b['title']}")
        n = sum(1 for r in cikti if r["durum"] == "denetlendi")
        print(f"\n{n} dosya denetlendi, "
              f"{sum(r.get('hata', 0) for r in cikti)} hata.")

    if okunamayan:
        return EXIT_USAGE
    return EXIT_FINDINGS if hata_var else EXIT_OK


def _sev(b: Any) -> str:
    return b.get("severity", "") if isinstance(b, dict) else getattr(
        b, "severity", "")


def _bulgu(b: Any) -> dict[str, str]:
    if isinstance(b, dict):
        return {"severity": b.get("severity", ""),
                "rule_id": b.get("rule_id", ""),
                "title": b.get("title", ""),
                "detail": b.get("detail", "")}
    return {"severity": b.severity, "rule_id": b.rule_id,
            "title": b.title, "detail": b.detail}


# ---------------------------------------------------------------------------
# primer
# ---------------------------------------------------------------------------


def cmd_primer(args: argparse.Namespace) -> int:
    from core import db, primers as pr

    db.init_db()
    yol = Path(args.dosya)
    if not yol.exists():
        print(f"dosya yok: {args.dosya}", file=sys.stderr)
        return EXIT_USAGE
    ad, dizi = read_sequence(yol)
    if not dizi:
        print(f"dizi okunamadi: {args.dosya}", file=sys.stderr)
        return EXIT_USAGE

    k = pr.colony_pcr(dizi, args.bas, args.son, circular=not args.dogrusal)
    capalar = pr.catalog_anchors(dizi, db.list_parts(),
                                 circular=not args.dogrusal)
    sonuc = {
        "komut": "primer", "ad": ad, "uzunluk": len(dizi),
        "koloni_pcr": {
            "forward": {"name": k.forward.name, "sequence": k.forward.sequence,
                        "tm": k.forward.tm, "gc": k.forward.gc},
            "reverse": {"name": k.reverse.name, "sequence": k.reverse.sequence,
                        "tm": k.reverse.tm, "gc": k.reverse.gc},
            "insert_varken": k.with_insert, "bos_vektorde": k.without_insert,
            "fark": k.fark,
        },
        "dizileme_capalari": capalar,
    }
    if args.json:
        print(json.dumps(sonuc, ensure_ascii=False, indent=1))
    else:
        print(f"{ad} — {len(dizi)} bp")
        print(f"\nKoloni PCR ({args.bas}-{args.son} arası insert):")
        for p in (k.forward, k.reverse):
            tm = "—" if p.tm is None else f"{p.tm:.1f} C"
            print(f"  {p.name:10} {p.sequence}  Tm {tm}  "
                  f"GC %{p.gc:.0f}")
        print(f"  insert varken {k.with_insert} bp, boş vektörde "
              f"{k.without_insert} bp, fark {k.fark} bp")
        print(f"\nDizileme çapaları ({len(capalar)}):")
        for c in capalar:
            print(f"  {c['part_id']:32} {c['primer']}  "
                  f"Tm {c['tm']:.1f} C  {c['start']}-{c['end']}")
    return EXIT_OK


# ---------------------------------------------------------------------------
# siparis
# ---------------------------------------------------------------------------


def cmd_siparis(args: argparse.Namespace) -> int:
    from core import db, primers as pr, procure

    db.init_db()
    yol = Path(args.dosya)
    if not yol.exists():
        print(f"dosya yok: {args.dosya}", file=sys.stderr)
        return EXIT_USAGE
    ad, dizi = read_sequence(yol)
    if not dizi:
        print(f"dizi okunamadi: {args.dosya}", file=sys.stderr)
        return EXIT_USAGE

    k = pr.colony_pcr(dizi, args.bas, args.son)
    capalar = pr.catalog_anchors(dizi, db.list_parts())
    bom = procure.build(ad, part_ids=args.parca or [],
                        enzymes=args.enzim or [],
                        primers=[k.forward, k.reverse] + capalar,
                        host_name=args.konak or "")
    if args.csv:
        print(procure.to_csv(bom), end="")
    elif args.oligo_csv:
        print(procure.oligo_order_csv([k.forward, k.reverse] + capalar),
              end="")
    elif args.json:
        print(json.dumps({
            "komut": "siparis", "ad": ad,
            "toplam_gun": bom.total_days, "oligo_baz": bom.oligo_bases,
            "kalemler": [{"tur": i.kind, "ad": i.name,
                          "miktar": i.quantity, "birim": i.unit,
                          "olcek": i.scale, "saflastirma": i.purification,
                          "dizi": i.sequence, "firmalar": i.suppliers,
                          "not": i.note} for i in bom.items],
            "takvim": [{"asama": p.name, "gun": p.days, "not": p.note}
                       for p in bom.phases],
            "varsayimlar": bom.assumptions,
        }, ensure_ascii=False, indent=1))
    else:
        print(procure.to_markdown(bom), end="")
    return EXIT_OK


# ---------------------------------------------------------------------------
# tasi (konak taşıma)
# ---------------------------------------------------------------------------


def cmd_tasi(args: argparse.Namespace) -> int:
    from core import db, porting

    db.init_db()
    if not db.get_host(args.konak):
        print(f"konak bulunamadi: {args.konak}", file=sys.stderr)
        print("katalogdaki konaklar: "
              + ", ".join(sorted(h["host_id"] for h in db.list_hosts())),
              file=sys.stderr)
        return EXIT_USAGE
    if not args.parca:
        print("en az bir --parca gerekli", file=sys.stderr)
        return EXIT_USAGE

    r = porting.port(args.parca, args.konak)
    if args.json:
        print(json.dumps({
            "komut": "tasi", "hedef_konak": r.target_host,
            "hedef_organizma": r.target_organism,
            "tasinabilir": r.portable,
            "degisiklikler": [{"parca": s.part_id, "kategori": s.category,
                               "eylem": s.action, "yerine": s.replacement,
                               "gerekce": s.reason} for s in r.swaps],
            "tasinan_kaset": r.ported_parts,
            "tikayanlar": r.blockers,
            "iskelet_adaylari": r.backbones,
            "kodon": r.codon_report,
            "notlar": r.notes,
        }, ensure_ascii=False, indent=1))
    else:
        print(porting.to_markdown(r), end="")
    # Taşıma tıkandıysa boru hattının bunu fark etmesi gerekir.
    return EXIT_OK if r.portable else EXIT_FINDINGS


# ---------------------------------------------------------------------------
# iz (denetim izi)
# ---------------------------------------------------------------------------


def cmd_iz(args: argparse.Namespace) -> int:
    from core import audit, db

    db.init_db()
    if not args.parca:
        print("en az bir --parca gerekli", file=sys.stderr)
        return EXIT_USAGE
    r = audit.record(args.ad or "tasarim", args.parca,
                     host=args.konak or "", purpose=args.amac)
    print(audit.to_json(r) if args.json else audit.to_markdown(r), end="")
    return EXIT_FINDINGS if r.errors else EXIT_OK


# ---------------------------------------------------------------------------
# kimlik
# ---------------------------------------------------------------------------


def cmd_kimlik(args: argparse.Namespace) -> int:
    from core import db, provenance

    db.init_db()
    durum = provenance.verify()
    bozuk = bool(durum["degisen"] or durum["eksik"])
    bilgi = {
        "surum": provenance.VERSION,
        "veri_parmak_izi": provenance.data_fingerprint(),
        "veri_dosyasi": len(provenance.compute_checksums()),
        "degisen": durum["degisen"], "eksik": durum["eksik"],
        "yeni": durum["yeni"],
        "parca": len(db.list_parts()), "konak": len(db.list_hosts()),
    }
    if args.json:
        print(json.dumps(bilgi, ensure_ascii=False, indent=1))
    else:
        print(provenance.describe())
        print(f"{bilgi['parca']} parça, {bilgi['konak']} konak")
        if bozuk:
            print("DİKKAT: izlenen veri dosyaları değişmiş; bu kurulumla "
                  "üretilen çıktı, depodaki veriyle aynı olmayabilir.")
    return EXIT_FINDINGS if bozuk else EXIT_OK


# ---------------------------------------------------------------------------
# Çözümleyici
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vectorlab-cli",
        description="VectorLab komut satırı: toplu denetim, primer, "
                    "sipariş, konak taşıma ve denetim izi.")
    p.add_argument("--json", action="store_true",
                   help="makine okunur JSON çıktı")
    alt = p.add_subparsers(dest="komut", required=True)

    d = alt.add_parser("denetle", help="GenBank/FASTA dosyalarını denetle")
    d.add_argument("dosya", nargs="+")
    d.add_argument("--konak", default="", help="konak kimliği (ör. BL21DE3)")
    d.add_argument("--amac", default="klonlama",
                   choices=["ifade", "klonlama"],
                   help="dışarıdan gelen plazmit için varsayılan "
                        "'klonlama': ifade öğeleri aranmaz")
    d.add_argument("--esik", type=float, default=92.0,
                   help="parça tanıma benzerlik eşiği (%%), varsayılan 92")
    d.set_defaults(func=cmd_denetle)

    pr_ = alt.add_parser("primer", help="koloni PCR ve dizileme çapaları")
    pr_.add_argument("dosya")
    pr_.add_argument("--bas", type=int, required=True,
                     help="insertin başlangıç konumu (1 tabanlı)")
    pr_.add_argument("--son", type=int, required=True,
                     help="insertin bitiş konumu")
    pr_.add_argument("--dogrusal", action="store_true",
                     help="dizi dairesel değil doğrusal")
    pr_.set_defaults(func=cmd_primer)

    s = alt.add_parser("siparis", help="malzeme listesi ve takvim")
    s.add_argument("dosya")
    s.add_argument("--bas", type=int, default=1)
    s.add_argument("--son", type=int, default=300)
    s.add_argument("--parca", action="append",
                   help="tasarımın parça kimliği (birden çok kez verilir)")
    s.add_argument("--enzim", action="append", help="kesim enzimi")
    s.add_argument("--konak", default="")
    s.add_argument("--csv", action="store_true",
                   help="malzeme listesini CSV olarak ver")
    s.add_argument("--oligo-csv", action="store_true",
                   help="yalnızca oligo sipariş sayfasını CSV olarak ver")
    s.set_defaults(func=cmd_siparis)

    t = alt.add_parser("tasi", help="kaseti başka bir konağa uyarla")
    t.add_argument("konak", help="hedef konak kimliği (ör. LLNZ9000)")
    t.add_argument("--parca", action="append", required=False,
                   help="kasetin parça kimliği (birden çok kez verilir)")
    t.set_defaults(func=cmd_tasi)

    i = alt.add_parser("iz", help="denetim izi / uyum dosyası")
    i.add_argument("--parca", action="append", required=False)
    i.add_argument("--ad", default="", help="tasarımın adı")
    i.add_argument("--konak", default="")
    i.add_argument("--amac", default="ifade",
                   choices=["ifade", "klonlama"])
    i.set_defaults(func=cmd_iz)

    k = alt.add_parser("kimlik", help="sürüm ve veri parmak izi")
    k.set_defaults(func=cmd_kimlik)
    return p


def main(argv: list[str] | None = None) -> int:
    p = build_parser()
    args = p.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:                           # pragma: no cover
        return EXIT_USAGE


if __name__ == "__main__":                              # pragma: no cover
    sys.exit(main())
