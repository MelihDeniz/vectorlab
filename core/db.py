"""SQLite veri katmanı: parçalar, yapılar, literatür kayıtları."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any

from core import provenance
from core.extracted_seed import extracted_parts
from core.disease_seed import disease_parts
from core.igem_seed import igem_parts
from core.signal_seed import signal_parts
from core.parts_seed import ALL_PARTS, HOSTS
from core.pku_seed import PKU_HOSTS, pku_parts

DB_PATH = Path(__file__).resolve().parent.parent / "data" / "vectorlab.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS parts (
    part_id          TEXT PRIMARY KEY,
    name             TEXT NOT NULL,
    category         TEXT NOT NULL,
    sequence         TEXT NOT NULL DEFAULT '',
    seq_status       TEXT NOT NULL DEFAULT 'full',
    length           INTEGER,
    strength         TEXT,
    induction        TEXT,
    requires         TEXT NOT NULL DEFAULT '[]',
    incompatible_with TEXT NOT NULL DEFAULT '[]',
    copy_number      TEXT,
    inc_group        TEXT,
    notes            TEXT NOT NULL DEFAULT '',
    accession        TEXT NOT NULL DEFAULT '',
    pmid             TEXT NOT NULL DEFAULT '',
    source           TEXT NOT NULL DEFAULT 'seed',
    created_at       TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS hosts (
    host_id       TEXT PRIMARY KEY,
    name          TEXT NOT NULL,
    organism      TEXT NOT NULL DEFAULT 'E. coli',
    use           TEXT,
    t7_polymerase INTEGER NOT NULL DEFAULT 0,
    laciq         INTEGER NOT NULL DEFAULT 0,
    lacz_dm15     INTEGER NOT NULL DEFAULT 0,
    recA          INTEGER NOT NULL DEFAULT 1,
    endA          INTEGER NOT NULL DEFAULT 1,
    rare_trna     INTEGER NOT NULL DEFAULT 0,
    pir           INTEGER NOT NULL DEFAULT 0,
    genotype      TEXT NOT NULL DEFAULT '',
    notes         TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS constructs (
    construct_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name         TEXT NOT NULL UNIQUE,
    host_id      TEXT,
    topology     TEXT NOT NULL DEFAULT 'circular',
    part_ids     TEXT NOT NULL DEFAULT '[]',
    insert_seq   TEXT NOT NULL DEFAULT '',
    insert_name  TEXT NOT NULL DEFAULT '',
    notes        TEXT NOT NULL DEFAULT '',
    provenance   TEXT NOT NULL DEFAULT '',
    created_at   TEXT NOT NULL DEFAULT (datetime('now')),
    updated_at   TEXT NOT NULL DEFAULT (datetime('now'))
);

-- Bir yapının her kaydedilişi burada da saklanır. `constructs` tablosu
-- yapının BUGÜNKÜ halini tutar; bu tablo nasıl o hale geldiğini.
--
-- Neden gerekli: bir tasarımı altı ay sonra açtığınızda aynı sonucu
-- alıp almadığınız, kodun sürümüne ve o sırada yüklü veri dosyalarına
-- bağlıdır. `provenance` alanı her sürüm için o damgayı saklar, yani
-- "bu diziyi hangi veriyle ürettim" sorusu sonradan cevaplanabilir.
-- Laboratuvar defteri mantığı budur: kayıt, sonucun kendisi kadar
-- önemlidir.
CREATE TABLE IF NOT EXISTS construct_versions (
    version_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    name         TEXT NOT NULL,
    host_id      TEXT,
    part_ids     TEXT NOT NULL DEFAULT '[]',
    insert_seq   TEXT NOT NULL DEFAULT '',
    insert_name  TEXT NOT NULL DEFAULT '',
    notes        TEXT NOT NULL DEFAULT '',
    provenance   TEXT NOT NULL DEFAULT '',
    saved_at     TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_versions_name
    ON construct_versions(name, version_id DESC);

CREATE TABLE IF NOT EXISTS literature (
    pmid      TEXT PRIMARY KEY,
    doi       TEXT NOT NULL DEFAULT '',
    title     TEXT NOT NULL DEFAULT '',
    authors   TEXT NOT NULL DEFAULT '',
    journal   TEXT NOT NULL DEFAULT '',
    year      TEXT NOT NULL DEFAULT '',
    abstract  TEXT NOT NULL DEFAULT '',
    is_oa     INTEGER NOT NULL DEFAULT 0,
    source_db TEXT NOT NULL DEFAULT 'europepmc',
    saved_at  TEXT NOT NULL DEFAULT (datetime('now'))
);

CREATE TABLE IF NOT EXISTS part_refs (
    part_id TEXT NOT NULL,
    pmid    TEXT NOT NULL,
    note    TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (part_id, pmid)
);

CREATE INDEX IF NOT EXISTS idx_parts_category ON parts(category);
"""

_LIST_FIELDS = ("requires", "incompatible_with")

# Tabloya yazılan alanlar. Parça kayıtları bunun dışında alan taşıyabilir
# (örneğin iGEM kayıtlarındaki `igem_id`); fazlalıklar yazmadan önce atılır.
_PART_COLUMNS = (
    "part_id", "name", "category", "sequence", "seq_status", "length",
    "strength", "induction", "requires", "incompatible_with", "copy_number",
    "inc_group", "notes", "accession", "pmid", "source",
)


# ---------------------------------------------------------------------------
# Bağlantı
# ---------------------------------------------------------------------------


def connect() -> sqlite3.Connection:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB_PATH, check_same_thread=False)
    con.row_factory = sqlite3.Row
    con.execute("PRAGMA foreign_keys = ON")
    return con


def init_db(force_reseed: bool = False) -> None:
    """Şemayı kur ve başlangıç verisini yükle."""
    con = connect()
    with con:
        con.executescript(SCHEMA)

        if force_reseed:
            con.execute(
                "DELETE FROM parts WHERE source IN "
                "('seed','igem','vector','signal','pku','disease')")
            con.execute("DELETE FROM hosts")

        n_parts = con.execute(
            "SELECT COUNT(*) FROM parts WHERE source = 'seed'"
        ).fetchone()[0]
        if n_parts == 0:
            for part in ALL_PARTS:
                _insert_part(con, part)

        # iGEM parçaları ayrı bir kaynaktan gelir. Sayı tutmuyorsa veri
        # dosyası yenilenmiş demektir: o kaynağı baştan yaz.
        igem = igem_parts()
        n_igem = con.execute(
            "SELECT COUNT(*) FROM parts WHERE source = 'igem'"
        ).fetchone()[0]
        if n_igem != len(igem):
            con.execute("DELETE FROM parts WHERE source = 'igem'")
            for part in igem:
                _insert_part(con, part)

        # Referans vektörlerden çıkarılmış parçalar üçüncü kaynak.
        cikarilan = extracted_parts()
        n_vektor = con.execute(
            "SELECT COUNT(*) FROM parts WHERE source = 'vector'"
        ).fetchone()[0]
        if n_vektor != len(cikarilan):
            con.execute("DELETE FROM parts WHERE source = 'vector'")
            for part in cikarilan:
                _insert_part(con, part)

        sinyaller = signal_parts()
        n_sinyal = con.execute(
            "SELECT COUNT(*) FROM parts WHERE source = 'signal'"
        ).fetchone()[0]
        if n_sinyal != len(sinyaller):
            con.execute("DELETE FROM parts WHERE source = 'signal'")
            for part in sinyaller:
                _insert_part(con, part)

        # PKU (terapötik canlı bakteri) parçaları. Dizisi sonradan
        # çekilen parçalar stub'dan full'a döneceği için sayıya değil
        # (durum, uzunluk, erişim) üçlüsüne bakılır.
        pku = pku_parts()
        var = {r[0]: (r[1], r[2], r[3]) for r in con.execute(
            "SELECT part_id, seq_status, length, accession FROM parts "
            "WHERE source = 'pku'")}
        istenen = {p["part_id"]: (p["seq_status"], p["length"],
                                  p["accession"]) for p in pku}
        if var != istenen:
            con.execute("DELETE FROM parts WHERE source = 'pku'")
            for part in pku:
                _insert_part(con, part)

        # Hastalık vektörlerinin terapötik genleri. Burada da sayıya
        # değil üçlüye bakılır: veri dosyası yenilendiğinde gen sayısı
        # aynı kalıp dizi değişebilir (yeni bir UniProt sürümü, başka
        # bir suşun alleli) ve sayıya bakan kontrol bunu kaçırır.
        terapotik = disease_parts()
        var_h = {r[0]: (r[1], r[2], r[3]) for r in con.execute(
            "SELECT part_id, seq_status, length, accession FROM parts "
            "WHERE source = 'disease'")}
        istenen_h = {p["part_id"]: (p["seq_status"], p["length"],
                                    p["accession"]) for p in terapotik}
        if var_h != istenen_h:
            con.execute("DELETE FROM parts WHERE source = 'disease'")
            for part in terapotik:
                _insert_part(con, part)

        n_hosts = con.execute("SELECT COUNT(*) FROM hosts").fetchone()[0]
        if n_hosts == 0:
            for host in HOSTS:
                con.execute(
                    """INSERT INTO hosts (host_id, name, organism, use,
                           t7_polymerase, laciq, lacz_dm15, recA, endA,
                           rare_trna, pir, genotype, notes)
                       VALUES (:host_id, :name, :organism, :use,
                           :t7_polymerase, :laciq, :lacz_dm15, :recA, :endA,
                           :rare_trna, :pir, :genotype, :notes)""",
                    {"organism": "E. coli", "pir": False, **host,
                     **{k: int(host.get(k, False)) for k in
                        ("t7_polymerase", "laciq", "lacz_dm15",
                         "recA", "endA", "rare_trna", "pir")}},
                )

        # Konak tablosu yalnızca boşken dolduruluyordu; mevcut bir
        # veritabanına sonradan eklenen konaklar bu yüzden eksik kalırdı.
        for host in PKU_HOSTS:
            con.execute(
                """INSERT OR IGNORE INTO hosts (host_id, name, organism, use,
                       t7_polymerase, laciq, lacz_dm15, recA, endA,
                       rare_trna, pir, genotype, notes)
                   VALUES (:host_id, :name, :organism, :use,
                       :t7_polymerase, :laciq, :lacz_dm15, :recA, :endA,
                       :rare_trna, :pir, :genotype, :notes)""",
                {"pir": 0, **host,
                 **{k: int(host.get(k, False)) for k in
                    ("t7_polymerase", "laciq", "lacz_dm15",
                     "recA", "endA", "rare_trna")}},
            )
    con.close()


def _insert_part(con: sqlite3.Connection, part: dict[str, Any]) -> None:
    row = {k: part.get(k) for k in _PART_COLUMNS}
    for f in _LIST_FIELDS:
        row[f] = json.dumps(part.get(f) or [], ensure_ascii=False)
    con.execute(
        """INSERT OR REPLACE INTO parts
           (part_id, name, category, sequence, seq_status, length, strength,
            induction, requires, incompatible_with, copy_number, inc_group,
            notes, accession, pmid, source)
           VALUES (:part_id, :name, :category, :sequence, :seq_status, :length,
            :strength, :induction, :requires, :incompatible_with, :copy_number,
            :inc_group, :notes, :accession, :pmid, :source)""",
        row,
    )


def _row_to_part(row: sqlite3.Row) -> dict[str, Any]:
    d = dict(row)
    for f in _LIST_FIELDS:
        try:
            d[f] = json.loads(d.get(f) or "[]")
        except (TypeError, json.JSONDecodeError):
            d[f] = []
    return d


# ---------------------------------------------------------------------------
# Parçalar
# ---------------------------------------------------------------------------


def list_parts(category: str | None = None, search: str = "") -> list[dict[str, Any]]:
    con = connect()
    sql = "SELECT * FROM parts"
    clauses, params = [], []
    if category:
        clauses.append("category = ?")
        params.append(category)
    if search:
        clauses.append("(name LIKE ? OR part_id LIKE ? OR notes LIKE ?)")
        like = f"%{search}%"
        params += [like, like, like]
    if clauses:
        sql += " WHERE " + " AND ".join(clauses)
    sql += " ORDER BY category, name"
    rows = con.execute(sql, params).fetchall()
    con.close()
    return [_row_to_part(r) for r in rows]


def get_part(part_id: str) -> dict[str, Any] | None:
    con = connect()
    row = con.execute("SELECT * FROM parts WHERE part_id = ?", (part_id,)).fetchone()
    con.close()
    return _row_to_part(row) if row else None


def get_parts(part_ids: list[str]) -> list[dict[str, Any]]:
    """Verilen sırayı KORUYARAK parçaları getir (montaj sırası önemli)."""
    if not part_ids:
        return []
    by_id = {p["part_id"]: p for p in list_parts()}
    return [by_id[pid] for pid in part_ids if pid in by_id]


def upsert_part(part: dict[str, Any]) -> None:
    part = dict(part)
    part.setdefault("source", "user")
    part.setdefault("seq_status", "full" if part.get("sequence") else "stub")
    seq = (part.get("sequence") or "").upper().replace(" ", "").replace("\n", "")
    part["sequence"] = seq
    if part["seq_status"] == "full":
        part["length"] = len(seq)
    for key in ("strength", "induction", "copy_number", "inc_group", "length"):
        part.setdefault(key, None)
    for key in ("notes", "accession", "pmid"):
        part.setdefault(key, "")
    for key in _LIST_FIELDS:
        part.setdefault(key, [])
    con = connect()
    with con:
        _insert_part(con, part)
    con.close()


def delete_part(part_id: str) -> None:
    con = connect()
    with con:
        con.execute("DELETE FROM parts WHERE part_id = ?", (part_id,))
        con.execute("DELETE FROM part_refs WHERE part_id = ?", (part_id,))
    con.close()


def category_counts() -> dict[str, int]:
    con = connect()
    rows = con.execute(
        "SELECT category, COUNT(*) AS n FROM parts GROUP BY category"
    ).fetchall()
    con.close()
    return {r["category"]: r["n"] for r in rows}


# ---------------------------------------------------------------------------
# Konaklar
# ---------------------------------------------------------------------------


def list_hosts() -> list[dict[str, Any]]:
    con = connect()
    rows = con.execute(
        "SELECT * FROM hosts ORDER BY organism != 'E. coli', use, name"
    ).fetchall()
    con.close()
    return [dict(r) for r in rows]


def get_host(host_id: str) -> dict[str, Any] | None:
    con = connect()
    row = con.execute("SELECT * FROM hosts WHERE host_id = ?", (host_id,)).fetchone()
    con.close()
    return dict(row) if row else None


# ---------------------------------------------------------------------------
# Yapılar (constructs)
# ---------------------------------------------------------------------------


def save_construct(
    name: str,
    host_id: str | None,
    part_ids: list[str],
    insert_seq: str = "",
    insert_name: str = "",
    notes: str = "",
) -> bool:
    """Yapıyı kaydet ve değiştiyse yeni bir sürüm yaz.

    Dönen değer yeni bir sürüm oluşturulup oluşturulmadığını söyler.
    Aynı tasarımı üst üste kaydetmek sürüm geçmişini şişirmesin diye,
    bir önceki sürümle birebir aynıysa yeni sürüm yazılmaz - geçmiş
    ancak gerçekten değişen şeyleri gösterirse işe yarar.
    """
    damga = provenance.stamp()
    parcalar = json.dumps(part_ids)
    con = connect()
    with con:
        con.execute(
            """INSERT INTO constructs (name, host_id, part_ids, insert_seq,
                                       insert_name, notes, provenance)
               VALUES (?, ?, ?, ?, ?, ?, ?)
               ON CONFLICT(name) DO UPDATE SET
                 host_id = excluded.host_id,
                 part_ids = excluded.part_ids,
                 insert_seq = excluded.insert_seq,
                 insert_name = excluded.insert_name,
                 notes = excluded.notes,
                 provenance = excluded.provenance,
                 updated_at = datetime('now')""",
            (name, host_id, parcalar, insert_seq, insert_name,
             notes, damga),
        )

        son = con.execute(
            """SELECT host_id, part_ids, insert_seq, insert_name, notes,
                      provenance
               FROM construct_versions WHERE name = ?
               ORDER BY version_id DESC LIMIT 1""", (name,)).fetchone()
        simdi = (host_id, parcalar, insert_seq, insert_name, notes, damga)
        yeni_surum = son is None or tuple(son) != simdi
        if yeni_surum:
            con.execute(
                """INSERT INTO construct_versions
                       (name, host_id, part_ids, insert_seq, insert_name,
                        notes, provenance)
                   VALUES (?, ?, ?, ?, ?, ?, ?)""",
                (name, host_id, parcalar, insert_seq, insert_name,
                 notes, damga),
            )
    con.close()
    return yeni_surum


def get_construct(name: str) -> dict[str, Any] | None:
    """Tek bir yapı, yoksa None."""
    con = connect()
    row = con.execute("SELECT * FROM constructs WHERE name = ?",
                      (name,)).fetchone()
    con.close()
    if row is None:
        return None
    d = dict(row)
    d["part_ids"] = _parcalari_coz(d.get("part_ids"))
    return d


def list_construct_versions(name: str) -> list[dict[str, Any]]:
    """Bir yapının sürüm geçmişi, en yenisi başta."""
    con = connect()
    rows = con.execute(
        """SELECT * FROM construct_versions WHERE name = ?
           ORDER BY version_id DESC""", (name,)).fetchall()
    con.close()
    out = []
    for r in rows:
        d = dict(r)
        d["part_ids"] = _parcalari_coz(d.get("part_ids"))
        out.append(d)
    return out


def rename_construct(eski: str, yeni: str) -> bool:
    """Yapıyı yeniden adlandır; sürüm geçmişi de taşınır.

    Yeni ad kullanımdaysa False döner ve hiçbir şey değişmez.
    """
    eski, yeni = eski.strip(), yeni.strip()
    if not yeni or eski == yeni:
        return False
    con = connect()
    try:
        if con.execute("SELECT 1 FROM constructs WHERE name = ?",
                       (yeni,)).fetchone():
            return False
        with con:
            con.execute(
                "UPDATE constructs SET name = ?, updated_at = datetime('now')"
                " WHERE name = ?", (yeni, eski))
            con.execute(
                "UPDATE construct_versions SET name = ? WHERE name = ?",
                (yeni, eski))
        return True
    finally:
        con.close()


def _parcalari_coz(ham: Any) -> list[str]:
    """`part_ids` alanını listeye çevir; bozuksa boş liste."""
    try:
        val = json.loads(ham)
    except (TypeError, json.JSONDecodeError):
        return []
    return val if isinstance(val, list) else []


def list_constructs() -> list[dict[str, Any]]:
    """Kayıtlı yapılar, en son güncellenen başta; sürüm sayısıyla."""
    con = connect()
    rows = con.execute(
        """SELECT c.*, (SELECT COUNT(*) FROM construct_versions v
                        WHERE v.name = c.name) AS version_count
           FROM constructs c ORDER BY c.updated_at DESC""").fetchall()
    con.close()
    out = []
    for r in rows:
        d = dict(r)
        d["part_ids"] = _parcalari_coz(d.get("part_ids"))
        out.append(d)
    return out


def delete_construct(name: str) -> None:
    """Yapıyı ve sürüm geçmişini sil."""
    con = connect()
    with con:
        con.execute("DELETE FROM constructs WHERE name = ?", (name,))
        con.execute("DELETE FROM construct_versions WHERE name = ?", (name,))
    con.close()


# ---------------------------------------------------------------------------
# Literatür
# ---------------------------------------------------------------------------


def save_reference(ref: dict[str, Any]) -> None:
    con = connect()
    with con:
        con.execute(
            """INSERT OR REPLACE INTO literature
               (pmid, doi, title, authors, journal, year, abstract, is_oa, source_db)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                str(ref.get("pmid") or ref.get("id") or ""),
                ref.get("doi", ""),
                ref.get("title", ""),
                ref.get("authors", ""),
                ref.get("journal", ""),
                str(ref.get("year", "")),
                ref.get("abstract", ""),
                int(bool(ref.get("is_oa"))),
                ref.get("source_db", "europepmc"),
            ),
        )
    con.close()


def list_references() -> list[dict[str, Any]]:
    con = connect()
    rows = con.execute("SELECT * FROM literature ORDER BY saved_at DESC").fetchall()
    con.close()
    return [dict(r) for r in rows]


def delete_reference(pmid: str) -> None:
    con = connect()
    with con:
        con.execute("DELETE FROM literature WHERE pmid = ?", (pmid,))
        con.execute("DELETE FROM part_refs WHERE pmid = ?", (pmid,))
    con.close()


def link_part_ref(part_id: str, pmid: str, note: str = "") -> None:
    con = connect()
    with con:
        con.execute(
            "INSERT OR REPLACE INTO part_refs (part_id, pmid, note) VALUES (?, ?, ?)",
            (part_id, pmid, note),
        )
    con.close()


def refs_for_part(part_id: str) -> list[dict[str, Any]]:
    con = connect()
    rows = con.execute(
        """SELECT l.*, pr.note AS link_note FROM part_refs pr
           JOIN literature l ON l.pmid = pr.pmid
           WHERE pr.part_id = ?""",
        (part_id,),
    ).fetchall()
    con.close()
    return [dict(r) for r in rows]
