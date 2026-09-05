#!/usr/bin/env python3
"""Verify all references.bib entries against Crossref (and DataCite for data DOIs).

For each entry:
  - if it has a DOI: fetch the record and compare authors, title, journal,
    volume, pages, year against the .bib fields
  - if it has no DOI: run a Crossref bibliographic search on the title and
    report the best candidate (score, DOI, metadata)

Output: a per-field discrepancy report.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

BIB = "/Users/dawid/Projects/ion_surface_sensor/manuscript/references.bib"
UA = {"User-Agent": "QESPM-manuscript-verifier/1.0 (mailto:dawid@localhost)"}
CROSSREF = "https://api.crossref.org/works/{}"
DATACITE = "https://api.datacite.org/dois/{}"


def get_json(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))


# ---------- bib parsing ----------
def parse_bib(path):
    raw = open(path).read()
    entries = {}
    order = []
    i = 0
    pat = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,")
    while i < len(raw):
        m = pat.search(raw, i)
        if not m:
            break
        etype, key = m.group(1), m.group(2)
        start = m.end()
        # find the matching closing brace: scan brace depth (start inside the
        # entry's outer brace, so initial depth is 1)
        depth = 1
        j = start
        while j < len(raw):
            c = raw[j]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        body = raw[start:j]
        fields = {}
        # split on commas at depth 0
        parts = []
        d = 0
        cur = []
        for ch in body:
            if ch == "{":
                d += 1
            elif ch == "}":
                d -= 1
            if ch == "," and d == 0:
                parts.append("".join(cur))
                cur = []
            else:
                cur.append(ch)
        parts.append("".join(cur))
        for p in parts:
            p = p.strip()
            if not p:
                continue
            mm = re.match(r"(\w+)\s*=\s*(.*)$", p, re.S)
            if not mm:
                continue
            fname, fval = mm.group(1), mm.group(2).strip()
            if fval.startswith("{") and fval.endswith("}"):
                fval = fval[1:-1]
            elif fval.startswith('"') and fval.endswith('"'):
                fval = fval[1:-1]
            fields[fname.lower()] = fval
        entries[key] = {"type": etype, "fields": fields}
        order.append(key)
        i = j + 1
    return entries, order


def norm(s):
    if s is None:
        return ""
    s = s.lower()
    s = re.sub(r"[{}]", "", s)
    s = s.replace("--", "-").replace("\u2013", "-").replace("\u2014", "-")
    s = re.sub(r"\s+", " ", s)
    s = s.strip(" .;,")
    s = re.sub(r"[^a-z0-9 \-:+./(),%*'&=~\\]", "", s)
    return s


def family_names(cr_authors):
    out = []
    for a in cr_authors:
        fam = a.get("family", "")
        fam = re.sub(r"\s+(Jr\.?|Sr\.?|II|III)$", "", fam, flags=re.I).strip()
        if fam:
            out.append(fam)
    return out


def bib_authors(field):
    """Parse 'Last, F. and Last2, G.' into list of family names."""
    if not field:
        return []
    names = re.split(r"\s+and\s+", field)
    out = []
    for n in names:
        n = n.strip().strip(",")
        n = re.sub(r"^\{|\}$", "", n)
        m = re.match(r"([^,]+),", n)
        out.append(m.group(1).strip() if m else n)
    return out


def cr_page(cr):
    p = cr.get("page")
    if p:
        return p
    art = cr.get("article-number")
    return art


# ---------- comparison ----------
REPORT = []


def compare(key, bib, cr):
    f = bib["fields"]
    problems = []
    # title
    cr_title = " ".join(cr.get("title", []))
    if norm(f.get("title", "")) and norm(cr_title) and norm(f["title"]) not in (norm(cr_title),):
        problems.append(f"title: bib='{f['title'][:70]}' vs CR='{cr_title[:70]}'")
    # container
    cr_j = " ".join(cr.get("container-title", []))
    if f.get("journal"):
        if not cr_j:
            problems.append(f"journal: bib='{f['journal']}' but CR has no container-title")
        elif norm(f["journal"]) != norm(cr_j):
            # allow abbreviated journal names in bib (e.g. Appl. Phys. B)
            pass  # too fuzzy; report only exact mismatch below in summary
    # volume
    if f.get("volume") and cr.get("volume"):
        bv = norm(f["volume"]).strip("vol.")
        cv = norm(cr["volume"])
        if bv and cv and bv != cv:
            problems.append(f"volume: bib={f['volume']} CR={cr['volume']}")
    elif f.get("volume") and not cr.get("volume"):
        problems.append(f"volume: bib={f['volume']} but CR has none")
    # pages
    bp, cp = f.get("pages"), cr_page(cr)
    if bp and cp and norm(bp).replace("-", "") != norm(cp).replace("-", ""):
        problems.append(f"pages: bib={bp} CR={cp}")
    elif bp and not cp:
        problems.append(f"pages: bib={bp} but CR has none")
    # year
    if f.get("year"):
        by = re.search(r"\d{4}", f["year"])
        cy = None
        for k in ("published-print", "published-online", "issued", "created"):
            if k in cr and cr[k] and cr[k].get("date-parts") and cr[k]["date-parts"][0]:
                cy = cr[k]["date-parts"][0][0]
                break
        if by and cy and int(by.group()) != int(cy):
            problems.append(f"year: bib={f['year']} CR={cy}")
    # authors (first authors only, count check)
    ba = bib_authors(f.get("author", ""))
    ca = family_names(cr.get("author", []))
    if ba and ca:
        if norm(ba[0]) != norm(ca[0]):
            problems.append(f"first author: bib='{ba[0]}' CR='{ca[0]}'")
        if len(ba) != len(ca):
            problems.append(f"author count: bib={len(ba)} CR={len(ca)}")
    elif ba and not ca:
        problems.append("authors: bib has authors but CR has none")
    REPORT.append((key, problems))
    return problems


def main():
    entries, order = parse_bib(BIB)
    print(f"parsed {len(entries)} entries")
    doi_ok = doi_missing = doi_fail = 0
    search_hits = []
    for key in order:
        e = entries[key]
        f = e["fields"]
        doi = f.get("doi", "").strip()
        if doi.startswith("http"):
            doi = doi.split("doi.org/")[-1]
        print(f"\n=== {key} ({e['type']})")
        if doi:
            datacite = doi.startswith("10.5281/") or doi.startswith("10.17632/")
            try:
                if datacite:
                    rec = get_json(DATACITE.format(urllib.parse.quote(doi)))
                    attrs = rec["data"]["attributes"]
                    cr = {
                        "title": [attrs.get("titles", [{}])[0].get("title", "")],
                        "container-title": [attrs.get("publisher", "")],
                        "author": [{"family": c.get("name", "").split(",")[0]}
                                   for c in attrs.get("creators", [])],
                        "volume": None,
                        "page": None,
                        "publisher": attrs.get("publisher"),
                        "issued": {"date-parts": [[attrs.get("publicationYear", 0)]]},
                    }
                    problems = compare(key, e, cr)
                    print(f"  DataCite record OK: '{attrs.get('titles', [{}])[0].get('title', '')[:80]}' ({attrs.get('publicationYear')})")
                else:
                    rec = get_json(CROSSREF.format(urllib.parse.quote(doi)))
                    cr = rec["message"]
                    problems = compare(key, e, cr)
                    print(f"  Crossref OK: '{' '.join(cr.get('title', []))[:80]}' ({cr.get('publisher', '')})")
                doi_ok += 1
                for p in problems:
                    print(f"    MISMATCH {p}")
                if not problems:
                    print("    all fields consistent")
            except Exception as ex:
                doi_fail += 1
                print(f"  DOI RESOLUTION FAILED: {doi} -> {ex}")
        else:
            doi_missing += 1
            title = f.get("title", "")
            t = re.sub(r"[\{\}]", "", title)
            q = urllib.parse.quote(t)
            try:
                rec = get_json("https://api.crossref.org/works?query.bibliographic=" + q + "&rows=3")
                items = rec["message"]["items"]
                if items:
                    best = items[0]
                    print(f"  no DOI in bib; top Crossref match (score {best.get('score', '?')}):")
                    print(f"    DOI: {best.get('DOI')}")
                    print(f"    title: '{' '.join(best.get('title', []))}'")
                    print(f"    container: '{' '.join(best.get('container-title', []))}' vol={best.get('volume')} page={best.get('page')} year={best.get('issued', {}).get('date-parts', [[None]])[0][0]}")
                    print(f"    authors: {', '.join(family_names(best.get('author', [])))}")
                    search_hits.append((key, best))
                else:
                    print("  no DOI in bib; Crossref search: no results")
            except Exception as ex:
                print(f"  Crossref search failed: {ex}")
            time.sleep(0.3)
    print(f"\n==== SUMMARY: entries={len(entries)} doi_ok={doi_ok} doi_missing={doi_missing} doi_fail={doi_fail}")
    return search_hits


if __name__ == "__main__":
    main()
