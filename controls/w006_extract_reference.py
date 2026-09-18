#!/usr/bin/env python3
"""Build controls/w006_reference.json from the hash-bound Goorley archive.

The archive (W006_PROTOCOL.md §1.1; SHA-256
511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62) is
copyrighted supplementary material and is deliberately NOT stored in this
repository. This script reads a locally held copy, verifies its hash, and
transcribes only the published depth-kerma-rate comparison values needed
for gates W6-G1a/b/c — the data the deposit itself describes as "proposed
for verification of NCT treatment planning software programs" — together
with the archive's MCTAL tally values and relative errors for the same
points. Every value is cited to its archive member.

Sources per beam:
  * published profile  = DepthDoseProfiles/DepthDose.xls (neutron beams)
                         and DepthDoseProfiles/DepthDosePhotons.xls (photon
                         beams), the sheet named for the beam; converted to
                         CSV with LibreOffice headless.
  * mctal              = MCNPoutput/<beam>.m, tallies 14 (boron per ppm),
                         24 (thermal-neutron brain kerma), 34 (fast-neutron
                         brain kerma), 44 (photon brain kerma), 46 axial bins
                         from depth -0.6 to 17.4 cm in 0.4 cm steps.
Findings recorded by this script (see ledger.md):
  * csse (epithermal) and cssp10: xls values equal the MCTAL values.
  * css1/css2/css10/css100/css1M and cssp02..cssp5: xls values differ from
    the deposited MCTAL by up to ~4 % (0.2 cm depth) and typically < 0.5 %
    elsewhere — different runs; the xls error columns are inconsistent. The
    published xls value is the G1 comparison target ("published profile");
    the MCTAL value and relative error are recorded alongside.
  * css253: the deposited css253.m has a 192-cell layout that does not map
    to the deposited deck; only the xls profile is available for that beam.

Usage: w006_extract_reference.py [path/to/supplementary_material_1_1428758-sup-0001.zip]
"""
import csv
import glob
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "w006_reference.json")
ARCHIVE_SHA256 = "511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62"
DEFAULT_ARCHIVE = os.path.expanduser("~/Downloads/supplementary_material_1_1428758-sup-0001.zip")

NEUTRON_SHEETS = {"epithermal": "csse", "0.025eV": "css253", "1kev": "css1", "2kev": "css2",
                  "10kev": "css10", "100kev": "css100", "1000kev": "css1M"}
PHOTON_SHEETS = {".2MeV": "cssp02", ".5MeV": "cssp05", "1MeV": "cssp1", "2MeV": "cssp2",
                 "5MeV": "cssp5", "10MeV": "cssp10"}
TALLY = {"boron": "14", "thermal": "24", "fast": "34", "gamma": "44"}
DEPTH0, DZ, NBINS = -0.6, 0.4, 46


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_mctal(text):
    lines = text.splitlines()
    hdr = lines[0].split()
    nps = int(hdr[5])
    tallies, i = {}, 0
    while i < len(lines):
        if lines[i].startswith("tally"):
            tnum = int(lines[i].split()[1])
            j = i
            while not lines[j].startswith("vals"):
                j += 1
            vals = []
            j += 1
            while not lines[j].startswith("tfc"):
                vals.extend(lines[j].split())
                j += 1
            nums = [float(v) for v in vals]
            tallies[str(tnum)] = [(nums[k], nums[k + 1]) for k in range(0, len(nums), 2)]
            i = j
        i += 1
    return nps, tallies


def fl(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def sheet_points(rows, photon):
    """Return list of (depth, comp, value, abs_err) from a converted sheet,
    reading only the analytic ('Curved Surface') block (columns 1..9)."""
    pts = []
    if photon:
        for r in rows[2:]:
            d, v = fl(r[1]) if len(r) > 1 else None, fl(r[2]) if len(r) > 2 else None
            e = fl(r[3]) if len(r) > 3 else None
            if d is not None and v is not None:
                pts.append((d, "gamma", v, e))
        return pts
    comp_cols = {}
    for j, h in enumerate(rows[1][:11]):
        hl = h.lower()
        for c in ("thermal", "fast", "gamma", "boron"):
            if c in hl and c not in comp_cols:
                comp_cols[c] = j
    for r in rows[3:]:
        d = fl(r[1]) if len(r) > 1 else None
        if d is None:
            continue
        for c, j in comp_cols.items():
            v = fl(r[j]) if j < len(r) else None
            e = fl(r[j + 1]) if j + 1 < len(r) else None
            if v is not None:
                pts.append((d, c, v, e))
    return pts


def main():
    archive = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_ARCHIVE
    h = sha256(archive)
    if h != ARCHIVE_SHA256:
        sys.exit(f"archive hash mismatch: {h}")
    tmp = tempfile.mkdtemp(prefix="w006ref_")
    try:
        with zipfile.ZipFile(archive) as z:
            z.extractall(tmp)
        csvdir = os.path.join(tmp, "csv")
        os.makedirs(csvdir)
        for xls in ("DepthDose.xls", "DepthDosePhotons.xls"):
            subprocess.run(["soffice", "--headless", "--convert-to",
                            "csv:Text - txt - csv (StarCalc):44,34,76,1,,0,false,true,false,false,false,-1",
                            os.path.join(tmp, "DepthDoseProfiles", xls), "--outdir", csvdir],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=300)
        ref = {"archive_sha256": ARCHIVE_SHA256, "archive_file": os.path.basename(archive),
               "doi": "10.1118/1.1428758", "deposit": "E-MPHYA6-29-009201",
               "depth_grid_cm": [round(DEPTH0 + DZ * k, 4) for k in range(NBINS)],
               "beams": {}}
        for sheets, photon, prefix in ((NEUTRON_SHEETS, False, "DepthDose-"),
                                       (PHOTON_SHEETS, True, "DepthDosePhotons-")):
            for sheet, beam in sheets.items():
                path = os.path.join(csvdir, f"{prefix}{sheet}.csv")
                rows = list(csv.reader(open(path, encoding="utf-8", errors="replace")))
                pts = sheet_points(rows, photon)
                mpath = os.path.join(tmp, "MCNPoutput", f"{beam}.m")
                nps, tallies = parse_mctal(open(mpath, encoding="latin-1").read())
                mctal_ok = all(len(tallies.get(TALLY[c], [])) == NBINS
                               for c in set(p[1] for p in pts))
                entries = []
                for d, c, v, e in pts:
                    k = int(round((d - DEPTH0) / DZ))
                    ent = {"depth_cm": d, "component": c, "published_Gy_per_min": v,
                           "published_abs_err": e, "bin": k}
                    if mctal_ok and 0 <= k < NBINS:
                        mv, mr = tallies[TALLY[c]][k]
                        ent["mctal_Gy_per_min"] = mv
                        ent["mctal_rel_err"] = mr
                    entries.append(ent)
                ref["beams"][beam] = {
                    "sheet": f"{'DepthDosePhotons' if photon else 'DepthDose'}.xls:{sheet}",
                    "mctal_file": f"MCNPoutput/{beam}.m", "mctal_nps": nps,
                    "mctal_layout_matches_deck": bool(mctal_ok),
                    "n_points": len(entries), "points": entries}
                print(f"{beam}: {len(entries)} published points; mctal layout ok={mctal_ok}")
        with open(OUT, "w") as f:
            json.dump(ref, f, indent=1, sort_keys=True)
        print("written", OUT, "sha256", sha256(OUT))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
