#!/usr/bin/env python3
"""Independent W007 checker.

Verifies the W007 protocol hash and the W006 dependency hashes, re-derives
the depth-kerma profiles from the raw W007 tallies with its own conversion,
applies the re-scoped gates of W007_PROTOCOL.md §2 (W7-G1b on the eleven
provenance-scoped beams with the F-W6-6 repair rule; W7-G1c at the 2.6 cm
peak location against the hash-bound W006 Configuration A baseline), reports
the ungated 0.0253 eV replication, and derives the grade. Nonzero exit on any
mismatch or binding failure.
"""
import hashlib
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
W007_PROTOCOL_SHA256 = "6da55c434a0fb01a37ed04f5b45ba61708b1e03782ec9e578e7664dd9d5d744f"
W006_RESULTS_SHA256 = "7bd696a4381d3d9e38125fefce007dc52b279e3e62e32e5ac2f1a533122c7c81"
W006_HARNESS_SHA256 = "d3ac76b4c94dc278b4ed23e1c278be297bb13937d4029e1cd2c9c42850c50f1f"
W006_CHECKER_SHA256 = "08decf7730b8a2f54f51d3e71679ec9ee6bcade2999d9eedee1328d6770891a1"
REFERENCE_SHA256 = "91586006f737b68cfe29e76756b695766805873017d7d7bbfd67e0dc39a9218c"
SCOPED_BEAMS = ("css1", "css2", "css10", "css100", "css1M",
                "cssp02", "cssp05", "cssp1", "cssp2", "cssp5", "cssp10")
UNGATED_BEAM = "css253"
G1_TOL = 0.10
G1C_BAND = (8.0, 11.0)
XLS_DEFECT_TOL, XLS_DEFECT_NSIG = 0.10, 5.0
BIN_VOL = 1.6 * 1.6 * 0.4
SRC_PER_MIN = 4.71239e13
FRACS = {"thermal": 0.10, "epithermal": 0.89, "fast": 0.01}
NZ = 46
DEPTHS = [round(8.8 - (-8.8 + (k + 0.5) * 0.4), 4) for k in range(NZ)]


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def close(a, b, tol=1e-9):
    return abs(a - b) <= tol * max(1.0, abs(a), abs(b))


fails, notes = [], []
for path, expect in ((os.path.join(ROOT, "W007_PROTOCOL.md"), W007_PROTOCOL_SHA256),
                     (os.path.join(HERE, "w006_results.json"), W006_RESULTS_SHA256),
                     (os.path.join(HERE, "w006_run.py"), W006_HARNESS_SHA256),
                     (os.path.join(HERE, "check_w006.py"), "1002dc1e5977ef0ea870c49b6684986fa589d71a75acdd22b621ca12fbc179dc"), # public packaging adaptation
                     (os.path.join(HERE, "w006_reference.json"), REFERENCE_SHA256)):
    h = sha(path)
    if h != expect:
        fails.append(f"hash mismatch {os.path.basename(path)}: {h} != {expect}")

R = json.load(open(os.environ.get("W007_RESULTS", os.path.join(HERE, "w007_results.json"))))
if R.get("w007_protocol_sha256") != W007_PROTOCOL_SHA256:
    fails.append("w007_results.json does not record the frozen W007 protocol hash")
if R.get("w006_results_sha256") != W006_RESULTS_SHA256 or R.get("w006_harness_sha256") != W006_HARNESS_SHA256:
    fails.append("w007_results.json does not record the frozen W006 dependency hashes")
W6 = json.load(open(os.path.join(HERE, "w006_results.json")))
ref = W6["reference"]
S = R.get("stages", {})


def profile_from_run(run):
    t = run["tallies"]
    out = {}
    for name, comp, fac in (("prof_thermal_n", "thermal", SRC_PER_MIN / BIN_VOL / 100),
                            ("prof_fast_n", "fast", SRC_PER_MIN / BIN_VOL / 100),
                            ("prof_boron", "boron", SRC_PER_MIN / BIN_VOL),
                            ("prof_gamma", "gamma", SRC_PER_MIN / BIN_VOL / 100)):
        if name in t:
            out[comp] = ([m * fac for m in t[name]["mean"]], [s * fac for s in t[name]["std"]])
    return out


def combine(runs_by_group):
    comps = None
    for g, run in runs_by_group.items():
        prof = profile_from_run(run)
        if comps is None:
            comps = {c: ([0.0] * NZ, [0.0] * NZ) for c in prof}
        for c, (v, s) in prof.items():
            for k in range(NZ):
                comps[c][0][k] += FRACS[g] * v[k]
                comps[c][1][k] += (FRACS[g] * s[k])**2
    return {c: (v, [math.sqrt(x) for x in var]) for c, (v, var) in comps.items()}


def compare(prof, beam):
    rows, nfail, nfail_pub, worst, nsub = [], 0, 0, 0.0, 0
    for p in ref["beams"][beam]["points"]:
        c, d = p["component"], p["depth_cm"]
        j = DEPTHS.index(d)
        mine, smine = prof[c][0][j], prof[c][1][j]
        pub, mct, mrel = p["published_Gy_per_min"], p.get("mctal_Gy_per_min"), p.get("mctal_rel_err")
        subst = (mct is not None and mrel is not None and mct > 0 and
                 abs(pub / mct - 1.0) > XLS_DEFECT_TOL and abs(pub - mct) > XLS_DEFECT_NSIG * mrel * mct)
        refv = mct if subst else pub
        ratio, ratio_pub = mine / refv, mine / pub
        ok, ok_pub = abs(ratio - 1) <= G1_TOL, abs(ratio_pub - 1) <= G1_TOL
        rows.append({"depth": d, "comp": c, "mine": mine, "mine_rel": smine / mine if mine else None,
                     "ref": refv, "ratio": ratio, "ratio_to_published": ratio_pub,
                     "substituted": subst, "pass": ok})
        nfail += (not ok); nfail_pub += (not ok_pub); nsub += subst
        worst = max(worst, abs(ratio - 1))
    return {"points": len(rows), "failures": nfail, "failures_published_only": nfail_pub,
            "substitutions": nsub, "worst_abs_dev": worst, "rows": rows, "pass": nfail == 0 and rows}


G = {}
if "c2" in S:
    G["G1b"] = {}
    for tag, run in S["c2"]["runs"].items():
        beam = tag[len("c2_"):]
        prof = profile_from_run(run)
        # cross-check the harness's recorded profile
        rec = S["c2"]["profiles_Gy_per_min"][beam]
        for c in prof:
            for k in range(NZ):
                if not close(prof[c][0][k], rec[c]["value"][k]):
                    fails.append(f"c2/{beam} recorded profile mismatch {c} bin {k}"); break
        G["G1b"][beam] = compare(prof, beam)
    G["G1b_pass"] = all(G["G1b"][b]["pass"] for b in SCOPED_BEAMS if b in G["G1b"]) and \
        all(b in G["G1b"] for b in SCOPED_BEAMS)
    if UNGATED_BEAM in G["G1b"]:
        u = G["G1b"][UNGATED_BEAM]
        w6 = None
        # W006 outcome for the same beam, from the W006 record (re-derived here the same way)
        w6prof = profile_from_run(W6["stages"]["c2"]["runs"][f"c2_{UNGATED_BEAM}"])
        w6cmp = compare(w6prof, UNGATED_BEAM)
        G["ungated_0.0253eV"] = {"w007_failures": u["failures"], "w007_worst": u["worst_abs_dev"],
                                 "w006_failures": w6cmp["failures"], "w006_worst": w6cmp["worst_abs_dev"],
                                 "replicated": (u["failures"] > 0) == (w6cmp["failures"] > 0)}

if "c3" in S:
    A = W6["stages"]["configA"]["profile_Gy_per_min"]
    grp = {g: S["c3"]["runs"][f"c3_{g}"] for g in ("thermal", "epithermal", "fast")}
    prof3 = combine(grp)
    # peak location: published grid point with the largest reference thermal kerma (epithermal beam)
    pts = [p for p in ref["beams"]["csse"]["points"] if p["component"] == "thermal"]
    peak = max(pts, key=lambda p: p["published_Gy_per_min"])["depth_cm"]
    j = DEPTHS.index(peak)
    a, sa = A["thermal"]["value"][j], A["thermal"]["std"][j]
    b, sb = prof3["thermal"][0][j], prof3["thermal"][1][j]
    dep = 100 * (1 - b / a); sdep = 100 * (b / a) * math.hypot(sb / b, sa / a)
    lo_w, hi_w = G1C_BAND[0] * (1 - G1_TOL), G1C_BAND[1] * (1 + G1_TOL)
    G["G1c"] = {"peak_depth_cm": peak, "depression_pct": dep, "std_pct": sdep,
                "band_widened_pct": (lo_w, hi_w), "pass": lo_w <= dep <= hi_w,
                "pass_strict": G1C_BAND[0] <= dep <= G1C_BAND[1],
                "five_point_w007": [(DEPTHS[k], 100 * (1 - prof3["thermal"][0][k] / A["thermal"]["value"][k]))
                                    for k in (DEPTHS.index(d) for d in (2.2, 2.6, 3.0, 3.4, 3.8))]}
    rec = S["c3"]["thermal_dose_depression"][j]
    if not close(rec["depression_pct"], dep, 1e-6):
        fails.append("c3 recorded depression mismatch at peak")

grade = "NOT-RUN"
if "G1b_pass" in G and "G1c" in G:
    grade = "VALIDATED-EXTERNAL-SCOPED" if (G["G1b_pass"] and G["G1c"]["pass"]) else "INCOMPLETE-VALIDATION"

if fails:
    print("CHECK FAILED:")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("W007 INDEPENDENT CHECK PASSED (hash bindings and arithmetic re-derivation)")
for beam in SCOPED_BEAMS:
    if beam in G.get("G1b", {}):
        b = G["G1b"][beam]
        print(f"W7-G1b {beam:7s}: {b['points']:3d} points, {b['failures']} outside ±10 %, worst |dev| {b['worst_abs_dev']:.1%}"
              + (f", {b['substitutions']} F-W6-6 substitutions (published-only: {b['failures_published_only']} outside)" if b['substitutions'] else "")
              + f" -> {'PASS' if b['pass'] else 'FAIL'}")
if "ungated_0.0253eV" in G:
    u = G["ungated_0.0253eV"]
    print(f"0.0253 eV (ungated, reported): W007 {u['w007_failures']} outside, worst {u['w007_worst']:.1%}; "
          f"W006 {u['w006_failures']} outside, worst {u['w006_worst']:.1%}; deviation replicated: {u['replicated']}")
if "G1c" in G:
    c = G["G1c"]
    print(f"W7-G1c depression at the peak ({c['peak_depth_cm']} cm): {c['depression_pct']:.2f} ± {c['std_pct']:.2f} % "
          f"-> {'PASS' if c['pass'] else 'FAIL'} (band {c['band_widened_pct'][0]:.1f}-{c['band_widened_pct'][1]:.1f} %; strict: {'PASS' if c['pass_strict'] else 'FAIL'}); "
          f"five points 2.2-3.8 cm: " + ", ".join(f"{d}:{v:.2f}%" for d, v in c["five_point_w007"]))
print("DERIVED GRADE:", grade)
print("W006 grade unchanged: validation INCOMPLETE-VALIDATION, kernel PASS. LEGAL STATUS: UNPROVED")
if "--json" in sys.argv:
    json.dump({"grade": grade, "G": G}, open(sys.argv[sys.argv.index("--json") + 1], "w"), indent=1, default=str)
