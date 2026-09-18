#!/usr/bin/env python3
"""Number audit: every quantitative statement in manuscript.md is checked
against quantities re-derived by the independent checkers
(controls/check_w006.py --json, controls/check_w007.py --json) or read from
the hash-bound results files; a few record-sourced statements (photon
estimator, transcription defects) are checked against the ledger text.

Usage: audit_numbers.py            (runs both checkers, then audits)
       audit_numbers.py --no-run   (uses cached checker JSON in scratch)
Exit 0 iff every expected string is present and every guard holds.
"""
import json
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
CTRL = os.path.join(ROOT, "controls")
PY = os.environ.get("W003ENV_PYTHON", sys.executable)
SCRATCH = os.environ.get("AUDIT_SCRATCH", os.path.join(HERE, ".audit"))
os.makedirs(SCRATCH, exist_ok=True)
J6P, J7P = os.path.join(SCRATCH, "w006_check.json"), os.path.join(SCRATCH, "w007_check.json")

if "--no-run" not in sys.argv:
    for script, out in (("check_w006.py", J6P), ("check_w007.py", J7P)):
        r = subprocess.run([PY, script, "--json", out], cwd=CTRL, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if r.returncode != 0:
            print(r.stdout); sys.exit(f"{script} exited {r.returncode}")
        print(f"{script}: exit 0 ({r.stdout.strip().splitlines()[0]})")

J6, J7 = json.load(open(J6P)), json.load(open(J7P))
R6 = json.load(open(os.path.join(CTRL, "w006_results.json")))
R7 = json.load(open(os.path.join(CTRL, "w007_results.json")))
LEDGER = open(os.path.join(ROOT, "ledger.md"), encoding="utf-8").read()
MD = open(os.path.join(HERE, "manuscript.md"), encoding="utf-8").read()
TXT = re.sub(r"\^\{([^}]*)\}", r"\1", MD)
TXT = re.sub(r"_\{([^}]*)\}", r"\1", TXT)
TXT = TXT.replace("*", "")
CONV = 1.602176634e-10 * R6["constants"]["SOURCE_PER_MIN"]

checks = []


def expect(label, s, where=TXT):
    ok = s in where
    checks.append((label, s, ok))


def guard(label, cond, detail=""):
    checks.append((label, f"guard: {detail}", bool(cond)))


def rng(vals, fmt="{:.3f}", sep="–"):
    return fmt.format(min(vals)) + sep + fmt.format(max(vals))


# ---------------------------------------------------------------- G1a
a = J6["G1"]["G1a"]; rows = a["rows"]
expect("G1a count", f"{a['points']}/{a['points']}")
expect("G1a worst", f"worst {100*a['worst_abs_dev']:.1f} %")
expect("G1a rel err", f"≤ {100*max(r['mine_rel'] for r in rows):.2f} %")
th14 = [r["ratio_to_published"] for r in rows if r["comp"] == "thermal" and r["depth"] <= 14.0]
guard("thermal ≤14 cm within 2.5 %", max(abs(x - 1) for x in th14) <= 0.025, f"max dev {100*max(abs(x-1) for x in th14):.2f} %")
expect("thermal ratio range ≤14 cm", rng(th14))
th_deep = [r["ratio_to_published"] for r in rows if r["comp"] == "thermal" and r["depth"] > 14.0]
guard("thermal deepest +4–6 %", 1.035 <= min(th_deep) and max(th_deep) <= 1.065, rng(th_deep))
bo = [r["ratio_to_published"] for r in rows if r["comp"] == "boron"]
guard("boron tracks thermal", abs(min(bo) - min(th14 + th_deep)) < 0.002, rng(bo))
fa = [r["ratio_to_published"] for r in rows if r["comp"] == "fast"]
expect("fast ratio range", rng(fa)); guard("fast within 3 %", max(abs(x - 1) for x in fa) <= 0.03)
ga = [r["ratio_to_published"] for r in rows if r["comp"] == "gamma"]
guard("photon +1 to +6 %", 1.01 <= min(ga) and max(ga) <= 1.065, rng(ga)); expect("photon phrase", "+1 to +6 %")

# ---------------------------------------------------------------- G1b
B6, B7 = J6["G1"]["G1b"], J7["G"]["G1b"]
phot = [b for b in B6 if b.startswith("cssp")]
pw = [B6[b]["worst_abs_dev"] for b in phot] + [B7[b]["worst_abs_dev"] for b in phot]
expect("photon beams worst range", f"{100*min(pw):.1f}–{100*max(pw):.1f} %")
expect("photon beams ≤", f"≤ {100*max(pw):.1f} %")
expect("100 keV both seeds", f"{100*min(B6['css100']['worst_abs_dev'], B7['css100']['worst_abs_dev']):.1f}–{100*max(B6['css100']['worst_abs_dev'], B7['css100']['worst_abs_dev']):.1f} %")
expect("1 MeV both seeds", f"{100*min(B6['css1M']['worst_abs_dev'], B7['css1M']['worst_abs_dev']):.1f}–{100*max(B6['css1M']['worst_abs_dev'], B7['css1M']['worst_abs_dev']):.1f} %")
guard("100 keV/1 MeV within 8 %", max(B6[b]["worst_abs_dev"] for b in ("css100", "css1M")) <= 0.08 and max(B7[b]["worst_abs_dev"] for b in ("css100", "css1M")) <= 0.08)
kev = ["css1", "css2", "css10"]
guard("1/2/10 keV first seed all pass", all(B6[b]["failures"] == 0 for b in kev))
fails7 = [(b, r) for b in kev for r in B7[b]["rows"] if not r["pass"]]
guard("second seed five failing points", len(fails7) == 5, str(len(fails7)))
expect("failing ratios", rng([r["ratio"] for _, r in fails7], "{:.2f}"))
ns, stats, sd = [], [], []
for b, r in fails7:
    r6 = next(x for x in B6[b]["rows"] if x["depth"] == r["depth"] and x["comp"] == r["comp"])
    mine, ref = r["mine"], r["ref"]; smine = r["mine_rel"] * mine; sref = r6["ref_rel"] * ref
    ns.append((mine - ref) / math.hypot(smine, sref)); stats += [r["mine_rel"], r6["ref_rel"]]
    sd.append(abs(r6["ratio_to_published"] / r["ratio_to_published"] - 1))
expect("failing nsigma", f"−{abs(max(ns)):.1f}σ to −{abs(min(ns)):.1f}σ")
expect("failing stats", f"{100*min(stats):.0f}–{100*max(stats):.0f} %")
expect("seed spread", f"{100*min(sd):.0f}–{100*max(sd):.0f} %")
guard("all failing points are fast-neutron", all(r["comp"] == "fast" for _, r in fails7))
c253 = B6["css253"]; f253 = [r for r in c253["rows"] if not r["pass"]]
guard("0.0253 eV: 10 outside", c253["failures"] == 10)
expect("0.0253 eV excess", f"+10–{100*c253['worst_abs_dev']:.0f} %")
expect("0.0253 eV depths", f"{min(r['depth'] for r in f253)}–{max(r['depth'] for r in f253)} cm")
th253 = {r["depth"]: r["ratio_to_published"] for r in c253["rows"] if r["comp"] == "thermal"}
expect("0.0253 surface", f"{th253[0.2]:.2f} at the surface"); expect("0.0253 deep", f"{th253[13.4]:.2f} at 13.4 cm")
th253_7 = {r["depth"]: r["ratio_to_published"] for r in B7["css253"]["rows"] if r["comp"] == "thermal"}
expect("0.0253 deep second seed", f"replicated: {th253_7[13.4]:.2f}")

# ---------------------------------------------------------------- G1c / depression
p6 = {x["depth"]: x["depression_pct"] for x in J6["G1"]["G1c"]["peak_points"]}
p7 = {x["depth"]: x["depression_pct"] for x in J7["G"]["G1c"]["five_points"]} if "five_points" in J7["G"]["G1c"] else None
if p7 is None:
    dep7 = {d["depth"]: d["depression_pct"] for d in R7["stages"]["c3"]["thermal_dose_depression"]}
    p7 = {d: dep7[d] for d in p6}
expect("G1c range", f"{min(p6.values()):.1f}–{max(p6.values()):.1f} %")
expect("G1c peak W006", f"{p6[2.6]:.1f} % at the 2.6 cm"); expect("G1c peak W007", f"{p7[2.6]:.1f} % with the second seed")
expect("2.2 cm value", f"{p6[2.2]:.2f} %")
mx = max(abs(p6[d] - p7[d]) for d in p6)
guard("seed replication ≤0.15 pp", mx <= 0.15, f"{mx:.3f}"); expect("replication phrase", "0.15 percentage points")
dep = {d["depth"]: d["depression_pct"] for d in R6["stages"]["c3"]["thermal_dose_depression"]}
expect("depression 0.2 cm", f"{dep[0.2]:.1f} % at 0.2 cm"); expect("depression 4.6 cm", f"{dep[4.6]:.0f} % at 4.6 cm")
deep = [dep[d] for d in dep if 13.0 <= d <= 14.2]
guard("~29 % at 13–14 cm", all(28.5 <= x <= 29.5 for x in deep), rng(deep, "{:.1f}"))
expect("abstract depression", f"{dep[0.2]:.0f} % at the surface to {dep[13.0]:.0f} % at 13 cm")

# ---------------------------------------------------------------- scored block
S = J6["scored"]
for name, tag in (("DS-1", "DS-1"), ("DS-2", "DS-2")):
    o = S[name]
    expect(f"{tag} tumour gap", f"{100*o['gaps']['tumour']:.1f} %")
    expect(f"{tag} brain gap", f"{100*o['gaps']['brain']:.1f} %")
    expect(f"{tag} scalp gap", f"{100*o['gaps']['scalp']:.1f} %")
    guard(f"{tag} soundness", not o["soundness"])
    w = o["G4"]["tumour"]["divergence_window_P_T_Gyw"]
    expect(f"{tag} ratio", f"{o['G4']['tumour']['ratio_certified_over_comparator']:.3f}")
    expect(f"{tag} window", f"[{w[0]:.1f}, {w[1]:.1f})")
d1, d2 = S["DS-1"], S["DS-2"]
expect("abstract ratio", f"{d2['G4']['tumour']['ratio_certified_over_comparator']:.2f}×")
w2 = d2["G4"]["tumour"]["divergence_window_P_T_Gyw"]; w1 = d1["G4"]["tumour"]["divergence_window_P_T_Gyw"]
expect("DS-2 window text", f"{w2[0]:.1f}–{w2[1]:.1f} Gy-w"); expect("DS-1 window text", f"{w1[0]:.1f}–{w1[1]:.1f} Gy-w")
expect("window multiples", f"{w2[0]/20:.1f} to {w2[1]/20:.1f} times"); expect("abstract multiples", f"{w2[0]/20:.1f}–{w2[1]/20:.1f}×")
oar = [x for n in ("DS-1", "DS-2") for roi in ("brain", "scalp") for x in S[n]["G4"][roi]["divergence_window_P_T_Gyw"]]
expect("OAR windows", f"{min(oar):.0f}–{max(oar):.0f} Gy-w")
br = d2["brackets"]
iv = lambda roi, lo_fmt, hi_fmt: f"[{lo_fmt.format((br[roi]['L']-3*br[roi]['sL'])*CONV)}, {hi_fmt.format((br[roi]['U']+3*br[roi]['sU'])*CONV)}]"
expect("DS-2 tumour interval", iv("tumour", "{:.1f}", "{:.1f}"))
expect("DS-2 brain interval", iv("brain", "{:.2f}", "{:.2f}"))
expect("DS-2 scalp interval", iv("scalp", "{:.2f}", "{:.2f}"))
nom = d2["direct"]["interior_6"]
expect("DS-2 nominal rates", f"{nom['tumour'][0]*CONV:.1f}, {nom['brain'][0]*CONV:.2f} and {nom['scalp'][0]*CONV:.2f}")
noncorner = lambda o: [k for k in o["direct"] if not k.startswith("corner_") and "AD" not in k]
guard("22 non-corner evaluations", len(noncorner(d1)) + len(noncorner(d2)) == 22)
ad = sum(1 for o in (d1, d2) for k in o["direct"] if k.endswith("_AD1"))
guard("nine AD re-runs", ad == 9, str(ad)); expect("nine phrase", "nine such re-runs")
# tumour minimum / scalp maximum locations
for name in ("DS-1", "DS-2"):
    o = S[name]; params = R6["stages"]["scored"]["sets"][name]["params"]
    lab = {"corner_lo": params["corner_lo"], "corner_hi": params["corner_hi"]}
    for i, p in enumerate(params["interior"]): lab[f"interior_{i}"] = p
    for i, p in enumerate(params["comparator"]): lab[f"comparator_{i}"] = p
    tv = {k: o["direct"][k]["tumour"] for k in o["direct"]}
    kmin = min(tv, key=lambda k: tv[k][0]); lo = tv["corner_lo"]
    guard(f"{name} tumour min ≈ all-low corner", abs(tv[kmin][0] - lo[0]) <= 3 * math.hypot(tv[kmin][1], lo[1]), f"min at {kmin} {lab.get(kmin.replace('_AD1',''))}")
    sv = {k: o["direct"][k]["scalp"][0] for k in o["direct"]}
    kmax = max(sv, key=sv.get); pm = lab[kmax.replace("_AD1", "")]; ds = R6["constants"]["DS"][name]
    guard(f"{name} scalp max at (low ρT, high ρS, high B)", pm[0] == min(ds["rt"]) and pm[1] == max(ds["rs"]) and pm[2] == max(ds["B"]) and sv[kmax] >= sv["corner_hi"], str(pm))
# corner run cost / histories / precision
walls = [R6["stages"]["scored"]["sets"][n]["runs"][c]["wall_s"] for n in ("DS-1", "DS-2") for c in ("corner_lo", "corner_hi")]
expect("corner cost", f"{min(walls)/60:.0f}–{max(walls)/60:.0f} min")
guard("corner histories 4e7", all(R6["stages"]["scored"]["sets"][n]["runs"][c]["histories"] == 4e7 for n in ("DS-1", "DS-2") for c in ("corner_lo", "corner_hi")))
guard("interior histories 2e7", all(R6["stages"]["scored"]["sets"][n]["runs"][k]["histories"] == 2e7 for n in ("DS-1", "DS-2") for k in R6["stages"]["scored"]["sets"][n]["runs"] if not k.startswith("corner_") and not k.endswith("_AD1")))
guard("corner tallies ≤2.5 %", all(S[n]["worst_relerr"][c] <= 0.025 for n in ("DS-1", "DS-2") for c in ("corner_lo", "corner_hi")), str({n: {c: round(S[n]["worst_relerr"][c], 4) for c in ("corner_lo", "corner_hi")} for n in ("DS-1", "DS-2")}))
# configuration A / C2 histories
A = R6["stages"]["configA"]["runs"]
h = {g: A[f"configA_{g}"]["histories"] for g in ("thermal", "epithermal", "fast")}
guard("epithermal group histories 3e7/6e8/1.5e8", h == {"thermal": 3e7, "epithermal": 6e8, "fast": 1.5e8}, str(h))
c2 = R6["stages"]["c2"]["runs"]
nh = [c2[k]["histories"] for k in c2 if not k.startswith("c2_cssp")]; ph = [c2[k]["histories"] for k in c2 if k.startswith("c2_cssp")]
guard("neutron beams 6–13e7", 6e7 <= min(nh) and max(nh) <= 1.3e8, f"{min(nh):.2g}–{max(nh):.2g}")
guard("photon beams 2e7", set(ph) == {2e7}, str(set(ph)))
# declared constants
C = R6["constants"]
guard("DS-1 set [lo,nom,hi]", C["DS"]["DS-1"] == {"rt": [2.5, 3.5, 4.5], "rs": [0.8, 1.0, 1.2], "B": [15.0, 20.0, 25.0]}, str(C["DS"]["DS-1"]))
guard("DS-2 set [lo,nom,hi]", C["DS"]["DS-2"] == {"rt": [1.58, 3.5, 5.86], "rs": [0.8, 1.0, 1.2], "B": [19.36, 25.0, 31.55]}, str(C["DS"]["DS-2"]))
guard("criteria", C["CRITERIA"] == {"tumour": [">=", 20.0], "brain": ["<=", 11.0], "scalp": ["<=", 10.0]}, str(C["CRITERIA"]))
guard("weights", C["W"] == {"tumour": [3.8, 3.2, 3.2, 1.0], "brain": [1.35, 3.2, 3.2, 1.0], "scalp": [2.5, 3.2, 3.2, 1.0]}, str(C["W"]))
guard("C3 30 ppm", C["C3_PPM"] == 30.0); guard("eps_f 2 %", C["EPS_F"] == 0.02); guard("tumour r=1.5 cm, centre 5.0 cm below z=8.8", C["TUMOUR"] == [1.5, round(8.8 - 5.0, 6)], str(C["TUMOUR"]))
guard("brain ratio 1.0", C["BRAIN_RATIO"] == 1.0); guard("G1c band 7.2–12.1", [round(7.2, 3), round(12.1, 3)] == [round(x, 3) for x in J6["G1"]["G1c"]["band_widened_pct"]], str(J6["G1"]["G1c"]["band_widened_pct"]))
guard("grade", J6["grade"] == {"validation": "INCOMPLETE-VALIDATION", "kernel": "PASS", "G4": "WINDOW-REPORTED-BINARY-UNDECLARED"} and J7["grade"] == "INCOMPLETE-VALIDATION", str((J6["grade"], J7["grade"])))
# ---------------------------------------------------------------- record-sourced (ledger)
expect("photon heating 4.7 eV (ledger 4.66)", "4.66 eV per", LEDGER); expect("5.6e5 eV (ledger)", "5.6e5 eV", LEDGER)
expect("Cl diagnostics 20–40 % (ledger)", "20–40 %", LEDGER); expect("Cl diagnostics 10–15 % (ledger)", "10–15 %", LEDGER)
expect("192-cell (ledger)", "192-cell", LEDGER); expect("2 keV repair rows (ledger)", "14.6 cm", LEDGER)
guard("F-W6-6 substitutions = 5", B6["css2"]["reference_substitutions"] == 5, str(B6["css2"]["reference_substitutions"]))
guard("F-W6-6 touched no other beam", all(B6[b]["reference_substitutions"] == 0 for b in B6 if b != "css2"))
# ---------------------------------------------------------------- report
bad = [c for c in checks if not c[2]]
for label, s, ok in checks:
    print(("PASS " if ok else "FAIL ") + f"{label:45s} {s}")
print(f"\n{len(checks)-len(bad)}/{len(checks)} checks passed")
sys.exit(1 if bad else 0)
