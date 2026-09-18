#!/usr/bin/env python3
"""W007 harness — re-scoped validation controls, replicated with new seeds.

Thin wrapper over the hash-bound W006 harness (`controls/w006_run.py`,
d3ac76b4…): identical physical model, beams, tally grid, kerma tables and
history counts; only the seeds, run tags and output files differ
(W007_PROTOCOL.md §1, §3). Stages:

    w007_run.py c2        twelve monoenergetic beams, seeds 700-711
    w007_run.py c3        flux-depression run, group-split, seeds 800-802
    w007_run.py assemble  controls/w007_results.json

Depression baseline: the W006 Configuration A profile taken from the
hash-bound controls/w006_results.json (verified here).
"""
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w006_run as H  # noqa: E402  (hash-bound dependency)

ROOT = os.path.normpath(os.path.join(HERE, ".."))
# Test-only overrides (never set for the record): W007_OUT, W007_RUNROOT.
RUNROOT = os.environ.get("W007_RUNROOT", os.path.join(ROOT, "w007_runs"))
OUTDIR = os.environ.get("W007_OUT", HERE)
W006_RESULTS = os.path.join(HERE, "w006_results.json")
W006_RESULTS_SHA256 = "7bd696a4381d3d9e38125fefce007dc52b279e3e62e32e5ac2f1a533122c7c81"
W006_HARNESS_SHA256 = "d3ac76b4c94dc278b4ed23e1c278be297bb13937d4029e1cd2c9c42850c50f1f"
SEEDS_C2 = {tag: 700 + i for i, tag in enumerate(list(H.MONO_NEUTRON_EV) + list(H.MONO_PHOTON_EV))}
SEEDS_C3 = {"thermal": 800, "epithermal": 801, "fast": 802}


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def stage_file(stage):
    return os.path.join(OUTDIR, f"w007_{stage}.json")


def check_deps():
    h = sha(os.path.join(HERE, "w006_run.py"))
    assert h == W006_HARNESS_SHA256, f"w006_run.py hash {h} != frozen"
    h = sha(W006_RESULTS)
    assert h == W006_RESULTS_SHA256, f"w006_results.json hash {h} != frozen"


def dump(stage, payload):
    payload["w007_harness_sha256"] = sha(__file__)
    payload["w006_harness_sha256"] = sha(os.path.join(HERE, "w006_run.py"))
    payload["w007_protocol_sha256"] = sha(os.path.join(ROOT, "W007_PROTOCOL.md"))
    payload["threads"] = H.THREADS
    payload["history_scale"] = H.HISTORY_SCALE
    payload["record_location"] = "controls/" if OUTDIR == HERE else f"TEST:{OUTDIR}"
    payload["timestamp_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(stage_file(stage), "w") as f:
        json.dump(payload, f, indent=1, sort_keys=True)
    print(f"stage {stage} -> {stage_file(stage)}", flush=True)


def stage_c2():
    check_deps()
    os.makedirs(RUNROOT, exist_ok=True)
    runs, profiles = {}, {}
    for tag, E in H.MONO_NEUTRON_EV.items():
        r = H.run_case(f"w007_c2_{tag}", ("mono_n", E), {}, H.N_C2_N[tag], SEEDS_C2[tag],
                       profile=True, roi=False, root=RUNROOT)
        runs[f"c2_{tag}"] = r
        profiles[tag] = H.kerma_profile(r["tallies"])
    for tag, E in H.MONO_PHOTON_EV.items():
        r = H.run_case(f"w007_c2_{tag}", ("mono_p", E), {}, H.N_C2_P, SEEDS_C2[tag],
                       profile=True, roi=False, photon_only=True, root=RUNROOT)
        runs[f"c2_{tag}"] = r
        profiles[tag] = H.kerma_profile(r["tallies"])
    dump("c2", {"runs": runs, "profiles_Gy_per_min": profiles, "depths_cm": H.DEPTHS,
                "seeds": SEEDS_C2})


def stage_c3():
    check_deps()
    os.makedirs(RUNROOT, exist_ok=True)
    ppm = {"brain": H.C3_PPM, "scalp": H.C3_PPM, "cranium": 0.0}
    runs = {}
    for g in ("thermal", "epithermal", "fast"):
        runs[g] = H.run_case(f"w007_c3_{g}", ("group", g), ppm, H.N_C3[g], SEEDS_C3[g],
                             profile=True, roi=False, root=RUNROOT)
    combined = H.combine_groups(runs)
    A = json.load(open(W006_RESULTS))["stages"]["configA"]["profile_Gy_per_min"]
    dump("c3", {"runs": {f"c3_{g}": r for g, r in runs.items()}, "ppm": ppm,
                "profile_Gy_per_min": combined, "depths_cm": H.DEPTHS, "seeds": SEEDS_C3,
                "thermal_dose_depression": H.depression(A, combined),
                "baseline": "W006 Configuration A profile from controls/w006_results.json "
                            + W006_RESULTS_SHA256})


def stage_assemble():
    check_deps()
    out = {"w007_protocol_sha256": sha(os.path.join(ROOT, "W007_PROTOCOL.md")),
           "w006_results_sha256": W006_RESULTS_SHA256, "w006_harness_sha256": W006_HARNESS_SHA256,
           "w007_harness_sha256": sha(__file__), "stages": {}}
    for st in ("c2", "c3"):
        if os.path.exists(stage_file(st)):
            out["stages"][st] = json.load(open(stage_file(st)))
    with open(os.path.join(OUTDIR, "w007_results.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("assembled -> controls/w007_results.json", sorted(out["stages"]))


if __name__ == "__main__":
    {"c2": stage_c2, "c3": stage_c3, "assemble": stage_assemble}[sys.argv[1]]()
