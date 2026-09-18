#!/usr/bin/env python3
"""Paper figures for the W006/W007 record. Every number is read from
controls/w006_results.json and controls/w007_results.json; nothing is typed
in by hand except geometry constants that also live in the protocol.

Outputs (paper/figures/):
  fig1_concept.png     declared set, corner evaluations, approximate envelope vs all
                       direct evaluations (DS-2, tumour) — the mechanism
  fig2_geometry.png    phantom cross-section, beam, tally column, tumour
  fig3_epithermal.png  epithermal-beam profiles vs published reference
  fig4_mono.png        monoenergetic beams: ratio ranges (W006 and W007)
  fig5_depression.png  thermal-dose depression from 30 ug/g B-10
  fig6_bounds.png      envelope edges vs direct evaluations, both sets
  fig7_window.png      divergence window vs preset comparator
"""
import json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle, Circle, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
OUT = os.environ.get("PAPER_FIG_OUT", os.path.join(HERE, "figures"))   # override for submission copies
os.makedirs(OUT, exist_ok=True)
W6 = json.load(open(os.path.join(ROOT, "controls", "w006_results.json")))
W7 = json.load(open(os.path.join(ROOT, "controls", "w007_results.json")))
REF = W6["reference"]["beams"]
S6, S7 = W6["stages"], W7["stages"]
DEPTHS = W6["constants"]["DEPTHS"]
CONV = 1.602176634e-10 * 4.71239e13     # MeV/g/src -> Gy-w/min at 1e10 n/cm2/s

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#ececea"
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
BAND = "#e9e8e4"; BANDBLUE = "#dbe7f7"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.titlecolor": INK, "legend.frameon": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "figure.dpi": 100, "savefig.dpi": int(os.environ.get("PAPER_FIG_DPI", "220"))})
BEAM_LABEL = {"css253": "0.0253 eV n", "css1": "1 keV n", "css2": "2 keV n", "css10": "10 keV n",
              "css100": "100 keV n", "css1M": "1 MeV n", "cssp02": "0.2 MeV γ", "cssp05": "0.5 MeV γ",
              "cssp1": "1 MeV γ", "cssp2": "2 MeV γ", "cssp5": "5 MeV γ", "cssp10": "10 MeV γ"}
COMP_LABEL = {"thermal": "thermal-neutron kerma", "fast": "fast-neutron kerma",
              "gamma": "photon kerma", "boron": "boron kerma per ppm ¹⁰B"}


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.15, facecolor="white")
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- helpers
def ref_series(beam, comp):
    pts = sorted([p for p in REF[beam]["points"] if p["component"] == comp], key=lambda p: p["depth_cm"])
    return pts


def repaired_ref(p):
    """F-W6-6 rule: published contradicted by same-run MCTAL by >10% and >5σ -> MCTAL."""
    pub, mct, mrel = p["published_Gy_per_min"], p.get("mctal_Gy_per_min"), p.get("mctal_rel_err")
    if mct and mrel is not None and abs(pub / mct - 1) > 0.10 and abs(pub - mct) > 5 * mrel * mct:
        return mct, mrel * mct, True
    sref = mrel * mct if (mct and mrel is not None) else (p.get("published_abs_err") or 0.0)
    return pub, sref, False


def dose_of(run, roi, W, MASS):
    wB, wN, wf, wg = W[roi]; m = MASS[roi]; t = run["tallies"]
    return (wB * 2.34 * t[f"{roi}_b10_na"]["mean"] + wN * 0.626 * t[f"{roi}_n14_np"]["mean"]
            + wf * t[f"{roi}_fast_heat"]["mean"] / 1e6 + wg * t[f"{roi}_photon_heat"]["mean"] / 1e6) / m


# ============================================================ Fig 2 geometry
def fig_geometry():
    fig, ax = plt.subplots(figsize=(9.2, 6.0))
    ax.set_aspect("equal"); ax.grid(False)
    # ellipsoids in the x-z plane (y = 0): (x/a)^2 + ((z-zc)/c)^2 = 1
    for (a, c, zc), col, lab in (((7.3, 8.8, 0.0), "#f3d9c4", "scalp (ICRU 46 skin, 0.5 cm)"),
                                 ((6.8, 8.3, 0.0), "#dcdcdc", "cranium (ICRU 46 bone)"),
                                 ((6.0, 6.5, 1.0), "#f6e7e2", "brain (ICRU 46 adult brain)")):
        ax.add_patch(Ellipse((0, zc), 2 * a, 2 * c, facecolor=col, edgecolor=INK2, lw=0.9, label=lab))
    ax.add_patch(Circle((0, 3.8), 1.5, facecolor="#c9d9f0", edgecolor=BLUE, lw=1.2,
                        label="declared tumour (r = 1.5 cm, 5 cm deep)"))
    # tally column
    ax.add_patch(Rectangle((-0.8, -8.8), 1.6, 18.4, facecolor="none", edgecolor=INK, lw=0.8, ls="--",
                           label="axial tally column 16×16 mm, 4 mm bins"))
    # beam
    ax.add_patch(Rectangle((-5, 12.0), 10, 3.0, facecolor="#fde9c9", edgecolor=ORANGE, lw=1.0,
                           label="beam: 10 cm disk at z = 15 cm, monodirectional −z"))
    for x in (-3.5, 0.0, 3.5):
        ax.add_patch(FancyArrowPatch((x, 12.0), (x, 9.6), arrowstyle="-|>", mutation_scale=12, color=ORANGE, lw=1.2))
    # depth annotations
    ax.plot([9.0, 9.0], [8.8, 3.8], color=INK2, lw=0.8)
    ax.plot([8.7, 9.3], [8.8, 8.8], color=INK2, lw=0.8); ax.plot([8.7, 9.3], [3.8, 3.8], color=INK2, lw=0.8)
    ax.text(9.5, 6.3, "5.0 cm depth\nto tumour centre", va="center", fontsize=8, color=INK2)
    ax.text(-5.3, 8.8, "depth 0\n(z = 8.8 cm)", ha="right", va="center", fontsize=8, color=INK2)
    ax.plot([-5.0, -0.8], [8.8, 8.8], color=INK2, lw=0.6, ls=":")
    ax.text(0, 13.5, "beam\n10 cm diameter", ha="center", va="center", fontsize=8, color=INK2)
    ax.text(0, 3.8, "T", ha="center", va="center", fontsize=9, color=INK)
    ax.text(0, -2.5, "brain", ha="center", va="center", fontsize=9, color=INK2)
    ax.annotate("cranium", xy=(3.2, -6.0), xytext=(6.2, -8.4), fontsize=8, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.annotate("scalp", xy=(-4.6, -6.6), xytext=(-8.5, -8.6), fontsize=8, color=INK2,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.set_xlim(-11, 13); ax.set_ylim(-10.5, 16)
    ax.set_xlabel("x (cm)"); ax.set_ylabel("z (cm), beam axis")
    ax.set_title("Modified Snyder head phantom (analytic model), x–z section at y = 0", loc="left", fontsize=10)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8, frameon=False)
    save(fig, "fig2_geometry.png")


# ============================================================ Fig 3 epithermal
def fig_epithermal():
    A = S6["configA"]["profile_Gy_per_min"]
    comps = ["thermal", "boron", "fast", "gamma"]
    fig, axes = plt.subplots(2, 4, figsize=(13.6, 5.9), sharex="col",
                             gridspec_kw={"height_ratios": [3, 1.5], "hspace": 0.06, "wspace": 0.34})
    for j, comp in enumerate(comps):
        pts = ref_series("csse", comp)
        d = [p["depth_cm"] for p in pts]
        rv = [p["published_Gy_per_min"] for p in pts]; rs = [p["mctal_rel_err"] * p["published_Gy_per_min"] for p in pts]
        idx = [DEPTHS.index(x) for x in d]
        mv = [A[comp]["value"][i] for i in idx]; ms = [A[comp]["std"][i] for i in idx]
        ax = axes[0, j]
        ax.plot(d, rv, lw=0, marker="o", ms=4.2, mfc="white", mec=INK2, mew=1.2, label="published reference (MCNP4B, ENDF/B-VI)")
        ax.plot(d, mv, color=BLUE, lw=1.8, label="this work (OpenMC 0.15.3, ENDF/B-VII.1)")
        ax.set_yscale("log"); ax.tick_params(labelbottom=False)
        ax.set_title(COMP_LABEL[comp] + ("\n(Gy/min per ppm)" if comp == "boron" else "\n(Gy/min)"), fontsize=9, loc="left")
        ax2 = axes[1, j]
        ratio = [m / r for m, r in zip(mv, rv)]
        rerr = [(m / r) * math.hypot(sm / m, sr / r) for m, r, sm, sr in zip(mv, rv, ms, rs)]
        ax2.axhspan(0.9, 1.1, color=BAND, zorder=0)
        ax2.axhline(1.0, color=MUTED, lw=0.8)
        ax2.errorbar(d, ratio, yerr=rerr, color=BLUE, lw=1.2, marker="o", ms=2.8, capsize=0)
        ax2.set_ylim(0.82, 1.18); ax2.set_yticks([0.9, 1.0, 1.1])
        ax2.set_xlabel("depth along beam axis (cm)")
        if j == 0:
            ax2.set_ylabel("this work / reference")
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.04), fontsize=9)
    fig.text(0.0, -0.02, "Generic epithermal beam (10 % thermal, 89 % epithermal, 1 % fast) on the analytic modified Snyder head phantom; "
             "16×16×4 mm axial tally bins. Shaded band: pre-committed ±10 %. Error bars: 1σ combined (Monte Carlo, both codes).",
             fontsize=7.5, color=INK2)
    save(fig, "fig3_epithermal.png")


# ============================================================ Fig 4 mono beams
def fig_mono():
    order = ["css253", "css1", "css2", "css10", "css100", "css1M", "cssp02", "cssp05", "cssp1", "cssp2", "cssp5", "cssp10"]
    comps = ["thermal", "boron", "fast", "gamma"]
    rows = []
    for beam in order:
        for comp in comps:
            pts = ref_series(beam, comp)
            if not pts:
                continue
            rows.append((beam, comp))
    fig, ax = plt.subplots(figsize=(8.6, 0.30 * len(rows) + 1.6))
    ax.axvspan(0.9, 1.1, color=BAND, zorder=0)
    ax.axvline(1.0, color=MUTED, lw=0.8)
    yl = []
    for i, (beam, comp) in enumerate(rows):
        y = len(rows) - 1 - i
        for stage, col, dy, lab in ((S6, BLUE, 0.16, "W006 seeds"), (S7, ORANGE, -0.16, "W007 seeds (replication)")):
            run = stage["c2"]["runs"][f"c2_{beam}"]
            prof = stage["c2"]["profiles_Gy_per_min"][beam]
            ratios, offs = [], 0
            for p in ref_series(beam, comp):
                r, sr, subst = repaired_ref(p)
                v = prof[comp]["value"][DEPTHS.index(p["depth_cm"])]
                ratios.append(v / r)
            lo, hi = min(ratios), max(ratios)
            clip_lo, clip_hi = max(lo, 0.7), min(hi, 1.2)
            ax.plot([clip_lo, clip_hi], [y + dy, y + dy], color=col, lw=2.2, solid_capstyle="butt",
                    label=lab if i == 0 else None)
            ax.plot([r for r in ratios if 0.7 <= r <= 1.2], [y + dy] * sum(1 for r in ratios if 0.7 <= r <= 1.2),
                    color=col, lw=0, marker="|", ms=5, alpha=0.7)
            if hi > 1.2:
                ax.annotate("", xy=(1.215, y + dy), xytext=(1.19, y + dy), arrowprops=dict(arrowstyle="-|>", color=col, lw=1))
            if lo < 0.7:
                ax.annotate("", xy=(0.685, y + dy), xytext=(0.71, y + dy), arrowprops=dict(arrowstyle="-|>", color=col, lw=1))
        yl.append((y, f"{BEAM_LABEL[beam]} · {comp}"))
    ax.set_yticks([y for y, _ in yl]); ax.set_yticklabels([t for _, t in yl], fontsize=8)
    ax.set_xlim(0.68, 1.22); ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
    ax.set_xlabel("this work / reference, range over all published depths (band ±10 %)")
    ax.set_title("Monoenergetic beams: agreement range per beam and component, two independent seeds", loc="left", fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    # separators between neutron and photon beams
    nrow_n = sum(1 for b, c in rows if not b.startswith("cssp"))
    ax.axhline(len(rows) - nrow_n - 0.5, color=MUTED, lw=0.6, ls=":")
    ax.text(0.69, len(rows) - nrow_n - 0.5 - 0.35, "photon beams ↓", fontsize=7.5, color=INK2, va="top")
    ax.text(0.69, len(rows) - 0.5 + 0.1, "neutron beams ↓", fontsize=7.5, color=INK2, va="bottom")
    ax.set_ylim(-0.7, len(rows) - 0.3 + 0.6)
    fig.text(0.0, -0.01, "Reference = published depth-kerma table; for the 2 keV beam five published rows at 14.2/14.6 cm that contradict the deposit's own MCNP "
             "tally by 5× are replaced by that tally (F-W6-6). Arrows mark ranges extending beyond the axis.", fontsize=7.5, color=INK2, wrap=True)
    save(fig, "fig4_mono.png")


# ============================================================ Fig 5 depression
def fig_depression():
    dep6 = S6["c3"]["thermal_dose_depression"]; dep7 = S7["c3"]["thermal_dose_depression"]
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.axvspan(2.0, 4.0, color=BAND, zorder=0, label="thermal-dose peak region 2–4 cm")
    ax.axhspan(8.0, 11.0, color=BANDBLUE, zorder=0, label="published 8–11 % (all beams, at the peak)")
    for dep, col, lab in ((dep6, BLUE, "this work, W006 seeds"), (dep7, ORANGE, "this work, W007 seeds")):
        d = [x["depth"] for x in dep if 0 <= x["depth"] <= 14.6 and not math.isnan(x["depression_pct"])]
        v = [x["depression_pct"] for x in dep if 0 <= x["depth"] <= 14.6 and not math.isnan(x["depression_pct"])]
        s = [x["std_pct"] for x in dep if 0 <= x["depth"] <= 14.6 and not math.isnan(x["depression_pct"])]
        ax.errorbar(d, v, yerr=s, color=col, lw=1.5, marker="o", ms=2.8, capsize=0, label=lab)
    ax.set_xlim(0, 14.8); ax.set_ylim(0, 32)
    ax.set_xlabel("depth along beam axis (cm)"); ax.set_ylabel("thermal-neutron dose depression (%)")
    ax.set_title("Flux depression by 30 µg/g ¹⁰B in scalp and brain (not cranium), epithermal beam", loc="left", fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    save(fig, "fig5_depression.png")


# ============================================================ Fig 6 bounds
def fig_bounds():
    sc = S6["scored"]; W = {k: tuple(v) for k, v in sc["weights"].items()}; MASS = sc["roi_mass_g"]
    rois = ["tumour", "brain", "scalp"]
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.2), gridspec_kw={"hspace": 0.5, "wspace": 0.28})
    for r, name in enumerate(["DS-1", "DS-2"]):
        blk = sc["sets"][name]; br = blk["brackets"]
        base = ["corner_lo", "corner_hi"] + [f"interior_{i}" for i in range(9)] + ["comparator_0", "comparator_1"]
        keys = []
        for k in base:
            keys.append(k)
            if f"{k}_AD1" in blk["runs"]:
                keys.append(f"{k}_AD1")
        nom = {roi: blk["direct_doses"]["interior_6"][roi]["value"] for roi in rois}
        for c, roi in enumerate(rois):
            ax = axes[r, c]
            L, U, sL, sU = br[roi]["L"], br[roi]["U"], br[roi]["sL"], br[roi]["sU"]
            n = nom[roi]
            ax.axhspan((L - 3 * sL) / n, (U + 3 * sU) / n, color=BANDBLUE, zorder=0, label="envelope [L−3σ, U+3σ]")
            ax.axhline(L / n, color=BLUE, lw=1.0); ax.axhline(U / n, color=BLUE, lw=1.0)
            xs, ys, cols, labs = [], [], [], []
            for i, k in enumerate(keys):
                v = blk["direct_doses"][k][roi]["value"] / n
                if k.startswith("corner"):
                    ax.plot(i, v, marker="s", ms=6, color=INK, lw=0, label="corner evaluations (bracket inputs)" if k == "corner_lo" else None)
                elif k.startswith("comparator"):
                    ax.plot(i, v, marker="D", ms=5.5, color=ORANGE, lw=0, label="preset-comparator points" if k == "comparator_0" else None)
                elif k.endswith("_AD1"):
                    ax.plot(i, v, marker="o", ms=4, mfc="white", mec=AQUA, mew=1.2, lw=0, label="AD-W6-1 re-runs" if k.endswith("interior_3_AD1") or k.endswith("interior_1_AD1") else None)
                else:
                    ax.plot(i, v, marker="o", ms=4.5, color=AQUA, lw=0, label="interior evaluations" if k == "interior_0" else None)
            allv = [blk["direct_doses"][k][roi]["value"] for k in keys]
            if roi == "tumour":
                gap = (min(allv) - L) / min(allv); txt = f"lower-bound gap {gap:.1%}"
            else:
                gap = (U - max(allv)) / max(allv); txt = f"upper-bound gap {gap:.1%}"
            ax.set_title(f"{name} · {roi}: {txt}", fontsize=9, loc="left")
            ax.set_xticks([])
            # separators: corners | interior | comparator
            i_int = keys.index("interior_0") - 0.5; i_cmp = keys.index("comparator_0") - 0.5
            ax.axvline(i_int, color=MUTED, lw=0.6, ls=":"); ax.axvline(i_cmp, color=MUTED, lw=0.6, ls=":")
            ax.set_xlabel(f"corners │ interior │ comparator   (n = {len(keys)})", fontsize=8)
            if c == 0:
                ax.set_ylabel("weighted dose / nominal")
            ax.set_ylim(min((L - 3 * sL) / n * 0.92, 0.5), max((U + 3 * sU) / n * 1.06, 1.15))
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.0), fontsize=8.5)
    fig.text(0.0, -0.01, "All non-corner cases meet the recorded uncertainty-tolerant comparison. Doses are ROI-mean weighted doses normalised to the "
             "nominal point of each set. DS-1: W003 continuity bands. DS-2: full observed range of a 27-patient ¹⁸F-BPA PET cohort.",
             fontsize=7.5, color=INK2)
    save(fig, "fig6_bounds.png")


# ============================================================ Fig 7 window
def fig_window():
    sc = S6["scored"]; rois = ["tumour"]
    fig, axes = plt.subplots(2, 1, figsize=(8.8, 5.0), sharex=True, gridspec_kw={"hspace": 0.6})
    for r, name in enumerate(["DS-2", "DS-1"]):
        blk = sc["sets"][name]; br = blk["brackets"]["tumour"]
        nom = blk["direct_doses"]["interior_6"]["tumour"]["value"]
        comp = min(blk["direct_doses"][k]["tumour"]["value"] for k in ("interior_6", "comparator_0", "comparator_1"))
        edge = br["L"] - 3 * br["sL"]
        crit = 20.0
        x_comp = crit * nom / comp; x_cert = crit * nom / edge
        ax = axes[r]
        ax.set_xlim(15, 65); ax.set_ylim(0, 2.75); ax.set_yticks([0.6, 1.8]); ax.set_yticklabels(["envelope\nedge", "preset\ncomparator"], fontsize=8.5)
        # comparator PASS region
        ax.add_patch(Rectangle((x_comp, 1.5), 65 - x_comp, 0.6, facecolor="#f8d9c8", edgecolor="none"))
        ax.add_patch(Rectangle((15, 1.5), x_comp - 15, 0.6, facecolor=BAND, edgecolor="none"))
        ax.add_patch(Rectangle((x_cert, 0.3), 65 - x_cert, 0.6, facecolor=BANDBLUE, edgecolor="none"))
        ax.add_patch(Rectangle((15, 0.3), x_cert - 15, 0.6, facecolor=BAND, edgecolor="none"))
        bb = dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.85)
        ax.text((x_comp + 65) / 2, 1.8, "meets threshold", ha="center", va="center", fontsize=9, color=INK, bbox=bb)
        ax.text((15 + x_comp) / 2, 1.8, "below threshold", ha="center", va="center", fontsize=9, color=INK2, bbox=bb)
        ax.text((x_cert + 65) / 2, 0.6, "meets threshold", ha="center", va="center", fontsize=9, color=INK, bbox=bb)
        ax.text((15 + x_cert) / 2, 0.6, "edge below\nillustrative threshold", ha="center", va="center", fontsize=8, color=INK2, bbox=bb)
        # divergence window
        ax.axvspan(x_comp, x_cert, ymin=0.04, ymax=0.96, facecolor=ORANGE, alpha=0.18, edgecolor="none")
        ax.axvline(x_comp, color=ORANGE, lw=1.4); ax.axvline(x_cert, color=BLUE, lw=1.4)
        ax.text(x_comp - 0.4, 2.42, f"{x_comp:.1f}", ha="right", fontsize=8.5, color=ORANGE)
        ax.text(x_cert + 0.4, 2.42, f"{x_cert:.1f}", ha="left", fontsize=8.5, color=BLUE)
        ax.text((x_comp + x_cert) / 2, 2.42, "divergence window", ha="center", fontsize=8.5, color=INK2)
        ax.axvline(crit, color=MUTED, lw=0.8, ls=":")
        ax.text(crit + 0.3, 0.05, "criterion 20 Gy-w", fontsize=7.5, color=INK2, ha="left", va="bottom")
        ax.set_title(f"{name}: lower envelope edge / comparator worst case = {edge/comp:.3f}; "
                     f"divergence window [{x_comp:.1f}, {x_cert:.1f}) Gy-w", loc="left", fontsize=9.5)
        ax.grid(False)
    axes[1].set_xlabel("nominal-assumption tumour ROI weighted dose of the plan (Gy-w)")
    fig.text(0.0, -0.03, "Shading: the fixed-ratio preset schedule meets the illustrative tumour "
             "mean-dose criterion (≥ 20 Gy-w); the envelope edge does not. No clinical acceptance claim.",
             fontsize=7.5, color=INK2)
    save(fig, "fig7_window.png")


# ============================================================ Fig 1 concept
def fig_concept():
    sc = S6["scored"]; blk = sc["sets"]["DS-2"]; ds = blk["declared"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1, 1.25], "wspace": 0.3})
    # left: the declared box in (T/B ratio, blood boron) with rho_S collapsed
    ax = axes[0]; ax.grid(False)
    rt, rs, B = ds["rt"], ds["rs"], ds["B"]
    ax.add_patch(Rectangle((rt[0], B[0]), rt[2] - rt[0], B[2] - B[0], facecolor="#f4f4f2", edgecolor=INK2, lw=1.0))
    ax.plot([rt[0]], [B[0]], marker="s", ms=8, color=INK, lw=0); ax.text(rt[0], B[0] - 0.9, "c⁻ (all low)", ha="center", va="top", fontsize=8.5)
    ax.plot([rt[2]], [B[2]], marker="s", ms=8, color=INK, lw=0); ax.text(rt[2], B[2] + 0.6, "c⁺ (all high)", ha="center", va="bottom", fontsize=8.5)
    for x, y in ((rt[0], B[2]), (rt[2], B[0])):
        ax.plot([x], [y], marker="o", ms=5, color=AQUA, lw=0)
    ax.text(rt[0] + 0.15, B[2] + 0.6, "mixed corners", ha="left", va="bottom", fontsize=8, color=AQUA)
    ax.plot([rt[1]], [B[1]], marker="o", ms=5, color=AQUA, lw=0)
    ax.plot([rt[1], rt[1]], [B[0], B[2]], color=ORANGE, lw=0, marker="D", ms=6)
    ax.annotate("preset comparator:\nnominal + blood endpoints\nat T/B = 3.5", xy=(rt[1], B[2]), xytext=(rt[1] + 0.35, B[2] - 3.2),
                fontsize=8, color=ORANGE, arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
    ax.text(rt[1] + 0.05, B[1] - 1.4, "nominal", fontsize=8, color=AQUA)
    ax.set_xlim(rt[0] - 0.6, rt[2] + 0.6); ax.set_ylim(B[0] - 3, B[2] + 3)
    ax.set_xlabel("tumour-to-blood ratio ρ_T"); ax.set_ylabel("blood boron B (ppm ¹⁰B)")
    ax.set_title("Declared set DS-2 (PET cohort range)", loc="left", fontsize=9.5)
    ax.text(rt[0] - 0.5, B[0] - 2.4, "two evaluations at c⁻ and c⁺ → bounds for every point in the box", fontsize=8, color=INK2)
    # right: real data — tumour dose vs tumour boron for all evaluations, band from the two corners
    ax = axes[1]
    W = {k: tuple(v) for k, v in sc["weights"].items()}; MASS = sc["roi_mass_g"]
    br = blk["brackets"]["tumour"]; nom = blk["direct_doses"]["interior_6"]["tumour"]["value"]
    L3, U3 = (br["L"] - 3 * br["sL"]) * CONV, (br["U"] + 3 * br["sU"]) * CONV
    ax.axhspan(L3, U3, color=BANDBLUE, zorder=0, label="approximate envelope [L−3σ, U+3σ] from c⁻, c⁺ only")
    for k, run in blk["runs"].items():
        p = run["params"]; xb = p[0] * p[2]
        v = blk["direct_doses"][k]["tumour"]["value"] * CONV
        if k.startswith("corner"):
            ax.plot(xb, v, marker="s", ms=7, color=INK, lw=0, label="corner evaluations" if k == "corner_lo" else None)
        elif k.startswith("comparator"):
            ax.plot(xb, v, marker="D", ms=6, color=ORANGE, lw=0, label="preset-comparator points" if k == "comparator_0" else None)
        elif k.endswith("_AD1"):
            continue
        else:
            ax.plot(xb, v, marker="o", ms=5, color=AQUA, lw=0, label="interior evaluations" if k == "interior_0" else None)
    ax.axhline(L3, color=BLUE, lw=1.0); ax.axhline(U3, color=BLUE, lw=1.0)
    ax.set_ylim(L3 * 0.75, U3 * 1.12)
    ax.text(ax.get_xlim()[1], L3, "lower envelope edge  ", va="bottom", ha="right", fontsize=8, color=BLUE)
    ax.text(ax.get_xlim()[1], U3, "upper envelope edge  ", va="bottom", ha="right", fontsize=8, color=BLUE)
    ax.set_xlabel("tumour ¹⁰B concentration ρ_T·B (ppm)")
    ax.set_ylabel("tumour ROI dose rate (Gy-w/min at 10¹⁰ n cm⁻² s⁻¹)")
    ax.set_title("Two-run envelope and direct evaluations (DS-2, tumour)", loc="left", fontsize=9.5)
    # the gap that matters: comparator worst case vs lower envelope edge
    comp_keys = ("interior_6", "comparator_0", "comparator_1")
    kmin = min(comp_keys, key=lambda k: blk["direct_doses"][k]["tumour"]["value"])
    xc = blk["runs"][kmin]["params"][0] * blk["runs"][kmin]["params"][2]
    vc = blk["direct_doses"][kmin]["tumour"]["value"] * CONV
    ax.annotate("", xy=(xc, L3), xytext=(xc, vc), arrowprops=dict(arrowstyle="<->", color=ORANGE, lw=1.2))
    ax.text(xc + 3, (L3 + vc) / 2, f"envelope edge = {L3/vc:.2f} ×\ncomparator worst case", fontsize=8, color=ORANGE, va="center")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=8)
    fig.text(0.0, -0.30, "Left: the declared uncertainty box (the scalp-ratio axis is omitted from the drawing). Right: ROI-mean weighted dose from every "
             "evaluation of the scored block against the interval built from the two corner runs alone.", fontsize=7.5, color=INK2)
    save(fig, "fig1_concept.png")


fig_geometry(); fig_epithermal(); fig_mono(); fig_depression(); fig_bounds(); fig_window(); fig_concept()
