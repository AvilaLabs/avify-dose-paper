#!/usr/bin/env python3
"""W006 harness — reference study on the modified Snyder head phantom.

Implements the FROZEN protocol `W006_PROTOCOL.md`
(sha256 c9558f2896d3e1330cfc8130d0a0722e6c08bd63af5334832931c09c8dd9153d):
the transcribed analytic Snyder phantom (§2.1), ICRU 46 materials (§2.2),
the generic epithermal beam and the monoenergetic control beams (§2.3), the
16 x 16 x 4 mm central-axis tally grid (§2.4), the reference kerma factors
(§2.5), the declared tumour (§2.6), declared sets (§2.7), weighting (§2.8),
frozen tallies (§3), the W003 bracket construction (§4), and the frozen run
schedule (§5). Gates and grade are NOT computed here in any authoritative
sense: `check_w006.py` re-derives everything from the raw tallies written
to `controls/w006_results.json` and derives the grade.

Stages (frozen order, §5):
    w006_run.py c1        C1 setup control (determinism, direction, zero-boron)
    w006_run.py c2        C2 monoenergetic beams (6 neutron + 6 photon)
    w006_run.py c3        C3 flux-depression control (30 ug/g B-10 scalp+brain)
    w006_run.py configA   Configuration A validation profiles (group-split)
    w006_run.py scored    Configuration B scored block + preset comparator
    w006_run.py assemble  merge stage files into controls/w006_results.json
    w006_run.py bench     timing only; writes to the scratchpad, never to controls/

Each stage writes controls/w006_<stage>.json (raw tallies, seeds, histories,
atom densities); `assemble` merges them. Run inputs and statepoints go to
w006_runs/<tag>/ (bulk artifacts, excluded from the manifest as in W003).

Transcription notes recorded at build time (F-W6-6, tooling scope; also in
ledger.md before any run):
  T1. The reference deck's fast source group (`si7`) spans 10 keV to 2 MeV
      and the paper's own prose says "10 keV to 2 MeV". The parenthetical
      "> 100 keV" in protocol §2.3 came from the deck's header comment; the
      deck's `si` tables are what the reference actually sampled and are what
      is built here. Physics is unchanged from the reference.
  T2. The reference deck's air material lists ZAID 20000 (natural calcium)
      for the 1.283 % "argon" fraction. Protocol §2.2 transcribes it as Ar,
      which is the physical composition; Ar is used here. At 1.293e-3 g/cm3
      the choice is immaterial to any tally.
  T3. MCNP DE/DF dose functions interpolate log-log by default and hold the
      endpoint value outside the table; the OpenMC EnergyFunctionFilter is
      run log-log and the neutron tables are extended flat to 1e-5 eV so the
      endpoint rule is reproduced. Photon tables start at 1 keV, which is
      also the transport cutoff.
  T4. The reference's source-energy biasing (sb2 0.1/0.2/0.7) is reproduced
      as three separately-run source groups (thermal 10 %, epithermal 89 %,
      fast 1 %) combined linearly with their physical fractions; this is the
      same estimator as biased sampling and yields group-resolved profiles.
      Configuration B and C1 use the analog composite source.
  T5. The OpenMC `heating` score for PHOTONS (collision/analog estimator;
      tracklength is refused) was measured at build time to report ~1e-4 of
      the physical photon energy deposition (4.7 eV per source neutron for
      the whole phantom against 5.6e5 eV of capture-photon energy produced,
      from `heating-local` minus `heating` for neutrons). The ROI photon
      dose (frozen tally T5) is therefore scored as photon flux x NIST
      photon kerma for brain — the estimator the reference itself uses
      (F44) and the one used for the profile tallies — converted to eV
      deposited in the ROI so the W003 dose formula is unchanged. W003 used
      the photon `heating` score; its photon component was numerically
      absent as a consequence (disclosed in the ledger; W003 gate logic is
      unaffected, its tightness gaps would only shrink with the photon
      dose present). Neutron `heating` (MT301) is unaffected and is used
      for the fast component as in W003.
  T6. Run order: C1, C2, Configuration A, C3, [scored]. C3's comparison
      baseline is the Configuration A no-boron profile, so A precedes C3;
      the protocol's constraint that component validation precede any
      composite-spectrum run is honoured (C2 precedes A and C3).
"""
import glob
import hashlib
import json
import math
import os
import shutil
import sys
import time

import numpy as np
import openmc

DATA = os.path.expanduser("~/nuclear-data/endfb-vii.1-hdf5/cross_sections.xml")
os.environ.setdefault("OPENMC_CROSS_SECTIONS", DATA)
openmc.config["cross_sections"] = DATA
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
# Record locations. The environment overrides exist ONLY for pipeline
# testing at reduced scale in a scratch location; the record is produced with
# none of them set (OUTDIR = controls/, RUNROOT = w006_runs/, SCALE = 1).
RUNROOT = os.environ.get("W006_RUNROOT", os.path.join(ROOT, "w006_runs"))
OUTDIR = os.environ.get("W006_OUT", HERE)
SCRATCH = os.environ.get(
    "W006_SCRATCH",
    "/tmp/claude-1000/-home-connoravila-Documents-Project-Waddle/"
    "8ef62f76-38f7-49e5-ab82-5eed9b527ee6/scratchpad/w006_bench")
RESULTS = os.path.join(OUTDIR, "w006_results.json")
REFERENCE = os.path.join(HERE, "w006_reference.json")
HISTORY_SCALE = float(os.environ.get("W006_HISTORY_SCALE", "1"))
PROTOCOL_SHA256 = "c9558f2896d3e1330cfc8130d0a0722e6c08bd63af5334832931c09c8dd9153d"
ARCHIVE_SHA256 = "511aa05c4253c77b2778a2bd7aa52334a1fd996cca3d08b7c1899f32db38fd62"
THREADS = int(os.environ.get("W006_THREADS", "8"))
OPENMC_EXEC = os.path.join(sys.prefix, "bin", "openmc")   # ~/.venvs/w003env/bin/openmc
if not os.path.exists(OPENMC_EXEC):
    OPENMC_EXEC = "openmc"


def stage_file(stage):
    return os.path.join(OUTDIR, f"w006_{stage}.json")


def scaled(n):
    """History count under the (test-only) scale factor; >= 100 per batch."""
    return max(int(n * HISTORY_SCALE), 100 * BATCHES)


# =====================================================================
# Frozen physical model — W006_PROTOCOL.md §2 (transcribed values)
# =====================================================================

# §2.1 MCNP `sq` coefficients (A, B, C) and z-centre; the coefficients are
# authoritative, semi-axes 1/sqrt(A,B,C) = (6, 9, 6.5), (6.8, 9.8, 8.3),
# (7.3, 10.3, 8.8) cm.
SQ = {"brain": (0.0277778, 0.0123456, 0.0236686, 1.0),
      "cranium": (0.0216263, 0.0104123, 0.0145159, 0.0),
      "scalp": (0.0187652, 0.0094260, 0.0129133, 0.0)}
SEMI = {"brain": (6.0, 9.0, 6.5), "cranium": (6.8, 9.8, 8.3), "scalp": (7.3, 10.3, 8.8)}
R_WORLD = 50.0            # deck surface 40: `so 50`, vacuum beyond
Z_SURFACE = 8.8           # scalp/air boundary on +z axis (depth zero)

# §2.2 ICRU 46 mass fractions; keys with a mass number are nuclides, bare
# symbols are natural elements expanded to their ENDF/B-VII.1 isotopes (deck
# ZAIDs 12000, 16000, 17000, 19000, 20000 -> natural Mg, S, Cl, K, Ca; deck
# 6000 -> the library's natural-carbon evaluation C0). Air argon: note T2.
DENSITY = {"brain": 1.040, "cranium": 1.610, "scalp": 1.090, "air": 1.293e-3}
COMPOSITION = {
    "brain": {"H1": 0.107, "C0": 0.145, "N14": 0.022, "O16": 0.712, "Na23": 0.002,
              "P31": 0.004, "S": 0.002, "Cl": 0.003, "K": 0.003},
    "cranium": {"H1": 0.050, "C0": 0.212, "N14": 0.040, "O16": 0.435, "Na23": 0.001,
                "Mg": 0.002, "P31": 0.081, "S": 0.003, "Ca": 0.176},
    "scalp": {"H1": 0.100, "C0": 0.204, "N14": 0.042, "O16": 0.645, "Na23": 0.002,
              "P31": 0.001, "S": 0.002, "Cl": 0.003, "K": 0.001},
    "air": {"C0": 0.00012, "N14": 0.75527, "O16": 0.23178, "Ar": 0.01283},
}
SAB = "c_H_in_H2O"        # deck `lwtr.01t` on hydrogen in the three tissues

# §2.3 beam: monodirectional disk, z = +15 cm, direction (0,0,-1), r in [0,5]
BEAM_Z, BEAM_R, BEAM_DIR = 15.0, 5.0, (0.0, 0.0, -1.0)
BEAM_AREA = math.pi * BEAM_R**2                    # 78.53982 cm2
GROUP_FRACTIONS = {"thermal": 0.10, "epithermal": 0.89, "fast": 0.01}
MONO_NEUTRON_EV = {"css253": 0.0253, "css1": 1.0e3, "css2": 2.0e3, "css10": 1.0e4,
                   "css100": 1.0e5, "css1M": 1.0e6}
MONO_PHOTON_EV = {"cssp02": 0.2e6, "cssp05": 0.5e6, "cssp1": 1.0e6, "cssp2": 2.0e6,
                  "cssp5": 5.0e6, "cssp10": 10.0e6}

# §2.4 tally grid: 16 x 16 mm cross-section, 4 mm thick, 46 bins z=-8.8..9.6
MESH_XY, MESH_Z0, MESH_Z1, MESH_NZ = 0.8, -8.8, 9.6, 46
BIN_VOLUME = (2 * MESH_XY) ** 2 * (MESH_Z1 - MESH_Z0) / MESH_NZ   # 1.024 cm3
DEPTHS = [round(Z_SURFACE - (MESH_Z0 + (k + 0.5) * 0.4), 4) for k in range(MESH_NZ)]
# DEPTHS[0] = 17.4 (deepest), DEPTHS[45] = -0.6 (air, upstream)

# §2.5 normalisation used by the reference tallies: 1e10 n/cm2/s over the
# 5 cm disk = 7.853982e11 n/s = 4.71239e13 n/min. Neutron/photon brain
# kermas are in cGy cm2 (divide by 100), boron kerma in Gy cm2 per 1 ppm.
SOURCE_PER_MIN = 4.71239e13
GY_PER_MEV_PER_G = 1.602176634e-10   # 1 MeV/g = 1.602e-13 J / 1e-3 kg

# §2.6 declared tumour (this project's extension; Configuration B only)
TUMOUR_R, TUMOUR_Z = 1.5, 3.8

# §2.7 declared sets. Axes (rho_T, rho_S, B): rho_T = tumour/blood ratio,
# rho_S = the W003 rho_M band applied to the scalp (the surface normal-tissue
# ROI, as mucosa was in W003), B = blood boron in ppm. Brain follows blood at
# BRAIN_RATIO (fixed 1.0, the W003 convention). The protocol states that
# normal-tissue ratios are "declared from cited literature" but records no
# values; this continuity reading was recorded in the ledger as an OPEN item
# and ADOPTED by the principal 2026-08-17. It affects Configuration B only.
BRAIN_RATIO = 1.0
DS = {
    "DS-1": {"rt": (2.5, 3.5, 4.5), "rs": (0.8, 1.0, 1.2), "B": (15.0, 20.0, 25.0)},
    "DS-2": {"rt": (1.58, 3.5, 5.86), "rs": (0.8, 1.0, 1.2), "B": (19.36, 25.0, 31.55)},
}

# §2.8 weighting (w_B, w_N, w_f, w_g)
W = {"tumour": (3.8, 3.2, 3.2, 1.0), "brain": (1.35, 3.2, 3.2, 1.0),
     "scalp": (2.5, 3.2, 3.2, 1.0)}
E_B, E_N = 2.34, 0.626    # MeV per 10B(n,a) (charged, local) and per 14N(n,p)
EPS_F = 0.02              # declared fast-component remainder (W003 §4)
BORON_GAMMA_BRANCH = 0.94

# §5 C3 control: 30 ug/g 10B in scalp and brain, not cranium
C3_PPM = 30.0

# §6.1 criteria (Gy-w) — evaluation requires a plan normalisation, see G4 note
CRITERIA = {"tumour": (">=", 20.0), "brain": ("<=", 11.0), "scalp": ("<=", 10.0)}
G1_TOL = 0.10
G1C_BAND_PCT = (8.0, 11.0)
G1C_PEAK_DEPTHS = (2.2, 2.6, 3.0, 3.4, 3.8)     # published grid points in 2-4 cm
GAP_KILL = 0.35
AD_W6_1_RELERR = 0.025

# History counts, set from the 2026-08-17 benchmark (~50k histories/s on
# 8 threads for neutron+photon runs) and recorded in the ledger before any
# run. Configuration A: chosen so that the combined relative error is
# <= 1 % at all 184 published points (§5 item 4); the deepest bins drive
# the epithermal and fast group counts. Configuration B: chosen so that
# every scored ROI tally is expected below the AD-W6-1 threshold (2.5 %) —
# the tumour fast-heat tally is the limiting one (~15 % at 4e5 histories).
N_PC = 2_000_000
N_C2_N = {"css253": 130_000_000, "css1": 60_000_000, "css2": 60_000_000,
          "css10": 60_000_000, "css100": 60_000_000, "css1M": 60_000_000}
N_C2_P = 20_000_000       # per monoenergetic photon beam
N_A = {"thermal": 30_000_000, "epithermal": 600_000_000, "fast": 150_000_000}
N_C3 = {"thermal": 10_000_000, "epithermal": 150_000_000, "fast": 40_000_000}
N_FULL, N_INT = 40_000_000, 20_000_000
BATCHES = 20                       # as W003; large profile runs use more
PARTICLES_PER_BATCH_MAX = 5_000_000   # caps the fixed-source bank memory


def n_batches(n):
    return max(BATCHES, int(math.ceil(n / PARTICLES_PER_BATCH_MAX)))

# ---------------------------------------------------------------------
# Kerma tables transcribed from the reference archive member
# MCNPinput/csse (de/df cards): boron per 1 ppm 10B (Gy cm2), ICRU 46 adult
# brain neutron kerma (cGy cm2; ICRU 63/ENDF/B-VI, S and Cl from JENDL-3.2,
# extrapolated 1/v below 0.0253 eV to 1e-4 eV as in the deck), and NIST
# photon kerma for brain (cGy cm2). Energies in MeV as in the deck.
# ---------------------------------------------------------------------
KERMA_B10_PER_PPM_GYCM2 = (
    (1e-11, 4.36e-12), (2.53e-08, 8.666e-14), (9.4e-06, 4.49e-15),
    (0.00015, 1.12e-15), (0.00025, 8.663e-16), (0.00035, 7.314e-16),
    (0.00045, 6.444e-16), (0.00055, 5.824e-16), (0.00065, 5.355e-16),
    (0.00075, 4.982e-16), (0.00085, 4.677e-16), (0.00095, 4.42e-16),
    (0.0015, 3.509e-16), (0.0025, 2.708e-16), (0.0035, 2.284e-16),
    (0.0045, 2.011e-16), (0.0055, 1.816e-16), (0.0065, 1.669e-16),
    (0.0075, 1.552e-16), (0.0085, 1.457e-16), (0.0095, 1.377e-16),
    (0.015, 1.095e-16), (0.02, 9.493e-17), (0.024, 8.681e-17),
    (0.03, 7.797e-17), (0.045, 6.462e-17), (0.055, 5.918e-17),
    (0.065, 5.519e-17), (0.075, 5.212e-17), (0.085, 4.968e-17),
    (0.095, 4.77e-17), (0.1, 4.68e-17), (0.12, 4.388e-17), (0.15, 4.039e-17),
    (0.17, 3.828e-17), (0.18, 3.728e-17), (0.19, 3.628e-17),
    (0.2, 3.532e-17), (0.21, 3.434e-17), (0.22, 3.339e-17),
    (0.23, 3.247e-17), (0.235, 3.202e-17), (0.24, 3.157e-17),
    (0.245, 3.113e-17), (0.25, 3.071e-17), (0.26, 2.986e-17),
    (0.27, 2.905e-17), (0.28, 2.827e-17), (0.3, 2.677e-17),
    (0.32, 2.547e-17), (0.34, 2.439e-17), (0.36, 2.359e-17),
    (0.38, 2.311e-17), (0.4, 2.295e-17), (0.42, 2.307e-17),
    (0.44, 2.332e-17), (0.46, 2.347e-17), (0.48, 2.325e-17),
    (0.5, 2.255e-17), (0.52, 2.145e-17), (0.54, 2.011e-17),
    (0.56, 1.873e-17), (0.58, 1.74e-17), (0.6, 1.62e-17), (0.65, 1.375e-17),
    (0.7, 1.199e-17), (0.75, 1.069e-17), (0.8, 9.721e-18), (0.85, 8.974e-18),
    (0.9, 8.393e-18), (0.95, 7.936e-18), (1.0, 7.577e-18), (1.05, 7.396e-18),
    (1.1, 7.409e-18), (1.15, 7.56e-18), (1.2, 7.81e-18), (1.25, 8.13e-18),
    (1.3, 8.557e-18), (1.35, 9.139e-18), (1.4, 9.864e-18), (1.45, 1.073e-17),
    (1.5, 1.172e-17), (1.55, 1.317e-17), (1.6, 1.527e-17), (1.65, 1.774e-17),
    (1.7, 2.015e-17), (1.75, 2.197e-17), (1.8, 2.267e-17), (1.85, 2.274e-17),
    (1.9, 2.231e-17), (1.95, 2.16e-17), (2.0, 2.069e-17), (2.05, 1.946e-17),
    (2.1, 1.817e-17), (2.15, 1.69e-17), (2.2, 1.582e-17), (2.25, 1.506e-17),
    (2.3, 1.451e-17), (2.35, 1.408e-17), (2.4, 1.376e-17), (2.45, 1.369e-17),
    (2.5, 1.415e-17), (2.55, 1.508e-17), (2.6, 1.607e-17), (2.65, 1.703e-17),
    (2.7, 1.802e-17), (2.75, 1.977e-17), (2.8, 2.095e-17), (2.85, 2.034e-17),
    (2.9, 1.948e-17), (2.95, 1.869e-17), (3.0, 1.795e-17), (3.05, 1.731e-17),
    (3.1, 1.678e-17), (3.15, 1.624e-17), (3.2, 1.585e-17), (3.3, 1.613e-17),
    (3.4, 1.721e-17), (3.6, 1.981e-17), (3.7, 2.119e-17), (3.8, 2.26e-17),
    (4.0, 2.23e-17), (4.2, 1.791e-17), (4.3, 1.724e-17), (4.4, 1.643e-17),
    (4.6, 1.377e-17), (4.8, 1.212e-17), (5.0, 1.076e-17), (5.2, 9.885e-18),
    (5.4, 8.972e-18), (5.5, 8.807e-18), (5.6, 8.892e-18), (5.8, 9.318e-18),
    (6.0, 9.685e-18), (6.2, 1.005e-17), (6.4, 1.03e-17), (6.6, 1.043e-17),
    (6.8, 1.045e-17), (7.0, 1.027e-17), (7.2, 9.887e-18), (7.4, 9.427e-18),
    (7.6, 8.954e-18), (7.8, 8.516e-18), (8.0, 8.143e-18), (8.2, 7.795e-18),
    (8.4, 7.495e-18), (8.6, 7.234e-18), (8.8, 7.004e-18), (9.0, 6.796e-18),
    (9.2, 6.609e-18), (9.4, 6.444e-18), (9.6, 6.303e-18), (9.8, 6.185e-18),
    (10.0, 6.091e-18), (10.5, 5.856e-18), (11.0, 5.74e-18),
    (11.5, 5.695e-18), (12.0, 5.764e-18), (12.5, 5.963e-18),
    (13.0, 6.259e-18), (13.5, 6.641e-18), (14.0, 7.032e-18),
    (14.5, 7.329e-18), (15.0, 7.492e-18), (15.5, 7.462e-18), (16.0, 7.4e-18),
    (16.5, 7.329e-18), (17.0, 7.25e-18), (17.5, 7.155e-18),
    (18.0, 7.047e-18), (18.5, 6.926e-18), (19.0, 6.792e-18),
    (19.5, 6.646e-18), (20.0, 6.489e-18),
)

# de24/df24 (E <= 0.5 eV) and de34/df34 (E >= 0.5 eV) merged; the deck's
# 1e-50 sentinels that implement the 0.5 eV cut are replaced by explicit
# EnergyFilter bins (thermal < 0.5 eV, fast > 0.5 eV) on the same table.
KERMA_BRAIN_N_CGYCM2 = (
    (1e-10, 2.84836e-10), (2.53e-08, 1.79156e-11), (3.6e-08, 1.51863e-11),
    (6.3e-08, 1.14605e-11), (1.1e-07, 8.68979e-12), (2e-07, 6.46317e-12),
    (3.6e-07, 4.80674e-12), (5e-07, 4.0973e-12), (6.3e-07, 3.65016e-12),
    (1.1e-06, 2.7692e-12), (2e-06, 2.06768e-12), (3.6e-06, 1.56446e-12),
    (6.3e-06, 1.22285e-12), (1.1e-05, 9.93368e-13), (2e-05, 8.67077e-13),
    (3.6e-05, 8.73031e-13), (6.3e-05, 1.04735e-12), (0.00011, 1.45912e-12),
    (0.0002, 2.35827e-12), (0.000236228, 2.7312e-12),
    (0.000285955, 3.24754e-12), (0.000317034, 3.57212e-12),
    (0.000341898, 3.83283e-12), (0.000357438, 3.99667e-12),
    (0.00036, 4.02383e-12), (0.000369869, 4.12962e-12),
    (0.000376085, 4.19786e-12), (0.000382301, 4.27008e-12),
    (0.000386186, 4.32166e-12), (0.000389294, 4.3751e-12),
    (0.000392013, 4.45426e-12), (0.000393567, 4.54948e-12),
    (0.000394927, 4.75157e-12), (0.000395704, 5.04322e-12),
    (0.000396384, 5.69169e-12), (0.000396772, 6.59789e-12),
    (0.000397161, 8.87168e-12), (0.000397404, 1.25184e-11),
    (0.000397598, 1.92702e-11), (0.000397841, 4.02373e-11),
    (0.00039792, 4.89985e-11), (0.00039796, 5.18963e-11),
    (0.000398, 5.29405e-11), (0.00039804, 5.18752e-11),
    (0.00039808, 4.89616e-11), (0.000398159, 4.01935e-11),
    (0.000398277, 2.78043e-11), (0.000398512, 1.48024e-11),
    (0.000398747, 9.92659e-12), (0.000398982, 7.77267e-12),
    (0.000399453, 6.02491e-12), (0.000399923, 5.3616e-12),
    (0.000400746, 4.90651e-12), (0.000401687, 4.71488e-12),
    (0.000403334, 4.60009e-12), (0.000405215, 4.56478e-12),
    (0.000408508, 4.56428e-12), (0.000412271, 4.58953e-12),
    (0.000416975, 4.6317e-12), (0.000424501, 4.70637e-12),
    (0.000435791, 4.82298e-12), (0.000450844, 4.98071e-12),
    (0.000473423, 5.21862e-12), (0.000503528, 5.5366e-12),
    (0.000548686, 6.01428e-12), (0.000608897, 6.65192e-12),
    (0.00063, 6.87555e-12), (0.0011, 1.18485e-11), (0.002, 2.12787e-11),
    (0.00220448, 2.34081e-11), (0.00292701, 3.09073e-11),
    (0.0034087, 3.58891e-11), (0.0036, 3.78642e-11),
    (0.00364954, 3.8369e-11), (0.00383018, 4.02082e-11),
    (0.0039506, 4.14327e-11), (0.00402586, 4.21975e-11),
    (0.00408607, 4.28091e-11), (0.00413123, 4.32676e-11),
    (0.00416134, 4.35733e-11), (0.00418768, 4.38408e-11),
    (0.00420273, 4.39938e-11), (0.00421778, 4.41469e-11),
    (0.00422531, 4.42245e-11), (0.00423284, 4.43023e-11),
    (0.00423848, 4.4363e-11), (0.00424224, 4.44084e-11),
    (0.00424507, 4.44503e-11), (0.00424695, 4.44963e-11),
    (0.00424836, 4.45582e-11), (0.0042493, 4.46597e-11),
    (0.00425024, 4.48206e-11), (0.00425071, 4.52448e-11),
    (0.00425118, 4.57924e-11), (0.00425165, 4.70518e-11),
    (0.00425183, 4.99682e-11), (0.00425191, 5.12363e-11),
    (0.004252, 5.16511e-11), (0.00425209, 5.18013e-11),
    (0.00425217, 5.16522e-11), (0.00425235, 5.12397e-11),
    (0.00425262, 4.99748e-11), (0.00425299, 4.80241e-11),
    (0.00425346, 4.64707e-11), (0.00425419, 4.55827e-11),
    (0.00425493, 4.50345e-11), (0.00425641, 4.48334e-11),
    (0.00425788, 4.46864e-11), (0.00426047, 4.46546e-11),
    (0.00426342, 4.46455e-11), (0.00426859, 4.46815e-11),
    (0.00427449, 4.47307e-11), (0.00428188, 4.48012e-11),
    (0.00429369, 4.49187e-11), (0.0043055, 4.5037e-11),
    (0.00432912, 4.52759e-11), (0.00435274, 4.55148e-11),
    (0.00439408, 4.59335e-11), (0.00444133, 4.64121e-11),
    (0.00451219, 4.71296e-11), (0.00460668, 4.8086e-11),
    (0.00477204, 4.97582e-11), (0.00496102, 5.16672e-11),
    (0.00529173, 5.50029e-11), (0.00576417, 5.97574e-11),
    (0.0063, 6.51359e-11), (0.011, 1.10066e-10), (0.02, 1.90608e-10),
    (0.036, 3.15707e-10), (0.063, 4.88476e-10), (0.082, 5.93968e-10),
    (0.086, 6.13173e-10), (0.09, 6.32378e-10), (0.094, 6.52636e-10),
    (0.098, 6.70771e-10), (0.105, 7.01356e-10), (0.115, 7.44549e-10),
    (0.125, 7.85529e-10), (0.135, 8.23295e-10), (0.145, 8.60107e-10),
    (0.155, 8.97e-10), (0.165, 9.30271e-10), (0.175, 9.62687e-10),
    (0.185, 9.95179e-10), (0.195, 1.02648e-09), (0.21, 1.06971e-09),
    (0.23, 1.12623e-09), (0.25, 1.17659e-09), (0.27, 1.225e-09),
    (0.29, 1.28468e-09), (0.31, 1.32372e-09), (0.33, 1.37525e-09),
    (0.35, 1.42928e-09), (0.37, 1.47882e-09), (0.39, 1.54217e-09),
    (0.42, 1.69417e-09), (0.46, 1.66227e-09), (0.5, 1.66729e-09),
    (0.54, 1.72979e-09), (0.58, 1.78877e-09), (0.62, 1.85065e-09),
    (0.66, 1.91171e-09), (0.7, 1.96912e-09), (0.74, 2.01813e-09),
    (0.78, 2.07804e-09), (0.82, 2.13004e-09), (0.86, 2.18651e-09),
    (0.9, 2.25552e-09), (0.94, 2.36364e-09), (0.98, 2.54054e-09),
    (1.05, 2.55888e-09), (1.15, 2.5495e-09), (1.25, 2.66336e-09),
    (1.35, 2.73723e-09), (1.45, 2.78737e-09), (1.55, 2.86468e-09),
    (1.65, 2.97877e-09), (1.75, 3.02627e-09), (1.85, 3.15611e-09),
    (1.95, 3.16622e-09), (2.1, 3.26881e-09), (2.3, 3.31828e-09),
    (2.5, 3.4701e-09), (2.7, 3.61684e-09), (2.9, 3.75627e-09),
    (3.1, 3.8782e-09), (3.3, 4.21016e-09), (3.5, 4.30374e-09),
    (3.7, 4.41548e-09), (3.9, 4.34663e-09), (4.2, 4.46838e-09),
    (4.6, 4.47802e-09), (5.0, 4.75189e-09), (5.4, 4.62536e-09),
    (5.8, 4.83295e-09), (6.2, 4.96441e-09), (6.6, 5.15965e-09),
    (7.0, 5.28666e-09), (7.4, 5.50396e-09), (7.8, 5.50521e-09),
    (8.2, 5.48272e-09), (8.6, 5.65167e-09), (9.0, 5.71392e-09),
    (9.4, 5.86699e-09), (9.8, 6.03676e-09), (10.5, 6.08986e-09),
    (11.5, 6.4657e-09), (12.5, 6.41512e-09), (13.5, 6.58102e-09),
    (14.5, 6.62826e-09), (16.0, 6.77503e-09), (18.0, 6.93483e-09),
    (20.0, 7.03169e-09),
)

KERMA_BRAIN_P_CGYCM2 = (
    (0.001, 5.9056e-08), (0.00103542, 5.579e-08), (0.0010721, 5.2682e-08),
    (0.0015, 2.9873e-08), (0.002, 1.7807e-08), (0.0021455, 1.564e-08),
    (0.00230297, 1.3973e-08), (0.002472, 1.2242e-08),
    (0.0026414, 1.0935e-08), (0.0028224, 9.6499e-09), (0.003, 8.7671e-09),
    (0.0036074, 6.1496e-09), (0.004, 5.1641e-09), (0.005, 3.3261e-09),
    (0.006, 2.3071e-09), (0.008, 1.2802e-09), (0.01, 8.0381e-10),
    (0.015, 3.3934e-10), (0.02, 1.8284e-10), (0.03, 7.8298e-11),
    (0.04, 4.6623e-11), (0.05, 3.52e-11), (0.06, 3.1598e-11),
    (0.08, 3.3735e-11), (0.1, 4.0984e-11), (0.15, 6.6258e-11),
    (0.2, 9.4689e-11), (0.3, 1.5275e-10), (0.4, 2.0912e-10), (0.5, 2.63e-10),
    (0.6, 3.1415e-10), (0.8, 4.09e-10), (1.0, 4.9475e-10), (1.25, 5.91e-10),
    (1.5, 6.7772e-10), (2.0, 8.3153e-10), (3.0, 1.0906e-09),
    (4.0, 1.3157e-09), (5.0, 1.5245e-09), (6.0, 1.7236e-09),
    (8.0, 2.1085e-09), (10.0, 2.4866e-09), (15.0, 3.4247e-09),
    (20.0, 4.3707e-09),
)

# Source energy grids (MeV) transcribed from the deck's si5/si6/si7 cards:
# 100 equal-lethargy bins per group, equal probability per bin, uniform in
# energy within a bin (MCNP `sp D 0 1 99r`).
SI_THERMAL_MEV = (
    1e-09, 1.06412e-09, 1.13235e-09, 1.20495e-09, 1.28221e-09, 1.36442e-09,
    1.4519e-09, 1.545e-09, 1.64406e-09, 1.74947e-09, 1.86165e-09,
    1.98101e-09, 2.10803e-09, 2.24319e-09, 2.38702e-09, 2.54007e-09,
    2.70293e-09, 2.87624e-09, 3.06066e-09, 3.2569e-09, 3.46572e-09,
    3.68794e-09, 3.9244e-09, 4.17603e-09, 4.44378e-09, 4.72871e-09,
    5.0319e-09, 5.35454e-09, 5.69786e-09, 6.06319e-09, 6.45195e-09,
    6.86563e-09, 7.30584e-09, 7.77428e-09, 8.27275e-09, 8.80318e-09,
    9.36762e-09, 9.96825e-09, 1.06074e-08, 1.12875e-08, 1.20112e-08,
    1.27814e-08, 1.36009e-08, 1.4473e-08, 1.54009e-08, 1.63884e-08,
    1.74392e-08, 1.85573e-08, 1.97472e-08, 2.10134e-08, 2.23607e-08,
    2.37944e-08, 2.532e-08, 2.69435e-08, 2.86711e-08, 3.05094e-08,
    3.24656e-08, 3.45472e-08, 3.67623e-08, 3.91194e-08, 4.16277e-08,
    4.42967e-08, 4.71369e-08, 5.01593e-08, 5.33754e-08, 5.67977e-08,
    6.04394e-08, 6.43146e-08, 6.84384e-08, 7.28265e-08, 7.74959e-08,
    8.24648e-08, 8.77523e-08, 9.33788e-08, 9.9366e-08, 1.05737e-07,
    1.12517e-07, 1.19731e-07, 1.27408e-07, 1.35577e-07, 1.4427e-07,
    1.5352e-07, 1.63364e-07, 1.73838e-07, 1.84984e-07, 1.96845e-07,
    2.09466e-07, 2.22897e-07, 2.37188e-07, 2.52396e-07, 2.6858e-07,
    2.858e-07, 3.04125e-07, 3.23625e-07, 3.44375e-07, 3.66456e-07,
    3.89952e-07, 4.14955e-07, 4.41561e-07, 4.69873e-07, 5e-07,
)
SI_EPI_MEV = (
    5e-07, 5.520524e-07, 6.0952371e-07, 6.7297806e-07, 7.4303831e-07,
    8.2039217e-07, 9.0579893e-07, 1.000097e-06, 1.1042119e-06, 1.2191656e-06,
    1.3460866e-06, 1.4862207e-06, 1.6409434e-06, 1.8117735e-06,
    2.0003878e-06, 2.2086378e-06, 2.4385676e-06, 2.6924342e-06,
    2.9727296e-06, 3.282205e-06, 3.6238983e-06, 4.0011636e-06, 4.4177039e-06,
    4.8776081e-06, 5.3853906e-06, 5.9460356e-06, 6.5650465e-06,
    7.2484994e-06, 8.003103e-06, 8.8362645e-06, 9.7561621e-06, 1.0771825e-05,
    1.1893224e-05, 1.3131366e-05, 1.4498404e-05, 1.6007758e-05,
    1.7674242e-05, 1.9514216e-05, 2.154574e-05, 2.3788755e-05, 2.6265278e-05,
    2.899962e-05, 3.201862e-05, 3.5351912e-05, 3.9032216e-05, 4.3095657e-05,
    4.7582122e-05, 5.253565e-05, 5.8004864e-05, 6.4043449e-05, 7.0710679e-05,
    7.8072001e-05, 8.6199671e-05, 9.5173471e-05, 0.00010508149,
    0.00011602097, 0.00012809932, 0.00014143507, 0.00015615914,
    0.00017241606, 0.0001903654, 0.00021018335, 0.00023206445, 0.00025622347,
    0.00028289757, 0.00031234856, 0.00034486555, 0.00038076771,
    0.00042040746, 0.00046417389, 0.00051249663, 0.00056584999,
    0.00062475769, 0.00068979797, 0.00076160925, 0.00084089644, 0.0009284378,
    0.0010250926, 0.0011318097, 0.0012496365, 0.0013797297, 0.0015233662,
    0.0016819559, 0.0018570556, 0.002050384, 0.0022638389, 0.0024995154,
    0.0027597269, 0.0030470278, 0.003364238, 0.0037144714, 0.0041011657,
    0.0045281167, 0.0049995154, 0.005519989, 0.0060946464, 0.0067291284,
    0.007429663, 0.0082031266, 0.0090571115, 0.01,
)
SI_FAST_MEV = (
    0.01, 0.010544119, 0.011117845, 0.011722788, 0.012360647, 0.013033213,
    0.013742375, 0.014490124, 0.015278559, 0.016109895, 0.016986465,
    0.01791073, 0.018885287, 0.019912872, 0.020996369, 0.022138822,
    0.023343437, 0.024613598, 0.02595287, 0.027365016, 0.028853998,
    0.030423999, 0.032079427, 0.033824929, 0.035665408, 0.037606031,
    0.039652247, 0.041809801, 0.044084752, 0.046483487, 0.049012742,
    0.051679618, 0.054491605, 0.057456597, 0.060582919, 0.063879351,
    0.067355148, 0.07102007, 0.074884407, 0.07895901, 0.08325532,
    0.087785401, 0.092561972, 0.097598445, 0.10290896, 0.10850843,
    0.11441258, 0.12063799, 0.12720213, 0.13412344, 0.14142136, 0.14911636,
    0.15723007, 0.16578525, 0.17480594, 0.18431747, 0.19434653, 0.2049213,
    0.21607145, 0.22782831, 0.24022489, 0.25329598, 0.2670783, 0.28161053,
    0.2969335, 0.31309022, 0.33012605, 0.34808884, 0.36702901, 0.38699976,
    0.40805715, 0.43026032, 0.4536716, 0.47835674, 0.50438504, 0.53182959,
    0.56076745, 0.59127987, 0.62345253, 0.65737577, 0.69314484, 0.73086017,
    0.77062766, 0.81255898, 0.85677186, 0.90339045, 0.95254564, 1.0043755,
    1.0590254, 1.116649, 1.177408, 1.241473, 1.309024, 1.3802504, 1.4553525,
    1.534541, 1.6180383, 1.7060788, 1.7989098, 1.8967919, 2.0,
)
SI_GROUPS = {"thermal": SI_THERMAL_MEV, "epithermal": SI_EPI_MEV, "fast": SI_FAST_MEV}


# =====================================================================
# Derived constants
# =====================================================================
def _ellipsoid_volume(a, b, c):
    return 4.0 / 3.0 * math.pi * a * b * c


V_TUMOUR = 4.0 / 3.0 * math.pi * TUMOUR_R**3
V_BRAIN_TOTAL = _ellipsoid_volume(*SEMI["brain"])
V_CRANIUM = _ellipsoid_volume(*SEMI["cranium"]) - V_BRAIN_TOTAL
V_SCALP = _ellipsoid_volume(*SEMI["scalp"]) - _ellipsoid_volume(*SEMI["cranium"])
ROI_VOLUME = {"tumour": V_TUMOUR, "brain": V_BRAIN_TOTAL - V_TUMOUR, "scalp": V_SCALP}
ROI_MASS = {"tumour": V_TUMOUR * DENSITY["brain"],
            "brain": (V_BRAIN_TOTAL - V_TUMOUR) * DENSITY["brain"],
            "scalp": V_SCALP * DENSITY["scalp"]}
ROIS = ("tumour", "brain", "scalp")


def _kerma_filter(table_mev, low_extend_ev=None):
    """EnergyFunctionFilter (log-log) from a (MeV, value) table; optionally
    extend flat to `low_extend_ev` to reproduce MCNP's endpoint rule."""
    e = [p[0] * 1.0e6 for p in table_mev]
    y = [p[1] for p in table_mev]
    if low_extend_ev is not None and low_extend_ev < e[0]:
        e = [low_extend_ev] + e
        y = [y[0]] + y
    f = openmc.EnergyFunctionFilter(e, y)
    f.interpolation = "log-log"
    return f


# =====================================================================
# Model construction
# =====================================================================
def make_material(region, b10_ppm=0.0, name=None):
    """ICRU 46 material for `region` with an optional B-10 mass loading in
    ppm (ug/g), taken out of the oxygen fraction as in W003."""
    m = openmc.Material(name=name or region)
    m.set_density("g/cm3", DENSITY[region])
    wo = b10_ppm * 1e-6
    for key, frac in COMPOSITION[region].items():
        f = frac - wo if key == "O16" else frac
        if any(ch.isdigit() for ch in key):
            m.add_nuclide(key, f, "wo")
        else:
            m.add_element(key, f, "wo")
    if wo > 0.0:
        m.add_nuclide("B10", wo, "wo")
    if region != "air":
        m.add_s_alpha_beta(SAB)
    return m


def quadric(region):
    A, B, C, zc = SQ[region]
    return openmc.Quadric(a=A, b=B, c=C, j=-2.0 * C * zc, k=C * zc * zc - 1.0)


def make_geometry(mats, tumour):
    s_br, s_cr, s_sc = quadric("brain"), quadric("cranium"), quadric("scalp")
    s_world = openmc.Sphere(r=R_WORLD, boundary_type="vacuum")
    cells = {}
    if tumour:
        s_t = openmc.Sphere(z0=TUMOUR_Z, r=TUMOUR_R)
        cells["tumour"] = openmc.Cell(name="tumour", fill=mats["tumour"], region=-s_t)
        cells["brain"] = openmc.Cell(name="brain", fill=mats["brain"], region=-s_br & +s_t)
    else:
        cells["brain"] = openmc.Cell(name="brain", fill=mats["brain"], region=-s_br)
    cells["cranium"] = openmc.Cell(name="cranium", fill=mats["cranium"], region=+s_br & -s_cr)
    cells["scalp"] = openmc.Cell(name="scalp", fill=mats["scalp"], region=+s_cr & -s_sc)
    cells["air"] = openmc.Cell(name="air", fill=mats["air"], region=+s_sc & -s_world)
    return openmc.Geometry(list(cells.values())), cells


def beam_space():
    return openmc.stats.CylindricalIndependent(
        r=openmc.stats.PowerLaw(0.0, BEAM_R, 1.0),
        phi=openmc.stats.Uniform(0.0, 2.0 * math.pi),
        z=openmc.stats.Discrete([BEAM_Z], [1.0]))


def group_energy(group):
    edges = [e * 1.0e6 for e in SI_GROUPS[group]]
    n = len(edges) - 1
    return openmc.stats.Mixture([1.0 / n] * n,
                                [openmc.stats.Uniform(edges[i], edges[i + 1]) for i in range(n)])


def make_sources(beam):
    """beam: 'epithermal' (analog composite), ('group', g), ('mono_n', E_eV),
    ('mono_p', E_eV)."""
    def src(particle, energy, strength=1.0):
        s = openmc.IndependentSource()
        s.space = beam_space()
        s.angle = openmc.stats.Monodirectional(BEAM_DIR)
        s.energy = energy
        s.particle = particle
        s.strength = strength
        return s
    if beam == "epithermal":
        return [src("neutron", group_energy(g), GROUP_FRACTIONS[g])
                for g in ("thermal", "epithermal", "fast")]
    kind, arg = beam
    if kind == "group":
        return [src("neutron", group_energy(arg))]
    if kind == "mono_n":
        return [src("neutron", openmc.stats.Discrete([arg], [1.0]))]
    if kind == "mono_p":
        return [src("photon", openmc.stats.Discrete([arg], [1.0]))]
    raise ValueError(beam)


def profile_tallies(photon_only=False):
    mesh = openmc.RegularMesh()
    mesh.dimension = (1, 1, MESH_NZ)
    mesh.lower_left = (-MESH_XY, -MESH_XY, MESH_Z0)
    mesh.upper_right = (MESH_XY, MESH_XY, MESH_Z1)
    mf = openmc.MeshFilter(mesh)
    pf_n, pf_p = openmc.ParticleFilter("neutron"), openmc.ParticleFilter("photon")
    out = []
    if not photon_only:
        kn_th = _kerma_filter(KERMA_BRAIN_N_CGYCM2, low_extend_ev=1e-5)
        kn_fa = _kerma_filter(KERMA_BRAIN_N_CGYCM2, low_extend_ev=1e-5)
        kb = _kerma_filter(KERMA_B10_PER_PPM_GYCM2, low_extend_ev=1e-6)
        t = openmc.Tally(name="prof_thermal_n")
        t.filters = [mf, pf_n, openmc.EnergyFilter([0.0, 0.5]), kn_th]
        t.scores = ["flux"]
        out.append(t)
        t = openmc.Tally(name="prof_fast_n")
        t.filters = [mf, pf_n, openmc.EnergyFilter([0.5, 2.0e7]), kn_fa]
        t.scores = ["flux"]
        out.append(t)
        t = openmc.Tally(name="prof_boron")
        t.filters = [mf, pf_n, kb]
        t.scores = ["flux"]
        out.append(t)
    kp = _kerma_filter(KERMA_BRAIN_P_CGYCM2)
    t = openmc.Tally(name="prof_gamma")
    t.filters = [mf, pf_p, kp]
    t.scores = ["flux"]
    out.append(t)
    return out


def roi_tallies(cells):
    out = []
    for roi in ROIS:
        if roi not in cells:
            continue
        cf = openmc.CellFilter([cells[roi]])
        t1 = openmc.Tally(name=f"{roi}_thermal_flux")
        t1.filters = [cf, openmc.ParticleFilter("neutron"), openmc.EnergyFilter([0.0, 0.5])]
        t1.scores = ["flux"]
        t2 = openmc.Tally(name=f"{roi}_b10_na")
        t2.filters = [cf]; t2.nuclides = ["B10"]; t2.scores = ["(n,a)"]
        t3 = openmc.Tally(name=f"{roi}_n14_np")
        t3.filters = [cf]; t3.nuclides = ["N14"]; t3.scores = ["(n,p)"]
        t4 = openmc.Tally(name=f"{roi}_fast_heat")
        t4.filters = [cf, openmc.ParticleFilter("neutron"), openmc.EnergyFilter([1.0e4, 2.1e7])]
        t4.scores = ["heating"]
        # T5 photon dose: photon flux x NIST photon kerma (brain), the same
        # estimator as the reference F44 tally. NOT the OpenMC `heating`
        # score for photons, which was found at build time to report ~1e-4
        # of the physical photon energy deposition in this build (see
        # ledger, tooling note T5). Raw tally units: cGy cm3 per source;
        # converted to eV deposited in the ROI in read_statepoint.
        t5 = openmc.Tally(name=f"{roi}_photon_fluxkerma_raw")
        t5.filters = [cf, openmc.ParticleFilter("photon"), _kerma_filter(KERMA_BRAIN_P_CGYCM2)]
        t5.scores = ["flux"]
        out.extend([t1, t2, t3, t4, t5])
    g1 = openmc.Tally(name="global_b10_na"); g1.nuclides = ["B10"]; g1.scores = ["(n,a)"]
    g2 = openmc.Tally(name="global_h1_ng"); g2.nuclides = ["H1"]; g2.scores = ["(n,gamma)"]
    out.extend([g1, g2])
    return out


def build(rundir, beam, ppm, n_particles, seed, tumour, profile, roi, photon_only=False):
    """ppm: dict region -> B-10 ppm for brain, cranium, scalp, tumour."""
    os.makedirs(rundir, exist_ok=True)
    mats = {r: make_material(r, ppm.get(r, 0.0)) for r in ("brain", "cranium", "scalp", "air")}
    if tumour:
        mats["tumour"] = make_material("brain", ppm.get("tumour", 0.0), name="tumour")
    materials = openmc.Materials(mats.values())
    materials.cross_sections = DATA
    geometry, cells = make_geometry(mats, tumour)
    st = openmc.Settings()
    st.run_mode = "fixed source"
    st.photon_transport = True
    st.source = make_sources(beam)
    st.batches = n_batches(n_particles)
    st.particles = n_particles // st.batches
    st.seed = seed
    st.output = {"tallies": False, "summary": True}
    tallies = openmc.Tallies()
    if profile:
        tallies.extend(profile_tallies(photon_only=photon_only))
    if roi:
        tallies.extend(roi_tallies(cells))
    openmc.Model(geometry, materials, st, tallies).export_to_model_xml(
        os.path.join(rundir, "model.xml"))
    dens = {}
    for r, m in mats.items():
        dens[r] = float(m.get_nuclide_atom_densities().get("B10", 0.0))
    return dens


EV_PER_CGY_PER_G = 1.0e-2 / 1.602176634e-19 * 1.0e-3   # 1 cGy = 6.2415e13 eV per gram
ROI_DENSITY = {"tumour": DENSITY["brain"], "brain": DENSITY["brain"], "scalp": DENSITY["scalp"]}


def read_statepoint(rundir, profile, roi, batches=BATCHES):
    sp = openmc.StatePoint(os.path.join(rundir, f"statepoint.{batches}.h5"))
    out = {}
    for t in sp.tallies.values():
        mean = t.mean.ravel(); std = t.std_dev.ravel()
        if t.name.startswith("prof_"):
            out[t.name] = {"mean": [float(x) for x in mean], "std": [float(x) for x in std]}
        else:
            out[t.name] = {"mean": float(mean[0]), "std": float(std[0])}
        if t.name.endswith("_photon_fluxkerma_raw"):
            # cGy cm3 per source -> eV deposited in the ROI per source:
            # (raw / V) [cGy] * EV_PER_CGY_PER_G * mass = raw * rho * EV_PER_CGY_PER_G
            r = t.name.split("_photon_fluxkerma_raw")[0]
            fac = ROI_DENSITY[r] * EV_PER_CGY_PER_G
            out[f"{r}_photon_heat"] = {"mean": float(mean[0]) * fac, "std": float(std[0]) * fac,
                                       "derived_from": t.name, "factor_eV_per_raw": fac}
    sp.close()
    return out


def run_case(tag, beam, ppm, n, seed, tumour=False, profile=False, roi=True,
             photon_only=False, root=None, params=None):
    root = root or RUNROOT
    rundir = os.path.join(root, tag)
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    n = scaled(n)
    dens = build(rundir, beam, ppm, n, seed, tumour, profile, roi, photon_only)
    t0 = time.time()
    openmc.run(cwd=rundir, output=False, threads=THREADS, openmc_exec=OPENMC_EXEC)
    wall = time.time() - t0
    tallies = read_statepoint(rundir, profile, roi, n_batches(n))
    beam_desc = beam if isinstance(beam, str) else list(beam)
    out = {"beam": beam_desc, "ppm": ppm, "params": params, "seed": seed,
           "histories": n, "batches": n_batches(n), "tumour": tumour, "wall_s": wall,
           "b10_atom_density": dens, "tallies": tallies}
    if roi:
        v_brain = V_BRAIN_TOTAL - V_TUMOUR if tumour else V_BRAIN_TOTAL
        n_tot = (dens["brain"] * v_brain + dens["cranium"] * V_CRANIUM
                 + dens["scalp"] * V_SCALP)
        if tumour:
            n_tot += dens["tumour"] * V_TUMOUR
        out["b10_atoms_total_bcm3"] = n_tot
    print(f"  [{tag}] {n:,} histories, {wall:.0f} s", flush=True)
    return out


# =====================================================================
# Post-processing (harness side; the checker re-derives independently)
# =====================================================================
def kerma_profile(tallies):
    """Convert profile tallies to Gy/min (boron: Gy/min per ppm) per bin,
    ordered by DEPTHS. Returns {comp: {'value': [...], 'std': [...]}}."""
    out = {}
    fac_n = SOURCE_PER_MIN / BIN_VOLUME / 100.0    # cGy cm2 -> Gy/min
    fac_b = SOURCE_PER_MIN / BIN_VOLUME             # Gy cm2 per ppm -> Gy/min/ppm
    for name, comp, fac in (("prof_thermal_n", "thermal", fac_n), ("prof_fast_n", "fast", fac_n),
                            ("prof_boron", "boron", fac_b), ("prof_gamma", "gamma", fac_n)):
        if name in tallies:
            out[comp] = {"value": [v * fac for v in tallies[name]["mean"]],
                         "std": [s * fac for s in tallies[name]["std"]]}
    return out


def combine_groups(runs):
    """Linear combination of group-split profile runs with the physical
    fractions; sigma in quadrature. runs: {group: run dict}."""
    comps = None
    for g, r in runs.items():
        prof = kerma_profile(r["tallies"])
        if comps is None:
            comps = {c: {"value": np.zeros(MESH_NZ), "var": np.zeros(MESH_NZ)} for c in prof}
        for c in prof:
            f = GROUP_FRACTIONS[g]
            comps[c]["value"] += f * np.array(prof[c]["value"])
            comps[c]["var"] += (f * np.array(prof[c]["std"]))**2
    return {c: {"value": [float(x) for x in d["value"]],
                "std": [float(math.sqrt(v)) for v in d["var"]]} for c, d in comps.items()}


def get(res, name):
    t = res["tallies"][name]
    return t["mean"], t["std"]


def dose(res, roi):
    """Direct weighted dose (MeV per g per source particle) with 1-sigma."""
    wB, wN, wf, wg = W[roi]
    m = ROI_MASS[roi]
    b, sb = get(res, f"{roi}_b10_na")
    nn, sn = get(res, f"{roi}_n14_np")
    fh, sf = get(res, f"{roi}_fast_heat")
    ph, sp_ = get(res, f"{roi}_photon_heat")
    val = (wB * E_B * b + wN * E_N * nn + wf * fh / 1e6 + wg * ph / 1e6) / m
    sig = math.sqrt((wB * E_B * sb)**2 + (wN * E_N * sn)**2 + (wf * sf / 1e6)**2
                    + (wg * sp_ / 1e6)**2) / m
    return val, sig


def brackets(lo, hi):
    """Certified [L,U] per ROI from the two corner runs only (W003 §4)."""
    out = {}
    nt_lo, nt_hi = lo["b10_atoms_total_bcm3"], hi["b10_atoms_total_bcm3"]
    G1lo, sG1lo = get(lo, "global_b10_na"); G1hi, sG1hi = get(hi, "global_b10_na")
    G2lo, sG2lo = get(lo, "global_h1_ng"); G2hi, sG2hi = get(hi, "global_h1_ng")
    KG_lo, KG_hi = G1lo / nt_lo, G1hi / nt_hi
    G1_lower, G1_upper = nt_lo * KG_hi, nt_hi * KG_lo
    for roi in ROIS:
        wB, wN, wf, wg = W[roi]; m = ROI_MASS[roi]
        nlo = lo["b10_atom_density"][roi]; nhi = hi["b10_atom_density"][roi]
        b_lo, sb_lo = get(lo, f"{roi}_b10_na"); b_hi, sb_hi = get(hi, f"{roi}_b10_na")
        K_lo, K_hi = b_lo / nlo, b_hi / nhi
        B_low, B_up = nlo * K_hi, nhi * K_lo
        sB_low, sB_up = nlo * sb_hi / nhi, nhi * sb_lo / nlo
        n_lo_, sn_lo = get(lo, f"{roi}_n14_np"); n_hi_, sn_hi = get(hi, f"{roi}_n14_np")
        f_lo, sf_lo = get(lo, f"{roi}_fast_heat"); f_hi, sf_hi = get(hi, f"{roi}_fast_heat")
        fmean = 0.5 * (f_lo + f_hi)
        fast_ok = abs(f_hi - f_lo) <= 3 * math.hypot(sf_lo, sf_hi) + EPS_F * fmean
        p_lo, sp_lo = get(lo, f"{roi}_photon_heat"); p_hi, sp_hi = get(hi, f"{roi}_photon_heat")
        bg = BORON_GAMMA_BRANCH
        det = G2lo * bg * G1hi - G2hi * bg * G1lo
        gH = (p_lo * bg * G1hi - p_hi * bg * G1lo) / det
        gB = (G2lo * p_hi - G2hi * p_lo) / det
        photon_valid = (gH > -3 * sp_lo / max(G2lo, 1e-30)) and \
                       (gB >= -3 * (sp_hi + sp_lo) / max(bg * G1lo, 1e-30))
        gH_c, gB_c = max(gH, 0.0), max(gB, 0.0)
        ph_low = gH_c * G2hi + gB_c * bg * G1_lower
        ph_up = gH_c * G2lo + gB_c * bg * G1_upper
        sph = math.hypot(sp_lo, sp_hi)
        L = (wB * E_B * B_low + wN * E_N * n_hi_ + wf * (fmean / 1e6) + wg * ph_low / 1e6) / m
        U = (wB * E_B * B_up + wN * E_N * n_lo_ + wf * (fmean / 1e6) + wg * ph_up / 1e6) / m
        sfast = (3 * math.hypot(sf_lo, sf_hi) + EPS_F * fmean) / 3 / 1e6
        sL = math.sqrt((wB * E_B * sB_low)**2 + (wN * E_N * sn_hi)**2 + (wf * sfast)**2
                       + (wg * sph / 1e6)**2) / m
        sU = math.sqrt((wB * E_B * sB_up)**2 + (wN * E_N * sn_lo)**2 + (wf * sfast)**2
                       + (wg * sph / 1e6)**2) / m
        out[roi] = {"L": L, "U": U, "sL": sL, "sU": sU,
                    "kernels": {"K_lo": K_lo, "K_hi": K_hi, "gammaH": gH, "gammaB": gB},
                    "photon_valid": bool(photon_valid), "fast_applicability_ok": bool(fast_ok)}
    return out


def ppm_map(rt, rs, B):
    return {"tumour": rt * B, "brain": BRAIN_RATIO * B, "scalp": rs * B, "cranium": 0.0}


def schedule(ds):
    """Corners, six mixed corners, nominal, two interior points (W003 rule:
    fractional positions (0.25,0.75,0.25) and (0.75,0.25,0.75)) and the two
    preset-comparator blood endpoints at nominal ratios."""
    rt, rs, B = ds["rt"], ds["rs"], ds["B"]
    lo, hi = (rt[0], rs[0], B[0]), (rt[2], rs[2], B[2])
    mixed = [(rt[i], rs[j], B[k]) for i in (0, 2) for j in (0, 2) for k in (0, 2)
             if (i, j, k) not in ((0, 0, 0), (2, 2, 2))]
    def frac(axis, f):
        return axis[0] + f * (axis[2] - axis[0])
    interior = mixed + [(rt[1], rs[1], B[1]),
                        (frac(rt, .25), frac(rs, .75), frac(B, .25)),
                        (frac(rt, .75), frac(rs, .25), frac(B, .75))]
    comparator = [(rt[1], rs[1], B[0]), (rt[1], rs[1], B[2])]
    return lo, hi, interior, comparator


def dump(stage, payload):
    payload["harness_sha256"] = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
    payload["protocol_sha256"] = PROTOCOL_SHA256
    payload["openmc_version"] = openmc.__version__
    payload["threads"] = THREADS
    payload["timestamp_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    payload["history_scale"] = HISTORY_SCALE
    payload["record_location"] = "controls/" if OUTDIR == HERE else f"TEST:{OUTDIR}"
    with open(stage_file(stage), "w") as f:
        json.dump(payload, f, indent=1, sort_keys=True)
    print(f"stage {stage} -> {stage_file(stage)}", flush=True)


# =====================================================================
# Stages
# =====================================================================
def stage_c1():
    """C1 setup control on both declared sets' corners: PC-A determinism,
    PC-B direction (thermal flux and per-atom boron kernel lower at c+ than
    c-, > 3 sigma), PC-C zero-boron dominance."""
    os.makedirs(RUNROOT, exist_ok=True)
    runs, checks = {}, {}
    lo1 = schedule(DS["DS-1"])[0]
    a1 = run_case("c1_det_1", "epithermal", ppm_map(*lo1), N_PC, 42, tumour=True)
    a2 = run_case("c1_det_2", "epithermal", ppm_map(*lo1), N_PC, 42, tumour=True)
    runs["c1_det_1"], runs["c1_det_2"] = a1, a2
    checks["PC-A_deterministic"] = all(
        abs(a1["tallies"][k]["mean"] - a2["tallies"][k]["mean"]) <= 1e-12 * max(1.0, abs(a1["tallies"][k]["mean"]))
        for k in a1["tallies"])
    zero = run_case("c1_zero_boron", "epithermal", ppm_map(0.0, 0.0, 0.0), N_PC, 103, tumour=True)
    runs["c1_zero_boron"] = zero
    seed = 110
    for name, ds in DS.items():
        lo, hi, _, _ = schedule(ds)
        rlo = run_case(f"c1_{name}_corner_lo", "epithermal", ppm_map(*lo), N_PC, seed, tumour=True)
        rhi = run_case(f"c1_{name}_corner_hi", "epithermal", ppm_map(*hi), N_PC, seed + 1, tumour=True)
        seed += 2
        runs[f"c1_{name}_corner_lo"], runs[f"c1_{name}_corner_hi"] = rlo, rhi
        for roi in ROIS:
            flo, slo = get(rlo, f"{roi}_thermal_flux"); fhi, shi = get(rhi, f"{roi}_thermal_flux")
            fz, sz = get(zero, f"{roi}_thermal_flux")
            checks[f"PC-B_flux_dir_{name}_{roi}"] = (flo - fhi) > 3 * math.hypot(slo, shi)
            klo = get(rlo, f"{roi}_b10_na")[0] / rlo["b10_atom_density"][roi]
            khi = get(rhi, f"{roi}_b10_na")[0] / rhi["b10_atom_density"][roi]
            sklo = get(rlo, f"{roi}_b10_na")[1] / rlo["b10_atom_density"][roi]
            skhi = get(rhi, f"{roi}_b10_na")[1] / rhi["b10_atom_density"][roi]
            checks[f"PC-B_kernel_dir_{name}_{roi}"] = (klo - khi) > 3 * math.hypot(sklo, skhi)
            checks[f"PC-C_zeroboron_{name}_{roi}"] = (fz - flo) > 3 * math.hypot(sz, slo)
    checks = {k: bool(v) for k, v in checks.items()}
    checks["C1_all"] = all(checks.values())
    print(json.dumps(checks, indent=1, sort_keys=True))
    dump("c1", {"checks": checks, "runs": runs})
    return checks["C1_all"]


def stage_c2():
    """C2 monoenergetic beams: profiles on the unmodified phantom."""
    os.makedirs(RUNROOT, exist_ok=True)
    runs, profiles = {}, {}
    seed = 200
    for tag, E in MONO_NEUTRON_EV.items():
        r = run_case(f"c2_{tag}", ("mono_n", E), {}, N_C2_N[tag], seed, profile=True, roi=False)
        runs[f"c2_{tag}"] = r
        profiles[tag] = kerma_profile(r["tallies"])
        seed += 1
    for tag, E in MONO_PHOTON_EV.items():
        r = run_case(f"c2_{tag}", ("mono_p", E), {}, N_C2_P, seed, profile=True, roi=False,
                     photon_only=True)
        runs[f"c2_{tag}"] = r
        profiles[tag] = kerma_profile(r["tallies"])
        seed += 1
    dump("c2", {"runs": runs, "profiles_Gy_per_min": profiles, "depths_cm": DEPTHS})


def stage_configA():
    """Configuration A: generic epithermal beam, unmodified phantom,
    group-split runs combined with the physical group fractions."""
    os.makedirs(RUNROOT, exist_ok=True)
    runs = {}
    seed = 300
    for g in ("thermal", "epithermal", "fast"):
        runs[g] = run_case(f"configA_{g}", ("group", g), {}, N_A[g], seed, profile=True, roi=False)
        seed += 1
    combined = combine_groups(runs)
    dump("configA", {"runs": {f"configA_{g}": r for g, r in runs.items()},
                     "group_fractions": GROUP_FRACTIONS,
                     "profile_Gy_per_min": combined, "depths_cm": DEPTHS,
                     "per_group_profiles_Gy_per_min": {g: kerma_profile(r["tallies"]) for g, r in runs.items()}})


def stage_c3():
    """C3 flux-depression control: 30 ug/g 10B in scalp and brain (not
    cranium), same beam and tally grid as Configuration A."""
    os.makedirs(RUNROOT, exist_ok=True)
    ppm = {"brain": C3_PPM, "scalp": C3_PPM, "cranium": 0.0}
    runs = {}
    seed = 400
    for g in ("thermal", "epithermal", "fast"):
        runs[g] = run_case(f"c3_{g}", ("group", g), ppm, N_C3[g], seed, profile=True, roi=False)
        seed += 1
    combined = combine_groups(runs)
    payload = {"runs": {f"c3_{g}": r for g, r in runs.items()}, "ppm": ppm,
               "profile_Gy_per_min": combined, "depths_cm": DEPTHS}
    if os.path.exists(stage_file("configA")):
        A = json.load(open(stage_file("configA")))["profile_Gy_per_min"]
        payload["thermal_dose_depression"] = depression(A, combined)
    dump("c3", payload)


def depression(no_boron, with_boron):
    """Fig. 17 definition: (D_noB - D_B) / D_noB per depth, thermal
    neutron kerma; sigma from independent runs."""
    dep = []
    for k in range(MESH_NZ):
        a, sa = no_boron["thermal"]["value"][k], no_boron["thermal"]["std"][k]
        b, sb = with_boron["thermal"]["value"][k], with_boron["thermal"]["std"][k]
        if a > 0 and b > 0:
            d = 1.0 - b / a
            sd = (b / a) * math.hypot(sb / b, sa / a)
        else:
            d, sd = float("nan"), float("nan")
        dep.append({"depth": DEPTHS[k], "depression_pct": 100 * d, "std_pct": 100 * sd})
    return dep


def stage_scored(force=False):
    """Configuration B scored block + preset comparator, both declared sets.
    Refuses to run until the OPEN declarations (normal-tissue ratios; plan
    normalisation for W6-G4) are recorded as adopted in the ledger and the
    ADOPTED flag below is set — protocol-before-evidence."""
    if not (DECLARATIONS_ADOPTED or force or os.environ.get("W006_FORCE_SCORED") == "1"):
        sys.exit("scored block blocked: OPEN declarations not adopted (see ledger.md)")
    if os.environ.get("W006_FORCE_SCORED") == "1" and OUTDIR == HERE:
        sys.exit("W006_FORCE_SCORED is a test-only override; refuse to write into controls/")
    c1 = json.load(open(stage_file("c1")))
    assert c1["checks"]["C1_all"], "C1 setup control has not passed"
    os.makedirs(RUNROOT, exist_ok=True)
    result = {"sets": {}}
    seed = 500
    for name, ds in DS.items():
        lo, hi, interior, comparator = schedule(ds)
        runs = {}
        runs["corner_lo"] = run_case(f"B_{name}_corner_lo", "epithermal", ppm_map(*lo), N_FULL,
                                     seed, tumour=True, params=lo)
        runs["corner_hi"] = run_case(f"B_{name}_corner_hi", "epithermal", ppm_map(*hi), N_FULL,
                                     seed + 1, tumour=True, params=hi)
        seed += 2
        for i, p in enumerate(interior):
            runs[f"interior_{i}"] = run_case(f"B_{name}_interior_{i}", "epithermal", ppm_map(*p),
                                             N_INT, seed, tumour=True, params=p)
            seed += 1
        for i, p in enumerate(comparator):
            runs[f"comparator_{i}"] = run_case(f"B_{name}_comparator_{i}", "epithermal",
                                               ppm_map(*p), N_INT, seed, tumour=True, params=p)
            seed += 1
        # AD-W6-1: any scored ROI tally > 2.5 % relative error -> double once;
        # both runs logged, neither discarded.
        for key in list(runs):
            r = runs[key]
            worst = max(t["std"] / abs(t["mean"]) for n_, t in r["tallies"].items()
                        if abs(t["mean"]) > 0 and not n_.startswith("global"))
            if worst > AD_W6_1_RELERR:
                p = tuple(r["params"])
                n0 = N_FULL if key.startswith("corner_") else N_INT
                runs[f"{key}_AD1"] = run_case(f"B_{name}_{key}_AD1", "epithermal", ppm_map(*p),
                                              2 * n0, seed, tumour=True, params=p)
                seed += 1
        params = {"corner_lo": lo, "corner_hi": hi,
                  "interior": interior, "comparator": comparator}
        br = brackets(runs["corner_lo"], runs["corner_hi"])
        direct = {k: {roi: dict(zip(("value", "std"), dose(r, roi))) for roi in ROIS}
                  for k, r in runs.items()}
        result["sets"][name] = {"declared": ds, "params": params, "runs": runs,
                                "brackets": br, "direct_doses": direct}
    result["brain_ratio"] = BRAIN_RATIO
    result["weights"] = W
    result["roi_mass_g"] = ROI_MASS
    result["plan_normalisation"] = PLAN_NORMALISATION
    dump("scored", result)


# ---------------------------------------------------------------------
# Declarations for Configuration B / W6-G4 that the frozen protocol left
# unstated, ADOPTED by the principal 2026-08-17 (ledger entry "W006 open
# declarations adopted") before any scored run:
#   1. normal-tissue ratios: continuity reading (BRAIN_RATIO = 1.0 fixed;
#      scalp takes the W003 rho_M band, DS["*"]["rs"]) — as coded above;
#   2. W6-G4 reported as the normalisation-free divergence window and the
#      certified/comparator ratio; no plan is declared, so the binary
#      DIVERGENCE/NULL is recorded as undeclared-by-adoption
#      (PLAN_NORMALISATION stays None; a later cited plan may be added by a
#      further ledger entry and the checker will then decide the binary);
#   3. ROI reading: whole-ROI tallies with min/max over the declared set.
# ---------------------------------------------------------------------
DECLARATIONS_ADOPTED = True
PLAN_NORMALISATION = None


def stage_assemble():
    out = {"protocol_sha256": PROTOCOL_SHA256, "archive_sha256": ARCHIVE_SHA256,
           "harness_sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
           "constants": {"W": W, "E_B": E_B, "E_N": E_N, "EPS_F": EPS_F, "DS": DS,
                         "BRAIN_RATIO": BRAIN_RATIO, "ROI_MASS": ROI_MASS,
                         "ROI_VOLUME": ROI_VOLUME, "BIN_VOLUME": BIN_VOLUME,
                         "SOURCE_PER_MIN": SOURCE_PER_MIN, "DEPTHS": DEPTHS,
                         "GROUP_FRACTIONS": GROUP_FRACTIONS, "CRITERIA": CRITERIA,
                         "G1_TOL": G1_TOL, "G1C_BAND_PCT": G1C_BAND_PCT,
                         "G1C_PEAK_DEPTHS": G1C_PEAK_DEPTHS, "GAP_KILL": GAP_KILL,
                         "C3_PPM": C3_PPM, "TUMOUR": [TUMOUR_R, TUMOUR_Z]},
           "stages": {}}
    for st in ("c1", "c2", "configA", "c3", "scored"):
        f = stage_file(st)
        if os.path.exists(f):
            out["stages"][st] = json.load(open(f))
    if "configA" in out["stages"] and "c3" in out["stages"]:
        out["stages"]["c3"]["thermal_dose_depression"] = depression(
            out["stages"]["configA"]["profile_Gy_per_min"], out["stages"]["c3"]["profile_Gy_per_min"])
    if os.path.exists(REFERENCE):
        out["reference"] = json.load(open(REFERENCE))
    with open(RESULTS, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("assembled ->", RESULTS, "stages:", sorted(out["stages"]))


def stage_bench():
    """Timing only. Writes to the scratchpad; produces no record."""
    os.makedirs(SCRATCH, exist_ok=True)
    n = int(os.environ.get("W006_BENCH_N", "400000"))
    for tag, beam, ppm, kw in (
            ("bench_A_epi", ("group", "epithermal"), {}, dict(profile=True, roi=False)),
            ("bench_A_fast", ("group", "fast"), {}, dict(profile=True, roi=False)),
            ("bench_A_th", ("group", "thermal"), {}, dict(profile=True, roi=False)),
            ("bench_B", "epithermal", ppm_map(3.5, 1.0, 25.0), dict(tumour=True, roi=True)),
            ("bench_p1", ("mono_p", 1.0e6), {}, dict(profile=True, roi=False, photon_only=True)),
            ("bench_n253", ("mono_n", 0.0253), {}, dict(profile=True, roi=False))):
        r = run_case(tag, beam, ppm, n, 7, root=SCRATCH, **kw)
        rate = n / r["wall_s"]
        print(f"{tag}: {rate:,.0f} hist/s")
        if kw.get("profile"):
            prof = kerma_profile(r["tallies"])
            for c in prof:
                v = prof[c]["value"]; s = prof[c]["std"]
                rel = [s[k] / v[k] if v[k] else float("nan") for k in range(MESH_NZ)]
                print(f"   {c}: rel err at depth 0.2 {rel[DEPTHS.index(0.2)]:.3f}, 5.0 {rel[DEPTHS.index(5.0)]:.3f}, "
                      f"13.4 {rel[DEPTHS.index(13.4)]:.3f}, 17.4 {rel[DEPTHS.index(17.4)]:.3f}; "
                      f"value at 3.0: {v[DEPTHS.index(3.0)]:.4g}")
        else:
            for roi in ROIS:
                d, s = dose(r, roi)
                worst = max(r["tallies"][f"{roi}_{k}"]["std"] / max(r["tallies"][f"{roi}_{k}"]["mean"], 1e-300)
                            for k in ("thermal_flux", "b10_na", "n14_np", "fast_heat", "photon_heat"))
                print(f"   {roi}: dose {d:.4g} MeV/g/src (+-{s/d:.3%}), worst tally rel err {worst:.3%}")


if __name__ == "__main__":
    os.environ.setdefault("OPENMC_CROSS_SECTIONS", DATA)
    stages = {"bench": stage_bench, "c1": stage_c1, "c2": stage_c2, "c3": stage_c3,
              "configA": stage_configA, "scored": stage_scored, "assemble": stage_assemble}
    if len(sys.argv) < 2 or sys.argv[1] not in stages:
        sys.exit(__doc__)
    stages[sys.argv[1]]()
