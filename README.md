# The Weak Memory of the Geodynamo: An RMT Verdict on Geomagnetic Reversals

Ruqing Chen, GUT Geoservice Inc., Montreal, Canada · June 2026

---

The **tenth** study in a unified Random Matrix Theory (RMT) program, and the first to
place a **deep-Earth, core-generated process** on the program's memory spectrum.

**The question:** does the geodynamo have memory? A reversal could be a charge-and-release
event (lose dipole stability → reverse → rebuild, imposing a recovery interval → level
repulsion), or it could be triggered randomly by core turbulence (→ Poisson). The
paleomagnetic literature, fitting gamma distributions to chron durations, has hovered
ambiguously near the memoryless value (k≈1.2–1.6 vs Poisson k=1) for decades.

**The answer:** the geodynamo has genuine but **weak memory**.

## Key result

Analyzing the **201 reversals** of the Cande & Kent (1995) C-sequence (0–83 Ma), with
local unfolding to remove the reversal-frequency trend:

| | ⟨r⟩ | CV | gamma k | position |
|---|---|---|---|---|
| Raw (no unfolding) | 0.456 | 1.49 | 0.91 | Poisson–GOE |
| **Unfolded (trend removed)** | **0.457** | **0.78** | **0.91** | **weak memory** |
| Poisson null (matched) | 0.389 | — | 1.0 | reference |
| GOE null (matched) | 0.531 | — | — | reference |

Reference values: Poisson 0.386, GOE 0.531, GUE 0.603.

- **+3.0σ above Poisson** (only 0.2% of Poisson realizations reach 0.457) → not memoryless
- **−3.3σ below GOE** → far from strong charge-and-release repulsion
- Robust to unfolding (CV drops 1.49→0.78, ⟨r⟩ unchanged)
- gamma k=0.91–1.14 reproduces the literature → method validated

### The cryptochron dependence (the key insight)

⟨r⟩ rises monotonically toward GOE as the shortest polarity intervals are removed:

| shortest chrons removed | ⟨r⟩ |
|---|---|
| 0% (all) | 0.457 |
| 10% | 0.487 |
| 20% | 0.530 (at GOE) |
| 25% | 0.534 |

The very brief intervals (cryptochrons/microchrons, some lasting only a few kyr) pull the
statistic toward Poisson. **The inferred memory strength of the geodynamo is a direct
function of how many brief intervals are accepted as true global reversals** — connecting
to an active paleomagnetic controversy (Laj & Channell 2007; the Mono Lake–Laschamp debate).
This turns a qualitative debate into a falsifiable measurement: as the global reality of the
shortest events is settled, ⟨r⟩ is pinned to a definite point on 0.46–0.53.

## Physical interpretation

The geodynamo falls at the **weak-memory end** of the program's spectrum:

- **Strong memory (GOE):** crustal charge-and-release — sealed ore veins, volcanic/geyser
  reservoirs, mantle plumes.
- **Weak memory (⟨r⟩≈0.46):** the geodynamo — a turbulent, chaotic MHD system where any
  reservoir memory is heavily diluted by core turbulence.
- **Memoryless (Poisson):** single solar active region; memoryless impact bombardment.

Core memory is real but diluted by turbulence — exactly as expected for a chaotic dynamo.

## Honest notes on scope

- **C-sequence only.** We use the well-calibrated CK95 0–83 Ma record. We do **not** include
  the M-sequence: a naive evenly-spaced rendering of M-sequence ages produced a spurious high
  ⟨r⟩ (a data artifact), which we decline to report.
- **Superchron** (Cretaceous Normal Superchron, 83–121 Ma) treated as a separate known anomaly.
- Events are **discrete and independently dated** (marine anomalies + Ar/Ar + astronomical
  tuning) — no peak-extraction from a continuous curve, so no extraction artifact applies.

## Reproduce

```bash
pip install -r requirements.txt
python code/geomagnetic_reversal_rmt_pipeline.py data/gpts_reversals_cseq.csv
```

## Repository Structure

```
geomagnetic-reversal-rmt/
├── README.md
├── LICENSE                       # MIT
├── requirements.txt
├── .zenodo.json
├── CITATION.cff
├── paper/
│   ├── paper.tex
│   ├── paper.pdf
│   └── figs/fig1_geodynamo.pdf
├── code/
│   └── geomagnetic_reversal_rmt_pipeline.py
├── data/
│   └── gpts_reversals_cseq.csv   # 201 reversal ages, CK95 C-sequence
├── figures/
└── results/
    └── geomag_rmt_results.json
```

## The unified program (10 systems)

| # | Domain | Memory | DOI |
|---|---|---|---|
| 1 | Stratigraphy | GOE | [20774581](https://zenodo.org/records/20774581) |
| 2 | Seismotectonics | scale-dependent | [20768130](https://zenodo.org/records/20768130) |
| 3 | Mantle plumes | GOE | [20768420](https://zenodo.org/records/20768420) |
| 4 | Metallogeny | GOE/GUE | [20768849](https://zenodo.org/records/20768849) |
| 5 | Evolution | GOE | [20783763](https://zenodo.org/records/20783763) |
| 6 | Hydrogeology | GOE/super-GUE | [20780389](https://zenodo.org/records/20780389) |
| 7 | Solar flares | Poisson (depletion limit) | [20784967](https://zenodo.org/records/20784967) |
| 8 | Orbital architecture | GUE (local-frame) | [20785613](https://zenodo.org/records/20785613) |
| 9 | Lunar craters | CSR / clustering (negative control) | [20787412](https://zenodo.org/records/20787412) |
| **10** | **Geomagnetic reversals** | **weak memory (Poisson–GOE)** | this work |

## License

MIT (code) · Reversal ages from Cande & Kent (1995).
