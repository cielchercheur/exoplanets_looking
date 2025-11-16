# Sorting and Ranking TOI List

`sorting_toi.py` reads an ExoFOP-style TESS Objects of Interest (TOI) CSV, 
filters out low‑quality or non‑M‑dwarf candidates, 
and assigns each remaining target a priority score for follow‑up transit observations. 

The result is a ranked list of promising M‑dwarf planet candidates.

Run:
```
python sorting_toi.py --toi_csv tois.csv
```

## Features

- Applies physical and observational cuts:
  - Stellar effective temperature, surface gravity, radius.
  - TESS magnitude, transit depth, period, and duration.
  - TESS disposition (`PC` only) and comment‑based vetting (e.g. binaries, false positives).
- Computes:
  - **Stellar merit score** (M‑dwarf friendly host).
  - **Observational score** (brightness, depth, period, duration).
  - **Combined priority score** used for ranking.
- Writes the ranked subset to CSV and prints a compact summary.

## Parameters
```python
TEFF_MIN, TEFF_MAX      = 2400.0, 4000.0        # K

LOGG_MIN                = 4.3

RSTAR_MIN, RSTAR_MAX    = 0.1, 0.7              # R_sun

TESS_DISP_OK            = {"PC"}

TESS_MAG_CUTOFF         = 14.0

DEPTH_PPM_MIN, DEPTH_PPM_MAX = 2000.0, 7000.0

PERIOD_PIVOT_D, PERIOD_MAX_D = 5.0, 15.0

DUR_MIN_H, DUR_PEAK_H, DUR_MAX_H = 0.5, 2.0, 5.0

# Brightness reference range used for normalization
TBRIGHT_MIN_REF, TBRIGHT_MAX_REF = 12.0, 15.0

COMMENT_FLAGS = [
    'v-shaped','v shaped','eb','eclips','odd-even','sb2','binary','fp','false',
    'retired','low snr','contamin','centroid offset'
]
