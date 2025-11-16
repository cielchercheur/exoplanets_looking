# Sorting and ranking toi list
# Run:
# python sorting_toi.py --toi_csv tois.csv

import re
import numpy as np
import pandas as pd
import argparse

# Defaults
TEFF_MIN, TEFF_MAX       = 2400.0, 4000.0   # K
LOGG_MIN                 = 4.3
RSTAR_MIN, RSTAR_MAX     = 0.1, 0.7         # R_sun

TESS_DISP_OK             = {"PC"}
TESS_MAG_CUTOFF          = 14.0

# Brightness reference range used for normalization
TBRIGHT_MIN_REF, TBRIGHT_MAX_REF = 12.0, 15.0

DEPTH_PPM_MIN, DEPTH_PPM_MAX = 2000.0, 7000.0
PERIOD_PIVOT_D, PERIOD_MAX_D = 5.0, 15.0

DUR_MIN_H, DUR_PEAK_H, DUR_MAX_H = 0.5, 2.0, 5.0

COMMENT_FLAGS = [
    'v-shaped','v shaped','eb','eclips','odd-even','sb2','binary',
    'fp','false positive','retired','low snr','contamin','centroid offset'
]

#COMMENT_FLAGS = []

# Regex: case-insensitive; handle v-shaped/v shaped; add word boundaries for key phrases
FLAGS_PATTERN = re.compile(
    r'(?:v[ -]?shaped|\beb\b|eclips|odd-even|\bsb2\b|\bbinary\b|fp|false positive|\blow snr\b|contamin|centroid offset)',
    flags=re.IGNORECASE
)

NUM_COLS = [
    'Stellar Eff Temp (K)', 'Stellar log(g) (cm/s^2)', 'Stellar Radius (R_Sun)',
    'TESS Mag', 'Depth (ppm)', 'Period (days)', 'Duration (hours)'
]

def filter_candidates(df: pd.DataFrame) -> pd.DataFrame:
    # Coerce numeric columns
    df = df.copy()
    for col in NUM_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Vectorized comment flag
    comments = df.get('Comments', pd.Series(index=df.index, dtype=object))
    flagged = comments.fillna('').str.contains(FLAGS_PATTERN, na=False)

    gate = (
        df['TESS Disposition'].isin(TESS_DISP_OK) &
        df['Stellar Eff Temp (K)'].between(TEFF_MIN, TEFF_MAX, inclusive='both') &
        (df['Stellar log(g) (cm/s^2)'] >= LOGG_MIN) &
        (df['Stellar Radius (R_Sun)'] <= RSTAR_MAX) &
        (df['TESS Mag'] <= TESS_MAG_CUTOFF) &
        (df['Depth (ppm)'] >= DEPTH_PPM_MIN) &
        (df['Period (days)'] <= PERIOD_MAX_D) &
        df['Duration (hours)'].between(DUR_MIN_H, DUR_MAX_H, inclusive='both') &
        (~flagged)
    )
    return df[gate].copy()

def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    # Convenience aliases
    teff   = df['Stellar Eff Temp (K)']
    rstar  = df['Stellar Radius (R_Sun)']
    tmag   = df['TESS Mag']
    depth  = df['Depth (ppm)']
    period = df['Period (days)']
    dur    = df['Duration (hours)']

    # Stellar merit (cooler, smaller better)
    m_teff = (TEFF_MAX - teff) / max(1e-9, (TEFF_MAX - TEFF_MIN))
    m_rad  = (RSTAR_MAX - rstar) / max(1e-9, (RSTAR_MAX - RSTAR_MIN))
    m_teff = m_teff.clip(0, 1).fillna(0.0)
    m_rad  = m_rad.clip(0, 1).fillna(0.0)
    score_m = 0.6 * m_teff + 0.4 * m_rad

    # Observability (brighter, deeper, shorter period, ~2 h duration)
    tbright = 1.0 - (tmag - TBRIGHT_MIN_REF) / max(1e-9, (TBRIGHT_MAX_REF - TBRIGHT_MIN_REF))
    tbright = tbright.clip(0, 1).fillna(0.0)

    depth_s = (depth - DEPTH_PPM_MIN) / max(1e-9, (DEPTH_PPM_MAX - DEPTH_PPM_MIN))
    depth_s = depth_s.clip(0, 1).fillna(0.0)

    period_s = 1.0 - (period - PERIOD_PIVOT_D) / max(1e-9, (PERIOD_MAX_D - PERIOD_PIVOT_D))
    period_s = period_s.clip(0, 1).fillna(0.0)

    # Vectorized triangular duration
    tri = np.zeros(len(df), dtype=float)
    mask = dur.notna() & (dur > DUR_MIN_H) & (dur < DUR_MAX_H)
    left = mask & (dur <= DUR_PEAK_H)
    right = mask & (dur > DUR_PEAK_H)
    tri[left]  = (dur[left]  - DUR_MIN_H) / max(1e-9, (DUR_PEAK_H - DUR_MIN_H))
    tri[right] = 1.0 - (dur[right] - DUR_PEAK_H) / max(1e-9, (DUR_MAX_H - DUR_PEAK_H))

    score_obs = (0.40 * tbright +
                 0.40 * depth_s +
                 0.15 * period_s +
                 0.05 * tri)

    df['score_m']        = score_m
    df['score_obs']      = score_obs
    df['priority_score'] = 0.5 * score_m + 0.5 * score_obs
    return df



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--toi_csv", required=True,
                   help="Path to full TOI CSV (ExoFOP-style).")

    # Fixed because of location and conditions [NOT USED]
    #p.add_argument("--site_lat", type=float, default=31.043416667, help="Observer latitude (deg).")
    #p.add_argument("--site_lon", type=float, default=-115.454763889, help="Observer longitude (deg, East positive).")
    #p.add_argument("--elevation_m", type=float, default=2780.0, help="Elevation (meters).")
    #p.add_argument("--min_alt_deg", type=float, default=30.0, help="Minimum altitude (deg) required through block.")
    #p.add_argument("--sun_limit", type=float, default=-18.0, help="Sun altitude threshold (deg) for 'dark'.")
    #p.add_argument("--dark_frac", type=float, default=0.8, help="Minimum fraction of block at/under sun_limit.")
    #p.add_argument("--cadence_min", type=float, default=5.0, help="Sampling cadence (minutes) within block.")
    #p.add_argument("--strict", action="store_true", help="Require 100% astronomical night (overrides --dark_frac to 1.0).")

    args = p.parse_args()

    # Read and keep original column names
    df = pd.read_csv(args.toi_csv)

    # Filter and rank
    cands = filter_candidates(df).copy()
    cands = rank_candidates(cands)
    cands = cands.sort_values('priority_score', ascending=False).reset_index(drop=True)

    cands.to_csv("m_dwarf_candidates_ranked.csv", index=False)

    print(f"Filtered + ranked candidates: {len(cands)}")
    cols_to_print = ["TOI", "TIC ID", "TESS Mag", "Period (days)", "priority_score", "Comments"]
    print(cands[cols_to_print].to_string(index=False))


if __name__ == "__main__":
    main()