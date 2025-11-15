# Sorting and ranking toi list
# Run:
# python sorting_toi.py --toi_csv tois.csv

import argparse
import pandas as pd


# Defaults
TEFF_MIN, TEFF_MAX = 2400, 4000        # K
LOGG_MIN = 4.3
RSTAR_MAX = 0.7                        # R_sun
TESS_DISP_OK = {"PC"}
TESS_MAG_MAX = 14.0
DEPTH_PPM_MIN = 2000.0
PERIOD_MAX_D = 15.0
DUR_MIN_H, DUR_MAX_H = 0.5, 5.0

COMMENT_FLAGS = [
    'v-shaped','v shaped','eb','eclips','odd-even','sb2','binary','fp','false',
    'retired','low snr','contamin','centroid offset'
]

def flagged_comment(txt: str) -> bool:
    if pd.isna(txt): return False
    s = str(txt).lower()
    return any(flag in s for flag in COMMENT_FLAGS)


def filter_candidates(df: pd.DataFrame) -> pd.DataFrame:
    gate = (
        df['TESS Disposition'].isin(TESS_DISP_OK) &
        df['Stellar Eff Temp (K)'].between(TEFF_MIN, TEFF_MAX, inclusive='both') &
        (df['Stellar log(g) (cm/s^2)'] >= LOGG_MIN) &
        (df['Stellar Radius (R_Sun)'] <= RSTAR_MAX) &
        (df['TESS Mag'] <= TESS_MAG_MAX) &
        (df['Depth (ppm)'] >= DEPTH_PPM_MIN) &
        (df['Period (days)'] <= PERIOD_MAX_D) &
        (df['Duration (hours)'].between(DUR_MIN_H, DUR_MAX_H, inclusive='both')) &
        (~df['Comments'].apply(flagged_comment))
    )
    return df[gate].copy()


def clamp01(x: float) -> float:
    return float(max(0.0, min(1.0, x))) if not pd.isna(x) else 0.0


def duration_tri(h: float) -> float:
    if pd.isna(h): return 0.0
    if h <= DUR_MIN_H or h >= DUR_MAX_H: return 0.0
    if h <= 2.0: return (h - DUR_MIN_H) / (2.0 - DUR_MIN_H)
    return 1.0 - (h - 2.0) / (DUR_MAX_H - 2.0)


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    m_teff_score = (TEFF_MAX - df['Stellar Eff Temp (K)']) / (TEFF_MAX - 2600)
    m_rad_score  = (RSTAR_MAX - df['Stellar Radius (R_Sun)']) / (RSTAR_MAX - 0.1)
    score_m = 0.6*m_teff_score.map(clamp01) + 0.4*m_rad_score.map(clamp01)

    tbright_score = 1 - (df['TESS Mag'] - 12) / (15 - 12)
    depth_score   = (df['Depth (ppm)'] - DEPTH_PPM_MIN) / (7000 - DEPTH_PPM_MIN)
    period_score  = 1 - (df['Period (days)'] - 5) / (PERIOD_MAX_D - 5)

    score_obs = (0.40 * tbright_score.map(clamp01) +
                 0.40 * depth_score.map(clamp01) +
                 0.15 * period_score.map(clamp01) +
                 0.05 * df['Duration (hours)'].map(duration_tri))

    df['score_m'] = score_m
    df['score_obs'] = score_obs
    df['priority_score'] = 0.5 * score_m + 0.5 * score_obs
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--toi_csv", required=True, help="Path to full TOI CSV (ExoFOP-style).")

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

    # Date Modified is now preserved in the output CSV
    cands.to_csv("m_dwarf_candidates_ranked.csv", index=False)

    print(f"Filtered + ranked candidates: {len(cands)}")
    cols_to_print = ["TOI", "TIC ID", "TESS Mag", "Period (days)", "priority_score", "Comments"]
    print(cands[cols_to_print].to_string(index=False))


if __name__ == "__main__":
    main()