# Sorting and ranking toi list
# Run
# python sorting_toi.py --toi_csv tois.csv --start 2025-11-17 --end 2025-11-30

import argparse
import math
import numpy as np
import pandas as pd
from typing import Tuple, List

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle
from astroplan import Observer
from astroplan.constraints import AltitudeConstraint
from astroplan.utils import time_grid_from_range

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


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {
        'TIC ID':'TIC','TESS Mag':'Tmag','Planet Num':'pl_num','TESS Disposition':'tess_disp',
        'TFOPWG Disposition':'tfop_disp','RA':'ra','Dec':'dec',
        'Epoch (BJD)':'epoch_bjd','Period (days)':'period','Duration (hours)':'duration_hr',
        'Depth (ppm)':'depth_ppm','Depth (mmag)':'depth_mmag','Planet Radius (R_Earth)':'rp_re',
        'Planet SNR':'tess_snr','Stellar Eff Temp (K)':'teff','Stellar log(g) (cm/s^2)':'logg',
        'Stellar Radius (R_Sun)':'rstar','Stellar Mass (M_Sun)':'mstar','Comments':'comments','Sectors':'sectors'
    }
    return df.rename(columns=rename_map)


def filter_candidates(df: pd.DataFrame) -> pd.DataFrame:
    gate = (
        df['tess_disp'].isin(TESS_DISP_OK) &
        df['teff'].between(TEFF_MIN, TEFF_MAX, inclusive='both') &
        (df['logg'] >= LOGG_MIN) &
        (df['rstar'] <= RSTAR_MAX) &
        (df['Tmag'] <= TESS_MAG_MAX) &
        (df['depth_ppm'] >= DEPTH_PPM_MIN) &
        (df['period'] <= PERIOD_MAX_D) &
        (df['duration_hr'].between(DUR_MIN_H, DUR_MAX_H, inclusive='both')) &
        (~df['comments'].apply(flagged_comment))
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
    m_teff_score = (TEFF_MAX - df['teff']) / (TEFF_MAX - 2600)
    m_rad_score  = (RSTAR_MAX - df['rstar']) / (RSTAR_MAX - 0.1)
    score_m = 0.6*m_teff_score.map(clamp01) + 0.4*m_rad_score.map(clamp01)

    tbright_score = 1 - (df['Tmag'] - 12) / (15 - 12)
    depth_score   = (df['depth_ppm'] - DEPTH_PPM_MIN) / (7000 - DEPTH_PPM_MIN)
    period_score  = 1 - (df['period'] - 5) / (PERIOD_MAX_D - 5)

    score_obs = (0.40 * tbright_score.map(clamp01) +
                 0.40 * depth_score.map(clamp01) +
                 0.15 * period_score.map(clamp01) +
                 0.05 * df['duration_hr'].map(duration_tri))

    df['score_m'] = score_m
    df['score_obs'] = score_obs
    df['priority_score'] = 0.45*score_m + 0.55*score_obs
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--toi_csv", required=True, help="Path to full TOI CSV (ExoFOP-style).")

    # not yet used
    p.add_argument("--start", required=True, help="UTC start date (YYYY-MM-DD).")
    p.add_argument("--end", required=True, help="UTC end date (YYYY-MM-DD).")

    # Fixed because of location and conditions
    p.add_argument("--site_lat", type=float, default=31.043416667, help="Observer latitude (deg).")
    p.add_argument("--site_lon", type=float, default=-115.454763889, help="Observer longitude (deg, East positive).")
    p.add_argument("--elevation_m", type=float, default=2780.0, help="Elevation (meters).")
    p.add_argument("--min_alt_deg", type=float, default=30.0, help="Minimum altitude (deg) required through block.")
    p.add_argument("--sun_limit", type=float, default=-18.0, help="Sun altitude threshold (deg) for 'dark'.")
    p.add_argument("--dark_frac", type=float, default=0.8, help="Minimum fraction of block at/under sun_limit.")
    p.add_argument("--cadence_min", type=float, default=5.0, help="Sampling cadence (minutes) within block.")
    p.add_argument("--strict", action="store_true", help="Require 100% astronomical night (overrides --dark_frac to 1.0).")

    args = p.parse_args()

    # Read, normalize, keep essentials
    raw = pd.read_csv(args.toi_csv)
    df = normalize_columns(raw)

    cols = ['TOI','TIC','tess_disp','Tmag','ra','dec','epoch_bjd','period','duration_hr','depth_ppm','teff','logg','rstar','comments']
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns in CSV: {missing}")

    df = df[cols].dropna(subset=['TOI','TIC','ra','dec','epoch_bjd','period','duration_hr','Tmag','teff','logg','rstar','tess_disp'])

    # Filter and rank
    cands = filter_candidates(df).copy()
    cands = rank_candidates(cands)
    cands = cands.sort_values('priority_score', ascending=False).reset_index(drop=True)

    cands.to_csv("m_dwarf_candidates_ranked.csv", index=False)

    print(f"Filtered + ranked candidates: {len(cands)}")
    cols_to_print = ["TOI", "TIC", "Tmag", "period", "priority_score", "comments"]
    print(cands[cols_to_print].to_string(index=False))


if __name__ == "__main__":
    main()