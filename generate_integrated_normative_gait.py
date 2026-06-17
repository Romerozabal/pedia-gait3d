#!/usr/bin/env python3
"""Generate pediatric normative gait tables and figures.

This script is intentionally self-contained inside
`PediaGait3D/`.

It combines:
  - Van Hamme et al. age/speed discrete feature models;
  - Schwartz et al. speed-binned angular and spatiotemporal references;
  - a Python implementation of the D2W gait-generator style angular-to-spatial
    transform.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


G = 9.81
SPEED_GROUP_RANGES = {
    "Very Slow": (0.090, 0.226),
    "Slow": (0.234, 0.359),
    "Free": (0.372, 0.482),
    "Fast": (0.502, 0.636),
    "Very Fast": (0.638, 0.844),
}
SPEED_GROUP_MEANS = {
    "Very Slow": 0.172,
    "Slow": 0.290,
    "Free": 0.429,
    "Fast": 0.560,
    "Very Fast": 0.694,
}
SPEED_GROUP_ORDER = ["Very Slow", "Slow", "Free", "Fast", "Very Fast"]
PROFILE_LABELS = {
    "trunk_tilt": "Trunk tilt",
    "trunk_obliquity": "Trunk obliquity",
    "trunk_rotation": "Trunk rotation",
    "pelvis_tilt": "Pelvis tilt",
    "pelvis_obliquity": "Pelvis obliquity",
    "pelvis_rotation": "Pelvis rotation",
    "hip_flex": "Hip flexion",
    "hip_adduction": "Hip adduction",
    "hip_rotation": "Hip rotation",
    "knee_flex": "Knee flexion",
    "knee_adduction": "Knee adduction",
    "knee_rotation": "Knee rotation",
    "ankle_dorsi": "Ankle dorsiflexion",
    "foot_progression": "Foot progression",
}
SCHWARTZ_PROFILE_ORDER = [
    "trunk_tilt",
    "trunk_obliquity",
    "trunk_rotation",
    "pelvis_tilt",
    "pelvis_obliquity",
    "pelvis_rotation",
    "hip_flex",
    "hip_adduction",
    "hip_rotation",
    "knee_flex",
    "knee_adduction",
    "knee_rotation",
    "ankle_dorsi",
    "foot_progression",
]
KINEMATIC_ROTATION_PROFILES = {
    "pelvis_roll_x": "pelvis_obliquity",
    "pelvis_pitch_y": "pelvis_tilt",
    "pelvis_yaw_z": "pelvis_rotation",
    "hip_roll_x": "hip_adduction",
    "hip_pitch_y": "hip_flex",
    "hip_yaw_z": "hip_rotation",
    "knee_roll_x": "knee_adduction",
    "knee_pitch_y": "knee_flex",
    "knee_yaw_z": "knee_rotation",
    "ankle_pitch_y": "ankle_dorsi",
    "terminal_foot_yaw_z": "foot_progression",
}
KINEMATIC_ROTATION_SIGNS = {
    "pelvis_obliquity": 1.0,
    "pelvis_tilt": -1.0,
    "pelvis_rotation": 1.0,
    "hip_adduction": -1.0,
    "hip_flex": 1.0,
    "hip_rotation": 1.0,
    "knee_adduction": -1.0,
    "knee_flex": -1.0,
    "knee_rotation": 1.0,
    "ankle_dorsi": 1.0,
    "foot_progression": 1.0,
}
D2W_SEGMENT_PROPORTIONS = {
    "leg_origin_offset": (0.0, 0.0, -0.530),
    "pelvis_to_hip": (0.0, 0.191 / 2, -0.01),
    "thigh": (0.0, 0.0, -0.246),
    "shank": (0.0, 0.0, -0.246),
    "forefoot": (-0.122 / 1.5, 0.0, -0.05),
    "toe": (-0.122, 0.0, -0.045),
    "heel": (0.03 / 1.5, 0.0, -0.05),
}
SUBJECT_BLUE = "#0072B2"
OUTPUT_ROOT_NAME = "outputs"
ANIMATION_FPS = 20


@dataclass(frozen=True)
class PatientSpec:
    patient_id: str
    age_y: float
    height_m: float
    speed_m_s: float
    leg_length_m: float | None = None


@dataclass(frozen=True)
class ContactModelConfig:
    enabled: bool = True
    pattern: str = "heel_strike"
    stance_start_percent: float = 0.0
    toe_off_percent: float = 60.0
    heel_to_flat_percent: float = 12.0
    forefoot_loading_percent: float = 45.0
    min_swing_clearance_m: float = 0.025


@dataclass(frozen=True)
class GenerationOptions:
    export_world_spatial: bool = True
    export_d2w_centered: bool = True
    schwartz_speed_interpolation: str = "linear"


@dataclass(frozen=True)
class SchwartzSpeedInterpolationPlan:
    mode: str
    lower_group: str
    upper_group: str
    lower_v_star: float
    upper_v_star: float
    alpha: float
    status: str

    @property
    def group_label(self) -> str:
        if self.lower_group == self.upper_group:
            return self.lower_group
        return f"{self.lower_group}_to_{self.upper_group}"


def package_dir() -> Path:
    return Path(__file__).resolve().parent


def model_dir() -> Path:
    return package_dir() / "models"


def default_out_dir(patient_id: str) -> Path:
    safe_id = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in patient_id).strip("_")
    if not safe_id:
        safe_id = "subject"
    return package_dir() / OUTPUT_ROOT_NAME / safe_id


def resolve_package_out_dir(out_dir: Path | None, patient_id: str) -> Path:
    resolved = default_out_dir(patient_id) if out_dir is None else (package_dir() / out_dir if not out_dir.is_absolute() else out_dir).resolve()
    package_root = package_dir().resolve()
    if not resolved.is_relative_to(package_root):
        raise ValueError(f"Output directory must stay inside the package directory: {resolved}")
    return resolved


def dimensionless_speed(speed_m_s: float, leg_length_m: float) -> float:
    return speed_m_s / math.sqrt(G * leg_length_m)


def schwartz_speed_group(v_star: float) -> tuple[str, str]:
    for group, (lo, hi) in SPEED_GROUP_RANGES.items():
        if lo <= v_star <= hi:
            return group, "inside_reported_range"
    closest = min(SPEED_GROUP_MEANS, key=lambda g: abs(v_star - SPEED_GROUP_MEANS[g]))
    return closest, "outside_reported_range_closest_group_used"


def schwartz_speed_interpolation_plan(v_star: float, mode: str = "linear") -> SchwartzSpeedInterpolationPlan:
    mode = mode.lower().strip()
    if mode not in {"linear", "nearest"}:
        raise ValueError("schwartz_speed_interpolation must be 'linear' or 'nearest'.")

    nearest_group, nearest_status = schwartz_speed_group(v_star)
    if mode == "nearest":
        v_nearest = SPEED_GROUP_MEANS[nearest_group]
        return SchwartzSpeedInterpolationPlan(
            mode=mode,
            lower_group=nearest_group,
            upper_group=nearest_group,
            lower_v_star=v_nearest,
            upper_v_star=v_nearest,
            alpha=0.0,
            status=f"{nearest_status}_nearest_waveform",
        )

    for group in SPEED_GROUP_ORDER:
        group_v = SPEED_GROUP_MEANS[group]
        if abs(v_star - group_v) <= 1.0e-12:
            return SchwartzSpeedInterpolationPlan(
                mode=mode,
                lower_group=group,
                upper_group=group,
                lower_v_star=group_v,
                upper_v_star=group_v,
                alpha=0.0,
                status="exact_speed_group_mean",
            )

    first = SPEED_GROUP_ORDER[0]
    last = SPEED_GROUP_ORDER[-1]
    first_v = SPEED_GROUP_MEANS[first]
    last_v = SPEED_GROUP_MEANS[last]
    if v_star <= first_v:
        return SchwartzSpeedInterpolationPlan(
            mode=mode,
            lower_group=first,
            upper_group=first,
            lower_v_star=first_v,
            upper_v_star=first_v,
            alpha=0.0,
            status="below_reference_mean_clamped_to_very_slow",
        )
    if v_star >= last_v:
        return SchwartzSpeedInterpolationPlan(
            mode=mode,
            lower_group=last,
            upper_group=last,
            lower_v_star=last_v,
            upper_v_star=last_v,
            alpha=0.0,
            status="above_reference_mean_clamped_to_very_fast",
        )

    for lower, upper in zip(SPEED_GROUP_ORDER[:-1], SPEED_GROUP_ORDER[1:]):
        lower_v = SPEED_GROUP_MEANS[lower]
        upper_v = SPEED_GROUP_MEANS[upper]
        if lower_v <= v_star <= upper_v:
            alpha = (v_star - lower_v) / (upper_v - lower_v)
            return SchwartzSpeedInterpolationPlan(
                mode=mode,
                lower_group=lower,
                upper_group=upper,
                lower_v_star=lower_v,
                upper_v_star=upper_v,
                alpha=float(alpha),
                status="linear_interpolation_between_group_means",
            )

    # Floating point fall-through should not be reached, but keep a deterministic fallback.
    v_nearest = SPEED_GROUP_MEANS[nearest_group]
    return SchwartzSpeedInterpolationPlan(
        mode=mode,
        lower_group=nearest_group,
        upper_group=nearest_group,
        lower_v_star=v_nearest,
        upper_v_star=v_nearest,
        alpha=0.0,
        status="fallback_nearest_waveform",
    )


def rotation_matrix(roll: float, pitch: float, yaw: float) -> np.ndarray:
    """D2W convention: R = Rz * Ry * Rx, angles in radians."""
    cr, sr = math.cos(roll), math.sin(roll)
    cp, sp = math.cos(pitch), math.sin(pitch)
    cy, sy = math.cos(yaw), math.sin(yaw)
    rx = np.array([[1, 0, 0], [0, cr, -sr], [0, sr, cr]], dtype=float)
    ry = np.array([[cp, 0, sp], [0, 1, 0], [-sp, 0, cp]], dtype=float)
    rz = np.array([[cy, -sy, 0], [sy, cy, 0], [0, 0, 1]], dtype=float)
    return rz @ ry @ rx


def circular_gaussian(data: np.ndarray, window_size: int) -> np.ndarray:
    if len(data) == 0:
        return data
    ws = window_size if window_size % 2 else window_size + 1
    half = ws // 2
    sigma = ws / 3.0
    x = np.arange(ws) - half
    kernel = np.exp(-0.5 * (x * x) / (sigma * sigma))
    kernel = kernel / kernel.sum()
    out = np.zeros_like(data, dtype=float)
    for i in range(len(data)):
        acc = np.zeros(data.shape[1], dtype=float)
        for k, j in zip(kernel, range(-half, half + 1)):
            acc += k * data[(i + j) % len(data)]
        out[i] = acc
    return out


def interpolate_profile(df: pd.DataFrame, profile_name: str, stat: str, group: str, n_points: int) -> np.ndarray:
    g = df[
        df["profile_name"].eq(profile_name)
        & df["stat"].eq(stat)
        & df["schwartz_speed_group"].eq(group)
    ].sort_values("cycle_percent")
    if g.empty:
        raise ValueError(f"Missing Schwartz profile {profile_name} {stat} {group}.")
    source_x = g["cycle_percent"].to_numpy(float)
    source_y = g["value_deg"].to_numpy(float)
    target_x = np.linspace(0, 100, n_points, endpoint=False)
    return np.interp(target_x, source_x, source_y)


def interpolate_profile_for_speed(
    df: pd.DataFrame,
    profile_name: str,
    stat: str,
    speed_plan: SchwartzSpeedInterpolationPlan,
    n_points: int,
) -> np.ndarray:
    lower = interpolate_profile(df, profile_name, stat, speed_plan.lower_group, n_points)
    if speed_plan.lower_group == speed_plan.upper_group:
        return lower
    upper = interpolate_profile(df, profile_name, stat, speed_plan.upper_group, n_points)
    return (1.0 - speed_plan.alpha) * lower + speed_plan.alpha * upper


def build_schwartz_profiles(speed_plan: SchwartzSpeedInterpolationPlan, n_points: int) -> pd.DataFrame:
    src = pd.read_csv(model_dir() / "lam_schwartz_2008_joint_rotation_reference_long.csv")
    missing = [p for p in SCHWARTZ_PROFILE_ORDER if p not in set(src["profile_name"])]
    if missing:
        raise ValueError(f"Missing expected Schwartz joint-rotation profiles in model table: {missing}")
    target_percent = np.linspace(0, 100, n_points, endpoint=False)
    rows = []
    for profile in SCHWARTZ_PROFILE_ORDER:
        meta = src[src["profile_name"].eq(profile)].iloc[0].to_dict()
        for stat in ["-1sd", "mean", "+1sd"]:
            vals = interpolate_profile_for_speed(src, profile, stat, speed_plan, n_points)
            for pct, val in zip(target_percent, vals):
                rows.append(
                    {
                        "profile_name": profile,
                        "segment": meta.get("segment", ""),
                        "anatomical_motion": meta.get("anatomical_motion", ""),
                        "d2w_rotation_axis": meta.get("d2w_rotation_axis", ""),
                        "schwartz_angle_label": meta.get("schwartz_angle_label", ""),
                        "schwartz_speed_group": speed_plan.group_label,
                        "schwartz_speed_lower_group": speed_plan.lower_group,
                        "schwartz_speed_upper_group": speed_plan.upper_group,
                        "schwartz_speed_lower_v_star": speed_plan.lower_v_star,
                        "schwartz_speed_upper_v_star": speed_plan.upper_v_star,
                        "schwartz_speed_interpolation_alpha": speed_plan.alpha,
                        "schwartz_speed_interpolation_mode": speed_plan.mode,
                        "schwartz_speed_interpolation_status": speed_plan.status,
                        "stat": stat,
                        "cycle_percent": pct,
                        "value_deg": val,
                        "source": meta.get("source", "Schwartz et al. 2008 Excel addendum"),
                        "reference_note": meta.get(
                            "reference_note",
                            "Speed-binned waveform reference; not age-specific to 3-5 years.",
                        ),
                    }
                )
    return pd.DataFrame(rows)


def build_d2w_kinematic_model_proportions(height_m: float) -> pd.DataFrame:
    rows = []
    for segment, vec in D2W_SEGMENT_PROPORTIONS.items():
        x, y, z = vec
        rows.append(
            {
                "segment_vector": segment,
                "proportion_x_per_height": x,
                "proportion_y_per_height": y,
                "proportion_z_per_height": z,
                "x_m": x * height_m,
                "y_m": y * height_m,
                "z_m": z * height_m,
                "source": "D2W gait_generator.cpp segment proportions scaled by patient height",
            }
        )
    return pd.DataFrame(rows)


def van_hamme_prediction(row: pd.Series, age_y: float, v_star: float) -> float:
    return (
        float(row["a_age"]) * age_y
        + float(row["b_speed"]) * v_star
        + float(row["c_age_speed"]) * age_y * v_star
        + float(row["d_intercept"])
    )


def build_van_hamme_features(age_y: float, v_star: float) -> pd.DataFrame:
    coeff = pd.read_csv(model_dir() / "lam_van_hamme_2015_regression_coefficients.csv")
    base_rows = []
    values: dict[tuple[str, str], tuple[float, float, str]] = {}
    for _, row in coeff.iterrows():
        mean = van_hamme_prediction(row, age_y, v_star)
        sd = float(row["residual_sd"])
        key = (row["quantity"], row["code"])
        values[key] = (mean, sd, row["description"])
        base_rows.append(
            {
                "metric_name": f"{row['code']}_{row['quantity']}",
                "quantity": row["quantity"],
                "van_hamme_codes": row["code"],
                "unit": "deg" if row["quantity"] == "peak_value" else "%",
                "normal_mean": mean,
                "normal_sd": sd,
                "normal_low": mean - sd,
                "normal_high": mean + sd,
                "reference_feature_definition": row["description"],
                "source": "Van Hamme et al. 2015 supplementary regression coefficients",
                "normal_model": "Y = a*age + b*dimensionless_speed + c*age*dimensionless_speed + d",
            }
        )

    def get(quantity: str, code: str) -> tuple[float, float, str]:
        return values[(quantity, code)]

    derived = []
    derived_specs = [
        ("ankle_dorsi_peak_max_angle_deg", "deg", "peak_value", "A_A2", "Max dorsiflexion (26-50%)"),
        ("ankle_dorsi_peak_max_time_percent", "%", "time_percent", "A_A2", "Timing of max dorsiflexion (26-50%)"),
        ("knee_flex_peak_max_angle_deg", "deg", "peak_value", "A_K3", "Max knee flexion (56-100%)"),
        ("knee_flex_peak_max_time_percent", "%", "time_percent", "A_K3", "Timing of max knee flexion (56-100%)"),
        ("hip_flex_peak_max_angle_deg", "deg", "peak_value", "A_H3", "Late-cycle max hip flexion (70-100%)"),
        ("hip_flex_peak_max_time_percent", "%", "time_percent", "A_H3", "Timing of late-cycle max hip flexion (70-100%)"),
    ]
    for metric, unit, quantity, code, desc in derived_specs:
        mean, sd, _ = get(quantity, code)
        derived.append((metric, unit, code, mean, sd, desc))

    rom_specs = [
        ("ankle_dorsi_rom_deg", "A_A2", "A_A3", "Ankle sagittal ROM from dorsiflexion peak minus late plantarflexion peak"),
        ("knee_flex_rom_deg", "A_K3", "A_K2", "Knee sagittal ROM from swing flexion peak minus mid-cycle minimum"),
        ("hip_flex_rom_deg", "A_H3", "A_H2", "Hip sagittal ROM from late flexion peak minus mid-cycle minimum"),
    ]
    for metric, high_code, low_code, desc in rom_specs:
        high_mean, high_sd, _ = get("peak_value", high_code)
        low_mean, low_sd, _ = get("peak_value", low_code)
        mean = high_mean - low_mean
        sd = math.sqrt(high_sd * high_sd + low_sd * low_sd)
        derived.append((metric, "deg", f"{high_code}-{low_code}", mean, sd, desc))

    derived_rows = []
    for metric, unit, codes, mean, sd, desc in derived:
        derived_rows.append(
            {
                "metric_name": metric,
                "quantity": "derived_manuscript_feature",
                "van_hamme_codes": codes,
                "unit": unit,
                "normal_mean": mean,
                "normal_sd": sd,
                "normal_low": mean - sd,
                "normal_high": mean + sd,
                "reference_feature_definition": desc,
                "source": "Van Hamme et al. 2015 supplementary regression coefficients",
                "normal_model": "Y = a*age + b*dimensionless_speed + c*age*dimensionless_speed + d",
            }
        )
    return pd.DataFrame(base_rows + derived_rows)


def build_spatiotemporal_reference(group: str, leg_length_m: float) -> pd.DataFrame:
    src = pd.read_csv(model_dir() / "lam_schwartz_2008_anova_reference_long.csv")
    src = src[src["schwartz_speed_group"].eq(group)].copy()
    rows = []
    for _, row in src.iterrows():
        metric = row["metric_name"]
        mean = float(row["mean"])
        sd = float(row["sd"])
        rows.append(row.to_dict() | {"converted_metric_name": "", "converted_mean": np.nan, "converted_sd": np.nan, "converted_unit": ""})
        if metric == "Dimensionless Step Length":
            rows.append(
                row.to_dict()
                | {
                    "converted_metric_name": "step_length_m",
                    "converted_mean": mean * leg_length_m,
                    "converted_sd": sd * leg_length_m,
                    "converted_unit": "m",
                }
            )
        elif metric == "Dimensionless Cadence":
            scale = math.sqrt(leg_length_m / G)
            cadence_steps_s = mean / scale
            cadence_sd_steps_s = sd / scale
            rows.append(
                row.to_dict()
                | {
                    "converted_metric_name": "cadence_steps_min",
                    "converted_mean": cadence_steps_s * 60.0,
                    "converted_sd": cadence_sd_steps_s * 60.0,
                    "converted_unit": "steps/min",
                }
            )
            rows.append(
                row.to_dict()
                | {
                    "converted_metric_name": "step_time_s",
                    "converted_mean": 1.0 / cadence_steps_s,
                    "converted_sd": np.nan,
                    "converted_unit": "s",
                }
            )
            rows.append(
                row.to_dict()
                | {
                    "converted_metric_name": "stride_time_s",
                    "converted_mean": 2.0 / cadence_steps_s,
                    "converted_sd": np.nan,
                    "converted_unit": "s",
                }
            )
    return pd.DataFrame(rows)


def smoothstep(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def periodic_interpolate(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return np.zeros_like(values, dtype=float)
    if len(idx) == len(values):
        return values.astype(float)
    n = len(values)
    x = np.arange(n)
    xp = np.r_[idx - n, idx, idx + n]
    fp = np.r_[values[idx], values[idx], values[idx]]
    return np.interp(x, xp, fp)


def contact_weights(cycle_percent: np.ndarray, cfg: ContactModelConfig) -> pd.DataFrame:
    pct = np.asarray(cycle_percent, dtype=float)
    w_heel = np.zeros_like(pct)
    w_forefoot = np.zeros_like(pct)
    w_toe = np.zeros_like(pct)
    stance = (pct >= cfg.stance_start_percent) & (pct <= cfg.toe_off_percent)

    if not cfg.enabled:
        w_heel[:] = 0.5
        w_forefoot[:] = 0.5
        stance[:] = True
    elif cfg.pattern == "heel_strike":
        idx = (pct >= cfg.stance_start_percent) & (pct <= cfg.heel_to_flat_percent)
        s1 = smoothstep((pct[idx] - cfg.stance_start_percent) / max(1.0e-6, cfg.heel_to_flat_percent - cfg.stance_start_percent))
        w_heel[idx] = 1.0 - 0.5 * s1
        w_forefoot[idx] = 0.5 * s1

        idx = (pct > cfg.heel_to_flat_percent) & (pct <= cfg.forefoot_loading_percent)
        w_heel[idx] = 0.5
        w_forefoot[idx] = 0.5

        idx = (pct > cfg.forefoot_loading_percent) & (pct <= cfg.toe_off_percent)
        s2 = smoothstep((pct[idx] - cfg.forefoot_loading_percent) / max(1.0e-6, cfg.toe_off_percent - cfg.forefoot_loading_percent))
        w_heel[idx] = 0.5 * (1.0 - s2)
        w_forefoot[idx] = 0.5 * (1.0 - s2)
        w_toe[idx] = s2
    elif cfg.pattern == "midfoot_strike":
        w_heel[stance] = 0.4
        w_forefoot[stance] = 0.6
    elif cfg.pattern == "forefoot_strike":
        w_forefoot[stance] = 0.7
        w_toe[stance] = 0.3
    else:
        raise ValueError(f"Unknown contact pattern: {cfg.pattern}")

    total = w_heel + w_forefoot + w_toe
    valid = total > 0.0
    w_heel[valid] /= total[valid]
    w_forefoot[valid] /= total[valid]
    w_toe[valid] /= total[valid]
    return pd.DataFrame(
        {
            "cycle_percent": pct,
            "is_stance": stance,
            "w_heel": w_heel,
            "w_forefoot": w_forefoot,
            "w_toe": w_toe,
            "contact_pattern": cfg.pattern,
        }
    )


def stride_time_from_spatiotemporal(spatiotemporal: pd.DataFrame) -> float | None:
    stride = spatiotemporal[spatiotemporal["converted_metric_name"].eq("stride_time_s")]
    if stride.empty:
        return None
    val = float(stride["converted_mean"].iloc[0])
    return val if np.isfinite(val) and val > 0 else None


def converted_spatiotemporal_value(spatiotemporal: pd.DataFrame, metric_name: str) -> float:
    row = spatiotemporal[spatiotemporal["converted_metric_name"].eq(metric_name)]
    if row.empty:
        return np.nan
    val = float(row["converted_mean"].iloc[0])
    return val if np.isfinite(val) else np.nan


def compute_raw_fk_points(
    schwartz_profiles: pd.DataFrame,
    height_m: float,
    n_points: int,
    angle_sign: float = 1.0,
) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    """Compute pelvis-frame FK points before ground/contact correction."""
    pct = np.linspace(0, 100, n_points, endpoint=False)
    mean = schwartz_profiles[schwartz_profiles["stat"].eq("mean")]
    missing = [profile for profile in KINEMATIC_ROTATION_PROFILES.values() if profile not in set(mean["profile_name"])]
    if missing:
        raise ValueError(f"Cannot build 3D trajectory; missing required kinematic rotation profiles: {missing}")

    def profile_rad(profile_name: str) -> np.ndarray:
        vals = mean[mean["profile_name"].eq(profile_name)].sort_values("cycle_percent")["value_deg"].to_numpy(float)
        if len(vals) == 0:
            return np.zeros(n_points)
        if len(vals) != n_points:
            raise ValueError(f"Profile {profile_name} has {len(vals)} points; expected {n_points}.")
        return np.deg2rad(angle_sign * KINEMATIC_ROTATION_SIGNS.get(profile_name, 1.0) * vals)

    zeros = np.zeros(n_points)
    pelvis_roll = profile_rad("pelvis_obliquity")
    pelvis_pitch = profile_rad("pelvis_tilt")
    pelvis_yaw = profile_rad("pelvis_rotation")
    hip_roll = profile_rad("hip_adduction")
    hip_pitch = profile_rad("hip_flex")
    hip_yaw = profile_rad("hip_rotation")
    knee_roll = profile_rad("knee_adduction")
    knee_pitch = profile_rad("knee_flex")
    knee_yaw = profile_rad("knee_rotation")
    ankle_roll = zeros
    ankle_pitch = profile_rad("ankle_dorsi")
    ankle_yaw = profile_rad("foot_progression")

    leg = np.array(D2W_SEGMENT_PROPORTIONS["leg_origin_offset"]) * height_m
    pelvis_dim = np.array(D2W_SEGMENT_PROPORTIONS["pelvis_to_hip"]) * height_m
    shank = np.array(D2W_SEGMENT_PROPORTIONS["shank"]) * height_m
    thigh = np.array(D2W_SEGMENT_PROPORTIONS["thigh"]) * height_m
    foot = np.array(D2W_SEGMENT_PROPORTIONS["forefoot"]) * height_m
    toe = np.array(D2W_SEGMENT_PROPORTIONS["toe"]) * height_m
    heel = np.array(D2W_SEGMENT_PROPORTIONS["heel"]) * height_m

    k_p = -leg
    pelvis_origin = np.zeros((n_points, 3))
    hip = np.zeros((n_points, 3))
    knee = np.zeros((n_points, 3))
    ankle = np.zeros((n_points, 3))
    forefoot = np.zeros((n_points, 3))
    toe_path = np.zeros((n_points, 3))
    heel_path = np.zeros((n_points, 3))

    for i in range(n_points):
        r_pelvis = rotation_matrix(pelvis_roll[i], pelvis_pitch[i], pelvis_yaw[i])
        r_hip = rotation_matrix(hip_roll[i], hip_pitch[i], hip_yaw[i])
        r_knee = rotation_matrix(knee_roll[i], knee_pitch[i], knee_yaw[i])
        r_ankle = rotation_matrix(ankle_roll[i], ankle_pitch[i], ankle_yaw[i])
        hip_rel = r_pelvis @ pelvis_dim
        knee_rel = r_pelvis @ r_hip @ thigh
        ankle_rel = r_pelvis @ r_hip @ r_knee @ shank
        foot_rel = r_pelvis @ r_hip @ r_knee @ r_ankle @ foot
        toe_rel = r_pelvis @ r_hip @ r_knee @ r_ankle @ toe
        heel_rel = r_pelvis @ r_hip @ r_knee @ r_ankle @ heel
        pelvis_origin[i] = np.zeros(3)
        hip[i] = k_p + hip_rel
        knee[i] = hip[i] + knee_rel
        ankle[i] = knee[i] + ankle_rel
        forefoot[i] = ankle[i] + foot_rel
        toe_path[i] = ankle[i] + toe_rel
        heel_path[i] = ankle[i] + heel_rel

    return (
        pct,
        {
            "pelvis_origin": pelvis_origin,
            "hip": hip,
            "knee": knee,
            "ankle": ankle,
            "forefoot": forefoot,
            "toe": toe_path,
            "heel": heel_path,
        },
        leg,
    )


def schwartz_mean_angles_rad(schwartz_profiles: pd.DataFrame, n_points: int, angle_sign: float = 1.0) -> dict[str, np.ndarray]:
    """Return the mean Schwartz kinematic inputs in radians."""
    mean = schwartz_profiles[schwartz_profiles["stat"].eq("mean")]
    missing = [profile for profile in KINEMATIC_ROTATION_PROFILES.values() if profile not in set(mean["profile_name"])]
    if missing:
        raise ValueError(f"Cannot build 3D trajectory; missing required kinematic rotation profiles: {missing}")

    def profile_rad(profile_name: str) -> np.ndarray:
        vals = mean[mean["profile_name"].eq(profile_name)].sort_values("cycle_percent")["value_deg"].to_numpy(float)
        if len(vals) == 0:
            return np.zeros(n_points)
        if len(vals) != n_points:
            raise ValueError(f"Profile {profile_name} has {len(vals)} points; expected {n_points}.")
        return np.deg2rad(angle_sign * KINEMATIC_ROTATION_SIGNS.get(profile_name, 1.0) * vals)

    zeros = np.zeros(n_points)
    return {
        "pelvis_roll": profile_rad("pelvis_obliquity"),
        "pelvis_pitch": profile_rad("pelvis_tilt"),
        "pelvis_yaw": profile_rad("pelvis_rotation"),
        "hip_roll": profile_rad("hip_adduction"),
        "hip_pitch": profile_rad("hip_flex"),
        "hip_yaw": profile_rad("hip_rotation"),
        "knee_roll": profile_rad("knee_adduction"),
        "knee_pitch": profile_rad("knee_flex"),
        "knee_yaw": profile_rad("knee_rotation"),
        "ankle_roll": zeros,
        "ankle_pitch": profile_rad("ankle_dorsi"),
        "ankle_yaw": profile_rad("foot_progression"),
    }


def compute_raw_bilateral_fk_points(
    schwartz_profiles: pd.DataFrame,
    height_m: float,
    n_points: int,
    angle_sign: float = 1.0,
) -> tuple[np.ndarray, dict[str, dict[str, np.ndarray] | np.ndarray], np.ndarray]:
    """Compute bilateral pelvis-frame FK before common-pelvis translation.

    The pelvis orientation is evaluated once at the global gait phase. The left
    limb uses the contralateral 50% limb phase, mirrored in the pelvis local
    frame before applying the global pelvis orientation.
    """
    pct = np.linspace(0, 100, n_points, endpoint=False)
    angles = schwartz_mean_angles_rad(schwartz_profiles, n_points, angle_sign)

    leg = np.array(D2W_SEGMENT_PROPORTIONS["leg_origin_offset"]) * height_m
    pelvis_dim = np.array(D2W_SEGMENT_PROPORTIONS["pelvis_to_hip"]) * height_m
    shank = np.array(D2W_SEGMENT_PROPORTIONS["shank"]) * height_m
    thigh = np.array(D2W_SEGMENT_PROPORTIONS["thigh"]) * height_m
    forefoot_vec = np.array(D2W_SEGMENT_PROPORTIONS["forefoot"]) * height_m
    toe_vec = np.array(D2W_SEGMENT_PROPORTIONS["toe"]) * height_m
    heel_vec = np.array(D2W_SEGMENT_PROPORTIONS["heel"]) * height_m

    k_p = -leg
    shift = n_points // 2
    raw: dict[str, dict[str, np.ndarray] | np.ndarray] = {
        "pelvis_origin": np.zeros((n_points, 3)),
        "R": {joint: np.zeros((n_points, 3)) for joint in ["hip", "knee", "ankle", "forefoot", "toe", "heel"]},
        "L": {joint: np.zeros((n_points, 3)) for joint in ["hip", "knee", "ankle", "forefoot", "toe", "heel"]},
    }

    for i in range(n_points):
        r_pelvis = rotation_matrix(angles["pelvis_roll"][i], angles["pelvis_pitch"][i], angles["pelvis_yaw"][i])
        for side, phase_shift, y_sign in [("R", 0, 1.0), ("L", shift, -1.0)]:
            j = (i + phase_shift) % n_points
            mirror = np.diag([1.0, y_sign, 1.0])
            r_hip = rotation_matrix(angles["hip_roll"][j], angles["hip_pitch"][j], angles["hip_yaw"][j])
            r_knee = rotation_matrix(angles["knee_roll"][j], angles["knee_pitch"][j], angles["knee_yaw"][j])
            r_ankle = rotation_matrix(angles["ankle_roll"][j], angles["ankle_pitch"][j], angles["ankle_yaw"][j])
            hip_rel = r_pelvis @ (mirror @ pelvis_dim)
            knee_rel = r_pelvis @ (mirror @ (r_hip @ thigh))
            ankle_rel = r_pelvis @ (mirror @ (r_hip @ r_knee @ shank))
            foot_chain = r_hip @ r_knee @ r_ankle
            forefoot_rel = r_pelvis @ (mirror @ (foot_chain @ forefoot_vec))
            toe_rel = r_pelvis @ (mirror @ (foot_chain @ toe_vec))
            heel_rel = r_pelvis @ (mirror @ (foot_chain @ heel_vec))
            side_arrays = raw[side]
            assert isinstance(side_arrays, dict)
            side_arrays["hip"][i] = k_p + hip_rel
            side_arrays["knee"][i] = side_arrays["hip"][i] + knee_rel
            side_arrays["ankle"][i] = side_arrays["knee"][i] + ankle_rel
            side_arrays["forefoot"][i] = side_arrays["ankle"][i] + forefoot_rel
            side_arrays["toe"][i] = side_arrays["ankle"][i] + toe_rel
            side_arrays["heel"][i] = side_arrays["ankle"][i] + heel_rel

    return pct, raw, leg


def bilateral_support_weights(cycle_percent: np.ndarray, toe_off_percent: float = 60.0) -> pd.DataFrame:
    """Support load-sharing schedule for a symmetric 50% phase-offset gait."""
    pct = np.asarray(cycle_percent, dtype=float)
    ds = max(0.0, toe_off_percent - 50.0)
    lambda_r = np.zeros_like(pct)
    lambda_l = np.zeros_like(pct)
    if ds <= 1.0e-6:
        lambda_r[pct < 50.0] = 1.0
        lambda_l[pct >= 50.0] = 1.0
    else:
        idx = pct < ds
        s = smoothstep(pct[idx] / ds)
        lambda_r[idx] = s
        lambda_l[idx] = 1.0 - s

        idx = (pct >= ds) & (pct < 50.0)
        lambda_r[idx] = 1.0

        idx = (pct >= 50.0) & (pct < 50.0 + ds)
        s = smoothstep((pct[idx] - 50.0) / ds)
        lambda_r[idx] = 1.0 - s
        lambda_l[idx] = s

        idx = pct >= 50.0 + ds
        lambda_l[idx] = 1.0

    return pd.DataFrame({"cycle_percent": pct, "lambda_R": lambda_r, "lambda_L": lambda_l})


def contact_point_from_weights(side_arrays: dict[str, np.ndarray], contact_df: pd.DataFrame) -> np.ndarray:
    w_heel = contact_df["w_heel"].to_numpy(float)[:, None]
    w_forefoot = contact_df["w_forefoot"].to_numpy(float)[:, None]
    w_toe = contact_df["w_toe"].to_numpy(float)[:, None]
    return w_heel * side_arrays["heel"] + w_forefoot * side_arrays["forefoot"] + w_toe * side_arrays["toe"]


def contact_anchor_from_weights(base_heel_anchor: np.ndarray, contact_df: pd.DataFrame, height_m: float) -> np.ndarray:
    """Build a ground rocker anchor; toe-off uses the forefoot anchor initially."""
    heel_vec = np.array(D2W_SEGMENT_PROPORTIONS["heel"]) * height_m
    forefoot_vec = np.array(D2W_SEGMENT_PROPORTIONS["forefoot"]) * height_m
    forefoot_offset = np.array([forefoot_vec[0] - heel_vec[0], 0.0, 0.0])

    w_heel = contact_df["w_heel"].to_numpy(float)[:, None]
    w_forefoot = contact_df["w_forefoot"].to_numpy(float)[:, None]
    w_toe = contact_df["w_toe"].to_numpy(float)[:, None]
    w_forefoot_anchor = w_forefoot + w_toe
    return (
        w_heel * base_heel_anchor
        + w_forefoot_anchor * (base_heel_anchor + forefoot_offset)
    )


def solve_overground_pelvis_translation(
    raw: dict[str, dict[str, np.ndarray] | np.ndarray],
    pct: np.ndarray,
    contact_cfg: ContactModelConfig,
    step_length_m: float,
    height_m: float,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Estimate one common pelvis translation from bilateral stance anchors."""
    raw_r = raw["R"]
    raw_l = raw["L"]
    assert isinstance(raw_r, dict) and isinstance(raw_l, dict)
    phase_r = pct
    phase_l = (pct + 50.0) % 100.0
    contact_r = contact_weights(phase_r, contact_cfg)
    contact_l = contact_weights(phase_l, contact_cfg)
    support = bilateral_support_weights(pct, contact_cfg.toe_off_percent)

    c_r = contact_point_from_weights(raw_r, contact_r)
    c_l = contact_point_from_weights(raw_l, contact_l)
    lambda_r = support["lambda_R"].to_numpy(float)[:, None]
    lambda_l = support["lambda_L"].to_numpy(float)[:, None]
    ds = max(0.0, contact_cfg.toe_off_percent - 50.0)
    forward_sign = -1.0 if c_r[0, 0] < 0.0 else 1.0

    heel_anchor_r = np.zeros_like(c_r)
    heel_anchor_l = np.zeros_like(c_l)
    heel_anchor_l[:, 0] = np.where(pct < ds, -forward_sign * step_length_m, forward_sign * step_length_m)

    r_single = (lambda_r.ravel() > 0.999) & contact_r["is_stance"].to_numpy(bool)
    l_single = (lambda_l.ravel() > 0.999) & contact_l["is_stance"].to_numpy(bool)
    heel_anchor_r[:, 1] = float(np.median(c_r[r_single, 1])) if np.any(r_single) else float(np.median(c_r[:, 1]))
    heel_anchor_l[:, 1] = float(np.median(c_l[l_single, 1])) if np.any(l_single) else float(np.median(c_l[:, 1]))

    anchor_r = contact_anchor_from_weights(heel_anchor_r, contact_r, height_m)
    anchor_l = contact_anchor_from_weights(heel_anchor_l, contact_l, height_m)

    p_r = anchor_r - c_r
    p_l = anchor_l - c_l
    total = np.maximum(lambda_r + lambda_l, 1.0e-9)
    pelvis_origin_world = (lambda_r * p_r + lambda_l * p_l) / total
    pelvis_origin_world[:, 1] = circular_gaussian(pelvis_origin_world[:, [1]], 31).ravel() * 0.4

    rows = []
    for side, phase, contact_df, contact_point, anchor, heel_anchor, candidate, lam in [
        ("R", phase_r, contact_r, c_r, anchor_r, heel_anchor_r, p_r, lambda_r.ravel()),
        ("L", phase_l, contact_l, c_l, anchor_l, heel_anchor_l, p_l, lambda_l.ravel()),
    ]:
        is_stance = contact_df["is_stance"].to_numpy(bool)
        anchor_error = np.linalg.norm(pelvis_origin_world + contact_point - anchor, axis=1)
        for i, cycle_percent in enumerate(pct):
            rows.append(
                {
                    "side": side,
                    "cycle_percent": cycle_percent,
                    "source_phase_percent": phase[i],
                    "is_stance": bool(is_stance[i]),
                    "lambda_support": lam[i],
                    "w_heel": float(contact_df["w_heel"].iloc[i]),
                    "w_forefoot": float(contact_df["w_forefoot"].iloc[i]),
                    "w_toe": float(contact_df["w_toe"].iloc[i]),
                    "heel_anchor_x_m": heel_anchor[i, 0],
                    "heel_anchor_y_m": heel_anchor[i, 1],
                    "heel_anchor_z_m": heel_anchor[i, 2],
                    "anchor_x_m": anchor[i, 0],
                    "anchor_y_m": anchor[i, 1],
                    "anchor_z_m": anchor[i, 2],
                    "contact_x_m": contact_point[i, 0],
                    "contact_y_m": contact_point[i, 1],
                    "contact_z_m": contact_point[i, 2],
                    "pelvis_candidate_x_m": candidate[i, 0],
                    "pelvis_candidate_y_m": candidate[i, 1],
                    "pelvis_candidate_z_m": candidate[i, 2],
                    "pelvis_world_origin_x_m": pelvis_origin_world[i, 0],
                    "pelvis_world_origin_y_m": pelvis_origin_world[i, 1],
                    "pelvis_world_origin_z_m": pelvis_origin_world[i, 2],
                    "contact_anchor_error_m": anchor_error[i] if lam[i] > 0.0 else np.nan,
                    "active_contact_anchor_error_m": lam[i] * anchor_error[i] if lam[i] > 0.0 else np.nan,
                    "contact_pattern": contact_cfg.pattern,
                }
            )
    return pelvis_origin_world, pd.DataFrame(rows)


def apply_bilateral_world_translation(
    raw: dict[str, dict[str, np.ndarray] | np.ndarray],
    pelvis_origin_world: np.ndarray,
    leg: np.ndarray,
) -> dict[str, dict[str, np.ndarray] | np.ndarray]:
    world: dict[str, dict[str, np.ndarray] | np.ndarray] = {
        "pelvis": pelvis_origin_world + np.array([0.0, 0.0, -leg[2] * 1.1]),
        "R": {},
        "L": {},
    }
    for side in ["R", "L"]:
        side_raw = raw[side]
        side_world = world[side]
        assert isinstance(side_raw, dict) and isinstance(side_world, dict)
        for joint, arr in side_raw.items():
            side_world[joint] = arr + pelvis_origin_world
    return world


def project_world_to_pelvis_fixed_treadmill(world: dict[str, dict[str, np.ndarray] | np.ndarray]) -> dict[str, dict[str, np.ndarray] | np.ndarray]:
    pelvis = world["pelvis"]
    assert isinstance(pelvis, np.ndarray)
    pelvis_x = pelvis[:, 0].copy()
    treadmill: dict[str, dict[str, np.ndarray] | np.ndarray] = {"pelvis": pelvis.copy(), "R": {}, "L": {}}
    treadmill_pelvis = treadmill["pelvis"]
    assert isinstance(treadmill_pelvis, np.ndarray)
    treadmill_pelvis[:, 0] = 0.0
    for side in ["R", "L"]:
        side_world = world[side]
        side_treadmill = treadmill[side]
        assert isinstance(side_world, dict) and isinstance(side_treadmill, dict)
        for joint, arr in side_world.items():
            shifted = arr.copy()
            shifted[:, 0] -= pelvis_x
            side_treadmill[joint] = shifted
    return treadmill


def convert_treadmill_to_d2w_convention(treadmill: dict[str, dict[str, np.ndarray] | np.ndarray]) -> dict[str, dict[str, np.ndarray] | np.ndarray]:
    d2w: dict[str, dict[str, np.ndarray] | np.ndarray] = {"pelvis": treadmill["pelvis"].copy(), "R": {}, "L": {}}
    pelvis = d2w["pelvis"]
    assert isinstance(pelvis, np.ndarray)
    pelvis[:, 0] *= -1.0
    pelvis[:, 1] *= -0.6
    for side in ["R", "L"]:
        side_treadmill = treadmill[side]
        side_d2w = d2w[side]
        assert isinstance(side_treadmill, dict) and isinstance(side_d2w, dict)
        for joint, arr in side_treadmill.items():
            converted = arr.copy()
            converted[:, 0] *= -1.0
            converted[:, 1] *= -1.0
            side_d2w[joint] = converted
    return d2w


def bilateral_arrays_to_spatial_df(
    arrays: dict[str, dict[str, np.ndarray] | np.ndarray],
    pct: np.ndarray,
    source: str,
) -> pd.DataFrame:
    rows = []
    pelvis = arrays["pelvis"]
    assert isinstance(pelvis, np.ndarray)
    for i, cycle_percent in enumerate(pct):
        point = pelvis[i]
        rows.append(
            {
                "side": "C",
                "joint": "pelvis",
                "cycle_percent": cycle_percent,
                "x_m": point[0],
                "y_m": point[1],
                "z_m": point[2],
                "rotation_inputs": ";".join(
                    f"{axis}={profile},sign={KINEMATIC_ROTATION_SIGNS.get(profile, 1.0):g}"
                    for axis, profile in KINEMATIC_ROTATION_PROFILES.items()
                ),
                "segment_proportions": "D2W height-scaled vectors in d2w_kinematic_model_proportions.csv",
                "source": source,
            }
        )
    for side in ["R", "L"]:
        side_arrays = arrays[side]
        assert isinstance(side_arrays, dict)
        for i, cycle_percent in enumerate(pct):
            for joint, arr in side_arrays.items():
                point = arr[i]
                rows.append(
                    {
                        "side": side,
                        "joint": joint,
                        "cycle_percent": cycle_percent,
                        "x_m": point[0],
                        "y_m": point[1],
                        "z_m": point[2],
                        "rotation_inputs": ";".join(
                            f"{axis}={profile},sign={KINEMATIC_ROTATION_SIGNS.get(profile, 1.0):g}"
                            for axis, profile in KINEMATIC_ROTATION_PROFILES.items()
                        ),
                        "segment_proportions": "D2W height-scaled vectors in d2w_kinematic_model_proportions.csv",
                        "source": source,
                    }
                )
    return pd.DataFrame(rows)


def build_d2w_ankle_command_trajectories(spatial: pd.DataFrame) -> pd.DataFrame:
    """Robot-oriented ankle command table with centered ankle X only."""
    ankle = spatial[spatial["joint"].eq("ankle") & spatial["side"].isin(["R", "L"])].copy()
    if ankle.empty:
        return ankle
    ankle_x_offset = 0.5 * (float(ankle["x_m"].max()) + float(ankle["x_m"].min()))
    ankle["x_m"] = ankle["x_m"] - ankle_x_offset
    ankle["command_note"] = "Ankle-only centered X command; do not use as full-body geometric trajectory."
    return ankle


def apply_phase_contact_model(
    arrays: dict[str, np.ndarray],
    contact_df: pd.DataFrame,
    min_swing_clearance_m: float,
    walking_speed_m_s: float,
    stride_time_s: float | None,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    heel = arrays["heel"]
    forefoot = arrays["forefoot"]
    toe = arrays["toe"]
    w_heel = contact_df["w_heel"].to_numpy(float)[:, None]
    w_forefoot = contact_df["w_forefoot"].to_numpy(float)[:, None]
    w_toe = contact_df["w_toe"].to_numpy(float)[:, None]
    is_stance = contact_df["is_stance"].to_numpy(bool)
    pct = contact_df["cycle_percent"].to_numpy(float)

    contact_point = w_heel * heel + w_forefoot * forefoot + w_toe * toe
    translation = np.zeros_like(contact_point)
    translation[is_stance, 1] = -contact_point[is_stance, 1]
    translation[is_stance, 2] = -contact_point[is_stance, 2]
    if stride_time_s:
        time_s = pct / 100.0 * stride_time_s
        desired_contact_x = -walking_speed_m_s * time_s
        translation[is_stance, 0] = desired_contact_x[is_stance] - contact_point[is_stance, 0]

    for axis in range(3):
        translation[:, axis] = periodic_interpolate(translation[:, axis], is_stance)

    swing = ~is_stance
    if np.any(swing):
        swing_floor = np.minimum.reduce([heel[swing, 2], forefoot[swing, 2], toe[swing, 2]])
        correction = np.maximum(0.0, min_swing_clearance_m - (swing_floor + translation[swing, 2]))
        translation[swing, 2] += correction

    corrected = {name: arr + translation for name, arr in arrays.items()}
    quality = contact_df.copy()
    quality["contact_x_m"] = contact_point[:, 0]
    quality["contact_y_m"] = contact_point[:, 1]
    quality["contact_z_m"] = contact_point[:, 2]
    quality["translation_x_m"] = translation[:, 0]
    quality["translation_y_m"] = translation[:, 1]
    quality["translation_z_m"] = translation[:, 2]
    return corrected, quality


def finalize_world_arrays(arrays: dict[str, np.ndarray], leg: np.ndarray, contact_df: pd.DataFrame, min_swing_clearance_m: float) -> dict[str, np.ndarray]:
    pelvis_path = arrays["pelvis_origin"].copy()
    pelvis_path[:, 2] += -leg[2] * 1.1
    smoothed = {
        "pelvis": circular_gaussian(pelvis_path, 20),
        "hip": circular_gaussian(arrays["hip"], 20),
        "knee": circular_gaussian(arrays["knee"], 20),
        "ankle": circular_gaussian(arrays["ankle"], 20),
        "forefoot": circular_gaussian(arrays["forefoot"], 20),
        "toe": circular_gaussian(arrays["toe"], 20),
        "heel": circular_gaussian(arrays["heel"], 20),
    }
    heel = smoothed["heel"]
    forefoot = smoothed["forefoot"]
    toe = smoothed["toe"]
    w_heel = contact_df["w_heel"].to_numpy(float)[:, None]
    w_forefoot = contact_df["w_forefoot"].to_numpy(float)[:, None]
    w_toe = contact_df["w_toe"].to_numpy(float)[:, None]
    is_stance = contact_df["is_stance"].to_numpy(bool)
    contact_point = w_heel * heel + w_forefoot * forefoot + w_toe * toe
    z_correction = np.zeros(len(contact_df), dtype=float)
    z_correction[is_stance] = -contact_point[is_stance, 2].ravel()
    swing = ~is_stance
    corrected_floor = np.minimum.reduce(
        [
            heel[:, 2] + z_correction,
            forefoot[:, 2] + z_correction,
            toe[:, 2] + z_correction,
        ]
    )
    z_correction += np.maximum(0.0, -corrected_floor)
    if np.any(swing):
        swing_floor = np.minimum.reduce(
            [
                heel[swing, 2] + z_correction[swing],
                forefoot[swing, 2] + z_correction[swing],
                toe[swing, 2] + z_correction[swing],
            ]
        )
        z_correction[swing] += np.maximum(0.0, min_swing_clearance_m - swing_floor)
    for arr in smoothed.values():
        arr[:, 2] += z_correction
    return smoothed


def arrays_to_spatial_df(
    arrays: dict[str, np.ndarray],
    pct: np.ndarray,
    source: str,
    d2w_centered: bool,
) -> pd.DataFrame:
    n_points = len(pct)
    shift = n_points // 2
    working = {name: arr.copy() for name, arr in arrays.items()}
    if d2w_centered:
        for arr in working.values():
            arr[:, 0] *= -1
            arr[:, 1] *= -1
        working["pelvis"][:, 1] *= 0.6
        ankle_x_offset = 0.5 * (working["ankle"][:, 0].max() + working["ankle"][:, 0].min())
    else:
        ankle_x_offset = 0.0

    rows = []
    for side, phase_shift, y_sign in [("R", 0, 1), ("L", shift, -1)]:
        for i, cycle_percent in enumerate(pct):
            j = (i + phase_shift) % n_points
            for joint, arr in working.items():
                point = arr[j].copy()
                if joint != "pelvis":
                    point[1] *= y_sign
                if d2w_centered and joint == "ankle":
                    point[0] -= ankle_x_offset
                rows.append(
                    {
                        "side": side,
                        "joint": joint,
                        "cycle_percent": cycle_percent,
                        "x_m": point[0],
                        "y_m": point[1],
                        "z_m": point[2],
                        "rotation_inputs": ";".join(
                            f"{axis}={profile},sign={KINEMATIC_ROTATION_SIGNS.get(profile, 1.0):g}"
                            for axis, profile in KINEMATIC_ROTATION_PROFILES.items()
                        ),
                        "segment_proportions": "D2W height-scaled vectors in d2w_kinematic_model_proportions.csv",
                        "source": source,
                    }
                )
    return pd.DataFrame(rows)


def reconstruct_spatial_layers(
    schwartz_profiles: pd.DataFrame,
    height_m: float,
    n_points: int,
    contact_cfg: ContactModelConfig,
    step_length_m: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    pct, raw_arrays, leg = compute_raw_bilateral_fk_points(schwartz_profiles, height_m, n_points)
    pelvis_origin_world, contact_quality = solve_overground_pelvis_translation(raw_arrays, pct, contact_cfg, step_length_m, height_m)
    world_arrays = apply_bilateral_world_translation(raw_arrays, pelvis_origin_world, leg)
    treadmill_arrays = project_world_to_pelvis_fixed_treadmill(world_arrays)
    d2w_arrays = convert_treadmill_to_d2w_convention(treadmill_arrays)
    world = bilateral_arrays_to_spatial_df(
        world_arrays,
        pct,
        "Schwartz mean angles transformed through bilateral FK and overground bilateral contact anchoring",
    )
    d2w = bilateral_arrays_to_spatial_df(
        d2w_arrays,
        pct,
        "Overground bilateral trajectory projected to pelvis-fixed treadmill frame and remapped to D2W signs",
    )
    return world, d2w, contact_quality


def build_spatial_reconstruction_quality(
    spatial: pd.DataFrame,
    spatiotemporal: pd.DataFrame,
    contact_quality: pd.DataFrame,
    toe_off_percent: float = 60.0,
) -> pd.DataFrame:
    rows = []
    target_step = spatiotemporal[spatiotemporal["converted_metric_name"].eq("step_length_m")]
    target_cadence = spatiotemporal[spatiotemporal["converted_metric_name"].eq("cadence_steps_min")]
    step_length_target = float(target_step["converted_mean"].iloc[0]) if not target_step.empty else np.nan
    cadence_target = float(target_cadence["converted_mean"].iloc[0]) if not target_cadence.empty else np.nan

    def metric(name: str, value: float, unit: str, description: str) -> None:
        rows.append(
            {
                "metric_name": name,
                "value": value,
                "unit": unit,
                "description": description,
            }
        )

    foot = spatial[spatial["joint"].isin(["forefoot", "toe", "heel"])]
    min_foot_z = float(foot["z_m"].min()) if not foot.empty else np.nan
    metric("min_foot_z_m", min_foot_z, "m", "Minimum final heel/forefoot/toe height after bilateral contact anchoring.")
    metric("max_ground_penetration_m", max(0.0, -min_foot_z), "m", "Positive diagnostic value if final heel/forefoot/toe penetrates below z=0.")

    foot_with_phase = foot.copy()
    foot_with_phase["source_phase_percent"] = np.where(
        foot_with_phase["side"].eq("L"),
        (foot_with_phase["cycle_percent"] + 50.0) % 100.0,
        foot_with_phase["cycle_percent"],
    )
    if not contact_quality.empty:
        contact_state = contact_quality[
            [
                "side",
                "cycle_percent",
                "lambda_support",
                "w_heel",
                "w_forefoot",
                "w_toe",
            ]
        ].copy()
        foot_with_phase = foot_with_phase.merge(contact_state, on=["side", "cycle_percent"], how="left")
        joint_weight_name = {
            "heel": "w_heel",
            "forefoot": "w_forefoot",
            "toe": "w_toe",
        }
        foot_with_phase["contact_point_weight"] = [
            row.get(joint_weight_name.get(row["joint"], ""), 0.0)
            for _, row in foot_with_phase.iterrows()
        ]
        foot_with_phase["is_active_contact_point"] = (
            foot_with_phase["lambda_support"].fillna(0.0).gt(0.0)
            & foot_with_phase["contact_point_weight"].fillna(0.0).gt(0.0)
        )
        non_active_foot = foot_with_phase[~foot_with_phase["is_active_contact_point"]]
        metric(
            "min_non_active_foot_z_m",
            float(non_active_foot["z_m"].min()) if not non_active_foot.empty else np.nan,
            "m",
            "Minimum heel/forefoot/toe Z for points not used as active stance contacts; separates non-active foot penetration from active contact anchoring.",
        )
    else:
        foot_with_phase["is_active_contact_point"] = False
        metric("min_non_active_foot_z_m", np.nan, "m", "Unavailable because contact quality table is empty.")
    swing = foot_with_phase[foot_with_phase["source_phase_percent"].gt(toe_off_percent)]
    metric("swing_foot_clearance_min_m", float(swing["z_m"].min()) if not swing.empty else np.nan, "m", "Diagnostic minimum heel/forefoot/toe height during side-specific swing phase; swing clearance is not imposed by moving the common pelvis.")

    stance = foot_with_phase[foot_with_phase["source_phase_percent"].le(toe_off_percent)]
    if not stance.empty:
        contact = stance.sort_values(["side", "joint", "cycle_percent"]).copy()
        contact[["dx", "dy", "dz"]] = contact.groupby(["side", "joint"])[["x_m", "y_m", "z_m"]].diff()
        speed_like = np.sqrt(contact["dx"] ** 2 + contact["dy"] ** 2 + contact["dz"] ** 2).dropna()
        stance_foot_velocity_rms = float(np.sqrt(np.mean(speed_like**2))) if len(speed_like) else np.nan
    else:
        stance_foot_velocity_rms = np.nan
    metric("stance_foot_displacement_rms_per_sample_m", stance_foot_velocity_rms, "m/sample", "RMS sample-to-sample heel/forefoot/toe displacement during nominal stance; lower values indicate a more fixed contact reference.")

    ankle = spatial[spatial["joint"].eq("ankle")]
    step_length_reconstructed = np.nan
    if not ankle.empty:
        per_side_range = ankle.groupby("side")["x_m"].agg(lambda s: float(s.max() - s.min()))
        if len(per_side_range):
            step_length_reconstructed = float(per_side_range.mean())
    metric("ankle_x_excursion_m", step_length_reconstructed, "m", "Mean ankle X excursion in the pelvis-fixed geometric D2W frame; not a true overground step length.")
    metric("step_length_target_m", step_length_target, "m", "Schwartz dimensionless step length converted using subject leg length.")
    metric("cadence_target_steps_min", cadence_target, "steps/min", "Schwartz dimensionless cadence converted using subject leg length.")

    pelvis = spatial[spatial["joint"].eq("pelvis")]
    pelvis_z_rom = float(pelvis["z_m"].max() - pelvis["z_m"].min()) if not pelvis.empty else np.nan
    pelvis_x_rom = float(pelvis["x_m"].max() - pelvis["x_m"].min()) if not pelvis.empty else np.nan
    metric("pelvis_x_rom_treadmill_m", pelvis_x_rom, "m", "Longitudinal pelvis range after projection to the pelvis-fixed D2W/treadmill frame.")
    metric("pelvis_vertical_rom_m", pelvis_z_rom, "m", "Vertical range of the exported D2W pelvis path.")

    if not contact_quality.empty:
        active = contact_quality[contact_quality["lambda_support"].gt(0.0)]
        err = active["active_contact_anchor_error_m"].dropna().to_numpy(float)
        metric("active_contact_anchor_error_rms_m", float(np.sqrt(np.mean(err**2))) if len(err) else np.nan, "m", "RMS active stance-contact anchor error from the bilateral pelvis solve.")
        active_contact_world_z = (
            active["pelvis_world_origin_z_m"].to_numpy(float)
            + active["contact_z_m"].to_numpy(float)
        )
        active_contact_world_z = active_contact_world_z[np.isfinite(active_contact_world_z)]
        metric(
            "active_contact_world_z_rms_m",
            float(np.sqrt(np.mean(active_contact_world_z**2))) if len(active_contact_world_z) else np.nan,
            "m",
            "RMS world Z of the active weighted contact point, computed as pelvis_world_origin_z_m + contact_z_m; isolates active ground anchoring from non-active foot penetration.",
        )
        for side in ["R", "L"]:
            side_err = active[active["side"].eq(side)]["active_contact_anchor_error_m"].dropna().to_numpy(float)
            metric(
                f"{side.lower()}_contact_anchor_error_rms_m",
                float(np.sqrt(np.mean(side_err**2))) if len(side_err) else np.nan,
                "m",
                f"RMS active contact-anchor residual for the {side} side.",
            )
        pelvis_origin_x = contact_quality.groupby("cycle_percent")["pelvis_world_origin_x_m"].first()
        pelvis_origin_y = contact_quality.groupby("cycle_percent")["pelvis_world_origin_y_m"].first()
        pelvis_origin_z = contact_quality.groupby("cycle_percent")["pelvis_world_origin_z_m"].first()
        metric(
            "pelvis_x_rom_world_m",
            float(pelvis_origin_x.max() - pelvis_origin_x.min()) if len(pelvis_origin_x) else np.nan,
            "m",
            "Longitudinal pelvis-origin range in the overground bilateral contact layer.",
        )
        metric(
            "pelvis_y_rom_world_m",
            float(pelvis_origin_y.max() - pelvis_origin_y.min()) if len(pelvis_origin_y) else np.nan,
            "m",
            "Mediolateral pelvis-origin range in the overground bilateral contact layer.",
        )
        ds = max(0.0, toe_off_percent - 50.0)
        z_by_phase = pelvis_origin_z.sort_index()
        z_jumps = []
        boundary_jumps = []
        origin_by_phase = (
            contact_quality.groupby("cycle_percent")[
                ["pelvis_world_origin_x_m", "pelvis_world_origin_y_m", "pelvis_world_origin_z_m"]
            ]
            .first()
            .sort_index()
        )
        for phase in [ds, 50.0, 50.0 + ds]:
            idx = int(np.argmin(np.abs(z_by_phase.index.to_numpy(float) - phase)))
            prev_idx = (idx - 1) % len(z_by_phase)
            z_jumps.append(abs(float(z_by_phase.iloc[idx] - z_by_phase.iloc[prev_idx])))
            if len(origin_by_phase):
                coords = origin_by_phase.to_numpy(float)
                idx = int(np.argmin(np.abs(origin_by_phase.index.to_numpy(float) - phase)))
                prev_idx = (idx - 1) % len(origin_by_phase)
                boundary_jumps.append(float(np.linalg.norm(coords[idx] - coords[prev_idx])))
        metric(
            "pelvis_z_jump_at_ds_m",
            float(max(z_jumps)) if z_jumps else np.nan,
            "m",
            "Maximum pelvis-origin Z sample jump around double-support transition boundaries.",
        )
        metric(
            "boundary_transition_jump_m",
            float(max(boundary_jumps)) if boundary_jumps else np.nan,
            "m",
            "Maximum pelvis-origin 3D sample jump at the configured support transition boundaries.",
        )
        anchors = contact_quality[contact_quality["lambda_support"].gt(0.5)]
        if not anchors.empty:
            right_heel_x = anchors[anchors["side"].eq("R")]["heel_anchor_x_m"]
            left_heel_x = anchors[anchors["side"].eq("L")]["heel_anchor_x_m"]
            if not right_heel_x.empty and not left_heel_x.empty:
                step_reconstructed = float(np.max(np.abs(left_heel_x.to_numpy(float) - float(right_heel_x.median()))))
            else:
                step_reconstructed = np.nan
        else:
            step_reconstructed = np.nan
        metric("step_length_reconstructed_m", step_reconstructed, "m", "Step length implied by the bilateral stance anchors.")

        candidate = contact_quality.pivot(index="cycle_percent", columns="side", values=["pelvis_candidate_x_m", "pelvis_candidate_y_m", "pelvis_candidate_z_m", "lambda_support"])
        if not candidate.empty and ("lambda_support", "R") in candidate.columns and ("lambda_support", "L") in candidate.columns:
            double_support = candidate[("lambda_support", "R")].gt(0.0) & candidate[("lambda_support", "L")].gt(0.0)
            for axis, col in [("x", "pelvis_candidate_x_m"), ("y", "pelvis_candidate_y_m"), ("z", "pelvis_candidate_z_m")]:
                if (col, "R") in candidate.columns and (col, "L") in candidate.columns:
                    gap = (candidate.loc[double_support, (col, "R")] - candidate.loc[double_support, (col, "L")]).dropna().to_numpy(float)
                    metric(
                        f"double_support_candidate_{axis}_gap_rms_m",
                        float(np.sqrt(np.mean(gap**2))) if len(gap) else np.nan,
                        "m",
                        f"RMS {axis.upper()} gap between right and left pelvis candidates during double support.",
                    )

    jump_errors = []
    for (_, _), sub in spatial.sort_values("cycle_percent").groupby(["side", "joint"]):
        coords = sub[["x_m", "y_m", "z_m"]].to_numpy(float)
        if len(coords) > 1:
            diffs = np.diff(np.vstack([coords, coords[:1]]), axis=0)
            jump_errors.append(float(np.max(np.linalg.norm(diffs, axis=1))))
    metric(
        "max_sample_jump_treadmill_m",
        float(max(jump_errors)) if jump_errors else np.nan,
        "m/sample",
        "Maximum circular sample-to-sample jump in the exported D2W/treadmill trajectory.",
    )

    symmetry_errors = []
    for joint in ["hip", "knee", "ankle", "forefoot", "toe", "heel"]:
        right = spatial[spatial["side"].eq("R") & spatial["joint"].eq(joint)].sort_values("cycle_percent")[["x_m", "y_m", "z_m"]].to_numpy(float)
        left = spatial[spatial["side"].eq("L") & spatial["joint"].eq(joint)].sort_values("cycle_percent")[["x_m", "y_m", "z_m"]].to_numpy(float)
        if len(right) and len(right) == len(left):
            left_compare = left.copy()
            left_compare[:, 1] *= -1.0
            shift = len(right) // 2
            left_compare = np.roll(left_compare, -shift, axis=0)
            symmetry_errors.append(np.sqrt(np.mean((right - left_compare) ** 2)))
    metric("left_right_symmetry_error_m", float(np.mean(symmetry_errors)) if symmetry_errors else np.nan, "m", "RMS difference after undoing mediolateral mirroring and 50% phase shift.")

    return pd.DataFrame(rows)


def write_table(df: pd.DataFrame, tables_dir: Path, name: str) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(tables_dir / name, index=False)


def save_fig(fig: plt.Figure, figures_dir: Path, stem: str) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_angular_profiles(df: pd.DataFrame, figures_dir: Path, patient_id: str) -> None:
    grid = [
        ["trunk_tilt", "trunk_obliquity", "trunk_rotation"],
        ["pelvis_tilt", "pelvis_obliquity", "pelvis_rotation"],
        ["hip_flex", "hip_adduction", "hip_rotation"],
        ["knee_flex", "knee_adduction", "knee_rotation"],
        ["ankle_dorsi", None, "foot_progression"],
    ]
    fig, axes = plt.subplots(5, 3, figsize=(10.8, 10.2), sharex=True)
    for row_axes, row_profiles in zip(axes, grid):
        for ax, profile in zip(row_axes, row_profiles):
            if profile is None:
                ax.axis("off")
                continue
            sub = df[df["profile_name"].eq(profile)]
            mean = sub[sub["stat"].eq("mean")].sort_values("cycle_percent")
            low = sub[sub["stat"].eq("-1sd")].sort_values("cycle_percent")
            high = sub[sub["stat"].eq("+1sd")].sort_values("cycle_percent")
            if mean.empty:
                ax.axis("off")
                continue
            ax.fill_between(mean["cycle_percent"], low["value_deg"], high["value_deg"], color=SUBJECT_BLUE, alpha=0.16, lw=0)
            ax.plot(mean["cycle_percent"], mean["value_deg"], color=SUBJECT_BLUE, lw=1.8)
            ax.set_title(PROFILE_LABELS[profile], fontsize=10)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Gait cycle (%)", fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("Angle (deg)", fontsize=9)
    fig.suptitle(f"{patient_id}: angular trajectories vs normalized gait cycle", fontsize=12, y=0.98)
    fig.subplots_adjust(top=0.93, hspace=0.58, wspace=0.32)
    save_fig(fig, figures_dir, "schwartz_speed_matched_angular_profiles")


def plot_spatial_views(df: pd.DataFrame, figures_dir: Path, patient_id: str) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    colors = {"pelvis": "#666666", "hip": "#0072B2", "knee": "#009E73", "ankle": "#D55E00", "forefoot": "#CC79A7", "toe": "#E69F00", "heel": "#000000"}
    views = [
        ("Sagittal", "x_m", "z_m", "X (m)", "Z (m)"),
        ("Frontal", "y_m", "z_m", "Y (m)", "Z (m)"),
        ("Transverse", "x_m", "y_m", "X (m)", "Y (m)"),
    ]
    side_style = {"L": "--", "R": "-"}
    for ax, (title, x_col, y_col, x_lab, y_lab) in zip(axes, views):
        for joint in ["pelvis", "hip", "knee", "ankle", "forefoot", "toe", "heel"]:
            sides = ["C"] if joint == "pelvis" else ["L", "R"]
            for side in sides:
                sub = df[df["side"].eq(side) & df["joint"].eq(joint)].sort_values("cycle_percent")
                ax.plot(
                    sub[x_col],
                    sub[y_col],
                    color=colors[joint],
                    lw=1.7,
                    ls=side_style.get(side, "-"),
                    alpha=0.95,
                    label=f"{joint} {side}",
                )
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(x_lab, fontsize=10)
        ax.set_ylabel(y_lab, fontsize=10)
        ax.axis("equal")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    keep = []
    seen = set()
    for h, lab in zip(handles, labels):
        joint = lab.split()[0]
        if joint not in seen:
            keep.append((h, joint))
            seen.add(joint)
    side_handles = [
        plt.Line2D([0], [0], color="0.2", lw=1.8, ls="-", label="R"),
        plt.Line2D([0], [0], color="0.2", lw=1.8, ls="--", label="L"),
    ]
    fig.legend(
        [h for h, _ in keep] + side_handles,
        [lab for _, lab in keep] + ["R side", "L side"],
        frameon=False,
        fontsize=8,
        loc="lower center",
        ncol=8,
    )
    fig.suptitle(f"{patient_id}: spatial trajectories in sagittal, frontal and transverse views", fontsize=12, y=0.98)
    fig.subplots_adjust(top=0.78, bottom=0.24, wspace=0.35)
    save_fig(fig, figures_dir, "spatial_trajectory_three_plane_views")


def plot_spatial_cycle_profiles(df: pd.DataFrame, figures_dir: Path, patient_id: str) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 6.2), sharex=True)
    colors = {"pelvis": "#666666", "hip": "#0072B2", "knee": "#009E73", "ankle": "#D55E00"}
    coord_specs = [("x_m", "X (m)"), ("y_m", "Y (m)"), ("z_m", "Z (m)")]
    for ax, (coord, ylabel) in zip(axes, coord_specs):
        for joint in ["pelvis", "hip", "knee", "ankle"]:
            side_specs = [("C", "-")] if joint == "pelvis" else [("R", "-"), ("L", "--")]
            for side, ls in side_specs:
                sub = df[df["side"].eq(side) & df["joint"].eq(joint)].sort_values("cycle_percent")
                ax.plot(sub["cycle_percent"], sub[coord], color=colors[joint], ls=ls, lw=1.7, label=f"{joint} {side}")
        ax.set_ylabel(ylabel, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].set_xlabel("Normalized gait cycle (%)", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, loc="lower center", ncol=5)
    fig.suptitle(f"{patient_id}: spatial trajectories vs normalized gait cycle", fontsize=12)
    fig.subplots_adjust(bottom=0.18, hspace=0.28)
    save_fig(fig, figures_dir, "spatial_trajectories_vs_normalized_gait_cycle")


def plot_spatial_segment_sequence(df: pd.DataFrame, figures_dir: Path, patient_id: str) -> None:
    side_df = df[df["side"].eq("R")].copy()
    pelvis_df = df[df["side"].eq("C") & df["joint"].eq("pelvis")].copy()
    if side_df.empty or pelvis_df.empty:
        return
    side_df = pd.concat([side_df, pelvis_df], ignore_index=True)
    views = [
        ("Sagittal", "x_m", "z_m", "X (m)", "Z (m)"),
        ("Frontal", "y_m", "z_m", "Y (m)", "Z (m)"),
        ("Transverse", "x_m", "y_m", "X (m)", "Y (m)"),
    ]
    links = [
        ("pelvis", "hip"),
        ("hip", "knee"),
        ("knee", "ankle"),
        ("ankle", "forefoot"),
        ("forefoot", "toe"),
        ("ankle", "heel"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    cmap = plt.cm.Greys
    norm = plt.Normalize(0, 100)
    for ax, (title, x_col, y_col, x_lab, y_lab) in zip(axes, views):
        line_segments = []
        colors = []
        for pct in sorted(side_df["cycle_percent"].unique()):
            phase = side_df[side_df["cycle_percent"].eq(pct)].set_index("joint")
            for a, b in links:
                if a not in phase.index or b not in phase.index:
                    continue
                line_segments.append(
                    [
                        [phase.loc[a, x_col], phase.loc[a, y_col]],
                        [phase.loc[b, x_col], phase.loc[b, y_col]],
                    ]
                )
                colors.append(pct)
        lc = LineCollection(line_segments, cmap=cmap, norm=norm, linewidths=1.0, alpha=0.85)
        lc.set_array(np.asarray(colors, dtype=float))
        ax.add_collection(lc)
        latest = side_df[side_df["cycle_percent"].eq(side_df["cycle_percent"].max())].set_index("joint")
        for a, b in links:
            if a in latest.index and b in latest.index:
                ax.plot([latest.loc[a, x_col], latest.loc[b, x_col]], [latest.loc[a, y_col], latest.loc[b, y_col]], color="0.05", lw=1.8)
        ax.autoscale()
        ax.axis("equal")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(x_lab, fontsize=10)
        ax.set_ylabel(y_lab, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.03)
    cbar.set_label("Gait cycle (%)", fontsize=9)
    fig.suptitle(f"{patient_id}: connected spatial segments across gait cycle", fontsize=12, y=0.98)
    fig.subplots_adjust(top=0.80, bottom=0.18, wspace=0.35)
    save_fig(fig, figures_dir, "spatial_connected_segments_gait_cycle_gradient")


def animate_spatial_gait_3d(
    df: pd.DataFrame,
    figures_dir: Path,
    patient_id: str,
    stride_time_s: float | None,
    animation_fps: int = ANIMATION_FPS,
) -> None:
    from matplotlib.animation import FuncAnimation, PillowWriter

    figures_dir.mkdir(parents=True, exist_ok=True)
    pelvis = df[df["side"].eq("C") & df["joint"].eq("pelvis")].sort_values("cycle_percent")
    if pelvis.empty:
        return
    cycles = pelvis["cycle_percent"].to_numpy(float)
    joint_paths: dict[str, dict[str, np.ndarray]] = {"C": {}, "R": {}, "L": {}}
    for side in ["C", "R", "L"]:
        side_df = df[df["side"].eq(side)]
        for joint in side_df["joint"].unique():
            sub = side_df[side_df["joint"].eq(joint)].sort_values("cycle_percent")
            if len(sub) == len(cycles):
                joint_paths[side][joint] = sub[["x_m", "y_m", "z_m"]].to_numpy(float)
    if "pelvis" not in joint_paths["C"] or not joint_paths["R"] or not joint_paths["L"]:
        return
    real_time_duration_s = stride_time_s if stride_time_s and np.isfinite(stride_time_s) and stride_time_s > 0.0 else len(cycles) / 10.0
    frame_count = max(12, int(round(real_time_duration_s * animation_fps)))
    sample_frames = np.floor(np.linspace(0, len(cycles), frame_count, endpoint=False)).astype(int)
    sample_frames = np.clip(sample_frames, 0, len(cycles) - 1)

    links = [
        ("pelvis", "hip"),
        ("hip", "knee"),
        ("knee", "ankle"),
        ("ankle", "forefoot"),
        ("forefoot", "toe"),
        ("ankle", "heel"),
    ]
    side_colors = {"R": "#0072B2", "L": "#D55E00"}
    all_points = df[["x_m", "y_m", "z_m"]].to_numpy(float)
    mins = all_points.min(axis=0)
    maxs = all_points.max(axis=0)
    span = np.maximum(maxs - mins, 0.05)
    pad = 0.12 * span
    mins -= pad
    maxs += pad

    fig = plt.figure(figsize=(9.5, 4.8))
    axes = [
        fig.add_subplot(1, 2, 1, projection="3d"),
        fig.add_subplot(1, 2, 2, projection="3d"),
    ]
    view_specs = [
        ("Lateral 3D view", 8, -90),
        ("Frontal 3D view", 8, 0),
    ]
    artists = []
    for ax, (title, elev, azim) in zip(axes, view_specs):
        ax.set_title(title, fontsize=11)
        ax.view_init(elev=elev, azim=azim)
        ax.set_proj_type("ortho")
        ax.set_xlim(mins[0], maxs[0])
        ax.set_ylim(mins[1], maxs[1])
        ax.set_zlim(mins[2], maxs[2])
        ax.set_xlabel("X (m)", fontsize=8)
        ax.set_ylabel("Y (m)", fontsize=8)
        ax.set_zlabel("Z (m)", fontsize=8)
        ax.set_box_aspect((span[0], span[1], span[2]))
        ax.grid(True, alpha=0.2)
        view_lines = []
        view_points = []
        for side in ["R", "L"]:
            for _ in links:
                (line,) = ax.plot([], [], [], color=side_colors[side], lw=2.2, alpha=0.95)
                view_lines.append((side, line))
            point = ax.scatter([], [], [], color=side_colors[side], s=18, depthshade=False)
            view_points.append((side, point))
        pelvis_point = ax.scatter([], [], [], color="#222222", s=26, depthshade=False)
        phase_text = ax.text2D(0.03, 0.94, "", transform=ax.transAxes, fontsize=9)
        artists.append((view_lines, view_points, pelvis_point, phase_text))

    def point_for(side: str, joint: str, frame: int) -> np.ndarray | None:
        if joint == "pelvis":
            return joint_paths["C"]["pelvis"][frame]
        if joint not in joint_paths[side]:
            return None
        return joint_paths[side][joint][frame]

    def update(frame_no: int):
        frame = int(sample_frames[frame_no])
        changed = []
        for view_lines, view_points, pelvis_point, phase_text in artists:
            line_i = 0
            for side in ["R", "L"]:
                side_joint_points = []
                for a, b in links:
                    pa = point_for(side, a, frame)
                    pb = point_for(side, b, frame)
                    _, line = view_lines[line_i]
                    line_i += 1
                    if pa is None or pb is None:
                        line.set_data_3d([], [], [])
                    else:
                        line.set_data_3d([pa[0], pb[0]], [pa[1], pb[1]], [pa[2], pb[2]])
                        changed.append(line)
                for joint in ["hip", "knee", "ankle", "forefoot", "toe", "heel"]:
                    p = point_for(side, joint, frame)
                    if p is not None:
                        side_joint_points.append(p)
                _, point = next(item for item in view_points if item[0] == side)
                if side_joint_points:
                    pts = np.vstack(side_joint_points)
                    point._offsets3d = (pts[:, 0], pts[:, 1], pts[:, 2])
                    changed.append(point)
            pp = joint_paths["C"]["pelvis"][frame]
            pelvis_point._offsets3d = ([pp[0]], [pp[1]], [pp[2]])
            phase_text.set_text(f"{cycles[frame]:.0f}% gait cycle | {real_time_duration_s:.2f} s/stride")
            changed.extend([pelvis_point, phase_text])
        return changed

    animation = FuncAnimation(fig, update, frames=frame_count, interval=1000.0 / animation_fps, blit=False)
    fig.suptitle(f"{patient_id}: 3D gait animation", fontsize=12)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.90, wspace=0.05)
    animation.save(figures_dir / "spatial_gait_3d_frontal_lateral_animation.gif", writer=PillowWriter(fps=animation_fps), dpi=120)
    plt.close(fig)


def write_run_readme(
    spec: PatientSpec,
    out_dir: Path,
    group: str,
    status: str,
    v_star: float,
    speed_plan: SchwartzSpeedInterpolationPlan,
) -> None:
    text = f"""# PediaGait3D output: {spec.patient_id}

Input:

- age_y: {spec.age_y}
- height_m: {spec.height_m}
- speed_m_s: {spec.speed_m_s}
- leg_length_m: {spec.leg_length_m if spec.leg_length_m is not None else 'estimated as 0.530 * height'}
- dimensionless_speed: {v_star:.6f}
- Schwartz spatiotemporal speed group: {group} ({status})
- Schwartz angular waveform speed interpolation: {speed_plan.group_label}
- Schwartz angular waveform interpolation mode: {speed_plan.mode}
- Schwartz angular waveform interpolation alpha: {speed_plan.alpha:.6f}
- Schwartz angular waveform interpolation status: {speed_plan.status}

Folders:

- `tables/`: CSV outputs.
- `figures/`: PNG summaries.
- `tables/spatial_reconstruction_quality.csv`: automatic geometric quality checks.
- `tables/contact_phase_model.csv`: bilateral contact weights, support sharing, anchors and residuals.
- `tables/world_spatial_joint_trajectories.csv`: overground bilateral-contact world-layer trajectory.
- `tables/d2w_style_spatial_joint_trajectories.csv`: pelvis-fixed D2W-sign geometric trajectory.
- `tables/d2w_ankle_command_trajectories.csv`: ankle-only centered X command table for robot workspace compatibility.
- `figures/spatial_gait_3d_frontal_lateral_animation.gif`: 3D gait animation with lateral and frontal views, timed to the Schwartz-derived stride time.

Method notes:

- Schwartz profiles are speed-matched waveform references, not age-matched targets.
- Schwartz angular profiles are generated from the requested dimensionless speed using the configured interpolation mode; `linear` interpolates between adjacent Schwartz speed-group mean curves and clamps outside the published mean range.
- Van Hamme features are age/speed discrete references for young children.
- D2W spatial trajectories are pelvis-fixed geometric reconstructions, not a BTS/Helen-Hayes clinical model.
- Spatial reconstruction is split into bilateral raw FK, bilateral support/contact anchoring, pelvis-fixed treadmill projection, and D2W-sign geometric remapping.
- Pelvis, hip and knee use available Schwartz three-plane rotations in the D2W `R = Rz * Ry * Rx` chain.
- Schwartz pelvis tilt, hip adduction, knee adduction and knee flexion are sign-converted for the D2W FK convention; hip and knee axial rotations keep their original signs.
- Right and left limbs are solved together with a common pelvis translation; the left limb uses a 50% phase shift while pelvis orientation remains at the global gait phase.
- Ground anchors use the raw-FK forward sign, side-specific mediolateral foot tracks, and a heel/forefoot/toe rocker anchor matched to the active contact weights.
- Swing clearance and foot penetration are reported as diagnostics; the common pelvis is not lifted to hide swing-foot errors because doing so would break stance contact.
- The full-body D2W CSV preserves segment geometry; ankle-only centering is separated into `d2w_ankle_command_trajectories.csv`.
- Ankle uses Schwartz dorsi/plantarflexion; foot progression is used as terminal foot/ankle yaw because independent ankle frontal/axial rotations are not provided in the Schwartz supplement.
- The height-scaled D2W segment vectors used for reconstruction are written to `tables/d2w_kinematic_model_proportions.csv`.
- Trunk curves are included in the angular table and figure, but are not used in the lower-limb D2W spatial chain.
"""
    (out_dir / "README_OUTPUT.md").write_text(text)


def run_one(
    spec: PatientSpec,
    out_dir: Path,
    n_points: int,
    contact_cfg: ContactModelConfig | None = None,
    options: GenerationOptions | None = None,
) -> None:
    contact_cfg = contact_cfg or ContactModelConfig()
    options = options or GenerationOptions()
    out_dir.mkdir(parents=True, exist_ok=True)
    tables_dir = out_dir / "tables"
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for old_svg in figures_dir.glob("*.svg"):
        old_svg.unlink()
    leg_length_m = spec.leg_length_m if spec.leg_length_m else 0.530 * spec.height_m
    v_star = dimensionless_speed(spec.speed_m_s, leg_length_m)
    group, status = schwartz_speed_group(v_star)
    speed_plan = schwartz_speed_interpolation_plan(v_star, options.schwartz_speed_interpolation)

    context = pd.DataFrame(
        [
            {
                "patient_id": spec.patient_id,
                "age_y": spec.age_y,
                "height_m": spec.height_m,
                "leg_length_m": leg_length_m,
                "walking_speed_m_s": spec.speed_m_s,
                "dimensionless_speed": v_star,
                "schwartz_speed_group": group,
                "schwartz_group_assignment_status": status,
                "schwartz_angular_speed_group": speed_plan.group_label,
                "schwartz_angular_speed_lower_group": speed_plan.lower_group,
                "schwartz_angular_speed_upper_group": speed_plan.upper_group,
                "schwartz_angular_speed_lower_v_star": speed_plan.lower_v_star,
                "schwartz_angular_speed_upper_v_star": speed_plan.upper_v_star,
                "schwartz_angular_speed_interpolation_alpha": speed_plan.alpha,
                "schwartz_angular_speed_interpolation_mode": speed_plan.mode,
                "schwartz_angular_speed_interpolation_status": speed_plan.status,
                "leg_length_note": "estimated as 0.530*height unless leg_length_m is provided",
            }
        ]
    )
    write_table(context, tables_dir, "patient_context.csv")

    spatiotemporal = build_spatiotemporal_reference(group, leg_length_m)
    spatiotemporal.insert(0, "patient_id", spec.patient_id)
    spatiotemporal.insert(1, "dimensionless_speed_input", v_star)
    write_table(spatiotemporal, tables_dir, "normative_spatiotemporal_reference.csv")

    vh = build_van_hamme_features(spec.age_y, v_star)
    vh.insert(0, "patient_id", spec.patient_id)
    vh.insert(1, "age_y", spec.age_y)
    vh.insert(2, "dimensionless_speed", v_star)
    write_table(vh, tables_dir, "van_hamme_age_speed_discrete_reference.csv")

    schwartz_profiles = build_schwartz_profiles(speed_plan, n_points)
    schwartz_profiles.insert(0, "patient_id", spec.patient_id)
    schwartz_profiles.insert(1, "dimensionless_speed_input", v_star)
    write_table(schwartz_profiles, tables_dir, "schwartz_speed_matched_angular_profiles.csv")

    proportions = build_d2w_kinematic_model_proportions(spec.height_m)
    proportions.insert(0, "patient_id", spec.patient_id)
    proportions.insert(1, "height_m", spec.height_m)
    write_table(proportions, tables_dir, "d2w_kinematic_model_proportions.csv")

    stride_time_s = stride_time_from_spatiotemporal(spatiotemporal)
    step_length_m = converted_spatiotemporal_value(spatiotemporal, "step_length_m")
    if not np.isfinite(step_length_m) and stride_time_s:
        step_length_m = spec.speed_m_s * stride_time_s / 2.0
    if not np.isfinite(step_length_m):
        step_length_m = 0.35 * spec.height_m
    world_spatial, spatial, contact_quality = reconstruct_spatial_layers(
        schwartz_profiles,
        spec.height_m,
        n_points,
        contact_cfg,
        step_length_m,
    )

    contact_quality.insert(0, "patient_id", spec.patient_id)
    write_table(contact_quality, tables_dir, "contact_phase_model.csv")

    world_spatial.insert(0, "patient_id", spec.patient_id)
    world_spatial.insert(1, "age_y", spec.age_y)
    world_spatial.insert(2, "height_m", spec.height_m)
    world_spatial.insert(3, "walking_speed_m_s", spec.speed_m_s)
    world_spatial.insert(4, "dimensionless_speed", v_star)
    world_spatial.insert(5, "schwartz_speed_group", group)
    world_spatial.insert(6, "schwartz_angular_speed_group", speed_plan.group_label)
    world_spatial.insert(7, "schwartz_angular_speed_interpolation_alpha", speed_plan.alpha)
    if options.export_world_spatial:
        write_table(world_spatial, tables_dir, "world_spatial_joint_trajectories.csv")

    spatial.insert(0, "patient_id", spec.patient_id)
    spatial.insert(1, "age_y", spec.age_y)
    spatial.insert(2, "height_m", spec.height_m)
    spatial.insert(3, "walking_speed_m_s", spec.speed_m_s)
    spatial.insert(4, "dimensionless_speed", v_star)
    spatial.insert(5, "schwartz_speed_group", group)
    spatial.insert(6, "schwartz_angular_speed_group", speed_plan.group_label)
    spatial.insert(7, "schwartz_angular_speed_interpolation_alpha", speed_plan.alpha)
    if options.export_d2w_centered:
        write_table(spatial, tables_dir, "d2w_style_spatial_joint_trajectories.csv")
        write_table(build_d2w_ankle_command_trajectories(spatial), tables_dir, "d2w_ankle_command_trajectories.csv")

    spatial_quality = build_spatial_reconstruction_quality(spatial, spatiotemporal, contact_quality, contact_cfg.toe_off_percent)
    spatial_quality.insert(0, "patient_id", spec.patient_id)
    write_table(spatial_quality, tables_dir, "spatial_reconstruction_quality.csv")

    plot_spatial_views(spatial, figures_dir, spec.patient_id)
    plot_spatial_cycle_profiles(spatial, figures_dir, spec.patient_id)
    plot_spatial_segment_sequence(spatial, figures_dir, spec.patient_id)
    plot_angular_profiles(schwartz_profiles, figures_dir, spec.patient_id)
    animate_spatial_gait_3d(spatial, figures_dir, spec.patient_id, stride_time_s)
    write_run_readme(spec, out_dir, group, status, v_star, speed_plan)

    print(f"patient_id={spec.patient_id}")
    print(f"dimensionless_speed={v_star:.4f} -> Schwartz group: {group} ({status})")
    print(
        "angular_waveforms="
        f"{speed_plan.group_label} mode={speed_plan.mode} alpha={speed_plan.alpha:.4f} "
        f"({speed_plan.status})"
    )
    print(f"wrote tables and figures to {out_dir}")


def specs_from_csv(path: Path) -> list[PatientSpec]:
    df = pd.read_csv(path)
    required = {"patient_id", "age_y", "height_m", "speed_m_s"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    specs = []
    for _, row in df.iterrows():
        leg_length = row.get("leg_length_m", np.nan)
        specs.append(
            PatientSpec(
                patient_id=str(row["patient_id"]),
                age_y=float(row["age_y"]),
                height_m=float(row["height_m"]),
                speed_m_s=float(row["speed_m_s"]),
                leg_length_m=None if pd.isna(leg_length) else float(leg_length),
            )
        )
    return specs


def contact_config_from_dict(data: dict | None) -> ContactModelConfig:
    data = data or {}
    return ContactModelConfig(
        enabled=bool(data.get("enabled", True)),
        pattern=str(data.get("pattern", "heel_strike")),
        stance_start_percent=float(data.get("stance_start_percent", 0.0)),
        toe_off_percent=float(data.get("toe_off_percent", 60.0)),
        heel_to_flat_percent=float(data.get("heel_to_flat_percent", 12.0)),
        forefoot_loading_percent=float(data.get("forefoot_loading_percent", 45.0)),
        min_swing_clearance_m=float(data.get("min_swing_clearance_m", 0.025)),
    )


def spec_from_config(path: Path) -> tuple[PatientSpec, Path, int, ContactModelConfig, GenerationOptions]:
    with path.open() as f:
        cfg = json.load(f)
    required = ["patient_id", "age_y", "height_m", "speed_m_s"]
    missing = [name for name in required if name not in cfg]
    if missing:
        raise ValueError(f"{path} is missing required config keys: {missing}")
    n_points = int(cfg.get("n_points", 100))
    if n_points < 20:
        raise ValueError("n_points must be at least 20.")
    leg_length = cfg.get("leg_length_m")
    spec = PatientSpec(
        patient_id=str(cfg["patient_id"]),
        age_y=float(cfg["age_y"]),
        height_m=float(cfg["height_m"]),
        speed_m_s=float(cfg["speed_m_s"]),
        leg_length_m=None if leg_length is None else float(leg_length),
    )
    out_dir = resolve_package_out_dir(Path(cfg["out_dir"]) if cfg.get("out_dir") else None, spec.patient_id)
    contact_cfg = contact_config_from_dict(cfg.get("contact_model"))
    options = GenerationOptions(
        export_world_spatial=bool(cfg.get("export_world_spatial", True)),
        export_d2w_centered=bool(cfg.get("export_d2w_centered", True)),
        schwartz_speed_interpolation=str(cfg.get("schwartz_speed_interpolation", "linear")),
    )
    return spec, out_dir, n_points, contact_cfg, options


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=None, help="JSON config file with subject parameters and output folder.")
    parser.add_argument("--patient-id", default=None)
    parser.add_argument("--age-y", type=float, default=None)
    parser.add_argument("--height-m", type=float, default=None)
    parser.add_argument("--speed-m-s", type=float, default=None)
    parser.add_argument("--leg-length-m", type=float, default=None)
    parser.add_argument("--n-points", type=int, default=100)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--examples-csv", type=Path, default=None, help="Batch CSV with patient_id, age_y, height_m, speed_m_s, optional leg_length_m.")
    parser.add_argument(
        "--schwartz-speed-interpolation",
        choices=["linear", "nearest"],
        default="linear",
        help="Angular waveform speed handling: linear interpolation between Schwartz speed-group means, or nearest discrete group.",
    )
    args = parser.parse_args()

    if args.config:
        spec, out_dir, n_points, contact_cfg, options = spec_from_config(args.config)
        run_one(spec, out_dir, n_points, contact_cfg, options)
        return

    if args.examples_csv:
        batch_root = resolve_package_out_dir(args.out_dir or Path(OUTPUT_ROOT_NAME) / "batch_outputs", "batch_outputs")
        specs = specs_from_csv(args.examples_csv)
        batch_options = GenerationOptions(schwartz_speed_interpolation=args.schwartz_speed_interpolation)
        for spec in specs:
            run_one(spec, batch_root / spec.patient_id, args.n_points, ContactModelConfig(), batch_options)
        return

    missing = [name for name, value in [("patient-id", args.patient_id), ("age-y", args.age_y), ("height-m", args.height_m), ("speed-m-s", args.speed_m_s)] if value is None]
    if missing:
        raise SystemExit(f"Missing required arguments. Use --config or provide: {', '.join(missing)}")
    spec = PatientSpec(
        patient_id=args.patient_id,
        age_y=args.age_y,
        height_m=args.height_m,
        speed_m_s=args.speed_m_s,
        leg_length_m=args.leg_length_m,
    )
    run_one(
        spec,
        resolve_package_out_dir(args.out_dir, spec.patient_id),
        args.n_points,
        ContactModelConfig(),
        GenerationOptions(schwartz_speed_interpolation=args.schwartz_speed_interpolation),
    )


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
