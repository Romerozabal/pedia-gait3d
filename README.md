# PediaGait3D

Normative pediatric gait trajectory reconstruction from subject parameters,
published reference models, and a D2W-compatible kinematic transform.

This repository contains a standalone Python workflow that estimates
spatiotemporal reference values, discrete angular reference parameters,
speed-matched angular waveforms, and approximate 3D spatial gait trajectories
from a small subject configuration file.

The package is intended as a reproducible reference generator. It is not a
validated clinical inverse-kinematics pipeline, not a BTS/Helen-Hayes
implementation, and not a force/contact dynamics model.

## Repository Contents

```text
.
├── generate_integrated_normative_gait.py
├── configs/
│   └── subject_gait_config.example.json
├── models/
│   ├── lam_schwartz_2008_anova_reference_long.csv
│   ├── lam_schwartz_2008_joint_rotation_reference_long.csv
│   └── lam_van_hamme_2015_regression_coefficients.csv
├── outputs/
│   └── .gitkeep
├── PACKAGE_MANIFEST.csv
└── requirements.txt
```

Generated outputs are written under `outputs/<patient_id>/` and are ignored by
Git. The repository keeps `outputs/.gitkeep` only so that the output root exists
after cloning.

## Installation

The code is intentionally small and uses standard scientific Python packages.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Current `requirements.txt`:

```text
numpy
pandas
matplotlib
```

## Basic Use

Edit:

```text
configs/subject_gait_config.example.json
```

Example:

```json
{
  "patient_id": "subject_001",
  "age_y": 9.0,
  "height_m": 1.25,
  "speed_m_s": 0.2,
  "leg_length_m": null,
  "n_points": 100,
  "schwartz_speed_interpolation": "linear",
  "out_dir": "outputs/subject_001",
  "contact_model": {
    "enabled": true,
    "pattern": "heel_strike",
    "stance_start_percent": 0.0,
    "toe_off_percent": 60.0,
    "heel_to_flat_percent": 12.0,
    "forefoot_loading_percent": 45.0,
    "min_swing_clearance_m": 0.025
  },
  "export_world_spatial": true,
  "export_d2w_centered": true
}
```

Run:

```bash
python3 generate_integrated_normative_gait.py \
  --config configs/subject_gait_config.example.json
```

The output folder is resolved relative to the repository root. Absolute output
paths outside the repository are rejected to keep generated files contained.

## Inputs

The required subject parameters are:

```text
patient_id
age_y
height_m
speed_m_s
```

Optional:

```text
leg_length_m
n_points
schwartz_speed_interpolation
contact_model
out_dir
```

If `leg_length_m` is not provided, the script estimates:

```text
leg_length_m = 0.530 * height_m
```

Dimensionless speed is computed as:

```text
v* = speed_m_s / sqrt(g * leg_length_m)
g  = 9.81 m/s^2
```

## Output Files

Each execution writes:

```text
outputs/<patient_id>/
├── README_OUTPUT.md
├── tables/
│   ├── patient_context.csv
│   ├── normative_spatiotemporal_reference.csv
│   ├── van_hamme_age_speed_discrete_reference.csv
│   ├── schwartz_speed_matched_angular_profiles.csv
│   ├── d2w_kinematic_model_proportions.csv
│   ├── contact_phase_model.csv
│   ├── world_spatial_joint_trajectories.csv
│   ├── d2w_style_spatial_joint_trajectories.csv
│   ├── d2w_ankle_command_trajectories.csv
│   └── spatial_reconstruction_quality.csv
└── figures/
    ├── spatial_trajectory_three_plane_views.png
    ├── spatial_trajectories_vs_normalized_gait_cycle.png
    ├── spatial_connected_segments_gait_cycle_gradient.png
    ├── schwartz_speed_matched_angular_profiles.png
    └── spatial_gait_3d_frontal_lateral_animation.gif
```

The full-body spatial tables include pelvis, hip, knee, ankle, forefoot, toe,
and heel trajectories. The ankle-only command table is exported separately for
D2W/robot workspace compatibility and should not be interpreted as a complete
body skeleton.

## Source Models And Local Tables

The repository distributes only the local CSV model tables consumed by the
generator. Original article PDFs and supplementary files are not included.
Users should obtain the cited source documents from the publishers when needed.

### Schwartz et al. 2008

Schwartz et al. 2008 provides the speed-binned typically developing pediatric
reference used here for:

- full-cycle angular waveforms;
- dimensionless step length;
- dimensionless cadence;
- dimensionless speed bins.

Local model tables:

```text
models/lam_schwartz_2008_joint_rotation_reference_long.csv
models/lam_schwartz_2008_anova_reference_long.csv
```

The angular model table was derived from the Schwartz supplementary joint
rotation data. It includes:

```text
trunk_tilt
trunk_obliquity
trunk_rotation
pelvis_tilt
pelvis_obliquity
pelvis_rotation
hip_flex
hip_adduction
hip_rotation
knee_flex
knee_adduction
knee_rotation
ankle_dorsi
foot_progression
```

For each profile, the table stores:

```text
schwartz_speed_group
stat = mean, -1sd, +1sd
cycle_percent
value_deg
source metadata
```

The Schwartz supplement does not provide independent ankle frontal-plane or
ankle axial-rotation curves. This implementation therefore uses ankle
dorsi/plantarflexion as ankle pitch and foot progression as the terminal
foot/ankle yaw approximation.

### Van Hamme et al. 2015

Van Hamme et al. 2015 is used for discrete age/speed angular reference
features in young children. The local model table stores regression
coefficients in:

```text
models/lam_van_hamme_2015_regression_coefficients.csv
```

The implemented regression form is:

```text
Y = a * age
  + b * v*
  + c * age * v*
  + d
```

where `v*` is dimensionless walking speed.

These values are exported as discrete reference parameters. They do not
currently warp or rescale the Schwartz full-cycle angular waveforms.

### D2W-Compatible Transform Convention

The former D2W reference snapshot is not required at runtime and is not
distributed in this public repository. The relevant transform assumptions are
implemented directly in `generate_integrated_normative_gait.py`: height-scaled
segment vectors, the `R = Rz * Ry * Rx` rotation convention, pelvis-fixed
treadmill projection, D2W Cartesian sign remapping, and the separate ankle-only
command export.

## Angular Waveform Construction

The script first computes the subject dimensionless speed:

```text
v* = speed_m_s / sqrt(g * leg_length_m)
```

The Schwartz speed bins are:

```text
Very Slow: [0.090, 0.226]
Slow:      [0.234, 0.359]
Free:      [0.372, 0.482]
Fast:      [0.502, 0.636]
Very Fast: [0.638, 0.844]
```

The spatiotemporal values are selected from the closest applicable Schwartz
speed group.

For angular waveforms, `schwartz_speed_interpolation` controls whether the
script uses a single speed group or interpolates between groups:

```text
nearest: use the closest Schwartz group waveform
linear:  interpolate between adjacent Schwartz group mean-speed waveforms
```

The temporal interpolation to the requested gait-cycle grid is:

```text
t_i = 100 * i / N,  i = 0, 1, ..., N-1
theta_g(t_i) = interp(theta_g_source, t_i)
```

For `linear` speed interpolation:

```text
alpha = (v* - v_low) / (v_high - v_low)

theta(t_i, v*) =
    (1 - alpha) * theta_low(t_i)
  + alpha       * theta_high(t_i)
```

where `v_low` and `v_high` are adjacent Schwartz group mean dimensionless
speeds. Outside the published mean-speed range, the angular waveform is clamped
to `Very Slow` or `Very Fast`; no extrapolation is performed.

The spatial reconstruction uses only the `mean` angular curves. The `-1sd` and
`+1sd` curves are exported and plotted for reference.

## Kinematic Convention

The transform uses the D2W rotation order:

```text
R = Rz(yaw) * Ry(pitch) * Rx(roll)
```

The angular profile mapping used by the forward kinematic chain is:

```text
pelvis_roll  Rx = + pelvis_obliquity
pelvis_pitch Ry = - pelvis_tilt
pelvis_yaw   Rz = + pelvis_rotation

hip_roll     Rx = - hip_adduction
hip_pitch    Ry = + hip_flex
hip_yaw      Rz = + hip_rotation

knee_roll    Rx = - knee_adduction
knee_pitch   Ry = - knee_flex
knee_yaw     Rz = + knee_rotation

ankle_roll   Rx = 0
ankle_pitch  Ry = + ankle_dorsi
ankle_yaw    Rz = + foot_progression
```

The sign conversions are explicit in `KINEMATIC_ROTATION_SIGNS` inside
`generate_integrated_normative_gait.py`. They represent the current mapping
between Schwartz clinical angle conventions and the D2W rotation axes used by
this FK implementation.

Trunk curves are exported in the angular table and plotted, but they are not
used in the lower-limb spatial reconstruction.

## Segment Geometry

The spatial model uses height-scaled segment vectors derived from the D2W
reference geometry:

```text
segment_vector_m = height_m * segment_proportion
```

The current proportional vectors are:

```text
leg_origin_offset = (0.0, 0.0, -0.530)
pelvis_to_hip     = (0.0, 0.191 / 2, -0.01)
thigh             = (0.0, 0.0, -0.246)
shank             = (0.0, 0.0, -0.246)
forefoot          = (-0.122 / 1.5, 0.0, -0.05)
toe               = (-0.122, 0.0, -0.045)
heel              = (0.03 / 1.5, 0.0, -0.05)
```

The exact scaled vectors used in a run are exported to:

```text
tables/d2w_kinematic_model_proportions.csv
```

## Spatial Reconstruction Workflow

The spatial reconstruction is intentionally separated into layers:

```text
Schwartz angular waveforms
  -> bilateral raw FK
  -> bilateral support/contact anchoring
  -> overground common pelvis trajectory
  -> pelvis-fixed treadmill projection
  -> D2W sign remapping
  -> ankle-only command export
```

### 1. Bilateral Forward Kinematics

The pelvis orientation is evaluated at the global gait-cycle phase. The right
limb uses that phase directly:

```text
phi_R = phi
```

The left limb uses a contralateral 50% shift:

```text
phi_L = (phi + 50) mod 100
```

The left limb is mirrored in the pelvis-local mediolateral axis before applying
the global pelvis orientation.

For each side, the FK chain is:

```text
pelvis
  -> hip
  -> knee
  -> ankle
  -> forefoot, toe, heel
```

with rotations:

```text
R_pelvis = Rz(pelvis_yaw) * Ry(pelvis_pitch) * Rx(pelvis_roll)
R_hip    = Rz(hip_yaw)    * Ry(hip_pitch)    * Rx(hip_roll)
R_knee   = Rz(knee_yaw)   * Ry(knee_pitch)   * Rx(knee_roll)
R_ankle  = Rz(ankle_yaw)  * Ry(ankle_pitch)  * Rx(ankle_roll)
```

### 2. Contact Weights

The default contact model is phase-kinematic and uses heel, forefoot, and toe
weights during stance:

```text
c_side(phi) =
    w_heel(phi)     * heel_side(phi)
  + w_forefoot(phi) * forefoot_side(phi)
  + w_toe(phi)      * toe_side(phi)
```

For the default `heel_strike` pattern:

```text
0-12%    heel-to-flat transition
12-45%   flat-foot support
45-60%   forefoot/toe transition
>60%     swing
```

Other configured patterns are `midfoot_strike` and `forefoot_strike`.

### 3. Bilateral Support Sharing

The right and left supports are combined with phase-dependent support weights:

```text
lambda_R(phi), lambda_L(phi)
```

The default toe-off is 60%, so there is a 10% double-support interval around
each limb transition.

### 4. Overground Common Pelvis

For each side, the model computes the pelvis translation candidate needed to
place the active contact point on its ground anchor:

```text
P_R_candidate(phi) = anchor_R(phi) - c_R(phi)
P_L_candidate(phi) = anchor_L(phi) - c_L(phi)
```

The common pelvis origin is the support-weighted bilateral solution:

```text
P_pelvis(phi) =
    [lambda_R(phi) * P_R_candidate(phi)
   + lambda_L(phi) * P_L_candidate(phi)]
  / [lambda_R(phi) + lambda_L(phi)]
```

The resulting overground trajectories are exported to:

```text
tables/world_spatial_joint_trajectories.csv
```

### 5. Pelvis-Fixed Treadmill Projection

The D2W geometric layer removes the overground longitudinal pelvis progression:

```text
x_treadmill_joint(phi) =
    x_world_joint(phi) - x_world_pelvis(phi)

x_treadmill_pelvis(phi) = 0
```

The vertical and mediolateral pelvis motion are retained from the overground
solution.

### 6. D2W Cartesian Remapping

The Cartesian sign remapping is:

```text
x_D2W = -x_treadmill
y_D2W = -y_treadmill
```

For the pelvis lateral coordinate, the inherited D2W scaling is kept:

```text
y_D2W_pelvis = -0.6 * y_treadmill_pelvis
```

The full-body D2W-style trajectory is exported to:

```text
tables/d2w_style_spatial_joint_trajectories.csv
```

The ankle-only centered command table is exported separately:

```text
tables/d2w_ankle_command_trajectories.csv
```

## Quality Metrics

The file:

```text
tables/spatial_reconstruction_quality.csv
```

contains diagnostics for geometry, contact anchoring, and continuity. These are
not clinical validation metrics.

Key contact metrics:

```text
active_contact_anchor_error_rms_m
active_contact_world_z_rms_m
min_non_active_foot_z_m
max_ground_penetration_m
double_support_candidate_z_gap_rms_m
boundary_transition_jump_m
```

The active contact height is computed as:

```text
z_active_contact(phi) =
    pelvis_world_origin_z(phi) + contact_z(phi)
```

`active_contact_world_z_rms_m` is therefore the preferred ground-height
diagnostic for the currently weighted contact point.

`min_non_active_foot_z_m` reports the lowest heel, forefoot, or toe point that
is not part of the active stance contact at that sample. This separates active
contact anchoring from penetration caused by simplified foot geometry.

`stance_foot_displacement_rms_per_sample_m` is retained as a legacy/global
stance-foot motion diagnostic. It is not the preferred rocker-contact quality
metric because heel, forefoot, and toe are not all expected to remain fixed
through the whole stance phase.

## Methodological Scope

This repository distinguishes the model components as follows:

- Schwartz et al. 2008 supplies speed-matched waveform and spatiotemporal
  reference context.
- Van Hamme et al. 2015 supplies discrete age/speed angular reference values.
- The D2W transform supplies a robot-compatible geometric reconstruction
  convention.
- The implemented contact model is kinematic and phase-defined; it does not
  estimate contact force, pressure, center of pressure, or unilateral dynamics.
- The resulting trajectories are generated references, not measured subject
  kinematics.

## Batch Mode

Batch generation is supported with a CSV file containing:

```csv
patient_id,age_y,height_m,speed_m_s,leg_length_m
S01,5,1.16,1.0,
S02,4,1.05,0.8,
```

Run:

```bash
python3 generate_integrated_normative_gait.py \
  --examples-csv path/to/subjects.csv \
  --out-dir outputs/batch_outputs
```

Each subject is written into its own subfolder.

## Manifest Check

The repository includes `PACKAGE_MANIFEST.csv` with SHA-256 checksums for
versioned package files. To verify it:

```bash
python3 - <<'PY'
import csv
import hashlib
from pathlib import Path

errors = []
root = Path(".")

with open("PACKAGE_MANIFEST.csv", newline="") as f:
    for row in csv.DictReader(f):
        path = root / row["package_path"]
        if not path.exists():
            errors.append((row["package_path"], "missing"))
            continue
        data = path.read_bytes()
        sha = hashlib.sha256(data).hexdigest()
        size = str(len(data))
        if sha != row["sha256"] or size != row["size_bytes"]:
            errors.append((row["package_path"], sha, size))

if errors:
    print("PACKAGE_MANIFEST.csv: FAIL")
    for item in errors:
        print(item)
    raise SystemExit(1)

print("PACKAGE_MANIFEST.csv: OK")
PY
```

## Minimal Validation Run

```bash
python3 -m py_compile generate_integrated_normative_gait.py
python3 generate_integrated_normative_gait.py \
  --config configs/subject_gait_config.example.json
```

The generated run should create CSV tables, PNG figures, and one GIF animation
inside `outputs/<patient_id>/`.
