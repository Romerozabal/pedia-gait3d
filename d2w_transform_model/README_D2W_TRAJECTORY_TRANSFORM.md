# D2W Robot Trajectory Transform Reference

This folder captures the D2W robot-side implementation that converts normalized
joint-angle gait references into spatial Cartesian trajectories. It is copied
from the ROS2 package at:

`/home/biorobotics/ros2_ws/src/d2w`

The source repository was not modified.

## Files Included

- `src/gait_generator.cpp`
  - Main implementation of angle loading, forward kinematics, pelvis/world-frame
    reconstruction, Cartesian remapping, smoothing, and interpolation.
- `include/gait_generator.hpp`
  - Function declarations and interpolator state used by the gait generator.
- `reference_angles/P_XYZ.txt`
- `reference_angles/H_XYZ.txt`
- `reference_angles/K_XYZ.txt`
- `reference_angles/A_XYZ.txt`
  - Pelvis, hip, knee, and ankle reference rotations loaded by the generator.
- `srv/GaitTrajectoryParameters.srv`
  - Related trajectory-parameter service definition with X/Z bounds, velocity,
    and patient height fields.
- `web/graficasEnDirecto.js`
  - Web plotting code that receives synchronized reference setpoints and displays
    spatial trajectories by gait-cycle phase.
- `metrics_context/metrics_node.py`
  - Metrics-node context showing that downstream robot metrics use published
    reference ankle trajectory points as per-cycle geometry.

## Active Transform Chain

The active transform is not exposed as one single `to_spatial()` function. It is
implemented as a sequence in `gait_generator.cpp`:

1. `loadAnglesFromFiles()` loads normalized 100-point reference rotations from
   `P_XYZ.txt`, `H_XYZ.txt`, `K_XYZ.txt`, and `A_XYZ.txt`.
2. `updateBodySegments(patient_height)` scales segment vectors from fixed body
   proportions.
3. `updateGaitGainsXYZ(...)` optionally rescales sagittal reference rotations
   and calls `applySagittalScaling()`.
4. `updateDirectGaitKinematics()` applies chained rotation matrices to compute
   pelvis-frame hip, knee, ankle, heel, and forefoot Cartesian positions.
5. `updateAnklePelvisPathsXYZ()` converts those pelvis-frame trajectories into
   a world/external-frame trajectory, constraining pelvis translation so that the
   pelvis does not carry longitudinal treadmill progression.
6. `applyCartesianScaling()` remaps ankle X and Z trajectory shape to the
   requested clinical trajectory bounds when amplitude control is active.
7. `updateSpatioTemporalGaitData()` derives step length and timing from the
   final ankle path.
8. `updatePaths()` creates gait-cycle interpolators for position, velocity, and
   acceleration.
9. `updateGaitPercentage()` samples those interpolators at the current right and
   left gait phases and publishes `GaitSetpoints`.

## Key Source Locations

In `src/gait_generator.cpp`:

- Lines 427-460: direct kinematics from rotations to joint/foot Cartesian
  positions.
- Lines 463-524: spatiotemporal timing from the final ankle path.
- Lines 706-818: sagittal rotation scaling with `parametrizeAndRemap`.
- Lines 821-900: Cartesian X/Z remapping of the ankle path.
- Lines 902-1001: generic affine shape remapping preserving the original
  normalized waveform shape.
- Lines 1003-1017: angle-reference file loading.
- Lines 1019-1030: high-level parameter update sequence.
- Lines 1065-1132: phase sampling and X-centering of the ankle path during
  publication.
- Lines 1158-1167: segment scaling from patient height.
- Lines 1169-1256: pelvis/world-frame trajectory reconstruction.
- Lines 1415-1453: circular Gaussian smoothing, preserving cycle continuity.
- Lines 1457-1599: output paths, derivatives, smoothing, and interpolators.

In `web/graficasEnDirecto.js`:

- Lines 823-914 receive synchronized reference setpoints and write them into
  gait-cycle phase buffers.

In `metrics_context/metrics_node.py`:

- Lines 186-203 show that the metrics node stores published reference ankle
  phase and position samples as per-cycle reference geometry.

## Pelvis-Frame To External-Frame Correction

The important correction is in `updateAnklePelvisPathsXYZ()`:

1. Direct kinematics first produces `ankle_absolute_xyz_` relative to the pelvis
   model.
2. A pelvis translation vector is built from the heel/ankle trajectory:
   `pelvis_trans.y = -heel1.y` and `pelvis_trans.z = -heel1.z`.
3. The lateral component is centered and symmetrized across the two half-cycles.
4. This translation is added to pelvis, hip, knee, ankle, forefoot, and heel
   points to shift the kinematic chain into an external/world-like frame.
5. The X and Y signs are corrected for the robot coordinate convention.
6. The final trajectories are smoothed with circular boundary handling.

For the longitudinal axis, pelvis translation is not made to follow treadmill
progression. In the current implementation the pelvis translation vector has no
active X component, and the published ankle X coordinate is centered by:

`x_offset = (max(ankle_path_x) + min(ankle_path_x)) / 2`

then:

`ankle_ref.x = interp_ankle_x(phase) - x_offset`

This is the part that makes the ankle trajectory interpretable in an external
cycle frame while keeping the pelvis longitudinal reference effectively static.

## Cartesian Scaling Convention

When amplitude control is active, `applyCartesianScaling()` extracts the current
ankle X and Z trajectories and remaps them to:

- `x_new_min = -hip_amplitude / 2`
- `x_new_max =  hip_amplitude / 2`
- `z_new_min =  knee_amplitude`
- `z_new_max =  ankle_amplitude`

The remapping preserves the original normalized waveform morphology:

`f = (source - source_min) / (source_max - source_min)`

`target = target_min + f * (target_max - target_min)`

So this is not a new fitted gait curve; it is an affine amplitude remapping of
the existing normalized reference trajectory.

## Practical Interpretation For LAM/Trajectory Work

If this method is reused outside ROS, the minimum ingredients are:

1. 100-point pelvis/hip/knee/ankle rotation profiles.
2. Patient height or segment lengths.
3. Segment vectors equivalent to the D2W body proportions.
4. Rotation convention `R = Rz * Ry * Rx`.
5. Pelvis/world-frame correction from `updateAnklePelvisPathsXYZ()`.
6. Circular smoothing and centered ankle X.

The method is appropriate for reconstructing the D2W reference geometry. It
should not be presented as a BTS/Helen Hayes clinical reconstruction method,
because its segment definitions, pelvis correction, and trajectory scaling are
robot-controller conventions.
