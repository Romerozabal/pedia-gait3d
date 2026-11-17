#!/usr/bin/env python3
"""Regenerate the animations and figures embedded in the documentation.

Run from the repository root:

    python3 docs/make_readme_media.py

Everything is computed with the functions of
`generate_integrated_normative_gait.py`; no trajectory is drawn by hand.
Outputs are written to `docs/media/`.
"""

from __future__ import annotations

import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image, ImageSequence
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import generate_integrated_normative_gait as g  # noqa: E402

MEDIA_DIR = REPO_ROOT / "docs" / "media"
EXAMPLE_CONFIG = REPO_ROOT / "configs" / "subject_gait_config.comfortable_speed.example.json"

RIGHT = "#0072B2"
LEFT = "#D55E00"
INK = "#222222"
MUTED = "#6B6B6B"
GRID = "#DDDDDD"
GROUND = "#9A9A9A"
SIDE_COLORS = {"R": RIGHT, "L": LEFT}
SIDE_NAMES = {"R": "Right", "L": "Left"}
LIMB_CHAIN = ["hip", "knee", "ankle"]
FOOT_CHAIN = ["heel", "ankle", "forefoot", "toe"]

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
    }
)


def reconstruct(spec: g.PatientSpec, n_points: int = 100, mode: str = "linear") -> dict:
    """Run the generator pipeline in memory for one subject."""
    leg_length_m = spec.leg_length_m if spec.leg_length_m else 0.530 * spec.height_m
    v_star = g.dimensionless_speed(spec.speed_m_s, leg_length_m)
    group, _ = g.schwartz_speed_group(v_star)
    plan = g.schwartz_speed_interpolation_plan(v_star, mode)
    profiles = g.build_schwartz_profiles(plan, n_points)
    spatiotemporal = g.build_spatiotemporal_reference(group, leg_length_m)
    step_length_m = g.converted_spatiotemporal_value(spatiotemporal, "step_length_m")
    cadence = g.converted_spatiotemporal_value(spatiotemporal, "cadence_steps_min")
    stride_time_s = g.stride_time_from_spatiotemporal(spatiotemporal)
    contact_cfg = g.ContactModelConfig()
    world, d2w, contact = g.reconstruct_spatial_layers(profiles, spec.height_m, n_points, contact_cfg, step_length_m)
    return {
        "spec": spec,
        "v_star": v_star,
        "group": group,
        "plan": plan,
        "profiles": profiles,
        "step_length_m": step_length_m,
        "cadence_steps_min": cadence,
        "stride_time_s": stride_time_s,
        "world": world,
        "d2w": d2w,
        "contact": contact,
        "contact_cfg": contact_cfg,
    }


def joint_paths(spatial: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    paths = {}
    for (side, joint), sub in spatial.groupby(["side", "joint"]):
        paths[(side, joint)] = sub.sort_values("cycle_percent")[["x_m", "y_m", "z_m"]].to_numpy(float)
    return paths


def profile_curve(profiles: pd.DataFrame, name: str, stat: str) -> np.ndarray:
    sub = profiles[profiles["profile_name"].eq(name) & profiles["stat"].eq(stat)].sort_values("cycle_percent")
    return sub["value_deg"].to_numpy(float)


def shrink_gif(path: Path, fps: int, colors: int = 48) -> None:
    """Re-encode with one shared palette so only changed pixels are stored per frame."""
    with Image.open(path) as im:
        frames = [frame.convert("RGB") for frame in ImageSequence.Iterator(im)]
    palette = frames[len(frames) // 2].quantize(colors=colors, method=Image.Quantize.MEDIANCUT)
    quantized = [frame.quantize(palette=palette, dither=Image.Dither.NONE) for frame in frames]
    quantized[0].save(path, save_all=True, append_images=quantized[1:], duration=int(round(1000 / fps)), loop=0, optimize=True)


def save_gif(animation: FuncAnimation, name: str, fps: int, dpi: int) -> None:
    MEDIA_DIR.mkdir(parents=True, exist_ok=True)
    path = MEDIA_DIR / name
    animation.save(path, writer=PillowWriter(fps=fps), dpi=dpi)
    shrink_gif(path, fps)
    print(f"wrote {path.relative_to(REPO_ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


def make_pipeline_overview(run: dict) -> None:
    """Angles -> bilateral 3D skeleton -> ankle trajectory, one gait cycle."""
    spec = run["spec"]
    profiles = run["profiles"]
    paths = joint_paths(run["d2w"])
    command = joint_paths(g.build_d2w_ankle_command_trajectories(run["d2w"]))
    n = len(paths[("C", "pelvis")])
    pct = np.linspace(0, 100, n, endpoint=False)
    half = n // 2
    toe_off = run["contact_cfg"].toe_off_percent

    fig = plt.figure(figsize=(11.0, 5.2))
    grid = fig.add_gridspec(
        3, 3, width_ratios=[1.0, 1.15, 1.0], height_ratios=[1, 1, 1],
        left=0.065, right=0.985, top=0.83, bottom=0.17, wspace=0.22, hspace=0.42,
    )
    fig.suptitle(
        f"PediaGait3D: {spec.age_y:g}-year-old, {spec.height_m:.2f} m, {spec.speed_m_s:.2f} m/s "
        f"(dimensionless speed {run['v_star']:.2f})",
        fontsize=11, x=0.065, ha="left", y=0.975,
    )
    fig.text(0.065, 0.895, "1  Speed-matched joint angles", fontsize=10, fontweight="bold")
    fig.text(0.385, 0.895, "2  Bilateral 3D forward kinematics", fontsize=10, fontweight="bold")
    fig.text(0.715, 0.895, "3  Ankle trajectory for the robot", fontsize=10, fontweight="bold")

    # 1. Angular waveforms with right/left phase cursors.
    angle_specs = [("hip_flex", "Hip flexion"), ("knee_flex", "Knee flexion"), ("ankle_dorsi", "Ankle dorsiflexion")]
    angle_dots = []
    for row, (name, label) in enumerate(angle_specs):
        ax = fig.add_subplot(grid[row, 0])
        mean = profile_curve(profiles, name, "mean")
        low = profile_curve(profiles, name, "-1sd")
        high = profile_curve(profiles, name, "+1sd")
        ax.fill_between(pct, low, high, color=MUTED, alpha=0.18, lw=0)
        ax.plot(pct, mean, color=INK, lw=1.6)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_ylabel("deg", fontsize=8)
        ax.set_title(label, fontsize=9, loc="left", pad=3)
        ax.tick_params(labelsize=8)
        if row < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Gait cycle (%)", fontsize=8)
        (dot_r,) = ax.plot([], [], "o", color=RIGHT, ms=7, mec="white", mew=1.2, zorder=5)
        (dot_l,) = ax.plot([], [], "o", color=LEFT, ms=7, mec="white", mew=1.2, zorder=5)
        angle_dots.append((mean, dot_r, dot_l))
        if row == 0:
            handles = [
                plt.Line2D([], [], marker="o", ls="", color=RIGHT, ms=6, label="Right limb"),
                plt.Line2D([], [], marker="o", ls="", color=LEFT, ms=6, label="Left limb (same curve, 50% later)"),
                plt.Line2D([], [], color=INK, lw=1.6, label="Reference mean (band: ±1 SD)"),
            ]
            fig.legend(handles=handles, fontsize=8, frameon=False, ncol=3, loc="lower left",
                       bbox_to_anchor=(0.06, 0.0), handlelength=1.4, columnspacing=1.6, handletextpad=0.5)

    # 2. 3D skeleton in the pelvis-fixed frame.
    ax3d = fig.add_subplot(grid[0:2, 1], projection="3d")
    pts = np.vstack([arr for (side, _), arr in paths.items() if side in ("R", "L")])
    mins, maxs = pts.min(axis=0), pts.max(axis=0)
    mins[1] -= 0.12
    maxs[1] += 0.12
    mins[2] = min(mins[2], 0.0)
    span = maxs - mins
    ax3d.set_xlim(mins[0], maxs[0])
    ax3d.set_ylim(mins[1], maxs[1])
    ax3d.set_zlim(mins[2], maxs[2])
    ax3d.set_box_aspect(tuple(span), zoom=1.45)
    ax3d.view_init(elev=13, azim=-62)
    ax3d.set_proj_type("ortho")
    ax3d.set_axis_off()
    ground = [[(mins[0], mins[1], 0.0), (maxs[0], mins[1], 0.0), (maxs[0], maxs[1], 0.0), (mins[0], maxs[1], 0.0)]]
    ax3d.add_collection3d(Poly3DCollection(ground, facecolor="#EFEFEF", edgecolor=GROUND, lw=0.8, alpha=0.9))
    for side in ["R", "L"]:
        ankle = paths[(side, "ankle")]
        ax3d.plot(ankle[:, 0], ankle[:, 1], ankle[:, 2], color=SIDE_COLORS[side], lw=0.8, alpha=0.35)
    skeleton_lines = {}
    for side in ["L", "R"]:
        (limb,) = ax3d.plot([], [], [], color=SIDE_COLORS[side], lw=2.6, marker="o", ms=4, solid_capstyle="round")
        (foot,) = ax3d.plot([], [], [], color=SIDE_COLORS[side], lw=2.6, solid_capstyle="round")
        skeleton_lines[side] = (limb, foot)
    (pelvis_line,) = ax3d.plot([], [], [], color=INK, lw=2.6, solid_capstyle="round")
    ax3d.text2D(0.5, -0.07, "Pelvis-fixed frame, forward = +X", transform=ax3d.transAxes, ha="center", fontsize=8, color=MUTED)

    # Stance bars under the skeleton.
    ax_phase = fig.add_subplot(grid[2, 1])
    stance_spans = {"R": [(0.0, toe_off)], "L": [(50.0, 100.0), (0.0, toe_off - 50.0)]}
    for row, side in enumerate(["R", "L"]):
        ax_phase.barh(row, 100, left=0, height=0.52, color="#F0F0F0", lw=0)
        for start, end in stance_spans[side]:
            ax_phase.barh(row, end - start, left=start, height=0.52, color=SIDE_COLORS[side], lw=0)
    ax_phase.set_yticks([0, 1])
    ax_phase.set_yticklabels(["Right", "Left"], fontsize=8)
    ax_phase.set_ylim(1.6, -0.6)
    ax_phase.set_xlim(0, 100)
    ax_phase.set_xticks([0, 25, 50, 75, 100])
    ax_phase.tick_params(labelsize=8, length=0, axis="y")
    ax_phase.set_xlabel("Gait cycle (%)", fontsize=8)
    ax_phase.set_title("Stance (coloured) and swing (grey)", fontsize=9, loc="left", pad=3)
    ax_phase.spines["left"].set_visible(False)
    phase_cursor = ax_phase.axvline(0, color=INK, lw=1.4)

    # 3. Ankle command: sagittal path and the two time series a robot would track.
    ax_path = fig.add_subplot(grid[0, 2])
    ax_x = fig.add_subplot(grid[1, 2])
    ax_z = fig.add_subplot(grid[2, 2], sharex=ax_x)
    right_cmd = command[("R", "ankle")]
    closed = np.vstack([right_cmd, right_cmd[:1]])
    ax_path.plot(closed[:, 0], closed[:, 2], color=MUTED, lw=1.4)
    pad_x = 0.06 * (right_cmd[:, 0].max() - right_cmd[:, 0].min())
    ax_path.set_xlim(right_cmd[:, 0].min() - pad_x, right_cmd[:, 0].max() + pad_x)
    ax_path.set_ylim(right_cmd[:, 2].min() - 0.03, right_cmd[:, 2].max() + 0.03)
    ax_path.set_aspect("equal", adjustable="box", anchor="W")
    ax_path.set_title("Sagittal path, X vs Z (m), same for both sides", fontsize=9, loc="left", pad=3)
    ax_path.tick_params(labelsize=8)
    ax_path.grid(True, color=GRID, lw=0.6)
    ax_path.set_axisbelow(True)
    command_dots = {}
    for side in ["R", "L"]:
        arr = command[(side, "ankle")]
        ax_x.plot(pct, arr[:, 0], color=SIDE_COLORS[side], lw=1.6)
        ax_z.plot(pct, arr[:, 2], color=SIDE_COLORS[side], lw=1.6)
        (dot_path,) = ax_path.plot([], [], "o", color=SIDE_COLORS[side], ms=7, mec="white", mew=1.2, zorder=5)
        (dot_x,) = ax_x.plot([], [], "o", color=SIDE_COLORS[side], ms=7, mec="white", mew=1.2, zorder=5)
        (dot_z,) = ax_z.plot([], [], "o", color=SIDE_COLORS[side], ms=7, mec="white", mew=1.2, zorder=5)
        command_dots[side] = (arr, dot_path, dot_x, dot_z)
    for ax, label in [(ax_x, "Ankle X, centred (m)"), (ax_z, "Ankle Z (m)")]:
        ax.set_title(label, fontsize=9, loc="left", pad=3)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8)
    ax_x.tick_params(labelbottom=False)
    ax_z.set_xlabel("Gait cycle (%)", fontsize=8)
    phase_text = fig.text(0.985, 0.975, "", fontsize=10, ha="right", va="top", color=INK)

    frames = np.arange(0, n, 2)

    def update(frame_no: int):
        i = int(frames[frame_no])
        j = (i + half) % n
        for mean, dot_r, dot_l in angle_dots:
            dot_r.set_data([pct[i]], [mean[i]])
            dot_l.set_data([pct[j]], [mean[j]])
        for side in ["R", "L"]:
            limb, foot = skeleton_lines[side]
            limb_pts = np.vstack([paths[(side, joint)][i] for joint in LIMB_CHAIN])
            foot_pts = np.vstack([paths[(side, joint)][i] for joint in FOOT_CHAIN])
            limb.set_data_3d(limb_pts[:, 0], limb_pts[:, 1], limb_pts[:, 2])
            foot.set_data_3d(foot_pts[:, 0], foot_pts[:, 1], foot_pts[:, 2])
            arr, dot_path, dot_x, dot_z = command_dots[side]
            dot_path.set_data([arr[i, 0]], [arr[i, 2]])
            dot_x.set_data([pct[i]], [arr[i, 0]])
            dot_z.set_data([pct[i]], [arr[i, 2]])
        hips = np.vstack([paths[("R", "hip")][i], paths[("L", "hip")][i]])
        pelvis_line.set_data_3d(hips[:, 0], hips[:, 1], hips[:, 2])
        phase_cursor.set_xdata([pct[i], pct[i]])
        phase_text.set_text(f"{pct[i]:3.0f}% of the right gait cycle")
        return []

    animation = FuncAnimation(fig, update, frames=len(frames), interval=50, blit=False)
    save_gif(animation, "pipeline_overview.gif", fps=20, dpi=90)
    plt.close(fig)


def make_speed_comparison(age_y: float, height_m: float) -> None:
    """Overground walking of the same child at three reference speeds, in real time."""
    leg_length_m = 0.530 * height_m
    duration_s = 3.6
    fps = 20
    lanes = []
    for label in ["Slow", "Free", "Fast"]:
        v_star = g.SPEED_GROUP_MEANS[label]
        speed = v_star * math.sqrt(g.G * leg_length_m)
        run = reconstruct(g.PatientSpec(f"speed_{label.lower()}", age_y, height_m, speed))
        paths = joint_paths(run["world"])
        n = len(paths[("C", "pelvis")])
        stride_m = 2.0 * run["step_length_m"]
        cycles = int(math.ceil(duration_s / run["stride_time_s"])) + 1
        tiled = {}
        for key, arr in paths.items():
            forward = -arr[:, 0]
            xs = np.concatenate([forward + k * stride_m for k in range(cycles)])
            zs = np.tile(arr[:, 2], cycles)
            tiled[key] = (xs, zs)
        offset = 0.35 - tiled[("R", "hip")][0][0]
        times = np.arange(cycles * n) * run["stride_time_s"] / n
        lanes.append({"label": label, "run": run, "tiled": tiled, "times": times, "offset": offset, "speed": speed})

    walkway_m = max(lane["speed"] for lane in lanes) * duration_s + 1.0
    fig, axes = plt.subplots(3, 1, figsize=(9.2, 5.6), sharex=True)
    fig.subplots_adjust(left=0.045, right=0.985, top=0.86, bottom=0.09, hspace=0.42)
    fig.suptitle(
        f"Same child ({age_y:g} y, {height_m:.2f} m), three reference speeds: overground reconstruction in real time",
        fontsize=11, x=0.045, ha="left", y=0.975,
    )
    fig.text(0.045, 0.915, "Speed changes the joint angles, the step length and the cadence, not just the playback rate.",
             fontsize=9, color=MUTED)
    time_text = fig.text(0.985, 0.975, "", fontsize=10, ha="right", va="top")

    def sample(lane: dict, key: tuple[str, str], t: float) -> tuple[float, float]:
        xs, zs = lane["tiled"][key]
        return float(np.interp(t, lane["times"], xs)) + lane["offset"], float(np.interp(t, lane["times"], zs))

    artists = []
    for ax, lane in zip(axes, lanes):
        run = lane["run"]
        ax.set_xlim(0, walkway_m)
        ax.set_ylim(-0.07, 0.80)
        ax.set_aspect("equal", adjustable="box")
        ax.axhline(0, color=GROUND, lw=1.0)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="x", length=3, labelsize=8)
        ax.set_title(
            f"{lane['label']}:  {lane['speed']:.2f} m/s   |   step length {run['step_length_m']:.2f} m   |   "
            f"cadence {run['cadence_steps_min']:.0f} steps/min",
            fontsize=9, loc="left", pad=4,
        )
        lines = {}
        for side in ["L", "R"]:
            (limb,) = ax.plot([], [], color=SIDE_COLORS[side], lw=2.4, marker="o", ms=3.5, solid_capstyle="round",
                              zorder=3 if side == "R" else 2)
            (foot,) = ax.plot([], [], color=SIDE_COLORS[side], lw=2.4, solid_capstyle="round", zorder=3 if side == "R" else 2)
            (prints,) = ax.plot([], [], "^", color=SIDE_COLORS[side], ms=6, mec="white", mew=0.8, zorder=1, clip_on=False)
            lines[side] = (limb, foot, prints)
        stride_s = run["stride_time_s"]
        strikes = {}
        for side, phase0 in [("R", 0.0), ("L", 0.5)]:
            strike_times = [(k + phase0) * stride_s for k in range(int(duration_s / stride_s) + 2)]
            strikes[side] = [(t, sample(lane, (side, "heel"), t)[0]) for t in strike_times if t <= duration_s]
        artists.append((lane, lines, strikes))
    axes[-1].set_xlabel("Distance walked (m)   ▲ heel strikes: blue = right, orange = left", fontsize=8)

    frame_count = int(round(duration_s * fps))

    def update(frame_no: int):
        t = frame_no / fps
        for lane, lines, strikes in artists:
            for side in ["R", "L"]:
                limb, foot, prints = lines[side]
                limb_pts = np.array([sample(lane, (side, joint), t) for joint in LIMB_CHAIN])
                foot_pts = np.array([sample(lane, (side, joint), t) for joint in FOOT_CHAIN])
                limb.set_data(limb_pts[:, 0], limb_pts[:, 1])
                foot.set_data(foot_pts[:, 0], foot_pts[:, 1])
                done = [x for (ts, x) in strikes[side] if ts <= t + 1e-9]
                prints.set_data(done, [-0.045] * len(done))
        time_text.set_text(f"t = {t:3.1f} s")
        return []

    animation = FuncAnimation(fig, update, frames=frame_count, interval=1000 / fps, blit=False)
    save_gif(animation, "speed_comparison_overground.gif", fps=fps, dpi=85)
    plt.close(fig)


def make_speed_interpolation_figure(age_y: float, height_m: float, speed_m_s: float = 0.92) -> None:
    """Static figure: five Schwartz speed groups and a waveform interpolated between two of them."""
    run = reconstruct(g.PatientSpec("speed_interpolation", age_y, height_m, speed_m_s))
    plan = run["plan"]
    src = pd.read_csv(g.model_dir() / "lam_schwartz_2008_joint_rotation_reference_long.csv")
    n = 100
    pct = np.linspace(0, 100, n, endpoint=False)
    shades = ["#C9D7E8", "#9DB8D9", "#6E97C6", "#3F74AE", "#1B4F8A"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), sharex=True)
    fig.subplots_adjust(left=0.06, right=0.985, top=0.80, bottom=0.16, wspace=0.26)
    for ax, (name, label) in zip(axes, [("hip_flex", "Hip flexion"), ("knee_flex", "Knee flexion"), ("ankle_dorsi", "Ankle dorsiflexion")]):
        for shade, group in zip(shades, g.SPEED_GROUP_ORDER):
            ax.plot(pct, g.interpolate_profile(src, name, "mean", group, n), color=shade, lw=1.5,
                    label=f"{group} (v* = {g.SPEED_GROUP_MEANS[group]:.2f})")
        ax.plot(pct, profile_curve(run["profiles"], name, "mean"), color=LEFT, lw=2.0, ls=(0, (4, 2)),
                label=f"{speed_m_s:.2f} m/s (v* = {run['v_star']:.2f}): {1 - plan.alpha:.0%} {plan.lower_group} + {plan.alpha:.0%} {plan.upper_group}")
        ax.set_title(label, fontsize=9, loc="left")
        ax.set_xlabel("Gait cycle (%)", fontsize=8)
        ax.set_xlim(0, 100)
        ax.grid(True, color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("Angle (deg)", fontsize=8)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8, frameon=False, ncol=6, loc="upper left", bbox_to_anchor=(0.05, 1.0),
               handlelength=1.6, columnspacing=1.1)
    path = MEDIA_DIR / "speed_interpolation.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path.relative_to(REPO_ROOT)} ({path.stat().st_size / 1024:.0f} kB)")


def copy_generator_outputs() -> None:
    """Run the real CLI entry point and keep two of its figures as example outputs."""
    spec, out_dir, n_points, contact_cfg, options = g.spec_from_config(EXAMPLE_CONFIG)
    g.run_one(spec, out_dir, n_points, contact_cfg, options)
    for name in ["spatial_gait_3d_frontal_lateral_animation.gif", "schwartz_speed_matched_angular_profiles.png"]:
        target = MEDIA_DIR / f"example_output_{name}"
        shutil.copyfile(out_dir / "figures" / name, target)
        print(f"wrote {target.relative_to(REPO_ROOT)} ({target.stat().st_size / 1024:.0f} kB)")


def main() -> None:
    MEDIA_DIR.mkdir(parents=True, exist_ok=True)
    spec, _, n_points, _, options = g.spec_from_config(EXAMPLE_CONFIG)
    run = reconstruct(spec, n_points, options.schwartz_speed_interpolation)
    make_pipeline_overview(run)
    make_speed_comparison(spec.age_y, spec.height_m)
    make_speed_interpolation_figure(spec.age_y, spec.height_m)
    copy_generator_outputs()


if __name__ == "__main__":
    main()
