import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fps = 50
begin_frames = [list(i for _ in range(fps // 4)) for i in range(1, 11, 1)]
flattened = [i for sublist in begin_frames for i in sublist]
frames_range = np.array(
    list(flattened) + list(range(1, 11, 1)) + list(range(12, 2502, 2))
)
param_labels = ["Density", "Pressure", "Temperature", "Velocity", "Flux j", "n"]

plt.style.use("dark_background")


def get_limits(katalog_idx, frames_range):
    results = {
        "density": [float("inf"), float("-inf")],
        "pressure": [float("inf"), float("-inf")],
        "temperature": [float("inf"), float("-inf")],
        "velocity": [float("inf"), float("-inf")],
        "flux_j": [float("inf"), float("-inf")],
        "n": [float("inf"), float("-inf")],
    }

    for frame_idx in frames_range:
        katalog = f"data/wyniki{katalog_idx}"
        # === NPTV ===
        try:
            data = np.loadtxt(f"{katalog}/nptv_{frame_idx}.dat")
            for i, key in enumerate(
                ["density", "pressure", "temperature", "velocity", "flux_j", "n"]
            ):
                col = data[:, i + 1]  # kolumny 1–6 (0 to x)
                results[key][0] = min(results[key][0], col.min())
                results[key][1] = max(results[key][1], col.max())
        except Exception:
            pass

    return results


def make_hist_update(ax, katalog_idx, xlim=None, ylim=None):
    def update(frame_idx):
        ax.clear()
        katalog = f"data/wyniki{katalog_idx}"
        filename = f"{katalog}/hist_{frame_idx}.dat"
        try:
            data = np.loadtxt(filename)
            ax.scatter(
                data[:, 0], data[:, 1], s=10, label="Numerical", alpha=0.7, zorder=2
            )
            ax.plot(
                data[:, 0], data[:, 2], color="#ff7f0e", label="Analytical", zorder=1
            )
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
            ax.set_xlabel("Velocity")
            ax.set_ylabel("Number of occurrences")
            ax.legend()
            ax.set_title(f"Histogram - Frame {frame_idx} - Zadanie {katalog_idx}")
        except Exception as e:
            print(f"Błąd w {filename}: {e}")

    return update


def make_xy_update(ax, katalog_idx, xlim=None, ylim=None):
    def update(frame_idx):
        ax.clear()
        katalog = f"data/wyniki{katalog_idx}"
        filename = f"{katalog}/xy_{frame_idx}.dat"
        try:
            data = np.loadtxt(filename)
            ax.scatter(data[:, 0], data[:, 1], s=5, alpha=0.6)
            ax.set_title(f"XY - Frame {frame_idx} - Zadanie {katalog_idx}")
            ax.set_xlabel("X")
            ax.set_ylabel("Y")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)
        except Exception as e:
            print(f"Błąd w {filename}: {e}")

    return update


def make_nptv_update(axes, fig, katalog_idx, lims):
    def update(frame_idx):
        katalog = f"data/wyniki{katalog_idx}"
        filename = f"{katalog}/nptv_{frame_idx}.dat"
        try:
            data = np.loadtxt(filename)
            x = data[:, 0]
            for i in range(6):
                axes[i].clear()
                axes[i].plot(x, data[:, i + 1])
                axes[i].set_title(f"{param_labels[i]}")
                axes[i].set_xlabel("x")
                axes[i].set_ylabel("value")
                if lims:
                    key = param_labels[i].lower().replace(" ", "_")
                    if key in lims:
                        axes[i].set_ylim(lims[key][0], lims[key][1])
            fig.suptitle(
                f"NPTV - Frame {frame_idx} - Zadanie {katalog_idx}", fontsize=12
            )
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        except Exception as e:
            print(f"Błąd w {filename}: {e}")

    return update


# Główna pętla po katalogach
for katalog_idx in range(1, 7):
    print(f"Przetwarzanie katalogu {katalog_idx}...")
    limits = get_limits(katalog_idx, frames_range)
    print(f"Tworzenie animacji dla katalogu {katalog_idx}...")

    # HIST
    fig_hist, ax_hist = plt.subplots()
    ani_hist = FuncAnimation(
        fig_hist,
        make_hist_update(ax_hist, katalog_idx),  # type: ignore
        frames=frames_range,
        interval=1000 / fps,
    )
    ani_hist.save(f"plots/hist_animation_{katalog_idx}.mp4", writer="ffmpeg", fps=fps)

    # XY
    fig_xy, ax_xy = plt.subplots()
    ani_xy = FuncAnimation(
        fig_xy, make_xy_update(ax_xy, katalog_idx, ylim=[0.0, 0.5] if katalog_idx == 6 else None), frames=frames_range, interval=1000 / fps  # type: ignore
    )
    ani_xy.save(f"plots/xy_animation_{katalog_idx}.mp4", writer="ffmpeg", fps=fps)

    # NPTV
    fig_nptv, axes_nptv = plt.subplots(2, 3, figsize=(12, 6))
    axes_nptv = axes_nptv.flatten()
    ani_nptv = FuncAnimation(
        fig_nptv,
        make_nptv_update(axes_nptv, fig_nptv, katalog_idx, limits),  # type: ignore
        frames=frames_range,
        interval=1000 / fps,
    )
    ani_nptv.save(f"plots/nptv_animation_{katalog_idx}.mp4", writer="ffmpeg", fps=fps)
