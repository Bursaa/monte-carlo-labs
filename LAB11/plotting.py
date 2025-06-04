import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

plt.style.use("dark_background")

NX = 100
NY = 100
N = 1e6


def load_data(case_id):
    absorption = np.loadtxt(f"data/absorption_{case_id}.dat")
    reflectance = np.loadtxt(f"data/reflectance_{case_id}.dat")
    transmittance = np.loadtxt(f"data/transmittance_{case_id}.dat")
    return absorption, reflectance, transmittance


def calculate_art(absorption, reflectance, transmittance):
    A = np.sum(absorption)
    R = np.sum(reflectance)
    T = np.sum(transmittance)
    norm = (A + R + T) / N
    return A, R, T, norm


def plot_absorption(case_id, absorption, A, R, T, norm):
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(
        absorption.T,
        origin="lower",
        cmap="viridis",
        norm=LogNorm(vmin=1e-1, vmax=np.max(absorption)),
    )
    fig.colorbar(im, label="Absorpcja")
    ax.set_title(
        f"Przypadek #{case_id}\nA={A:.2e}, R={R:.2e}, T={T:.2e}, "
        rf"$\frac{{A+R+T}}{{N}}$={norm:.4f}",
        fontsize=10,
        loc="center",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.tight_layout()
    fig.savefig(f"plots/absorption_plot_{case_id}.pdf")
    plt.show()


def main():
    for case_id in range(1, 9):
        absorption, reflectance, transmittance = load_data(case_id)
        A, R, T, norm = calculate_art(absorption, reflectance, transmittance)
        plot_absorption(case_id, absorption, A, R, T, norm)


if __name__ == "__main__":
    main()
