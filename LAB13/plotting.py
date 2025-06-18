import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

plt.style.use("dark_background")

N = 10**8


def plot_energy_maps(filename="data/Energies.dat"):
    data = np.loadtxt(filename, dtype=float)
    a_vals = np.sort(np.unique(data[:, 0]))
    c_vals = np.sort(np.unique(data[:, 1]))

    E = data[:, 2].reshape(len(a_vals), len(c_vals))
    var = data[:, 3].reshape(len(a_vals), len(c_vals))
    sigma = data[:, 4].reshape(len(a_vals), len(c_vals)) + 1e-20

    extent = [a_vals[0], a_vals[-1], c_vals[0], c_vals[-1]]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    titles = [
        r"Energia $\langle E_L \rangle$",
        r"Wariancja $\sigma^2$",
        r"$\log_{10}(\sigma + 10^{-20})$",
    ]
    datasets = [E, var, np.log10(sigma)]

    # Ustaw ticki co 0.1
    xticks = np.round(np.arange(a_vals[0], a_vals[-1] + 0.001, 0.1), 2)
    yticks = np.round(np.arange(c_vals[0], c_vals[-1] + 0.001, 0.1), 2)

    for ax, title, data in zip(axes, titles, datasets):
        im = ax.imshow(
            data.T,
            extent=extent,
            origin="lower",
            aspect="auto",
            cmap="viridis",
        )
        ax.set_title(title)
        ax.set_xlabel("a")
        ax.set_ylabel("c")
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        fig.colorbar(im, ax=ax, label="wartość")

    fig.tight_layout()
    fig.savefig("plots/Energies.pdf")


def generate_log_intervals(max_n: int) -> np.ndarray:
    frames = []
    max_exp = int(np.log10(max_n))

    for exp in range(7, max_exp + 1):
        base = 10**exp
        for i in range(1, 101):
            val = i * base // 100
            if val <= max_n:
                frames.append(val)

    return np.sort(np.array(np.unique(frames)))


def animate_histograms(a, c, ylim):
    frames = generate_log_intervals(N)
    fig, ax = plt.subplots(figsize=(8, 5))
    (line_data,) = ax.plot(
        [], [], label="symulacja", marker="o", linestyle="-", markersize=3
    )
    (line_exact,) = ax.plot([], [], "--", label="dokładna", alpha=0.5)
    ax.set_xlim(0, 8)
    ax.set_ylim(0.0, 1.0)
    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel("r")
    ax.set_ylabel("gęstość prawdopodobieństwa")
    ax.legend()

    def init():
        line_data.set_data([], [])
        line_exact.set_data([], [])
        return line_data, line_exact

    def update(frame):
        data = np.loadtxt(f"data/hist_a={a}_c={c}_{frame}.dat")
        r, p, exact = data[:, 0], data[:, 1], data[:, 2]
        line_data.set_data(r, p)
        line_exact.set_data(r, exact)
        ax.set_title(f"Krok MC = {frame}")
        return line_data, line_exact

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=frames,
        init_func=init,
        interval=100,
    )
    ani.save(f"plots/histogram_a={a}_c={c}.mp4", writer="ffmpeg", fps=10)


print("Plotting Energy, Variance and sigma values")
plot_energy_maps()  # wykres mapy energii i wariancji

print("Animation of histogram 1 ongoing..")
animate_histograms(1, 0, ylim=(0.0, 0.6))


print("Animation of histogram 2 ongoing..")
animate_histograms(0.5, -0.5, ylim=(0.0, 0.3))
