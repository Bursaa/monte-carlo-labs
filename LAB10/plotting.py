import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

plt.style.use("dark_background")


def animate_uxt_with_exact(npath: int) -> None:
    # === WCZYTYWANIE DANYCH MONTE CARLO ===
    mc_filename = f"data/u(x,t)_npaths={npath}.txt"
    exact_filename = "data/u(x,t)_exact.txt"

    x_vals = np.arange(0.0, 2.0 + 1e-3, 0.01)
    mc_data = np.loadtxt(mc_filename)
    t_vals = mc_data[:, 0]
    u_mc = mc_data[:, 1:]

    # === WCZYTYWANIE DANYCH DOKŁADNYCH ===
    exact_data = np.loadtxt(exact_filename)
    u_exact = exact_data[:, 1:]

    # === WYKRES ===
    fig, ax = plt.subplots()
    scat = ax.scatter(x_vals, u_mc[0], s=20, label="MC")
    (line,) = ax.plot(x_vals, u_exact[0], color="#ff5555", label="Exact")
    time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)
    ax.set_xlim(x_vals.min(), x_vals.max())
    ax.set_ylim(min(u_mc.min(), u_exact.min()), max(u_mc.max(), u_exact.max()))
    ax.set_xlabel("x [m]")
    ax.set_ylabel("u(x, t) [V]")
    ax.set_title(f"u(x, t): MC vs Exact  npath={npath}")
    ax.legend(loc="upper right")

    # === AKTUALIZACJA KLATKI ===
    def update(frame):
        scat.set_offsets(np.c_[x_vals, u_mc[frame]])
        line.set_ydata(u_exact[frame])
        time_text.set_text(f"t = {t_vals[frame]*1e9:.1f} ns")
        return scat, line, time_text

    ani = animation.FuncAnimation(
        fig, update, frames=len(t_vals), interval=200, blit=True
    )
    plt.tight_layout()
    gif_filename = mc_filename.replace(".txt", ".gif").replace("data/", "plots/")
    mp4_filename = mc_filename.replace(".txt", ".mp4").replace("data/", "plots/")
    print(f"Zapisuję animację do {gif_filename} i {mp4_filename}...")
    ani.save(gif_filename, writer="pillow", fps=5)
    ani.save(mp4_filename, writer="ffmpeg", fps=5)
    print("Gotowe.")


for npath in [10, 100, 1000, 10000, 100000]:
    animate_uxt_with_exact(npath)
