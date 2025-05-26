import matplotlib.pyplot as plt
import numpy as np

plt.style.use("dark_background")

V_rel = np.loadtxt("data/Vrel.txt")
nx, ny = V_rel.shape
delta = 0.1
x = np.linspace(0, delta * (nx - 1), nx)
y = np.linspace(0, delta * (ny - 1), ny)
X, Y = np.meshgrid(x, y, indexing="ij")


def plot_map(
    data: np.ndarray, filename: str, title: str, N: int = None, B: int = None
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(X, Y, data, shading="auto", cmap="viridis")
    cb = fig.colorbar(c, ax=ax)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if N is not None and B is not None:
        ax.set_title(rf"{title}(x,y) - N={N} B={B}")
    else:
        ax.set_title(rf"{title}(x,y)")
    cb.set_label(title)
    fig.tight_layout()
    fig.savefig(f"plots/{filename}.pdf")


plot_map(data=V_rel, filename="V_rel", title=r"V_{rel}")
for N in [100, 300]:
    for B in [0, 1]:
        sigma_filename = f"SigmaV_N={N}_B={B}"
        V_mc_filename = f"Vmc_N={N}_B={B}"
        Absorbed_filename = f"Absorbed_N={N}_B={B}"
        sigma_data = np.loadtxt(f"data/{sigma_filename}.txt")
        V_mc_data = np.loadtxt(f"data/{V_mc_filename}.txt")
        Absorbed_data = np.loadtxt(f"data/{Absorbed_filename}.txt")

        plot_map(data=V_mc_data, filename=V_mc_filename, title=r"$V_{MC}$", N=N, B=B)
        plot_map(
            data=np.abs(V_mc_data - V_rel),
            filename=f"delta_V_N={N}_B={B}",
            title="$\Delta V$",
            N=N,
            B=B,
        )
        plot_map(
            data=sigma_data,
            filename=sigma_filename,
            title=r"$\sigma_{V_{MC}}$",
            N=N,
            B=B,
        )
        plot_map(data=Absorbed_data, filename=Absorbed_filename, title=r"$S$", N=N, B=B)


plt.show()
