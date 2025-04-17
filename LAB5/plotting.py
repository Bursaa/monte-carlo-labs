import matplotlib.pyplot as plt
import numpy as np
import os

plt.style.use("dark_background")
# Parametry
k_values = [2, 3, 4, 5, 6, 7]
Ns = [10**k for k in k_values]
function_names = ["C1", "C2", "C3"]
methods = ["elementary", "systematical", "layer"]

# Inicjalizacja słownika wyników
results = {
    f: {m: {"g": [], "sigma": [], "intervals": []} for m in methods}
    for f in function_names
}

# Wczytywanie danych g, sigma oraz intervals przy użyciu np.loadtxt
for fname in function_names:
    for k in k_values:
        filename = f"data/{fname}_values_k={k}.txt"
        intervals_filename = f"data/{fname}_intervals_k={k}.txt"
        data = np.loadtxt(filename)
        intervals_data = np.loadtxt(intervals_filename)

        for i, method in enumerate(methods):
            g, sigma = data[i]
            intervals = intervals_data[
                :, i
            ]  # Zakładając, że dane w 'intervals' mają formę [M, liczba metod]

            results[fname][method]["g"].append(g)
            results[fname][method]["sigma"].append(np.abs(sigma))
            results[fname][method]["intervals"].append(intervals)

# Wykresy g ± sigma, R oraz histogramy intervals
for fname in function_names:
    fig, ax = plt.subplots(2, 1, figsize=(8, 12), sharex=True)

    # g ± sigma
    for method in methods:
        g_vals = results[fname][method]["g"]
        sigma_vals = results[fname][method]["sigma"]
        ax[0].errorbar(
            Ns[: len(g_vals)],
            g_vals,
            yerr=sigma_vals,
            fmt="o-",
            capsize=4,
            label=method,
        )
    ax[0].set_xscale("log")
    ax[0].set_ylabel("g")
    ax[0].set_title(f"{fname}: g ± σ oraz R = (σ/g)·100%")
    ax[0].legend()

    # R = (sigma / g) * 100%
    for method in methods:
        g_vals = results[fname][method]["g"]
        sigma_vals = results[fname][method]["sigma"]
        R_vals = [
            100 * sigma / g if g != 0 else 0 for g, sigma in zip(g_vals, sigma_vals)
        ]
        ax[1].plot(Ns[: len(R_vals)], R_vals, marker="s", label=method)
    ax[1].set_xscale("log")
    ax[1].set_ylabel("R [%]")
    ax[1].set_xlabel("N")
    ax[1].legend()
    fig.savefig(f"plots/{fname}_g_R.pdf")
for fname in function_names:
    fig, ax = plt.subplots(3, 1, figsize=(6, 10))
    fig.suptitle(f"Histogramy dla {fname}")

    for i, method in enumerate(methods):
        ax[i].set_title(f"{method.capitalize()}")
        ax[i].set_xlabel("Numer przedziału")
        ax[i].set_ylabel("Częstotliwość")

        # Histogramy intervals dla danej metody
        for k in k_values:
            intervals = results[fname][method]["intervals"][
                k - 2
            ]  # Zakładając, że k_values zaczynają się od 2

            ax[i].bar(np.arange(1, len(intervals) + 1), intervals, label=f"k={k}")
        ax[i].legend()
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)
    fig.savefig(f"plots/{fname}_histograms.pdf")
    plt.show()
