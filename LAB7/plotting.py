import numpy as np
import matplotlib.pyplot as plt
import numba as nb
from numba.typed import List

plt.style.use("dark_background")

# Parametry
tmax = 200.0  # maksymalny czas
N = 50  # liczba przedziałów czasowych
delta_t = tmax / N


# Funkcja do obliczania statystyk
@nb.njit
def calculate_statistics(x3_data, Pmax, N):
    h0 = np.zeros(N, dtype=np.float64)  # Sumowanie x3
    h1 = np.zeros(N, dtype=np.float64)  # Sumowanie średnich x3
    h2 = np.zeros(N, dtype=np.float64)  # Sumowanie kwadratów średnich x3
    ncount = np.zeros(
        N, dtype=np.int64
    )  # Liczenie punktów w każdym przedziale czasowym

    # Obliczanie sum dla każdego p (każdego zestawu danych)
    for p in range(Pmax):
        for time_data in x3_data[p]:
            t = time_data[0]  # Czas
            x3 = time_data[3]  # Wartość x3
            l = int(t // delta_t)  # Indeks odpowiadający czasowi
            if l < N:  # Zapewniamy, że indeks mieści się w zakresie
                h0[l] += x3  # Sumowanie wartości x3
                ncount[l] += 1  # Liczymy ilość punktów w danym przedziale czasowym
        for l in range(N):
            x3t = h0[l] / ncount[l]
            h1[l] += x3t
            h2[l] += x3t * x3t

    x3t_sr_1 = h1 / Pmax
    x3t_sr_2 = h2 / Pmax
    sigma_t = np.sqrt((x3t_sr_2 - x3t_sr_1 * x3t_sr_1) / Pmax)

    return x3t_sr_1, sigma_t


# Funkcja do generowania wykresu
def plot_statistics(Pmax_values):
    fig, axs = plt.subplots(
        len(Pmax_values), 1, figsize=(7, 5 * len(Pmax_values))
    )  # Utwórz 3 wykresy obok siebie
    # Przechodzimy przez różne wartości Pmax i tworzymy wykresy
    for idx, Pmax in enumerate(Pmax_values):
        axs[idx].set_xlabel("t")
        axs[idx].set_ylabel(r"$\overline{x_3}$")
        x3_data = List()

        # Zbieranie danych z plików
        for p in range(Pmax):
            data = np.loadtxt(f"data/x(t)_p={p+1}.txt")
            x3_data.append(data)  # Każdy element to (t, x1, x2, x3) z pliku

        # Obliczanie średnich i odchyleń standardowych
        h1, sigma_t = calculate_statistics(x3_data, Pmax, N)

        # Rysowanie wykresu dla danej wartości Pmax
        time = (np.arange(N) + 0.5) * delta_t
        axs[idx].errorbar(
            time,
            h1,
            yerr=sigma_t,
            fmt="o",  # Punkty jako kółka
            capsize=3,
            label=f"$\overline{{x_3}}$ for $Pmax={Pmax}$",
            linestyle="-",  # Łączy punkty linią
            linewidth=0.5,  # Zmniejsza szerokość linii
            markersize=4,  # Zmniejsza wielkość punktów
        )
        axs[idx].legend()
        axs[idx].set_title(f"$\overline{{x_3}}$ for $Pmax={Pmax}$")

    fig.suptitle(r"$\overline{x_3}$ for different $Pmax$")
    fig.tight_layout(rect=[0, 0, 1, 0.96])  # Dostosowanie układu wykresów
    fig.savefig(r"plots/x3_sr_for_different_Pmax.pdf")


fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlabel("t")
ax.set_ylabel(r"$x_{1,2,3}$")
data = np.loadtxt(f"data/x(t)_p=1.txt")
ax.plot(data[:, 0], data[:, 1], color="#1f77b4", label=r"$x_1$")
ax.plot(data[:, 0], data[:, 2], color="#ff7f0e", label=r"$x_2$")
ax.plot(data[:, 0], data[:, 3], color="#2ca02c", label=r"$x_3$")
ax.legend()
fig.suptitle(r"$x_{1,2,3}(t)$ for $p = 1$")
fig.savefig(r"plots/x_n(t)_p=1.pdf")

fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlabel("t")
ax.set_ylabel(r"$x_{1,2,3}$")
for i in range(5):
    data = np.loadtxt(f"data/x(t)_p={i+1}.txt")
    ax.plot(data[:, 0], data[:, 1], color="#1f77b4", label=r"$x_1$" if i == 0 else None)
    ax.plot(data[:, 0], data[:, 2], color="#ff7f0e", label=r"$x_2$" if i == 0 else None)
    ax.plot(data[:, 0], data[:, 3], color="#2ca02c", label=r"$x_3$" if i == 0 else None)

ax.legend()
fig.suptitle(r"$x_{1,2,3}(t)$ for $p \in [1,5]$")
fig.savefig(r"plots/x_n(t)_p=5.pdf")

Pmax_values = [5, 10, 100]
plot_statistics(Pmax_values)
plt.show()
