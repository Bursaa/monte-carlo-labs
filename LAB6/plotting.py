import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.style.use("dark_background")

# Parametry
dt = 0.1
t_max = 100
D_arr = [r"$D_{xx}$", r"$D_{yy}$", r"$D_{xy}$"]
k_values = [2, 3, 4, 5]
t_arr = np.arange(dt, t_max - dt, dt)

figD, axD = plt.subplots(2, 2, figsize=(10, 10))
for i, k in enumerate(k_values):
    data = np.loadtxt(f"data/D_values_k={k}.txt")
    ax = axD[i % 2, i // 2]
    ax.set_title(f"$N=10^{k}$")
    ax.set_xlabel("t [s]")
    ax.set_ylabel(r"$D_{ii}$")
    for j, D in enumerate(D_arr):
        ax.plot(t_arr, data[:, j], label=D)
    ax.legend()
figD.suptitle(r"$D_{ii}$ time evolution")
figD.savefig("plots/D_ii_time_evolution.pdf")
plt.tight_layout()
plt.show()


def update(frame):
    scat.set_offsets(np.column_stack((x_data[frame], y_data[frame])))
    ax.set_title(f"t = {frame*dt:.1f}")
    return (scat,)


# s_values = [10, 5, 3, 1]
# for i, k in enumerate(k_values):
#     x_file = f"data/Zad1_x(i,t)_k={k}.txt"
#     y_file = f"data/Zad1_y(i,t)_k={k}.txt"

#     # Wczytanie danych
#     x_data = np.loadtxt(x_file)
#     y_data = np.loadtxt(y_file)

#     # Przygotowanie wykresu
#     fig, ax = plt.subplots()
#     fig.suptitle(f"Rozpraszanie cząstek (k={k})")
#     scat = ax.scatter([], [], s=s_values[i])
#     ax.set_xlim(np.min(x_data), np.max(x_data))
#     ax.set_ylim(np.min(y_data), np.max(y_data))
#     ax.set_title("")
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     # Funkcja aktualizująca ramkę
#     # Tworzenie animacji
#     ani = animation.FuncAnimation(
#         fig, update, frames=range(0, x_data.shape[0], 5), interval=10, blit=True
#     )
#     ani.save(f"plots/animacja_ruchy_czastek_k={k}.gif", writer="pillow", fps=30)

#     plt.close()  # opcjonalnie: nie pokazuj interaktywnie


Ra_arr = [0.1, 0.5]
omega_arr = [10, 50, 100]
t_max = 1000.0
t_arr = np.arange(0, t_max, dt)
figN, axN = plt.subplots(figsize=(15, 10))
axN.set_title(f"n(t)")
axN.set_xlabel("t[s]")
axN.set_ylabel("n")
for omega in omega_arr:
    for Ra in Ra_arr:
        data = np.loadtxt(f"data/Zad2_n(t)_omega={omega}_Ra={Ra:.1f}.txt")
        axN.plot(t_arr, data, label=f"$\omega = {omega}$, $R_a = {Ra:.1f}$ ")
axN.legend()
figN.suptitle(r"Active particles time evolution")
figN.savefig("plots/n(t).pdf")
plt.tight_layout()
plt.show()

xr, yr, Rr = 0.0, 0.0, 5.0
xa, ya = 3.0, 0.0
xs, ys = -4.5, 0.0

Ra_arr = [0.1, 0.5]
omega_arr = [10, 50, 100]
duration = 20  # sekundy

for omega in omega_arr:
    for Ra in Ra_arr:
        x_file = f"data/Zad2_x(i,t)_omega={omega}_Ra={Ra:.1f}.txt"
        y_file = f"data/Zad2_y(i,t)_omega={omega}_Ra={Ra:.1f}.txt"

        # Wczytanie danych
        x_data = np.loadtxt(x_file)
        y_data = np.loadtxt(y_file)
        max_rows = 2000
        step = 2
        # Przycinasz dane
        x_data = np.loadtxt(x_file)[:max_rows]
        y_data = np.loadtxt(y_file)[:max_rows]

        # Obliczenia dla animacji
        num_frames = max_rows // step  # = 1000
        fps = num_frames // duration  # = 50
        interval = 1000 / fps  # = 20 ms
        # Obliczenia dla animacji
        num_frames = max_rows // step  # = 1000
        fps = num_frames // duration  # = 50
        interval = 1000 / fps  # = 20 ms

        fig, ax = plt.subplots()
        fig.suptitle(f"Rozpraszanie cząstek ($\\omega={omega}$, $R_a = {Ra}$)")

        # Dodaj okręgi: główny obszar i absorber
        main_circle = plt.Circle((xr, yr), Rr, fill=False, linewidth=1)
        absorber_circle = plt.Circle((xa, ya), Ra, fill=False, linewidth=1)
        ax.add_patch(main_circle)
        ax.add_patch(absorber_circle)

        # Dodaj punkt źródła jako większy marker
        ax.plot(xs, ys, marker="o", markersize=8, label="Źródło")

        scat = ax.scatter([], [], s=1)

        ax.set_xlim(np.min(x_data), np.max(x_data))
        ax.set_ylim(np.min(y_data), np.max(y_data))
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("")
        ax.set_aspect("equal")
        ax.legend(loc="upper right")

        ani = animation.FuncAnimation(
            fig,
            update,
            frames=range(0, max_rows, step),
            interval=interval,
            blit=True,
        )

        ani.save(
            f"plots/Ruchy_czastek_omega={omega}_Ra={Ra:.1f}.gif",
            writer="pillow",
            fps=fps,
        )
        print(f"plots/Ruchy_czastek_omega={omega}_Ra={Ra:.1f}.gif")

        plt.close()
