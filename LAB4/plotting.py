import matplotlib.pyplot as plt
import numpy as np

plt.style.use("dark_background")

n_values = np.logspace(2, 6, num=5, base=10)
KApoints = np.loadtxt("data/KA_points.txt")
KBpoints = np.loadtxt("data/KB_points.txt")

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.scatter(KApoints[:, 0], KApoints[:, 1], s=0.2)
ax1.scatter(KBpoints[:, 0], KBpoints[:, 1], s=0.2)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_aspect("equal")  # Zachowanie proporcji
ax1.set_title("Circles with uniform distribution")
fig1.savefig("plots/circles.pdf")

alphas = ["A", "B"]
rAzeroindicators = ["!=", "="]
fig2, ax2 = plt.subplots(1, 2, figsize=(15, 10))

for i, rAZeroIndicator in enumerate(rAzeroindicators):
    for alpha in alphas:
        filename = f"data/mu_points_{alpha}_xA{rAZeroIndicator}0.txt"
        mu_points = np.loadtxt(filename)

        m1 = mu_points[:, 0]  # Średnie wartości
        sigma = mu_points[:, 1]  # Niepewności

        ax2[i].errorbar(
            n_values,
            m1,
            yerr=sigma,
            fmt="o",
            capsize=3,
            label=f"$\\alpha$ = {alpha}",
        )

    ax2[i].set_xlabel(r"$n$")
    ax2[i].set_ylabel("Common area")
    ax2[i].set_title(
        f"Plot for xA {'= RB + 0.5RA' if rAZeroIndicator == '!=' else '= 0'}"
    )
    ax2[i].set_xscale("log")
    ax2[i].legend()

fig2.tight_layout()
fig2.savefig("plots/mu_points_errorbars.pdf")
