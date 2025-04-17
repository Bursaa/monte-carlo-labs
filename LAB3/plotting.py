import matplotlib.pyplot as plt
import numpy as np

plt.style.use("dark_background")

points = np.loadtxt("data/points.txt")
macierz_norm = np.loadtxt("data/macierz_norm.txt")
macierz_jedn = np.loadtxt("data/macierz_jedn.txt")

r_xy_jedn = macierz_jedn[0, 1] / np.sqrt(macierz_jedn[0, 0] * macierz_jedn[1, 1])
r_xy_norm = macierz_norm[0, 1] / np.sqrt(macierz_norm[0, 0] * macierz_norm[1, 1])

fig1, ax1 = plt.subplots(figsize=(10, 10))
ax1.scatter(points[:, 0], points[:, 1], s=0.2)
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Gauss Distribution in 2D")
fig1.savefig("plots/gauss_distribution.pdf")


fig2, ax2 = plt.subplots(figsize=(10, 10))
ax2.scatter(points[:, 2], points[:, 3], s=0.2)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("Uniform Distribution in 2D")
fig2.savefig("plots/uniform_distribution.pdf")


fig3, ax3 = plt.subplots(figsize=(10, 10))
ax3.scatter(points[:, 4], points[:, 5], s=0.2)
ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_title(
    f"Afinical Transformated Uniform Distribution in 2D r_xy = {r_xy_jedn:.3f}"
)
fig3.savefig("plots/transformed_uniform_distribution.pdf")

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(points[:, 6], points[:, 7], s=0.2)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Afinical Transformated Gauss Distribution in 2D r_xy = {r_xy_norm:.3f}")
fig.savefig("plots/transformed_gauss_distribution.pdf")
plt.show()
