import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import chisquare

df = pd.read_csv("data/sampling_results.txt")
df = df.apply(pd.to_numeric, errors="coerce")
Bins = 10
alpha = 0.05
plt.style.use("dark_background")

complex_distribution_array = df["Complex Distribution"].to_numpy()

elimination_method_array = df["Elimination Method"].to_numpy()

marcov_chain_array_d5 = df["Markov Chain (delta=0.5)"].to_numpy()

marcov_chain_array_d05 = df["Markov Chain (delta=0.05)"].to_numpy()
f_vals = df["function"].to_numpy()

x_vals = np.linspace(0, 1, f_vals.shape[0])
fig1, ax1 = plt.subplots(figsize=(8, 6))
fig2, ax2 = plt.subplots(figsize=(8, 6))
fig3, ax3 = plt.subplots(figsize=(8, 6))
fig4, ax4 = plt.subplots(figsize=(8, 6))
ax1.hist(
    complex_distribution_array,
    bins=Bins,
    density=True,
    edgecolor="white",
    label="Complex Distribution",
)
ax1.plot(x_vals, f_vals, label="f(x)", color="red")
ax1.legend()
ax1.set_title("Complex Distribution Histogram")

ax2.hist(
    elimination_method_array,
    bins=Bins,
    density=True,
    edgecolor="white",
    label="Elimination Method",
)
ax2.plot(x_vals, f_vals, label="f(x)", color="red")
ax2.legend()
ax2.set_title("Elimination Method Histogram")

ax3.hist(
    marcov_chain_array_d5,
    bins=Bins,
    density=True,
    edgecolor="white",
    label="Markov Chain (delta=0.5)",
)
ax3.plot(x_vals, f_vals, label="f(x)", color="red")
ax3.legend()
ax3.set_title("Markov Chain (delta=0.5) Histogram")


ax4.hist(
    marcov_chain_array_d05,
    bins=Bins,
    density=True,
    edgecolor="white",
    label="Markov Chain (delta=0.05)",
)
ax4.plot(x_vals, f_vals, label="f(x)", color="red")
ax4.legend()
ax4.set_title("Markov Chain (delta=0.05) Histogram")

plt.tight_layout()
plt.show()

fig1.savefig("outputs/complex_distribution_histogram.pdf")
fig2.savefig("outputs/elimination_method_histogram.pdf")
fig3.savefig("outputs/marcov_chain_delta=0.5.pdf")
fig4.savefig("outputs/marcov_chain_delta=0.05.pdf")


# Chi-Square Test

for name in [
    "Complex Distribution",
    "Elimination Method",
    "Markov Chain (delta=0.5)",
    "Markov Chain (delta=0.05)",
]:
    data = df[name].to_numpy()
    observed, _ = np.histogram(data, bins=Bins)
    expected, _ = np.histogram(f_vals, bins=Bins)
    chi2_stat, p_value = chisquare(observed, expected)
    print(name)
    print(f"chi2 = {chi2_stat:.2f}, p = {p_value:.5f}")
    if p_value < alpha:
        print(
            " Odrzucamy hipotezę zerową (rozkład różni się istotnie od teoretycznego)"
        )
    else:
        print("Brak podstaw do odrzucenia hipotezy zerowej")
