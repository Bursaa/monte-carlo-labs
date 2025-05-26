import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np

plt.style.use("dark_background")


def mean_r(wielokaty):
    points = np.array([punkt for wielokat in wielokaty for punkt in wielokat])
    dist = np.linalg.norm(points, axis=1)
    return np.mean(dist)


def wczytaj_wielokaty(plik):
    wielokaty = []
    blok = []
    puste_linie = 0

    with open(plik, "r") as f:
        for linia in f:
            linia = linia.strip()
            if linia == "":
                puste_linie += 1
                if puste_linie == 2 and blok:
                    wielokaty.append(blok)
                    blok = []
                    puste_linie = 0
            else:
                puste_linie = 0
                try:
                    x, y, z = map(float, linia.split())
                    blok.append([x, y, z])
                except ValueError:
                    continue
        if blok:
            wielokaty.append(blok)
    return wielokaty


def rysuj_wielokaty_3d(wielokaty, name, index=None):
    if index is not None:
        wielokaty = [wielokaty[index]]

    fig = go.Figure()

    for wielokat in wielokaty:
        if len(wielokat) < 3:
            continue

        punkty = np.array(wielokat)
        # jeśli nie zamknięty, zamknij pętlę
        if not np.allclose(punkty[0], punkty[-1]):
            punkty = np.vstack([punkty, punkty[0]])

        x, y, z = punkty[:, 0], punkty[:, 1], punkty[:, 2]

        # linie wielokąta
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines+markers",
                line=dict(color="black", width=5),
                marker=dict(size=6, color="black"),
                name="Wielokąt",
            )
        )

        # powierzchnia - tylko jeśli wielokąt jest płaski (można tu rozwinąć dla ogólnych)
        fig.add_trace(
            go.Mesh3d(
                x=x[:-1],
                y=y[:-1],
                z=z[:-1],  # bez ostatniego punktu domknięcia
                color="gold",
                opacity=0.5,
                alphahull=0,  # wypukła powierzchnia
                showscale=False,
                name="Powierzchnia",
            )
        )

    fig.update_layout(
        title=f"Wielokąty 3D - {name}",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, b=0, t=30),
    )

    fig.write_html(f"plots/3D_fulleren_{name}.html")  # 🔁 Interaktywny plik HTML


name = "zad1"
wielokaty = wczytaj_wielokaty(f"data/polygons_{name}.txt")
rysuj_wielokaty_3d(wielokaty, name)

for i in range(1, 6):
    name = "zad" + str(i)
    wielokaty = wczytaj_wielokaty(f"data/polygons_{name}.txt")
    rysuj_wielokaty_3d(wielokaty, name)
    r_mean = mean_r(wielokaty)

    pcf = np.loadtxt(f"data/pcf_{name}.txt")
    x_values = np.array([m * 2.5 * r_mean / len(pcf) for m in range(len(pcf))])
    dx = x_values[1] - x_values[0]
    fig_pcf, ax_pcf = plt.subplots(figsize=(8, 8))
    ax_pcf.bar(x_values, pcf, width=dx, edgecolor="white")
    ax_pcf.set_xticks(np.arange(min(x_values), max(x_values) + 0.5, 0.5))
    ax_pcf.set_xlabel("Indeks")
    ax_pcf.set_ylabel("Wartość pcf")
    ax_pcf.set_title(f"Wykres słupkowy pcf - {name}")
    fig_pcf.savefig(f"plots/pcf_hist_{name}.pdf")

    if i > 1:
        values = np.loadtxt(f"data/values_{name}.txt")
        iterations = values[:, 0]
        V = values[:, 1]
        r_sr = values[:, 2]
        beta = values[:, 3]

        fig_V, ax_V = plt.subplots(figsize=(8, 8))
        ax_V.plot(iterations, V)
        ax_V.set_xlabel("iteration")
        ax_V.set_ylabel("Energy [eV]")
        ax_V.set_title(f"Energy - {name}")
        ax_V.set_ylim(top=0)
        fig_V.savefig(f"plots/Energy_{name}.pdf")

        fig_r, ax_r = plt.subplots(figsize=(8, 8))
        ax_r.plot(iterations, r_sr)
        ax_r.set_xlabel("iteration")
        ax_r.set_ylabel(r"$r_{sr}$ [A]")
        ax_r.set_title(rf"$r_{{sr}}$ - {name}")
        fig_r.savefig(f"plots/r_sr_{name}.pdf")

        fig_beta, ax_beta = plt.subplots(figsize=(8, 8))
        ax_beta.plot(iterations, beta)
        ax_beta.set_xlabel("iteration")
        ax_beta.set_ylabel(r"$\beta$")
        ax_beta.set_title(rf"$\beta$ - {name}")
        fig_beta.savefig(f"plots/beta_{name}.pdf")
    plt.close("all")

for i in range(30, 41):
    name = "zad6_n=" + str(i)

    wielokaty = wczytaj_wielokaty(f"data/polygons_{name}.txt")
    rysuj_wielokaty_3d(wielokaty, name)
    r_mean = mean_r(wielokaty)

    pcf = np.loadtxt(f"data/pcf_{name}.txt")
    x_values = np.array([m * 2.5 * r_mean / len(pcf) for m in range(len(pcf))])
    dx = x_values[1] - x_values[0]
    fig_pcf, ax_pcf = plt.subplots(figsize=(8, 8))
    ax_pcf.bar(x_values, pcf, width=dx, edgecolor="white")
    ax_pcf.set_xticks(np.arange(min(x_values), max(x_values) + 0.5, 0.5))
    ax_pcf.set_xlabel("Indeks")
    ax_pcf.set_ylabel("Wartość pcf")
    ax_pcf.set_title(f"Wykres słupkowy pcf - {name}")
    fig_pcf.savefig(f"plots/pcf_hist_{name}.pdf")

    values = np.loadtxt(f"data/values_{name}.txt")
    iterations = values[:, 0]
    V = values[:, 1]
    r_sr = values[:, 2]
    beta = values[:, 3]

    fig_V, ax_V = plt.subplots(figsize=(8, 8))
    ax_V.plot(iterations, V)
    ax_V.set_xlabel("iteration")
    ax_V.set_ylabel("Energy [eV]")
    ax_V.set_title(f"Energy - {name}")
    ax_V.set_ylim(top=0)
    fig_V.savefig(f"plots/Energy_{name}.pdf")

    fig_r, ax_r = plt.subplots(figsize=(8, 8))
    ax_r.plot(iterations, r_sr)
    ax_r.set_xlabel("iteration")
    ax_r.set_ylabel(r"$r_{sr}$ [A]")
    ax_r.set_title(rf"$r_{{sr}}$ - {name}")
    fig_r.savefig(f"plots/r_sr_{name}.pdf")

    fig_beta, ax_beta = plt.subplots(figsize=(8, 8))
    ax_beta.plot(iterations, beta)
    ax_beta.set_xlabel("iteration")
    ax_beta.set_ylabel(r"$\beta$")
    ax_beta.set_title(rf"$\beta$ - {name}")
    fig_beta.savefig(f"plots/beta_{name}.pdf")
    plt.close("all")
