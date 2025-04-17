import random
import matplotlib.pyplot as plt

plt.style.use("dark_background")
fig1, ax1 = plt.subplots(figsize=(9, 6))
fig2, ax2 = plt.subplots(figsize=(9, 6))
fig3, ax3 = plt.subplots(figsize=(9, 6))
fig4, ax4 = plt.subplots(figsize=(9, 6))
fig5, ax5 = plt.subplots(figsize=(9, 6))
fig6, ax6 = plt.subplots(figsize=(9, 6))

ax1.set_xscale("log")
ax1.set_yscale("log")
ax2.set_xscale("log")
ax2.set_yscale("log")
ax3.set_xscale("log")
ax3.set_yscale("log")
ax4.set_xscale("log")
ax4.set_yscale("log")
ax5.set_xscale("log")
ax5.set_yscale("log")
ax6.set_xscale("log")
ax6.set_yscale("log")

ax1.set_xlabel("N")
ax2.set_xlabel("N")
ax3.set_xlabel("N")
ax4.set_xlabel("N")
ax5.set_xlabel("N")
ax6.set_xlabel("N")

ax1.set_ylabel(r"$\widehat{X}$")
ax2.set_ylabel(r"$\widehat{X^2}$")
ax3.set_ylabel(r"$\Delta X$")
ax4.set_ylabel(r"$\sigma_{\widehat{X}}^2$")
ax5.set_ylabel(r"$\sigma_{Z}^2$")
ax6.set_ylabel(r"$\Delta \sigma^2$")


p_array = [0.1, 0.5, 0.9]
N = 10000001

for px in p_array:
    k = 2
    sumx1 = 0
    sumx2 = 0
    est_x = []
    est_x2 = []
    err_x = []
    var_num = []
    var_teo = []
    err_var = []
    x_axis = []
    for i in range(N):
        u = random.random()
        x = 0
        if u < px:
            x = 1

        sumx1 += x
        sumx2 += x**2
        if i == 10**k:
            x_axis.append(10**k)
            k += 1
            est_x.append(sumx1 / i)
            est_x2.append(sumx2 / i)
            err_x.append(abs((est_x[-1] - px) / px))
            var_num.append((est_x2[-1] - est_x[-1] ** 2) / i)
            var_teo.append((px - px**2) / i)
            err_var.append(abs((var_num[-1] - var_teo[-1]) / var_teo[-1]))
    ax1.plot(x_axis, est_x, label=f"p={px}")
    ax2.plot(x_axis, est_x2, label=f"p={px}")
    ax3.plot(x_axis, err_x, label=f"p={px}")
    ax4.plot(x_axis, var_num, label=f"p={px}")
    ax5.plot(x_axis, var_teo, label=f"p={px}")
    ax6.plot(x_axis, err_var, label=f"p={px}")

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
ax6.legend()

fig1.savefig("outputs/estimation_x.pdf")
fig2.savefig("outputs/estimation_x2.pdf")
fig3.savefig("outputs/error_x.pdf")
fig4.savefig("outputs/variance_numerical.pdf")
fig5.savefig("outputs/variance_theoretical.pdf")
fig6.savefig("outputs/error_variance.pdf")

plt.show()
