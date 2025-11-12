"""
I prefer the way that these streamplot plots look in Python compared to 
tikz/pgf (difficult because there are ode simulations to do) and Matlab. So 
this function generates the phase portraits for the pendulum example using 
the controller passed to it, which I compute in Matlab. 
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches
from scipy.integrate import odeint
from examples.example26ExportedControlFunction import u
from sympy import *

theta0 = -np.pi
xmax = 3.5
ymax = 4.11
x1_range = np.linspace(-xmax, xmax, 400)
x2_range = np.linspace(-ymax, ymax, 200)
X1, X2 = np.meshgrid(x1_range, x2_range)
U, V = np.zeros_like(X1), np.zeros_like(X2)


def f(x1, x2):
    return [x2, u(x1, x2) - np.sin(x1)]


def pendulum_derivatives(x, t):
    x1, x2 = x
    x1_dot, x2_dot = f(x1, x2)
    return [x1_dot, x2_dot]


def plot_function(x1_range, x2_range, f, xmax, theta0):
    U = np.zeros((len(x2_range), len(x1_range)))
    V = np.zeros((len(x2_range), len(x1_range)))
    for i in range(len(x1_range)):
        for j in range(len(x2_range)):
            x1 = x1_range[i]
            x2 = x2_range[j]
            U[j, i], V[j, i] = f(x1, x2)

    fig, ax = plt.subplots()
    ax.streamplot(
        X1,
        X2,
        U,
        V,
        color="lightgray",
        density=0.8,
        linewidth=1,
        cmap=plt.cm.gray,
        broken_streamlines=False,
    )
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-ymax, ymax)

    x0 = [theta0, np.sqrt(2 * (1 - np.cos(theta0)))]
    x_t = odeint(pendulum_derivatives, x0, np.linspace(0, 5, 1000))
    ax.plot(x_t[:, 0], x_t[:, 1], "r-", linewidth=2)

    ax.set_frame_on(False)
    plt.axis("off")  # this rows the rectangular frame
    ax.get_xaxis().set_visible(False)  # this removes the ticks and numbers for x axis
    ax.get_yaxis().set_visible(False)  # this removes the ticks and numbers for y axis

    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()  # Method 1


plot_function(x1_range, x2_range, f, xmax, theta0)
plt.savefig(
    f"plots/example26_phasePortrait_u{int(controlDegree)}.pdf",
    bbox_inches="tight",
    pad_inches=0,
    transparent=True,
)
plt.savefig(
    "plots/example26_phasePortrait.pdf",
    bbox_inches="tight",
    pad_inches=0,
    transparent=True,
)
plt.close("all")
x1, x2 = symbols("x1 x2")
print("u(x1,x2) = " + str(u(x1, x2)))
