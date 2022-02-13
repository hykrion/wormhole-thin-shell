# wormhole thin shell
This is faulty code in both C and Wolfram...

The main paper is:

- Adria Delhom, Caio F. B. Macedo, Gonzalo J. Olmo, Luís C. B.
Crispino (2019), Absorption by black hole remnants in metric-affine gravity.
arXiv:1906.06411v1 [gr-qc] DOI 10.1103/PhysRevD.100.024016

This software was developed using [GSL library](https://www.gnu.org/software/gsl/), [glib](https://docs.gtk.org/glib/) and [Wolfram](https://wolfram.com/) (using [jupyter](https://jupyter.org/) with the free Wolfram engine) I've done all the development in Windows so I've used the [MSYS2](https://www.msys2.org/) system. You can use [Chocolatey](https://chocolatey.org/) or install MSYS2 directly.

After that you just have to install the GSL and glib packages using [pacman](https://archlinux.org/pacman/pacman.8.html). pacman is very convenient...

# GUI
The GUI was made with [Tcl/Tk](https://www.tcl.tk/) You need one tcl distribution to use the GUI (but you don't need the GUI to use the program). Another useful tool is [gnuplot](http://www.gnuplot.info/), you can install it with pacman too.

# Partial results
I was not able to extend the Gauss hypergeometrical function, 2_F_1, to calculate equations (4) and (8) using GSL (maybe [Arb](https://arblib.org/) would be a better choice for this...) so I decided to use Wolfram.
Compare the wolfram code to calculate the tortoise coordinates
```
(* Funciones *)
z[x_] = Sqrt[(x^2 + Sqrt[x^2 + 4]) / 2];
zp[x_] = 1 + 1 / z^4 /. z -> z[x];
zm[x_] = 1 - 1 / z^4 /. z -> z[x];
h[x_] = -1 / dc + Sqrt[z^4 - 1] (Hypergeometric2F1[1 / 2, 3 / 4, 3 / 2, 1 - z^4] + Hypergeometric2F1[1 / 2, 7 / 4, 3 / 2, 1 - z^4]) / 2 /. z -> z[x];
a[x_] = (1 - rSbyrc (1 + dc h[x]) / (z Sqrt[zm[x]])) / zp[z] /. z -> z[x];
zyy = a[x] zp[x] D[a[x] zp[x] D[z[x], x], x];
tortoiseXY = NDSolveValue[{x'[y] == a[x[y]] zp[x[y]], x[yL] == yL + 14.}, x, {y, yL, yR}];
```
with the equivalent in C an GSL [tortoise](https://github.com/hykrion/black-bounce-double/blob/main/tortoise.c) There's no colour...

![tortoise](/img/tortoise.png)
![Veff](/img/veff.png)

unfortunately I was not able to calculate the coeffcients R & T to calculate the partial absorption cross section... so I decided to use wolfram to calculate the 2_F_1 and use the results in the C code to get sigma_l. But I was not able to do it... and this's the cause

![Black Bounce GUI](/img/coefficients.png)

as you can see R + T != 1 so I cannot go ahead and calculate sigma_l.

I'm very confident that the correct results are almost in the wolfram code... but some work is need to get to the goal.
In the C code the only flaw is the analytic extension of 2_F_1 as the code has good results in the [Schwarzschild](https://github.com/hykrion/schwarzschild-double) and [Black Bounce](https://github.com/hykrion/black-bounce-double) cases.

# Documentation
For  more information you can consult the pdf's in the *doc* folder (English version it's just an automatic translation)
