# Spin System Animation with Exchange and Dzyaloshinskii-Moriya Interactions

This project provides an animation of a spin system subject to both **exchange interactions** and **Dzyaloshinskii-Moriya (DM) interactions**.

Essentially, I used **gnuplot** to create a real-time animation of the evolution of a system under ferromagnetic exchange and antisymmetric (DM) interactions, using **Monte Carlo simulations** with the **Metropolis algorithm**. The code also implements an **over-relaxation algorithm** to accelerate the dynamics of the modulated structures that emerge from the competition between the two interactions: a **skyrmion**.

## Phases

Depending on the internal parameters (temperature `T`, external magnetic field `H`, and anisotropy `B`), different modulated phases can appear:

- **Helical phase**: characterized by stripe patterns in two dimensions.  
- **Skyrmion phase**: a "bubble" phase with antiparallel cores and edges, exhibiting a vortex-like pattern between them.

## Running the Code

The file is self-contained. To run it, simply download it, compile with:

```bash
gcc animacao.c -lm -O3
```
and execute the program, piping it to gnuplot:
```
./a.out | gnuplot
```
## Modifying the Plots

Plot settings can be modified in two functions inside the code:

initGnuplot() – Contains the plot characteristics: window size, axis labels, tick marks, colorbar display, etc.

printGnuplot() – Generates the actual plot, converting a 1D array into a 2D grid, producing the color map first and then overlaying the vector plot.

Note: This project is intended for educational and research purposes in the study of spin dynamics and skyrmion formation.

## Stripes pattern simulation
![stripes gif](gifs/stripes.gif)


## Skyrmion lattice simulation
![skyrmions gif](gifs/skyrmions.gif)
