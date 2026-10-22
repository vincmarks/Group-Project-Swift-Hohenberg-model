# Group-Project-Swift-Hohenberg model

This is a repo to work on the group project.

## Model

For the phase field crystal model:

```math
\partial_t \phi = M \nabla^2 \frac{\delta F}{\delta \phi},
\qquad
F = \int \mathrm{d}^d x \left[ \frac{1}{2}\, \phi \left( -\epsilon + \left( q_0^2 + \nabla^2 \right)^2 \right) \phi + \frac{\phi^4}{4} \right]
```

Treating the linear term implicitly and the nonlinear term $\phi^3$ explicitly in Fourier space, we have found the following update formula:

```math
\hat{\phi}_{i+1,m} = \frac{\hat{\phi}_{i,m} - \Delta t\, M\, k_m^2 \sum_{j} (\phi_{i,j})^3\, e^{-i k_m x_j}}{1 + \Delta t\, M\, k_m^2 \left( (q_0^2 - k_m^2)^2 - \epsilon \right)}
```

Here $i$ labels the time step, $m$ the wavenumber $k_m$ and $j$ the grid point $x_j$. The sum over $j$ is the discrete Fourier transform of $\phi^3$, computed with FFTW. In 2D, $k_m^2$ becomes $k_x^2 + k_y^2$.

## Reproducing the results

We use **Julia 1.12.6** 

```bash
git clone https://github.com/vincmarks/Group-Project-Swift-Hohenberg-model.git
cd Group-Project-Swift-Hohenberg-model
julia setup.jl
```

or, from a Julia REPL started in the repository folder:

```julia
include("setup.jl") # this may take a few minutes. Go get a coffee! ;)
```

`setup.jl` does everything in one go:

1. activates the project environment and installs the pinned package versions,
2. deletes the `Plots` folder and creates it again, so it only contains results from the current run,
3. runs `pfc.jl`(results obtained from the PFC simulation) and then `gsh.jl`(results obtained from the generalized Swift-Hohenberg simulation).


All figures and GIFs are saved in the `Plots` folder