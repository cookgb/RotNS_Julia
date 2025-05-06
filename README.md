# RotNS_Julia
A Julia based code to compute Rotating Neutron Star models

## Motivation
In the early 1990's, I wrote a `FORTRAN` code called `RotNS` to compute equilibrium models for rotating Neutron Stars based on General Relativity.
Four widely referenced papers resulted from this work:
  1. Gregory B. Cook, Stuart L. Shapiro, and Saul A. Teukolsky, “Spin-Up of a Rapidly Rotating Star by Angular Momentum Loss: Effects of General Relativity.” [*Ap. J.* **398** (1992), pp. 203—223](https://ui.adsabs.harvard.edu/abs/1992ApJ...398..203C/abstract).
  2. Gregory B. Cook, Stuart L. Shapiro, and Saul A. Teukolsky, “Rapidly Rotating Polytropes in General Relativity.” [*Ap. J.* **422** (1994), pp. 227—242](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..227C/abstract).
  3. Gregory B. Cook, Stuart L. Shapiro, and Saul A. Teukolsky, “Rapidly Rotating Neutron Stars in General Relativity: Realistic Equations of State.” [*Ap. J.* **424** (1994), pp. 823—845](https://ui.adsabs.harvard.edu/abs/1994ApJ...424..823C/abstract).
  4. Gregory B. Cook, Stuart L. Shapiro, and Saul A. Teukolsky, “Recycling Pulsars to Millisecond Periods in General Relativity.” [*Ap. J. Letters* **423** (1994), pp. 117—120](https://ui.adsabs.harvard.edu/abs/1994ApJ...423L.117C/abstract).

I was interested in learning more about the [Julia programming language](https://julialang.org/) and thought that a fun project would be to get some of my undergraduate students to help me create a new version of this historic code.  And, I also hope that the new, public code might serve as a useful tool for both researchers and educators.

## Status
### Equations of State (EOS)
Constructing Neutron-Star models requires that we have some EOS for nuclear matter.  A simple approximation is to simply use a polytropic equations of state which can be evaluated analytically, but many tabulated EOS exist which represent various theoretical models for how matter should behave down through nuclear densities.  All of the nuclear EOS used in [*Ap. J.* **424** (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...424..823C/abstract) can be found in my public [Neutron_Star_EOS](https://github.com/cookgb/Neutron_Star_EOS.git) GitHub repository.  This repository also contains a Julia code that can convert these individual ASCI text tables into a single [HDF5](https://www.hdfgroup.org/solutions/hdf5/) file.  A modern repository for various EOS, including cold Neutron Star EOS, can be found at the [CompOSE web site](https://compose.obspm.fr/home/).

#### `RotNS_EOS.jl`
`RotNS_EOS.jl` is a collection of Julia data structures and routines that can encapsulate all of the functionality needed to deal with either a polytropic or realistic equation of state.  To create a polytropic EOS opject, all you needs to do is to create a `PolytropicEOS` object based on a chosen polytropic index $n$.
```
eos = PolytropicEOS(1.5)
```
Similarly, to create a `RealisticEOS` object, you simply need to tell it where the `HDF5` file resides and which EOS to choose.
```
eos = RealisticEOS("./Tabulated_EOS.h5","FPS")
```
Given an `EOS` object, it is then straightforward to compute various quantities of interest.  For example, if you know the pressure $P$, then it the corresponding total energy density $\epsilon$ can be computed using
```
ϵ = EnergyDensityfromPressure(eos,P)
```
`RotNS_EOS.jl` provides the functionality to deal with pressure $P$, total energy density $\epsilon$, baryon-mass density $\rho_0$, and enthalpy $h$.

#### `RotNS_TOV.jl`
A first step in computing rotating neutron star models is to compute non-rotating models.  This is accomplished by solving the [Tolman–Oppenheimer–Volkoff(TOV) equations](https://en.wikipedia.org/wiki/Tolman%E2%80%93Oppenheimer%E2%80%93Volkoff_equation).  `RotNS_TOV.jl` provides the data structures and routines needed to solve the TOV equations in either standard Schwarzschild coordinates, or in the isotropic coordinates needed to work with rotating models.  Any non-rotating Neutron-Star model requries a choice of the EOS and the value of the central total energy density $\epsilon_c$.  The integration scheme also needs an idea of what the radius of the star will be, so we must provide an upper limit on the radius.  To construct a TOV model using Schwarzschild coordinates, we can do the following.
```
eos = RealisticEOS("./Tabulated_EOS.h5","FPS")
ϵ = 2.0
Rmax = 2.0
solution = TOVSchwarzschild(eos,ϵ,Rmax)
```
Alternatively, to construct a TOV model in isotropic coordinates, we would call `TOVModel`.
```
solution = TOVModel(eos,ϵ,Rmax)
```

#### The next step is to construct rotating models!
