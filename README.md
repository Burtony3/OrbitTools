# OrbitTools
 Assorted collection of tools and functions for orbital calculations and visualization.

## Code Groups

- `elm`: Code pertaining to Keplarian orbital elements
- `eom`: Functions to be passed to ODE solvers (i.e.: `ode45()`)
- `lambert`: Arrangement of solvers by different author's for Lambert's Problem
- `opt`: Functions to minimization of functions
- `prop`: Analytical propagation functions for orbits
- `spk`: Scripts for interacting, creating, or extracting ephemerides

## Usage

Code from each group can be accessed through dot indexing:

```matlab
odeFunc = @(t, X) eom.Kep(t, X, 398600);
[t, X] = ode45(odeFunc, X0, tspan, ...);
```