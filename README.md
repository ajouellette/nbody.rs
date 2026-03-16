nbody.rs

## Goal
Learn Rust by writing a simple N-body simulation code.

## Inspirations
- [GADGET](https://wwwmpa.mpa-garching.mpg.de/gadget4/) - Perhaps the most well known cosmological N-body simulation codes.
- [The Art of Computational Science](http://www.artcompsci.org/) by Piet Hut and Jun Makino
- [nbody-python](https://github.com/pmocz/nbody-python) - A very simple N-body simulation code in Python.

## Implemented features
- [x] Working simulation using leapfrog or Euler integration
- [x] System of units with conversions between simulation and physical units
- [x] Code is generic over number of spatial dimensions (meant for up to and including 3D)

## Todo
No particular order, no guarantee that anything will actually get implemented.
## Simulation features
- [ ] Periodic boundary conditions
- [ ] Cosmology
- [ ] Generate initial conditions (Plummer sphere, cosmological initial conditions, ...)
- [ ] Hermite integrator
### I/O
- [ ] Logging: log simulation progress, energy conservation
- [ ] Save/load simulation state, output particle snapshots
### Config / Simulation parameters
- [ ] Setup simulation by loading from a yaml file
- [ ] Command line arguments should override config file
### Performance
- [ ] Implement efficient approximation schemes for accelerations (trees)
- [ ] Parallelization?
### Post-processing
- [ ] Python interface to visualize outputs
