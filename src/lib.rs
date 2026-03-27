pub mod constants;
pub mod integrators;
pub mod units;
pub mod vec;

use crate::units::Units;
use crate::vec::*;

/// A simulation particle in `NDIM`-dimensional space.
#[derive(Debug)]
pub struct Particle<const NDIM: usize> {
    pub mass: f64,
    pub pos: NdVec<NDIM>,
    pub vel: NdVec<NDIM>,
    pub acc: NdVec<NDIM>,
}

impl<const NDIM: usize> Particle<NDIM> {
    /// Get the position vector to this particle from another.
    pub fn rvec(&self, other: &Self) -> NdVec<NDIM> {
        let mut res: NdVec<NDIM> = [0.0; NDIM];
        for d in 0..NDIM {
            res[d] = self.pos[d] - other.pos[d]
        }
        res
    }
}

/// Current simulation state.
/// Stores all information needed to start or resume a simulation.
#[derive(Debug)]
pub struct SimState<const NDIM: usize> {
    /// Current simulation time.
    pub curr_time: f64,
    /// Softening length for gravitational interactions.
    pub softening: f64,
    /// Simulation units.
    pub units: Units,
    /// List of particles in the simulation.
    pub particles: Vec<Particle<NDIM>>,
}

impl<const NDIM: usize> SimState<NDIM> {
    /// Setup a new simulation state with the given particles, softening length, and units.
    pub fn new(particles: Vec<Particle<NDIM>>, softening: f64, units: Units) -> Self {
        let curr_time = 0.0;
        Self {
            curr_time,
            softening,
            units,
            particles,
        }
    }

    /// Compute the current kinetic energy of the system.
    pub fn energy_kin(&self) -> f64 {
        0.5 * self
            .particles
            .iter()
            .map(|p| p.mass * norm2(p.vel))
            .sum::<f64>()
    }

    /// Compute the current potential energy of the system.
    pub fn energy_pot(&self) -> f64 {
        let mut tot_pot: f64 = 0.0;
        for i in 0..self.particles.len() {
            for j in i + 1..self.particles.len() {
                let dist = norm(self.particles[i].rvec(&self.particles[j]));
                tot_pot -=
                    self.units.grav_const() * self.particles[i].mass * self.particles[j].mass
                        / dist;
            }
        }
        tot_pot
    }

    /// Compute center of mass of the system.
    pub fn center_of_mass(&self) -> NdVec<NDIM> {
        let total_mass: f64 = self.particles.iter().map(|p| p.mass).sum();
        let mut com: NdVec<NDIM> = [0.0; NDIM];
        for p in self.particles.iter() {
            let frac_mass = p.mass / total_mass;
            for d in 0..NDIM {
                com[d] += frac_mass * p.pos[d];
            }
        }
        com
    }

    /// Brute-force update all of the particle accelerations in the system.
    pub fn accelerate(&mut self) {
        for p in self.particles.iter_mut() {
            p.acc = [0.0; NDIM];
        }
        let eps2 = self.softening * self.softening;
        let n = self.particles.len();
        for i in 0..n {
            for j in i + 1..n {
                let rvec = self.particles[i].rvec(&self.particles[j]);
                let dist3_2 = norm2(rvec) + eps2;
                let dist3_2 = dist3_2 * dist3_2.sqrt();
                let force = rvec.map(|x| -self.units.grav_const() * x / dist3_2);
                for d in 0..NDIM {
                    self.particles[i].acc[d] += force[d] * self.particles[j].mass;
                    self.particles[j].acc[d] -= force[d] * self.particles[i].mass;
                }
            }
        }
    }
}
