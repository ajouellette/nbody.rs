pub mod constants;
pub mod integrators;
pub mod units;
pub mod vec;

use crate::units::Units;
use crate::vec::*;

#[derive(Debug)]
pub struct Particle<const NDIM: usize> {
    pub mass: f64,
    pub pos: NdVec<NDIM>,
    pub vel: NdVec<NDIM>,
    pub acc: NdVec<NDIM>,
}

impl<const NDIM: usize> Particle<NDIM> {
    pub fn rvec(&self, other: &Self) -> NdVec<NDIM> {
        let mut res: NdVec<NDIM> = [0.0; NDIM];
        for d in 0..NDIM {
            res[d] = self.pos[d] - other.pos[d]
        }
        res
    }
}

#[derive(Debug)]
pub struct SimState<const NDIM: usize> {
    pub curr_time: f64,
    pub softening: f64,
    // units
    pub units: Units,
    // list of particles
    pub particles: Vec<Particle<NDIM>>,
}

impl<const NDIM: usize> SimState<NDIM> {
    pub fn new(particles: Vec<Particle<NDIM>>, softening: f64, units: Units) -> Self {
        let curr_time = 0.0;
        Self {
            curr_time,
            softening,
            units,
            particles,
        }
    }

    pub fn energy_kin(&self) -> f64 {
        0.5 * self
            .particles
            .iter()
            .map(|p| p.mass * norm2(p.vel))
            .sum::<f64>()
    }

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

    // fn center_of_mass(&self) -> NdVec<NDIM> {
    //     let total_mass: f64 = self.particles.iter().map(|p| p.mass).sum();
    // }

    pub fn accelerate(&mut self) {
        for p in self.particles.iter_mut() {
            p.acc = [0.0; NDIM];
        }
        let eps2 = self.softening * self.softening;
        for i in 0..self.particles.len() {
            for j in i + 1..self.particles.len() {
                let rvec = self.particles[i].rvec(&self.particles[j]);
                let dist2 = norm2(rvec) + eps2;
                let force = rvec.map(|x| -self.units.grav_const() * x / dist2 / dist2.sqrt());
                for d in 0..NDIM {
                    self.particles[i].acc[d] += force[d] * self.particles[j].mass;
                    self.particles[j].acc[d] -= force[d] * self.particles[i].mass;
                }
            }
        }
    }
}
