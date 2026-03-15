pub mod vec;

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
    // pub // todo: units
    pub particles: Vec<Particle<NDIM>>,
}

impl<const NDIM: usize> SimState<NDIM> {
    pub fn new(curr_time: f64, softening: f64, particles: Vec<Particle<NDIM>>) -> Self {
        Self {
            curr_time,
            softening,
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
                tot_pot -= self.particles[i].mass * self.particles[j].mass / dist
            }
        }
        tot_pot
    }

    // fn center_of_mass(&self) -> NdVec<NDIM> {
    //     let total_mass: f64 = self.particles.iter().map(|p| p.mass).sum();
    // }

    pub fn accelerate(&mut self) {
        for i in 0..self.particles.len() {
            self.particles[i].acc = [0.0; NDIM];
        }
        for i in 0..self.particles.len() {
            for j in i + 1..self.particles.len() {
                let rvec = self.particles[i].rvec(&self.particles[j]);
                let dist2 = norm2(rvec) + self.softening * self.softening;
                let force = rvec.map(|x| -x / dist2 / dist2.sqrt());
                for d in 0..NDIM {
                    self.particles[i].acc[d] += force[d] * self.particles[j].mass;
                    self.particles[j].acc[d] -= force[d] * self.particles[i].mass;
                }
            }
        }
    }
}

pub fn step_euler<const NDIM: usize>(sim_state: &mut SimState<NDIM>, dt: f64) {
    for i in 0..sim_state.particles.len() {
        for d in 0..NDIM {
            sim_state.particles[i].pos[d] += sim_state.particles[i].vel[d] * dt;
        }
    }
    sim_state.accelerate();
    for i in 0..sim_state.particles.len() {
        for d in 0..NDIM {
            sim_state.particles[i].vel[d] += sim_state.particles[i].acc[d] * dt;
        }
    }
}

pub fn step_leapfrog<const NDIM: usize>(sim_state: &mut SimState<NDIM>, dt: f64) {
    for i in 0..sim_state.particles.len() {
        // kick
        for d in 0..NDIM {
            sim_state.particles[i].vel[d] += sim_state.particles[i].acc[d] * dt / 2.0;
        }
        // drift
        for d in 0..NDIM {
            sim_state.particles[i].pos[d] += sim_state.particles[i].vel[d] * dt;
        }
    }
    sim_state.accelerate();
    for i in 0..sim_state.particles.len() {
        // kick
        for d in 0..NDIM {
            sim_state.particles[i].vel[d] += sim_state.particles[i].acc[d] * dt / 2.0;
        }
    }
}

pub fn integrate<const NDIM: usize>(
    sim_state: &mut SimState<NDIM>,
    total_time: f64,
    dt: f64,
    stepper_fn: &dyn Fn(&mut SimState<NDIM>, f64),
) {
    sim_state.accelerate();
    while sim_state.curr_time + dt <= total_time {
        stepper_fn(sim_state, dt);
        sim_state.curr_time += dt;
    }
    let dt_end = total_time - sim_state.curr_time;
    stepper_fn(sim_state, dt_end);
    sim_state.curr_time += dt_end;
}

// #[cfg(test)]
// mod test {
//     use super::*;

//     #[test]
//     fn test_kinetic() {
//         let parts = vec![
//             Particle {
//                 mass: 1.0,
//                 pos: [0.0; 2],
//                 vel: [1.0, -1.0],
//                 acc: [0.0; 2],
//             },
//             Particle {
//                 mass: 1.0,
//                 pos: [3.0, 4.0],
//                 vel: [-2.0, 1.0],
//                 acc: [0.0; 2],
//             },
//         ];
//         let sim = SimState::new(0.0, 0.0, parts);
//         assert_eq!(sim.energy_kin(), 3.5);
//     }
// }
