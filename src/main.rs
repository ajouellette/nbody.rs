use clap::Parser;
use rand::prelude::*;

use nbody::{integrators, integrators::integrate, units::Units, *};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    // periodic: bool,
    // boxsize: f64,
}

fn main() {
    let args = Args::parse();

    const NDIM: usize = 3;
    let nparts = 200;

    // random initial conditions
    let mut particles = Vec::<Particle<NDIM>>::with_capacity(nparts);
    let mut rng = SmallRng::seed_from_u64(42);
    for _i in 0..nparts {
        particles.push(Particle {
            mass: 1.0,
            pos: rng.random(),
            vel: [0.0; NDIM],
            acc: [0.0; NDIM],
        });
    }

    let units = Units::new_kpc_msun_km_s();
    println!(
        "gravitational constant in sim units: {:.5e}",
        units.grav_const()
    );
    println!(
        "sim time unit in years: {:.5e}",
        units.time_unit() / constants::YR_IN_S
    );

    let mut sim_state = SimState::new(particles, 1e-3, units);
    println!("total particles: {}", sim_state.particles.len());

    let ke = sim_state.energy_kin();
    let pe = sim_state.energy_pot();
    println!(
        "total / kinetic / potential: {:.5e} / {:.5e} / {:.5e}",
        ke + pe,
        ke,
        pe
    );

    integrate(&mut sim_state, 10.0, 1e-3, &integrators::step_leapfrog);

    let ke = sim_state.energy_kin();
    let pe = sim_state.energy_pot();
    println!(
        "total / kinetic / potential: {:.5e} / {:.5e} / {:.5e}",
        ke + pe,
        ke,
        pe
    );
}
