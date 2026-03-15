use clap::Parser;
use rand::prelude::*;

use nbody::*;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    // periodic: bool,
    // boxsize: f64,
}

fn main() {
    let args = Args::parse();

    const NDIM: usize = 3;
    let nparts = 20;

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

    let mut sim_state = SimState::new(particles, 1e-3);
    println!("total particles: {}", sim_state.particles.len());

    let ke = sim_state.energy_kin();
    let pe = sim_state.energy_pot();
    println!(
        "total / kinetic / potential: {:.5e} / {:.5e} / {:.5e}",
        ke + pe,
        ke,
        pe
    );

    integrate(&mut sim_state, 5.0, 1e-3, &step_leapfrog);

    let ke = sim_state.energy_kin();
    let pe = sim_state.energy_pot();
    println!(
        "total / kinetic / potential: {:.5e} / {:.5e} / {:.5e}",
        ke + pe,
        ke,
        pe
    );
}
