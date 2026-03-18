use crate::SimState;

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
