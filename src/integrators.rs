use crate::SimState;

/// Compute one simulation step using the symmetric Euler method.
///
/// $x^{n+1} = x^n + v^n \Delta t$
/// $v^{n+1} = v^n + a^{n+1} \Delta t$
/// This is a first-order, time-symmetric method.
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

/// Compute one simulation step using the leapfrog method.
///
/// 1/2-Kick: $v^{n+1/2} = v^n + a^n \Delta t/2$
/// Drift: $x^{n+1} = x^n + v^{n+1/2} \Delta t$
/// 1/2-Kick: $v^{n+1} = v^{n+1/2} + a^{n+1} \Delta t/2$
/// This is a second-order, time-symmetric method.
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

/// Time-integrate the system for given total time.
pub fn integrate<const NDIM: usize>(
    sim_state: &mut SimState<NDIM>,
    total_time: f64,
    dt: f64,
    stepper_fn: impl Fn(&mut SimState<NDIM>, f64),
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
