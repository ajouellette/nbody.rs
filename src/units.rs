use crate::constants::*;

/// A system of units defined by a unit length, mass, and velocity.
/// Base units are assumed to be SI (m, kg, m/s).
#[derive(Debug, Clone, Copy)]
pub struct Units {
    length_unit: f64,
    mass_unit: f64,
    vel_unit: f64,
    grav_const: f64,
}

impl Units {
    /// Create a new `Units` instance with the given length, mass, and velocity units.
    pub fn new(length_unit: f64, mass_unit: f64, vel_unit: f64) -> Self {
        let grav_const = G_NEWTON * mass_unit / (length_unit * vel_unit.powi(2));
        Self {
            length_unit,
            mass_unit,
            vel_unit,
            grav_const,
        }
    }

    /// Create a new `Units` instance with SI (m, kg, m/s) units.
    pub fn new_si() -> Self {
        Self::new(1.0, 1.0, 1.0)
    }

    /// Create a new `Units` instance with kpc, M_sun, and km/s units.
    pub fn new_kpc_msun_km_s() -> Self {
        Self::new(MPC_IN_M / 1e3, M_SUN, 1e3)
    }

    /// Create a new `Units` instance with GADGET default units (kpc, 10^10 M_sun, km/s).
    pub fn new_gadget_default() -> Self {
        Self::new(MPC_IN_M / 1e3, 1e10 * M_SUN, 1e3)
    }

    /// Get the length unit.
    pub fn length_unit(&self) -> f64 {
        self.length_unit
    }

    /// Get the mass unit.
    pub fn mass_unit(&self) -> f64 {
        self.mass_unit
    }

    /// Get the velocity unit.
    pub fn vel_unit(&self) -> f64 {
        self.vel_unit
    }

    /// Get the time unit.
    pub fn time_unit(&self) -> f64 {
        self.length_unit / self.vel_unit
    }

    /// Get the gravitational constant.
    pub fn grav_const(&self) -> f64 {
        self.grav_const
    }

    // functions to convert key quantities
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_units_si() {
        let units = Units::new_si();
        assert_eq!(units.length_unit(), 1.0);
        assert_eq!(units.mass_unit(), 1.0);
        assert_eq!(units.vel_unit(), 1.0);
        assert_eq!(units.time_unit(), 1.0);
        assert_eq!(units.grav_const(), G_NEWTON);
    }

    #[test]
    fn test_units_gadget() {
        let units = Units::new_gadget_default();
        assert_eq!(units.length_unit(), MPC_IN_M / 1e3);
        assert_eq!(units.mass_unit(), 1e10 * M_SUN);
        assert_eq!(units.vel_unit(), 1e3);
        let time_unit_gadget = 9.8e8 * YR_IN_S;
        assert!((units.time_unit() - time_unit_gadget).abs() / time_unit_gadget < 1e-2);
        let grav_const_gadget = 43007.1;
        assert!((units.grav_const() - grav_const_gadget).abs() / grav_const_gadget < 1e-3);
    }
}
