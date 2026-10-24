//! Analytical solution of the 1-D advection-dispersion-reaction equation over a
//! single pipe segment and a single time step.
//!
//! This is the same closed form as the original `ADE_solver.py`: a source at
//! x = 0 whose concentration varies linearly from `w1` to `w2` during the step,
//! an initial concentration `c0` in the segment, and first-order bulk decay `K`.
//! The difference is that every `exp(..) * erfc(..)` product is evaluated in
//! log space, so low velocities and long segments no longer overflow into NaN.

use crate::special::exp_erfc;

/// Empirical dispersion model `D = coef * Re^exponent`, with `Re = U d / nu`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dispersion {
    /// Kinematic viscosity of water, m²/s.
    pub nu: f64,
    pub coef: f64,
    pub exponent: f64,
}

impl Dispersion {
    pub const CHLORINE: Dispersion = Dispersion { nu: 1.26e-6, coef: 1e-4, exponent: 0.66 };
    pub const FLUORINE: Dispersion = Dispersion { nu: 1.0533e-6, coef: 2e-5, exponent: 0.83 };

    /// Dispersion coefficient (m²/s) for velocity `u` (m/s) in a pipe of diameter `d` (m).
    pub fn coefficient(&self, u: f64, d: f64) -> f64 {
        self.coef * (u.abs() * d / self.nu).powf(self.exponent)
    }
}

impl Default for Dispersion {
    fn default() -> Self {
        Self::CHLORINE
    }
}

/// Physical parameters of one segment for one time step.
#[derive(Debug, Clone, Copy)]
pub struct Segment {
    /// Flow velocity, m/s (> 0).
    pub u: f64,
    /// Dispersion coefficient, m²/s (> 0).
    pub d: f64,
    /// First-order bulk decay coefficient, 1/s.
    pub k: f64,
    /// Segment length, m.
    pub x: f64,
    /// Time step, s.
    pub dt: f64,
}

/// Result of one segment step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepResult {
    pub value: f64,
    /// True when the closed form was not finite and plug flow was used instead.
    pub fell_back: bool,
}

impl Segment {
    /// Concentration at the downstream end of the segment at the end of the step.
    ///
    /// `w1`, `w2`: source (upstream) concentration at the start and end of the step.
    /// `c0`: initial concentration in the segment.
    pub fn step(&self, w1: f64, w2: f64, c0: f64) -> StepResult {
        if self.x <= 0.0 {
            return StepResult { value: w2, fell_back: false };
        }
        let value = self.closed_form(w1, w2, c0);
        if value.is_finite() {
            StepResult { value, fell_back: false }
        } else {
            StepResult { value: self.plug_flow(w1, w2, c0), fell_back: true }
        }
    }

    fn closed_form(&self, w1: f64, w2: f64, c0: f64) -> f64 {
        let Segment { u, d, k, x, dt: tau } = *self;
        let w = (u * u + 4.0 * d * k).sqrt();
        let a = (w2 - w1) / tau;
        let sqrt_dt = 2.0 * (d * tau).sqrt();

        // Source terms. The original also subtracts the same terms evaluated at
        // the end of the step; for x > 0 those are exactly zero.
        let source = |m: f64| {
            let p = x / (2.0 * d) * (u + m * w);
            let z = (x + m * tau * w) / sqrt_dt;
            let ramp = a * (tau / 2.0 + m * x / (2.0 * w)) + w1 / 2.0;
            exp_erfc(p, z) * ramp
        };
        // Initial-condition terms.
        let initial = |m: f64| {
            let p = -k * tau + (1.0 + m) * u * x / (2.0 * d);
            let z = (x + m * u * tau) / sqrt_dt;
            c0 / 2.0 * exp_erfc(p, z)
        };

        source(-1.0) + source(1.0) + c0 * (-k * tau).exp() - initial(-1.0) - initial(1.0)
    }

    /// Pure advection with decay; used only when the closed form breaks down.
    fn plug_flow(&self, w1: f64, w2: f64, c0: f64) -> f64 {
        let travel = self.x / self.u;
        if travel <= self.dt {
            let entered_at = self.dt - travel;
            (w1 + (w2 - w1) * entered_at / self.dt) * (-self.k * travel).exp()
        } else {
            c0 * (-self.k * self.dt).exp()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(u: f64, x: f64, k: f64, dt: f64) -> Segment {
        Segment { u, d: Dispersion::CHLORINE.coefficient(u, 0.3), k, x, dt }
    }

    #[test]
    fn constant_source_without_decay_is_preserved() {
        // Fluid fully replaced within the step: outlet equals the source.
        let r = seg(1.0, 1.0, 0.0, 900.0).step(0.8, 0.8, 0.8);
        assert!((r.value - 0.8).abs() < 1e-9, "{r:?}");
        assert!(!r.fell_back);
    }

    #[test]
    fn decay_matches_first_order_kinetics() {
        let k = 1.5 / 86400.0;
        let s = seg(0.5, 1.0, k, 900.0);
        let r = s.step(1.0, 1.0, 1.0);
        // Uniform field: everything decays by exp(-k * dt) at most; advection
        // from a fresh source makes the segment end only slightly older.
        let lo = (-k * 900.0).exp();
        assert!(r.value >= lo - 1e-9 && r.value <= 1.0 + 1e-9, "{r:?}");
    }

    #[test]
    fn slow_flow_long_segment_stays_finite() {
        // The Python version produced inf * 0 = NaN here and silently returned 0.
        let s = seg(1e-4, 50.0, 1e-5, 900.0);
        let r = s.step(1.0, 0.5, 0.7);
        assert!(r.value.is_finite());
        // Front has not arrived: outlet is the decayed initial content.
        assert!((r.value - 0.7 * (-1e-5 * 900.0f64).exp()).abs() < 1e-6, "{r:?}");
    }

    #[test]
    fn matches_original_python_solver() {
        // (u, diameter, x, k, dt, w1, w2, c0, ADE_solver.solver output)
        let cases = [
            (0.5, 0.3, 1.0, 1.736111111111111e-05, 900.0, 1.0, 0.6, 0.8, 0.6008679977431729),
            (0.05, 0.15, 2.0, 1.1574074074074073e-05, 900.0, 0.4, 0.9, 0.2, 0.8773779233817107),
            (1.2, 0.5, 10.0, 3.472222222222222e-05, 600.0, 0.0, 1.0, 0.5, 0.9858261973824519),
            (0.01, 0.1, 5.0, 0.0, 900.0, 1.0, 1.0, 0.0, 0.9107904698529736),
        ];
        for (u, diameter, x, k, dt, w1, w2, c0, expected) in cases {
            let s = Segment { u, d: Dispersion::CHLORINE.coefficient(u, diameter), k, x, dt };
            let r = s.step(w1, w2, c0);
            assert!((r.value - expected).abs() < 1e-12, "{r:?} vs {expected}");
        }
    }

    #[test]
    fn source_ramp_reaches_outlet_with_delay() {
        // Short travel time: outlet sees close to the end-of-step source value.
        let r = seg(1.0, 10.0, 0.0, 900.0).step(0.0, 1.0, 0.0);
        let expected = (900.0 - 10.0) / 900.0;
        assert!((r.value - expected).abs() < 0.01, "{r:?} vs {expected}");
    }
}
