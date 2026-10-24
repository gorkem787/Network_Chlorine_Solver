//! Complementary error function, evaluated in log space so that products like
//! `exp(p) * erfc(z)` stay finite when `p` is huge and `erfc(z)` underflows.

/// Chebyshev coefficients for erfc on z >= 0 (Numerical Recipes, 3rd ed., §6.2.2).
const COF: [f64; 28] = [
    -1.3026537197817094,
    6.4196979235649026e-1,
    1.9476473204185836e-2,
    -9.561514786808631e-3,
    -9.46595344482036e-4,
    3.66839497852761e-4,
    4.2523324806907e-5,
    -2.0278578112534e-5,
    -1.624290004647e-6,
    1.303655835580e-6,
    1.5626441722e-8,
    -8.5238095915e-8,
    6.529054439e-9,
    5.059343495e-9,
    -9.91364156e-10,
    -2.27365122e-10,
    9.6467911e-11,
    2.394038e-12,
    -6.886027e-12,
    8.94487e-13,
    3.13092e-13,
    -1.12708e-13,
    3.81e-16,
    7.106e-15,
    -1.523e-15,
    -9.4e-17,
    1.21e-16,
    -2.8e-17,
];

/// ln(erfc(z)) for z >= 0, accurate even where erfc(z) underflows to 0.
fn ln_erfc_nonneg(z: f64) -> f64 {
    let t = 2.0 / (2.0 + z);
    let ty = 4.0 * t - 2.0;
    let (mut d, mut dd) = (0.0, 0.0);
    for &c in COF[1..].iter().rev() {
        let tmp = d;
        d = ty * d - dd + c;
        dd = tmp;
    }
    t.ln() - z * z + 0.5 * (COF[0] + ty * d) - dd
}

/// Complementary error function.
pub fn erfc(z: f64) -> f64 {
    if z >= 0.0 {
        ln_erfc_nonneg(z).exp()
    } else {
        2.0 - ln_erfc_nonneg(-z).exp()
    }
}

/// Computes `exp(p) * erfc(z)` without intermediate overflow or underflow.
pub fn exp_erfc(p: f64, z: f64) -> f64 {
    if z.is_infinite() {
        return if z > 0.0 { 0.0 } else { 2.0 * p.exp() };
    }
    if z >= 0.0 {
        (p + ln_erfc_nonneg(z)).exp()
    } else {
        p.exp() * erfc(z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rel(a: f64, b: f64) -> f64 {
        ((a - b) / b).abs()
    }

    #[test]
    fn erfc_matches_reference_values() {
        assert!(rel(erfc(0.0), 1.0) < 1e-14);
        assert!(rel(erfc(0.5), 0.4795001221869535) < 1e-14);
        assert!(rel(erfc(1.0), 0.15729920705028513) < 1e-14);
        assert!(rel(erfc(2.0), 0.004677734981047266) < 1e-13);
        assert!(rel(erfc(5.0), 1.5374597944280349e-12) < 1e-12);
        assert!(rel(erfc(-1.0), 1.8427007929497148) < 1e-14);
    }

    #[test]
    fn exp_erfc_survives_extremes() {
        // erfc(30) ~ 2.6e-393 underflows; exp(900) overflows. The product is ~1.
        let v = exp_erfc(900.0, 30.0);
        assert!(v.is_finite() && v > 0.0);
        let direct = (900.0 + ln_erfc_nonneg(30.0)).exp();
        assert!(rel(v, direct) < 1e-12);
        assert_eq!(exp_erfc(5.0, f64::INFINITY), 0.0);
    }
}
