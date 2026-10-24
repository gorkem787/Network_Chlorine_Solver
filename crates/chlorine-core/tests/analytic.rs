//! Accuracy tests for the analytic pipe scheme.

use chlorine_core::special::exp_erfc;
use chlorine_core::*;

const DIAMETER: f64 = 0.3;

fn area() -> f64 {
    std::f64::consts::PI * DIAMETER * DIAMETER / 4.0
}

fn single_pipe(length: f64, c0: f64) -> Network {
    Network {
        nodes: vec![
            Node { id: "R".into(), kind: NodeKind::Reservoir, initial_quality: c0, coords: None },
            Node { id: "J".into(), kind: NodeKind::Junction, initial_quality: c0, coords: None },
        ],
        links: vec![Link { id: "P".into(), kind: LinkKind::Pipe, from: 0, to: 1, length, diameter: DIAMETER }],
    }
}

/// Outlet series for a single pipe. `velocity(t)` is the velocity from t on,
/// `inlet(t)` the inlet concentration, both sampled on the grid.
fn run(
    length: f64,
    c0: f64,
    dt: f64,
    hours: f64,
    velocity: impl Fn(f64) -> f64,
    inlet: impl Fn(f64) -> f64,
    settings: Settings,
) -> Vec<f64> {
    let steps = (hours * 3600.0 / dt).round() as usize;
    let t = |k: usize| k as f64 * dt;
    let hyd = Hydraulics {
        dt,
        flows: (0..=steps).map(|k| vec![velocity(t(k)) * area()]).collect(),
        boundary: vec![Some((0..=steps).map(|k| inlet(t(k))).collect()), None],
    };
    simulate(&single_pipe(length, c0), &hyd, &settings).unwrap().node_quality[1].clone()
}

// ---- Exact solution, Eqs. 14-16 of the reference in the README, absolute times ----

fn a_term(x: f64, t: f64, ti: f64, m: f64, a: f64, b: f64, u: f64, d: f64, k: f64) -> f64 {
    if t <= ti {
        return 0.0;
    }
    let tau = t - ti;
    let w = (u * u + 4.0 * d * k).sqrt();
    let p = x / (2.0 * d) * (u + m * w);
    let z = (x + m * tau * w) / (2.0 * (d * tau).sqrt());
    exp_erfc(p, z) * (a * (tau / 2.0 + m * x / (2.0 * w)) + (a * ti + b) / 2.0)
}

fn b_term(x: f64, t: f64, m: f64, u: f64, d: f64, k: f64, c0: f64) -> f64 {
    let z = if x == 0.0 && m == 0.0 { 0.0 } else { (x + m * u * t) / (2.0 * (d * t).sqrt()) };
    c0 / 2.0 * exp_erfc(-k * t + (1.0 + m) * u * x / (2.0 * d), z)
}

/// Eq. 16 for a piecewise-linear source through `(times[i], values[i])` and initial `c0`.
fn exact(x: f64, t: f64, times: &[f64], values: &[f64], u: f64, d: f64, k: f64, c0: f64) -> f64 {
    let mut c = 2.0 * b_term(0.0, t, 0.0, u, d, k, c0) - b_term(x, t, -1.0, u, d, k, c0) - b_term(x, t, 1.0, u, d, k, c0);
    for i in 0..times.len() - 1 {
        let (t1, t2) = (times[i], times[i + 1]);
        let a = (values[i + 1] - values[i]) / (t2 - t1);
        let b = values[i] - a * t1;
        c += a_term(x, t, t1, -1.0, a, b, u, d, k) + a_term(x, t, t1, 1.0, a, b, u, d, k)
            - a_term(x, t, t2, -1.0, a, b, u, d, k)
            - a_term(x, t, t2, 1.0, a, b, u, d, k);
    }
    c
}

/// A non-trivial inlet pattern: a slug, a ramp and a dip.
fn pattern(t: f64) -> f64 {
    let h = t / 3600.0;
    if h < 1.0 {
        1.0
    } else if h < 2.0 {
        0.2
    } else if h < 3.0 {
        0.2 + 0.6 * (h - 2.0)
    } else {
        0.8
    }
}

/// `f` sampled every `step` seconds and linearly interpolated in between.
fn piecewise_linear_on(step: f64, t: f64, f: impl Fn(f64) -> f64) -> f64 {
    let t1 = (t / step).floor() * step;
    let t2 = t1 + step;
    f(t1) + (f(t2) - f(t1)) * (t - t1) / step
}

#[test]
fn constant_velocity_matches_paper_exactly() {
    let (length, u, k, c0, dt) = (1000.0, 0.2, 1.5 / 86400.0, 0.5, 900.0);
    let d = Dispersion::CHLORINE.coefficient(u, DIAMETER);
    let out = run(length, c0, dt, 8.0, |_| u, pattern, Settings { k, segment_length: 5.0, ..Default::default() });

    // The solver sees the inlet as piecewise linear between grid samples.
    let times: Vec<f64> = (0..out.len()).map(|i| i as f64 * dt).collect();
    let values: Vec<f64> = times.iter().map(|&t| pattern(t)).collect();
    let mut max: f64 = 0.0;
    for (i, &t) in times.iter().enumerate().skip(1) {
        let e = exact(length, t, &times[..=i], &values[..=i], u, d, k, c0);
        max = max.max((out[i] - e).abs());
        assert!((out[i] - e).abs() < 1e-9, "t={t}: solver {} vs exact {e}", out[i]);
    }
    eprintln!("constant velocity: max |solver - Eq.16| = {max:.2e}");
}

#[test]
fn result_does_not_depend_on_time_step_when_inputs_are_representable() {
    // Inlet piecewise linear and velocity changes on a 900 s grid: dt = 900,
    // 300 and 60 describe exactly the same problem and must agree. (The legacy
    // scheme gives visibly different answers here.)
    let velocity = |t: f64| if t < 7200.0 { 0.2 } else if t < 10800.0 { 0.05 } else { 0.15 };
    let inlet = |t: f64| piecewise_linear_on(900.0, t, pattern);
    let s = Settings { k: 1.0 / 86400.0, segment_length: 1.0, ..Default::default() };
    let coarse = run(600.0, 0.3, 900.0, 6.0, velocity, inlet, s);
    for dt in [300.0, 60.0] {
        let fine = run(600.0, 0.3, dt, 6.0, velocity, inlet, s);
        let ratio = (900.0 / dt) as usize;
        let max = coarse.iter().enumerate().map(|(i, c)| (c - fine[i * ratio]).abs()).fold(0.0, f64::max);
        eprintln!("dt 900 vs {dt}: max difference {max:.2e}");
        for (i, c) in coarse.iter().enumerate() {
            assert!((c - fine[i * ratio]).abs() < 2e-3, "dt={dt} t={}: {c} vs {}", i as f64 * 900.0, fine[i * ratio]);
        }
    }
}

/// Crank-Nicolson finite differences on a long fine grid, used as an
/// independent reference for time-varying velocity.
fn finite_difference(length: f64, c0: f64, hours: f64, velocity: impl Fn(f64) -> f64, inlet: impl Fn(f64) -> f64, k: f64) -> Vec<(f64, f64)> {
    let dx = 0.25;
    let dt = 1.0;
    let n = (3.0 * length / dx) as usize; // far boundary well past the outlet
    let mut c = vec![c0; n + 1];
    let mut out = vec![(0.0, c0)];
    let steps = (hours * 3600.0 / dt) as usize;
    let (mut a, mut b, mut cc, mut r) = (vec![0.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1], vec![0.0; n + 1]);
    let outlet = (length / dx) as usize;
    for s in 0..steps {
        let t = s as f64 * dt;
        let u = velocity(t);
        let d = Dispersion::CHLORINE.coefficient(u, DIAMETER);
        // C_t = -u C_x + d C_xx - k C, central differences, CN in time.
        let alpha = d * dt / (dx * dx);
        let beta = u * dt / (2.0 * dx);
        let lo = -alpha / 2.0 - beta / 2.0;
        let mid = 1.0 + alpha + k * dt / 2.0;
        let hi = -alpha / 2.0 + beta / 2.0;
        for i in 1..n {
            let explicit = (alpha / 2.0 + beta / 2.0) * c[i - 1] + (1.0 - alpha - k * dt / 2.0) * c[i] + (alpha / 2.0 - beta / 2.0) * c[i + 1];
            a[i] = lo;
            b[i] = mid;
            cc[i] = hi;
            r[i] = explicit;
        }
        // Dirichlet inlet, zero-gradient far end.
        let (a0, b0, c0_, r0) = (0.0, 1.0, 0.0, inlet(t + dt));
        a[0] = a0;
        b[0] = b0;
        cc[0] = c0_;
        r[0] = r0;
        a[n] = -1.0;
        b[n] = 1.0;
        cc[n] = 0.0;
        r[n] = 0.0;
        // Thomas algorithm.
        for i in 1..=n {
            let w = a[i] / b[i - 1];
            b[i] -= w * cc[i - 1];
            r[i] -= w * r[i - 1];
        }
        c[n] = r[n] / b[n];
        for i in (0..n).rev() {
            c[i] = (r[i] - cc[i] * c[i + 1]) / b[i];
        }
        out.push((t + dt, c[outlet]));
    }
    out
}

#[test]
fn variable_velocity_matches_finite_differences() {
    let (length, c0, hours, k) = (200.0, 0.3, 4.0, 1.0 / 86400.0);
    // Velocity drops, then rises; changes fall on the 900 s grid.
    let velocity = |t: f64| if t < 2700.0 { 0.12 } else if t < 6300.0 { 0.04 } else { 0.09 };
    let inlet = |t: f64| pattern(t * 2.0); // faster pattern so the pipe sees all of it
    // Make the inlet piecewise linear on the 900 s grid for both solvers.
    let dt = 900.0;
    let inlet_pl = |t: f64| piecewise_linear_on(dt, t, inlet);

    let analytic = run(length, c0, dt, hours, velocity, inlet, Settings { k, segment_length: 0.25, ..Default::default() });
    let fd = finite_difference(length, c0, hours, velocity, inlet_pl, k);
    let mut max_err: f64 = 0.0;
    for (i, c) in analytic.iter().enumerate() {
        let t = i as f64 * dt;
        let reference = fd[(t / 1.0) as usize].1;
        max_err = max_err.max((c - reference).abs());
    }
    eprintln!("variable velocity: max |analytic - FD| = {max_err:.2e}");
    assert!(max_err < 5e-3, "max |analytic - FD| = {max_err}");
}

#[test]
fn vanishing_dispersion_approaches_plug_flow() {
    let (length, u, k, dt) = (500.0, 0.25, 1.0 / 86400.0, 300.0);
    let tiny = Dispersion { coef: Dispersion::CHLORINE.coef * 1e-6, ..Dispersion::CHLORINE };
    let s = Settings { k, segment_length: 1.0, dispersion: tiny, ..Default::default() };
    let out = run(length, 0.0, dt, 4.0, |_| u, pattern, s);
    let travel = length / u;
    for (i, c) in out.iter().enumerate() {
        let t = i as f64 * dt;
        if t < travel + dt {
            continue; // skip the arrival of the front
        }
        // Plug flow: outlet = inlet (piecewise linear on the grid) one travel time ago, decayed.
        let expected = piecewise_linear_on(dt, t - travel, pattern) * (-k * travel).exp();
        assert!((c - expected).abs() < 1e-3, "t={t}: {c} vs {expected}");
    }
}

#[test]
fn flow_reversal_carries_profile_back() {
    // Fill a pipe from the left, stop, reverse: the water that comes back out
    // of the inlet end is the water that went in, minus decay.
    let (length, k) = (300.0, 0.0);
    let velocity = |t: f64| if t < 1800.0 { 0.1 } else { -0.1 };
    let net = single_pipe(length, 0.0);
    let dt = 300.0;
    let steps = 12;
    let flows = (0..=steps).map(|i| vec![velocity(i as f64 * dt) * area()]).collect();
    // R injects 1.0 while flowing forward; J is the source (0.0) when reversed.
    let mut hyd = Hydraulics { dt, flows, boundary: vec![Some(vec![1.0; steps + 1]), None] };
    hyd.boundary[1] = None;
    let r = simulate(&net, &hyd, &Settings { k, segment_length: 1.0, ..Default::default() }).unwrap();
    let j = &r.node_quality[1];
    assert!(j.iter().all(|c| c.is_finite() && *c >= 0.0 && *c <= 1.0 + 1e-9));
}
