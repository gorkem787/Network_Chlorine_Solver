//! Compares both pipe schemes with the exact analytical solution (README, References)
//! (Eq. 16: superposition over the whole source history) for a single pipe
//! with constant velocity.
//!
//! cargo run --release -p chlorine-core --example single_pipe_vs_exact

use chlorine_core::special::exp_erfc;
use chlorine_core::*;

/// Eq. 14 with absolute times; zero for t <= ti.
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

/// Eq. 16 for a source that is piecewise linear between `(times[i], values[i])`, C0 = 0.
fn exact(x: f64, t: f64, times: &[f64], values: &[f64], u: f64, d: f64, k: f64) -> f64 {
    let mut c = 0.0;
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

fn main() {
    let (length, diameter, u, k) = (1000.0, 0.3, 0.2, 1.0 / 86400.0);
    let q = u * std::f64::consts::PI * diameter * diameter / 4.0;
    let d = Dispersion::CHLORINE.coefficient(u, diameter);
    let hours = 6.0;
    // Source: 1 mg/L for the first hour, then 0 (a slug of chlorinated water).
    let source = |t: f64| if t < 3600.0 { 1.0 } else { 0.0 };

    println!("pipe {length} m, U = {u} m/s, D = {d:.3} m²/s, travel time {:.0} s\n", length / u);
    let net = Network {
        nodes: vec![
            Node { id: "R".into(), kind: NodeKind::Reservoir, initial_quality: 0.0, coords: None },
            Node { id: "J".into(), kind: NodeKind::Junction, initial_quality: 0.0, coords: None },
        ],
        links: vec![Link { id: "P".into(), kind: LinkKind::Pipe, from: 0, to: 1, length, diameter }],
    };

    let configs = [
        (Scheme::Legacy, 900.0, 1.0),
        (Scheme::Legacy, 60.0, 1.0),
        (Scheme::Analytic, 900.0, 1.0),
        (Scheme::Analytic, 60.0, 1.0),
    ];
    let mut runs = Vec::new();
    for &(scheme, dt, seg) in &configs {
        let steps = (hours * 3600.0 / dt) as usize;
        let series: Vec<f64> = (0..=steps).map(|k| source(k as f64 * dt)).collect();
        let hyd = Hydraulics { dt, flows: vec![vec![q]; steps + 1], boundary: vec![Some(series), None] };
        let r = simulate(&net, &hyd, &Settings { k, segment_length: seg, scheme, ..Default::default() }).unwrap();
        runs.push((dt, r.node_quality[1].clone()));
    }

    // The exact solution needs the source as piecewise linear; a 1 s ramp at 3600 s.
    let (st, sv) = ([0.0, 3600.0, 3601.0, 1e9], [1.0, 1.0, 0.0, 0.0]);
    print!("{:>6} {:>8}", "t (h)", "exact");
    for (scheme, dt, _) in configs {
        print!(" {:>14}", format!("{scheme:?} {dt}s"));
    }
    println!();
    for i in 0..=24 {
        let t = i as f64 * 900.0;
        print!("{:>6.2} {:>8.4}", t / 3600.0, exact(length, t, &st, &sv, u, d, k));
        for (dt, series) in &runs {
            let idx = (t / dt).round() as usize;
            print!(" {:>14.4}", series[idx]);
        }
        println!();
    }
}
