//! How much error do frequent period restarts add? A single pipe with constant
//! velocity, except that the flow is perturbed by a negligible amount every
//! step so that every step starts a new period. Compared with the unperturbed
//! (single-period, exact) run.
//!
//! cargo run --release -p chlorine-core --example restart_error

use chlorine_core::*;

fn main() {
    // A slow pipe: u*dt is below the cell length for dt = 60 s, the case
    // that froze water on a fixed grid.
    let (length, diameter, u, hours) = (300.0, 0.3, 0.01, 24.0);
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let net = Network {
        nodes: vec![
            Node { id: "R".into(), kind: NodeKind::Reservoir, initial_quality: 0.0, coords: None },
            Node { id: "J".into(), kind: NodeKind::Junction, initial_quality: 0.0, coords: None },
        ],
        links: vec![Link { id: "P".into(), kind: LinkKind::Pipe, from: 0, to: 1, length, diameter }],
    };
    let source = |t: f64| if (t / 3600.0) % 2.0 < 1.0 { 1.0 } else { 0.2 };

    println!("{:>6} {:>6} {:>8} {:>12}", "dt", "cell", "D scale", "max error");
    for dt in [900.0, 300.0, 60.0] {
        for cell in [1.0, 0.25] {
            for scale in [1.0, 1e-3] {
                let steps = (hours * 3600.0 / dt) as usize;
                let disp = Dispersion { coef: Dispersion::CHLORINE.coef * scale, ..Dispersion::CHLORINE };
                let s = Settings { k: 1.0 / 86400.0, segment_length: cell, dispersion: disp, ..Default::default() };
                let series: Vec<f64> = (0..=steps).map(|k| source(k as f64 * dt)).collect();
                let mut out = Vec::new();
                // 1e-5 is above the default period tolerance, so every step restarts.
                for perturb in [0.0, 1e-5] {
                    let flows = (0..=steps).map(|k| vec![u * area * (1.0 + perturb * (k % 2) as f64)]).collect();
                    let hyd = Hydraulics { dt, flows, boundary: vec![Some(series.clone()), None] };
                    out.push(simulate(&net, &hyd, &s).unwrap().node_quality[1].clone());
                }
                let err = out[0].iter().zip(&out[1]).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
                println!("{dt:>6} {cell:>6} {scale:>8} {err:>12.2e}");
            }
        }
    }
}
