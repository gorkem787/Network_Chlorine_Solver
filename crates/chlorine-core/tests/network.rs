use chlorine_core::*;

fn node(id: &str, kind: NodeKind, c0: f64) -> Node {
    Node { id: id.into(), kind, initial_quality: c0, coords: None }
}

fn pipe(id: &str, from: usize, to: usize, length: f64, diameter: f64) -> Link {
    Link { id: id.into(), kind: LinkKind::Pipe, from, to, length, diameter }
}

/// Constant flows and constant boundary values for `steps` steps.
fn steady(net: &Network, flows: &[f64], sources: &[(usize, f64)], dt: f64, steps: usize) -> Hydraulics {
    let mut boundary = vec![None; net.nodes.len()];
    for &(n, c) in sources {
        boundary[n] = Some(vec![c; steps + 1]);
    }
    Hydraulics { dt, flows: vec![flows.to_vec(); steps + 1], boundary }
}

fn last(r: &Results, node: usize) -> f64 {
    *r.node_quality[node].last().unwrap()
}

#[test]
fn single_pipe_reaches_source_concentration_without_decay() {
    let net = Network {
        nodes: vec![node("R", NodeKind::Reservoir, 1.0), node("J", NodeKind::Junction, 0.0)],
        links: vec![pipe("P", 0, 1, 500.0, 0.3)],
    };
    let hyd = steady(&net, &[0.035], &[(0, 1.0)], 900.0, 20);
    let s = Settings { k: 0.0, segment_length: 10.0, ..Default::default() };
    let r = simulate(&net, &hyd, &s).unwrap();
    assert!((last(&r, 1) - 1.0).abs() < 1e-3, "{}", last(&r, 1));
    assert_eq!(r.diagnostics, Diagnostics::default());
}

#[test]
fn steady_state_decay_follows_travel_time() {
    let (length, d, q) = (1000.0, 0.3, 0.035);
    let net = Network {
        nodes: vec![node("R", NodeKind::Reservoir, 1.0), node("J", NodeKind::Junction, 1.0)],
        links: vec![pipe("P", 0, 1, length, d)],
    };
    let k = 2.0 / 86400.0;
    let hyd = steady(&net, &[q], &[(0, 1.0)], 900.0, 40);
    let r = simulate(&net, &hyd, &Settings { k, segment_length: 5.0, ..Default::default() }).unwrap();
    let u = q / (std::f64::consts::PI * d * d / 4.0);
    let expected = (-k * length / u).exp();
    assert!((last(&r, 1) - expected).abs() / expected < 0.02, "{} vs {expected}", last(&r, 1));
}

#[test]
fn junction_mixes_by_flow() {
    let net = Network {
        nodes: vec![
            node("A", NodeKind::Reservoir, 1.0),
            node("B", NodeKind::Reservoir, 0.2),
            node("J", NodeKind::Junction, 0.0),
        ],
        links: vec![pipe("PA", 0, 2, 50.0, 0.2), pipe("PB", 1, 2, 50.0, 0.2)],
    };
    let hyd = steady(&net, &[0.01, 0.03], &[(0, 1.0), (1, 0.2)], 900.0, 10);
    let r = simulate(&net, &hyd, &Settings { k: 0.0, segment_length: 5.0, ..Default::default() }).unwrap();
    let expected = (0.01 * 1.0 + 0.03 * 0.2) / 0.04;
    assert!((last(&r, 2) - expected).abs() < 1e-3, "{}", last(&r, 2));
}

#[test]
fn stagnant_node_decays_exponentially() {
    let net = Network {
        nodes: vec![node("R", NodeKind::Reservoir, 1.0), node("J", NodeKind::Junction, 0.8)],
        links: vec![pipe("P", 0, 1, 100.0, 0.2)],
    };
    let k = 1.0 / 86400.0;
    let steps = 24;
    let hyd = steady(&net, &[0.0], &[(0, 1.0)], 3600.0, steps);
    let r = simulate(&net, &hyd, &Settings { k, ..Default::default() }).unwrap();
    let expected = 0.8 * (-k * 3600.0 * steps as f64).exp();
    assert!((last(&r, 1) - expected).abs() < 1e-12);
}

#[test]
fn flow_reversal_keeps_results_finite() {
    let net = Network {
        nodes: vec![
            node("A", NodeKind::Reservoir, 1.0),
            node("J", NodeKind::Junction, 0.5),
            node("B", NodeKind::Tank, 0.3),
        ],
        links: vec![pipe("P1", 0, 1, 200.0, 0.2), pipe("P2", 1, 2, 200.0, 0.2)],
    };
    let steps = 48;
    let flows: Vec<Vec<f64>> = (0..=steps)
        .map(|k| {
            // P2 alternates direction every 6 hours: tank fills, then drains.
            let q2 = if (k / 24) % 2 == 0 { 0.01 } else { -0.01 };
            vec![0.02, q2]
        })
        .collect();
    let mut boundary = vec![None; 3];
    boundary[0] = Some(vec![1.0; steps + 1]);
    boundary[2] = Some(vec![0.3; steps + 1]);
    let hyd = Hydraulics { dt: 900.0, flows, boundary };
    let r = simulate(&net, &hyd, &Settings { segment_length: 10.0, ..Default::default() }).unwrap();
    assert!(r.node_quality.iter().flatten().all(|c| c.is_finite() && *c >= 0.0));
    // While the tank drains into J, J must lie between the two source values.
    let j = &r.node_quality[1];
    assert!(j[steps] > 0.3 - 1e-9 && j[steps] < 1.0 + 1e-9, "{}", j[steps]);
}

#[test]
fn pump_loop_terminates_and_reports_cycle() {
    // R -> J1 -> J2 -> J3 -> (pump) -> J1 forms a directed loop.
    let mut links = vec![
        pipe("P0", 0, 1, 100.0, 0.2),
        pipe("P1", 1, 2, 100.0, 0.2),
        pipe("P2", 2, 3, 100.0, 0.2),
    ];
    links.push(Link { id: "PU".into(), kind: LinkKind::Pump, from: 3, to: 1, length: 0.0, diameter: 0.0 });
    let net = Network {
        nodes: vec![
            node("R", NodeKind::Reservoir, 1.0),
            node("J1", NodeKind::Junction, 0.0),
            node("J2", NodeKind::Junction, 0.0),
            node("J3", NodeKind::Junction, 0.0),
        ],
        links,
    };
    let hyd = steady(&net, &[0.01, 0.02, 0.02, 0.01], &[(0, 1.0)], 900.0, 30);
    let r = simulate(&net, &hyd, &Settings { k: 0.0, segment_length: 10.0, ..Default::default() }).unwrap();
    assert!(r.diagnostics.cycle_breaks > 0);
    // With no decay everything converges to the source concentration.
    for n in 1..4 {
        assert!((last(&r, n) - 1.0).abs() < 0.05, "node {n}: {}", last(&r, n));
    }
}

#[test]
fn legacy_scheme_still_runs() {
    let net = Network {
        nodes: vec![node("R", NodeKind::Reservoir, 1.0), node("J", NodeKind::Junction, 0.0)],
        links: vec![pipe("P", 0, 1, 500.0, 0.3)],
    };
    let hyd = steady(&net, &[0.035], &[(0, 1.0)], 900.0, 20);
    let s = Settings { k: 0.0, segment_length: 10.0, scheme: Scheme::Legacy, ..Default::default() };
    let r = simulate(&net, &hyd, &s).unwrap();
    assert!((last(&r, 1) - 1.0).abs() < 1e-3, "{}", last(&r, 1));
}

#[test]
fn rejects_network_without_source() {
    let net = Network {
        nodes: vec![node("A", NodeKind::Junction, 1.0), node("B", NodeKind::Junction, 0.0)],
        links: vec![pipe("P", 0, 1, 10.0, 0.1)],
    };
    let hyd = steady(&net, &[0.01], &[], 900.0, 2);
    assert!(matches!(simulate(&net, &hyd, &Settings::default()), Err(SimError::InvalidInput(_))));
}

#[test]
fn progress_callback_can_cancel() {
    let net = Network {
        nodes: vec![node("R", NodeKind::Reservoir, 1.0), node("J", NodeKind::Junction, 0.0)],
        links: vec![pipe("P", 0, 1, 100.0, 0.2)],
    };
    let hyd = steady(&net, &[0.01], &[(0, 1.0)], 900.0, 10);
    let r = simulate_with_progress(&net, &hyd, &Settings::default(), |k, _| k < 3);
    assert_eq!(r.unwrap_err(), SimError::Cancelled);
}
