//! Reads an EPANET input file, runs EPANET's hydraulic and water quality
//! simulation, and converts the results into `chlorine-core` inputs on a
//! uniform time grid.

mod project;
pub mod report;

use std::path::Path;

use chlorine_core::network::resample_step;
use chlorine_core::{Hydraulics, Link, LinkKind, Network, Node, NodeKind};
use epanet_sys as sys;

pub use project::{EpanetError, Project};

#[derive(Debug, Clone, Default)]
pub struct LoadOptions {
    /// Simulation time step, s.
    pub dt: f64,
    /// Overrides the input file's duration, s.
    pub duration: Option<f64>,
    /// When set (1/s), EPANET's own reaction coefficients are replaced with this
    /// first-order bulk decay (and wall reactions disabled), so that EPANET's
    /// results are a like-for-like reference.
    pub sync_decay: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct Loaded {
    pub network: Network,
    pub hydraulics: Hydraulics,
    /// EPANET's own water quality result, `reference[node][k]`, mg/L.
    pub reference: Vec<Vec<f64>>,
    pub flow_units: &'static str,
    pub warnings: Vec<String>,
}

struct Units {
    name: &'static str,
    /// m³/s per flow unit.
    flow: f64,
    /// m per length unit.
    length: f64,
    /// m per diameter unit.
    diameter: f64,
}

fn units(code: i32) -> Result<Units, EpanetError> {
    let us = |name, flow| Units { name, flow, length: 0.3048, diameter: 0.0254 };
    let si = |name, flow| Units { name, flow, length: 1.0, diameter: 0.001 };
    Ok(match code {
        sys::EN_CFS => us("CFS", 0.028316846592),
        sys::EN_GPM => us("GPM", 6.30901964e-5),
        sys::EN_MGD => us("MGD", 0.0438126364),
        sys::EN_IMGD => us("IMGD", 0.0526168042),
        sys::EN_AFD => us("AFD", 0.0142764101),
        sys::EN_LPS => si("LPS", 1e-3),
        sys::EN_LPM => si("LPM", 1e-3 / 60.0),
        sys::EN_MLD => si("MLD", 1e3 / 86400.0),
        sys::EN_CMH => si("CMH", 1.0 / 3600.0),
        sys::EN_CMD => si("CMD", 1.0 / 86400.0),
        other => return Err(EpanetError::other(format!("unknown flow unit code {other}"))),
    })
}

pub fn load(path: &Path, opts: &LoadOptions) -> Result<Loaded, EpanetError> {
    if !(opts.dt > 0.0) {
        return Err(EpanetError::other("time step must be positive"));
    }
    let mut p = Project::open(path)?;
    let mut warnings = Vec::new();
    let u = units(p.flow_units()?)?;
    let n_nodes = p.count(sys::EN_NODECOUNT)?;
    let n_links = p.count(sys::EN_LINKCOUNT)?;

    let mut nodes = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        let kind = match p.node_type(i)? {
            sys::EN_RESERVOIR => NodeKind::Reservoir,
            sys::EN_TANK => NodeKind::Tank,
            _ => NodeKind::Junction,
        };
        nodes.push(Node {
            id: p.node_id(i)?,
            kind,
            initial_quality: p.node_value(i, sys::EN_INITQUAL)?,
            coords: p.node_coords(i),
        });
    }

    let mut links = Vec::with_capacity(n_links);
    for i in 0..n_links {
        let kind = match p.link_type(i)? {
            sys::EN_CVPIPE | sys::EN_PIPE => LinkKind::Pipe,
            sys::EN_PUMP => LinkKind::Pump,
            _ => LinkKind::Valve,
        };
        let (from, to) = p.link_nodes(i)?;
        let (length, diameter) = if kind == LinkKind::Pipe {
            (p.link_value(i, sys::EN_LENGTH)? * u.length, p.link_value(i, sys::EN_DIAMETER)? * u.diameter)
        } else {
            (0.0, 0.0)
        };
        links.push(Link { id: p.link_id(i)?, kind, from, to, length, diameter });
    }

    let duration = match opts.duration {
        Some(d) => {
            p.set_time_param(sys::EN_DURATION, d.round() as i64)?;
            d
        }
        None => p.time_param(sys::EN_DURATION)? as f64,
    };
    if duration < opts.dt {
        return Err(EpanetError::other(format!(
            "simulation duration ({duration} s) is shorter than the time step ({} s); \
             the network needs an extended-period [TIMES] Duration",
            opts.dt
        )));
    }
    let steps = (duration / opts.dt + 1e-9).floor() as usize;

    // Quality keeps changing between hydraulic events, so EPANET must stop at
    // every grid time; it limits its hydraulic step to the next reporting time.
    // (The Python version did this with setTimeReportingStep.)
    let dt_whole = opts.dt.round();
    if (opts.dt - dt_whole).abs() > 1e-9 || dt_whole < 1.0 {
        return Err(EpanetError::other(format!("time step must be a whole number of seconds, got {}", opts.dt)));
    }
    p.set_time_param(sys::EN_REPORTSTEP, dt_whole as i64)?;
    p.set_time_param(sys::EN_REPORTSTART, 0)?;
    if p.time_param(sys::EN_HYDSTEP)? > dt_whole as i64 {
        p.set_time_param(sys::EN_HYDSTEP, dt_whole as i64)?;
    }

    if p.qual_type()? != sys::EN_CHEM {
        warnings.push("input file has no chemical quality analysis; enabled one as 'Chlorine' (mg/L)".into());
        p.set_chemical("Chlorine", "mg/L")?;
    }
    if let Some(k) = opts.sync_decay {
        let per_day = -k * 86400.0;
        for (i, link) in links.iter().enumerate() {
            if link.kind == LinkKind::Pipe {
                p.set_link_value(i, sys::EN_KBULK, per_day)?;
                p.set_link_value(i, sys::EN_KWALL, 0.0)?;
            }
        }
        for (i, node) in nodes.iter().enumerate() {
            if node.kind == NodeKind::Tank {
                p.set_node_value(i, sys::EN_TANK_KBULK, per_day)?;
            }
        }
    }

    // EPANET returns error 240 for nodes without a [SOURCES] entry.
    let source_strength: Vec<f64> = (0..n_nodes).map(|i| p.node_value(i, sys::EN_SOURCEQUAL).unwrap_or(0.0)).collect();
    let (mut times, mut flows, mut quality) = (Vec::new(), Vec::new(), Vec::new());
    let mut injects = vec![false; n_nodes];
    p.run(|p, t| {
        times.push(t as f64);
        flows.push((0..n_links).map(|i| Ok(p.link_value(i, sys::EN_FLOW)? * u.flow)).collect::<Result<Vec<_>, EpanetError>>()?);
        quality.push((0..n_nodes).map(|i| p.node_value(i, sys::EN_QUALITY)).collect::<Result<Vec<_>, _>>()?);
        for (i, inj) in injects.iter_mut().enumerate() {
            // Negative demand = water entering the network (the Python version's
            // way of modelling a pump inflow).
            *inj |= p.node_value(i, sys::EN_DEMAND)? < 0.0;
        }
        Ok(())
    })?;

    let missing = (0..=steps).filter(|k| !times.iter().any(|t| (t - *k as f64 * opts.dt).abs() < 0.5)).count();
    if missing > 0 {
        warnings.push(format!("EPANET did not stop at {missing} grid times; their values are held from the previous event"));
    }
    let flows = resample_step(&times, &flows, opts.dt, steps);
    let quality = resample_step(&times, &quality, opts.dt, steps);
    let reference: Vec<Vec<f64>> = (0..n_nodes).map(|n| quality.iter().map(|row| row[n]).collect()).collect();

    let boundary = nodes
        .iter()
        .enumerate()
        .map(|(n, node)| {
            let is_source = node.kind != NodeKind::Junction || injects[n] || source_strength[n] > 0.0;
            is_source.then(|| reference[n].clone())
        })
        .collect();

    Ok(Loaded {
        network: Network { nodes, links },
        hydraulics: Hydraulics { dt: opts.dt, flows, boundary },
        reference,
        flow_units: u.name,
        warnings,
    })
}
