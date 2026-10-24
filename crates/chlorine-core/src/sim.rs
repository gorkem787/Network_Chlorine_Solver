//! Network-wide water quality simulation.
//!
//! At every time step nodes are visited in flow order: a node's concentration
//! is the flow-weighted mix of what its inflowing links deliver, and each
//! outflowing pipe then carries that concentration to its far end using the
//! selected [`Scheme`].

use std::collections::VecDeque;
use std::fmt;

use crate::ade::{Dispersion, Segment};
use crate::analytic::AnalyticPipe;
use crate::network::{Hydraulics, Network};

/// How concentration is carried through a pipe.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Scheme {
    /// The analytical solution's Eq. 16 (README, References) applied to the whole pipe, restarted from
    /// the actual in-pipe profile only when the velocity changes.
    #[default]
    Analytic,
    /// The original Python scheme: Eq. 13 restarted every step on every
    /// segment, with the upstream point's previous value as the initial
    /// concentration. Kept to reproduce earlier results; it is not convergent.
    Legacy,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Settings {
    /// First-order bulk decay coefficient, 1/s.
    pub k: f64,
    /// Target segment length, m: the profile cell size for [`Scheme::Analytic`],
    /// the chained segment length for [`Scheme::Legacy`].
    pub segment_length: f64,
    pub dispersion: Dispersion,
    /// Velocities at or below this (m/s) are treated as stagnant water.
    pub stagnant_velocity: f64,
    pub scheme: Scheme,
    /// [`Scheme::Analytic`]: relative flow change below which a pipe keeps its
    /// current constant-velocity period. Absorbs the tiny differences EPANET
    /// produces when it re-solves unchanged hydraulics at every report time.
    pub period_tolerance: f64,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            k: 1.0 / 86400.0,
            segment_length: 1.0,
            dispersion: Dispersion::CHLORINE,
            stagnant_velocity: 1e-9,
            scheme: Scheme::Analytic,
            period_tolerance: 1e-4,
        }
    }
}

/// Things the solver had to work around. All zero for a well-behaved run.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Diagnostics {
    /// Steps × nodes where a directed flow cycle had to be broken by using the
    /// previous step's value for one inflow.
    pub cycle_breaks: usize,
    /// Segment steps where the closed form was not finite and plug flow was used.
    pub plug_flow_fallbacks: usize,
    /// Segment steps that produced a (slightly) negative concentration, clamped to 0.
    pub negative_clamped: usize,
    /// [`Scheme::Analytic`]: pipe periods restarted because the flow changed.
    pub period_restarts: usize,
}

#[derive(Debug, Clone)]
pub struct Results {
    pub times: Vec<f64>,
    /// `node_quality[node][k]`, mg/L.
    pub node_quality: Vec<Vec<f64>>,
    pub diagnostics: Diagnostics,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SimError {
    InvalidInput(String),
    Cancelled,
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SimError::InvalidInput(msg) => write!(f, "invalid input: {msg}"),
            SimError::Cancelled => write!(f, "simulation cancelled"),
        }
    }
}

impl std::error::Error for SimError {}

pub fn simulate(net: &Network, hyd: &Hydraulics, settings: &Settings) -> Result<Results, SimError> {
    simulate_with_progress(net, hyd, settings, |_, _| true)
}

/// Like [`simulate`], calling `progress(step, total_steps)` after every step.
/// Returning `false` from `progress` cancels the run.
pub fn simulate_with_progress(
    net: &Network,
    hyd: &Hydraulics,
    settings: &Settings,
    mut progress: impl FnMut(usize, usize) -> bool,
) -> Result<Results, SimError> {
    validate(net, hyd, settings)?;

    let steps = hyd.steps();
    let dt = hyd.dt;
    let decay = (-settings.k * dt).exp();
    let n_nodes = net.nodes.len();
    let n_links = net.links.len();

    let mut node_quality: Vec<Vec<f64>> = net
        .nodes
        .iter()
        .zip(&hyd.boundary)
        .map(|(node, b)| {
            let mut series = Vec::with_capacity(steps + 1);
            series.push(b.as_ref().map_or(node.initial_quality, |b| b[0]));
            series
        })
        .collect();

    let mut pipes: Vec<Transport> = net
        .links
        .iter()
        .map(|link| {
            if link.is_pass_through() {
                return Transport::PassThrough;
            }
            match settings.scheme {
                // EPANET fills each pipe with its end node's initial quality.
                Scheme::Analytic => {
                    let c0 = net.nodes[link.to].initial_quality;
                    Transport::Analytic(AnalyticPipe::new(link.length, settings.segment_length, c0))
                }
                // The Python version used the start node.
                Scheme::Legacy => {
                    let c0 = net.nodes[link.from].initial_quality;
                    let n = ((link.length / settings.segment_length).round() as usize).max(1);
                    Transport::Legacy(PipeState { segment: link.length / n as f64, cur: vec![c0; n + 1], prev: vec![c0; n + 1] })
                }
            }
        })
        .collect();

    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (i, link) in net.links.iter().enumerate() {
        incident[link.from].push(i);
        if link.to != link.from {
            incident[link.to].push(i);
        }
    }

    let mut diag = Diagnostics::default();
    let mut outlet = vec![0.0; n_links];
    let mut prev_outlet = vec![0.0; n_links];
    let mut done_link = vec![false; n_links];
    let mut upstream = vec![usize::MAX; n_links]; // usize::MAX = inactive this step
    let mut pending = vec![0usize; n_nodes];
    let mut resolved = vec![false; n_nodes];
    let mut queue = VecDeque::with_capacity(n_nodes);
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get()).min(n_links.max(1));

    for k in 1..=steps {
        // The flow during [t_{k-1}, t_k] is the one EPANET reports at t_{k-1}.
        // The legacy scheme used the value at t_k, as the Python code did.
        let flows = match settings.scheme {
            Scheme::Analytic => &hyd.flows[k - 1],
            Scheme::Legacy => &hyd.flows[k],
        };

        let mut restarting = 0;
        for (i, link) in net.links.iter().enumerate() {
            let q = flows[i];
            upstream[i] = match (is_active(link, q, settings), q > 0.0) {
                (false, _) => usize::MAX,
                (true, true) => link.from,
                (true, false) => link.to,
            };
            done_link[i] = false;
            if let Transport::Analytic(pipe) = &pipes[i] {
                restarting += pipe.will_restart(q, settings.period_tolerance) as usize;
            }
        }

        // Pipes are independent here. Restarting a period re-evaluates the
        // whole in-pipe profile, which dominates the run time, so spread
        // those steps over threads when there are several.
        let prepare = |pipes: &mut [Transport], links: &[crate::network::Link], flows: &[f64]| -> usize {
            pipes.iter_mut().zip(links).zip(flows).map(|((p, l), &q)| p.prepare(l, q, settings, dt, decay)).sum()
        };
        diag.period_restarts += if restarting >= 2 && threads > 1 {
            let chunk = n_links.div_ceil(threads);
            std::thread::scope(|scope| {
                let handles: Vec<_> = pipes
                    .chunks_mut(chunk)
                    .zip(net.links.chunks(chunk))
                    .zip(flows.chunks(chunk))
                    .map(|((p, l), f)| scope.spawn(move || prepare(p, l, f)))
                    .collect();
                handles.into_iter().map(|h| h.join().expect("pipe worker panicked")).sum::<usize>()
            })
        } else {
            prepare(&mut pipes, &net.links, flows)
        };

        let downstream = |i: usize, up: &[usize]| {
            let l = &net.links[i];
            if up[i] == l.from { l.to } else { l.from }
        };

        queue.clear();
        for n in 0..n_nodes {
            resolved[n] = false;
            pending[n] = incident[n]
                .iter()
                .filter(|&&i| upstream[i] != usize::MAX && downstream(i, &upstream) == n)
                .count();
            if hyd.boundary[n].is_some() || pending[n] == 0 {
                queue.push_back(n);
                resolved[n] = true;
            }
        }

        let mut n_resolved = 0;
        while n_resolved < n_nodes {
            let node = match queue.pop_front() {
                Some(n) => n,
                None => {
                    // Directed flow cycle (only possible through a pump). Break
                    // it at the node with the fewest unresolved inflows.
                    let n = (0..n_nodes)
                        .filter(|&n| !resolved[n])
                        .min_by_key(|&n| pending[n])
                        .expect("unresolved node must exist");
                    resolved[n] = true;
                    diag.cycle_breaks += 1;
                    n
                }
            };
            n_resolved += 1;

            let value = match &hyd.boundary[node] {
                Some(series) => series[k],
                None => {
                    let (mut mass, mut flow) = (0.0, 0.0);
                    for &i in &incident[node] {
                        if upstream[i] != usize::MAX && downstream(i, &upstream) == node {
                            let q = flows[i].abs();
                            let c = if done_link[i] { outlet[i] } else { prev_outlet[i] };
                            mass += q * c;
                            flow += q;
                        }
                    }
                    if flow > 0.0 { mass / flow } else { node_quality[node][k - 1] * decay }
                }
            };
            node_quality[node].push(value);

            for &i in &incident[node] {
                if upstream[i] != node {
                    continue;
                }
                let link = &net.links[i];
                let before = node_quality[node][k - 1];
                outlet[i] = match &mut pipes[i] {
                    Transport::PassThrough => value,
                    Transport::Analytic(pipe) => pipe.advance(before, value, settings.k, dt).max(0.0),
                    Transport::Legacy(pipe) => {
                        let u = flows[i].abs() / link.area();
                        let seg = Segment {
                            u,
                            d: settings.dispersion.coefficient(u, link.diameter),
                            k: settings.k,
                            x: pipe.segment,
                            dt,
                        };
                        pipe.advance(&seg, flows[i] > 0.0, before, value, &mut diag)
                    }
                };
                done_link[i] = true;

                let down = downstream(i, &upstream);
                pending[down] = pending[down].saturating_sub(1);
                if pending[down] == 0 && !resolved[down] {
                    resolved[down] = true;
                    queue.push_back(down);
                }
            }
        }

        for i in 0..n_links {
            if done_link[i] {
                prev_outlet[i] = outlet[i];
            }
        }

        if !progress(k, steps) {
            return Err(SimError::Cancelled);
        }
    }

    Ok(Results { times: hyd.times(), node_quality, diagnostics: diag })
}

fn is_active(link: &crate::network::Link, q: f64, settings: &Settings) -> bool {
    if link.is_pass_through() {
        q.abs() > 0.0
    } else {
        q.abs() / link.area() > settings.stagnant_velocity
    }
}

enum Transport {
    PassThrough,
    Analytic(AnalyticPipe),
    Legacy(PipeState),
}

impl Transport {
    /// Start-of-step bookkeeping for flow `q`. Returns the number of period restarts (0 or 1).
    fn prepare(&mut self, link: &crate::network::Link, q: f64, settings: &Settings, dt: f64, decay: f64) -> usize {
        let active = is_active(link, q, settings);
        match self {
            Transport::PassThrough => 0,
            Transport::Legacy(pipe) => {
                std::mem::swap(&mut pipe.prev, &mut pipe.cur);
                if !active {
                    // Stagnant water just decays in place. (The Python version
                    // used (1 - K) per step, which ignores the step length.)
                    for (c, p) in pipe.cur.iter_mut().zip(&pipe.prev) {
                        *c = p * decay;
                    }
                }
                0
            }
            Transport::Analytic(pipe) => {
                let tol = settings.period_tolerance;
                if active {
                    let u = q.abs() / link.area();
                    pipe.begin_step(q, u, settings.dispersion.coefficient(u, link.diameter), settings.k, tol) as usize
                } else {
                    let restarted = pipe.begin_step(q, 0.0, 0.0, settings.k, tol);
                    pipe.idle(dt);
                    restarted as usize
                }
            }
        }
    }
}

/// Concentrations at the segment boundaries of one pipe, ordered from
/// `Link::from` (index 0) to `Link::to` (last index) regardless of flow
/// direction, so a flow reversal does not scramble the stored profile.
#[derive(Debug, Clone, Default)]
struct PipeState {
    segment: f64,
    prev: Vec<f64>,
    cur: Vec<f64>,
}

impl PipeState {
    /// Advances the profile by one step and returns the outlet concentration.
    fn advance(&mut self, seg: &Segment, forward: bool, inlet_before: f64, inlet_now: f64, diag: &mut Diagnostics) -> f64 {
        let last = self.cur.len() - 1;
        let at = |j: usize| if forward { j } else { last - j };
        self.prev[at(0)] = inlet_before;
        self.cur[at(0)] = inlet_now;
        for j in 0..last {
            let (a, b) = (at(j), at(j + 1));
            let r = seg.step(self.prev[a], self.cur[a], self.prev[a]);
            diag.plug_flow_fallbacks += r.fell_back as usize;
            self.cur[b] = if r.value < 0.0 {
                diag.negative_clamped += 1;
                0.0
            } else {
                r.value
            };
        }
        self.cur[at(last)]
    }
}

fn validate(net: &Network, hyd: &Hydraulics, s: &Settings) -> Result<(), SimError> {
    let err = |m: String| Err(SimError::InvalidInput(m));
    if !(hyd.dt > 0.0) {
        return err(format!("time step must be positive, got {}", hyd.dt));
    }
    if hyd.flows.len() < 2 {
        return err("need at least one time step of hydraulic results".into());
    }
    if !(s.segment_length > 0.0) {
        return err(format!("segment length must be positive, got {}", s.segment_length));
    }
    if !(s.k >= 0.0) || !s.k.is_finite() {
        return err(format!("decay coefficient must be a non-negative number, got {}", s.k));
    }
    if hyd.boundary.len() != net.nodes.len() {
        return err(format!("{} boundary entries for {} nodes", hyd.boundary.len(), net.nodes.len()));
    }
    for (k, row) in hyd.flows.iter().enumerate() {
        if row.len() != net.links.len() {
            return err(format!("step {k}: {} flows for {} links", row.len(), net.links.len()));
        }
        if let Some(i) = row.iter().position(|q| !q.is_finite()) {
            return err(format!("step {k}: flow in link '{}' is not finite", net.links[i].id));
        }
    }
    for (n, b) in hyd.boundary.iter().enumerate() {
        if let Some(b) = b {
            if b.len() != hyd.flows.len() {
                return err(format!("boundary series for node '{}' has {} values, expected {}", net.nodes[n].id, b.len(), hyd.flows.len()));
            }
        }
    }
    for l in &net.links {
        if l.from >= net.nodes.len() || l.to >= net.nodes.len() {
            return err(format!("link '{}' refers to a missing node", l.id));
        }
    }
    if hyd.boundary.iter().all(Option::is_none) {
        return err("network has no source (reservoir, tank or injection node)".into());
    }
    Ok(())
}
