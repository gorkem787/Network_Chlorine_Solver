//! Network description and hydraulic results, independent of where they come
//! from (EPANET, a test fixture, ...). Everything here is in SI units.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeKind {
    Junction,
    Reservoir,
    Tank,
}

#[derive(Debug, Clone)]
pub struct Node {
    pub id: String,
    pub kind: NodeKind,
    /// Initial concentration, mg/L (= g/m³).
    pub initial_quality: f64,
    /// Map coordinates, if the input file has them.
    pub coords: Option<(f64, f64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinkKind {
    Pipe,
    Pump,
    Valve,
}

#[derive(Debug, Clone)]
pub struct Link {
    pub id: String,
    pub kind: LinkKind,
    /// Index into `Network::nodes` of the start node.
    pub from: usize,
    /// Index into `Network::nodes` of the end node.
    pub to: usize,
    /// Length, m. Pumps and valves are treated as zero-length links.
    pub length: f64,
    /// Diameter, m.
    pub diameter: f64,
}

impl Link {
    /// Pumps, valves and zero-length pipes pass concentration through instantly.
    pub fn is_pass_through(&self) -> bool {
        self.kind != LinkKind::Pipe || self.length <= 0.0 || self.diameter <= 0.0
    }

    pub fn area(&self) -> f64 {
        std::f64::consts::PI * self.diameter * self.diameter / 4.0
    }
}

#[derive(Debug, Clone, Default)]
pub struct Network {
    pub nodes: Vec<Node>,
    pub links: Vec<Link>,
}

/// Hydraulic and boundary data on the simulation's uniform time grid
/// `t_k = k * dt`, `k = 0..=steps`.
#[derive(Debug, Clone)]
pub struct Hydraulics {
    pub dt: f64,
    /// `flows[k][link]`, m³/s, positive from `Link::from` to `Link::to`.
    pub flows: Vec<Vec<f64>>,
    /// `boundary[node]`: imposed concentration series (one value per grid time)
    /// for reservoirs, tanks and injection points; `None` for nodes the solver
    /// computes by mixing.
    pub boundary: Vec<Option<Vec<f64>>>,
}

impl Hydraulics {
    pub fn steps(&self) -> usize {
        self.flows.len().saturating_sub(1)
    }

    pub fn times(&self) -> Vec<f64> {
        (0..self.flows.len()).map(|k| k as f64 * self.dt).collect()
    }
}

/// Samples a piecewise-constant series (value `values[i]` holds from `times[i]`
/// until the next time) onto the grid `0, dt, 2dt, ..., steps*dt`.
///
/// EPANET reports results at every hydraulic event (tank fill/empty, control
/// actions, ...), not only at the reporting step. The Python version indexed
/// those results as if they were evenly spaced, which shifted flows in time
/// whenever an intermediate event occurred.
pub fn resample_step<T: Clone>(times: &[f64], values: &[T], dt: f64, steps: usize) -> Vec<T> {
    assert_eq!(times.len(), values.len());
    assert!(!times.is_empty(), "no samples to resample");
    let mut out = Vec::with_capacity(steps + 1);
    let mut i = 0;
    for k in 0..=steps {
        let t = k as f64 * dt;
        // Tolerance absorbs EPANET's integer-second clock.
        while i + 1 < times.len() && times[i + 1] <= t + 1e-6 {
            i += 1;
        }
        out.push(values[i].clone());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resample_holds_last_value_and_skips_intermediate_events() {
        // Events at 0, 900, 1234 (tank fills), 1800.
        let times = [0.0, 900.0, 1234.0, 1800.0];
        let values = [1, 2, 3, 4];
        let r = resample_step(&times, &values, 900.0, 3);
        assert_eq!(r, vec![1, 2, 4, 4]);
    }

    #[test]
    fn resample_handles_coarser_input() {
        let r = resample_step(&[0.0, 3600.0], &[1, 2], 900.0, 5);
        assert_eq!(r, vec![1, 1, 1, 1, 2, 2]);
    }
}
