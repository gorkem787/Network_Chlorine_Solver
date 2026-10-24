//! Analytical pipe transport (see README, References) with time-varying velocity.
//!
//! Velocity in a pipe is piecewise constant in time (one value per hydraulic
//! period). Within a period the concentration anywhere in the pipe is, exactly,
//!
//! * Eq. 16 of the paper: superposition of the linear source units seen at the
//!   inlet since the period began, plus
//! * the evolution of the in-pipe profile present when the period began, with a
//!   zero-concentration inlet. This is the paper's `B` term generalised from a
//!   uniform initial concentration to a piecewise-constant profile (method of
//!   images on the semi-infinite domain; for a single uniform cell it reduces
//!   to `2B(0,t,0) - B(x,t,-1) - B(x,t,1)`).
//!
//! When the velocity changes, the profile is evaluated on cells and the next
//! period starts from it. Because the equation is linear this restart is exact
//! apart from representing the profile on cells. The cells move with the flow
//! (their edges shift by `u * tau` at each restart, new cells appear at the
//! inlet and small neighbours are merged), so slow flow is not lost to
//! rounding the shift to a fixed grid.

use crate::special::{erfc, exp_erfc};

/// Number of erfc length scales beyond which contributions are dropped.
const CUTOFF: f64 = 8.0;

#[derive(Debug, Clone)]
pub(crate) struct AnalyticPipe {
    length: f64,
    /// Target cell length, m.
    h: f64,
    /// Cell edges (ascending, from 0 to `length`) and values, measured from
    /// the inlet of the current period: from `Link::from` when `forward`.
    edges: Vec<f64>,
    vals: Vec<f64>,
    forward: bool,
    period: Option<Period>,
}

#[derive(Debug, Clone)]
struct Period {
    /// Signed flow defining this period; a different value starts a new one.
    q: f64,
    /// Velocity (m/s, >= 0); zero means stagnant.
    u: f64,
    d: f64,
    /// Inlet concentration at the grid times of this period (first = period start).
    inlet: Vec<f64>,
    /// Elapsed time since the period started, s.
    tau: f64,
    /// Source units before this index no longer reach x <= length.
    first_unit: usize,
}

impl AnalyticPipe {
    pub fn new(length: f64, cell_length: f64, c0: f64) -> Self {
        let n = ((length / cell_length).round() as usize).max(1);
        let h = length / n as f64;
        AnalyticPipe {
            length,
            h,
            edges: (0..=n).map(|j| if j == n { length } else { j as f64 * h }).collect(),
            vals: vec![c0; n],
            forward: true,
            period: None,
        }
    }

    /// Whether [`begin_step`](Self::begin_step) with this flow would restart an existing period.
    pub fn will_restart(&self, q: f64, rel_tol: f64) -> bool {
        self.period.as_ref().is_some_and(|p| Self::differs(p.q, q, rel_tol))
    }

    fn differs(old: f64, new: f64, rel_tol: f64) -> bool {
        (old - new).abs() > rel_tol * old.abs().max(new.abs()) || (old >= 0.0) != (new >= 0.0)
    }

    /// Prepares the step ending at `t_k`. `q` is the flow during the step; a
    /// new period begins whenever it differs from the current one by more than
    /// the relative tolerance `rel_tol`. Returns true when an existing period was restarted.
    pub fn begin_step(&mut self, q: f64, u: f64, d: f64, k: f64, rel_tol: f64) -> bool {
        let changed = match &self.period {
            None => true,
            Some(p) => Self::differs(p.q, q, rel_tol),
        };
        let mut restarted = false;
        if changed {
            if let Some(p) = self.period.take() {
                self.restart(&p, k);
                restarted = true;
            }
            let forward = q >= 0.0;
            if forward != self.forward {
                self.mirror();
            }
            self.period = Some(Period { q, u, d, inlet: Vec::new(), tau: 0.0, first_unit: 0 });
        }
        restarted
    }

    /// Advances a flowing pipe by `dt` and returns the outlet concentration.
    /// `inlet_before`/`inlet_now`: upstream node concentration at the start/end of the step.
    pub fn advance(&mut self, inlet_before: f64, inlet_now: f64, k: f64, dt: f64) -> f64 {
        let mut p = self.period.take().expect("begin_step must be called first");
        if p.inlet.is_empty() {
            p.inlet.push(inlet_before);
        }
        p.inlet.push(inlet_now);
        p.tau += dt;

        // Drop source units whose water has entirely passed the outlet.
        let w = (p.u * p.u + 4.0 * p.d * k).sqrt();
        let units = p.inlet.len() - 1;
        while p.first_unit < units {
            let tau_end = p.tau - (p.first_unit + 1) as f64 * dt;
            if tau_end > 0.0 && tau_end * w - self.length > CUTOFF * 2.0 * (p.d * tau_end).sqrt() {
                p.first_unit += 1;
            } else {
                break;
            }
        }

        let c = self.concentration(&p, self.length, k, dt);
        self.period = Some(p);
        c
    }

    /// Stagnant pipe: only the elapsed time is tracked; decay is applied when
    /// the period ends.
    pub fn idle(&mut self, dt: f64) {
        if let Some(p) = &mut self.period {
            p.tau += dt;
        }
    }

    fn mirror(&mut self) {
        // Cells kept beyond the old outlet would land upstream of the new inlet.
        let l = self.length;
        let n = self.edges[1..].partition_point(|&e| e < l) + 1;
        self.edges.truncate(n + 1);
        self.vals.truncate(n);
        self.edges[n] = l;
        self.edges = self.edges.iter().rev().map(|e| l - e).collect();
        self.vals.reverse();
        self.forward = !self.forward;
    }

    /// Replaces the profile with its state at the end of period `p`.
    fn restart(&mut self, p: &Period, k: f64) {
        if p.u <= 0.0 || p.tau <= 0.0 || p.inlet.len() < 2 {
            let decay = (-k * p.tau).exp();
            self.vals.iter_mut().for_each(|c| *c *= decay);
            return;
        }
        let dt = p.tau / (p.inlet.len() - 1) as f64;
        let shift = p.u * p.tau;
        let l = self.length;
        // The solution is for a semi-infinite pipe: water that has left still
        // disperses back towards the outlet. Keep it for a few dispersion lengths.
        let keep = l + (CUTOFF * 2.0 * (p.d * p.tau).sqrt()).max(self.h);

        // New cell edges: water that entered during the period, then the old
        // cells carried downstream.
        let mut edges = vec![0.0];
        let entered = shift.min(keep);
        let n_in = (entered / self.h).ceil().max(1.0) as usize;
        for i in 1..=n_in {
            edges.push(entered * i as f64 / n_in as f64);
        }
        for &e in &self.edges[1..] {
            if edges.last().is_some_and(|&last| last >= keep) {
                break;
            }
            edges.push(e + shift);
        }
        let end = *edges.last().unwrap();
        if end < l {
            edges.push(l);
        }
        let end = *edges.last().unwrap();
        let mut vals: Vec<f64> = edges
            .windows(2)
            .map(|w| self.concentration(p, 0.5 * (w[0] + w[1]), k, dt).max(0.0))
            .collect();

        // Merge neighbours that together are no longer than one cell.
        let (mut me, mut mv) = (vec![0.0], Vec::with_capacity(vals.len()));
        let mut start = 0.0;
        let mut acc = 0.0; // length-weighted sum of the open cell
        for (j, c) in vals.drain(..).enumerate() {
            let (a, b) = (edges[j], edges[j + 1]);
            if b - start > self.h * (1.0 + 1e-9) && a > start {
                mv.push(acc / (a - start));
                me.push(a);
                start = a;
                acc = 0.0;
            }
            acc += c * (b - a);
        }
        mv.push(acc / (end - start));
        me.push(end);
        self.edges = me;
        self.vals = mv;
    }

    /// Concentration at distance `x` from the inlet at the period's current time.
    fn concentration(&self, p: &Period, x: f64, k: f64, dt: f64) -> f64 {
        let (u, d, tau) = (p.u, p.d, p.tau);
        let w = (u * u + 4.0 * d * k).sqrt();

        // Eq. 14 term. `tau_i` = time since the unit boundary, `c_i` = inlet value there.
        let a_term = |tau_i: f64, m: f64, slope: f64, c_i: f64| {
            if tau_i <= 0.0 {
                return 0.0;
            }
            let p_exp = x / (2.0 * d) * (u + m * w);
            let z = (x + m * tau_i * w) / (2.0 * (d * tau_i).sqrt());
            exp_erfc(p_exp, z) * (slope * (tau_i / 2.0 + m * x / (2.0 * w)) + c_i / 2.0)
        };

        let mut source = 0.0;
        for i in p.first_unit..p.inlet.len() - 1 {
            let (w1, w2) = (p.inlet[i], p.inlet[i + 1]);
            let slope = (w2 - w1) / dt;
            let tau_start = tau - i as f64 * dt;
            let tau_end = tau_start - dt;
            source += a_term(tau_start, -1.0, slope, w1) + a_term(tau_start, 1.0, slope, w1)
                - a_term(tau_end, -1.0, slope, w2)
                - a_term(tau_end, 1.0, slope, w2);
        }

        source + (-k * tau).exp() * self.initial_term(p, x)
    }

    /// Evolution of the period's starting profile with a zero inlet, without decay.
    fn initial_term(&self, p: &Period, x: f64) -> f64 {
        let (u, d, tau) = (p.u, p.d, p.tau);
        let n = self.vals.len();
        let s = 2.0 * (d * tau).sqrt();
        let first_ending_after = |pos: f64| self.edges[1..].partition_point(|&e| e <= pos);

        // Direct term: the profile advected by u*tau and spread with variance
        // 2*d*tau. Cells farther than CUTOFF*s from y contribute erfc
        // differences of 0 - 0 or 2 - 2.
        //
        // Summed by parts over cell edges, Σ_j c_j/2 (E(e_j+1) - E(e_j)) with
        // E(e) = erfc((y - e)/s) becomes Σ_edges E(e) (c_left - c_right)/2,
        // so edges between equal cells cost nothing. Beyond the last cell the
        // profile continues at the last value (the paper's dC/dx -> 0 as
        // x -> infinity), which cancels the last edge's term.
        let y = x - u * tau;
        let (lo, hi) = (first_ending_after(y - CUTOFF * s), (first_ending_after(y + CUTOFF * s) + 1).min(n));
        let mut sum = if hi == n { self.vals[n - 1] } else { 0.0 };
        for j in lo..=hi.min(n - 1) {
            let left = if j > lo { self.vals[j - 1] } else { 0.0 };
            let right = if j < hi { self.vals[j] } else { 0.0 };
            let jump = left - right;
            if jump != 0.0 {
                sum += jump / 2.0 * erfc((y - self.edges[j]) / s);
            }
        }
        let (end, last) = (self.edges[n], self.vals[n - 1]);

        // Image term enforcing a zero inlet; only cells near the inlet matter.
        let bound = u * x / d - (x + u * tau).powi(2) / (s * s);
        if bound > -50.0 {
            let hi = (first_ending_after(CUTOFF * s) + 1).min(n);
            let ux_d = u * x / d;
            for j in 0..hi {
                let (x1, x2) = (self.edges[j], self.edges[j + 1]);
                sum -= self.vals[j] / 2.0 * (exp_erfc(ux_d, (x + u * tau + x1) / s) - exp_erfc(ux_d, (x + u * tau + x2) / s));
            }
            sum -= last / 2.0 * exp_erfc(ux_d, (x + u * tau + end) / s);
        }
        sum
    }
}
