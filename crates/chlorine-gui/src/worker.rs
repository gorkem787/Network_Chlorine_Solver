//! Runs EPANET and the solver on a background thread.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::mpsc::{self, Receiver};
use std::sync::Arc;
use std::time::{Duration, Instant};

use chlorine_core::{simulate_with_progress, Dispersion, Results, Scheme, Settings, SimError};
use chlorine_epanet::report::{self, NodeStats};
use chlorine_epanet::{load, LoadOptions, Loaded};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Species {
    Chlorine,
    Fluorine,
}

#[derive(Debug, Clone)]
pub struct Params {
    /// Bulk decay coefficient, 1/day.
    pub k_per_day: f64,
    /// Time step, s.
    pub dt: u32,
    /// `None` uses the input file's duration.
    pub hours: Option<f64>,
    pub segment: f64,
    pub species: Species,
    pub scheme: Scheme,
    pub dispersion_scale: f64,
    pub period_tolerance: f64,
    pub keep_epanet_reactions: bool,
}

impl Default for Params {
    fn default() -> Self {
        let s = Settings::default();
        Params {
            k_per_day: 1.0,
            dt: 900,
            hours: None,
            segment: s.segment_length,
            species: Species::Chlorine,
            scheme: s.scheme,
            dispersion_scale: 1.0,
            period_tolerance: s.period_tolerance,
            keep_epanet_reactions: false,
        }
    }
}

pub struct Output {
    pub path: PathBuf,
    pub params: Params,
    pub loaded: Loaded,
    pub results: Results,
    pub stats: Vec<NodeStats>,
    pub overall_rmse: Option<f64>,
    pub load_time: Duration,
    pub solve_time: Duration,
}

pub enum Stage {
    Epanet,
    Solving,
}

pub struct Job {
    pub stage: Arc<AtomicUsize>, // 0 = EPANET, 1 = solving
    pub step: Arc<AtomicUsize>,
    pub total: Arc<AtomicUsize>,
    pub cancel: Arc<AtomicBool>,
    pub done: Receiver<Result<Output, String>>,
}

impl Job {
    pub fn stage(&self) -> Stage {
        if self.stage.load(Ordering::Relaxed) == 0 { Stage::Epanet } else { Stage::Solving }
    }

    pub fn fraction(&self) -> f32 {
        let total = self.total.load(Ordering::Relaxed);
        if total == 0 { 0.0 } else { self.step.load(Ordering::Relaxed) as f32 / total as f32 }
    }
}

pub fn start(path: PathBuf, params: Params, ctx: eframe::egui::Context) -> Job {
    let (tx, rx) = mpsc::channel();
    let job = Job {
        stage: Arc::default(),
        step: Arc::default(),
        total: Arc::default(),
        cancel: Arc::default(),
        done: rx,
    };
    let (stage, step, total, cancel) = (job.stage.clone(), job.step.clone(), job.total.clone(), job.cancel.clone());
    std::thread::spawn(move || {
        let result = run(path, params, &stage, &step, &total, &cancel, &ctx);
        let _ = tx.send(result);
        ctx.request_repaint();
    });
    job
}

fn run(
    path: PathBuf,
    params: Params,
    stage: &AtomicUsize,
    step: &AtomicUsize,
    total: &AtomicUsize,
    cancel: &AtomicBool,
    ctx: &eframe::egui::Context,
) -> Result<Output, String> {
    let k = params.k_per_day / 86400.0;
    let started = Instant::now();
    let loaded = load(
        &path,
        &LoadOptions {
            dt: params.dt as f64,
            duration: params.hours.map(|h| h * 3600.0),
            sync_decay: (!params.keep_epanet_reactions).then_some(k),
        },
    )
    .map_err(|e| e.to_string())?;
    let load_time = started.elapsed();
    if cancel.load(Ordering::Relaxed) {
        return Err(SimError::Cancelled.to_string());
    }

    let base = match params.species {
        Species::Chlorine => Dispersion::CHLORINE,
        Species::Fluorine => Dispersion::FLUORINE,
    };
    let settings = Settings {
        k,
        segment_length: params.segment,
        dispersion: Dispersion { coef: base.coef * params.dispersion_scale, ..base },
        scheme: params.scheme,
        period_tolerance: params.period_tolerance,
        ..Default::default()
    };
    stage.store(1, Ordering::Relaxed);
    let started = Instant::now();
    let mut last_repaint = Instant::now();
    let results = simulate_with_progress(&loaded.network, &loaded.hydraulics, &settings, |k, n| {
        step.store(k, Ordering::Relaxed);
        total.store(n, Ordering::Relaxed);
        if last_repaint.elapsed() > Duration::from_millis(50) {
            ctx.request_repaint();
            last_repaint = Instant::now();
        }
        !cancel.load(Ordering::Relaxed)
    })
    .map_err(|e| e.to_string())?;
    let solve_time = started.elapsed();

    let stats = report::node_stats(&loaded, &results);
    let overall_rmse = report::overall_rmse(&stats);
    Ok(Output { path, params, loaded, results, stats, overall_rmse, load_time, solve_time })
}
