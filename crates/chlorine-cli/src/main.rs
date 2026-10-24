use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use chlorine_core::{simulate, Dispersion, NodeKind, Scheme, Settings};
use chlorine_epanet::{load, LoadOptions};
use clap::{Parser, ValueEnum};

/// Simulates chlorine transport, dispersion and decay in an EPANET network and
/// compares the result with EPANET's own water quality simulation.
#[derive(Parser, Debug)]
#[command(name = "chlorine", version)]
struct Args {
    /// EPANET input file (.inp).
    input: PathBuf,

    /// Bulk decay coefficient K, 1/day.
    #[arg(short, long, default_value_t = 1.0)]
    k: f64,

    /// Time step, seconds.
    #[arg(long, default_value_t = 900.0)]
    dt: f64,

    /// Simulation duration in hours (default: the input file's duration).
    #[arg(long)]
    hours: Option<f64>,

    /// Target pipe segment length, metres.
    #[arg(long, default_value_t = 1.0)]
    segment: f64,

    /// Dispersion model.
    #[arg(long, value_enum, default_value_t = Species::Chlorine)]
    species: Species,

    /// Pipe transport scheme. `legacy` reproduces the original Python code.
    #[arg(long, value_enum, default_value_t = SchemeArg::Analytic)]
    scheme: SchemeArg,

    /// Multiplies the dispersion coefficient. Values near 0 approach plug flow,
    /// which is what EPANET models; useful to check transport without dispersion.
    #[arg(long, default_value_t = 1.0)]
    dispersion_scale: f64,

    /// Relative flow change that starts a new constant-velocity period in a
    /// pipe (analytic scheme). Smaller is more exact and slower.
    #[arg(long, default_value_t = Settings::default().period_tolerance)]
    period_tolerance: f64,

    /// Keep the input file's reaction coefficients in the EPANET reference run
    /// instead of replacing them with K.
    #[arg(long)]
    keep_epanet_reactions: bool,

    /// Write per-node time series (model and EPANET) to this CSV file.
    #[arg(short, long)]
    out: Option<PathBuf>,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SchemeArg {
    Analytic,
    Legacy,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum Species {
    Chlorine,
    Fluorine,
}

fn main() -> ExitCode {
    match run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: Args) -> Result<(), Box<dyn std::error::Error>> {
    let k = args.k / 86400.0;
    let started = Instant::now();
    let loaded = load(
        &args.input,
        &LoadOptions {
            dt: args.dt,
            duration: args.hours.map(|h| h * 3600.0),
            sync_decay: (!args.keep_epanet_reactions).then_some(k),
        },
    )?;
    let loaded_in = started.elapsed();
    for w in &loaded.warnings {
        eprintln!("warning: {w}");
    }

    let net = &loaded.network;
    let settings = Settings {
        k,
        segment_length: args.segment,
        dispersion: {
            let d = match args.species {
                Species::Chlorine => Dispersion::CHLORINE,
                Species::Fluorine => Dispersion::FLUORINE,
            };
            Dispersion { coef: d.coef * args.dispersion_scale, ..d }
        },
        period_tolerance: args.period_tolerance,
        scheme: match args.scheme {
            SchemeArg::Analytic => Scheme::Analytic,
            SchemeArg::Legacy => Scheme::Legacy,
        },
        ..Default::default()
    };
    let started = Instant::now();
    let results = simulate(net, &loaded.hydraulics, &settings)?;
    let solved_in = started.elapsed();

    println!(
        "{}: {} nodes, {} links, {} steps of {} s, flow units {}",
        args.input.display(),
        net.nodes.len(),
        net.links.len(),
        loaded.hydraulics.steps(),
        args.dt,
        loaded.flow_units
    );
    println!("EPANET load + run: {:.2?}   solver: {:.2?}", loaded_in, solved_in);
    let d = results.diagnostics;
    if d.period_restarts > 0 {
        println!("pipe period restarts: {}", d.period_restarts);
    }
    if d.cycle_breaks + d.plug_flow_fallbacks + d.negative_clamped > 0 {
        println!(
            "diagnostics: {} cycle breaks, {} plug-flow fallbacks, {} negative values clamped",
            d.cycle_breaks, d.plug_flow_fallbacks, d.negative_clamped
        );
    }

    println!();
    println!("{:<12} {:<10} {:>10} {:>10} {:>10} {:>10}", "node", "kind", "final", "EPANET", "RMSE", "max |Δ|");
    let (mut sq_sum, mut count) = (0.0, 0usize);
    for (n, node) in net.nodes.iter().enumerate() {
        let model = &results.node_quality[n];
        let reference = &loaded.reference[n];
        let diffs = model.iter().zip(reference).map(|(a, b)| a - b);
        let (sq, max) = diffs.fold((0.0, 0.0f64), |(s, m), d| (s + d * d, m.max(d.abs())));
        let rmse = (sq / model.len() as f64).sqrt();
        let is_boundary = loaded.hydraulics.boundary[n].is_some();
        if !is_boundary {
            sq_sum += sq;
            count += model.len();
        }
        let kind = match (node.kind, is_boundary) {
            (NodeKind::Reservoir, _) => "reservoir",
            (NodeKind::Tank, _) => "tank",
            (NodeKind::Junction, true) => "source",
            (NodeKind::Junction, false) => "junction",
        };
        println!(
            "{:<12} {:<10} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
            node.id,
            kind,
            model.last().unwrap(),
            reference.last().unwrap(),
            rmse,
            max
        );
    }
    if count > 0 {
        println!("\noverall RMSE over computed nodes: {:.4} mg/L", (sq_sum / count as f64).sqrt());
    }

    if let Some(path) = &args.out {
        let mut w = BufWriter::new(File::create(path)?);
        writeln!(w, "node_id,time_h,model,epanet")?;
        for (n, node) in net.nodes.iter().enumerate() {
            for (k, t) in results.times.iter().enumerate() {
                writeln!(w, "{},{},{},{}", node.id, t / 3600.0, results.node_quality[n][k], loaded.reference[n][k])?;
            }
        }
        w.flush()?;
        println!("wrote {}", path.display());
    }
    Ok(())
}
