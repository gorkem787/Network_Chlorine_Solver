use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use chlorine_core::{simulate, Dispersion, Scheme, Settings};
use chlorine_epanet::{load, report, LoadOptions};
use clap::{Parser, ValueEnum};

/// Simulates chlorine transport, dispersion and decay in an EPANET network and
/// compares the result with EPANET's own water quality simulation.
#[derive(Parser, Debug)]
#[command(
    name = "chlorine",
    version,
    after_help = "Copyright (C) 2024-2026 Gorkem Bora Candan. License GPL-3.0-or-later: <https://gnu.org/licenses/gpl.html>.\n\
                  This is free software with ABSOLUTELY NO WARRANTY. Third-party notices: THIRD_PARTY_LICENSES.txt."
)]
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
    let stats = report::node_stats(&loaded, &results);
    for (node, s) in net.nodes.iter().zip(&stats) {
        println!(
            "{:<12} {:<10} {:>10.4} {:>10.4} {:>10.4} {:>10.4}",
            node.id,
            s.role.label(),
            s.final_model,
            s.final_epanet,
            s.rmse,
            s.max_abs_diff
        );
    }
    if let Some(rmse) = report::overall_rmse(&stats) {
        println!("\noverall RMSE over computed nodes: {rmse:.4} mg/L");
    }

    if let Some(path) = &args.out {
        report::write_csv(path, &loaded, &results)?;
        println!("wrote {}", path.display());
    }
    Ok(())
}
