// No console window in release builds on Windows.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod about;
mod app;
mod map;
mod worker;

use eframe::egui;

fn main() -> eframe::Result {
    // An .inp path on the command line (or dropped on the exe) is opened at
    // start; `--run` also starts the simulation.
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    let autorun = args.iter().any(|a| a == "--run");
    let initial = args.into_iter().find(|a| a != "--run").map(std::path::PathBuf::from);
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Network Chlorine Solver")
            .with_inner_size([1360.0, 860.0])
            .with_min_inner_size([900.0, 600.0])
            .with_drag_and_drop(true),
        ..Default::default()
    };
    eframe::run_native(
        "Network Chlorine Solver",
        options,
        Box::new(move |_cc| Ok(Box::new(app::App::new(initial, autorun)))),
    )
}
