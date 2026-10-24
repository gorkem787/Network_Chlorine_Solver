use std::path::PathBuf;
use std::sync::atomic::Ordering;
use std::time::{Duration, Instant};

use eframe::egui::{self, RichText};
use egui_plot::{Legend, Line, LineStyle, Plot, PlotPoints, VLine};

use chlorine_core::Scheme;
use chlorine_epanet::report::{self, NodeRole};

use crate::about::About;
use crate::map::{self, MapView};
use crate::worker::{self, Job, Output, Params, Species, Stage};

pub struct App {
    path: Option<PathBuf>,
    params: Params,
    use_file_duration: bool,
    hours: f64,
    job: Option<Job>,
    output: Option<Output>,
    error: Option<String>,
    notice: Option<String>,
    time_index: usize,
    selected: Option<usize>,
    view: MapView,
    filter: String,
    playing: bool,
    last_tick: Instant,
    autorun: bool,
    about: About,
}

impl App {
    pub fn new(initial: Option<PathBuf>, autorun: bool) -> Self {
        App {
            autorun,
            about: About::default(),
            path: initial,
            params: Params::default(),
            use_file_duration: true,
            hours: 24.0,
            job: None,
            output: None,
            error: None,
            notice: None,
            time_index: 0,
            selected: None,
            view: MapView::default(),
            filter: String::new(),
            playing: false,
            last_tick: Instant::now(),
        }
    }

    fn pick_file(&mut self) {
        let mut dialog = rfd::FileDialog::new().add_filter("EPANET (*.inp)", &["inp"]);
        if let Some(dir) = self.path.as_ref().and_then(|p| p.parent()) {
            dialog = dialog.set_directory(dir);
        }
        if let Some(p) = dialog.pick_file() {
            self.path = Some(p);
            self.error = None;
        }
    }

    fn run(&mut self, ctx: &egui::Context) {
        let Some(path) = self.path.clone() else { return };
        let mut params = self.params.clone();
        params.hours = (!self.use_file_duration).then_some(self.hours);
        self.error = None;
        self.notice = None;
        self.playing = false;
        self.job = Some(worker::start(path, params, ctx.clone()));
    }

    fn export_csv(&mut self) {
        let Some(out) = &self.output else { return };
        let stem = out.path.file_stem().map(|s| s.to_string_lossy().into_owned()).unwrap_or_else(|| "results".into());
        let mut dialog = rfd::FileDialog::new()
            .add_filter("CSV", &["csv"])
            .set_file_name(format!("{stem}_results.csv"));
        if let Some(dir) = out.path.parent() {
            dialog = dialog.set_directory(dir);
        }
        if let Some(p) = dialog.save_file() {
            match report::write_csv(&p, &out.loaded, &out.results) {
                Ok(()) => self.notice = Some(format!("Saved {}", p.display())),
                Err(e) => self.error = Some(format!("Could not write CSV: {e}")),
            }
        }
    }

    fn poll_job(&mut self) {
        let Some(job) = &self.job else { return };
        if let Ok(result) = job.done.try_recv() {
            self.job = None;
            match result {
                Ok(out) => {
                    let steps = out.results.times.len().saturating_sub(1);
                    self.time_index = self.time_index.min(steps);
                    if self.time_index == 0 {
                        self.time_index = steps;
                    }
                    // Keep the selection if the same network was re-run.
                    let same_network = self.output.as_ref().is_some_and(|o| o.loaded.network.nodes.len() == out.loaded.network.nodes.len());
                    if !same_network {
                        self.selected = out.stats.iter().position(|s| s.role == NodeRole::Junction);
                        self.view.reset();
                    }
                    self.output = Some(out);
                }
                Err(e) => self.error = Some(e),
            }
        }
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.heading("Network Chlorine Solver");
            ui.separator();
            if ui.button("📂 Open network…").clicked() {
                self.pick_file();
            }
            let can_run = self.path.is_some() && self.job.is_none();
            if ui.add_enabled(can_run, egui::Button::new("▶ Run")).clicked() {
                self.run(ui.ctx());
            }
            if ui.add_enabled(self.output.is_some(), egui::Button::new("💾 Export CSV…")).clicked() {
                self.export_csv();
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                egui::widgets::global_theme_preference_buttons(ui);
                if ui.button("ℹ About").clicked() {
                    self.about.open = true;
                }
            });
        });
    }

    fn parameters(&mut self, ui: &mut egui::Ui) {
        ui.heading("Network");
        match &self.path {
            Some(p) => {
                ui.label(RichText::new(p.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default()).strong());
                ui.label(RichText::new(p.parent().map(|d| d.display().to_string()).unwrap_or_default()).small().weak());
            }
            None => {
                ui.label("No file selected.");
                ui.label(RichText::new("You can also drop an .inp file onto the window.").small().weak());
            }
        }

        ui.add_space(8.0);
        ui.heading("Parameters");
        egui::Grid::new("params").num_columns(2).spacing([8.0, 6.0]).show(ui, |ui| {
            ui.label("Decay coefficient K");
            ui.add(egui::DragValue::new(&mut self.params.k_per_day).speed(0.05).range(0.0..=100.0).suffix(" 1/day"));
            ui.end_row();

            ui.label("Time step");
            ui.add(egui::DragValue::new(&mut self.params.dt).speed(10).range(1..=86400).suffix(" s"));
            ui.end_row();

            ui.label("Duration");
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.use_file_duration, "from file");
                ui.add_enabled(
                    !self.use_file_duration,
                    egui::DragValue::new(&mut self.hours).speed(1.0).range(0.1..=10000.0).suffix(" h"),
                );
            });
            ui.end_row();

            ui.label("Substance");
            egui::ComboBox::from_id_salt("species")
                .selected_text(match self.params.species {
                    Species::Chlorine => "Chlorine",
                    Species::Fluorine => "Fluoride",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.params.species, Species::Chlorine, "Chlorine");
                    ui.selectable_value(&mut self.params.species, Species::Fluorine, "Fluoride");
                });
            ui.end_row();

            ui.label("Scheme");
            egui::ComboBox::from_id_salt("scheme")
                .selected_text(match self.params.scheme {
                    Scheme::Analytic => "Analytic",
                    Scheme::Legacy => "Legacy (Python)",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.params.scheme, Scheme::Analytic, "Analytic");
                    ui.selectable_value(&mut self.params.scheme, Scheme::Legacy, "Legacy (Python)");
                });
            ui.end_row();

            ui.label("Cell length");
            ui.add(egui::DragValue::new(&mut self.params.segment).speed(0.1).range(0.05..=1000.0).suffix(" m"));
            ui.end_row();
        });

        egui::CollapsingHeader::new("Advanced").show(ui, |ui| {
            egui::Grid::new("advanced").num_columns(2).spacing([8.0, 6.0]).show(ui, |ui| {
                ui.label("Dispersion factor").on_hover_text("1 = empirical D. Values towards 0 approach plug flow (EPANET).");
                ui.add(egui::DragValue::new(&mut self.params.dispersion_scale).speed(0.01).range(0.0..=100.0));
                ui.end_row();

                ui.label("Period tolerance").on_hover_text("Relative flow changes below this do not start a new velocity period.");
                ui.add(egui::DragValue::new(&mut self.params.period_tolerance).speed(1e-5).range(0.0..=0.1));
                ui.end_row();
            });
            ui.checkbox(&mut self.params.keep_epanet_reactions, "Use the .inp reaction coefficients in the EPANET reference");
        });

        ui.add_space(8.0);
        if let Some(job) = &self.job {
            let label = match job.stage() {
                Stage::Epanet => "Running EPANET hydraulics…".to_string(),
                Stage::Solving => format!(
                    "Solving: step {} / {}",
                    job.step.load(Ordering::Relaxed),
                    job.total.load(Ordering::Relaxed)
                ),
            };
            ui.label(label);
            ui.add(egui::ProgressBar::new(job.fraction()).show_percentage());
            if ui.button("■ Cancel").clicked() {
                job.cancel.store(true, Ordering::Relaxed);
            }
        } else {
            let run = ui.add_enabled(self.path.is_some(), egui::Button::new(RichText::new("▶ Run").size(16.0)).min_size(egui::vec2(ui.available_width(), 32.0)));
            if run.clicked() {
                self.run(ui.ctx());
            }
        }

        if let Some(e) = &self.error {
            ui.add_space(6.0);
            ui.colored_label(ui.visuals().error_fg_color, e);
        }
        if let Some(n) = &self.notice {
            ui.add_space(6.0);
            ui.label(RichText::new(n).small());
        }

        if let Some(out) = &self.output {
            ui.add_space(10.0);
            ui.separator();
            ui.heading("Last run");
            let net = &out.loaded.network;
            ui.label(format!("{} nodes, {} links", net.nodes.len(), net.links.len()));
            ui.label(format!("{} steps × {} s, flow units {}", out.loaded.hydraulics.steps(), out.params.dt, out.loaded.flow_units));
            ui.label(format!("EPANET {:.2?}, solver {:.2?}", out.load_time, out.solve_time));
            if let Some(r) = out.overall_rmse {
                ui.label(format!("RMSE vs EPANET: {r:.4} mg/L"));
            }
            let d = out.results.diagnostics;
            if d.cycle_breaks + d.plug_flow_fallbacks + d.negative_clamped > 0 {
                ui.label(
                    RichText::new(format!(
                        "Warning: {} cycle breaks, {} plug-flow fallbacks, {} negative values clamped",
                        d.cycle_breaks, d.plug_flow_fallbacks, d.negative_clamped
                    ))
                    .color(ui.visuals().warn_fg_color),
                );
            }
            for w in &out.loaded.warnings {
                ui.label(RichText::new(w).small().color(ui.visuals().warn_fg_color));
            }
        }
    }

    fn node_table(&mut self, ui: &mut egui::Ui) {
        let Some(out) = &self.output else { return };
        ui.heading("Nodes");
        ui.add(egui::TextEdit::singleline(&mut self.filter).hint_text("Search…"));
        ui.add_space(4.0);
        let filter = self.filter.to_lowercase();
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
            egui::Grid::new("nodes").num_columns(4).striped(true).spacing([10.0, 4.0]).show(ui, |ui| {
                ui.label(RichText::new("ID").strong());
                ui.label(RichText::new("Type").strong());
                ui.label(RichText::new("RMSE").strong());
                ui.label(RichText::new("max |Δ|").strong());
                ui.end_row();
                for (n, (node, s)) in out.loaded.network.nodes.iter().zip(&out.stats).enumerate() {
                    if !filter.is_empty() && !node.id.to_lowercase().contains(&filter) {
                        continue;
                    }
                    if ui.selectable_label(self.selected == Some(n), &node.id).clicked() {
                        self.selected = Some(n);
                    }
                    ui.label(s.role.label());
                    ui.label(format!("{:.4}", s.rmse));
                    ui.label(format!("{:.4}", s.max_abs_diff));
                    ui.end_row();
                }
            });
        });
    }

    fn time_controls(&mut self, ui: &mut egui::Ui, steps: usize, dt: f64) {
        ui.horizontal(|ui| {
            if ui.button(if self.playing { "⏸" } else { "▶" }).on_hover_text("Play / pause").clicked() {
                self.playing = !self.playing;
                if self.playing && self.time_index >= steps {
                    self.time_index = 0;
                }
                self.last_tick = Instant::now();
            }
            let hours = self.time_index as f64 * dt / 3600.0;
            ui.label(format!("t = {hours:6.2} h"));
            ui.checkbox(&mut self.view.show_ids, "IDs");
            ui.spacing_mut().slider_width = (ui.available_width() - 20.0).max(100.0);
            ui.add(egui::Slider::new(&mut self.time_index, 0..=steps).show_value(false));
        });
        if self.playing {
            if self.last_tick.elapsed() >= Duration::from_millis(80) {
                self.last_tick = Instant::now();
                if self.time_index < steps {
                    self.time_index += 1;
                } else {
                    self.playing = false;
                }
            }
            ui.ctx().request_repaint_after(Duration::from_millis(80));
        }
    }

    fn central(&mut self, ui: &mut egui::Ui) {
        let Some(out) = &self.output else {
            ui.centered_and_justified(|ui| {
                ui.label(RichText::new("Open an EPANET network and press ▶ Run.").size(18.0).weak());
            });
            return;
        };
        let steps = out.results.times.len() - 1;
        let dt = out.loaded.hydraulics.dt;
        self.time_controls(ui, steps, dt);
        let out = self.output.as_ref().unwrap();

        let scale_max = out
            .results
            .node_quality
            .iter()
            .flatten()
            .chain(out.loaded.reference.iter().flatten())
            .fold(0.0f64, |m, &c| m.max(c));

        let total = ui.available_height();
        let map_height = (total * 0.58).max(200.0);
        let mut clicked = None;
        ui.allocate_ui(egui::vec2(ui.available_width(), map_height), |ui| {
            clicked = map::show(ui, out, self.time_index, self.selected, &mut self.view, scale_max);
        });
        ui.label(RichText::new("Drag: pan · Wheel: zoom · Double-click: fit · Click: select node").small().weak());
        if clicked.is_some() {
            self.selected = clicked;
        }

        ui.add_space(4.0);
        let Some(n) = self.selected else {
            ui.label("Select a node on the map or in the list to plot it.");
            return;
        };
        let node = &out.loaded.network.nodes[n];
        ui.label(RichText::new(format!("Node {} ({})", node.id, out.stats[n].role.label())).strong());
        let to_points = |series: &[f64]| -> PlotPoints<'static> {
            out.results.times.iter().zip(series).map(|(t, c)| [t / 3600.0, *c]).collect::<Vec<_>>().into()
        };
        let model = to_points(&out.results.node_quality[n]);
        let epanet = to_points(&out.loaded.reference[n]);
        let now = self.time_index as f64 * dt / 3600.0;
        Plot::new("series")
            .legend(Legend::default())
            .x_axis_label("Time (h)")
            .y_axis_label("Concentration (mg/L)")
            .show(ui, |p| {
                p.line(Line::new("Model", model).width(2.0));
                p.line(Line::new("EPANET", epanet).style(LineStyle::dashed_dense()).width(1.5));
                p.vline(VLine::new("t", now).color(egui::Color32::GRAY));
            });
    }
}


impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        if std::mem::take(&mut self.autorun) {
            self.run(ui.ctx());
        }
        self.poll_job();

        let dropped = ui.ctx().input(|i| {
            i.raw.dropped_files.iter().map(|f| f.path()).find(|p| !p.as_os_str().is_empty()).map(|p| p.to_path_buf())
        });
        if let Some(p) = dropped {
            if p.extension().is_some_and(|e| e.eq_ignore_ascii_case("inp")) {
                self.path = Some(p);
                self.error = None;
            } else {
                self.error = Some("Only EPANET .inp files can be opened.".into());
            }
        }

        egui::Panel::top("top").show(ui, |ui| {
            ui.add_space(4.0);
            self.top_bar(ui);
            ui.add_space(4.0);
        });
        egui::Panel::left("parameters").resizable(true).default_size(290.0).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| self.parameters(ui));
        });
        if self.output.is_some() {
            egui::Panel::right("nodes").resizable(true).default_size(260.0).show(ui, |ui| self.node_table(ui));
        }
        egui::CentralPanel::default().show(ui, |ui| self.central(ui));
        self.about.show(ui.ctx());
    }
}
