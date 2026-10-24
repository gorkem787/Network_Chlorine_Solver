//! About window: version, license and third-party notices.

use eframe::egui::{self, RichText};

const LICENSE: &str = include_str!("../../../LICENSE");
const THIRD_PARTY: &str = include_str!("../../../THIRD_PARTY_LICENSES.txt");

pub const COPYRIGHT: &str = "Copyright (C) 2024-2026 Gorkem Bora Candan";

#[derive(Default, PartialEq, Clone, Copy)]
enum Tab {
    #[default]
    About,
    License,
    ThirdParty,
}

#[derive(Default)]
pub struct About {
    pub open: bool,
    tab: Tab,
}

impl About {
    pub fn show(&mut self, ctx: &egui::Context) {
        let mut open = self.open;
        egui::Window::new("About Network Chlorine Solver")
            .open(&mut open)
            .default_size([620.0, 480.0])
            .collapsible(false)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.selectable_value(&mut self.tab, Tab::About, "About");
                    ui.selectable_value(&mut self.tab, Tab::License, "License (GPL-3.0)");
                    ui.selectable_value(&mut self.tab, Tab::ThirdParty, "Third-party notices");
                });
                ui.separator();
                match self.tab {
                    Tab::About => about_tab(ui),
                    Tab::License => text_tab(ui, LICENSE),
                    Tab::ThirdParty => text_tab(ui, THIRD_PARTY),
                }
            });
        self.open = open;
    }
}

fn about_tab(ui: &mut egui::Ui) {
    ui.heading(format!("Network Chlorine Solver {}", env!("CARGO_PKG_VERSION")));
    ui.label("Transport, dispersion and first-order decay of chlorine in water distribution networks.");
    ui.add_space(8.0);
    ui.label(COPYRIGHT);
    ui.label(
        "This program is free software: you can redistribute it and/or modify it under the terms of the \
         GNU General Public License as published by the Free Software Foundation, either version 3 of the \
         License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY; see the \
         License tab for details.",
    );
    ui.hyperlink_to("Source code", "https://github.com/gorkem787/Network_Chlorine_Solver");

    ui.add_space(10.0);
    ui.label(RichText::new("Method").strong());
    ui.label("Pipe transport uses an analytical solution of the one-dimensional advection-dispersion-reaction equation.");
    ui.add_space(6.0);
    ui.label(RichText::new("Hydraulics").strong());
    ui.label("Hydraulic and reference water quality results are computed with EPANET 2.2 (MIT License).");
    ui.hyperlink_to("OpenWaterAnalytics/EPANET", "https://github.com/OpenWaterAnalytics/EPANET");
}

fn text_tab(ui: &mut egui::Ui, text: &str) {
    egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
        ui.add(egui::Label::new(RichText::new(text).monospace().size(11.0)).wrap());
    });
}
