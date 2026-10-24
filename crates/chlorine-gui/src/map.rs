//! Network map: nodes and links coloured by concentration at the selected time.

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};

use chlorine_core::{LinkKind, NodeKind};

use crate::worker::Output;

pub struct MapView {
    zoom: f32,
    pan: Vec2,
    pub show_ids: bool,
}

impl Default for MapView {
    fn default() -> Self {
        MapView { zoom: 1.0, pan: Vec2::ZERO, show_ids: false }
    }
}

impl MapView {
    pub fn reset(&mut self) {
        self.zoom = 1.0;
        self.pan = Vec2::ZERO;
    }
}

/// Five-stop blue → teal → green → yellow → red ramp for `v` in [0, 1].
pub fn ramp(v: f64) -> Color32 {
    const STOPS: [(f32, f32, f32); 5] = [
        (49.0, 54.0, 149.0),
        (33.0, 145.0, 140.0),
        (94.0, 201.0, 98.0),
        (253.0, 231.0, 37.0),
        (215.0, 48.0, 39.0),
    ];
    let v = v.clamp(0.0, 1.0) as f32 * (STOPS.len() - 1) as f32;
    let i = (v.floor() as usize).min(STOPS.len() - 2);
    let f = v - i as f32;
    let (a, b) = (STOPS[i], STOPS[i + 1]);
    let mix = |x: f32, y: f32| (x + (y - x) * f).round() as u8;
    Color32::from_rgb(mix(a.0, b.0), mix(a.1, b.1), mix(a.2, b.2))
}

/// World positions of the nodes: file coordinates, or a circle if any are missing.
fn positions(out: &Output) -> Vec<Pos2> {
    let nodes = &out.loaded.network.nodes;
    if nodes.iter().all(|n| n.coords.is_some()) {
        nodes.iter().map(|n| n.coords.map(|(x, y)| Pos2::new(x as f32, y as f32)).unwrap()).collect()
    } else {
        let n = nodes.len().max(1) as f32;
        (0..nodes.len())
            .map(|i| {
                let a = i as f32 / n * std::f32::consts::TAU;
                Pos2::new(a.cos() * 100.0, a.sin() * 100.0)
            })
            .collect()
    }
}

/// Draws the map and returns the node the user clicked, if any.
pub fn show(ui: &mut egui::Ui, out: &Output, k: usize, selected: Option<usize>, view: &mut MapView, scale_max: f64) -> Option<usize> {
    let (response, painter) = ui.allocate_painter(ui.available_size(), Sense::click_and_drag());
    let rect = response.rect;
    painter.rect_filled(rect, 4.0, ui.visuals().extreme_bg_color);

    let world = positions(out);
    let (mut min, mut max) = (Pos2::new(f32::MAX, f32::MAX), Pos2::new(f32::MIN, f32::MIN));
    for p in &world {
        min = min.min(*p);
        max = max.max(*p);
    }
    let size = (max - min).max(Vec2::splat(1e-3));
    let fit = ((rect.width() - 40.0) / size.x).min((rect.height() - 40.0) / size.y).max(1e-6);
    let world_center = min + size / 2.0;

    // Pan by dragging, zoom with the wheel around the pointer, reset on double click.
    if response.dragged() {
        view.pan += response.drag_delta();
    }
    if response.double_clicked() {
        view.reset();
    }
    if let (true, Some(pointer)) = (response.hovered(), response.hover_pos()) {
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            let factor = (scroll * 0.002).exp();
            let offset = pointer - rect.center() - view.pan;
            view.zoom = (view.zoom * factor).clamp(0.05, 200.0);
            view.pan = pointer - rect.center() - offset * factor;
        }
    }
    let to_screen = |p: Pos2| {
        let d = (p - world_center) * fit * view.zoom;
        rect.center() + view.pan + Vec2::new(d.x, -d.y)
    };
    let screen: Vec<Pos2> = world.iter().map(|p| to_screen(*p)).collect();
    let painter = painter.with_clip_rect(rect);

    let q = &out.results.node_quality;
    let norm = |c: f64| if scale_max > 0.0 { c / scale_max } else { 0.0 };
    for link in &out.loaded.network.links {
        let c = 0.5 * (q[link.from][k] + q[link.to][k]);
        let width = if link.kind == LinkKind::Pipe { 2.5 } else { 4.0 };
        painter.line_segment([screen[link.from], screen[link.to]], Stroke::new(width, ramp(norm(c))));
    }
    for (n, node) in out.loaded.network.nodes.iter().enumerate() {
        let color = ramp(norm(q[n][k]));
        let p = screen[n];
        match node.kind {
            NodeKind::Junction => {
                painter.circle_filled(p, 4.0, color);
            }
            NodeKind::Reservoir => {
                painter.rect_filled(Rect::from_center_size(p, Vec2::splat(12.0)), 1.0, color);
            }
            NodeKind::Tank => {
                painter.rect_filled(Rect::from_center_size(p, Vec2::new(16.0, 11.0)), 3.0, color);
            }
        }
        if Some(n) == selected {
            painter.circle_stroke(p, 9.0, Stroke::new(2.0, ui.visuals().strong_text_color()));
        }
        if view.show_ids {
            painter.text(p + Vec2::new(7.0, -7.0), Align2::LEFT_BOTTOM, &node.id, FontId::proportional(11.0), ui.visuals().weak_text_color());
        }
    }

    // Nearest node under the pointer.
    let hovered = response.hover_pos().and_then(|pointer| {
        screen
            .iter()
            .enumerate()
            .map(|(i, p)| (i, p.distance(pointer)))
            .filter(|(_, d)| *d < 10.0)
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _)| i)
    });
    if let Some(n) = hovered {
        let node = &out.loaded.network.nodes[n];
        let text = format!(
            "{}\nmodel  {:.4} mg/L\nEPANET {:.4} mg/L",
            node.id, q[n][k], out.loaded.reference[n][k]
        );
        let galley = painter.layout_no_wrap(text, FontId::monospace(12.0), ui.visuals().text_color());
        let pos = screen[n] + Vec2::new(12.0, -galley.size().y - 6.0);
        let bg = Rect::from_min_size(pos, galley.size()).expand(5.0);
        painter.rect_filled(bg, 4.0, ui.visuals().window_fill);
        painter.rect_stroke(bg, 4.0, ui.visuals().window_stroke, egui::StrokeKind::Outside);
        painter.galley(pos, galley, ui.visuals().text_color());
    }

    // Colour legend.
    let legend = Rect::from_min_size(rect.right_bottom() - Vec2::new(190.0, 34.0), Vec2::new(170.0, 10.0));
    let steps = 40;
    for i in 0..steps {
        let x0 = legend.left() + legend.width() * i as f32 / steps as f32;
        let x1 = legend.left() + legend.width() * (i + 1) as f32 / steps as f32;
        let r = Rect::from_min_max(Pos2::new(x0, legend.top()), Pos2::new(x1, legend.bottom()));
        painter.rect_filled(r, 0.0, ramp(i as f64 / (steps - 1) as f64));
    }
    let font = FontId::proportional(11.0);
    let color = ui.visuals().text_color();
    painter.text(legend.left_bottom() + Vec2::new(0.0, 3.0), Align2::LEFT_TOP, "0", font.clone(), color);
    painter.text(legend.right_bottom() + Vec2::new(0.0, 3.0), Align2::RIGHT_TOP, format!("{scale_max:.2} mg/L"), font, color);

    (response.clicked()).then_some(hovered).flatten()
}
