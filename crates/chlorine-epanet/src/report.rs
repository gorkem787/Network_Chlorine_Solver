//! Comparison with EPANET's own result and CSV export, shared by the CLI and GUI.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use chlorine_core::{NodeKind, Results};

use crate::Loaded;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeRole {
    Reservoir,
    Tank,
    /// Junction whose concentration is imposed from EPANET (injection or negative demand).
    Source,
    /// Junction computed by the solver.
    Junction,
}

impl NodeRole {
    pub fn label(self) -> &'static str {
        match self {
            NodeRole::Reservoir => "reservoir",
            NodeRole::Tank => "tank",
            NodeRole::Source => "source",
            NodeRole::Junction => "junction",
        }
    }
}

#[derive(Debug, Clone)]
pub struct NodeStats {
    pub role: NodeRole,
    pub final_model: f64,
    pub final_epanet: f64,
    pub rmse: f64,
    pub max_abs_diff: f64,
}

pub fn node_stats(loaded: &Loaded, results: &Results) -> Vec<NodeStats> {
    loaded
        .network
        .nodes
        .iter()
        .enumerate()
        .map(|(n, node)| {
            let model = &results.node_quality[n];
            let reference = &loaded.reference[n];
            let (sq, max) = model
                .iter()
                .zip(reference)
                .map(|(a, b)| a - b)
                .fold((0.0, 0.0f64), |(s, m), d| (s + d * d, m.max(d.abs())));
            let role = match (node.kind, loaded.hydraulics.boundary[n].is_some()) {
                (NodeKind::Reservoir, _) => NodeRole::Reservoir,
                (NodeKind::Tank, _) => NodeRole::Tank,
                (NodeKind::Junction, true) => NodeRole::Source,
                (NodeKind::Junction, false) => NodeRole::Junction,
            };
            NodeStats {
                role,
                final_model: *model.last().unwrap(),
                final_epanet: *reference.last().unwrap(),
                rmse: (sq / model.len() as f64).sqrt(),
                max_abs_diff: max,
            }
        })
        .collect()
}

/// RMSE over all time steps of the nodes the solver computes (not sources).
pub fn overall_rmse(stats: &[NodeStats]) -> Option<f64> {
    let computed: Vec<f64> = stats.iter().filter(|s| s.role == NodeRole::Junction).map(|s| s.rmse * s.rmse).collect();
    (!computed.is_empty()).then(|| (computed.iter().sum::<f64>() / computed.len() as f64).sqrt())
}

/// Long-format CSV: `node_id,time_h,model,epanet`.
pub fn write_csv(path: &Path, loaded: &Loaded, results: &Results) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "node_id,time_h,model,epanet")?;
    for (n, node) in loaded.network.nodes.iter().enumerate() {
        for (k, t) in results.times.iter().enumerate() {
            writeln!(w, "{},{},{},{}", node.id, t / 3600.0, results.node_quality[n][k], loaded.reference[n][k])?;
        }
    }
    w.flush()
}
