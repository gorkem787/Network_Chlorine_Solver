//! Chlorine (or any first-order decaying substance) transport in water
//! distribution networks, using the analytical advection-dispersion-reaction
//! solution of the original Python `ADE_solver.py`.

pub mod ade;
mod analytic;
pub mod network;
pub mod sim;
pub mod special;

pub use ade::Dispersion;
pub use network::{Hydraulics, Link, LinkKind, Network, Node, NodeKind};
pub use sim::{simulate, simulate_with_progress, Diagnostics, Results, Scheme, Settings, SimError};
