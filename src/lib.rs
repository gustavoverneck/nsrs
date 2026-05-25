pub mod core;
pub use core::constants;

// Reexports no crate root para facilitar imports:
// use nsrs::{GM1, GM3, HadronsMatter, NlemModel, EngineMode, Solver, Artist, generate_mr_curve, write_eos_with_mr};
pub use core::io_utils::{read_eos_file, write_eos_with_mr};
pub use core::model::{GM1, GM3};
pub use core::physics::{HadronsMatter, MagneticTopology, NlemModel};
pub use core::darkphotons::DarkPhotonsMatter;
pub use core::quarks::QuarksMatter;
pub use core::hybrid::HybridMatter;
pub use core::plotting::Artist;
pub use core::solver::{EngineMode, Solver};
pub use core::tov_solver::generate_mr_curve;
