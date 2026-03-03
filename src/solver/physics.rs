// solver/physics.rs

use crate::solver::model::{GM1, GM3, ModelParams};

struct PhysicsEngine {
    model_params: ModelParams,   // Nuclear model params: GM1, GM3 ()
    bg: f64,                     // Magnetic Field in Gauss
}


impl PhysicsEngine {
    pub fn new(model_params: ModelParams, bg: f64) -> Self {
        
        PhysicsEngine {
            model_params,
            bg,
        }
    }
}
