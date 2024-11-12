use crate::{
    // task1::RitzData,
    // task1_concrete::{System13, WFunctionSystem},
    task6::EllipticEquationSolver,
    task6_concrete::EE3,
};

pub struct Storage {
    pub is_busy: bool,
    // pub task1data: Option<RitzData<System13, WFunctionSystem>>,
    pub task6data: Option<EllipticEquationSolver<EE3>>,
}

impl Storage {
    pub fn default() -> Storage {
        Storage {
            is_busy: false,
            // task1data: None,
            task6data: None,
        }
    }
}
