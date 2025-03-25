use crate::task6::EllipticEquationSolver;

pub struct Storage {
    pub is_busy: bool,
    pub task6data: Option<EllipticEquationSolver>,
}

impl Storage {
    pub fn default() -> Storage {
        Storage {
            is_busy: false,
            task6data: None,
        }
    }
}
