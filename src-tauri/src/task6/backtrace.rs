#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct SolverBacktraceData {
    pub step: usize,
    pub discrepancy_norm: f64,     // ||F - AUk|| = max|Lh(uij) + fij|
    pub relative_discrepancy: f64, // disrepancy norm / ||F - AU0||
    pub norm_adj_difference: f64,  // ||Uk - Uk-1||
    pub apposterior_est: f64,      // p(H) / (1 - p(H)) * norm adj difference
    // p(H) for optimal simple iteration (and simple iteration)
    // = (Delta - delta) / (Delta + delta)
    pub spectral_radius_approach: f64, // pk = norm adj difference [cur] / norm adj difference [prev]

    pub u: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct SolverBacktrace {
    pub datas: Vec<SolverBacktraceData>,

    pub tau: f64,
    pub cheb_params: Vec<usize>,
    pub taus: Vec<f64>,
}

impl SolverBacktrace {
    pub fn new() -> SolverBacktrace {
        SolverBacktrace {
            datas: Vec::new(),
            tau: 0_f64,
            cheb_params: Vec::<usize>::new(),
            taus: Vec::<f64>::new(),
        }
    }

    pub fn add_data(&mut self, sbd: SolverBacktraceData) {
        self.datas.push(sbd);
    }

    pub fn clear(&mut self) {
        self.datas.clear();
        self.tau = 0_f64;
        self.cheb_params.clear();
        self.taus.clear();
    }
}
