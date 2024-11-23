#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct SolverBacktraceStep {
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

impl SolverBacktraceStep {
    pub fn new(
        step: usize,
        discrepancy_norm: f64,
        u0_norm: f64,
        norm_adj_difference: f64,
        norm_adj_difference_prev: f64,
        ph1ph: f64,
        u: Vec<Vec<f64>>,
    ) -> Self {
        SolverBacktraceStep {
            step,
            discrepancy_norm,
            relative_discrepancy: discrepancy_norm / u0_norm,
            norm_adj_difference,
            apposterior_est: ph1ph * norm_adj_difference,
            spectral_radius_approach: norm_adj_difference / norm_adj_difference_prev,

            u,
        }
    }
}

#[derive(Debug, Clone, serde::Deserialize, serde::Serialize)]
pub struct SolverBacktraceData {
    pub solver_steps: Vec<SolverBacktraceStep>,

    pub tau: f64,
    pub cheb_params: Vec<usize>,
    pub taus: Vec<f64>,
}

impl SolverBacktraceData {
    pub fn new() -> SolverBacktraceData {
        SolverBacktraceData {
            solver_steps: Vec::new(),

            tau: 0_f64,
            cheb_params: Vec::<usize>::new(),
            taus: Vec::<f64>::new(),
        }
    }

    pub fn clear(&mut self) {
        self.solver_steps.clear();
        self.tau = 0_f64;
        self.cheb_params.clear();
        self.taus.clear();
    }
}

#[derive(PartialEq, PartialOrd, Debug, Clone, serde::Deserialize, serde::Serialize)]
pub enum BacktraceLevel {
    Everything,
    NoSteps,
}

#[derive(Debug, Clone)]
pub struct SolverBacktrace {
    level: BacktraceLevel,
    pub backtrace_data: SolverBacktraceData,

    steps_count: usize,
    step: usize,

    u0_norm: f64,
    norm_adj_diff_prev: f64,

    preserving_start_end: usize,
    division_middle: usize,

    counter_start_end: usize,
    counter_middle: usize,
}

pub struct SolverBacktraceCfg {
    pub level: BacktraceLevel,
    pub steps_count: usize,
    pub u0_norm: f64,
    pub preserving_start_end: usize,
    pub division_middle: usize,
}

impl SolverBacktrace {
    pub fn new() -> SolverBacktrace {
        SolverBacktrace {
            level: BacktraceLevel::Everything,
            backtrace_data: SolverBacktraceData::new(),

            steps_count: 0_usize,
            step: 0_usize,

            u0_norm: 0_f64,
            norm_adj_diff_prev: 0_f64,

            preserving_start_end: 0_usize,
            division_middle: 0_usize,

            counter_start_end: 0_usize,
            counter_middle: 0_usize,
        }
    }

    pub fn start(&mut self, cfg: SolverBacktraceCfg, tau: f64, cheb_params: Vec<usize>) {
        self.level = cfg.level;
        self.steps_count = cfg.steps_count;

        self.u0_norm = cfg.u0_norm;

        self.preserving_start_end = cfg.preserving_start_end;
        self.division_middle = cfg.division_middle;

        // Inner state
        self.counter_start_end = self.preserving_start_end;

        self.backtrace_data.tau = tau;
        self.backtrace_data.cheb_params = cheb_params;
    }

    pub fn check_current_step(&mut self) -> bool {
        self.step += 1;
        if self.counter_middle == 0 {
            if self.counter_start_end == 0 {
                self.counter_middle = self.division_middle - 1;
            } else {
                self.counter_start_end -= 1;
            }
            true
        } else {
            self.counter_middle -= 1;

            if self.steps_count - self.step < self.preserving_start_end + 2 {
                self.counter_start_end = self.preserving_start_end;
                self.counter_middle = 0;
            }
            false
        }
    }

    pub fn add_data(
        &mut self,
        cur_tau: f64,
        discrepancy_norm: f64,
        norm_adj_difference: f64,
        ph1ph: f64,
        u: Vec<Vec<f64>>,
    ) {
        self.backtrace_data.taus.push(cur_tau);
        if self.level == BacktraceLevel::Everything {
            self.backtrace_data
                .solver_steps
                .push(SolverBacktraceStep::new(
                    self.step,
                    discrepancy_norm,
                    self.u0_norm,
                    norm_adj_difference,
                    self.norm_adj_diff_prev,
                    ph1ph,
                    u,
                ));
        }
        self.norm_adj_diff_prev = norm_adj_difference;
    }

    pub fn clear(&mut self) {
        self.backtrace_data.clear();

        self.step = 0_usize;
        self.u0_norm = 1_f64;
        self.norm_adj_diff_prev = 1_f64;

        self.preserving_start_end = 0_usize;
        self.division_middle = 0_usize;

        self.counter_start_end = 0_usize;
        self.counter_middle = 0_usize;
    }
}
