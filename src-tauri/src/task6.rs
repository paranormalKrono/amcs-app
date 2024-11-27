use std::{
    collections::VecDeque,
    f64::consts::{PI, SQRT_2},
    fmt,
    sync::mpsc::{self, Receiver, Sender},
    thread::{self},
    vec,
};

pub mod backtrace;
pub mod elliptic_equation;

use backtrace::*;
use elliptic_equation::*;

type Matrix = Vec<Vec<f64>>;

#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum SolvingMethod {
    SimpleIteration,
    OptimalSimpleIteration,
    Zeidel,
    UpperRelaxation,
    IterationChebyshevsky,
    AlternatingTriangular,
    AlternatingTriangularChebyshevsky,
    AlternatingDirections,
}

#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum OptimizationType {
    Nothing,
    Parallelism,
    Multithreading,
}

#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum CalculationType {
    MaxSteps,
    Calculated,
}

#[derive(Debug)]
pub struct EllipticEquationSolver<T: EllipticEquation> {
    elliptic_equation: T,

    area: (f64, f64),
    p_boundaries: (f64, f64),
    q_boundaries: (f64, f64),
    mu: fn(f64, f64) -> f64,

    n: usize,
    m: usize,

    hx: f64,
    hy: f64,

    pub u: Vec<Vec<f64>>,
    steps_count: usize,

    pub backtrace: SolverBacktrace,
    backtrace_level: BacktraceLevel,
    preserving_start_end: usize,
    division_middle: usize,

    t: f64,
}

impl<T: EllipticEquation + fmt::Display> fmt::Display for EllipticEquationSolver<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let buf: String = format!(
            "Problem: {}\n
            Area: {:?}\n
            P boundaries: {:?}\n
            Q boundaries: {:?}\n
            N = {}\n
            M = {}\n
            hx = {}\n
            hy = {}\n
            Steps count = {}\n
            t = {}\n",
            self.elliptic_equation,
            self.area,
            self.p_boundaries,
            self.q_boundaries,
            self.n,
            self.m,
            self.hx,
            self.hy,
            self.steps_count,
            self.t,
        );

        write!(f, "{}", buf)
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum ThreadID {
    First,
    Last,
    Other,
}

struct RawTrisf64Mut {
    id: ThreadID,
    n: usize,
    l: *mut f64,
    d: *mut f64,
    r: *mut f64,
    b: *mut f64,
}

struct RawTrisf64 {
    l: *const f64,
    d: *mut f64,
    r: *const f64,
    b: *mut f64,
}

unsafe impl Send for RawTrisf64Mut {}
unsafe impl Send for RawTrisf64 {}

impl<T: EllipticEquation> EllipticEquationSolver<T> {
    pub fn new(
        elliptic_equation: T,
        area: (f64, f64),
        p_boundaries: (f64, f64),
        q_boundaries: (f64, f64),
        mu: fn(f64, f64) -> f64,
        (n, m): (usize, usize),
    ) -> EllipticEquationSolver<T> {
        EllipticEquationSolver {
            elliptic_equation,

            area,
            p_boundaries,
            q_boundaries,
            mu,

            n,
            m,

            hx: area.0 / ((n - 1) as f64),
            hy: area.1 / ((m - 1) as f64),

            u: vec![vec![0_f64; m]; n],
            steps_count: 0_usize,
            backtrace: SolverBacktrace::new(),
            backtrace_level: BacktraceLevel::None,
            preserving_start_end: 8,
            division_middle: 5,

            t: 0_f64,
        }
    }

    pub fn reset_process(&mut self) {
        self.u = vec![vec![0_f64; self.m]; self.n];
        self.steps_count = 0_usize;
        self.backtrace.clear();
    }

    pub fn set_size(&mut self, n: usize, m: usize) {
        if (self.n != n || self.m != m) && n >= 3 && m >= 3 {
            self.hx = self.area.0 / ((n - 1) as f64);
            self.hy = self.area.1 / ((m - 1) as f64);
            self.n = n;
            self.m = m;

            self.reset_process();
        }
    }

    pub fn set_backtrace(
        &mut self,
        backtrace_level: BacktraceLevel,
        first_end_preserving: usize,
        middle_division: usize,
    ) {
        self.backtrace_level = backtrace_level;
        self.preserving_start_end = first_end_preserving;
        self.division_middle = middle_division;
    }

    pub fn get_size(&self) -> (usize, usize) {
        (self.n, self.m)
    }

    pub fn get_h(&self) -> (f64, f64) {
        (self.hx, self.hy)
    }

    pub fn get_h2(&self) -> (f64, f64) {
        (self.hx * self.hx, self.hy * self.hy)
    }

    // fn check_delta_abs_equality(&self) -> bool {
    //     let (c1, c2) = self.p_boundaries;
    //     let (d1, d2) = self.q_boundaries;
    //     false
    // }

    fn get_delta(&self, solving_method: &SolvingMethod) -> (f64, f64) {
        let (hx, hy) = self.get_h();
        let (lx, ly) = self.area;

        let (c1, c2) = self.p_boundaries;
        let (d1, d2) = self.q_boundaries;

        let h1 = 4_f64 / hx.powi(2);
        let h2 = 4_f64 / hy.powi(2);

        let hl1 = PI * hx / 2_f64 / lx;
        let hl2 = PI * hy / 2_f64 / ly;
        match solving_method {
            SolvingMethod::SimpleIteration
            | SolvingMethod::OptimalSimpleIteration
            | SolvingMethod::Zeidel
            | SolvingMethod::IterationChebyshevsky => {
                let small_delta = c1 * h1 * hl1.sin().powi(2) + d1 * h2 * hl2.sin().powi(2);
                let big_delta = c2 * h1 * hl1.cos().powi(2) + d2 * h2 * hl2.cos().powi(2);

                (small_delta, big_delta)
            }
            SolvingMethod::AlternatingTriangular
            | SolvingMethod::AlternatingTriangularChebyshevsky => {
                let small_delta = c1 * h1 * hl1.sin().powi(2) + d1 * h2 * hl2.sin().powi(2);
                let big_delta = c2 * h1 + d2 * h2;

                (small_delta, big_delta)
            }
            SolvingMethod::AlternatingDirections => {
                let small_delta =
                    f64::min(c1 * h1 * hl1.sin().powi(2), d1 * h2 * hl2.sin().powi(2));
                let big_delta = f64::max(c2 * h1 * hl1.cos().powi(2), d2 * h2 * hl2.cos().powi(2));

                (small_delta, big_delta)
            }
            _ => (0_f64, 0_f64),
        }
    }

    fn calculate_boundaries(&mut self) {
        let (lx, ly) = self.area;

        let mut fi: f64;
        for i in 0..self.n {
            fi = i as f64;
            self.u[i][0] = (self.mu)(fi * self.hx, 0_f64);
            self.u[i][self.m - 1] = (self.mu)(fi * self.hx, ly);
        }

        let mut fj: f64;
        for j in 1..(self.m - 1) {
            fj = j as f64;
            self.u[0][j] = (self.mu)(0_f64, fj * self.hy);
            self.u[self.n - 1][j] = (self.mu)(lx, fj * self.hy);
        }
    }

    pub fn get_discrepancy_matrix(
        &self,
        mf: &[Vec<f64>],
        mp: &[Vec<f64>],
        mq: &[Vec<f64>],
        u: &[Vec<f64>],
    ) -> Vec<Vec<f64>> {
        let mut discr: Vec<Vec<f64>> = vec![vec![0_f64; self.m]; self.n];

        let (hx2, hy2) = self.get_h2();

        for i in 1..(self.n - 1) {
            for j in 1..(self.m - 1) {
                discr[i][j] = mf[i - 1][j - 1] + mp[i][j - 1] * (u[i + 1][j] - u[i][j]) / hx2
                    - mp[i - 1][j - 1] * (u[i][j] - u[i - 1][j]) / hx2
                    + mq[i - 1][j] * (u[i][j + 1] - u[i][j]) / hy2
                    - mq[i - 1][j - 1] * (u[i][j] - u[i][j - 1]) / hy2;
            }
        }

        discr
    }

    fn get_discrepancy_norm(
        &self,
        mf: &[Vec<f64>],
        mp: &[Vec<f64>],
        mq: &[Vec<f64>],
        u: &[Vec<f64>],
    ) -> f64 {
        let mut max_mod = 0_f64;

        let (hx2, hy2) = self.get_h2();

        for i in 1..(self.n - 1) {
            for j in 1..(self.m - 1) {
                max_mod = f64::max(
                    max_mod,
                    (mf[i - 1][j - 1] + mp[i][j - 1] * (u[i + 1][j] - u[i][j]) / hx2
                        - mp[i - 1][j - 1] * (u[i][j] - u[i - 1][j]) / hx2
                        + mq[i - 1][j] * (u[i][j + 1] - u[i][j]) / hy2
                        - mq[i - 1][j - 1] * (u[i][j] - u[i][j - 1]) / hy2)
                        .abs(),
                );
            }
        }

        max_mod
    }

    pub fn set_param(&mut self, t: f64) {
        self.t = t;
    }

    fn get_discrepancy_zero_norm(&self, mf: &[Vec<f64>]) -> f64 {
        let mut max_mod = 0_f64;

        for i in 0..mf.len() {
            for j in 0..mf[0].len() {
                max_mod = f64::max(max_mod, mf[i][j].abs());
            }
        }

        max_mod
    }

    fn get_norm_adj_difference(&self, uk_prev: &[Vec<f64>]) -> f64 {
        let mut max_mod = 0_f64;

        for (i, ui) in uk_prev.iter().enumerate().take(self.n - 1).skip(1) {
            for (j, uij) in ui.iter().enumerate().take(self.m - 1).skip(1) {
                max_mod = f64::max(max_mod, (uij - self.u[i][j]).abs())
            }
        }

        max_mod
    }

    pub fn get_interior_values(&self) -> (Matrix, Matrix, Matrix) {
        let mut mf = vec![vec![0_f64; self.m - 2]; self.n - 2];
        let mut mp = vec![vec![0_f64; self.m - 2]; self.n - 1];
        let mut mq = vec![vec![0_f64; self.m - 1]; self.n - 2];

        let (hx, hy) = self.get_h();

        // Calculating functions results
        for i in 0..mf.len() {
            for j in 0..mf[0].len() {
                mf[i][j] = self
                    .elliptic_equation
                    .f((i + 1) as f64 * hx, (j + 1) as f64 * hy);
            }
        }

        for i in 0..mp.len() {
            for j in 0..mp[0].len() {
                mp[i][j] = self
                    .elliptic_equation
                    .p((i as f64 + 0.5_f64) * hx, (j + 1) as f64 * hy);
            }
        }

        for i in 0..mq.len() {
            for j in 0..mq[0].len() {
                mq[i][j] = self
                    .elliptic_equation
                    .q((i + 1) as f64 * hx, (j as f64 + 0.5_f64) * hy);
            }
        }

        (mf, mp, mq)
    }

    pub fn solve_ees(
        &mut self,
        solving_method: SolvingMethod,
        calculation: CalculationType,
        eps: f64,
        max_steps: usize,
        optimization: OptimizationType,
        threads_count: usize,
    ) {
        if self.steps_count > 0 {
            self.reset_process();
        }

        self.calculate_boundaries();

        let (n, m) = (self.n, self.m);
        let (hx2, hy2) = self.get_h2();
        // let &mut u = self.u.borrow_mut();

        let (small_delta, big_delta): (f64, f64);

        let mut calculated_steps_count: usize = max_steps;

        // For backtrace
        let mut uk_prev: Vec<Vec<f64>> = Vec::new();
        let mut ph1ph: f64 = 0_f64;

        // Parameters
        let mut tau: f64 = 1_f64;
        let mut cheb_params: Vec<usize> = Vec::new();
        let mut cheb_coeffs: Vec<f64> = Vec::new();
        let mut cheb_coeffs_iter: std::iter::Cycle<std::slice::Iter<'_, f64>> =
            cheb_coeffs.iter().cycle();

        // Functions values
        let (mf, mp, mq) = self.get_interior_values();

        let mut ksi: f64 = 1_f64;

        // Temporary values for solving
        let mut a1: f64;
        let mut av0: Vec<f64> = Vec::new();
        let mut av01: Vec<f64> = Vec::new();
        let mut av02: Vec<f64> = Vec::new();
        let mut av1: Vec<f64> = Vec::new();
        let mut av11: Vec<f64> = Vec::new();
        let mut av12: Vec<f64> = Vec::new();
        let mut avs0: Vec<f64> = Vec::new();
        let mut avs1: Vec<f64> = Vec::new();
        let mut avt0: Vec<f64> = Vec::new();
        let mut avt1: Vec<f64> = Vec::new();
        let mut am: Vec<Vec<f64>> = Vec::new();

        // Parallel
        let (tx01, rx10) = mpsc::channel::<usize>();
        let (tx10, rx01) = mpsc::channel::<usize>();
        let (txr01, rxr10) = mpsc::channel::<RawTrisf64>();

        // Multithread
        let mut txmv0e: Vec<Sender<RawTrisf64Mut>> = Vec::with_capacity(threads_count);
        let (txme0, rxm0e) = mpsc::channel::<usize>();
        let mut txmvu0e: Vec<Sender<usize>> = Vec::with_capacity(threads_count);

        match solving_method {
            SolvingMethod::SimpleIteration => {
                (small_delta, big_delta) = self.get_delta(&solving_method);
                ksi = small_delta / big_delta;
                calculated_steps_count = ((1_f64 / eps).ln() / 2_f64 / ksi) as usize + 1;
            }
            SolvingMethod::OptimalSimpleIteration => {
                (small_delta, big_delta) = self.get_delta(&solving_method);
                ksi = small_delta / big_delta;
                calculated_steps_count = ((1_f64 / eps).ln() / 2_f64 / ksi) as usize + 1;

                tau = 2_f64 / (big_delta + small_delta);
            }
            SolvingMethod::Zeidel => {
                (small_delta, big_delta) = self.get_delta(&solving_method);
                ksi = small_delta / big_delta;
                calculated_steps_count = ((1_f64 / eps).ln() / 4_f64 / ksi) as usize + 1;
            }
            SolvingMethod::UpperRelaxation => {
                let h = (hx2 + hy2).sqrt();
                tau = 2_f64 / (1_f64 + (h * PI).sin());
                tau += self.t;
            }
            SolvingMethod::IterationChebyshevsky => {
                (small_delta, big_delta) = self.get_delta(&solving_method);

                cheb_params = get_optimal_chebyshev_params(max_steps);
                cheb_coeffs =
                    calculate_chebyshev_coefficients(&cheb_params, small_delta, big_delta);
                cheb_coeffs_iter = cheb_coeffs.iter().cycle();
            }
            SolvingMethod::AlternatingTriangular => {
                (small_delta, big_delta) = self.get_delta(&solving_method);
                ksi = small_delta / big_delta;

                let h = (hx2 + hy2).sqrt();
                calculated_steps_count = ((1_f64 / eps).ln() / (2_f64 * PI * h)) as usize + 1;

                let w = 2_f64 / (small_delta * big_delta).sqrt();
                let g1 = small_delta / (2_f64 + 2_f64 * ksi.sqrt());
                let g2 = small_delta / (4_f64 * ksi.sqrt());

                tau = 2_f64 / (g1 + g2);

                av0 = vec![
                    w / hx2, // 0 - k1
                    w / hy2, // 1 - k2
                ];
                am = vec![vec![0_f64; m]; n];
            }
            SolvingMethod::AlternatingTriangularChebyshevsky => {
                (small_delta, big_delta) = self.get_delta(&solving_method);
                ksi = small_delta / big_delta;
                calculated_steps_count =
                    ((2_f64 / eps).ln() / (2_f64 * SQRT_2 * ksi.powf(0.25_f64))) as usize + 1;

                let w = 2_f64 / (small_delta * big_delta).sqrt();

                cheb_params = get_optimal_chebyshev_params(max_steps);

                av0 = vec![
                    w / hx2, // 0 - k1
                    w / hy2, // 1 - k2
                ];

                let g1 = small_delta / (2_f64 + 2_f64 * ksi.sqrt());
                let g2 = small_delta / (4_f64 * ksi.sqrt());

                cheb_coeffs = calculate_chebyshev_coefficients(&cheb_params, g1, g2);
                cheb_coeffs_iter = cheb_coeffs.iter().cycle();

                am = vec![vec![0_f64; m]; n];
            }
            SolvingMethod::AlternatingDirections => {
                (small_delta, big_delta) = self.get_delta(&solving_method);

                tau = 2_f64 / (small_delta * big_delta).sqrt();
                // let h = (hx2 + hy2).sqrt();
                // tau = h.powi(2) / (PI * h).sin();
                tau += self.t;

                calculated_steps_count = ((1_f64 / eps).ln() / (2_f64 * PI)) as usize + 1;

                avs0 = vec![0_f64; n];
                av0 = vec![0_f64; n];
                av01 = vec![0_f64; n];
                av02 = vec![0_f64; n];

                avs1 = vec![0_f64; m];
                av1 = vec![0_f64; m];
                av11 = vec![0_f64; m];
                av12 = vec![0_f64; m];

                avt0 = vec![0_f64; n - 2];
                avt1 = vec![0_f64; m - 2];

                // println!("mf\n{:?}", mf);
                // println!("mp\n{:?}", mp);
                // println!("mq\n{:?}", mq);

                match optimization {
                    OptimizationType::Parallelism => {
                        thread::spawn(move || {
                            parallel_raw_run_method_second(tx10, rx10, rxr10);
                        });
                    }
                    OptimizationType::Multithreading => {
                        for _ in 0..threads_count {
                            let (tx, rx) = mpsc::channel::<RawTrisf64Mut>();
                            let txue0 = txme0.clone();
                            let (txu, rxu) = mpsc::channel::<usize>();
                            thread::spawn(move || {
                                multithread_run_method(rx, txue0, rxu);
                            });
                            txmv0e.push(tx);
                            txmvu0e.push(txu);
                        }
                    }
                    _ => {}
                }
            }
        }
        calculated_steps_count *= 5;

        match calculation {
            CalculationType::MaxSteps => {
                self.steps_count = max_steps;
            }
            CalculationType::Calculated => {
                self.steps_count = calculated_steps_count;
            }
        }

        // Backtrace initializing
        if self.backtrace_level != BacktraceLevel::None {
            self.backtrace.start(
                SolverBacktraceCfg {
                    level: self.backtrace_level.clone(),
                    steps_count: self.steps_count,
                    u0_norm: self.get_discrepancy_zero_norm(&mf),
                    preserving_start_end: self.preserving_start_end,
                    division_middle: self.division_middle,
                },
                tau,
                cheb_params,
            );
            uk_prev = self.u.clone();
            ph1ph = (1_f64 - ksi) / (1_f64 + ksi);
            ph1ph = ph1ph / (1_f64 - ph1ph);
        }

        // _______________ Solving values
        let mut fij: f64;

        // tmp
        let mut ut: f64;
        let mut uq1v: VecDeque<f64>;

        // For implicit methods we use shifting
        if solving_method == SolvingMethod::SimpleIteration
            || solving_method == SolvingMethod::OptimalSimpleIteration
            || solving_method == SolvingMethod::IterationChebyshevsky
        {
            uq1v = VecDeque::from(vec![0_f64; n - 2]);
        } else {
            uq1v = VecDeque::new();
        }

        // [i][j] : i - row, j - column

        let mut p1: f64; // p row top
        let mut q1: f64; // q column left
        let mut p2: f64; // p row bottom
        let mut q2: f64; // q column right

        let mut up1: f64; // u row top
        let mut uq1: f64; // u column left
        let mut up2: f64; // u row bottom
        let mut uq2: f64; // u column right

        let mut cur_u;
        let mut cur_tau = tau;

        let backtrace_update = |eq: &mut Self,
                                cur_tau: f64,
                                mf: &[Vec<f64>],
                                mp: &[Vec<f64>],
                                mq: &[Vec<f64>],
                                ph1ph: f64,
                                uk_prev: &[Vec<f64>]| {
            if eq.backtrace_level != BacktraceLevel::None && eq.backtrace.check_current_step() {
                eq.backtrace.add_data(
                    cur_tau,
                    eq.get_discrepancy_norm(mf, mp, mq, &eq.u),
                    eq.get_norm_adj_difference(uk_prev),
                    ph1ph,
                    eq.u.clone(),
                );
            };
        };

        for _ in 0..self.steps_count {
            match solving_method {
                SolvingMethod::SimpleIteration | SolvingMethod::OptimalSimpleIteration => {
                    for i in 1..(n - 1) {
                        uq1v[i - 1] = self.u[i][0];
                    }
                }
                SolvingMethod::IterationChebyshevsky => {
                    for i in 1..(n - 1) {
                        uq1v[i - 1] = self.u[i][0];
                    }

                    // Iterator is cycled for this method and len > 0, so we can unwrap it safely.
                    cur_tau = *cheb_coeffs_iter.next().unwrap();
                }
                SolvingMethod::AlternatingTriangularChebyshevsky => {
                    // Iterator is cycled for this method and len > 0, so we can unwrap it safely.
                    cur_tau = *cheb_coeffs_iter.next().unwrap();
                }
                _ => {}
            }

            match solving_method {
                SolvingMethod::SimpleIteration => {
                    for j in 1..(m - 1) {
                        ut = self.u[0][j];

                        for i in 1..(n - 1) {
                            // Setting values
                            fij = mf[i - 1][j - 1];

                            p1 = mp[i - 1][j - 1];
                            q1 = mq[i - 1][j - 1];
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            up2 = self.u[i + 1][j];
                            uq2 = self.u[i][j + 1];

                            // Shifting
                            up1 = ut;
                            uq1 = uq1v.pop_front().unwrap();

                            // Preparing for shift
                            ut = self.u[i][j];
                            uq1v.push_back(self.u[i][j]);

                            // Calculating
                            cur_u =
                                ((p1 * up1 + p2 * up2) / hx2 + (q1 * uq1 + q2 * uq2) / hy2 + fij)
                                    / (p1 / hx2 + p2 / hx2 + q1 / hy2 + q2 / hy2);

                            // Writing
                            self.u[i][j] = cur_u;
                        }
                    }
                }
                SolvingMethod::OptimalSimpleIteration | SolvingMethod::IterationChebyshevsky => {
                    // Difference in optimal parameter
                    for j in 1..(m - 1) {
                        ut = self.u[0][j];

                        for i in 1..(n - 1) {
                            // Setting values
                            fij = mf[i - 1][j - 1];

                            p1 = mp[i - 1][j - 1];
                            q1 = mq[i - 1][j - 1];
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            up2 = self.u[i + 1][j];
                            uq2 = self.u[i][j + 1];

                            // Shifting
                            up1 = ut;
                            uq1 = uq1v.pop_front().unwrap();

                            // Preparing for shift
                            ut = self.u[i][j];
                            uq1v.push_back(self.u[i][j]);

                            // Calculating
                            cur_u = ut
                                + cur_tau
                                    * (p2 * (up2 - ut) / hx2 - p1 * (ut - up1) / hx2
                                        + q2 * (uq2 - ut) / hy2
                                        - q1 * (ut - uq1) / hy2
                                        + fij);

                            // Writing
                            self.u[i][j] = cur_u;
                        }
                    }
                }
                SolvingMethod::Zeidel => {
                    for j in 1..(m - 1) {
                        for i in 1..(n - 1) {
                            // Setting values
                            fij = mf[i - 1][j - 1];

                            p1 = mp[i - 1][j - 1];
                            q1 = mq[i - 1][j - 1];
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            up1 = self.u[i - 1][j];
                            uq1 = self.u[i][j - 1];
                            up2 = self.u[i + 1][j];
                            uq2 = self.u[i][j + 1];

                            // Calculating
                            cur_u =
                                ((p1 * up1 + p2 * up2) / hx2 + (q1 * uq1 + q2 * uq2) / hy2 + fij)
                                    / (p1 / hx2 + p2 / hx2 + q1 / hy2 + q2 / hy2);

                            // Writing
                            self.u[i][j] = cur_u;
                        }
                    }
                }
                SolvingMethod::UpperRelaxation => {
                    // Like Zeidel
                    for j in 1..(m - 1) {
                        for i in 1..(n - 1) {
                            // Setting values
                            fij = mf[i - 1][j - 1];

                            p1 = mp[i - 1][j - 1];
                            q1 = mq[i - 1][j - 1];
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            up2 = self.u[i + 1][j];
                            uq2 = self.u[i][j + 1];

                            // Shifting
                            up1 = self.u[i - 1][j];
                            uq1 = self.u[i][j - 1];

                            // Preparing for shift
                            ut = self.u[i][j];

                            // Calculating
                            cur_u = ut
                                + tau
                                    * (fij + p2 * (up2 - ut) / hx2 - p1 * (ut - up1) / hx2
                                        + q2 * (uq2 - ut) / hy2
                                        - q1 * (ut - uq1) / hy2)
                                    / (p1 / hx2 + p2 / hx2 + q1 / hy2 + q2 / hy2);

                            // Writing
                            self.u[i][j] = cur_u;
                        }
                    }
                }
                SolvingMethod::AlternatingTriangular
                | SolvingMethod::AlternatingTriangularChebyshevsky => {
                    let (k1, k2) = (av0[0], av0[1]);
                    for j in 1..(m - 1) {
                        for i in 1..(n - 1) {
                            // Setting values
                            fij = mf[i - 1][j - 1];

                            p1 = mp[i - 1][j - 1];
                            q1 = mq[i - 1][j - 1];
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            up2 = self.u[i + 1][j];
                            uq2 = self.u[i][j + 1];

                            // Shifting
                            up1 = self.u[i - 1][j];
                            uq1 = self.u[i][j - 1];

                            // Preparing for shift
                            ut = self.u[i][j];

                            // Calculating
                            cur_u = p2 * (up2 - ut) / hx2 - p1 * (ut - up1) / hx2
                                + q2 * (uq2 - ut) / hy2
                                - q1 * (ut - uq1) / hy2;
                            am[i][j] =
                                (k1 * p1 * am[i - 1][j] + k2 * q1 * am[i][j - 1] + cur_u + fij)
                                    / (1_f64 + k1 * p1 + k2 * q1);
                        }
                    }

                    for j in (1..(m - 1)).rev() {
                        for i in (1..(n - 1)).rev() {
                            // Setting values
                            p2 = mp[i][j - 1];
                            q2 = mq[i - 1][j];

                            // Calculating
                            am[i][j] = (k1 * p2 * am[i + 1][j] + k2 * q2 * am[i][j + 1] + am[i][j])
                                / (1_f64 + k1 * p2 + k2 * q2);
                            self.u[i][j] += cur_tau * am[i][j];
                        }
                    }
                }
                SolvingMethod::AlternatingDirections => {
                    for j in 0..(m - 2) {
                        // Filling vector b
                        for i in 0..(n - 2) {
                            uq1 = self.u[i + 1][j];
                            ut = self.u[i + 1][j + 1];
                            uq2 = self.u[i + 1][j + 2];
                            avs0[i + 1] = ut
                                + tau / 2_f64
                                    * (mf[i][j] + mq[i][j + 1] * (uq2 - ut) / hy2
                                        - mq[i][j] * (ut - uq1) / hy2);
                        }

                        avs0[0] = self.u[0][j + 1];
                        avs0[n - 1] = self.u[n - 1][j + 1];

                        // Filling tridiagonal matrix as 3 vectors
                        a1 = -tau / 2_f64 / hx2;
                        for t in 1..(n - 1) {
                            av0[t] = a1 * mp[t - 1][j];
                            av02[t] = a1 * mp[t][j];
                            av01[t] = 1_f64 - av0[t] - av02[t];
                        }

                        av01[0] = 1_f64;
                        av01[n - 1] = 1_f64;

                        // Calculating solution by passing method
                        match optimization {
                            OptimizationType::Nothing => {
                                run_method(&av0, &mut av01, &av02, &mut avs0);
                            }
                            OptimizationType::Parallelism => {
                                parallel_raw_run_method_main(
                                    &tx01, &rx01, &txr01, &mut av0, &mut av01, &mut av02, &mut avs0,
                                );
                            }
                            OptimizationType::Multithreading => {
                                multithread_run_method_main(
                                    &txmv0e, &rxm0e, &txmvu0e, &mut av0, &mut av01, &mut av02,
                                    &mut avs0,
                                );
                            }
                        }

                        // Using temporary vector to save memory
                        if j > 0 {
                            // Writing result for previous column
                            for (i, &avt0i) in avt0.iter().enumerate() {
                                self.u[i + 1][j] = avt0i;
                            }
                        }

                        if j < m - 3 {
                            // Storing result for this column
                            avt0.copy_from_slice(&avs0[1..(n - 1)]);
                        } else {
                            // Writing result for last column
                            for (i, &avs0i) in avs0.iter().enumerate().take(n - 1).skip(1) {
                                self.u[i][j + 1] = avs0i;
                            }
                        }
                    }

                    for i in 0..(n - 2) {
                        for j in 0..(m - 2) {
                            up1 = self.u[i][j + 1];
                            ut = self.u[i + 1][j + 1];
                            up2 = self.u[i + 2][j + 1];
                            avs1[j + 1] = ut
                                + tau / 2_f64
                                    * (mf[i][j] + mp[i + 1][j] * (up2 - ut) / hx2
                                        - mp[i][j] * (ut - up1) / hx2);
                        }

                        avs1[0] = self.u[i + 1][0];
                        avs1[m - 1] = self.u[i + 1][m - 1];

                        a1 = -tau / 2_f64 / hy2;
                        for t in 1..(m - 1) {
                            av1[t] = a1 * mq[i][t - 1];
                            av12[t] = a1 * mq[i][t];
                            av11[t] = 1_f64 - av1[t] - av12[t];
                        }

                        av11[0] = 1_f64;
                        av11[m - 1] = 1_f64;

                        match optimization {
                            OptimizationType::Nothing => {
                                run_method(&av1, &mut av11, &av12, &mut avs1);
                            }
                            OptimizationType::Parallelism => {
                                parallel_raw_run_method_main(
                                    &tx01, &rx01, &txr01, &mut av1, &mut av11, &mut av12, &mut avs1,
                                );
                            }
                            OptimizationType::Multithreading => {
                                multithread_run_method_main(
                                    &txmv0e, &rxm0e, &txmvu0e, &mut av1, &mut av11, &mut av12,
                                    &mut avs1,
                                );
                            }
                        }

                        if i > 0 {
                            self.u[i][1..(m - 1)].copy_from_slice(&avt1);
                        }

                        if i < n - 3 {
                            avt1.copy_from_slice(&avs1[1..(m - 1)]);
                        } else {
                            self.u[i + 1][1..(m - 1)].copy_from_slice(&avs1[1..(m - 1)]);
                        }
                    }
                }
            }

            backtrace_update(self, cur_tau, &mf, &mp, &mq, ph1ph, &uk_prev);
            uk_prev = self.u.clone();
        }
    }
}

fn multithread_run_method_main(
    txmv0e: &[Sender<RawTrisf64Mut>],
    rxm0e: &Receiver<usize>,
    txmvu0e: &[Sender<usize>],
    l: &mut [f64],
    d: &mut [f64],
    r: &mut [f64],
    b: &mut [f64],
) {
    let threads_count = txmv0e.len();
    let n = b.len();

    let part = n / threads_count;
    let mut u = n % part;
    let last_part: usize = if u == 0 {
        part
    } else {
        n - (threads_count - 1) * part
    };

    let mut tn: usize;
    let mut id: ThreadID;
    for (i, tx) in txmv0e.iter().enumerate() {
        match i {
            0 => {
                id = ThreadID::First;
                tn = part;
                u = 0;
            }
            i if i == threads_count - 1 => {
                id = ThreadID::Last;
                tn = last_part;
                u = (threads_count - 1) * part;
            }
            _ => {
                id = ThreadID::Other;
                tn = part;
                u = i * part;
            }
        }
        // println!("pointers for {} at {} with part {}", i, u, tn);
        unsafe {
            tx.send(RawTrisf64Mut {
                id, // 1 is first, 0 is last
                n: tn,
                l: l.as_mut_ptr().add(u),
                d: d.as_mut_ptr().add(u),
                r: r.as_mut_ptr().add(u),
                b: b.as_mut_ptr().add(u),
            })
            .unwrap(); // 1s
        }
    }

    for _ in 0..threads_count {
        rxm0e.recv().unwrap(); // 2r
    }

    for txu in txmvu0e {
        txu.send(n).unwrap(); // 3s
    }

    for _ in 0..threads_count {
        rxm0e.recv().unwrap(); // 4r
    }

    // Solving internal matrix using passing method
    let mut g: f64;

    // i = 1..n
    for i in ((2 * part - 1)..(n - last_part)).step_by(part) {
        g = l[i] / d[i - part]; // g = l[i] / d[i - 1]
        d[i] -= g * r[i - part]; // d[i] -= g * r[i - 1]
        b[i] -= g * b[i - part]; // b[i] -= g * b[i - 1]
    }

    g = l[n - 1] / d[n - 1 - last_part]; // g = l[i] / d[i - 1]
    d[n - 1] -= g * r[n - 1 - last_part]; // d[i] -= g * r[i - 1]
    b[n - 1] -= g * b[n - 1 - last_part]; // b[i] -= g * b[i - 1]

    b[n - 1] /= d[n - 1]; // b[n - 1] /= d[n - 1]

    // i = 0..(n-1)
    b[n - 1 - last_part] =
        (b[n - 1 - last_part] - r[n - 1 - last_part] * b[n - 1]) / d[n - 1 - last_part];

    for i in ((part - 1)..(n - last_part - part)).step_by(part).rev() {
        // b[i] = (b[i] - r[i] * b[i + 1]) / d[i];
        b[i] = (b[i] - r[i] * b[i + part]) / d[i];
    }

    for txu in txmvu0e {
        txu.send(n).unwrap(); // 5s
    }

    for _ in 0..threads_count {
        rxm0e.recv().unwrap(); // 6r
    }
}

fn multithread_run_method(rx: Receiver<RawTrisf64Mut>, txu: Sender<usize>, rxu: Receiver<usize>) {
    loop {
        // Receiving part of matrix described as thread id, size and 4 raw pointers
        let RawTrisf64Mut { id, n, l, d, r, b } = {
            let rec = rx.recv(); // 1r
            match rec {
                Ok(value) => value,
                Err(_) => {
                    break;
                }
            }
        };

        let mut g: f64;
        let h: f64;

        for i in 1..n {
            unsafe {
                g = *l.add(i) / *d.add(i - 1); // g = l[i] / d[i-1]
                *l.add(i) = -g * *l.add(i - 1); // l[i] = -g * l[i-1]
                *d.add(i) -= g * *r.add(i - 1); // d[i] -= g * r[i-1]
                *b.add(i) -= g * *b.add(i - 1); // b[i] -= g * b[i-1]
            }
        }

        for i in (1..(n - 1)).rev() {
            unsafe {
                g = *r.add(i - 1) / *d.add(i); // g = r[i-1] / d[i]
                *r.add(i - 1) = -g * *r.add(i); // r[i-1] = -g * r[i]
                *l.add(i - 1) -= g * *l.add(i); // l[i-1] -= g * l[i]
                *b.add(i - 1) -= g * *b.add(i); // b[i-1] -= g * b[i]
            }
        }

        txu.send(n).unwrap(); // 2s
        rxu.recv().unwrap(); // 3r

        if id != ThreadID::Last {
            unsafe {
                g = *r.add(n - 1) / *d.add(n); // g = r[n-1] / d[n]
                *d.add(n - 1) -= g * *l.add(n); // d[n-1] -= g * l[n]
                *r.add(n - 1) = -g * *r.add(n); // r[n-1] = -g * r[n]
                *b.add(n - 1) -= g * *b.add(n); // b[n-1] -= g * b[n]
            }
        }

        txu.send(n).unwrap(); // 4s
        rxu.recv().unwrap(); // 5r

        // Calculating result
        if id == ThreadID::First {
            unsafe {
                g = *b.add(n - 1);
                for i in 0..(n - 1) {
                    //println!("A{} n{}", i, n);
                    // b[i] = (b[i] - r[i] * b[n-1]) / d[i]
                    *b.add(i) = (*b.add(i) - *r.add(i) * g) / *d.add(i);
                }
            }
        } else {
            unsafe {
                h = *b.sub(1);
                g = *b.add(n - 1);
                for i in 0..(n - 1) {
                    //println!("B{} n{}", i, n);
                    // b[i] = (b[i] - l[i] * b[-1] - r[i] * b[n-1]) / d[i]
                    *b.add(i) = (*b.add(i) - *l.add(i) * h - *r.add(i) * g) / *d.add(i);
                }
            }
        }

        // Sending message that thread is ended work
        txu.send(n).unwrap(); // 6s
    }
}

fn parallel_raw_run_method_main(
    tx01: &Sender<usize>,
    rx01: &Receiver<usize>,
    txr01: &Sender<RawTrisf64>,
    l: &mut [f64],
    d: &mut [f64],
    r: &mut [f64],
    b: &mut [f64],
) {
    let n = b.len();
    let n2 = n / 2;

    tx01.send(n2).unwrap(); // 1
    txr01
        .send(RawTrisf64 {
            l: l.as_ptr(),
            d: d.as_mut_ptr(),
            r: r.as_ptr(),
            b: b.as_mut_ptr(),
        })
        .unwrap(); // 2

    let mut c;
    for i in ((n2 + 1)..n).rev() {
        c = r[i - 1] / d[i];
        d[i - 1] -= c * l[i];
        b[i - 1] -= c * b[i];
    }

    tx01.send(n).unwrap(); // 3
    rx01.recv().unwrap(); // 4
    d[n2] = 1_f64;

    for i in (n2 + 1)..n {
        b[i] = (b[i] - l[i] * b[i - 1]) / d[i];
    }

    rx01.recv().unwrap(); // 5
}

fn parallel_raw_run_method_second(
    tx10: Sender<usize>,
    rx10: Receiver<usize>,
    rxr10: Receiver<RawTrisf64>,
) {
    // Checking for new data
    loop {
        // Checking if main thread is closed
        let n2 = {
            let rec = rx10.recv(); // 1
            match rec {
                Ok(value) => value,
                Err(_) => {
                    break;
                }
            }
        };
        // println!("{}", n2);

        let RawTrisf64 { l, d, r, mut b } = rxr10.recv().unwrap(); // 2

        let mut c;

        // pointers at 0
        // c1 = l[i] / d[i - 1];
        // d[i] -= c1 * r[i - 1];
        // b[i] -= c1 * b[i - 1];
        for i in 1..n2 {
            unsafe {
                c = *l.add(i) / *d.add(i - 1);
                *d.add(i) -= c * *r.add(i - 1);
                c *= *b;
                b = b.add(1);
                *b -= c;
            }
        }
        // b n2-1

        rx10.recv().unwrap(); // 3

        // let c = l[n2] / d[n2 - 1];
        // d20 -= c * r[n2 - 1];
        // b20 = (b20 - c * b[n2 - 1]) / d20;
        // b[n2 - 1] = (b[n2 - 1] - r[n2 - 1] * b20) / d[n2 - 1];
        unsafe {
            c = *l.add(n2) / *d.add(n2 - 1);
            *d.add(n2) -= c * *r.add(n2 - 1);
            c *= *b;
            b = b.add(1);
            *b = (*b - c) / *d.add(n2);
            c = *r.add(n2 - 1) * *b;
        }

        tx10.send(n2).unwrap(); // 4

        unsafe {
            b = b.sub(1);
            *b = (*b - c) / *d.add(n2 - 1);
        }

        // b[i - 1] = (b[i - 1] - r[i - 1] * b[i]) / d[i - 1];
        for i in (1..n2).rev() {
            unsafe {
                c = *r.add(i - 1) * *b;
                b = b.sub(1);
                *b = (*b - c) / *d.add(i - 1);
            }
        }

        tx10.send(n2).unwrap(); // 5
    }
}

fn run_method(l: &[f64], d: &mut [f64], r: &[f64], b: &mut [f64]) {
    let n = b.len();

    let mut g: f64;
    for i in 1..n {
        g = l[i] / d[i - 1];
        d[i] -= g * r[i - 1];
        b[i] -= g * b[i - 1];
    }

    b[n - 1] /= d[n - 1];

    for i in (0..(n - 1)).rev() {
        b[i] = (b[i] - r[i] * b[i + 1]) / d[i];
    }

    //println!("a0\n{:?}\na1\n{:?}\na2\n{:?}\n", a0, a1, a2);
}

fn calculate_chebyshev_coefficients(
    cheb_params: &[usize],
    small_delta: f64,
    big_delta: f64,
) -> Vec<f64> {
    let mut cur_cheb = 0_usize;
    let mut result = Vec::<f64>::new();
    result.resize_with(cheb_params.len(), || {
        let value = 2_f64
            / (big_delta
                + small_delta
                + (big_delta - small_delta)
                    * (cheb_params[cur_cheb] as f64 / (2 * cheb_params.len()) as f64 * PI).cos());
        cur_cheb += 1;
        value
    });

    result
}

fn get_optimal_chebyshev_params(n: usize) -> Vec<usize> {
    match n {
        n if n % 16 == 0 => {
            vec![1, 31, 15, 17, 7, 25, 9, 23, 3, 29, 13, 19, 5, 27, 11, 21]
        }
        n if n % 8 == 0 => {
            vec![1, 15, 7, 9, 3, 13, 5, 11]
        }
        n if n % 4 == 0 => {
            vec![1, 7, 3, 5]
        }
        n if n % 2 == 0 => {
            vec![1, 3]
        }
        _ => {
            vec![1]
        }
    }
}
