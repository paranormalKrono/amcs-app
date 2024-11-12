use std::{fmt, vec};

use peroxide::fuga::*;

//
pub trait RitzProblem {
    // fn boundary_conditions(&self) -> [f64; 4];
    fn fn_p(&self, t: f64) -> f64;
    fn fn_r(&self, t: f64) -> f64;
    fn fn_f(&self, t: f64) -> f64;
}

pub trait LinearlyIndependentSystem {
    fn set_reserving(&mut self, is_reserving: bool);
    fn set_size(&mut self, n: usize);
    fn calculate(&mut self, size: usize);
    fn fn_w(&self, i: usize, t: f64) -> f64;
    fn fn_wd(&self, i: usize, t: f64) -> f64;
    fn fn_ww(&self, i: usize, j: usize, t: f64) -> f64;
    fn fn_wwd(&self, i: usize, j: usize, t: f64) -> f64;
}

#[derive(PartialEq, Debug, serde::Serialize, serde::Deserialize)]
pub enum SolvingMethod {
    Ritz,
    Relaxation,
}

pub struct RitzData<T: RitzProblem, F: LinearlyIndependentSystem> {
    pub problem: T,
    pub functions_system: F,
    t_span: (f64, f64),
    is_reserving: bool,
    solving_method: SolvingMethod,

    a: Matrix,
    f: Vec<f64>,
    calculated_size: usize,
    cur_size: usize,
}

impl<T: RitzProblem + std::fmt::Display, F: LinearlyIndependentSystem + std::fmt::Display>
    fmt::Display for RitzData<T, F>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let buf: String = format!(
            "Problem:\n\t{}\n
        LIS:\n\t{}\n
        t_span:\n\t{:?}\n
        Is reserving:\n\t{}\n
        Solving method:\n\t{:?}\n
        Matrix a:\n\n{}\n\n
        Vector f:\n\t{:?}\n
        Calculated size:\n\t{}\n
        Current size:\n\t{}\n",
            self.problem,
            self.functions_system,
            self.t_span,
            self.is_reserving,
            self.solving_method,
            self.a,
            self.f,
            self.calculated_size,
            self.cur_size,
        );

        write!(f, "{}", buf)
    }
}

impl<T: RitzProblem, F: LinearlyIndependentSystem> RitzData<T, F> {
    pub fn new(
        ritz_problem: T,
        mut lis_functions_system: F,
        is_reserving: bool,
        is_reserving_lis: bool,
    ) -> Self {
        let mut v: Vec<Vec<f64>> = Vec::<Vec<f64>>::new();
        let v0 = vec![0_f64; 1];
        v.resize(1, v0.clone());
        let a = py_matrix(v);

        lis_functions_system.set_reserving(is_reserving_lis);

        RitzData {
            problem: ritz_problem,
            functions_system: lis_functions_system,
            t_span: (0_f64, 0_f64),
            is_reserving,
            solving_method: SolvingMethod::Ritz,
            a,
            f: v0,
            calculated_size: 0,
            cur_size: 1,
        }
    }

    pub fn get_matrix_a(&self) -> Matrix {
        self.a.clone()
    }

    pub fn get_span(&self) -> (f64, f64) {
        self.t_span
    }

    pub fn set_span(&mut self, t_span: (f64, f64)) {
        if t_span != self.t_span {
            self.t_span = t_span;
            self.calculated_size = 0;
        }
    }

    pub fn set_size(&mut self, n: usize) {
        if n < 1 {
            return;
        }

        match n {
            n if n > self.cur_size => {
                let mut new_m = py_matrix(vec![vec![0_f64; n]; n]);
                for i in 0..self.cur_size {
                    for j in 0..self.cur_size {
                        new_m[(i, j)] = self.a[(i, j)];
                    }
                }
                self.a = new_m;

                self.f.resize(n, 0_f64);

                self.cur_size = n;
            }
            n if !self.is_reserving && n < self.cur_size => {
                self.a = self.a.submat((0, 0), (n - 1, n - 1));

                self.f.truncate(n);

                self.cur_size = n;
                if self.calculated_size > n {
                    self.calculated_size = n;
                }
            }
            _ => {}
        }

        self.functions_system.set_size(n);
    }

    pub fn set_reserving(&mut self, is_reserving: bool, is_reserving_lis: bool) {
        self.is_reserving = is_reserving;
        self.functions_system.set_reserving(is_reserving_lis);

        if !is_reserving {
            self.set_size(self.cur_size);
        }
    }

    pub fn calculate(&mut self, n: usize) {
        self.set_size(n);
        self.functions_system.calculate(n);

        let integrator = Integral::G30K61(1e-3, 30);

        // Using A matrix symmetry
        // Upper right block
        for i in 0..self.calculated_size {
            for j in self.calculated_size..n {
                self.a[(i, j)] = integrate(
                    |t| {
                        self.problem.fn_p(t) * self.functions_system.fn_wwd(i, j, t)
                            + self.problem.fn_r(t) * self.functions_system.fn_ww(i, j, t)
                    },
                    self.t_span,
                    integrator,
                );

                // A matrix symmetry
                self.a[(j, i)] = self.a[(i, j)];
            }
        }

        // Down right part
        for i in self.calculated_size..n {
            for j in (i + 1)..n {
                self.a[(i, j)] = integrate(
                    |t| {
                        self.problem.fn_p(t) * self.functions_system.fn_wwd(i, j, t)
                            + self.problem.fn_r(t) * self.functions_system.fn_ww(i, j, t)
                    },
                    self.t_span,
                    integrator,
                );

                // A matrix symmetry
                self.a[(j, i)] = self.a[(i, j)];
            }
        }

        // Diagonal
        for i in self.calculated_size..n {
            self.a[(i, i)] = integrate(
                |t| {
                    self.problem.fn_p(t) * self.functions_system.fn_wd(i, t).powi(2)
                        + self.problem.fn_r(t) * self.functions_system.fn_w(i, t).powi(2)
                },
                self.t_span,
                integrator,
            );
        }

        // Calculating f
        for i in self.calculated_size..n {
            self.f[i] = integrate(
                |t| self.problem.fn_f(t) * self.functions_system.fn_w(i, t),
                self.t_span,
                integrator,
            );
        }

        self.calculated_size = n;
    }

    pub fn calculate_solution(&mut self, n: usize) -> Vec<f64> {
        match n {
            size if size > self.calculated_size => {
                self.calculate(size);
                self.a.solve(&self.f, LU)
            }
            size if size == self.calculated_size => self.a.solve(&self.f, LU),
            _ => {
                let fcalc = &self.f[..n].to_vec();
                let acalc = self.a.submat((0, 0), (n - 1, n - 1));
                acalc.solve(fcalc, LU)
            }
        }
    }

    pub fn calculate_solution_result(&mut self, sol: &[f64], t: f64) -> f64 {
        let mut res = 0_f64;

        for (i, &soli) in sol.iter().enumerate() {
            res += self.functions_system.fn_w(i, t) * soli;
        }

        res
    }
}
