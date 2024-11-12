use crate::task1::{LinearlyIndependentSystem, RitzData, RitzProblem};

use peroxide::fuga::*;
use peroxide::structure::polynomial::Calculus;
use std::{fmt, vec};

#[derive(Debug)]
pub struct System13 {
    // a1, a2, b1, b2
    pub boundary_conditions: [f64; 4],
}

impl fmt::Display for System13 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "system 13 boundary conditions {:?}",
            self.boundary_conditions
        )
    }
}

impl RitzProblem for System13 {
    // fn boundary_conditions(&self) -> [f64; 4] {
    //     self.boundary_conditions
    // }

    fn fn_p(&self, t: f64) -> f64 {
        1_f64 / (2_f64 * t + 3_f64)
    }

    fn fn_r(&self, t: f64) -> f64 {
        1_f64 + t.cos()
    }

    fn fn_f(&self, t: f64) -> f64 {
        1_f64 + t
    }
}

impl System13 {
    fn fn_pd(&self, t: f64) -> f64 {
        -2_f64 / (2_f64 * t + 3_f64).powi(2)
    }
}

pub fn get_task1() -> Option<RitzData<System13, WFunctionSystem>> {
    let ritz_data = RitzData::new(
        System13 {
            boundary_conditions: [1., 0., 0., 1.],
        },
        WFunctionSystem::new(),
        true,
        true,
    );
    Some(ritz_data)
}

pub struct WFunctionSystem {
    wfunctions: Vec<Polynomial>,
    wfunctions_derivative: Vec<Polynomial>,
    is_reserving: bool,
    calculated_size: usize,
    cur_size: usize,
}

impl WFunctionSystem {
    pub fn new() -> WFunctionSystem {
        WFunctionSystem {
            wfunctions: vec![
                Polynomial {
                    coef: vec![1_f64, 1_f64],
                },
                Polynomial {
                    coef: vec![-1_f64, 0_f64, 1_f64],
                },
            ],
            wfunctions_derivative: vec![
                Polynomial { coef: vec![1_f64] },
                Polynomial {
                    coef: vec![-2_f64, 0_f64],
                },
            ],
            calculated_size: 2,
            cur_size: 2,
            is_reserving: true,
        }
    }
}

impl fmt::Display for WFunctionSystem {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut buf: String = "W Functions system:\n".to_string();
        // Write strictly the first element into the supplied output
        // stream: `f`. Returns `fmt::Result` which indicates whether the
        // operation succeeded or failed. Note that `write!` uses syntax which
        // is very similar to `println!`.
        for i in 0..self.wfunctions.len() {
            buf = format!(
                "{buf}W{i}:\n\t{}\n\t{}\n\n",
                self.wfunctions[i], self.wfunctions_derivative[i]
            );
        }
        write!(f, "{}", buf)
    }
}

fn next_jacobi_polynomial(jp0: Polynomial, jp1: &Polynomial, n: f64) -> Polynomial {
    let mut jp2 = (2_f64 * n + 5_f64) / (n + 2_f64) * jp1.clone();
    jp2.coef.push(0_f64);
    jp2 = (n + 3_f64) / (n + 4_f64) * (jp2 - jp0);
    jp2
}

impl LinearlyIndependentSystem for WFunctionSystem {
    fn set_size(&mut self, size: usize) {
        if size < 2 {
            self.cur_size = 2;
        } else {
            self.cur_size = size;
        }

        if !self.is_reserving && self.cur_size < self.calculated_size {
            self.wfunctions.truncate(self.cur_size);
            self.wfunctions_derivative.truncate(self.cur_size);
            self.calculated_size = self.cur_size;
        }
    }

    fn set_reserving(&mut self, is_reserving: bool) {
        self.is_reserving = is_reserving;

        self.set_size(self.cur_size);
    }

    fn calculate(&mut self, size: usize) {
        if self.calculated_size > size {
            return;
        }
        self.set_size(size);

        // Calculating Jacobi poliniomials from 2 to calculated_n - 1
        let mut jp0 = Polynomial { coef: vec![1_f64] };
        let mut jp1 = Polynomial {
            coef: vec![2_f64, 0_f64],
        };
        let p = Polynomial {
            coef: vec![-1_f64, 0_f64, 1_f64],
        };

        for i in 0..(self.calculated_size - 2) {
            let jp2 = next_jacobi_polynomial(jp0, &jp1, i as f64);
            jp0 = jp1;
            jp1 = jp2;
        }

        // Calculating additional W functions and derivatives from cur n to n - 1
        self.wfunctions
            .reserve(self.cur_size - self.calculated_size);
        self.wfunctions_derivative
            .reserve(self.cur_size - self.calculated_size);

        for i in (self.calculated_size - 2)..(self.cur_size - 2) {
            let jp2 = next_jacobi_polynomial(jp0, &jp1, i as f64);
            jp0 = jp1;
            jp1 = jp2;

            let wp = jp0.clone() * p.clone();
            self.wfunctions.push(wp.clone());
            self.wfunctions_derivative.push(wp.derivative());
        }

        self.calculated_size = self.cur_size;
    }

    fn fn_w(&self, i: usize, t: f64) -> f64 {
        self.wfunctions[i].eval(t)
    }

    fn fn_wd(&self, i: usize, t: f64) -> f64 {
        self.wfunctions_derivative[i].eval(t)
    }

    fn fn_ww(&self, i: usize, j: usize, t: f64) -> f64 {
        self.wfunctions[i].eval(t) * self.wfunctions[j].eval(t)
    }

    fn fn_wwd(&self, i: usize, j: usize, t: f64) -> f64 {
        self.wfunctions_derivative[i].eval(t) * self.wfunctions_derivative[j].eval(t)
    }
}

pub fn get_real_solution(solver: &RitzData<System13, WFunctionSystem>, n: usize) -> Vec<f64> {
    let (start, h) = {
        let (a, b) = solver.get_span();
        (a, (b - a) / ((n - 1) as f64))
    };

    let mut a = py_matrix(vec![vec![0_f64; n]; n]);
    let mut fk = vec![0_f64; n];
    let mut pk = vec![0_f64; n];
    let mut rk = vec![0_f64; n];
    let mut pdk = vec![0_f64; n];

    let mut t = start + h;

    for i in 1..(n - 1) {
        fk[i] = solver.problem.fn_f(t);
        pk[i] = solver.problem.fn_p(t);
        rk[i] = solver.problem.fn_r(t);
        pdk[i] = solver.problem.fn_pd(t);
        t += h;
    }

    fk[0] = 0_f64;
    fk[n - 1] = 0.9307495259772893_f64; // 0_f64;

    a[(0, 0)] = 1_f64;
    // a[(n - 1, n - 2)] = -1_f64 / h;
    a[(n - 1, n - 1)] = 1_f64; // = 1_f64 / h;

    let h2 = h.powi(2);
    for i in 1..(n - 1) {
        a[(i, i - 1)] = -pk[i] / h2 + pdk[i] / h;
        a[(i, i)] = 2_f64 * pk[i] / h2 + rk[i];
        a[(i, i + 1)] = -pk[i] / h2 - pdk[i] / h;
    }

    a.solve(&fk, SolveKind::LU)
}

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TableRow {
    pub condition_number: f64,
    pub solution: Vec<f64>,
    pub discrepancy: Vec<f64>,
}

impl TableRow {
    pub fn new(divisions_count: usize) -> TableRow {
        TableRow {
            condition_number: 0_f64,
            solution: vec![0_f64; divisions_count],
            discrepancy: vec![0_f64; divisions_count],
        }
    }
}

#[derive(Debug, serde::Serialize, serde::Deserialize)]
pub struct Table {
    pub data: Vec<TableRow>,
}

impl Table {
    pub fn new(count: usize, divisions_count: usize) -> Table {
        Table {
            data: vec![TableRow::new(divisions_count); count],
        }
    }
}

pub fn calculate_table_rows(
    solver: &mut RitzData<System13, WFunctionSystem>,
    n: usize,
    divisions_count: usize,
) -> Vec<TableRow> {
    solver.calculate(n);

    let mut table = Table::new(n, divisions_count);

    let data: &mut Vec<TableRow> = table.data.as_mut();

    let a = solver.get_matrix_a();

    let (start, hx) = {
        let (a, b) = solver.get_span();
        (a, (b - a) / ((divisions_count - 1) as f64))
    };

    let mut t: f64;
    let real_sol = get_real_solution(solver, divisions_count);
    println!("Real solution:\n\t{:?}", real_sol);

    let mut cur_a: Matrix;
    let mut cur_ai: Matrix;
    let mut ritz_solution: Vec<f64>;
    let mut cur_data: f64;
    for (i, datai) in data.iter_mut().enumerate().take(divisions_count) {
        cur_a = a.submat((0, 0), (i, i));
        cur_ai = cur_a.inv();
        datai.condition_number = cur_a.norm(Norm::L2) * cur_ai.norm(Norm::L2);

        ritz_solution = solver.calculate_solution(i + 1);

        t = start;
        for (j, &real_solj) in real_sol.iter().enumerate().take(divisions_count) {
            cur_data = solver.calculate_solution_result(&ritz_solution, t);
            datai.solution[j] = cur_data;
            datai.discrepancy[j] = (cur_data - real_solj).abs();

            t += hx;
        }
    }

    table.data
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn solve_ritz_test() -> anyhow::Result<()> {
        let n = 7;
        let boundary_conditions = [1., 0., 0., 1.];
        let span = (-1_f64, 1_f64);

        let system = System13 {
            boundary_conditions,
        };
        let functions_system = WFunctionSystem::new();
        let mut ritz_data = RitzData::new(system, functions_system, true, true);

        ritz_data.set_span(span);
        ritz_data.set_size(n);
        let sol = ritz_data.calculate_solution(n);

        let solution_in_m05 = ritz_data.calculate_solution_result(&sol, -0.5_f64);
        let solution_in_0 = ritz_data.calculate_solution_result(&sol, 0_f64);
        let solution_in_05 = ritz_data.calculate_solution_result(&sol, 0.5_f64);

        println!(
            "Result:\n -0.5 \t {}\n 0 \t {}\n 0.5 \t {}",
            solution_in_m05, solution_in_0, solution_in_05
        );

        println!(
            "Check boundary condition: u({}) = {} = 0",
            span.0,
            ritz_data.calculate_solution_result(&sol, span.0)
        );

        Ok(())
    }
}
