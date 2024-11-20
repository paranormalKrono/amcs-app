use core::fmt;
use std::{f64::consts::PI, vec};

use crate::task6::*;
use backtrace::SolverBacktraceStep;
use elliptic_equation::EllipticEquation;

pub struct EE3;

impl EllipticEquation for EE3 {
    fn f(&self, x: f64, y: f64) -> f64 {
        2_f64 * x.sin() * y.cos()
    }

    fn p(&self, _x: f64, _y: f64) -> f64 {
        1_f64
    }
    fn q(&self, _x: f64, _y: f64) -> f64 {
        1_f64
    }
}

impl fmt::Display for EE3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let buf: String = "EE3: p = 1, q = 1, f = 2sin(x)cos(y)".to_string();

        write!(f, "{}", buf)
    }
}

fn mu(x: f64, y: f64) -> f64 {
    x.sin() * y.cos()
}

fn sol(x: f64, y: f64) -> f64 {
    x.sin() * y.cos()
}

pub fn get_task6() -> Option<EllipticEquationSolver<EE3>> {
    Some(EllipticEquationSolver::new(
        EE3 {},
        (PI, 1_f64),
        (1_f64, 1_f64),
        (1_f64, 1_f64),
        mu,
        (5, 5),
        true,
    ))
}

pub fn get_real_solution(solver: &EllipticEquationSolver<EE3>) -> Vec<Vec<f64>> {
    let (hx, hy) = solver.get_h();
    let (n, m) = solver.get_size();

    let mut solution = vec![vec![0_f64; m]; n];

    let mut fi = 0_f64;
    for (_, solutioni) in solution.iter_mut().enumerate().take(n) {
        let mut fj = 0_f64;
        for (_, solutionij) in solutioni.iter_mut().enumerate().take(m) {
            *solutionij = sol(fi * hx, fj * hy);
            fj += 1_f64;
        }
        fi += 1_f64;
    }

    solution
}

pub fn get_absolute_relative_errors(solver: &EllipticEquationSolver<EE3>) -> (Vec<f64>, Vec<f64>) {
    let backtrace_datas: &Vec<SolverBacktraceStep> =
        solver.backtrace.backtrace_data.solver_steps.as_ref();

    let mut absolute_errors = vec![0_f64; backtrace_datas.len()];
    let mut relative_errors = vec![0_f64; backtrace_datas.len()];

    let (hx, hy) = solver.get_h();
    let (n, m) = solver.get_size();

    let mut max_mod = 0_f64;

    let mut fi: f64;
    let mut fj: f64;
    for i in 1..(n - 1) {
        fi = i as f64;
        for j in 1..(m - 1) {
            fj = j as f64;

            max_mod = f64::max(max_mod, sol(fi * hx, fj * hy).abs());
        }
    }

    let u0_error = max_mod;

    for k in 0..backtrace_datas.len() {
        max_mod = 0_f64;
        for i in 1..(n - 1) {
            fi = i as f64;
            for j in 1..(m - 1) {
                fj = j as f64;

                max_mod = f64::max(
                    max_mod,
                    (backtrace_datas[k].u[i][j] - sol(fi * hx, fj * hy)).abs(),
                );
            }
        }
        absolute_errors[k] = max_mod;
        relative_errors[k] = max_mod / u0_error;
    }

    (absolute_errors, relative_errors)
}
