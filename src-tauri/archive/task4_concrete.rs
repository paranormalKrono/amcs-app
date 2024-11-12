use serde::{Deserialize, Serialize};
use std::vec;

#[derive(Debug, Clone, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum SolvingMethod {
    Explicit,
    PureImplicit,
    Implicit,
    Symmetric,
}

#[derive(Debug, Clone)]
pub struct AllSolutions {
    pub sols: Vec<Solution>,
    pub kappas: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Solution {
    pub solving_method: SolvingMethod,
    pub u: Vec<Vec<f64>>,
}

pub fn calculate_all(a: f64, b: f64, n: usize, m: usize, mu: f64) -> AllSolutions {
    AllSolutions {
        sols: vec![
            calculate(a, b, n, m, mu, SolvingMethod::Explicit),
            calculate(a, b, n, m, mu, SolvingMethod::PureImplicit),
            calculate(a, b, n, m, mu, SolvingMethod::Implicit),
            calculate(a, b, n, m, mu, SolvingMethod::Symmetric),
        ],
        kappas: calculate_kappas(a, b, n, m, mu),
    }
}

fn calculate_kappas(a: f64, b: f64, n: usize, m: usize, mu: f64) -> Vec<Vec<f64>> {
    let (h, tau) = (a / ((n - 1) as f64), b / ((m - 1) as f64));
    let th = tau / h;

    // it's transposed, don't ask me why or I'll go insane
    let mut kappas = vec![vec![0_f64; n]; m];

    let mut cur_x: f64;
    for kappasj in kappas.iter_mut().take(m) {
        cur_x = 1_f64;
        for kappasji in kappasj.iter_mut().take(n) {
            *kappasji = -mu * (1_f64 - cur_x) * th;
            cur_x += h;
        }
    }

    kappas
}

fn calculate(a: f64, b: f64, n: usize, m: usize, mu: f64, method: SolvingMethod) -> Solution {
    let mut u = vec![vec![0_f64; m]; n];

    let (h, tau) = (a / ((n - 1) as f64), b / ((m - 1) as f64));
    let th = tau / h;

    let mut kappa;

    // u(x,0) = x^2
    let mut cur_x = 1_f64;
    for ui in u.iter_mut().take(n) {
        ui[0] = cur_x.powi(2);
        cur_x += h;
    }

    // u(1,t) = 1
    for j in 0..m {
        u[0][j] = 1_f64;
    }

    // Calculating
    match method {
        SolvingMethod::Explicit => {
            for j in 1..m {
                cur_x = 1_f64;
                for i in 1..n {
                    cur_x += h;
                    kappa = -mu * (1_f64 - cur_x) * th;

                    u[i][j] = (1_f64 - kappa) * u[i][j - 1] + kappa * u[i - 1][j - 1];
                }
            }
        }
        SolvingMethod::PureImplicit => {
            for j in 1..m {
                cur_x = 1_f64;
                for i in 1..n {
                    cur_x += h;
                    kappa = -1_f64 / th / (1_f64 - cur_x) / mu;

                    u[i][j] = (u[i][j - 1] * kappa + u[i - 1][j]) / (1_f64 + kappa);
                }
            }
        }
        SolvingMethod::Implicit => {
            for j in 1..m {
                cur_x = 1_f64;
                for i in 1..n {
                    cur_x += h;
                    kappa = -1_f64 / th / (1_f64 - cur_x) / mu;

                    u[i][j] = (1_f64 - kappa) * u[i - 1][j] + kappa * u[i - 1][j - 1];
                }
            }
        }
        SolvingMethod::Symmetric => {
            for j in 1..m {
                cur_x = 1_f64;
                for i in 1..n {
                    cur_x += h;
                    kappa = -mu * (1_f64 - cur_x) * th;

                    u[i][j] = (-u[i - 1][j]
                        + u[i][j - 1]
                        + u[i - 1][j - 1]
                        + kappa * (-u[i][j - 1] + u[i - 1][j] + u[i - 1][j - 1]))
                        / (1_f64 + kappa);
                }
            }
        }
    }

    Solution {
        solving_method: method,
        u,
    }
}
