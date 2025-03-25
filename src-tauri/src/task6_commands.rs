use crate::{commands::Storage, error::Error, task6::*, task6_concrete::*};
use backtrace::{BacktraceLevel, SolverBacktraceData};
use equations::{MevalEquation, MevalEquationCfg};
use std::sync::Mutex;
use std::time::Instant;
use tauri::State;

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task6CalculationRequest {
    function_p: String,
    function_q: String,
    function_c: String,
    function_f: String,
    function_mu: String,
    function_sol: String,

    area_x: f64,
    area_y: f64,

    function_p_bound_left: f64,
    function_p_bound_right: f64,
    function_q_bound_left: f64,
    function_q_bound_right: f64,

    solving_method: SolvingMethod,
    n: usize,
    m: usize,
    t: f64,

    calculation: CalculationType,
    max_steps: usize,
    eps: f64,

    optimization: OptimizationType,
    threads_count: usize,

    backtrace_level: BacktraceLevel,
    first_end_preserving: usize,
    middle_division: usize,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task6RequestAnswer {
    real_solution: Vec<Vec<f64>>,
    solver: String,
    solution: Vec<Vec<f64>>,
    backtrace_data: SolverBacktraceData,
    delta_time: f64,

    steps_absolute_errors: Vec<f64>,
    steps_relative_errors: Vec<f64>,
    absolute_errors: Vec<Vec<f64>>,
    relative_errors: Vec<Vec<f64>>,
}

#[tauri::command(rename_all = "snake_case")]
pub async fn task6_request(
    request: Task6CalculationRequest,
    state_storage: State<'_, Mutex<Storage>>,
) -> Result<Task6RequestAnswer, Error> {
    let mut state = state_storage.lock()?;

    // println!("Task6 request: {:?}\nBusy: {}", request, state.is_busy);
    if state.is_busy {
        return Ok(Task6RequestAnswer {
            real_solution: Vec::<Vec<f64>>::new(),
            solver: String::new(),
            solution: Vec::new(),
            backtrace_data: SolverBacktraceData::new(),
            delta_time: 0_f64,
            steps_absolute_errors: Vec::<f64>::new(),
            steps_relative_errors: Vec::<f64>::new(),
            absolute_errors: Vec::<Vec<f64>>::new(),
            relative_errors: Vec::<Vec<f64>>::new(),
        });
    }
    state.is_busy = true;

    if state.task6data.is_none() {
        state.task6data = Some(EllipticEquationSolver::new());
    }

    let request_answer = {
        let solver = state.task6data.as_mut().unwrap();

        solver.set_param(request.t);
        solver.set_boundaries(
            request.function_p_bound_left,
            request.function_p_bound_right,
            request.function_q_bound_left,
            request.function_q_bound_right,
        );
        solver.set_size(request.n, request.m);
        solver.set_area(request.area_x, request.area_y);
        solver.set_backtrace(
            request.backtrace_level.clone(),
            request.first_end_preserving,
            request.middle_division,
        );

        let eeq = MevalEquation::new(MevalEquationCfg {
            p: request.function_p,
            q: request.function_q,
            c: request.function_c,
            f: request.function_f,
            mu: request.function_mu,
            sol: request.function_sol,
        });

        let solving_cfg = SolvingCfg {
            solving_method: request.solving_method,
            calculation: request.calculation,
            eps: request.eps,
            max_steps: request.max_steps,
            optimization: request.optimization,
            threads_count: request.threads_count,
        };

        let instant = Instant::now();
        solver.solve_ees(&eeq, solving_cfg);
        let delta_time = instant.elapsed().as_secs_f64();

        let mut u = Vec::new();
        let mut real_solution: Vec<Vec<f64>> = Vec::new();

        let mut steps_absolute_errors: Vec<f64> = Vec::new();
        let mut steps_relative_errors: Vec<f64> = Vec::new();
        let mut absolute_errors: Vec<Vec<f64>> = Vec::new();
        let mut relative_errors: Vec<Vec<f64>> = Vec::new();

        if request.backtrace_level != BacktraceLevel::None {
            u = solver.u.clone();
            let sol = |x, y| eeq.sol(x, y);
            real_solution = get_real_solution(solver, sol);
        }
        if request.backtrace_level == BacktraceLevel::Full
            || request.backtrace_level == BacktraceLevel::Graphs
        {
            (steps_absolute_errors, steps_relative_errors) =
                get_steps_absolute_relative_errors(solver, &real_solution);

            (absolute_errors, relative_errors) =
                get_absolute_relative_errors(&eeq, solver, &real_solution);
        }

        let str = solver.to_string();
        Task6RequestAnswer {
            real_solution,
            solver: str,
            solution: u,
            backtrace_data: solver.backtrace.backtrace_data.clone(),
            delta_time,
            steps_absolute_errors,
            steps_relative_errors,
            absolute_errors,
            relative_errors,
        }
    };

    state.is_busy = false;
    Ok(request_answer)
}
