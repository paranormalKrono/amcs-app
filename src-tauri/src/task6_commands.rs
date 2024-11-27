use crate::{commands::Storage, error::Error, task6::*, task6_concrete::*};
use backtrace::{BacktraceLevel, SolverBacktraceData};
use std::sync::Mutex;
use std::time::Instant;
use tauri::State;

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task6CalculationRequest {
    solving_method: SolvingMethod,
    threads_count: usize,
    t: f64,
    eps: f64,
    calculation: CalculationType,
    max_steps: usize,
    n: usize,
    m: usize,
    optimization: OptimizationType,
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
        state.task6data = get_task6();
    }

    let request_answer = {
        let solver = state.task6data.as_mut().unwrap();

        solver.set_param(request.t);
        solver.set_size(request.n, request.m);
        solver.set_backtrace(
            request.backtrace_level.clone(),
            request.first_end_preserving,
            request.middle_division,
        );

        let instant = Instant::now();

        solver.solve_ees(
            request.solving_method,
            request.calculation,
            request.eps,
            request.max_steps,
            request.optimization,
            request.threads_count,
        );

        let delta_time = instant.elapsed().as_secs_f64();

        let mut real_solution: Vec<Vec<f64>> = Vec::new();

        let mut steps_absolute_errors: Vec<f64> = Vec::new();
        let mut steps_relative_errors: Vec<f64> = Vec::new();
        let mut absolute_errors: Vec<Vec<f64>> = Vec::new();
        let mut relative_errors: Vec<Vec<f64>> = Vec::new();

        if request.backtrace_level != BacktraceLevel::None {
            real_solution = get_real_solution(solver);
        }
        if request.backtrace_level == BacktraceLevel::Full
            || request.backtrace_level == BacktraceLevel::Graphs
        {
            (steps_absolute_errors, steps_relative_errors) =
                get_steps_absolute_relative_errors(solver, &real_solution);

            (absolute_errors, relative_errors) =
                get_absolute_relative_errors(solver, &real_solution);
        }

        let u = solver.u.clone();
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
