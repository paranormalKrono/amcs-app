use crate::{commands::Storage, error::Error, task6::*, task6_concrete::*};
use backtrace::{BacktraceLevel, SolverBacktraceData};
use std::sync::Mutex;
use std::time::Instant;
use tauri::http::request;
use tauri::State;

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task6CalculationRequest {
    solving_method: SolvingMethod,
    optimization: OptimizationType,
    t: f64,
    eps: f64,
    calculation: CalculationType,
    max_steps: usize,
    n: usize,
    m: usize,
    is_backtrace: bool,
    first_end_preserving: usize,
    middle_division: usize,
    backtrace_level: BacktraceLevel,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task6RequestAnswer {
    absolute_errors: Vec<f64>,
    relative_errors: Vec<f64>,
    real_solution: Vec<Vec<f64>>,
    solver: String,
    backtrace_data: SolverBacktraceData,
    delta_time: f64,
}

#[tauri::command(rename_all = "snake_case")]
pub async fn task6_request(
    request: Task6CalculationRequest,
    state_storage: State<'_, Mutex<Storage>>,
) -> Result<Task6RequestAnswer, Error> {
    let mut state = state_storage.lock()?;

    println!("Task6 request: {:?}\nBusy: {}", request, state.is_busy);
    if state.is_busy {
        return Ok(Task6RequestAnswer {
            absolute_errors: Vec::<f64>::new(),
            relative_errors: Vec::<f64>::new(),
            real_solution: Vec::<Vec<f64>>::new(),
            solver: String::new(),
            backtrace_data: SolverBacktraceData::new(),
            delta_time: 0_f64,
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
            request.is_backtrace,
            request.first_end_preserving,
            request.middle_division,
        );

        let instant = Instant::now();

        solver.solve_ees(
            request.solving_method,
            request.optimization,
            request.calculation,
            request.eps,
            request.max_steps,
            request.backtrace_level,
        );

        let delta_time = instant.elapsed().as_secs_f64();

        let mut absolute_errors: Vec<f64> = Vec::new();
        let mut relative_errors: Vec<f64> = Vec::new();
        let mut real_solution: Vec<Vec<f64>> = Vec::new();

        if request.is_backtrace {
            (absolute_errors, relative_errors) = get_absolute_relative_errors(solver);
            real_solution = get_real_solution(solver);
        }

        let str = solver.to_string();
        Task6RequestAnswer {
            absolute_errors,
            relative_errors,
            real_solution,
            solver: str,
            backtrace_data: solver.backtrace.backtrace_data.clone(),
            delta_time,
        }
    };

    state.is_busy = false;
    Ok(request_answer)
}
