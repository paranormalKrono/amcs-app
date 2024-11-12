use crate::{commands::Storage, error::Error, task1::*, task1_concrete::*};
use std::sync::Mutex;
use tauri::State;

#[derive(Debug, PartialEq, serde::Deserialize, serde::Serialize)]
pub enum Task1RequestType {
    Calculate,
    CalculateTable,
    CalculateSolution,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task1CalculationRequest {
    request_type: Task1RequestType,
    solving_method: SolvingMethod,
    a: f64,
    b: f64,
    divisions_count: usize,
    n: usize,
    t: f64,
    is_reserving: bool,
    is_reserving_lis: bool,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task1RequestAnswer {
    table_rows: Vec<TableRow>,
    sol_t: f64,
    state: String,
}

#[tauri::command(rename_all = "snake_case")]
pub async fn task1_request(
    request: Task1CalculationRequest,
    state_storage: State<'_, Mutex<Storage>>,
) -> Result<Task1RequestAnswer, Error> {
    let mut state = state_storage.lock()?;

    println!("Task1 request: {:?}\nBusy: {}", request, state.is_busy);

    let mut answer = Task1RequestAnswer {
        state: "Unknown".to_string(),
        table_rows: vec![TableRow::new(request.divisions_count)],
        sol_t: 0_f64,
    };

    if state.is_busy {
        return Ok(answer);
    }

    state.is_busy = true;

    if state.task1data.is_none() {
        state.task1data = get_task1();
    }

    let data = state.task1data.as_mut().unwrap();

    data.set_reserving(request.is_reserving, request.is_reserving_lis);
    data.set_span((request.a, request.b));
    data.set_size(request.n);

    match request.request_type {
        Task1RequestType::Calculate => {
            data.calculate(request.n);
        }
        Task1RequestType::CalculateTable => {
            answer.table_rows = calculate_table_rows(data, request.n, request.divisions_count);
        }
        Task1RequestType::CalculateSolution => {
            let sol = data.calculate_solution(request.n);
            answer.sol_t = data.calculate_solution_result(&sol, request.t);
        }
    }

    answer.state = data.to_string();

    state.is_busy = false;
    Ok(answer)
}
