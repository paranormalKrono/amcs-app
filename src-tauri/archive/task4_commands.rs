use crate::{commands::Storage, error::Error, task4_concrete::*};
use std::sync::Mutex;
use tauri::State;

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Task4CalculationRequest {
    a: f64,
    b: f64,
    n: usize,
    m: usize,
    mu: f64,
}

#[derive(Debug, Default, serde::Deserialize, serde::Serialize)]
pub struct Task4RequestAnswer {
    data_file_name: String,
    sols: Vec<Solution>,
    kappas: Vec<Vec<f64>>,
}

#[tauri::command(rename_all = "snake_case")]
pub async fn task4_request(
    request: Task4CalculationRequest,
    state_storage: State<'_, Mutex<Storage>>,
) -> Result<Task4RequestAnswer, Error> {
    let mut state = state_storage.lock()?;

    println!("Task4 request: {:?}\nBusy: {}", request, state.is_busy);
    if state.is_busy {
        return Ok(Default::default());
    }
    state.is_busy = true;

    let sols = calculate_all(request.a, request.b, request.n, request.m, request.mu);

    // println!("{:?}", sols);

    let data_file_name = "data".to_owned();
    let request_answer = Task4RequestAnswer {
        data_file_name: data_file_name.clone(),
        sols: sols.sols.clone(),
        kappas: sols.kappas.clone(),
    };

    state.is_busy = false;
    Ok(request_answer)
}
