// Learn more about Tauri commands at https://tauri.app/v1/guides/features/command

use std::sync::Mutex;

mod commands;
mod error;
// mod task1;
// mod task1_commands;
// mod task1_concrete;
// mod task4_commands;
// mod task4_concrete;
mod task6;
mod task6_commands;
mod task6_concrete;

use tauri::generate_handler;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    println!("Hello from Rust!");
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .manage(Mutex::new(commands::Storage::default()))
        .invoke_handler(generate_handler![
            // task1_commands::task1_request,
            // task4_commands::task4_request,
            task6_commands::task6_request
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
