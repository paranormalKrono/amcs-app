


import { invoke } from '@tauri-apps/api/core';

import Plotly from 'plotly.js-dist-min';


const header = document.getElementById("header-element") as HTMLElement;
const main_form = document.getElementById("main_form") as HTMLFormElement;
const data_text1 = document.getElementById("data_text1") as HTMLElement;
const data_text2 = document.getElementById("data_text2") as HTMLElement;
const data_output = document.getElementById("data_output") as HTMLElement;

const button_change_task = document.getElementById("button_change_task") as HTMLButtonElement;
const button_send_request = document.getElementById("button_send_request") as HTMLButtonElement;


enum BacktraceLevel {
    None = "None",
    Tables = "Tables",
    Graphs = "Graphs",
    Full = "Full",
}

type request6 = {
    solving_method: string,
    t: number,
    eps: number,
    calculation: string,
    max_steps: number,
    n: number,
    m: number,
    optimization: string,
    threads_count: number,
    backtrace_level: BacktraceLevel,
    first_end_preserving: number,
    middle_division: number,
}

type RequestType = request6;
let request: RequestType;





type data6backtracestep = {
    step: number,
    discrepancy_norm: number,
    relative_discrepancy: number,
    norm_adj_difference: number,
    apposterior_est: number,

    spectral_radius_approach: number,
    u: number[][],
}

type data6solbacktracedata = {
    solver_steps: data6backtracestep[],

    tau: number,
    cheb_params: number[],
    taus: number[],
}

type data6 = {
    solver: string,
    real_solution: number[][],
    solution: number[][],
    backtrace_data: data6solbacktracedata,
    delta_time: number,
    steps_absolute_errors: number[],
    steps_relative_errors: number[],
    absolute_errors: number[][],
    relative_errors: number[][],
}

type DataType = data6;





const available_menus: number[] = [
    6  // 6 - current term 
]
const header_menu_texts: string[] = [
    "Эллиптические уравнения"
]
const MENUS_COUNT = available_menus.length;

let cur_menu_n = 6;

let forms: HTMLFormElement[] = [];
let data_menus: HTMLElement[] = [];


for (let i = 0; i < MENUS_COUNT; ++i) {
    let n = available_menus[i];

    forms[i] = document.getElementById('t' + n + '_form') as HTMLFormElement;
    data_menus[i] = document.getElementById('t' + n + '_data_menu') as HTMLElement;
}

const forms_display = forms[0].style.display;
const data_display = data_menus[0].style.display;

let cur_menu_id: number;
for (let i = 0; i < MENUS_COUNT; ++i) {
    forms[i].style.display = 'none';
    data_menus[i].style.display = 'none';

    if (cur_menu_n == available_menus[i]) {
        cur_menu_id = i;
        break;
    }
}

type response_function = ((request: RequestType, data: DataType) => void);
const response_fuctions: readonly response_function[] = [
    task6_response,
]

let cur_response_function: response_function;

function enable_current_menu() {
    header.innerText = header_menu_texts[cur_menu_id];
    cur_response_function = response_fuctions[cur_menu_id];
    forms[cur_menu_id].style.display = forms_display;
    data_menus[cur_menu_id].style.display = data_display;
}

enable_current_menu();

let output_text = "";

function choose_menu(value: number) {
    forms[cur_menu_id].style.display = 'none';
    data_menus[cur_menu_id].style.display = 'none';

    cur_menu_n = Number(value);
    for (let i = 0; i < available_menus.length; ++i) {
        if (available_menus[i] == cur_menu_n) {
            cur_menu_id = i;
        }
    }

    enable_current_menu();
}

function add_to_output(data: string) {
    output_text += data;
}

let previous_time: number;
let new_time: number;
function start_timer() {
    let cur_date = new Date();
    previous_time = cur_date.getTime();
    new_time = previous_time;
}

function stop_timer() {
    let cur_date = new Date();
    new_time = cur_date.getTime();
}

function get_time_difference(): number {
    // alert("Prev t " + previous_time);
    // alert("New t " + new_time);
    // alert("Diff " + (new_time - previous_time));
    return (new_time - previous_time) / 1000;
}

function print_output() {
    var current_time = new Date();

    data_output.innerText = `Date: `
        + `${current_time.getDate().toString().padStart(2, '0')}/`
        + `${(current_time.getMonth() + 1).toString().padStart(2, '0')}/`
        + `${current_time.getFullYear().toString().padStart(2, '0')}`
        + `\nTime: `
        + `${current_time.getHours().toString().padStart(2, '0')}:`
        + `${current_time.getMinutes().toString().padStart(2, '0')}:`
        + `${current_time.getSeconds().toString().padStart(2, '0')}`
        + `\nInner ETA: ${get_time_difference().toString()}`
        + '\n-------------------\n'
        + output_text
        + '\n-------------------\n' +
        data_output.innerText;

    output_text = "";
}

function get_data_text(data: DataType): string {
    let text: string = '';

    switch (cur_menu_n) {
        case 6:
            let data6 = data as data6;

            text = `${data6.solver}\n\n\n`
                // + `Real solution:\n${data6.real_solution.toString()}\n\n`
                // + `Absolute errors:\n${data6.absolute_errors.toString()}\n\n`
                // + `Relative errors:\n${data6.relative_errors.toString()}\n\n`
                // + `Tau:\n\t${data6.backtrace.tau}\n\n`
                // + `Taus:\n\t${data6.backtrace.taus.toString()}\n\n`
                + `Chebyshev parameters:\n\t${data6.backtrace_data.cheb_params.toString()}\n\n\n`
                + `ETA:\n\t${data6.delta_time}`;
            break;
        default:
    }

    return text;
}

function get_request_text(request: RequestType) {
    let text: string = '';

    switch (cur_menu_n) {
        case 6:
            let request6 = request as request6;

            text = `Solving method: ${request6.solving_method.toString()}\n`
                + `t: ${request6.t.toString()}\n`
                + `eps: ${request6.eps.toString()}\n`
                + `calculation: ${request6.calculation.toString()}\n`
                + `max_steps: ${request6.max_steps.toString()}\n`
                + `n: ${request6.n.toString()}\n`
                + `m: ${request6.m.toString()}\n`
                + `optimization: ${request6.optimization.toString()}\n`
                + `threads_count: ${request6.threads_count.toString()}\n`;
            + `backtrace_level: ${request6.backtrace_level.toString()}\n`
                + `first_end_preserving: ${request6.first_end_preserving.toString()}\n`
                + `middle_division: ${request6.middle_division.toString()}\n`
            break;
        default:
    }

    return text;
}

async function send_request_and_get_asnwer(request: RequestType) {
    let request_text = get_request_text(request);
    add_to_output(request_text);
    data_text1.innerText = "Processing request:\n\t" + request_text;
    let response_text: string = '';

    let request_method: string = 'task' + cur_menu_n.toString() + '_request';

    start_timer();

    await invoke(request_method, { request: request })
        .then((response) => {
            stop_timer();
            let data = response as DataType;

            data_text1.innerText = get_data_text(data);
            data_text2.innerText = `${data.delta_time.toString().padEnd(6, "0")} ${request.optimization.toString()} ${request.solving_method.toString()}\n${data_text2.innerText}`;

            cur_response_function(request, data);
        })
        .catch((error) => {
            stop_timer();
            data_text1.innerText = error + "\n\n\nTried to process request:\n\t" + JSON.stringify(request);
            add_to_output(response_text);
        });
    print_output();
}


function change_task() {
    let value = Number(main_form.cur_task.value);
    choose_menu(value);
}

function send_request_button_click() {
    fill_request_data();

    send_request_and_get_asnwer(request);
}

button_change_task.addEventListener('click', change_task);
button_send_request.addEventListener('click', send_request_button_click);

const t6_table_body = document.getElementById("t6_table_body") as HTMLTableElement;
const t6_grid1_body = document.getElementById("t6_grid1_body") as HTMLTableElement;
const t6_grid2_body = document.getElementById("t6_grid2_body") as HTMLTableElement;
const t6_grid1_header = document.getElementById("t6_grid1_header") as HTMLElement;
const t6_grid2_header = document.getElementById("t6_grid2_header") as HTMLElement;
const t6_surface_plot1 = document.getElementById("t6_surface_plot1") as HTMLElement;
const t6_surface_plot2 = document.getElementById("t6_surface_plot2") as HTMLElement;
const t6_surface_plot3 = document.getElementById("t6_surface_plot3") as HTMLElement;
const t6_surface_plot4 = document.getElementById("t6_surface_plot4") as HTMLElement;


function fill_request_data() {
    let f = forms[cur_menu_id];

    switch (cur_menu_n) {
        case 6:
            request = {
                solving_method: f.solving_method.value,
                t: parseFloat(f.t.value),
                eps: parseFloat(f.eps.value),
                calculation: f.calculation.value,
                max_steps: Number(f.max_steps.value),
                n: Number(f.n.value),
                m: Number(f.m.value),
                optimization: f.optimization.value,
                threads_count: Number(f.threads_count.value),
                backtrace_level: f.backtrace_level.value,
                first_end_preserving: Number(f.first_end_preserving.value),
                middle_division: Number(f.middle_division.value),
            }
            break;
        default:
    }
}

function td_cell(data: string | number) {
    return '<td>' + data + '</td>';
}

function clear_children_except(element: HTMLElement, except: number) {
    let childrens = element.children;
    if (childrens !== null) {
        for (let i = childrens.length - 1; i > except - 1; --i) {
            element.removeChild(childrens[i]);
        }
    }
}


function task6_response(request: RequestType, data: DataType) {
    request = request as request6;
    data = data as data6;


    clear_children_except(t6_table_body, 1);
    clear_children_except(t6_grid1_header, 1);
    clear_children_except(t6_grid2_header, 1);
    clear_children_except(t6_grid1_body, 1);
    clear_children_except(t6_grid2_body, 1);
    clear_children_except(t6_surface_plot1, 0);
    clear_children_except(t6_surface_plot2, 0);
    clear_children_except(t6_surface_plot3, 0);
    clear_children_except(t6_surface_plot4, 0);

    let u = data.solution;

    let cur_text;
    if (request.backtrace_level == BacktraceLevel.Full || request.backtrace_level == BacktraceLevel.Tables) {
        const steps_absolute_errors = data.steps_absolute_errors;
        const steps_relative_errors = data.steps_relative_errors;
        const taus = data.backtrace_data.taus;

        for (let i = 0; i < data.backtrace_data.solver_steps.length; ++i) {
            let solver_step = data.backtrace_data.solver_steps[i];
            cur_text = '<tr class="content">';
            cur_text += '<td class="sticky_cell">' + solver_step.step + '</td>';
            cur_text += td_cell(taus[i]);
            cur_text += td_cell(solver_step.discrepancy_norm);
            cur_text += td_cell(solver_step.relative_discrepancy);
            cur_text += td_cell(steps_absolute_errors[i]);
            cur_text += td_cell(steps_relative_errors[i]);
            cur_text += td_cell(solver_step.norm_adj_difference);
            cur_text += td_cell(solver_step.apposterior_est);
            cur_text += td_cell(solver_step.spectral_radius_approach);
            cur_text += '</tr>';
            t6_table_body.insertAdjacentHTML("beforeend", cur_text);
        }


        let a = 0;
        let hx = Math.PI / (request.n - 1);

        let b = 0;
        let hy = 1 / (request.m - 1);

        cur_text = '';
        for (let j = 0; j < request.m; ++j) {
            cur_text += '<td class="sticky_cell">' + (b + hy * j).toString() + '</td>';
        }
        t6_grid1_header.insertAdjacentHTML('beforeend', cur_text);
        t6_grid2_header.insertAdjacentHTML('beforeend', cur_text);


        // add_to_output("Solution:\n\t" + sol);
        // add_to_output("Calculated solution\n\t" + u);
        // add_to_output("Tau\n\t" + data.backtrace.tau);
        // add_to_output("Cheb params\n\t" + data.backtrace.cheb_params);


        for (let i = 0; i < request.n; ++i) {

            cur_text = '<tr class="content">';
            cur_text += '<td class="sticky_cell">' + (a + hx * i).toString() + '</td>';
            for (let j = 0; j < request.m; ++j) {
                cur_text += td_cell(u[i][j].toString());
            }
            cur_text += '</tr>';
            t6_grid1_body.insertAdjacentHTML('beforeend', cur_text);

            cur_text = '<tr class="content">';
            cur_text += '<td class="sticky_cell">' + (a + hx * i).toString() + '</td>';
            for (let j = 0; j < request.m; ++j) {
                cur_text += td_cell(data.real_solution[i][j].toString());
            }
            cur_text += '</tr>';
            t6_grid2_body.insertAdjacentHTML('beforeend', cur_text);
        }

    }

    if (request.backtrace_level == BacktraceLevel.Full || request.backtrace_level == BacktraceLevel.Graphs) {
        let layout = {
            scene: { camera: { eye: { x: 1.97, y: -2.2, z: 0.3 } } },
            autosize: false,
            width: 700,
            height: 700,
            margin: {
                l: 15,
                r: 30,
                b: 45,
                t: 80,
            },
            paper_bgcolor: "#1b1b1b",
        };

        let graph_real_solution = [{
            z: data.real_solution,
            type: 'surface',
            contours: {
                z: {
                    show: true,
                    usecolormap: true,
                    highlightcolor: "#42f462",
                    project: { z: true }
                }
            }
        }];

        let graph_calculated_solution = [{
            z: u,
            type: 'surface',
            contours: {
                z: {
                    show: true,
                    usecolormap: true,
                    highlightcolor: "#42f462",
                    project: { z: true }
                }
            }
        }];

        let graph_absolute_errors = [{
            z: data.absolute_errors,
            type: 'surface',
            contours: {
                z: {
                    show: true,
                    usecolormap: true,
                    highlightcolor: "#42f462",
                    project: { z: true }
                }
            }
        }];

        let graph_relative_errors = [{
            z: data.relative_errors,
            type: 'surface',
            contours: {
                z: {
                    show: true,
                    usecolormap: true,
                    highlightcolor: "#42f462",
                    project: { z: true }
                }
            }
        }];

        Plotly.newPlot(t6_surface_plot1, graph_real_solution as Plotly.Data[], layout);
        Plotly.newPlot(t6_surface_plot2, graph_calculated_solution as Plotly.Data[], layout);
        Plotly.newPlot(t6_surface_plot3, graph_absolute_errors as Plotly.Data[], layout);
        Plotly.newPlot(t6_surface_plot4, graph_relative_errors as Plotly.Data[], layout);
    }
}