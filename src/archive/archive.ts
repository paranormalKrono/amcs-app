type request1 = {
    request_type: string,
    solving_method: string,
    a: number,
    b: number,
    divisions_count: number,
    n: number,
    t: number,
    is_reserving: boolean,
    is_reserving_lis: boolean,
}
type request4 = {
    a: number,
    b: number,
    n: number,
    m: number,
    mu: number,
}


type RequestType = request1 | request4;


type data1tablerow = {
    condition_number: number,
    solution: number[],
    discrepancy: number[],
}

type data1 = {
    table_rows: data1tablerow[],
    sol_t: number,
    state: string,
}

type data4sol = {
    solving_method: string,
    u: number[][],
}

type data4 = {
    data_file_name: string,
    sols: data4sol[],
    kappas: number[][],
}


type DataType = data1 | data4;


const available_menus: number[] = [
    1, // 1 - completed
    4, // 4 - completed
    6  // 6 - current term 
]

// const response_fuctions: readonly response_function[] = [
//     task1_response,
//     task4_response,
//     task6_response,
// ]


let text_for_output = "";

const header_menu_texts: string[] = [
    "Ритц",
    "Разностные схемы",
    "Эллиптические уравнения"
]

const t1_table_body = document.getElementById("t1_table_body") as HTMLTableElement;
const t1_table_tdyn = document.getElementById('td_yn') as HTMLTableCellElement;
const t1_table_tdyxyn = document.getElementById('td_yxyn') as HTMLTableCellElement;

const t4_surface_plot_name = "t4_surface_plot";
const t4_table_name = "t4_table";
const t4_surface_plots_count = 4;
const t4_surface_plots: HTMLElement[] = [];
for (let i = 0; i < t4_surface_plots_count; ++i) {
    t4_surface_plots[i] = document.getElementById(t4_surface_plot_name + (i + 1)) as HTMLTableElement;
}

//function fill_request_data() {
//    let f = forms[cur_menu_id];
//
//    switch (cur_menu_n) {
//        case 1:
//            request = {
//                request_type: f.request_type.value,
//                solving_method: f.solving_method.value,
//                a: parseFloat(f.a.value),
//                b: parseFloat(f.b.value),
//                divisions_count: Number(f.divisions_count.value),
//                n: Number(f.n.value),
//                t: Number(f.t.value),
//                is_reserving: f.is_reserving.checked,
//                is_reserving_lis: f.is_reserving_lis.checked,
//            }
//            break;
//        case 4:
//            request = {
//                a: parseFloat(f.a.value),
//                b: parseFloat(f.b.value),
//                n: Number(f.n.value),
//                m: Number(f.m.value),
//                mu: parseFloat(f.mu.value),
//            }
//            break;
//        case 6:
//            request = {
//                solving_method: f.solving_method.value,
//                optimization: f.optimization.value,
//                t: parseFloat(f.t.value),
//                eps: parseFloat(f.eps.value),
//                calculation: f.calculation.value,
//                max_steps: Number(f.max_steps.value),
//                n: Number(f.n.value),
//                m: Number(f.m.value),
//                is_backtrace: f.is_backtrace.checked,
//                first_end_preserving: Number(f.first_end_preserving.value),
//                middle_division: Number(f.middle_division.value),
//            }
//            break;
//        default:
//    }
//}


// function task1_response(request: RequestType, data: DataType) {
//     request = request as request1;
//     data = data as data1;
//
//     if (request.request_type == "CalculateTable") {
//
//         clear_children_except(t1_table_body, 1);
//
//         t1_table_tdyn.colSpan = request.divisions_count;
//         t1_table_tdyxyn.colSpan = request.divisions_count;
//
//         let a = request.a;
//         let h = (request.b - a) / (request.divisions_count - 1);
//
//         let cur_text = '<tr>';
//         for (let i = 0; i < request.divisions_count; ++i) {
//             cur_text += '<td>x=' + a + '</td>';
//             a += h;
//         }
//         a = request.a;
//         for (let i = 0; i < request.divisions_count; ++i) {
//             cur_text += '<td>x=' + a + '</td>';
//             a += h;
//         }
//         cur_text += '</tr>';
//
//         t1_table_body.insertAdjacentHTML("beforeend", cur_text);
//
//         for (let i = 0; i < data.table_rows.length; ++i) {
//             let table_row = data.table_rows[i];
//             cur_text = '<tr class="content">';
//             cur_text += '<td>' + (i + 1) + '</td>';
//             cur_text += '<td>' + table_row.condition_number + '</td>';
//             for (let j = 0; j < request.divisions_count; ++j) {
//                 cur_text += '<td>' + table_row.solution[j] + '</td>';
//             }
//             for (let j = 0; j < request.divisions_count; ++j) {
//                 cur_text += '<td>' + table_row.discrepancy[j] + '</td>';
//             }
//             cur_text += '</tr>';
//             t1_table_body.insertAdjacentHTML("beforeend", cur_text);
//         }
//     }
//     if (request.request_type == "CalculateSolution") {
//         text_for_output += "\n\n Sol t = " + data.sol_t;
//     }
//     data_text.innerText = data.state;
//     text_for_output += '\n\n' + data.state;
//
//     add_to_output(JSON.stringify(request));
//     add_to_output(text_for_output);
// }
//
// function task4_response(request: RequestType, data: DataType) {
//     request = request as request4;
//     data = data as data4;
//
//     let l = t4_surface_plots_count;
//
//     // add_to_output(JSON.stringify(data));
//
//     let layout = {
//         title: "",
//         scene: { camera: { eye: { x: -1.2, y: -1.0, z: 0.2 } } },
//         autosize: false,
//         width: 600,
//         height: 400,
//         margin: {
//             l: 15,
//             r: 30,
//             b: 20,
//             t: 50,
//         },
//         paper_bgcolor: "#1b1b1b",
//     };
//
//     for (let i = 0; i < l; ++i) {
//         let sol = data.sols[i];
//
//         let data_u = [{
//             z: sol.u,
//             type: 'surface',
//             contours: {
//                 z: {
//                     show: true,
//                     usecolormap: true,
//                     highlightcolor: "#42f462",
//                     project: { z: true }
//                 }
//             }
//         }];
//
//         layout.title = sol.solving_method;
//
//         Plotly.newPlot(t4_surface_plot_name + cur_menu_n, data_u as Plotly.Data[], layout);
//     }
//
//     let table_header = [['Kappas in x|t']]
//     let table_side = []
//
//     let cur_x = 1;
//     let cur_t = 0;
//     let h = request.a / (request.n - 1);
//     let tau = request.b / (request.m - 1);
//
//     for (let i = 0; i < request.n; i++) {
//         table_side.push([String(cur_x)]);
//         cur_x += h;
//     }
//     for (let i = 0; i < request.m; i++) {
//         table_header.push([String(cur_t)]);
//         cur_t += tau;
//     }
//
//     let table_layout = {
//         autosize: true,
//         width: 1200,
//         height: 800,
//         margin: {
//             l: 10,
//             r: 10,
//             b: 10,
//             t: 20,
//         },
//         paper_bgcolor: "#1b1b1b"
//     };
//
//     // data.kappas.unshift(table_side);
//     let data_kappa = [{
//         type: 'table',
//         columnwidth: [40],
//         header: {
//             values: table_header,
//             align: ["left"],
//             height: 30,
//             line: { color: '#686868', width: 1 },
//             fill: { color: '#1b1b1b' },
//             font: { family: "CascadiaCode_VTT", size: 16, color: "white" },
//             format: "3f"
//         },
//         cells: {
//             values: [table_side, data.kappas],
//             align: ["left"],
//             height: 30,
//             line: { color: "#686868", width: 1 },
//             fill: { color: ['#1b1b1b', '#1b1b1b'] },
//             font: { family: "CascadiaCode_VTT", size: 16, color: "white" },
//             format: "3f"
//         }
//     }]
//
//     Plotly.newPlot(t4_table_name, data_kappa as Plotly.Data[], table_layout);
//
//     add_to_output(JSON.stringify(request));
// }