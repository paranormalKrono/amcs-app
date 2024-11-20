// fn parallel4_run_method_secondary(rx: &Receiver<Vec<Vec<f64>>>, txr: &Sender<Vec<Vec<f64>>>) {
//     loop {
//         let (l, mut d, r, mut b) = {
//             let rec = rx.recv();
//             let mut val = match rec {
//                 Ok(value) => value,
//                 Err(_) => {
//                     println!("Terminating.");
//                     break;
//                 }
//             };
//             (val.remove(0), val.remove(0), val.remove(0), val.remove(0))
//         };
//
//         let n2 = b.len() / 2;
//         let mut c1;
//
//         for i in 1..n2 {
//             c1 = l[i] / d[i - 1];
//             d[i] -= c1 * r[i - 1];
//             b[i] -= c1 * b[i - 1];
//         }
//
//         let m = rx.recv().unwrap();
//         let mut d20 = m[0][0];
//         let mut b20 = m[0][1];
//
//         let c = l[n2] / d[n2 - 1];
//         d20 -= c * r[n2 - 1];
//         b20 = (b20 - c * b[n2 - 1]) / d20;
//         b[n2 - 1] = (b[n2 - 1] - r[n2 - 1] * b20) / d[n2 - 1];
//
//         txr.send(vec![vec![b20]]).unwrap();
//
//         for i in (1..n2).rev() {
//             b[i - 1] = (b[i - 1] - r[i - 1] * b[i]) / d[i - 1];
//         }
//
//         txr.send(vec![b[..n2].to_vec()]).unwrap();
//     }
// }
//
// fn parallel4_run_method_main(
//     tx: &Sender<Vec<Vec<f64>>>,
//     rxr: &Receiver<Vec<Vec<f64>>>,
//     l: &[f64],
//     d: &mut [f64],
//     r: &[f64],
//     b: &mut [f64],
// ) {
//     let n2 = b.len() / 2;
//     let mut c2;
//
//     tx.send(vec![l.to_owned(), d.to_owned(), r.to_owned(), b.to_owned()])
//         .unwrap();
//
//     for i in (n2..(b.len() - 1)).rev() {
//         c2 = r[i] / d[i + 1];
//         d[i] -= c2 * l[i + 1];
//         b[i] -= c2 * b[i + 1];
//     }
//
//     tx.send(vec![vec![d[n2], b[n2]]]).unwrap();
//     b[n2] = rxr.recv().unwrap()[0][0];
//     d[n2] = 1_f64;
//
//     for i in n2..(b.len() - 1) {
//         b[i + 1] = (b[i + 1] - l[i + 1] * b[i]) / d[i + 1];
//     }
//
//     b[..n2].clone_from_slice(&rxr.recv().unwrap()[0]);
// }

// fn parallel2_run_method(al: &[f64], ad: &mut [f64], ar: &[f64], b: &mut [f64]) {
//     let n2 = b.len() / 2;
//
//     let (tx, rx) = mpsc::channel();
//     let (txr, rxr) = mpsc::channel();
//
//     let (ar1, ar2) = ar.split_at(n2);
//     let (ad1, ad2) = ad.split_at_mut(n2);
//     let (al1, al2) = al.split_at(n2);
//     let (b1, b2) = b.split_at_mut(n2);
//
//     thread::scope(|s| {
//         s.spawn(move || {
//             let mut c2;
//             for i in (0..(b2.len() - 1)).rev() {
//                 c2 = ar2[i] / ad2[i + 1];
//                 ad2[i] -= c2 * al2[i + 1];
//                 b2[i] -= c2 * b2[i + 1];
//             }
//
//             tx.send((al2[0], ad2[0], b2[0])).unwrap();
//             b2[0] = rxr.recv().unwrap();
//             ad2[0] = 1_f64;
//
//             for i in 0..(b2.len() - 1) {
//                 b2[i + 1] = (b2[i + 1] - al2[i + 1] * b2[i]) / ad2[i + 1];
//             }
//         });
//
//         s.spawn(move || {
//             let mut c1;
//             for i in 1..n2 {
//                 c1 = al1[i] / ad1[i - 1];
//                 ad1[i] -= c1 * ar1[i - 1];
//                 b1[i] -= c1 * b1[i - 1];
//             }
//
//             let (al20, mut ad20, mut b20) = rx.recv().unwrap();
//
//             let c = al20 / ad1[n2 - 1];
//             ad20 -= c * ar1[n2 - 1];
//             b20 = (b20 - c * b1[n2 - 1]) / ad20;
//             b1[n2 - 1] = (b1[n2 - 1] - ar1[n2 - 1] * b20) / ad1[n2 - 1];
//
//             txr.send(b20).unwrap();
//
//             for i in (1..n2).rev() {
//                 b1[i - 1] = (b1[i - 1] - ar1[i - 1] * b1[i]) / ad1[i - 1];
//             }
//         });
//     });
// }
//
// fn parallel_run_method(l: &[f64], d: &mut [f64], r: &[f64], b: &mut [f64]) {
//     let n2 = b.len() / 2;
//
//     thread::scope(|s| {
//         let (ad1, ad2) = d.split_at_mut(n2);
//         let (b1, b2) = b.split_at_mut(n2);
//
//         s.spawn(move || {
//             let mut c1;
//             for i in 1..n2 {
//                 // println!("A1_{}", i);
//                 // thread::sleep(time::Duration::from_millis(1));
//
//                 c1 = l[i] / ad1[i - 1];
//                 ad1[i] -= c1 * r[i - 1];
//                 b1[i] -= c1 * b1[i - 1];
//             }
//         });
//
//         s.spawn(move || {
//             let mut c2;
//             for i in (0..(b2.len() - 1)).rev() {
//                 // println!("B1_{}", i + n2);
//                 // thread::sleep(time::Duration::from_millis(1));
//
//                 c2 = r[i + n2] / ad2[i + 1];
//                 ad2[i] -= c2 * l[i + 1 + n2];
//                 b2[i] -= c2 * b2[i + 1];
//             }
//         });
//     });
//
//     let c = l[n2] / d[n2 - 1];
//     d[n2] -= c * r[n2 - 1];
//     b[n2] = (b[n2] - c * b[n2 - 1]) / d[n2];
//     b[n2 - 1] = (b[n2 - 1] - r[n2 - 1] * b[n2]) / d[n2 - 1];
//
//     thread::scope(|s| {
//         let (ad1, ad2) = d.split_at_mut(n2);
//         let (b1, b2) = b.split_at_mut(n2);
//
//         s.spawn(move || {
//             for i in (1..n2).rev() {
//                 b1[i - 1] = (b1[i - 1] - r[i - 1] * b1[i]) / ad1[i - 1];
//             }
//         });
//         s.spawn(move || {
//             for i in 0..(b2.len() - 1) {
//                 b2[i + 1] = (b2[i + 1] - l[i + 1 + n2] * b2[i]) / ad2[i + 1];
//             }
//         });
//     })
// }
//
// fn with_comments_parallel_run_method(l: &[f64], d: &mut [f64], r: &[f64], b: &mut [f64]) {
//     let n2 = b.len() / 2;
//
//     thread::scope(|s| {
//         let (ad1, ad2) = d.split_at_mut(n2);
//         let (b1, b2) = b.split_at_mut(n2);
//
//         s.spawn(move || {
//             let mut c1;
//             for i in 1..n2 {
//                 // println!("A1_{}", i);
//                 // thread::sleep(time::Duration::from_millis(1));
//
//                 // B C 0  i - 1
//                 // A B C  i
//                 c1 = l[i] / ad1[i - 1];
//                 ad1[i] -= c1 * r[i - 1];
//                 b1[i] -= c1 * b1[i - 1];
//                 // B C 0  i - 1
//                 // 0 K C  i
//             }
//         });
//
//         s.spawn(move || {
//             let mut c2;
//             for i in (0..(b2.len() - 1)).rev() {
//                 // println!("B1_{}", i + n2);
//                 // thread::sleep(time::Duration::from_millis(1));
//
//                 // A B C  i
//                 // 0 A B  i + 1
//                 c2 = r[i + n2] / ad2[i + 1];
//                 ad2[i] -= c2 * l[i + 1 + n2];
//                 b2[i] -= c2 * b2[i + 1];
//                 // A K 0  i
//                 // 0 A B  i + 1
//             }
//         });
//     });
//
//     // K C  n2 - 1
//     // A K  n2
//     let c = l[n2] / d[n2 - 1];
//     d[n2] -= c * r[n2 - 1];
//     b[n2] = (b[n2] - c * b[n2 - 1]) / d[n2];
//     // K C  n2 - 1
//     // 0 1  n2
//     b[n2 - 1] = (b[n2 - 1] - r[n2 - 1] * b[n2]) / d[n2 - 1];
//     // 1 0  n2 - 1
//     // 0 1  n2
//
//     thread::scope(|s| {
//         let (ad1, ad2) = d.split_at_mut(n2);
//         let (b1, b2) = b.split_at_mut(n2);
//
//         s.spawn(move || {
//             for i in (1..n2).rev() {
//                 // K C  i - 1
//                 // 0 1  i
//                 b1[i - 1] = (b1[i - 1] - r[i - 1] * b1[i]) / ad1[i - 1];
//                 // 1 0  i - 1
//                 // 0 1  i
//             }
//         });
//         s.spawn(move || {
//             for i in 0..(b2.len() - 1) {
//                 // 1 0  i
//                 // A K  i + 1
//                 b2[i + 1] = (b2[i + 1] - l[i + 1 + n2] * b2[i]) / ad2[i + 1];
//                 // 1 0  i
//                 // 0 1  i + 1
//             }
//         });
//     })
// }

//
// // #[cfg(test)]
// mod tests {
//     use crate::task6::*;
//
//     #[test]
//     fn check_chebyshev_12() {
//         assert_eq!(
//             get_optimal_chebyshev_params(12),
//             vec![1, 31, 15, 17, 7, 25, 9, 23, 3, 29, 13, 19]
//         );
//     }
//
//     #[test]
//     fn check_chebyshev_13() {
//         assert_eq!(
//             get_optimal_chebyshev_params(13),
//             vec![1, 31, 15, 17, 7, 25, 9, 23, 3, 29, 13, 19, 5]
//         );
//     }
//
//     #[test]
//     fn check_chebyshev_16() {
//         assert_eq!(
//             get_optimal_chebyshev_params(16),
//             vec![1, 31, 15, 17, 7, 25, 9, 23, 3, 29, 13, 19, 5, 27, 11, 21]
//         );
//     }
//
//     #[test]
//     fn check_chebyshev_7() {
//         assert_eq!(get_optimal_chebyshev_params(7), vec![1, 15, 7, 9, 3, 13, 5]);
//     }
//
//     #[test]
//     fn check_chebyshev_4() {
//         assert_eq!(get_optimal_chebyshev_params(4), vec![1, 7, 3, 5]);
//     }
//
//     #[test]
//     fn check_chebyshev_2() {
//         assert_eq!(get_optimal_chebyshev_params(2), vec![1, 3]);
//     }
// }

// Super Chebyshev. It's not how it should've work.
//
// let mut tmp;
// let mut calc_diff;
// let (opt, degree) = {
//     let mut degree: usize = n.ilog2() as usize;
//
//     let mut n2 = 1_usize;
//     for _ in 0..degree {
//         n2 <<= 1;
//     }
//     if n2 != n {
//         degree += 1;
//         n2 <<= 1;
//     }
//
//     let mut p2 = 1_usize;
//     calc_diff = 1_usize;
//
//     // finding optimal size for calculations
//     tmp = n & 1;
//     calc_diff += tmp;
//     let mut opt = n + tmp;
//     while opt != n2 {
//         opt >>= 1;
//         n2 >>= 1;
//         tmp = opt & 1;
//         opt += tmp;
//         calc_diff += p2 * tmp;
//         p2 *= 2;
//     }
//
//     (opt, degree)
// };
//
// let mut res = vec![0_usize; n + calc_diff];
// res[0] = 1;
// res[1] = 3;
//
// let mut m = 1_usize;
// let mut om = m;
// let mut p;
//
// for _ in 2..=degree {
//     m *= 2;
//     om *= 2;
//     if m == opt {
//         om -= 1;
//     }
//     p = om * 2;
//
//     for j in (0..om).rev() {
//         p -= 2;
//         res[p] = res[j];
//         res[p + 1] = 4 * m - res[p];
//     }
// }
//
// res.truncate(n);
// res
//

// Hell no!
//
// let mut fi: f64;
// let mut fj: f64;
//
// let mut fij: f64;
//
// fi = 1_f64;
// fj = 1_f64;
// fij = self.EQ.f(fi * self.hx, fj * self.hy);
// max_mod = (fij
//     + self.EQ.p(0.5_f64 * self.hx, self.hy) * self.u[0][1] / self.hx.powi(2)
//     + self.EQ.q(self.hx, 0.5_f64 * self.hy) * self.u[1][0] / self.hy.powi(2))
// .abs();
//
// fi = (self.N - 2) as f64;
// fj = 1_f64;
// fij = self.EQ.f(fi * self.hx, fj * self.hy);
// max_mod = f64::max(
//     max_mod,
//     (fij + self.EQ.p((fi + 0.5_f64) * self.hx, self.hy) * self.u[self.N - 1][1]
//         / self.hx.powi(2)
//         + self.EQ.q(fi * self.hx, (fj - 0.5_f64) * self.hy) * self.u[self.N - 2][0]
//             / self.hy.powi(2))
//     .abs(),
// );
//
// fi = 1_f64;
// fj = (self.M - 2) as f64;
// fij = self.EQ.f(fi * self.hx, fj * self.hy);
// max_mod = f64::max(
//     max_mod,
//     (fij + self.EQ.p((fi - 0.5_f64) * self.hx, self.hy) * self.u[0][self.M - 2]
//         / self.hx.powi(2)
//         + self.EQ.q(fi * self.hx, (fj + 0.5_f64) * self.hy) * self.u[1][self.M - 1]
//             / self.hy.powi(2))
//     .abs(),
// );
//
// fi = (self.N - 2) as f64;
// fj = (self.M - 2) as f64;
// fij = self.EQ.f(fi * self.hx, fj * self.hy);
// max_mod = f64::max(
//     max_mod,
//     (fij + self.EQ.p((fi + 0.5_f64) * self.hx, self.hy) * self.u[self.N - 1][self.M - 2]
//         / self.hx.powi(2)
//         + self.EQ.q(fi * self.hx, (fj + 0.5_f64) * self.hy)
//             * self.u[self.N - 2][self.M - 1]
//             / self.hy.powi(2))
//     .abs(),
// );
//
// for i in 2..(self.N - 2) {
//     fi = i as f64;
//     fj = 1_f64;
//
//     fij = self.EQ.f(fi * self.hx, fj * self.hy);
//
//     max_mod = f64::max(
//         max_mod,
//         (fij + self.EQ.q(fi * self.hx, (fj - 0.5_f64) * self.hy) * self.u[i][0]
//             / self.hy.powi(2))
//         .abs(),
//     );
//
//     fj = (self.M - 2) as f64;
//
//     fij = self.EQ.f(fi * self.hx, fj * self.hy);
//
//     max_mod = f64::max(
//         max_mod,
//         (fij + self.EQ.q(fi * self.hx, (fj + 0.5_f64) * self.hy) * self.u[i][self.M - 1]
//             / self.hy.powi(2))
//         .abs(),
//     );
// }
//
// for j in 2..(self.M - 2) {
//     fj = j as f64;
//     fi = 1_f64;
//
//     fij = self.EQ.f(fi * self.hx, fj * self.hy);
//
//     max_mod = f64::max(
//         max_mod,
//         (fij + self.EQ.p(fi * self.hx, (fj - 0.5_f64) * self.hy) * self.u[0][j]
//             / self.hx.powi(2))
//         .abs(),
//     );
//
//     fi = (self.N - 2) as f64;
//
//     fij = self.EQ.f(fi * self.hx, fj * self.hy);
//
//     max_mod = f64::max(
//         max_mod,
//         (fij + self.EQ.p(fi * self.hx, (fj + 0.5_f64) * self.hy) * self.u[self.N - 1][j]
//             / self.hx.powi(2))
//         .abs(),
//     );
// }
//
// for i in 3..(self.N - 3) {
//     fi = i as f64;
//     for j in 3..(self.N - 3) {
//         fj = j as f64;
//         fij = self.EQ.f(fi * self.hx, fj * self.hy);
//         max_mod = f64::max(max_mod, fij.abs());
//     }
// }
//
