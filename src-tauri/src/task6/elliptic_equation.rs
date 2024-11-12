pub trait EllipticEquation {
    fn f(&self, x: f64, y: f64) -> f64;

    fn p(&self, x: f64, y: f64) -> f64;
    fn q(&self, x: f64, y: f64) -> f64;
}
