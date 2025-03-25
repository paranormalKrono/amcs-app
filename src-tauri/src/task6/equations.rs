use meval::Expr;

pub trait EllipticEquation {
    fn p(&self, x: f64, y: f64) -> f64;
    fn q(&self, x: f64, y: f64) -> f64;

    fn f(&self, x: f64, y: f64) -> f64;

    fn mu(&self, x: f64, y: f64) -> f64;
}

pub trait HyperbolicEquation {
    fn p(&self, x: f64, y: f64) -> f64;
    fn q(&self, x: f64, y: f64) -> f64;
    fn c(&self, x: f64, y: f64) -> f64;

    fn f(&self, x: f64, y: f64) -> f64;

    fn mu(&self, x: f64, y: f64) -> f64;
}

pub struct EE3;

impl EllipticEquation for EE3 {
    fn f(&self, x: f64, y: f64) -> f64 {
        2_f64 * x.sin() * y.cos()
    }

    fn p(&self, _x: f64, _y: f64) -> f64 {
        1_f64
    }
    fn q(&self, _x: f64, _y: f64) -> f64 {
        1_f64
    }

    fn mu(&self, x: f64, y: f64) -> f64 {
        x.sin() * y.cos()
    }
}

impl std::fmt::Display for EE3 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let buf: String = "EE3: p = 1, q = 1, f = 2sin(x)cos(y)".to_string();

        write!(f, "{}", buf)
    }
}

// pub fn ee3_sol(x: f64, y: f64) -> f64 {
//     x.sin() * y.cos()
// }

type EEQFunctionsTuple = (
    Box<dyn Fn(f64, f64) -> f64>,
    Box<dyn Fn(f64, f64) -> f64>,
    Box<dyn Fn(f64, f64) -> f64>,
    Box<dyn Fn(f64, f64) -> f64>,
    Box<dyn Fn(f64, f64) -> f64>,
    Box<dyn Fn(f64, f64) -> f64>,
);

pub struct MevalEquation {
    pub meval_p: Box<dyn Fn(f64, f64) -> f64>,
    pub meval_q: Box<dyn Fn(f64, f64) -> f64>,
    pub meval_c: Box<dyn Fn(f64, f64) -> f64>,
    pub meval_f: Box<dyn Fn(f64, f64) -> f64>,
    pub meval_mu: Box<dyn Fn(f64, f64) -> f64>,
    pub meval_sol: Box<dyn Fn(f64, f64) -> f64>,
}

pub struct MevalEquationCfg {
    pub p: String,
    pub q: String,
    pub c: String,
    pub f: String,
    pub mu: String,
    pub sol: String,
}

impl MevalEquation {
    fn process_config(cfg: MevalEquationCfg) -> EEQFunctionsTuple {
        let mut ar: Vec<String> = vec![cfg.p, cfg.q, cfg.c, cfg.f, cfg.mu, cfg.sol];

        let mut arb = ar.iter_mut().map(|x| {
            let parsed = x.parse::<Expr>().unwrap();
            let binded = parsed.bind2("x", "y").unwrap();
            let b: Box<dyn Fn(f64, f64) -> f64> = Box::new(binded);
            b
        });

        (
            arb.next().unwrap(),
            arb.next().unwrap(),
            arb.next().unwrap(),
            arb.next().unwrap(),
            arb.next().unwrap(),
            arb.next().unwrap(),
        )
    }

    pub fn new(cfg: MevalEquationCfg) -> Self {
        let (p, q, c, f, mu, sol) = MevalEquation::process_config(cfg);
        MevalEquation {
            meval_p: p,
            meval_q: q,
            meval_c: c,
            meval_f: f,
            meval_mu: mu,
            meval_sol: sol,
        }
    }

    pub fn sol(&self, x: f64, y: f64) -> f64 {
        (self.meval_sol)(x, y)
    }
}

impl EllipticEquation for MevalEquation {
    fn f(&self, x: f64, y: f64) -> f64 {
        (self.meval_f)(x, y)
    }

    fn p(&self, x: f64, y: f64) -> f64 {
        (self.meval_p)(x, y)
    }

    fn q(&self, x: f64, y: f64) -> f64 {
        (self.meval_q)(x, y)
    }

    fn mu(&self, x: f64, y: f64) -> f64 {
        (self.meval_mu)(x, y)
    }
}

impl HyperbolicEquation for MevalEquation {
    fn p(&self, x: f64, y: f64) -> f64 {
        (self.meval_p)(x, y)
    }

    fn q(&self, x: f64, y: f64) -> f64 {
        (self.meval_q)(x, y)
    }

    fn c(&self, x: f64, y: f64) -> f64 {
        (self.meval_c)(x, y)
    }

    fn f(&self, x: f64, y: f64) -> f64 {
        (self.meval_f)(x, y)
    }

    fn mu(&self, x: f64, y: f64) -> f64 {
        (self.meval_mu)(x, y)
    }
}
