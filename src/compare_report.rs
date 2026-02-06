#[derive(Debug, Default, Clone)]
pub struct CompareReport {
    pub n: usize,
    pub n_over: u64,
    pub max_abs: f64,
    pub sum_abs: f64,
    pub sum_sq: f64,

    // for Pearson
    sum_x: f64,
    sum_y: f64,
    sum_x2: f64,
    sum_y2: f64,
    sum_xy: f64,

}    

impl CompareReport {
    pub fn update_from_bins(&mut self, a: &[f64], b: &[f64], eps: f64) {
        debug_assert_eq!(a.len(), b.len());

        for (&x, &y) in a.iter().zip(b.iter()) {
            let d = (x - y).abs();
            self.n += 1;
            self.sum_abs += d;
            self.sum_sq += d * d;

            if d > self.max_abs {
                self.max_abs = d;
            }
            if d > eps {
                self.n_over += 1;
            }
            // Pearson accumulators
            self.sum_x += x;
            self.sum_y += y;
            self.sum_x2 += x * x;
            self.sum_y2 += y * y;
            self.sum_xy += x * y;
        }
    }

    pub fn merge(&mut self, other: &CompareReport) {
        self.n += other.n;
        self.n_over += other.n_over;
        self.sum_abs += other.sum_abs;
        self.sum_sq += other.sum_sq;
        if other.max_abs > self.max_abs {
            self.max_abs = other.max_abs;
        }
        // Pearson accumulators
        self.sum_x += other.sum_x;
        self.sum_y += other.sum_y;
        self.sum_x2 += other.sum_x2;
        self.sum_y2 += other.sum_y2;
        self.sum_xy += other.sum_xy;
    }

    pub fn finish(&self) -> (usize, f64, f64, f64, f64, f64) {
        // (
        //   n_over_eps,
        //   frac_n_over_eps,
        //   mean_abs,
        //   var_abs,
        //   rmse,
        //   max_abs
        // )
        if self.n == 0 {
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        } else {
            let n = self.n as f64;

            let frac = self.n_over_eps as f64 / n;
            let mean_abs = self.sum_abs / n;

            // variance of absolute difference
            // E[d²] − (E[|d|])²
            let var_abs = (self.sum_sq / n) - (mean_abs * mean_abs);

            let rmse = (self.sum_sq / n).sqrt();

            (
                self.n_over_eps,
                frac,
                mean_abs,
                var_abs.max(0.0), // numerical safety
                rmse,
                self.max_abs,
            )
        }
    }

    pub fn finish_names() -> (&'static str, &'static str, &'static str, &'static str, &'static str) {
        (
            "n_over_eps",
            "frac_n_over_eps",
            "mean_abs",
            "var_abs",
            "rmse",
            "max_abs",
        )
    }

    pub fn pearson(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }

        let n = self.n as f64;
        let num = self.sum_xy - (self.sum_x * self.sum_y) / n;
        let den_x = self.sum_x2 - (self.sum_x * self.sum_x) / n;
        let den_y = self.sum_y2 - (self.sum_y * self.sum_y) / n;

        if den_x <= 0.0 || den_y <= 0.0 {
            0.0
        } else {
            num / (den_x * den_y).sqrt()
        }
    }
}