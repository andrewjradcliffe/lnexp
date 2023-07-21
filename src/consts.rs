pub mod f64 {
    // Thresholds for `ln_1p_exp` based on
    // http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    // with improvements from
    // https://github.com/JuliaStats/LogExpFunctions.jl/files/8218470/log1pexp.pdf
    use std::f64::consts::LN_2;
    const PRECISION: f64 = f64::MANTISSA_DIGITS as f64;
    const LN_OF_P_MUL_LN_2: f64 = 3.6037789929704576; // p_ln2.ln()
    pub const X0: f64 = -PRECISION * LN_2;
    pub const X1: f64 = (PRECISION - 1.0) * LN_2 / 2.0;
    pub const X2: f64 = {
        let p_ln2 = PRECISION * LN_2;
        p_ln2 + LN_OF_P_MUL_LN_2 * (1.0 / p_ln2 - 1.0)
    };

    // inverse logit bounds
    pub const ILOGIT_LOWER: f64 = -744.4400719213812; // logit(5.0e-324)
    pub const ILOGIT_UPPER: f64 = 36.7368005696771; // logit(1.0 - f64::EPSILON / 2.0)
}

pub mod f32 {
    use std::f32::consts::LN_2;
    const PRECISION: f32 = f32::MANTISSA_DIGITS as f32;
    const LN_OF_P_MUL_LN_2: f32 = 2.8115408; // p_ln2.ln()
    pub const X0: f32 = -PRECISION * LN_2;
    pub const X1: f32 = (PRECISION - 1.0) * LN_2 / 2.0;
    pub const X2: f32 = {
        let p_ln2 = PRECISION * LN_2;
        p_ln2 + LN_OF_P_MUL_LN_2 * (1.0 / p_ln2 - 1.0)
    };

    pub const ILOGIT_LOWER: f32 = -103.27893; // logit(1.0e-45_f32)
    pub const ILOGIT_UPPER: f32 = 16.635532; // logit(1.0 - f32::EPSILON / 2.0)
}
