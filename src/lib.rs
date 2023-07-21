mod consts;

pub trait LnExp {
    /// Returns `ln(1 - exp(x))`, computed as described in
    /// Martin Maechler (2012), Accurately Computing log(1 − exp(− |a|))
    /// http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
    fn ln_1m_exp(&self) -> Self;

    /// Returns `ln(1 + exp(x))`, computed as described in Maechler (2012).
    fn ln_1p_exp(&self) -> Self;

    /// Returns `ln(exp(x) - 1)`, the inverse of `ln_1p_exp`, computed by inverting
    /// the case analysis described in Maechler (2012).
    fn ln_exp_m1(&self) -> Self;

    /// Returns the logit, mapping from the closed interval [0,1] to a real number.
    fn logit(&self) -> Self;

    /// Returns the inverse-logit mapping from a real number to the closed interval [0,1].
    fn inv_logit(&self) -> Self;

    /// Returns the natural logarithm of the inv_logit function, computed
    /// more carefully than the composition of functions `x.inv_logit().ln()`.
    fn ln_inv_logit(&self) -> Self;
}

macro_rules! impl_lnexp {
    // N.B. the comments below are written for the `f64` implementation;
    // identical logic is applied to the `f32` case.
    ( $f:ident ) => {
        impl LnExp for $f {
            // See: Martin Maechler (2012), Accurately Computing log(1 − exp(− |a|))
            // http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
            fn ln_1m_exp(&self) -> $f {
                if *self < -std::$f::consts::LN_2 {
                    (-self.exp()).ln_1p()
                } else {
                    (-self.exp_m1()).ln()
                }
            }

            // See Section 3 of Maechler (2012).
            fn ln_1p_exp(&self) -> $f {
                if *self <= crate::consts::$f::X0 {
                    self.exp()
                } else if *self <= crate::consts::$f::X1 {
                    self.exp().ln_1p()
                } else if *self <= crate::consts::$f::X2 {
                    *self + (-*self).exp()
                } else {
                    *self
                }
            }

            fn ln_exp_m1(&self) -> $f {
                // can skip inversion of *self < -37.0 case as log1p(exp(x))
                // yields the same result as log(x).
                if *self <= crate::consts::$f::X1 {
                    self.exp_m1().ln()
                } else if *self <= crate::consts::$f::X2 {
                    *self - (-*self).exp()
                } else {
                    *self
                }
            }

            fn logit(&self) -> $f {
                (*self / (1.0 - *self)).ln()
            }

            fn inv_logit(&self) -> $f {
                // -744.4400719213812 is the logit of 5.0e-324, the smallest
                // value which can be represented as a 64-bit floating point
                // number.  Thus, anything less than logit(5.0e-324) must
                // necessarily be zero in order to preserve the identity
                // logit(inv_logit(x)) = x, and, alternately,
                // inv_logit(logit(5.0e-324)) = 5.0e-324.
                if *self < crate::consts::$f::ILOGIT_LOWER {
                    0.0
                } else {
                    // The next floating point number, 36.73680056967711,
                    // results in a value which when exponentiated is greater
                    // than ldexp(1.0, 53). When 1.0 is added to the
                    // exponentiated value, rounding causes the denominator to
                    // increase, leading to a value which is 2^-53 smaller
                    // than that computed from 36.7368005696771. Several of
                    // theThe function is no longer monotonic, thus, this is a
                    // suitable cutpoint.  In other words, logit(1.0 -
                    // f64::EPSILON / 2.0) = 36.7368005696771, thus, as the
                    // next linear-scale input (1.0) to `logit` maps to +inf,
                    // 1.0 must be returned for any value greater than
                    // logit(1.0 - f64::EPSILON / 2.0) when transforming from
                    // logit-scale to linear.
                    if *self > crate::consts::$f::ILOGIT_UPPER {
                        1.0
                    } else {
                        let e = self.exp();
                        e / (1.0 + e)
                    }
                }
            }

            fn ln_inv_logit(&self) -> $f {
                -(-*self).ln_1p_exp()
            }
        }
    };
}

impl_lnexp!(f64);
impl_lnexp!(f32);

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(test)]
    mod f64_impl {
        use super::*;

        #[test]
        fn ln_1m_exp_works() {
            let x: f64 = -1e-20;
            assert!(x.ln_1m_exp().is_finite());

            let x: f64 = -1.0;
            assert_eq!(x.ln_1m_exp(), (-x.exp()).ln_1p());

            let x: f64 = -0.1;
            assert_eq!(x.ln_1m_exp(), (-x.exp_m1()).ln());
        }

        #[test]
        fn ln_1p_exp_identity() {
            let xs: Vec<f64> = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                0.0,
                1.0,
                -50.0,
                -37.0,
                18.0,
                33.3,
                50.0,
            ];
            for x in xs {
                assert_eq!(x.ln_1p_exp().ln_exp_m1(), x);
            }
        }

        #[test]
        fn ln_1p_exp_misc() {
            let ln2 = std::f64::consts::LN_2;
            let eps = f64::EPSILON;
            let x: f64 = 0.0;
            assert!((x.ln_1p_exp() - ln2).abs() < eps);

            let e = std::f64::consts::E;
            let x: f64 = 1.0;
            assert!((x.ln_1p_exp() - e.ln_1p()).abs() < eps);

            let x: f64 = -1.0;
            assert!((x.ln_1p_exp() - (e.ln_1p() - 1.0)).abs() < eps);

            let xs: Vec<f64> = vec![1e2, 1e3, 1e4, 1e5];
            for x in xs {
                assert!((x.ln_1p_exp() - x).abs() < eps);
            }

            let xs: Vec<f64> = vec![-1e3, -1e4, -1e5];
            for x in xs {
                assert_eq!(x.ln_1p_exp(), 0.0);
            }
        }

        #[test]
        fn logit_identity() {
            let xs: Vec<f64> = vec![5e-324, 0.0, 0.5, 1.0];
            for x in xs {
                assert_eq!(x.logit().inv_logit(), x);
            }

            // Approximate identity
            let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![
                eps,
                eps / 2.0,
                1.0 - eps,
                1.0 - eps / 2.0,
                0.5 - eps,
                0.5 + eps,
            ];
            for x in xs {
                assert!((x.logit().inv_logit() - x).abs() < eps);
            }
        }

        #[test]
        fn inv_logit_misc() {
            let x: f64 = f64::INFINITY;
            assert_eq!(x.inv_logit(), 1.0);

            let x: f64 = f64::NEG_INFINITY;
            assert_eq!(x.inv_logit(), 0.0);

            let x: f64 = -744.0;
            assert!(x.inv_logit() > 0.0);

            let x: f64 = 36.0;
            assert!(x.inv_logit() < 1.0);

            let x: f64 = 0.5;
            assert_eq!(x.logit(), 0.0);

            let eps = f64::EPSILON;
            let x: f64 = 2.0;
            assert!((x.inv_logit() - 1.0 / (1.0 + (-2.0_f64).exp())).abs() < eps);
            assert!((x.inv_logit().logit() - 2.0) < eps);
        }

        #[test]
        fn logit_smooth_upper_bound() {
            let x: f64 = 36.7368005696771;
            assert_eq!(x.inv_logit(), 1.0 - f64::EPSILON / 2.0);

            let x: f64 = 36.73680056967711;
            assert_eq!(x.inv_logit(), 1.0);
        }

        #[test]
        fn logit_smooth_lower_bound() {
            let x: f64 = 5e-324;
            let rhs: f64 = -744.4400719213812;
            assert_eq!(x.logit(), rhs);

            let x: f64 = -744.4400719213813;
            assert_eq!(x.inv_logit(), 0.0);
        }

        #[test]
        fn ln_inv_logit_works() {
            // let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                // f64::EPSILON,
                0.0,
                // 1.0,
                -50.0,
                -37.0,
                // 18.0,
                // 33.3,
                // 50.0,
            ];

            for x in xs {
                let lhs = x.ln_inv_logit();
                let rhs = x.inv_logit().ln();
                assert_eq!(lhs, rhs);
                // assert!((lhs - rhs).abs() < eps);
            }
        }
    }
}
