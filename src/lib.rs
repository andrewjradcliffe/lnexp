mod consts;

/// The `LnExp` trait extends the functionality afforded by a `libm`, providing
/// careful evaluation of compositions of the `ln`, `exp`, `ln_1p` and `exp_m1`
/// functions so as to avoid underflow/overflow. Furthermore, in some cases, the
/// methods provided by this trait result in fewer calls than the equivalent
/// function composition.
pub trait LnExp {
    /// Returns `ln(1 - exp(x))`, computed as described in
    /// [Martin Maechler (2012), Accurately Computing log(1 − exp(− |a|))](http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf).
    /// For x > 0, the result is undefined, hence, the return value is `nan`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = -1e-20;
    /// assert!(x.ln_1m_exp().is_finite());
    /// // compare with the naive computation
    /// assert!((-x.exp()).ln_1p().is_infinite());
    ///
    /// let x: f64 = -50.0;
    /// assert_ne!(x.ln_1m_exp(), 0.0);
    /// // the result if we blindly applied `ln(-exp_m1(x))` to the entire range
    /// assert_eq!((-x.exp_m1()).ln(), 0.0);
    ///
    /// // An interesting property: the inverse of `ln_1m_exp` is `ln_1m_exp`.
    /// let x: f64 = -2.0;
    /// assert_eq!(x.ln_1m_exp().ln_1m_exp(), x);
    /// ```
    fn ln_1m_exp(&self) -> Self;

    /// Returns `ln(1 + exp(x))`, computed as described in Maechler (2012),
    /// for which a theoretical basis has been provided by [cossio (2022), untitled](https://github.com/JuliaStats/LogExpFunctions.jl/files/8218470/log1pexp.pdf).
    /// The inverse of this function is `ln_exp_m1`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = (1023.0 * (2.0_f64).ln()) + 1.0;
    /// assert!(x.ln_1p_exp().is_finite());
    /// // compare with naive computation
    /// assert_eq!(x.exp().ln_1p(), f64::INFINITY);
    /// ```
    fn ln_1p_exp(&self) -> Self;

    /// Returns `ln(exp(x) - 1)`, the inverse of `ln_1p_exp`, computed by inverting
    /// the case analysis described in Maechler (2012).
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = (1023.0 * (2.0_f64).ln()) + 1.0;
    /// assert!(x.ln_exp_m1().is_finite());
    /// // compare with naive computation
    /// assert_eq!(x.exp_m1().ln(), f64::INFINITY);
    /// ```
    fn ln_exp_m1(&self) -> Self;

    /// Returns the [logit](https://en.wikipedia.org/wiki/Logit), mapping from the closed interval \[0,1\] to a real number.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = 1.0;
    /// assert_eq!(x.logit(), f64::INFINITY);
    ///
    /// let x: f64 = 0.5;
    /// assert_eq!(x.logit(), 0.0);
    ///
    /// let x: f64 = 0.0;
    /// assert_eq!(x.logit(), f64::NEG_INFINITY);
    /// ```
    fn logit(&self) -> Self;

    /// Returns the inverse-logit mapping from a real number to the closed interval \[0,1\].
    /// This is the inverse of `logit`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = f64::INFINITY;
    /// assert_eq!(x.inv_logit(), 1.0);
    ///
    /// let x: f64 = 0.0;
    /// assert_eq!(x.inv_logit(), 0.5);
    ///
    /// let x: f64 = f64::NEG_INFINITY;
    /// assert_eq!(x.inv_logit(), 0.0);
    ///
    /// // Smooth lower bound by enforcing that values less than
    /// // `logit(5.0e-324)` (`logit(1.0e-45)` if `f32`) map to 0.0
    /// // when the inverse-logit is applied.
    /// let x: f64 = -745.0; // Slightly less than lower bound
    /// assert_eq!(x.inv_logit(), 0.0);
    /// // compare with naive computation
    /// let e = x.exp();
    /// assert_ne!(e / (1.0 + e), 0.0);
    ///
    /// // Smooth upper bound by enforcing that values greater than
    /// // `logit(1.0 - epsilon / 2.0)` map to 1.0 when the inverse-logit is applied.
    /// let x: f64 = 36.8; // Slightly larger than upper bound
    /// assert_eq!(x.inv_logit(), 1.0);
    /// // compare with naive computation
    /// let e = x.exp();
    /// assert_ne!(e / (1.0 + e), 1.0);
    /// ```
    fn inv_logit(&self) -> Self;

    /// Returns the natural logarithm of the inverse-logit function, computed
    /// more carefully than the composition of functions `x.inv_logit().ln()`.
    /// This utilizes the identity `ln(inv_logit(x)) = -ln(1 + exp(-x))`
    /// to compute the result via `ln_1p_exp`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = -745.0;
    /// assert!(x.ln_inv_logit().is_finite());
    /// // compare with naive computation
    /// assert_eq!(x.inv_logit().ln(), f64::NEG_INFINITY);
    /// ```
    fn ln_inv_logit(&self) -> Self;

    /// Returns the logit of the exponential of `self`, i.e. `logit(exp(x))`,
    /// computed more carefully and with fewer function calls than the composition.
    /// This is the inverse of `ln_inv_logit`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = 50.0;
    /// assert_eq!(x.ln_inv_logit().logit_exp(), x);
    /// let x: f64 = 743.0;
    /// assert!(x.ln_inv_logit().logit_exp().is_finite());
    /// // compare with naive computation
    /// assert_eq!(x.ln_inv_logit().exp().logit(), f64::INFINITY);
    /// ```
    fn logit_exp(&self) -> Self;

    /// Returns the natural logarithm of the 1 minus the inverse logit function,
    /// computed more carefully and with fewer function calls than the
    /// composition of functions `(1.0 - x.inv_logit()).ln()`. This exploits negation
    /// in the log-odds domain to compute the result via `ln_1p_exp`,
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = 50.0;
    /// assert!(x.ln_1m_inv_logit().is_finite());
    /// // compare with naive computation
    /// assert_eq!((1.0 - x.inv_logit()).ln(), f64::NEG_INFINITY);
    /// ```
    fn ln_1m_inv_logit(&self) -> Self;

    /// Returns the logit of 1 minus the exponential of `self`, i.e.
    /// `logit(1 - exp(x))`, computed more carefully and with fewer function calls
    /// than the composition. This is the inverse of `ln_1m_inv_logit`.
    ///
    /// # Examples
    /// ```
    /// use lnexp::LnExp;
    ///
    /// let x: f64 = 50.0;
    /// assert_eq!(x.ln_1m_inv_logit().logit_1m_exp(), x);
    /// let x: f64 = 743.0;
    /// assert!(x.ln_1m_inv_logit().is_finite());
    /// // compare with naive computation
    /// assert_eq!((1.0 - x.ln_1m_inv_logit().exp()).logit(), f64::INFINITY);
    /// ```
    fn logit_1m_exp(&self) -> Self;
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
                    // than that computed from 36.7368005696771.
                    // The function is no longer monotonic, thus, this is a
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

            fn logit_exp(&self) -> $f {
                // This uses the following identity:
                // logit(exp(x)) = log(exp(x) / (1 + exp(x))) = log(exp(x)) - log(1 - exp(x))
                *self - self.ln_1m_exp()
            }

            fn ln_1m_inv_logit(&self) -> $f {
                -self.ln_1p_exp()
            }

            fn logit_1m_exp(&self) -> $f {
                // This uses the same identity as `logit_exp`, followed by negation
                // on the log-odds scale.
                // That is, -logit(exp(x)) = log(1 - exp(x)) - log(exp(x))
                self.ln_1m_exp() - *self
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

            // limits
            let x: f64 = f64::NEG_INFINITY;
            assert_eq!(x.ln_1m_exp(), -0.0); // -0.0 due to log1p(-0.0)

            let x: f64 = 0.0;
            assert_eq!(x.ln_1m_exp(), f64::NEG_INFINITY);

            let x: f64 = -0.0;
            assert_eq!(x.ln_1m_exp(), f64::NEG_INFINITY);

            let x: f64 = 5.0e-324;
            assert!(x.ln_1m_exp().is_nan());
            let x: f64 = f64::MIN_POSITIVE;
            assert!(x.ln_1m_exp().is_nan());

            let x: f64 = -5.0e-324;
            assert!(x.ln_1m_exp().is_finite());
            let x: f64 = -f64::MIN_POSITIVE;
            assert!(x.ln_1m_exp().is_finite());

            // its own inverse
            let x: f64 = -0.5;
            assert_eq!(x.ln_1m_exp().ln_1m_exp(), x);

            let x: f64 = -f64::EPSILON;
            assert!((x.ln_1m_exp().ln_1m_exp() - x).abs() < 6.0e-31);
        }

        #[test]
        fn ln_1p_exp_identity() {
            let xs: Vec<f64> = vec![
                f64::INFINITY,
                f64::NEG_INFINITY,
                0.0,
                1.0,
                3.0,
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

            let x: f64 = 2.0;
            assert!((x.ln_1p_exp().ln_exp_m1() - x).abs() <= 2.0 * eps);

            // When naive method overflows
            let x: f64 = (1023.0 * (2.0_f64).ln()) + 1.0;
            assert!(x.ln_1p_exp().is_finite());

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
        fn ln_exp_m1_identity() {
            let xs: Vec<f64> = vec![f64::INFINITY, 0.0, 1.0, 18.0, 33.3, 50.0];
            for x in xs {
                assert_eq!(x.ln_exp_m1().ln_1p_exp(), x);
            }
        }

        #[test]
        fn ln_exp_m1_misc() {
            // limits
            let x: f64 = 0.0;
            assert_eq!(x.ln_exp_m1(), f64::NEG_INFINITY);
            let x: f64 = -0.0;
            assert_eq!(x.ln_exp_m1(), f64::NEG_INFINITY);

            let x: f64 = -5.0e-324;
            assert!(x.ln_exp_m1().is_nan());
            let x: f64 = -f64::MIN_POSITIVE;
            assert!(x.ln_exp_m1().is_nan());

            let x: f64 = f64::INFINITY;
            assert_eq!(x.ln_exp_m1(), f64::INFINITY);

            let x: f64 = (1023.0 * (2.0_f64).ln()) + 1.0;
            assert!(x.ln_exp_m1().is_finite());
        }

        #[test]
        fn logit_identity() {
            let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![5.0e-324, 0.0, 0.5, 1.0 - eps, 1.0 - eps / 2.0, 1.0];
            for x in xs {
                assert_eq!(x.logit().inv_logit(), x);
            }

            // Approximate identity
            let xs: Vec<f64> = vec![f64::MIN_POSITIVE, eps / 2.0, eps, 0.5 - eps, 0.5 + eps];
            let epsilons: Vec<f64> = vec![
                f64::MIN_POSITIVE,
                // eps^(3/2) is a loose but acceptable alternative to eps^2
                eps * eps.sqrt() / 2.0,
                eps * eps.sqrt(),
                eps / 2.0,
                eps / 2.0,
            ];
            for (x, eps) in xs.into_iter().zip(epsilons.into_iter()) {
                assert!((x.logit().inv_logit() - x).abs() < eps);
            }
        }

        #[test]
        fn inv_logit_misc() {
            let x: f64 = f64::INFINITY;
            assert_eq!(x.inv_logit(), 1.0);

            let x: f64 = f64::NEG_INFINITY;
            assert_eq!(x.inv_logit(), 0.0);

            let x: f64 = -745.0;
            assert_eq!(x.inv_logit(), 0.0);

            let x: f64 = -744.0;
            assert!(x.inv_logit() > 0.0);

            let x: f64 = 36.0;
            assert!(x.inv_logit() < 1.0);

            let x: f64 = 37.0;
            assert_eq!(x.inv_logit(), 1.0);

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
            let x: f64 = 5.0e-324;
            let rhs: f64 = -744.4400719213812;
            assert_eq!(x.logit(), rhs);

            let x: f64 = -744.4400719213813;
            assert_eq!(x.inv_logit(), 0.0);
        }

        #[test]
        fn ln_inv_logit_works() {
            let xs: Vec<f64> = vec![f64::INFINITY, f64::NEG_INFINITY, 0.0, -50.0, -37.0];

            for x in xs {
                let lhs = x.ln_inv_logit();
                let rhs = x.inv_logit().ln();
                assert_eq!(lhs, rhs);
            }

            let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![eps, 1.0, 18.0, 33.3, 50.0];

            for x in xs {
                let lhs = x.ln_inv_logit();
                let rhs = x.inv_logit().ln();
                assert!((lhs - rhs).abs() < eps);
            }
        }

        #[test]
        fn ln_inv_logit_misc() {
            let x: f64 = f64::INFINITY;
            assert_eq!(x.ln_inv_logit(), -0.0);
            let x: f64 = f64::NEG_INFINITY;
            assert_eq!(x.ln_inv_logit(), f64::NEG_INFINITY);

            let x: f64 = -745.0;
            assert!(x.ln_inv_logit().is_finite());
            let x: f64 = 50.0;
            assert!(x.ln_inv_logit().is_finite());
            let x: f64 = 745.0;
            assert!(x.ln_inv_logit().is_finite());
        }

        #[test]
        fn logit_exp_works() {
            let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![eps, eps.sqrt(), 0.2, 0.4, 0.8, 1.0 - eps.sqrt(), 1.0 - eps];
            let neg_xs: Vec<f64> = xs.iter().map(|&x| -x).collect();
            for x in xs {
                assert!((x.ln_inv_logit().logit_exp() - x).abs() < eps);
            }
            for x in neg_xs {
                assert!((x.ln_inv_logit().logit_exp() - x).abs() < 2.0 * eps);
            }

            let xs: Vec<f64> = vec![-f64::NEG_INFINITY, 0.0, f64::INFINITY];
            for x in xs {
                assert_eq!(x.ln_inv_logit().logit_exp(), x);
            }
        }

        #[test]
        fn ln_1m_inv_logit_works() {
            let x: f64 = f64::INFINITY;
            assert_eq!(x.ln_1m_inv_logit(), f64::NEG_INFINITY);
            let x: f64 = f64::NEG_INFINITY;
            assert_eq!(x.ln_1m_inv_logit(), -0.0);

            let x: f64 = -745.0;
            assert!(x.ln_1m_inv_logit().is_finite());
            let x: f64 = 50.0;
            assert!(x.ln_1m_inv_logit().is_finite());
            let x: f64 = 745.0;
            assert!(x.ln_1m_inv_logit().is_finite());
        }

        #[test]
        fn logit_1m_exp_works() {
            let eps = f64::EPSILON;
            let xs: Vec<f64> = vec![eps, eps.sqrt(), 0.2, 0.4, 0.8, 1.0 - eps.sqrt(), 1.0 - eps];
            let neg_xs: Vec<f64> = xs.iter().map(|&x| -x).collect();
            for x in xs {
                assert!((x.ln_1m_inv_logit().logit_1m_exp() - x).abs() < 2.0 * eps);
            }
            for x in neg_xs {
                assert!((x.ln_1m_inv_logit().logit_1m_exp() - x).abs() < eps);
            }

            let xs: Vec<f64> = vec![-f64::NEG_INFINITY, 0.0, f64::INFINITY];
            for x in xs {
                assert_eq!(x.ln_1m_inv_logit().logit_1m_exp(), x);
            }
        }
    }

    #[cfg(test)]
    mod f32_impl {
        use super::*;

        #[test]
        fn ln_1m_exp_works() {
            let x: f32 = -1e-20;
            assert!(x.ln_1m_exp().is_finite());

            let x: f32 = -1.0;
            assert_eq!(x.ln_1m_exp(), (-x.exp()).ln_1p());

            let x: f32 = -0.1;
            assert_eq!(x.ln_1m_exp(), (-x.exp_m1()).ln());

            // limits
            let x: f32 = f32::NEG_INFINITY;
            assert_eq!(x.ln_1m_exp(), -0.0); // -0.0 due to log1p(-0.0)

            let x: f32 = 0.0;
            assert_eq!(x.ln_1m_exp(), f32::NEG_INFINITY);

            let x: f32 = -0.0;
            assert_eq!(x.ln_1m_exp(), f32::NEG_INFINITY);

            let x: f32 = 1.0e-45;
            assert!(x.ln_1m_exp().is_nan());
            let x: f32 = f32::MIN_POSITIVE;
            assert!(x.ln_1m_exp().is_nan());

            let x: f32 = -1.0e-45;
            assert!(x.ln_1m_exp().is_finite());
            let x: f32 = -f32::MIN_POSITIVE;
            assert!(x.ln_1m_exp().is_finite());

            // its own inverse
            let x: f32 = -0.5;
            assert_eq!(x.ln_1m_exp().ln_1m_exp(), x);

            let x: f32 = -f32::EPSILON;
            assert!((x.ln_1m_exp().ln_1m_exp() - x).abs() < 6.0e-14);
        }

        #[test]
        fn ln_1p_exp_identity() {
            let xs: Vec<f32> = vec![
                f32::INFINITY,
                f32::NEG_INFINITY,
                0.0,
                2.0,
                3.0,
                -50.0,
                -37.0,
                18.0,
                33.3,
                50.0,
            ];
            for x in xs {
                assert_eq!(x.ln_1p_exp().ln_exp_m1(), x);
            }
            let x: f32 = 1.0;
            assert_eq!(x.ln_1p_exp().ln_exp_m1(), x - f32::EPSILON);
        }

        #[test]
        fn ln_1p_exp_misc() {
            let ln2 = std::f32::consts::LN_2;
            let eps = f32::EPSILON;
            let x: f32 = 0.0;
            assert!((x.ln_1p_exp() - ln2).abs() < eps);

            let e = std::f32::consts::E;
            let x: f32 = 1.0;
            assert!((x.ln_1p_exp() - e.ln_1p()).abs() < eps);

            let x: f32 = -1.0;
            assert!((x.ln_1p_exp() - (e.ln_1p() - 1.0)).abs() < eps);

            // When naive method overflows
            let x: f32 = (127.0 * (2.0_f32).ln()) + 1.0;
            assert!(x.ln_1p_exp().is_finite());

            let xs: Vec<f32> = vec![1e2, 1e3, 1e4, 1e5];
            for x in xs {
                assert!((x.ln_1p_exp() - x).abs() < eps);
            }

            let xs: Vec<f32> = vec![-1e3, -1e4, -1e5];
            for x in xs {
                assert_eq!(x.ln_1p_exp(), 0.0);
            }
        }

        #[test]
        fn ln_exp_m1_identity() {
            let xs: Vec<f32> = vec![f32::INFINITY, 0.0, 1.0, 18.0, 33.3, 50.0];
            for x in xs {
                assert_eq!(x.ln_exp_m1().ln_1p_exp(), x);
            }
        }

        #[test]
        fn ln_exp_m1_misc() {
            // limits
            let x: f32 = 0.0;
            assert_eq!(x.ln_exp_m1(), f32::NEG_INFINITY);
            let x: f32 = -0.0;
            assert_eq!(x.ln_exp_m1(), f32::NEG_INFINITY);

            let x: f32 = -1.0e-45;
            assert!(x.ln_exp_m1().is_nan());
            let x: f32 = -f32::MIN_POSITIVE;
            assert!(x.ln_exp_m1().is_nan());

            let x: f32 = f32::INFINITY;
            assert_eq!(x.ln_exp_m1(), f32::INFINITY);

            // When naive method overflows
            let x: f32 = (127.0 * (2.0_f32).ln()) + 1.0;
            assert!(x.ln_exp_m1().is_finite());
        }

        #[test]
        fn logit_identity() {
            let eps = f32::EPSILON;
            let xs: Vec<f32> = vec![1.0e-45, 0.0, 0.5, 1.0 - eps, 1.0];
            for x in xs {
                assert_eq!(x.logit().inv_logit(), x);
            }

            // Approximate identity
            let xs: Vec<f32> = vec![
                f32::MIN_POSITIVE,
                eps / 2.0,
                eps,
                0.5 - eps,
                0.5 + eps,
                1.0 - eps / 2.0,
            ];
            for x in xs {
                assert!((x.logit().inv_logit() - x).abs() < eps);
            }
        }

        #[test]
        fn inv_logit_misc() {
            let x: f32 = f32::INFINITY;
            assert_eq!(x.inv_logit(), 1.0);

            let x: f32 = f32::NEG_INFINITY;
            assert_eq!(x.inv_logit(), 0.0);

            let x: f32 = -104.0;
            assert_eq!(x.inv_logit(), 0.0);

            let x: f32 = -103.0;
            assert!(x.inv_logit() > 0.0);

            let x: f32 = 16.0;
            assert!(x.inv_logit() < 1.0);

            let x: f32 = 17.0;
            assert_eq!(x.inv_logit(), 1.0);

            let x: f32 = 0.5;
            assert_eq!(x.logit(), 0.0);

            let eps = f32::EPSILON;
            let x: f32 = 2.0;
            assert!((x.inv_logit() - 1.0 / (1.0 + (-2.0_f32).exp())).abs() < eps);
            assert!((x.inv_logit().logit() - 2.0) < eps);
        }

        #[test]
        fn logit_smooth_upper_bound() {
            let x: f32 = 16.63553;
            assert_eq!(x.inv_logit(), 1.0 - f32::EPSILON / 2.0);

            // This is the bound, but, exp(x) / (1.0 + exp(x)) rounds to 1.0
            let x: f32 = 16.635532;
            assert_eq!(x.inv_logit(), 1.0);

            let x: f32 = 16.635534;
            assert_eq!(x.inv_logit(), 1.0);
        }

        #[test]
        fn logit_smooth_lower_bound() {
            let x: f32 = 1.0e-45;
            let rhs: f32 = -103.27893;
            assert_eq!(x.logit(), rhs);

            let x: f32 = -103.27894;
            assert_eq!(x.inv_logit(), 0.0);
        }

        #[test]
        fn ln_inv_logit_works() {
            let xs: Vec<f32> = vec![f32::INFINITY, f32::NEG_INFINITY, 0.0, -50.0, -37.0];

            for x in xs {
                let lhs = x.ln_inv_logit();
                let rhs = x.inv_logit().ln();
                assert_eq!(lhs, rhs);
            }

            let eps = f32::EPSILON;
            let xs: Vec<f32> = vec![eps, 1.0, 18.0, 33.3, 50.0];

            for x in xs {
                let lhs = x.ln_inv_logit();
                let rhs = x.inv_logit().ln();
                assert!((lhs - rhs).abs() < eps);
            }
        }

        #[test]
        fn ln_inv_logit_misc() {
            let x: f32 = f32::INFINITY;
            assert_eq!(x.ln_inv_logit(), -0.0);
            let x: f32 = f32::NEG_INFINITY;
            assert_eq!(x.ln_inv_logit(), f32::NEG_INFINITY);

            let x: f32 = -103.0;
            assert!(x.ln_inv_logit().is_finite());
            let x: f32 = 35.0;
            assert!(x.ln_inv_logit().is_finite());
            let x: f32 = 103.0;
            assert!(x.ln_1m_inv_logit().is_finite());
        }

        #[test]
        fn logit_exp_works() {
            let eps = f32::EPSILON;
            let xs: Vec<f32> = vec![eps, eps.sqrt(), 0.2, 0.4, 0.8, 1.0 - eps.sqrt(), 1.0 - eps];
            let neg_xs: Vec<f32> = xs.iter().map(|&x| -x).collect();
            for x in xs {
                assert!((x.ln_inv_logit().logit_exp() - x).abs() < 2.0 * eps);
            }
            for x in neg_xs {
                assert!((x.ln_inv_logit().logit_exp() - x).abs() < 2.0 * eps);
            }

            let xs: Vec<f32> = vec![-f32::NEG_INFINITY, 0.0, f32::INFINITY];
            for x in xs {
                assert_eq!(x.ln_inv_logit().logit_exp(), x);
            }
        }

        #[test]
        fn ln_1m_inv_logit_works() {
            let x: f32 = f32::INFINITY;
            assert_eq!(x.ln_1m_inv_logit(), f32::NEG_INFINITY);
            let x: f32 = f32::NEG_INFINITY;
            assert_eq!(x.ln_1m_inv_logit(), -0.0);

            let x: f32 = -103.0;
            assert!(x.ln_1m_inv_logit().is_finite());
            let x: f32 = 35.0;
            assert!(x.ln_1m_inv_logit().is_finite());
            let x: f32 = 103.0;
            assert!(x.ln_1m_inv_logit().is_finite());
        }

        #[test]
        fn logit_1m_exp_works() {
            let eps = f32::EPSILON;
            let xs: Vec<f32> = vec![eps, eps.sqrt(), 0.2, 0.4, 0.8, 1.0 - eps.sqrt(), 1.0 - eps];
            let neg_xs: Vec<f32> = xs.iter().map(|&x| -x).collect();
            for x in xs {
                assert!((x.ln_1m_inv_logit().logit_1m_exp() - x).abs() < 2.0 * eps);
            }
            for x in neg_xs {
                assert!((x.ln_1m_inv_logit().logit_1m_exp() - x).abs() < 2.0 * eps);
            }

            let xs: Vec<f32> = vec![-f32::NEG_INFINITY, 0.0, f32::INFINITY];
            for x in xs {
                assert_eq!(x.ln_1m_inv_logit().logit_1m_exp(), x);
            }
        }
    }
}
