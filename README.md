# lnexp

## Description

Provides a trait (`LnExp`) for floating-point types to perform careful
evaluation of compositions of `ln`, `ln_1p`, `exp` and `exp_m1`.
Implementations are provided for `f64` and `f32`; see the
documentation for details.

## License

Licensed under either of

  * [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)
  * [MIT license](http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.

## Citations
- [Martin Maechler (2012), Accurately Computing log(1 − exp(− |a|))](http://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf)
- [cossio (2022), untitled](https://github.com/JuliaStats/LogExpFunctions.jl/files/8218470/log1pexp.pdf)
