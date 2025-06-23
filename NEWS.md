# geeCRT 1.1.4

## Major changes
1. Improved Numerical Stability: Addressed a potential numerical instability that caused a LAPACK error during checks on a specific CRAN platform (R-devel with OpenBLAS). While the package passed checks with standard libraries, we have proactively enhanced the functions.

2. Enhanced Generalized Inverse: To ensure robustness, the internal generalized inverse (ginv) calculations have been improved with optimized stability checks and matrix rescaling to better handle ill-conditioned matrices.