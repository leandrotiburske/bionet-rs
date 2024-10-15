use ndarray::Array2;
use statrs::distribution::{StudentsT, ContinuousCDF};
use rayon::prelude::*;
use std::sync::{Arc, Mutex};

pub fn pearson_to_pvalue(correlations: &Array2<f64>, n: usize) -> Array2<f64> {
    let (rows, cols) = correlations.dim();
    let p_values = Arc::new(Mutex::new(Array2::zeros((rows, cols))));

    // Parallel iteration over rows
    (0..rows).into_par_iter().for_each(|i| {

        let mut row_p_values = Array2::zeros((1, cols)); 

        for j in 0..cols {
            let r = correlations[[i, j]];

            if r.is_nan() || r < -1.0 || r > 1.0 {
                row_p_values[[0, j]] = 1.0;
                continue;
            }

            // Calculate t statistic
            let t = r * ((n as f64 - 2.0).sqrt() / (1.0 - r * r).sqrt());

            // Degrees of freedom
            let degrees_freedom = (n as f64 - 2.0) as usize;

            // Create t-distribution
            let t_dist = StudentsT::new(0.0, 1.0, degrees_freedom as f64).unwrap();

            // Calculate p-value (two-tailed)
            let p_value = 2.0 * t_dist.cdf(-t.abs());
            row_p_values[[0, j]] = p_value as f64;

        }
        
        // Lock and store the result for this row
        
        let mut p_values = p_values.lock().unwrap();
        for j in 0..cols { 
            p_values[[i, j]] = row_p_values[[0, j]]; 
        }
        
    });

    Arc::try_unwrap(p_values).unwrap().into_inner().unwrap()
}

pub fn adjust_pvalues(p_values: &mut Array2<f64>) -> Array2<f64> {

    let (rows, cols) = p_values.dim();
    for i in 0..rows {
        for j in 0..cols {

            if ( p_values[[i, j]] * rows as f64 ) <= 0.05 {
                p_values[[i, j]] = p_values[[i, j]] * rows as f64;
            }

            else {
                p_values[[i, j]] = 1.0;
            }
        }
    }
    p_values.clone()
}

pub fn filter_corr(mut correlation: Array2<f64>, p_values: &Array2<f64>) -> Array2<f64> {

    let (rows, cols) = correlation.dim();
    for i in 0..rows {
        for j in 0..cols {

            if i == j {
                correlation[[i, j]] = 0.0;

            } else if p_values[[i, j]] < 0.05 && ( correlation[[i, j]] < -0.8 ||  correlation[[i, j]] > 0.8 ) {
                continue;

            } else {
                correlation[[i, j]] = 0.0;
            }
        }
    }
    correlation
}