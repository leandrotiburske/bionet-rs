use bionet_rs::*;
use ndarray_stats::CorrelationExt;

fn main() {
    let abundance = read_array().unwrap();
    let correlation = abundance.pearson_correlation().unwrap();
    let n = abundance.shape()[1]; // Number of samples = number of columns
    let mut p_values = pearson_to_pvalue(&correlation, n);
    p_values = adjust_pvalues(&mut p_values);
    let corrected_corr = filter_corr(correlation, &p_values);
}