use bionet_rs::*;
use ndarray_stats::CorrelationExt;

fn main() {
    let abundance = read_abundance("/home/leandro/Documents/abundance.csv".to_string()).unwrap();
    let gene_names = read_gene_names("/home/leandro/Documents/abundance.csv".to_string()).unwrap();

    let correlation = abundance.pearson_correlation().unwrap();

    let n = abundance.shape()[1]; // Number of samples = number of columns
    let mut p_values = pearson_to_pvalue(&correlation, n);
    p_values = adjust_pvalues(&mut p_values);

    let corrected_corr = filter_corr(correlation, &p_values);

    let mut network = CoexpressionNetwork::new();
    network.create_graph(corrected_corr, gene_names);
    println!("{:?}", network.network);

    println!("{:?}", network.find_communities())

}