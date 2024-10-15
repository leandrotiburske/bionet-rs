use ndarray::{Array1, Array2};
use petgraph::prelude::*;
use std::collections::HashMap;

pub struct CoexpressionNetwork {
    pub network: UnGraph<String, f64>,
}

impl CoexpressionNetwork {

    pub fn new() -> Self {
        CoexpressionNetwork {
            network: UnGraph::<String, f64>::new_undirected(),
        }    
    }

    pub fn create_graph(&mut self, correlation: Array2<f64>, gene_names: Array1<String>) {

        let mut node_indices = HashMap::new(); // Map to keep track of node indices
        let (rows, cols) = correlation.dim();
    
        for i in 0..rows {
    
            let node_index = self.network.add_node(gene_names[i].clone());
            node_indices.insert(gene_names[i].clone(), node_index); // Store the index
    
            for j in 0..cols {
                if correlation[[i, j]] != 0.0 {
                    if let Some(&target_index) = node_indices.get(&gene_names[j]) {
                        // Add the edge with the correlation value as the weight
                        self.network.add_edge(node_index, target_index, correlation[[i, j]]);
                    }
                }
            }
        }    
    }

}