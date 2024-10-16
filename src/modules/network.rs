// Importing necessary modules from external crates
use ndarray::{Array1, Array2};
use petgraph::graph::Node;
use petgraph::prelude::*;
use std::collections::HashMap;
use petgraph::algo::kosaraju_scc;
use petgraph::visit::NodeFiltered;

// CoexpressionNetwork struct holds an undirected graph where nodes represent genes and edges represent correlations

pub struct CoexpressionNetwork {
    pub network: UnGraph<String, f64>,
}

impl CoexpressionNetwork {
    ///
    /// Creates a new, empty CoexpressionNetwork.
    ///
    /// # Returns
    /// A `CoexpressionNetwork` with an empty undirected graph.

    pub fn new() -> Self {
        CoexpressionNetwork {
            network: UnGraph::<String, f64>::new_undirected(),
        }    
    }

    /// Populates the graph with nodes representing genes and edges representing correlations between them.
    ///
    /// # Arguments
    /// * `correlation` - A 2D array of floating-point numbers representing the correlation matrix.
    /// * `gene_names` - A 1D array of gene names corresponding to the rows/columns of the correlation matrix.
    ///
    /// The function adds nodes for each gene and edges between genes whose correlation is non-zero.

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

    /// Find strongly connected components by applying the Korasaju's algorithm
    ///
    ///
    /// # Returns
    /// Vector (Vec<Vec<NodeIndex>>) of gene communities that are not composed of a single gene

    pub fn find_communities(&self) -> Vec<Vec<NodeIndex>> {

        let mut communities = kosaraju_scc(&self.network);

        let filtered_graph = NodeFiltered::from_fn(&self.network, |node| {
            communities.iter().any(|community| community.contains(&node))
        });

        communities.retain(|v| v.len() > 1); 
        communities 
    }

}