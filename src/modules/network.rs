// Importing necessary modules from external crates

use ndarray::{Array1, Array2};
use petgraph::prelude::*;
use webgestalt_lib::{methods::ora::{get_ora, ORAConfig}, readers::read_gmt_file};
use std::collections::HashMap;
use petgraph::algo::kosaraju_scc;
use ahash::{AHashMap, AHashSet};


///
/// CoexpressionNetwork struct holds an undirected graph where nodes represent genes and edges represent correlations
/// 

pub struct CoexpressionNetwork {
    pub network: UnGraph<String, f64>,
    pub communities: Option<Vec<Vec<NodeIndex>>>,
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
            communities: None,
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

        let mut node_indices: HashMap<String, NodeIndex> = HashMap::new(); // Map to keep track of node indices
        let (rows, cols) = correlation.dim();
    
        for i in 0..rows {
    
            let node_index: NodeIndex = self.network.add_node(gene_names[i].clone());
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

    pub fn find_communities(&mut self) {

        let mut communities: Vec<Vec<NodeIndex>>  = kosaraju_scc(&self.network);

        communities.retain(|v: &Vec<NodeIndex>| v.len() > 1);
        self.communities = Some(communities);

    }

    /// Perform over-representation analysis for each gene community
    /// 

    fn perform_ora(&mut self) {

        // Perform ORA with webgestalt_lib

        let communities: Vec<Vec<NodeIndex>> = self.communities.clone().unwrap();
        let communities: AHashSet<Vec<NodeIndex>> = communities.into_iter().collect();

        let gmt = match read_gmt_file("/home/leandro/Documents/teste.gmt".to_string()) {
            Ok(gmt) => gmt,
            Err(e) => {
                eprintln!("Failed to read GMT file: {:?}", e);
                return;
            }
        };
        println!("{:?}", gmt);

        let communities: AHashSet<Vec<String>> = communities
        .into_iter()
        .map(|vec_nodeindex| {
            vec_nodeindex
                .into_iter()
                .map(|node| self.network[node].clone())  // Convert NodeIndex to String
                .collect::<Vec<String>>()                                                        // Collect into Vec<String>
        })
        .collect();  // Collect into AHashSet<Vec<String>>

        println!("{:?}", communities);

        let teste: Vec<String> = vec!["teste".to_string()];
        let reference: &AHashSet<String> = &teste.into_iter().collect();

        for community in &communities {

            let community = AHashSet::from_iter(community.to_owned());

            let ora_config: ORAConfig = ORAConfig {
                min_overlap: 0,
                max_set_size: 500,
                min_set_size: 2,
                fdr_method: webgestalt_lib::stat::AdjustmentMethod::BH,
            };

            let ora = get_ora(&community, reference, gmt.clone(), ora_config);
            
            // Check if there are any results
            if ora.is_empty() {
                eprintln!("ORA returned no results for this community");
            } else {
                // Print the ORA results
                for result in ora {
                    println!("ORA result: {:?}", result);
                }
            }

        }

    }

    pub fn ora(&mut self) {

        match self.communities {

            None => println!("No communities found. Please run `find_communities()` first."),
            Some(_) => self.perform_ora(),

        }

    }

}