use csv::ReaderBuilder;
use ndarray::{Array1, Array2};
use std::error::Error;
use std::fs::File;

pub fn read_abundance(filepath: String) -> Result<Array2<f64>, Box<dyn Error>> {
    
    //
    // Open the CSV file
    let file = File::open(filepath)?;
    //

    //
    // Create a CSV reader
    //
    let mut reader = ReaderBuilder::new().from_reader(file);
    let mut values: Vec<Vec<f64>> = Vec::new();

    //
    // Read and skip the header
    //
    for result in reader.records() {
        let record = result?;
        // Parse the values and collect them into a row vector
        let row: Vec<f64> = record.iter().skip(1)
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        values.push(row);
    }

    //
    // Convert Vec<Vec<f64>> to Array2<f64>
    //
    let num_rows = values.len();
    let num_cols = values[0].len();
    let array2 = Array2::from_shape_vec((num_rows, num_cols), values
        .into_iter()
        .flat_map(|v| v)
        .collect())?;

    // Return Array2
    Ok(array2)
}


pub fn read_gene_names(filepath: String) -> Result<Array1<String>, Box<dyn Error>> {

    // Open the CSV file
    let file = File::open(filepath)?;
    //

    //
    // Create a CSV reader
    //
    let mut reader = ReaderBuilder::new().from_reader(file);
    let mut values: Vec<Vec<String>> = Vec::new();

    //
    // Read and skip the header
    //
    for result in reader.records() {
        let record = result?;
        // Parse the values and collect them into a row vector
        let row: Vec<String> = vec![record.iter().next().unwrap().to_owned()];
        values.push(row);
    }

    //
    // Convert Vec<Vec<f64>> to Array2<f64>
    //
    let num_rows = values.len();
    let array1 = Array1::from_shape_vec(num_rows, values
        .into_iter()
        .flat_map(|v| v)
        .collect())?;

    // Return Array2
    Ok(array1)

}