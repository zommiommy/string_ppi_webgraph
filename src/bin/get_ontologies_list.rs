//! Document to get the unique list of ontology prefixes from the STRING enrichment terms file.
//!
//! These ontology prefixes will look like:
//! * GO: Gene Ontology
//! * CL: Cell Ontology
//!
//! The list of ontology prefixes will be saved in a JSON file.

use anyhow::Result;
use dsi_progress_logger::*;
use flate2::read::GzDecoder;
use std::collections::BTreeSet;
use std::fs;
use std::io;
use std::io::prelude::*;

fn parse_string_enrichment_terms(ontologie_codes: &mut BTreeSet<String>) -> Result<()> {
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on protein.enrichment.terms.v12.0.txt.gz");
    let file = fs::File::open("../protein.enrichment.terms.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        // We get the ontology node name
        let ontology_node_name = line.split('\t').take(3).collect::<Vec<_>>()[2];

        // We get the prefix of the ontology code
        let ontology_code = if ontology_node_name.contains('-') {
            ontology_node_name.split('-').collect::<Vec<_>>()[0]
        } else if ontology_node_name.starts_with("IPR") {
            "IPR"
        } else if ontology_node_name.starts_with("PF") {
            "PF"
        } else if ontology_node_name.starts_with("SM") {
            "SM"
        } else if ontology_node_name.starts_with("WP") {
            "WP"
        } else if ontology_node_name.contains(':') {
            ontology_node_name.split(':').collect::<Vec<_>>()[0]
        } else {
            unreachable!("how to split {}", ontology_node_name)
        };

        // We insert the ontology code in the set
        ontologie_codes.insert(ontology_code.to_string());
        pl.light_update();
    }
    pl.done();
    Ok(())
}

pub fn main() -> Result<()> {
    stderrlog::new()
        .verbosity(2)
        .timestamp(stderrlog::Timestamp::Second)
        .init()
        .unwrap();

    // Initialize the ontology codes set
    let mut ontologie_codes = BTreeSet::new();

    // We parse the STRING enrichment terms file
    parse_string_enrichment_terms(&mut ontologie_codes)?;

    // We save the ontology codes as a JSON list
    let mut ontology_codes_file = fs::File::create("ontology_codes.json")?;
    ontology_codes_file.write_all(serde_json::to_string(&ontologie_codes)?.as_bytes())?;

    Ok(())
}
