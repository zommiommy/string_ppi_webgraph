
use std::io::prelude::*;
use std::io;
use std::fs;
use dsi_progress_logger::*;
use flate2::read::GzDecoder;
use anyhow::Result;
use std::collections::BTreeSet;


/// Check that all OMA entry ids in the oma-groups file are prefixed by codes
/// in the species file
fn check_oma_groups_prefixes() -> Result<()> {
    // load the OMA species codes
    let file = fs::File::open("../oma-species.txt")?;
    let gz = io::BufReader::new(file);
    let mut oma_codes = BTreeSet::new();

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let oma_code = vals[0];
        oma_codes.insert(oma_code.to_string());
    }

    // check that all OMA groups are in the species file
    let file = fs::File::open("../oma-groups.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for candidate_oma_entry in line.split('\t').skip(2) {
            if !oma_codes.contains(&candidate_oma_entry[..5]) {
                panic!("{}", &candidate_oma_entry[..5]);
            }
        }
    }
    
    Ok(())
}

/// check that in the eggnogg ortholog tsv, all the comma separated ids in the 
/// last column are valid STRING ppi ids
fn check_eggnog() -> Result<()> {
    let file = fs::File::open("../protein.info.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));
    
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Loading protein.info.v12.0.txt.gz");
    let mut ppi_ids = BTreeSet::new();
    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let ppi_id = vals[0];
        ppi_ids.insert(ppi_id.to_string());
        pl.light_update();
    }
    pl.done();
    
    let file = fs::File::open("../e6.og2seqs_and_species.tsv")?;
    let gz = io::BufReader::new(file);

    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Checking e6.og2seqs_and_species.tsv");
    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let ppi_vals = vals.last().unwrap();
        
        for maybe_string_id in ppi_vals.split(',') {
            if !ppi_ids.contains(maybe_string_id) {
                panic!("{}", maybe_string_id);
            }
        }
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

    check_oma_groups_prefixes()?;
    check_eggnog()?;
    Ok(())
}