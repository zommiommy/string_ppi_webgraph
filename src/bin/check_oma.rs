
use std::io::prelude::*;
use std::io;
use std::fs;
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

pub fn main() -> Result<()> {
    check_oma_groups_prefixes()?;
    Ok(())
}