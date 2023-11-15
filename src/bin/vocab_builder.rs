//! Scan all the files and build an unique vocabulary that maps all node names
//! to node ids.

use std::io::prelude::*;
use std::io;
use std::fs;
use flate2::read::GzDecoder;
use indicatif::ProgressIterator;
use std::collections::BTreeMap;
use anyhow::Result;
use indicatif::ProgressBar;

fn parse_oma_groups(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    println!("Working on oma-groups.txt.gz");
    let file = fs::File::open("../oma-groups.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for candidate_oma_entry in line.split('\t').skip(2) {
            let node_id = vocab.len();
            vocab.entry(candidate_oma_entry.to_string()).or_insert(node_id);
        }
    }
    Ok(())
}

fn parse_oma_species(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    println!("Working on oma-species.txt");
    let file = fs::File::open("../oma-species.txt")?;
    let gz = io::BufReader::new(file);

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let oma_code = vals[0];
        let ncbi_code = format!("NCBI:{}", vals[2]);
        
        let node_id = vocab.len();
        vocab.entry(oma_code.to_string()).or_insert(node_id);
        let node_id = vocab.len();
        vocab.entry(ncbi_code).or_insert(node_id);
    }

    Ok(())
}

fn parse_oma_uniprot(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    println!("Working on oma-uniprot.txt.gz");
    let file = fs::File::open("../oma-uniprot.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for node_name in line.split('\t') {
            let node_id = vocab.len();
            vocab.entry(node_name.to_string()).or_insert(node_id);
        }
    }
    Ok(())
}

fn parse_string_aliases(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    println!("Working on protein.aliases.v12.0.txt.gz");
    let file = fs::File::open("../protein.aliases.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let source = vals[2];
        if source != "UniProt_AC" {
            continue;
        }
        for node_name in line.split('\t').take(2) {
            let node_id = vocab.len();
            vocab.entry(node_name.to_string()).or_insert(node_id);
        }
    }
    Ok(())
}

fn parse_string_enrichment_terms(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    println!("Working on protein.enrichment.terms.v12.0.txt.gz");
    let file = fs::File::open("../protein.enrichment.terms.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();

        let string_protein_id = vals[0];
        let node_id = vocab.len();
        vocab.entry(string_protein_id.to_string()).or_insert(node_id);

        
        let go_term = vals[2];
        let node_id = vocab.len();
        vocab.entry(go_term.to_string()).or_insert(node_id);
    }
    Ok(())
}

fn parse_string_links(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    println!("Working on protein.links.full.v12.0.txt.gz");
    let file = fs::File::open("../protein.links.full.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().skip(1).progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split(' ').collect::<Vec<_>>();
        let src = vals[0];
        let node_id = vocab.len();
        vocab.entry(src.to_string()).or_insert(node_id);

        
        let dst = vals[1];
        let node_id = vocab.len();
        vocab.entry(dst.to_string()).or_insert(node_id);
    }
    Ok(())
}

fn dump_vocab(vocab: &BTreeMap<String, usize>) -> Result<()> {
    let mut vocab_file = io::BufWriter::new(fs::File::create("../vocab.sorted.tsv")?);
    for (node_name, node_id) in vocab.iter().progress() {
        writeln!(vocab_file, "{}\t{}", node_name, node_id)?;
    }

    let mut vocabs = vocab.iter().map(|(node_name, node_id)| {
        (node_id, node_name)
    }).collect::<Vec<_>>();
    vocabs.sort_by_key(|(node_id, _)| *node_id);

    let mut vocab_file = io::BufWriter::new(fs::File::create("../vocab.tsv")?);
    for (_, node_name) in vocabs.iter().progress() {
        writeln!(vocab_file, "{}", node_name)?;
    }
    Ok(())
}

pub fn main() -> Result<()> {
    let mut vocab = BTreeMap::new();
    parse_oma_groups(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_oma_species(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_oma_uniprot(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_string_aliases(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_string_enrichment_terms(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_string_links(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    dump_vocab(&vocab)?;

    Ok(())
}