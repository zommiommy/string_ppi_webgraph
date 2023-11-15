//! Scan all the files and build an unique vocabulary that maps all node names
//! to node ids.

use std::io::prelude::*;
use std::io;
use std::fs;
use dsi_progress_logger::*;
use flate2::read::GzDecoder;
use indicatif::ProgressIterator;
use std::collections::BTreeMap;
use anyhow::Result;
use indicatif::ProgressBar;

fn parse_oma_groups(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-groups.txt.gz");
    let file = fs::File::open("../oma-groups.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for candidate_oma_entry in line.split('\t').skip(2) {
            let node_id = vocab.len();
            vocab.entry(candidate_oma_entry.to_uppercase()).or_insert(node_id);
        }
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_oma_species(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-species.txt");
    let file = fs::File::open("../oma-species.txt")?;
    let gz = io::BufReader::new(file);

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let oma_code = vals[0];
        let ncbi_code = format!("NCBITAXON:{}", vals[2]);
        
        let node_id = vocab.len();
        vocab.entry(oma_code.to_uppercase()).or_insert(node_id);
        let node_id = vocab.len();
        vocab.entry(ncbi_code.to_uppercase()).or_insert(node_id);
        pl.light_update();
    }
    pl.done();

    Ok(())
}

fn parse_oma_uniprot(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-uniprot.txt.gz");
    let file = fs::File::open("../oma-uniprot.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for node_name in line.split('\t') {
            let node_id = vocab.len();
            vocab.entry(node_name.to_uppercase()).or_insert(node_id);
        }
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_string_aliases(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on protein.aliases.v12.0.txt.gz");
    let file = fs::File::open("../protein.aliases.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().progress_with(ProgressBar::new_spinner()) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let source = vals[2].to_uppercase();
        if source != "UNIPROT_AC" {
            continue;
        }
        for node_name in line.split('\t').take(2) {
            let node_id = vocab.len();
            vocab.entry(node_name.to_uppercase()).or_insert(node_id);
        }
        pl.light_update();
    }
    pl.done();
    Ok(())
}

const ENRICHMENT_FILTER: &[&str] = &[
    "BTO","CL","DOID","FBCV","GO","HP","MP","ZP"
];

fn parse_string_enrichment_terms(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on enrichment.terms.v12.0.txt.gz");
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
        vocab.entry(string_protein_id.to_uppercase()).or_insert(node_id);

        
        let term = vals[2].to_uppercase();

        for term_filter in ENRICHMENT_FILTER {
            if term.starts_with(term_filter) {
                let node_id = vocab.len();
                vocab.entry(term.to_string()).or_insert(node_id);
                break;
            }
        }
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_string_links(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on protein.links.full.v12.0.txt.gz");
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
        vocab.entry(src.to_uppercase()).or_insert(node_id);

        
        let dst = vals[1];
        let node_id = vocab.len();
        vocab.entry(dst.to_uppercase()).or_insert(node_id);
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_kgx_nodelist(vocab: &mut BTreeMap<String, usize>, file: &str) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start(format!("Working on {}", file));
    assert!(file.ends_with(".tsv"));
    let file = fs::File::open(format!("../{}", file))?;
    let gz = io::BufReader::new(file);

    let mut lines_iter = gz.lines();

    let header = lines_iter.next().unwrap()?;
    let vals: Vec<&str> = header.split('\t').collect::<Vec<_>>();
    assert_eq!(vals[0], "id");

    for line in lines_iter {
        let line = line?;

        let node_name = line.split('\t').next().unwrap();
        let node_id = vocab.len();
        vocab.entry(node_name.to_uppercase()).or_insert(node_id);

        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_eggnog_groups(vocab: &mut BTreeMap<String, usize>) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on e6.og2seqs_and_species.tsv");
    let file = fs::File::open("../e6.og2seqs_and_species.tsv")?;
    let gz = io::BufReader::new(file);

    // skip header
    for line in gz.lines() {
        let line = line?;

        let vals = line.split('\t').collect::<Vec<_>>();
        let node_name = format!("EGG:{}", vals[1]);
        let node_id = vocab.len();
        vocab.entry(node_name.to_uppercase()).or_insert(node_id);

        pl.light_update();
    }
    pl.done();
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

const KGX_FILES: &[&str] = &[
    "ncbitaxon_kgx_tsv_nodes.tsv",
    "go_kgx_tsv_nodes.tsv",
    "bto_kgx_tsv_nodes.tsv",
    "cl_kgx_tsv_nodes.tsv",
    "doid_kgx_tsv_nodes.tsv",
    "fbcv_kgx_tsv_nodes.tsv",
    "hp_kgx_tsv_nodes.tsv",
    "mp_kgx_tsv_nodes.tsv",
    "zp_kgx_tsv_nodes.tsv",
];

pub fn main() -> Result<()> {
    let mut vocab = BTreeMap::new();
    parse_oma_groups(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_oma_species(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_oma_uniprot(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_eggnog_groups(&mut vocab)?;
    println!("vocab size: {}", vocab.len());

    for file in KGX_FILES {
        parse_kgx_nodelist(&mut vocab, file)?;
        println!("vocab size: {}", vocab.len());
    }

    parse_string_aliases(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_string_enrichment_terms(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    parse_string_links(&mut vocab)?;
    println!("vocab size: {}", vocab.len());
    dump_vocab(&vocab)?;

    Ok(())
}