//! Using the pre-built vocab.tsv, merge all files into a graph
//!
//! you may want to increase the maximum number of open files, expecially if 
//! your batch_size is small, i.e. you have little RAM
//! ```
//! sudo sysctl -w fs.file-max=10000000
//! ```


use std::collections::BTreeMap;
use std::io::prelude::*;
use std::io;
use std::fs;
use anyhow::Result;
use flate2::read::GzDecoder;
use dsi_progress_logger::*;

use webgraph::prelude::*;
use webgraph::graph::arc_list_graph::ArcListGraph;
use webgraph::graph::bvgraph::parallel_compress_sequential_iter;

use rand::Rng;
use std::path::Path;
use itertools::{Dedup, Itertools};

/// Create a new random dir inside the given folder
pub fn temp_dir<P: AsRef<Path>>(base: P) -> String {
    let mut base = base.as_ref().to_owned();
    const ALPHABET: &[u8] = b"0123456789abcdef";
    let mut rnd = rand::thread_rng();
    let mut random_str = String::new();
    loop {
        random_str.clear();
        for _ in 0..16 {
            let idx = rnd.gen_range(0..ALPHABET.len());
            random_str.push(ALPHABET[idx] as char);
        }
        base.push(&random_str);

        if !base.exists() {
            std::fs::create_dir(&base).unwrap();
            return base.to_string_lossy().to_string();
        }
        base.pop();
    }
}

fn parse_oma_groups(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-groups.txt.gz");
    let file = fs::File::open("../oma-groups.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for src in line.split('\t').skip(2) {

            let src_id = vocab.get(src).unwrap();
            let src_prefix = vocab.get(&src[..5]).unwrap();
            sorted.push(*src_prefix, *src_id)?;
            pl.light_update();

            for dst in line.split('\t').skip(2) {
                if src == dst {
                    continue;
                }
                let dst_id = vocab.get(dst).unwrap();
                sorted.push(*src_id, *dst_id)?;
                pl.light_update();
            }
        }
    }
    pl.done();
    Ok(())
}

fn parse_oma_species(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-species.txt");
    let file = fs::File::open("../oma-species.txt")?;
    let gz = io::BufReader::new(file);

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let oma_code = vals[0];
        let ncbi_code = format!("NCBI:{}", vals[2]);
        
        let oma_code = vocab.get(oma_code).unwrap();
        let ncbi_code = vocab.get(&ncbi_code).unwrap();
        sorted.push(*oma_code, *ncbi_code)?;
        sorted.push(*ncbi_code, *oma_code)?;
        pl.light_update();
    }
    pl.done();

    Ok(())
}

fn parse_oma_uniprot(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on oma-uniprot.txt.gz");
    let file = fs::File::open("../oma-uniprot.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let oma_code = vals[0];
        let uniprot_code = vals[1];
        let oma_code_id = vocab.get(oma_code).unwrap();
        let uniprot_code_id = vocab.get(uniprot_code).unwrap();

        sorted.push(*oma_code_id, *uniprot_code_id)?;
        sorted.push(*uniprot_code_id, *oma_code_id)?;
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_string_aliases(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on protein.aliases.v12.0.txt.gz");
    let file = fs::File::open("../protein.aliases.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split('\t').collect::<Vec<_>>();
        let source = vals[2];
        if source != "UniProt_AC" {
            continue;
        }
        let string_code = vals[0];
        let string_code_id = vocab.get(string_code).unwrap();
        let uniprot_code = vals[1];
        let uniprot_code_id = vocab.get(uniprot_code).unwrap();
        sorted.push(*string_code_id, *uniprot_code_id)?;
        pl.light_update();
        sorted.push(*uniprot_code_id, *string_code_id)?;
        pl.light_update();
    }
    pl.done();
    Ok(())
}


fn parse_string_enrichment_terms(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    // check that all OMA groups are in the species file
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
        let vals = line.split('\t').collect::<Vec<_>>();

        let string_protein = vals[0];
        let string_protein_id = vocab.get(string_protein).unwrap();

        let go_term = vals[2];
        let go_term_id = vocab.get(go_term).unwrap();

        sorted.push(*string_protein_id, *go_term_id)?;
        pl.light_update();
        sorted.push(*go_term_id, *string_protein_id)?;
        pl.light_update();
    }
    pl.done();
    Ok(())
}

fn parse_string_links(vocab: &BTreeMap<String, usize>, sorted: &mut SortPairs) -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on protein.links.full.v12.0.txt.gz");
    let file = fs::File::open("../protein.links.full.v12.0.txt.gz")?;
    let gz = io::BufReader::new(GzDecoder::new(io::BufReader::new(file)));

    for line in gz.lines().skip(1) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let vals = line.split(' ').collect::<Vec<_>>();

        let src = vals[0];
        let src_id = vocab.get(src).unwrap();
        
        let dst = vals[1];
        let dst_id = vocab.get(dst).unwrap();
        
        // this file ***SHOULD*** be already undirected
        sorted.push(*src_id, *dst_id)?;
        pl.light_update();
        //sorted.push(*dst_id, *src_id)?;
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
    
    // load the vocab
    let mut vocab = BTreeMap::new();
    let f = std::io::BufReader::new(std::fs::File::open("../vocab.tsv")?);
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Loading vocab");
    for (node_id, num_node) in f.lines().enumerate() {
        let num_node = num_node?;
        vocab.insert(num_node, node_id);
        pl.light_update();
    }
    pl.done();
    let num_nodes = vocab.len();
    // a batch is 16GBs
    let mut sorted = SortPairs::new(1_000_000_000, temp_dir("/dfd/tmp"))?;

    parse_oma_uniprot(&vocab, &mut sorted)?;
    parse_oma_species(&vocab, &mut sorted)?;
    parse_oma_groups(&vocab, &mut sorted)?;
    parse_string_aliases(&vocab, &mut sorted)?;
    parse_string_enrichment_terms(&vocab, &mut sorted)?;
    parse_string_links(&vocab, &mut sorted)?;

    // conver the iter to a graph
    let g = ArcListGraph::new(
        num_nodes,
        sorted.iter().unwrap().map(|(src, dst, _)| (src, dst)).dedup(),
    );
    // compress it
    parallel_compress_sequential_iter::<&ArcListGraph<Dedup<std::iter::Map<KMergeIters<_>, _>>>, _>(
        "../res",
        &g,
        num_nodes,
        CompFlags::default(),
        1,
        temp_dir("/dfd/tmp"),
    )
    .unwrap();

    Ok(())
}