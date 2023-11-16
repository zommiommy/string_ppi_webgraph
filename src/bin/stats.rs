use anyhow::Result;
use dsi_progress_logger::*;
use flate2::read::GzDecoder;
use std::collections::BTreeMap;
use std::fs;
use std::io;
use std::io::prelude::*;

fn eggnogg_stats() -> Result<()> {
    // check that all OMA groups are in the species file
    let mut pl = ProgressLogger::default();
    pl.display_memory(true);
    pl.start("Working on e6.og2seqs_and_species.tsv");
    let mut dist_raw: BTreeMap<usize, usize> = BTreeMap::new();
    let mut dist_edges: BTreeMap<usize, usize> = BTreeMap::new();

    let file = fs::File::open("../e6.og2seqs_and_species.tsv")?;
    let gz = io::BufReader::new(file);

    for line in gz.lines() {
        let line = line?;
        let vals: Vec<&str> = line.split('\t').collect::<Vec<_>>();

        let group_len = vals[3].parse::<usize>()?;
        let num_edges = group_len + (group_len * (group_len - 1)) + 1;

        *dist_raw.entry(group_len).or_insert(0) += 1;
        *dist_edges.entry(num_edges).or_insert(0) += 1;
    }
    pl.done();

    let mut file = io::BufWriter::new(fs::File::create("egg_stats.tsv")?);
    for (k, v) in dist_raw.iter() {
        writeln!(file, "{}\t{}", k, v)?;
    }
    let mut file = io::BufWriter::new(fs::File::create("egg_stats_edges.tsv")?);
    for (k, v) in dist_edges.iter() {
        writeln!(file, "{}\t{}", k, v)?;
    }
    Ok(())
}

fn oma_stats() -> Result<()> {
    let mut dist_raw: BTreeMap<usize, usize> = BTreeMap::new();
    let mut dist_edges: BTreeMap<usize, usize> = BTreeMap::new();

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
        let group_len = line.split('\t').skip(2).count();
        let num_edges = group_len + (group_len * (group_len - 1)) + 1;

        *dist_raw.entry(group_len).or_insert(0) += 1;
        *dist_edges.entry(num_edges).or_insert(0) += 1;
    }
    pl.done();

    let mut file = io::BufWriter::new(fs::File::create("oma_stats.tsv")?);
    for (k, v) in dist_raw.iter() {
        writeln!(file, "{}\t{}", k, v)?;
    }
    let mut file = io::BufWriter::new(fs::File::create("oma_stats_edges.tsv")?);
    for (k, v) in dist_edges.iter() {
        writeln!(file, "{}\t{}", k, v)?;
    }
    Ok(())
}

pub fn main() -> Result<()> {
    stderrlog::new()
        .verbosity(2)
        .timestamp(stderrlog::Timestamp::Second)
        .init()?;

    oma_stats()?;
    eggnogg_stats()?;
    Ok(())
}
