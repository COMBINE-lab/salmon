use schema::{CellEQClass, CellExp, ProcessedExp, SCEQClass};

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io;
use std::io::Write;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::SplitWhitespace;

use bitvec::*;
use hashers::fnv::FNV1aHasher64;
use levenshtein::levenshtein;
use petgraph;
use petgraph::prelude::*;

pub fn read_eqclass_from_iter(elem_iter: &mut SplitWhitespace) -> Option<SCEQClass> {
    let parse_u32 = |x: &mut SplitWhitespace| x.next().expect("parse 32 error").parse::<u32>().expect("parse 32 error");
    let parse_u64 = |x: &mut SplitWhitespace| x.next().expect("parse 64 error").parse::<u64>().expect("parse 64 error");

    // first element is group size
    let gsize: u64 = parse_u64(elem_iter);
    let mut tvec: Vec<u32> = Vec::with_capacity(gsize as usize);
    for _ in 1..=gsize {
        tvec.push(parse_u32(elem_iter));
    }
    // The total number of reads in the eq class?
    let _count = parse_u64(elem_iter);
    let bsize = parse_u64(elem_iter);
    let mut cells: Vec<u32> = Vec::with_capacity(bsize as usize);
    let mut cell_count: Vec<u32> = Vec::with_capacity(bsize as usize);
    let mut umis: Vec<(String, u32)> = Vec::with_capacity(bsize as usize);
    cell_count.push(0);

    for _ in 1..=bsize {
        // barcode ID
        let bc = parse_u32(elem_iter);
        cells.push(bc);

        // size of ugroup
        let ugsize = parse_u32(elem_iter);
        let offset = ugsize + cell_count.last().unwrap_or(&0);
        cell_count.push(offset);

        for _ in 1..=ugsize {
            let umi = elem_iter.next().expect("no next error").to_string();
            let umi_count = parse_u32(elem_iter);
            umis.push((umi, umi_count));
        }
    }
    Some(SCEQClass::new(tvec, cells, cell_count, umis))
}

pub fn read_eq(file_reader: &mut BufReader<File>, buf: &mut String) -> Option<SCEQClass> {
    let _ = file_reader.read_line(buf);
    let eqc: Option<SCEQClass>;
    {
        let mut elem_iter = buf.trim_right().split_whitespace();
        eqc = read_eqclass_from_iter(&mut elem_iter);
    }
    buf.clear();
    eqc
}

pub fn parse_bfh(p: &Path) -> Option<ProcessedExp> {
    let file = File::open(&p);
    let mut file_reader = BufReader::new(file.expect("file unwrap error"));
    let mut buf = String::new();

    // read # of transcripts
    let _ = file_reader.read_line(&mut buf);
    buf.pop();
    let num_txp = buf.parse::<u64>().expect("parse 64 error numtxps");
    buf.clear();

    // read # of barcodes
    let _ = file_reader.read_line(&mut buf);
    buf.pop();
    let num_bc = buf.parse::<u64>().expect("parse 64 error numbcs");
    buf.clear();

    // read # of eq classes
    let _ = file_reader.read_line(&mut buf);
    buf.pop();
    let num_eq = buf.parse::<u64>().expect("parse 64 error numeqs");
    buf.clear();

    let mut txp_vec: Vec<String> = Vec::with_capacity(num_txp as usize);
    // read txps
    for _ in 1..=num_txp {
        let _ = file_reader.read_line(&mut buf);
        txp_vec.push(buf.trim_right().to_string());
        buf.clear();
    }

    let mut bc_vec: Vec<String> = Vec::with_capacity(num_bc as usize);
    // read barcodes
    for _ in 1..=num_bc {
        let _ = file_reader.read_line(&mut buf);
        bc_vec.push(buf.trim_right().to_string());
        buf.clear();
    }
    info!("txp vec size = {:?}", txp_vec.len());
    info!("bc vec size = {:?}", bc_vec.len());

    // now read each equivalence class
    let mut eqclasses: Vec<SCEQClass> = Vec::with_capacity(num_eq as usize);
    for n in 1..=num_eq {
        let eqc = read_eq(&mut file_reader, &mut buf);
        eqclasses.push(eqc.expect("can't push to vector"));
        if n % 10000 == 0 {
            print!("\rparsed {:?} eq classes", n);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
    }
    println!();
    let proc_exp = ProcessedExp {
        transcripts: txp_vec,
        num_cells: num_bc,
        barcodes: bc_vec,
        eq_classes: eqclasses,
    };
    Some(proc_exp)
}

pub fn extract_cell(fexp: &ProcessedExp, cn: u32) -> Option<(String, CellExp)> {
    let mut ceq: Vec<CellEQClass> = Vec::new();
    for eqc in fexp.eq_classes.iter() {
        let p = eqc.cells.iter().position(|&x| x == cn);
        match p {
            Some(i) => {
                let start = eqc.cell_counts[i] as usize;
                let stop = eqc.cell_counts[i + 1] as usize;
                let mut v: Vec<(String, u32)> = Vec::with_capacity(stop - start);
                v.extend(eqc.umis[start..stop].iter().cloned());
                ceq.push(CellEQClass {
                    transcripts: &eqc.transcripts,
                    umis: v,
                });
            }
            None => {}
        }
    }
    Some((
        fexp.barcodes.get(cn as usize).expect("can't get barcodes").to_string(),
        CellExp {
            transcripts: &fexp.transcripts,
            eq_classes: ceq,
        },
    ))
}

pub fn graph_from_cell(
    ce: &CellExp,
    tgmap: &HashMap<String, String, BuildHasherDefault<FNV1aHasher64>>,
    verbose: bool,
) -> (
    petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    BitVec<LittleEndian, u64>,
    f64,
    f64,
) {
    // build the txp -> eqclass vector map
    fn get_map() -> HashMap<u32, Vec<u32>, BuildHasherDefault<FNV1aHasher64>> {
        HashMap::default()
    }
    // build the txp -> eqclass vector map
    fn get_set() -> HashSet<u32, BuildHasherDefault<FNV1aHasher64>> {
        HashSet::default()
    }
    fn get_string_set() -> HashSet<String, BuildHasherDefault<FNV1aHasher64>> {
        HashSet::default()
    }

    // Get a map from each transcript to it's list of eq classes
    let mut tid_map = get_map();
    let mut multi_gene = 0u64;
    let mut multi_gene_vec: BitVec<LittleEndian, u64> = BitVec::with_capacity(ce.eq_classes.len());
    for (eqid, eq) in ce.eq_classes.iter().enumerate() {
        let mut gset = get_string_set();
        for t in eq.transcripts.iter() {
            tid_map.entry(*t).or_insert(Vec::new()).push(eqid as u32);
            //gset.insert(tgmap.get(&ce.transcripts[*t as usize]).expect("can't extract tname").to_string());
            let gene = tgmap.get(&ce.transcripts[*t as usize]);
            match gene {
                Some(gname) => gset.insert(gname.to_string()),
                None => {
                    println!("{:?}, {:?}", *t as usize, ce.transcripts[*t as usize]);
                    panic!("can't extract gene name from eqclass")
                },
            };
        }
        if gset.len() > 1 {
            multi_gene += 1;
            multi_gene_vec.push(true);
        } else {
            multi_gene_vec.push(false);
        }
    }

    if verbose {
        info!(
            "Of {:?} eq classes {:?} are mulit-gene",
            ce.eq_classes.len(),
            multi_gene
        );
    }

    enum EdgeType {
        NoEdge,
        BiDirected,
        XToY,
        YToX,
    }

    let has_edge = |x: &(String, u32), y: &(String, u32)| -> EdgeType {
        if &x.0 == &y.0 {
            return EdgeType::BiDirected;
        }
        if x.1 > (2 * y.1 - 1) {
            if levenshtein(&x.0, &y.0) < 2 {
                return EdgeType::XToY;
            } else {
                return EdgeType::NoEdge;
            }
        } else if y.1 > (2 * x.1 - 1) {
            if levenshtein(&x.0, &y.0) < 2 {
                return EdgeType::YToX;
            } else {
                return EdgeType::NoEdge;
            }
        }
        EdgeType::NoEdge
    };

    let mut bidirected = 0u64;
    let mut unidirected = 0u64;

    let mut bidirected_in_multigene = 0u64;
    let mut unidirected_in_multigene = 0u64;

    let mut graph = DiGraphMap::<(u32, u32), ()>::new();
    for eqid in 0..ce.eq_classes.len() {
        if verbose && eqid % 1000 == 0 {
            print!("\rprocessed {:?} eq classes", eqid);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        // for each equivalence class
        let eq = &ce.eq_classes[eqid as usize];
        let u1 = &eq.umis;
        for (xi, x) in u1.iter().enumerate() {
            graph.add_node((eqid as u32, xi as u32));
            for xi2 in (xi + 1)..u1.len() {
                let x2 = &u1[xi2];
                graph.add_node((eqid as u32, xi2 as u32));
                let et = has_edge(&x, &x2);
                match et {
                    EdgeType::BiDirected => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        bidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            bidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::XToY => {
                        graph.add_edge((eqid as u32, xi as u32), (eqid as u32, xi2 as u32), ());
                        unidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            unidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::YToX => {
                        graph.add_edge((eqid as u32, xi2 as u32), (eqid as u32, xi as u32), ());
                        unidirected += 1;
                        if multi_gene_vec[eqid] == true {
                            unidirected_in_multigene += 1;
                        }
                    }
                    EdgeType::NoEdge => {}
                }
            }
        }

        let mut hset = get_set();
        // for every transcript
        for t in eq.transcripts.iter() {
            // find the equivalence classes sharing this transcript
            for eq2id in tid_map[t].iter() {
                if (*eq2id as usize) <= eqid {
                    continue;
                }
                if hset.contains(eq2id) {
                    continue;
                }
                hset.insert(*eq2id);
                let eq2 = &ce.eq_classes[*eq2id as usize];
                // compare all the umis
                let u2 = &eq2.umis;
                for (xi, x) in u1.iter().enumerate() {
                    // Node for equiv : eqid and umi : xi
                    graph.add_node((eqid as u32, xi as u32));
                    for (yi, y) in u2.iter().enumerate() {
                        // Node for equiv : eq2id and umi : yi
                        graph.add_node((*eq2id as u32, yi as u32));
                        let et = has_edge(&x, &y);
                        match et {
                            EdgeType::BiDirected => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                bidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    bidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::XToY => {
                                graph.add_edge((eqid as u32, xi as u32), (*eq2id, yi as u32), ());
                                unidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    unidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::YToX => {
                                graph.add_edge((*eq2id, yi as u32), (eqid as u32, xi as u32), ());
                                unidirected += 1;
                                if multi_gene_vec[eqid] == true
                                    || multi_gene_vec[*eq2id as usize] == true
                                {
                                    unidirected_in_multigene += 1;
                                }
                            }
                            EdgeType::NoEdge => {}
                        }
                    }
                }
            }
        }
    }
    if verbose {
        info!("tid_map of size {:?}", tid_map.len());
        info!(
            "size of graph ({:?}, {:?})",
            graph.node_count(),
            graph.edge_count()
        );
    }

    let ambig_bi: f64 = if bidirected > 0 {
        100.0f64 * ((bidirected_in_multigene as f64) / bidirected as f64)
    } else {
        0.0f64
    };
    let ambig_uni: f64 = if unidirected > 0 {
        100.0f64 * ((unidirected_in_multigene as f64) / unidirected as f64)
    } else {
        0.0f64
    };

    if verbose {
        info!(
            "bi-directed edges: {:?}, uni-directed edges : {:?}",
            bidirected, unidirected
        );
        info!("bi-directed in multigene : {:?}%", ambig_bi);
        info!("uni-directed in multigene : {:?}%", ambig_uni);
    }

    (graph, multi_gene_vec, ambig_bi, ambig_uni)
}

pub fn parse_tgmap(p: &Path) -> Option<HashMap<String, String, BuildHasherDefault<FNV1aHasher64>>> {
    // build the txp -> eqclass vector map
    fn get_map() -> HashMap<String, String, BuildHasherDefault<FNV1aHasher64>> {
        HashMap::default()
    }
    let mut m = get_map();
    let file = File::open(&p);
    let file_reader = BufReader::new(file.expect("can't read file"));

    for line in file_reader.lines() {
        let uline = line.expect("can't extract line");
        let mut l = uline.trim_right().split_whitespace();
        m.insert(l.next().expect("can't extract line").to_string(), l.next().expect("can't extract next").to_string());
    }
    return Some(m);
}
