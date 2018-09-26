extern crate bitvec;
extern crate hashers;
extern crate levenshtein;
extern crate pbr;
extern crate petgraph;
extern crate rayon;
extern crate pretty_env_logger;
extern crate statrs;
extern crate rand;
extern crate libc;

#[macro_use]
extern crate log;

mod analyze_cell;
mod optimize;
mod parse;
mod schema;
mod bootstrap;

use std::collections::HashMap;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc::channel;
use std::thread;

use libc::{c_char};
use std::ffi::CStr;
use hashers::fnv::FNV1aHasher64;
use pbr::ProgressBar;
use rayon::prelude::*;

#[no_mangle]
pub extern "C" fn do_parse_bfh(bfh_str: *const c_char,
                               tgmap_str: *const c_char,
                               out_dir_str: *const c_char,
                               num_threads: u64) {
    pretty_env_logger::init_timed();
    // TODO: add the below features later
    let num_bootstraps = 0u32;
    let nthread = num_threads as usize;
    let quant_only_cell_with_name = None;
    let quant_only_gene_with_name = None;
    let is_only_cell = !quant_only_cell_with_name.is_none();
    let only_unique = false;
    ///////////////////
    let tgmap_file_path = unsafe {
        assert!(!tgmap_str.is_null());
        CStr::from_ptr(tgmap_str).to_str()
    };

    let bfh_file_path = unsafe {
        assert!(!bfh_str.is_null());
        CStr::from_ptr(bfh_str).to_str()
    };
    let out_file_path = unsafe {
        assert!(!out_dir_str.is_null());
        CStr::from_ptr(out_dir_str).to_str()
    };

    println!("Using: {:?}", tgmap_file_path);
    println!("Using: {:?}", bfh_file_path);
    println!("Using: {:?}", out_file_path);

    let tm = Path::new(&tgmap_file_path.unwrap()).canonicalize();
    let tgmap = match tm {
        Ok(p) => {
            if p.is_file() {
                parse::parse_tgmap(&p)
            } else {
                None
            }
        }
        Err(_) => None,
    };
    let m = tgmap.unwrap();
    let iput = Path::new(&bfh_file_path.unwrap()).canonicalize();
    let rpath = match iput {
        Ok(p) => {
                parse::parse_bfh(&p)
        }
        Err(_) => None,
    };
    if rpath.is_none() {
        eprintln!("Could not access the input file");
        std::process::exit(1);
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(nthread as usize)
        .build_global()
        .unwrap();
    let x = rpath.unwrap();
    info!("num cells {:?}", x.num_cells);
    fn get_map() -> HashMap<String, u32, BuildHasherDefault<FNV1aHasher64>> {
        HashMap::default()
    }
    let mut gid_map = get_map();
    let mut gctr = 0u32;
    let mut gnames: Vec<String> = Vec::new();

    for (_k, v) in &m {
        if !gid_map.contains_key(v) {
            gid_map.insert(v.to_string(), gctr);
            gctr += 1;
            gnames.push(v.to_string());
        }
    }
    let (sender, receiver) = channel();
    let mut cell_id = 0u64;
    let quant_file = out_file_path.unwrap().to_owned();
    let mut buffer = File::create(quant_file+"/cell_quants.tsv").unwrap();
    let num_cells = match is_only_cell {
        true => 1,
        false => {
            write!(
                buffer,
                "cell_id\t{}\n",
                gnames
                    .iter()
                    .fold(String::new(), |acc, ref arg| acc + &arg.to_string() + "\t")
            ).unwrap();
            x.num_cells
        }
    };
    thread::spawn(move || {
        let mut pb = ProgressBar::new(num_cells);
        pb.format("╢▌▌░╟");
        for (bcode, countstr) in receiver.iter() {
            pb.inc();
            cell_id += 1;
            if !is_only_cell {
                write!(buffer,
                        "{}\t{}\n",
                        bcode,
                        countstr)
                    .expect("cannot write to output file!");
            }
        }
        pb.finish_print("done");
    });
    use std::sync::atomic::{AtomicBool, Ordering};
    let early_break = AtomicBool::new(false);

    (0..x.num_cells)
        .into_par_iter()
        .for_each_with(sender, |s, cell_id| {
            if !early_break.load(Ordering::Relaxed) {
                let (bcode, cell_exp) = parse::extract_cell(&x, cell_id as u32)
                    .expect("can't extract cell");
                if !is_only_cell
                    || Some(bcode.clone()) == quant_only_cell_with_name
                {
                    let (g, _mgv, _ambig_bi, _ambig_uni) =
                        parse::graph_from_cell(&cell_exp, &m, false);
                    let gcounts: Option<Vec<f32>> = analyze_cell::get_num_molecules(
                        &g,
                        &cell_exp,
                        &x,
                        &gid_map,
                        &m,
                        num_bootstraps,
                        quant_only_gene_with_name.clone(),
                        is_only_cell,
                        only_unique,
                    );
                    let ostr = match is_only_cell {
                        true => "".to_string(),
                        false => gcounts.unwrap()
                                    .iter()
                                    .fold(String::new(), |acc, &arg| acc + &arg.to_string() + "\t"),
                    };

                    s.send((bcode, ostr)).expect("can't send through mpsc");
                    if is_only_cell {
                        early_break.store(true, Ordering::Relaxed);
                    }
                }
            }
        });
    //Ok(())
}
