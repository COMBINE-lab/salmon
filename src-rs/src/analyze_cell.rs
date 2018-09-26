use petgraph::prelude::*;
use petgraph::unionfind::*;
use petgraph::visit::NodeIndexable;
use hashers::fnv::FNV1aHasher64;
use petgraph;
//use rand;
//use rand::Rng;
use std;
use std::collections::{HashMap, HashSet, VecDeque};
use std::hash::BuildHasherDefault;

use bootstrap::{do_bootstrapping, write_bootstraps};
use optimize::{optimize, write_quants, SalmonEQClass};
use schema::{CellExp, ProcessedExp};

pub fn weakly_connected_components<G>(
    g: G,
) -> HashMap<u32, Vec<u32>, BuildHasherDefault<FNV1aHasher64>>
where
    G: petgraph::visit::NodeCompactIndexable + petgraph::visit::IntoEdgeReferences,
{
    let mut vertex_sets = UnionFind::new(g.node_bound());
    for edge in g.edge_references() {
        let (a, b) = (edge.source(), edge.target());

        // union the two vertices of the edge
        vertex_sets.union(g.to_index(a), g.to_index(b));
    }
    let labels = vertex_sets.into_labeling();
    fn get_map() -> HashMap<u32, Vec<u32>, BuildHasherDefault<FNV1aHasher64>> {
        HashMap::default()
    }

    let mut components = get_map();
    for (i, v) in labels.iter().enumerate() {
        let mut ve = components.entry(*v as u32).or_insert(Vec::new());
        ve.push(i as u32);
    }
    components
}

pub fn collapse_vertices(
    v: u32,
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    ce: &CellExp,
) -> (Vec<u32>, u32) {
    fn get_set(cap: u32) -> HashSet<u32, BuildHasherDefault<FNV1aHasher64>> {
        let mut h = HashSet::default();
        h.reserve(cap as usize);
        h
    }

    let mut largest_mcc: Vec<u32> = Vec::new();
    let mut chosen_txp = 0u32;
    let vert = g.from_index(v as usize);

    unsafe {
        for txp in ce
            .eq_classes
            .get_unchecked(vert.0 as usize)
            .transcripts
            .iter()
        {
            let mut bfs_list = VecDeque::new();
            bfs_list.push_back(v);

            let mut visited_set = get_set(16);
            visited_set.insert(v);

            let mut current_mcc = Vec::new();

            while let Some(cv) = bfs_list.pop_front() {
                current_mcc.push(cv);

                for nv in g.neighbors_directed(g.from_index(cv as usize), Outgoing) {
                    let n = g.to_index(nv) as u32;
                    if visited_set.contains(&n) {
                        continue;
                    } else {
                        visited_set.insert(n);
                    }
                    let n_labels = ce.eq_classes.get_unchecked(nv.0 as usize).transcripts;
                    if n_labels.contains(&txp) {
                        bfs_list.push_back(n);
                    }
                }
            }

            if largest_mcc.len() < current_mcc.len() {
                largest_mcc = current_mcc;
                chosen_txp = *txp;
            }
        }
    }

    (largest_mcc, chosen_txp)
}

pub fn get_num_molecules(
    g: &petgraph::graphmap::GraphMap<(u32, u32), (), petgraph::Directed>,
    ce: &CellExp,
    exp: &ProcessedExp,
    gid_map: &HashMap<String, u32, BuildHasherDefault<FNV1aHasher64>>,
    tgmap: &HashMap<String, String, BuildHasherDefault<FNV1aHasher64>>,
    num_bootstraps: u32,
    with_only_gene: Option<String>,
    is_only_cell: bool,
    only_unique: bool,
) -> Option<Vec<f32>> {
    fn get_set(cap: u32) -> HashSet<u32, BuildHasherDefault<FNV1aHasher64>> {
        let mut h = HashSet::default();
        h.reserve(cap as usize);
        h
    }

    let mut enable_gene_sampling = false;
    let mut subsample_gene_idx: u32 = <u32>::max_value();
    if !with_only_gene.clone().is_none() {
        enable_gene_sampling = true;
        subsample_gene_idx = *gid_map.get(&with_only_gene.expect("can't unwrap gene name"))
            .expect("can't found index of the given gene");

    }

    let comps = weakly_connected_components(g);
    let mut one_vertex_components: Vec<usize> = vec![0, 0];
    //let mut identified_txps: Vec<u32> = Vec::with_capacity(comps.len());

    // code for dumping pre collapse graph for sanity check.
    let dump_no_collapse = false;
    if dump_no_collapse {
        unsafe {
            for (_comp_label, comp_verts) in comps.iter() {
                for vin in comp_verts.iter() {
                    let vert = g.from_index(*vin as usize);
                    for txp in ce
                        .eq_classes
                        .get_unchecked(vert.0 as usize)
                        .transcripts
                        .iter()
                    {
                        eprint!("{},", txp);
                    }
                }
                eprintln!();
            }
        }
    }

    // Make salmon eqclasses for EM based counting
    let mut salmon_eqclass_hash: HashMap<Vec<u32>, u32, BuildHasherDefault<FNV1aHasher64>> =
        HashMap::default();

    for (_comp_label, comp_verts) in comps.iter() {
        if comp_verts.len() > 1 {
            let mut vset = get_set(comp_verts.len() as u32);
            for vin in comp_verts.iter() {
                vset.insert(*vin);
            }
            while !vset.is_empty() {
                let mut best_mcc: Vec<u32> = Vec::new();
                let mut best_covering_txp = std::u32::MAX;
                for v in vset.iter() {
                    let (new_mcc, covering_txp) = collapse_vertices(*v, g, ce);

                    if best_mcc.len() < new_mcc.len() {
                        best_mcc = new_mcc;
                        best_covering_txp = covering_txp;
                    }
                }

                // get gene_id of best covering transcript
                let tname = exp.transcripts.get(best_covering_txp as usize)
                    .expect("can't get transcript name");
                let gname = tgmap.get(tname)
                    .expect("can't get transcript to gene map");
                let best_covering_gene = gid_map.get(gname)
                    .expect("can't get transcript to gene map");

                if best_covering_txp == std::u32::MAX {
                    println!("Could not find a covering transcript");
                    std::process::exit(1);
                }

                let mut global_genes: HashSet<u32> = HashSet::new();
                unsafe {
                    for (index, vertex) in best_mcc.iter().enumerate() {
                        let vert = g.from_index((*vertex) as usize);
                        let mut local_genes = HashSet::<u32>::new();
                        for txp in ce
                            .eq_classes
                            .get_unchecked(vert.0 as usize)
                            .transcripts
                            .iter()
                        {
                            //eprint!("{},", txp);

                            let tname = exp.transcripts.get(*txp as usize)
                                .expect("can't get transcript name");
                            let gname = tgmap.get(tname)
                                .expect("can't get transcript to gene name");
                            let gene_idx = gid_map.get(gname)
                                .expect("can't get gene name");

                            // make a list of representative gene for this molecule
                            local_genes.insert(*gene_idx);
                        }

                        if index == 0 {
                            global_genes = local_genes;
                        }
                        else {
                            global_genes = global_genes
                                .intersection(&local_genes)
                                .cloned()
                                .collect();
                        }
                    }
                }
                //eprint!(";");

                // extract gene level salmon eqclass and increment count by 1
                assert!(
                    global_genes.len() > 0,
                    "can't find representative gene(s) for a molecule"
                );

                // assert the best covering gene in the global gene set
                assert!(global_genes.contains(best_covering_gene), 
                    "best gene not in covering set, shouldn't be possible");

                let mut eq_label: Vec<u32> = global_genes.iter().cloned().collect();
                eq_label.sort();
                
                let mut has_gene: bool = false;
                if enable_gene_sampling {
                    for gene_idx in global_genes {
                        if gene_idx == subsample_gene_idx {
                            has_gene = true;
                            break;
                        }
                    }
                }
                else{ has_gene = true; }

                if has_gene {
                    // incrementing the count of thie eqclass label by 1
                    let counter = salmon_eqclass_hash.entry(eq_label).or_insert(0);
                    *counter += 1;
                }

                //identified_txps.push(best_covering_txp);
                for rv in best_mcc.iter() {
                    vset.remove(rv);
                }
            } //end-while
        } else {
            let tv = comp_verts.first()
                .expect("can't extract first vertex");
            let tl = unsafe {
                ce.eq_classes
                    .get_unchecked(g.from_index(*tv as usize).0 as usize)
                    .transcripts
            };
            if tl.len() == 1 {
                one_vertex_components[0] += 1;
            } else {
                one_vertex_components[1] += 1;
            }

            let mut genes = HashSet::<u32>::new();

            for txp in tl {
                //eprint!("{},", txp);

                // extracting gene id
                // TODO: redundant code ********* START
                let tname = exp.transcripts.get(*txp as usize)
                    .expect("can't get transcript name");
                let gname = tgmap.get(tname)
                    .expect("can't get transcript to gene map");
                let gene_idx = gid_map.get(gname)
                    .expect("can't get gene name");

                // make a list of representative gene for this molecule
                genes.insert(*gene_idx);
            }

            // extract gene level salmon eqclass and increment count by 1
            assert!(
                genes.len() > 0,
                "can't find representative gene(s) for a molecule"
            );

            let mut eq_label: Vec<u32> = genes.iter().cloned().collect();
            eq_label.sort();

            let mut has_gene: bool = false;
            if enable_gene_sampling {
                for gene_idx in genes {
                    if gene_idx == subsample_gene_idx {
                        has_gene = true;
                        break;
                    }
                }
            }
            else{ has_gene = true; }
            
            if has_gene {
                // incrementing the count of thie eqclass label by 1
                let counter = salmon_eqclass_hash.entry(eq_label).or_insert(0);
                *counter += 1;
            }

            // TODO: redundant code ********* END
        }

        //let rand_cover = rand::thread_rng().choose(&tl)
        //    .expect("can;t get random cover");
        //identified_txps.push(*rand_cover as u32);
    }

    let mut salmon_eqclasses = Vec::<SalmonEQClass>::new();
    for (key, val) in salmon_eqclass_hash {
        salmon_eqclasses.push(SalmonEQClass {
            labels: key,
            counts: val,
        });
    }

    let mut unique_evidence: Vec<bool> = vec![false; gid_map.len()];
    let mut no_ambiguity: Vec<bool> = vec![true; gid_map.len()];
    if is_only_cell {
        info!("Total Networks: {}", comps.len());
        let num_txp_unique_networks: usize = one_vertex_components.iter().sum();
        let num_txp_ambiguous_networks: usize = comps.len() - num_txp_unique_networks;
        info!(
            ">1 vertices Network: {}, {}%",
            num_txp_ambiguous_networks,
            num_txp_ambiguous_networks as f32 * 100.0 / comps.len() as f32
        );
        info!(
            "1 vertex Networks w/ 1 txp: {}, {}%",
            one_vertex_components[0],
            one_vertex_components[0] as f32 * 100.0 / comps.len() as f32
        );
        info!(
            "1 vertex Networks w/ >1 txp: {}, {}%",
            one_vertex_components[1],
            one_vertex_components[1] as f32 * 100.0 / comps.len() as f32
        );

        //info!("Total Predicted Molecules {}", identified_txps.len());

        // iterate and extract gene names
        let mut gene_names: Vec<String> = vec!["".to_string(); gid_map.len()];
        for (gene_name, gene_idx) in gid_map {
            gene_names[*gene_idx as usize] = gene_name.clone();
        }

        if num_bootstraps > 0 {
            //entry point for bootstrapping
            let gene_counts: Vec<Vec<f32>> = do_bootstrapping(salmon_eqclasses, 
                &mut unique_evidence,
                &mut no_ambiguity,
                &num_bootstraps, 
                gid_map.len(),
                only_unique);
        
            write_bootstraps(gene_names, gene_counts, unique_evidence, 
                no_ambiguity, num_bootstraps);
            return None;
        }
        else{
            //entry point for EM
            //println!("{:?}", subsample_gene_idx);
            //println!("{:?}", &salmon_eqclasses);
            let gene_counts: Vec<f32> = optimize(salmon_eqclasses, &mut unique_evidence, 
                &mut no_ambiguity, gid_map.len(), only_unique);

            write_quants(gene_names, gene_counts, unique_evidence, no_ambiguity);
            return None;
        } // end-else
    }
    else {
        Some(optimize(salmon_eqclasses, 
            &mut unique_evidence,
            &mut no_ambiguity,
            gid_map.len(),
            only_unique
        ))
    }
    //identified_txps
}
