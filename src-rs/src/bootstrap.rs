use std::fs::File;
use std::io::Write;

use rand::StdRng;
use statrs::distribution::{Categorical, Distribution};

use optimize::{optimize, SalmonEQClass};

pub fn write_bootstraps(
    names: Vec<String>,
    alphas: Vec<Vec<f32>>,
    unique: Vec<bool>,
    not_ambiguous: Vec<bool>,
    num_bootstraps: u32,
) {
    let mut buffer = File::create("bootstrap_quants.sf").expect("can't open file");

    assert!(alphas.len() == num_bootstraps as usize, "different number of bootstraps");

    write!(buffer,
        "gene\talevin\thas_unique"
    ).expect("can't write quants to file");
    for boot_id in 0..num_bootstraps {
        write!(
            buffer,
            "\t{}",
            boot_id,
        ).expect("can't write bootstraps to file");
    }
    write!(buffer, "\n").expect("can't write bootstraps to file");

    for index in 0..names.len(){
        write!(buffer, 
            "{}\t{}\t{}", 
            names[index],
            not_ambiguous[index],
            unique[index])
            .expect("can't write quants to file");
        //TODO: might have to change based on the matrix-major format
        for boot_id in 0..num_bootstraps {
            write!(
                buffer,
                "\t{}",
                alphas[boot_id as usize][index],
            ).expect("can't write bootstraps to file");   
        }
        write!(buffer, "\n").expect("can't write bootstraps to file");
    }
}

pub fn do_bootstrapping(
    eqclasses: Vec<SalmonEQClass>, 
    unique_evidence: &mut Vec<bool>,
    no_ambiguity: &mut Vec<bool>,
    num_bootstraps: &u32,
    num_alphas: usize,
    only_unique: bool
) -> Vec<Vec<f32>> {
    let mut boot_quants: Vec<Vec<f32>> = Vec::new();
    let num_eqclasses = eqclasses.len();

    //extract the probability mass of each eqclass
    let mut prob_mass: Vec<f64> = vec![0.0; num_eqclasses];
    for (index, eqclass) in eqclasses.iter().enumerate() {
        prob_mass[index] += eqclass.counts as f64;
    }

    // calculate the multivariate categorical/discrete distribution
    let mut rng = StdRng::new().unwrap();
    let eq_distribution = Categorical::new(&prob_mass).unwrap();

    let get_total_num_molecules = |eqclasses: &Vec<SalmonEQClass>| -> usize {
        let mut num_molecules: usize = 0;
        eqclasses.iter().for_each(|eqclass| {
            num_molecules += eqclass.counts as usize;
        });
        num_molecules
    };

    // initially total num of molecules
    let total_molecules = get_total_num_molecules(&eqclasses);

    // iterate multiple times for bootstraps
    // TODO: current implementation can be memory heavy
    // Alternatively we can dump each round of bootstrap
    // insteead storing all in memory
    for _ in 0..*num_bootstraps {
        let mut boot_eqclasses = eqclasses.clone();

        // clear the counts of all eqclasses
        boot_eqclasses.iter_mut().for_each(|eqclass| {
            eqclass.counts = 0;
        });

        for _ in 0..total_molecules {
            // resample frim the multivariate distribution 
            // NOTE: should it ind_sample or sample ??
            let eq_id = eq_distribution.sample::<StdRng>(&mut rng);

            // increment the count for relevant eqclass
            boot_eqclasses[eq_id as usize].counts += 1;
        }

        assert!(total_molecules == get_total_num_molecules(&boot_eqclasses), 
            "The number of molecules not same between two rounds");


        //let mut unique_evidence: Vec<bool> = vec![false; num_alphas];
        let boot_alphas = optimize(boot_eqclasses, 
            unique_evidence,
            no_ambiguity,
            num_alphas.clone(),
            only_unique);

        boot_quants.push(boot_alphas);
    }

    return boot_quants;
}