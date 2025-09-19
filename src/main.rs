mod cyc_and_lie;
mod brown_polys;
mod densematrix;
mod auxiliary;

use cyc_and_lie::*;
use rustc_hash::FxHashMap;
use std::time::Instant;
use std::cmp::min;
use clap::{Command, Arg};
use brown_polys::*;

fn main() {
    let matches = Command::new("lkv")
        .version("0.1.0")
        .author("Florian Naef, Thomas Willwacher")
        .about("Linearized KV dimension computation")   
        .arg(
            Arg::new("start_len")
                .help("Starting length")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                .value_name("START_LEN"),
        )
        .arg(
            Arg::new("stop_len")
                .help("Stopping length")
                .required(true)
                .value_parser(clap::value_parser!(usize))
                .value_name("STOP_LEN"),
        )
        .arg(
            Arg::new("num_threads")
                .short('t')
                .long("threads")
                .help("Size of shared threadpool to use (excluding the worker threads)")
                .value_parser(clap::value_parser!(usize))
                .value_name("NUM")
                .default_value("4"),
        )
        .arg(
            Arg::new("prime")
                .short('p')
                .long("prime")
                .help("The prime number to use for modular arithmetic")
                .value_parser(clap::value_parser!(u32))
                .value_name("PRIME")
                //.default_value("27644437")
                ,
        )
        .arg(
            Arg::new("start_depth")
                .short('d')
                .help("The lowest depth to compute")
                .value_parser(clap::value_parser!(usize))
                .value_name("START_DEPTH")
                .default_value("1")
                ,
        )
        .arg(
            Arg::new("stop_depth")
                .short('D')
                .help("The highest depth to compute")
                .value_parser(clap::value_parser!(usize))
                .value_name("STOP_DEPTH")
                .default_value("99")
                ,
        )
        .arg(
            Arg::new("lower")
                .long("lower")
                .help("Compute a lower dimension bound as the number of linearly independent brackets of generators.")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("brown_path")
                .long("brown")
                .help("Path to the Brown generators")
                .value_parser(clap::value_parser!(String))
                .value_name("BROWN_PATH")
                .default_value("brown_polys")
        )
        .get_matches();

    let num_threads = *matches.get_one::<usize>("num_threads").unwrap_or(&0);
    let start_len = *matches.get_one::<usize>("start_len").unwrap();
    let stop_len = *matches.get_one::<usize>("stop_len").unwrap();
    // let p = *matches.get_one::<u32>("prime").unwrap_or(&127);
    let p = *matches.get_one::<u32>("prime").unwrap_or(&3323);
    let start_depth = *matches.get_one::<usize>("start_depth").unwrap_or(&1);
    let stop_depth = *matches.get_one::<usize>("stop_depth").unwrap_or(&99);
    let lower = matches.get_flag("lower");
    let brown_path: String = matches
        .get_one::<String>("brown_path")
        .cloned()
        .unwrap_or_else(|| String::from("brown_polys"));

    if num_threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .expect("Failed to create thread pool");
    }

    let p = MyUInt::from(p);


    if lower {
        println!("Computing lower bound for hat lkv dimensions. Using prime {} for modular arithmetic.", p);
        let mut results: FxHashMap<(usize, usize), usize> = FxHashMap::default();
        let gen_lists = get_generator_list_2(stop_len, min(stop_depth, stop_len), &brown_path);
        // let gen_lists = get_generator_list(stop_len, min(stop_depth, stop_len));
        for l in start_len..=stop_len {
            for d in start_depth..=stop_depth {
                let start = Instant::now();
                let dim = get_hat_lkv_dim_lowerbound2(l, d, p, &gen_lists);
                let duration = start.elapsed();
                println!("Lower bound for hat lkv dimension for length {} and depth {}: {} in {:?}", l, d, dim, duration);
                results.insert((l, d), dim);
            }
        }

        // print all results
        println!("\n\n ***** Results Summary *****");
        let mut sorted_results: Vec<_> = results.iter().collect();
        sorted_results.sort_by_key(|&(l, d)| (*l, *d));
        for (&(l, d), &dim) in sorted_results {
            if dim>0 {
                println!("Length {} and depth {}: Lower bound hat lkv dimension {}", l, d, dim);
            }
        }
    } else {
        println!("Computing KV dimensions. Using prime {} for modular arithmetic.", p);

        let mut results : FxHashMap<usize, FxHashMap<usize,usize>> = FxHashMap::default();
        let depths_list = (start_depth..=stop_depth).collect::<Vec<usize>>();
        for l in start_len..=stop_len {
            let start = Instant::now();
            let (total, kds) = get_lkv_dim(l, &depths_list, p);
            let duration = start.elapsed();
            println!("lkv dimensions computed in {:?}", duration);
            // println!("** make_div_matrix for length {} took {:?}. Size: {}x{}", l, duration, matrix.n_rows, matrix.n_cols);  
            results.insert(l, kds.clone());

            for (depth, kv_dim) in {
                let mut sorted: Vec<_> = kds.iter().collect();
                sorted.sort_by_key(|&(k, _)| k);
                sorted
            } {
                println!("depth: {}, kvdim: {}", depth, kv_dim);
            }
            println!("Total KV dimension for length {}: {}", l, total);
        }
        // print all results
        println!("\n\n ***** Results Summary *****");
        let mut sorted_results: Vec<_> = results.iter().collect();
        sorted_results.sort_by_key(|&(l, _)| *l);
        for (l, kds) in sorted_results {
            println!("Length {}: ", l);
            let mut total = 0;
            let mut sorted_kds: Vec<_> = kds.iter().collect();
            sorted_kds.sort_by_key(|&(depth, _)| *depth);
            for (depth, kv_dim) in sorted_kds {
                println!("  Depth {}: KV dimension {}", depth, kv_dim);
                total += kv_dim;
            }
            println!("Total KV dimension for length {}: {}", l, total);
        }
    }
}
