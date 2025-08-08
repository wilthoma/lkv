
use std::fs::File;
use std::io::BufReader;
use crate::cyc_and_lie::*;
use crate::densematrix::*;

use rustc_hash::FxHashMap;
use rayon::prelude::*;
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle};
use std::str::FromStr;


fn read_pairs_from_json_str(path: &str) -> std::io::Result<Vec<(CycInt, String)>> {
    // Open the file
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    // Deserialize the JSON into Vec<(CycInt, i64)>
    let data: Vec<(CycInt, String)> = serde_json::from_reader(reader)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;

    Ok(data)
}

pub fn json_to_vector_str(path: &str, word_length: usize) -> CycVector {
    // word length is the lie word length... the cyclic word length and hence output word length is one more
    let pairs = read_pairs_from_json_str(path,).unwrap();
    let mut vec = CycVector::new(word_length+1);
    for (key, value) in pairs {
        if let Err(e) = MyInt::from_str(&value) {
            println!("Failed to parse value '{}' as MyInt: {}", value, e);
        }
        vec.add_to_entry(key, MyInt::from_str(&value).unwrap());
    }
    vec
}



pub fn get_brown_generators(max_weight : usize, path : &str) -> FxHashMap<usize, Vec<CycVector>> {
    println!("Loading Brown generators from JSON files.");
    let start = Instant::now();
    let mut generators = FxHashMap::default();
    
    for i in (12..=max_weight).step_by(2) {
        let mut vecs = Vec::new();
        let numvecs = i / 12 - (if i % 12 == 2 {1} else {0});
        for j in 1..=numvecs {
            let file_path = std::path::Path::new(path)
                .join(format!("g{}_{}.json", i, j))
                .to_string_lossy()
                .to_string();
            if let Ok(_) = std::fs::metadata(&file_path) {
                vecs.push(json_to_vector_str(&file_path, i));
            } else {
                vecs.push(json_to_vector_str(&file_path, i));
            }
        }
        if !vecs.is_empty() {
            generators.insert(i, vecs);
        }
    }

    let duration = start.elapsed();
    println!("Finished loading Brown generators in {:?}", duration);
    generators
}

pub fn sigma_generator(k : usize) -> CycVector {
    let mut v = MyVec::new(1);
    v.add_to_entry(1, MyInt::from(1));
    let mut vx = MyVec::new(1);
    vx.add_to_entry(0, MyInt::from(1));
    for _ in 0..2*k {
        v = vx.lie_bracket(&v);
    }
    let mut vec = CycVector::new(2*k+2);
    for (i, val) in v.into_iter() {
        if val != MyInt::from(0) {
            // add one "y" to make a cyclic word
            vec.add_to_entry(1 + (i << 1), val);
        }
    }
    vec
}


pub fn get_generator_list_2(max_word_length: usize, max_depth : usize, brown_path : &str) -> FxHashMap<(usize, usize), Vec<Vec<CycVector>>> {
    let brown_generators = get_brown_generators(max_word_length, brown_path);
    let mut generators: FxHashMap<(usize, usize), Vec<Vec<CycVector>>> = FxHashMap::default();

    for word_length in 2..=max_word_length {
        for depth in 1..=max_depth {
            println!("Computing generators for word length {} and depth {}", word_length, depth);
            let mut veclist: Vec<Vec<CycVector>> = Vec::new();
            // in depth 1, add sigma generators
            if depth == 1 && word_length % 2 == 1 && word_length >= 3 {
                veclist.push(vec![sigma_generator((word_length - 1) / 2)]);
            }
            // in depth 4, need to add brown generators
            if depth == 4 {
                if let Some(brown_vecs) = brown_generators.get(&word_length) {
                    veclist.extend(brown_vecs.iter().map(|vec| vec![vec.clone()]));
                }
            }
            // brackets with sigma generators
            if depth >= 1 {
                let max_k = (word_length - 2) / 2;
                for k in 1..=max_k {
                    if let Some(vecs) = generators.get(&(word_length-2*k-1, depth - 1)) {
                        // we can only add sigma generators if the previous depth has been computed
                        if vecs.len() > 0 {
                            let sg = sigma_generator(k);

                            let new_vecs: Vec<_> = vecs.par_iter()
                                .map({
                                    move |vec| {
                                        let mut new_vec = vec.clone();
                                        new_vec.push(sg.clone());
                                        new_vec
                                    }
                                })
                                .collect();
                            veclist.extend(new_vecs);
                        }
                    }
                }
            }
            // brackets with brown generators
            if depth >= 5 {
                let max_k = (word_length - 2) / 2;
                for k in 1..=max_k {
                    if let Some(brown_vecs) = brown_generators.get(&(2*k)) {
                        if let Some(vecs) = generators.get(&(word_length-2*k, depth - 4)) {
                            let new_vecs: Vec<_> = brown_vecs.par_iter()
                                .flat_map_iter(|vec| {
                                    vecs.iter().map(move |vec2| 
                                        {
                                            let mut new_vec = vec2.clone();
                                            new_vec.push(vec.clone());
                                            new_vec
                                        }
                                    )
                                })
                                .collect();
                            veclist.extend(new_vecs);
                        }
                    }
                }
            }
            if !veclist.is_empty() {
                generators.insert((word_length, depth), veclist);
            }
        }
    }
    generators
}

pub fn iterated_half_lie_bracket(
    vecs: &[CycVector],
) -> CycVector {
    if vecs.is_empty() {
        return CycVector::new(0);
    }
    let mut result = vecs[0].clone();
    for vec in &vecs[1..] {
        result = result.half_lie_bracket(vec);
    }
    result
}

pub fn get_hat_lkv_dim_lowerbound2(word_length: usize, depth: usize, p: MyUInt, gen_lists : &FxHashMap<(usize, usize), Vec<Vec<CycVector>>>) -> usize {
    println!("Computing lower bound for hat lkv dimension for word length {} and depth {}", word_length, depth);
    if let Some(vecs) = gen_lists.get(&(word_length, depth)) {
        println!("Found {} generators for word length {} and depth {}", vecs.len(), word_length, depth);
        // println!("{}", vecs[0]);
        // make matrix of the vectors in vecs using a vector compressor
        let start = Instant::now();
        let cur_outwords = all_cyclic_words_ones2(word_length+1, depth+1, depth+1).into_iter().collect::<Vec<CycInt>>();
        // let out_words = all_cyclic_words(word_length+1).into_iter().collect::<Vec<CycInt>>();
        let duration = start.elapsed();
        println!("all_cyclic_words took {:?}", duration);

        let start = Instant::now();
        let target_size = sder2_dim_xy(word_length-depth, depth+1);
        if target_size == 0 {
            println!("No generators for word length {} and depth {}", word_length, depth);
            return 0;
        }
        let col_compressor = VectorCompressor::new_random(&cur_outwords, target_size, 2, p, word_length+1);
        let duration = start.elapsed();
        println!("VectorCompressor took {:?}", duration);
        let pb = ProgressBar::new(vecs.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"));

        let row_results: Vec<Vec<MyUInt>> = vecs.par_iter()
            .map(|v| {
            let vv = iterated_half_lie_bracket(v);
            let compressed = col_compressor.compress(&vv.data);
            pb.inc(1);
            compressed
            })
            .collect::<Vec<_>>();
        pb.finish_with_message("Done");

        let mut matrix = DenseMatrix{ 
            rows: vecs.len(),
            cols: target_size,
            data: row_results,
        };

        // compress the matrix
        // matrix.compress_rows(target_size, 5, p);

        let cur_rank = matrix.rank(p as MyUInt);
        cur_rank
    } else {
        0 // no generators for this word length and depth
    }
}



