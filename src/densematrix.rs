use rayon::prelude::*;
use rand::Rng;
// use num_traits::{Zero, One};
use std::fmt::Debug;
use indicatif::{ProgressBar, ProgressStyle};
use crate::auxiliary::{
    mod_inv, mod_add, mod_mul, mod_sub, GoodInteger,
};

#[derive(Debug, Clone)]
pub struct DenseMatrix<T: GoodInteger> {
    pub data: Vec<Vec<T>>,
    pub rows: usize,
    pub cols: usize,
}

impl<T: GoodInteger> DenseMatrix<T> {
    // pub fn new(data: Vec<Vec<T>>) -> Self {
    //     let rows = data.len();
    //     let cols = if rows > 0 { data[0].len() } else { 0 };
    //     Self { data, rows, cols }
    // }

    #[inline(always)]
    pub fn add_row_multiple_to_row(
        &mut self,
        target_row: usize,
        source_row: usize,
        multiple: T,
        p: T,
    ) {
        if target_row < self.rows && source_row < self.rows {
            for c in 0..self.cols {
                self.data[target_row][c] = mod_add(
                    self.data[target_row][c],
                    mod_mul(self.data[source_row][c], multiple, p),
                    p,
                );
            }
        }
    }
    // #[inline(always)]
    // pub fn add_col_multiple_to_col(
    //     &mut self,
    //     target_col: usize,
    //     source_col: usize,
    //     multiple: T,
    //     p: T,
    // ) {
    //     if target_col < self.cols && source_col < self.cols {
    //         for r in 0..self.rows {
    //             self.data[r][target_col] = mod_add(
    //                 self.data[r][target_col],
    //                 mod_mul(self.data[r][source_col], multiple, p),
    //                 p,
    //             );
    //         }
    //     }
    // }
    pub fn compress_rows(&mut self, to_rows: usize, spread: usize, p: T) {
        if to_rows >= self.rows {
            return; // No compression needed
        }

        // for each row from to_rows to rows, add random multiples of the row to spread random rows in the range [0, to_rows)
        // then delete the row
        let mut rng = rand::rng();
        for r in to_rows..self.rows {
            for _ in 0..spread {
                let target_row = rng.random_range(0..to_rows);
                if target_row != r {
                    let multiple = rng.random_range(T::one()..p); // Random multiple
                    self.add_row_multiple_to_row(target_row, r, multiple, p);
                }
            }
        }
        self.data.truncate(to_rows);
        self.rows = to_rows;
    }
    // pub fn compress_cols(&mut self, to_cols: usize, spread: usize, p: T) {
    //     if to_cols >= self.cols {
    //         return; // No compression needed
    //     }

    //     // for each column from to_cols to cols, add random multiples of the column to spread random columns in the range [0, to_cols)
    //     // then delete the column
    //     use rand::Rng;
    //     let mut rng = rand::rng();
    //     for c in to_cols..self.cols {
    //         for _ in 0..spread {
    //             let target_col = rng.random_range(0..to_cols);
    //             if target_col != c {
    //                 let multiple = rng.random_range(T::one()..p); // Random multiple
    //                 self.add_col_multiple_to_col(target_col, c, multiple, p);
    //             }
    //         }
    //     }
    //     for r in 0..self.rows {
    //         self.data[r].truncate(to_cols);
    //     }
    //     self.cols = to_cols;
    // }
    // pub fn compress(&mut self, to_rows:usize, to_cols:usize, spread:usize,p: T) {
    //     let start = std::time::Instant::now();
    //     let old_rows = self.rows;
    //     let old_cols = self.cols;
    //     self.compress_rows(to_rows, spread, p);
    //     self.compress_cols(to_cols, spread, p);
    //     let elapsed = start.elapsed();
    //     println!(
    //         "Compressed matrix from {} to {} rows and from {} to {} cols in {:?}",
    //         old_rows, self.rows, old_cols, self.cols, elapsed
    //     );
    // }

    pub fn rank(&mut self, p: T) -> usize {
        let mut rank = 0;
        let mut row = 0;

        let total_steps = self.rows.min(self.cols) as u64;
        let bar = ProgressBar::new(total_steps);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} rows, ETA: {eta} {msg}")
                .unwrap()
                .progress_chars("#>-"),
        );

        while row < self.rows && rank < self.cols {
            let pivot_row = (row..self.rows).find(|&r| self.data[r][rank] != T::zero());

            if let Some(pivot_idx) = pivot_row {
                if pivot_idx != row {
                    self.data.swap(row, pivot_idx);
                }

                let pivot_val = self.data[row][rank];
                let pivot_inv = mod_inv(pivot_val, p);

                for c in rank..self.cols {
                    self.data[row][c] = mod_mul(self.data[row][c], pivot_inv, p);
                }

                let pivot_row_data = self.data[row].clone();
                let current_rows = &mut self.data;

                current_rows
                    .par_iter_mut()
                    .enumerate()
                    .filter(|(r, _)| *r != row)
                    .for_each(|(_, ref mut current_row)| {
                        let factor = current_row[rank];
                        if factor != T::zero() {
                            for c in rank..self.cols {
                                current_row[c] = mod_sub(
                                    current_row[c],
                                    mod_mul(factor, pivot_row_data[c], p),
                                    p,
                                );
                            }
                        }
                    });

                row += 1;
                bar.inc(1);
            }

            rank += 1;
        }

        bar.finish_with_message("Done");

        self.data
            .iter()
            .filter(|r| r.iter().any(|&x| x != T::zero()))
            .count()
    }

    // pub fn horiz_stack(&self, other: &Self) -> Self {
    //     if self.rows != other.rows {
    //         panic!("Cannot horizontally stack matrices with different number of rows");
    //     }
    //     if self.rows == 0 {
    //         return DenseMatrix::new(vec![vec![]; 0]);
    //     }
    //     let mut new_data = Vec::with_capacity(self.rows);
    //     for (row_self, row_other) in self.data.iter().zip(other.data.iter()) {
    //         let mut new_row = row_self.clone();
    //         new_row.extend(row_other.iter().cloned());
    //         new_data.push(new_row);
    //     }
    //     assert_eq!(new_data.len(), self.rows);
    //     assert_eq!(new_data[0].len(), self.cols + other.cols);
    //     DenseMatrix::new(new_data)
    // }
    // pub fn vert_stack(&self, other: &Self) -> Self {
    //     if self.cols != other.cols {
    //         panic!("Cannot vertically stack matrices with different number of columns");
    //     }
    //     let mut new_data = Vec::with_capacity(self.rows + other.rows);
    //     new_data.extend(self.data.clone());
    //     new_data.extend(other.data.clone());
    //     DenseMatrix::new(new_data)
    // }
}
