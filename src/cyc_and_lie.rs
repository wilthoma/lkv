

use rayon::prelude::*;
use rustc_hash::{FxHashSet, FxHashMap};
use core::fmt;
use rand::Rng;
use crate::auxiliary::{ mod_add, mod_mul};
use crate::densematrix::DenseMatrix;
use std::time::Instant;
use indicatif::{ProgressBar, ProgressStyle};
use std::ops::{Add, AddAssign, Mul, MulAssign};

use ::i256::*;//{i256, i512};
// use i256::{i256, i512};
// use i512::i512;


// pub type MyInt = i128;  // for the vector construction
pub type MyInt = I512;  // for the vector construction
// pub type MyInt = i256;  // for the vector construction
// pub type MyUInt = u16; // matrix data type
pub type MyUInt = u32; // matrix data type
pub type CycInt = u128; // integer type for storing the cyclic words
// pub type MyInt = i16;

#[derive(Clone)]
pub struct MyVec {
    pub data : FxHashMap<CycInt, MyInt>,
    pub word_length : usize,
}

impl MyVec {
    #[inline(always)]
    pub fn get(&self, key : CycInt) -> MyInt {
        *self.data.get(&key).unwrap_or(&MyInt::from(0))
    }
    #[inline(always)]
    pub fn add_to_entry(&mut self, key : CycInt, val : MyInt) {
        let entry = self.data.entry(key).or_insert(MyInt::from(0));
        *entry += val;
        if *entry == MyInt::from(0) {
            self.data.remove(&key);
        }
    }
    pub fn new(word_length:usize) -> Self {
        MyVec {
            data: FxHashMap::default(),
            word_length : word_length
        }
    }
    pub fn into_iter(&self) -> impl Iterator<Item = (CycInt, MyInt)> + '_ {
        self.data.iter().map(|(&k, &v)| (k, v))
    }
    pub fn print(&self) {
        for (&x,val) in self.data.iter() {
            println!("{}: {}", word_to_bitstr(x, self.word_length), val)
        }
    }

    pub fn ass_product(&self, other: &MyVec) -> MyVec {
        let mut ret = MyVec::new(self.word_length + other.word_length);
        for (&x1, &v1) in self.data.iter() {
            for (&x2, &v2) in other.data.iter() {
                let y = x2 | (x1 << other.word_length);
                ret.add_to_entry(y, v1 * v2 );
            }
        }
        ret
    }

    pub fn subtract_in(&mut self, other: &MyVec) {
        for (&x2, &v2) in other.data.iter() {
            self.add_to_entry(x2, -v2);
        }
    }
    pub fn add_in(&mut self, other: &MyVec) {
        for (&x2, &v2) in other.data.iter() {
            self.add_to_entry(x2, v2);
        }
    }
    pub fn addmul_in(&mut self, other: &MyVec, factor : MyInt) {
        for (&x2, &v2) in other.data.iter() {
            self.add_to_entry(x2, v2 * factor);
        }
    }


    pub fn lie_bracket(&self, other: &MyVec) -> MyVec {
        let mut p1 = self.ass_product(other);
        let p2 = other.ass_product(self);
        p1.subtract_in(&p2);
        p1
    }

    pub fn lie_bracket_odd(&self, other: &MyVec) -> MyVec {
        let mut p1 = self.ass_product(other);
        let p2 = other.ass_product(self);
        if (self.word_length %2 !=0) && (other.word_length %2 !=0) {
            p1.add_in(&p2);
        } else {
            p1.subtract_in(&p2);
        }
        p1
    }
}


#[derive(Clone)]
pub struct CycVector {
    pub data : FxHashMap<CycInt, MyInt>,
    word_length : usize,
}

impl CycVector {
    pub fn new(word_length:usize) -> Self {
        CycVector {
            data: FxHashMap::default(),
            word_length : word_length
        }
    }
    #[inline(always)]
    pub fn get(&self, key : CycInt) -> MyInt {
        let key = canon_word(key, self.word_length);
        *self.data.get(&key).unwrap_or(&MyInt::from(0))
    }
    #[inline(always)]
    pub fn add_to_entry(&mut self, key : CycInt, val : MyInt) {
        let key = canon_word(key, self.word_length);
        let entry = self.data.entry(key).or_insert(MyInt::from(0));
        *entry += val;
        if *entry == MyInt::from(0) {
            self.data.remove(&key);
        }
    }

    pub fn from_myvec(myvec: MyVec, pad_to_length : usize) -> Self {
        let mut ret = CycVector::new(pad_to_length);
        for (x, v) in myvec.into_iter() {
            ret.add_to_entry(x, v);
        }
        ret
    }

    // the output cyclic words have word length length
    // here x must be a word of length length-1
    pub fn from_word(x : CycInt, length : usize) -> Self{
        let v = lie_to_ass(x, length-1);
        CycVector::from_myvec(v, length)
    }

    pub fn from_lyndon_word(x : CycInt, length : usize, lyndon_cache : &FxHashMap<usize, FxHashSet<CycInt>>) -> Self{
        let v = lyndon_to_ass(x, length-1, lyndon_cache);
        CycVector::from_myvec(v, length)
    }

    pub fn from_words(x1 : CycInt, length1 : usize, x2 : CycInt, length2 : usize) -> CycVector {
        let v1 = lie_to_ass(x1, length1);
        let v2 = lie_to_ass(x2, length2);
        let mut ret = CycVector::new(length1+length2);
        for (y1, val1) in v1.into_iter() {
            for (y2,val2) in v2.into_iter() {
                let y = y1 | (y2 << length1);
                ret.add_to_entry(y, val1 * val2);
            }
        }
        ret
    }

    pub fn from_lyndon_words(x1 : CycInt, length1 : usize, x2 : CycInt, length2 : usize, lyndon_cache : &FxHashMap<usize, FxHashSet<CycInt>>) -> CycVector {
        let v1 = lyndon_to_ass(x1, length1, lyndon_cache);
        let v2 = lyndon_to_ass(x2, length2, lyndon_cache);
        let mut ret = CycVector::new(length1+length2);
        for (y1, val1) in v1.into_iter() {
            for (y2,val2) in v2.into_iter() {
                let y = y1 | (y2 << length1);
                ret.add_to_entry(y, val1 * val2);
            }
        }
        ret
    }

    pub fn laplace(&self) -> CycVector {
        let mut ret = CycVector::new(self.word_length-1);
        for (&x, &val) in self.data.iter() {
            let tmp = laplace_word(x, self.word_length);
            for (&y, &val2) in tmp.data.iter() {
                ret.add_to_entry(y, val2 * val);
            }
        }
        ret
    }

    pub fn half_laplace(&self) -> CycVector {
        let mut ret = CycVector::new(self.word_length-1);
        for (&x, &val) in self.data.iter() {
            let tmp = half_laplace_word(x, self.word_length);
            for (&y, &val2) in tmp.data.iter() {
                ret.add_to_entry(y, val2 * val);
            }
        }
        ret
    }

    pub fn into_iter(&self) -> impl Iterator<Item = (CycInt, MyInt)> + '_ {
        self.data.iter().map(|(&k, &v)| (k, v))
    }

    pub fn lie_bracket(&self, other: &Self) -> Self {
        let mut ret = Self::new(self.word_length + other.word_length-1);
        for (x, val) in self.into_iter() {
            for (y, val2) in other.into_iter() {
                let mut xr = x;
                for _ in 0..self.word_length {
                    let mut yr = y;
                    for _ in 0..other.word_length {
                        if xr&0b1 == yr&0b1 {
                            let z = (xr>>1) << other.word_length | yr;
                            ret.add_to_entry(z, val * val2);
                            let z = xr << (other.word_length-1) | (yr>>1);
                            ret.add_to_entry(z, -val * val2);

                        }
                        yr = cyc_rotate_left(yr, other.word_length);
                    }
                    xr = cyc_rotate_left(xr, self.word_length);
                }
            }
        }
        ret
    }

    pub fn half_lie_bracket(&self, other: &Self) -> Self {
        let mut ret = Self::new(self.word_length + other.word_length-1);
        for (x, val) in self.into_iter() {
            for (y, val2) in other.into_iter() {
                let mut xr = x;
                for _ in 0..self.word_length {
                    let mut yr = y;
                    for _ in 0..other.word_length {
                        if xr&0b1 == yr&0b1 && xr&0b1 == 1 {
                            let z = (xr>>1) << other.word_length | yr;
                            ret.add_to_entry(z, val * val2);
                            let z = xr << (other.word_length-1) | (yr>>1);
                            ret.add_to_entry(z, -val * val2);

                        }
                        yr = cyc_rotate_left(yr, other.word_length);
                    }
                    xr = cyc_rotate_left(xr, self.word_length);
                }
            }
        }
        ret
    }

}


impl<'a, 'b> Add<&'b CycVector> for &'a CycVector {
    type Output = CycVector;

    fn add(self, other: &'b CycVector) -> CycVector {
        assert_eq!(self.word_length, other.word_length, "CycVectors must have the same word length to be added");
        let mut ret = CycVector::new(self.word_length);
        for (x, v) in self.into_iter() {
            ret.add_to_entry(x, v);
        }
        for (x, v) in other.into_iter() {
            ret.add_to_entry(x, v);
        }
        ret
    }
}

// CycVector + &CycVector
impl<'b> Add<&'b CycVector> for CycVector {
    type Output = CycVector;

    fn add(self, other: &'b CycVector) -> CycVector {
        &self + other
    }
}

// &CycVector + CycVector
impl<'a> Add<CycVector> for &'a CycVector {
    type Output = CycVector;

    fn add(self, other: CycVector) -> CycVector {
        self + &other
    }
}

// CycVector + CycVector
impl Add<CycVector> for CycVector {
    type Output = CycVector;

    fn add(self, other: CycVector) -> CycVector {
        &self + &other
    }
}

impl<'a> AddAssign<&'a CycVector> for CycVector {
    fn add_assign(&mut self, other: &'a CycVector) {
        assert_eq!(self.word_length, other.word_length, "CycVectors must have the same word length to be added");
        for (x, v) in other.into_iter() {
            self.add_to_entry(x, v);
        }
    }
}
impl<'a> Mul<MyInt> for &'a CycVector {
    type Output = CycVector;

    fn mul(self, scalar: MyInt) -> CycVector {
        let mut ret = CycVector::new(self.word_length);
        for (x, v) in self.into_iter() {
            ret.add_to_entry(x, v * scalar);
        }
        ret
    }
}


impl Mul<MyInt> for CycVector {
    type Output = CycVector;

    fn mul(self, scalar: MyInt) -> CycVector {
        &self * scalar
    }
}

impl MulAssign<MyInt> for CycVector {
    fn mul_assign(&mut self, scalar: MyInt) {
        if scalar == MyInt::from(0) {
            self.data.clear(); // if scalar is 0, clear the vector
            return;
        }
        for (_, v) in self.data.iter_mut() {
            *v *= scalar;
        }
    }
}


impl fmt::Display for CycVector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut iter = self.into_iter().peekable();
        while let Some((x, val)) = iter.next() {
            write!(f, "({}){}", val, word_to_bitstr(x, self.word_length))?;
            if iter.peek().is_some() {
                write!(f, "+")?;
            }
        }  
        Ok(())
    }
}
impl fmt::Display for MyVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut iter = self.into_iter().peekable();
        while let Some((x, val)) = iter.next() {
            write!(f, "({}){}", val, word_to_bitstr(x, self.word_length))?;
            if iter.peek().is_some() {
                write!(f, "+")?;
            }
        }  
        Ok(())
    }
}

pub struct VectorCompressor {
    target_length: usize,
    spread_pattern: FxHashMap<CycInt, Vec<(usize,MyUInt)>>,
    word_length: usize,
    p: MyUInt,
}
impl  VectorCompressor {
    pub fn new_random(basis : &[CycInt], target_length: usize, spread : usize, p: MyUInt, word_length: usize) -> Self {
        // create some random spread pattern
        assert!(target_length <= basis.len(), 
            "Target length must be less than or equal to the basis length ({} <= {})",
            target_length, basis.len());
        let mut spread_pattern = FxHashMap::default();
        let mut rng = rand::rng();
        for (i,&x) in basis.iter().enumerate() {
            let idx1 = i % target_length;
            let v1 = rng.random_range(1..p) as MyUInt; // random value between 1 and p-1
            let mut res = vec![(idx1, v1)];

            for _ in 0..spread {
                let idx2 = rng.random_range(0..target_length);
                let v2 = rng.random_range(1..p) as MyUInt; // random value between 1 and p-1
                res.push((idx2, v2));
            }
            spread_pattern.insert(x, res);
        }
        VectorCompressor {
            target_length,
            spread_pattern,
            word_length,
            p,
        }
    }
    pub fn compress(&self, vec: &FxHashMap<CycInt,MyInt>) -> Vec<MyUInt> {
        let mut ret = vec![0; self.target_length];
        for (&x, &v) in vec.iter() {
            // println!("Compressing word {}", word_to_bitstr(x, self.word_length));
            let v : MyUInt = normalize_mod(v, self.p);
            if let Some(spread) = self.spread_pattern.get(&x) {
                for &(idx, val) in spread.iter() {
                    ret[idx] = mod_add(ret[idx], mod_mul(v, val, self.p), self.p);
                }
            } else {
                panic!("No spread pattern found for word {}", word_to_bitstr(x, self.word_length));
            }
        }
        ret
    }
}

#[inline(always)]
pub fn laplace_word(x : CycInt, length : usize) -> CycVector {
    let mut ret = CycVector::new(length-1);
    let mut x = x;
    for k in 0..length {
        // if (x & 0b11) == 0b11  {
        if (x & 0b11) == 0b11 || (x & 0b11) == 0b00 {
            ret.add_to_entry(x >> 1, MyInt::from(1));
        }
        if k<length-1 {
            x = cyc_rotate_left(x, length);
        }
    }
    ret
}
#[inline(always)]
pub fn half_laplace_word(x : CycInt, length : usize) -> CycVector {
    let mut ret = CycVector::new(length-1);
    let mut x = x;
    for k in 0..length {
        if (x & 0b11) == 0b11  {
        // if (x & 0b11) == 0b11 || (x & 0b11) == 0b00 {
            ret.add_to_entry(x >> 1, MyInt::from(1));  
        }
        if k<length-1 {
            x = cyc_rotate_left(x, length);
        }
    }
    ret
}

pub fn lie_to_ass(x : CycInt, length: usize) -> MyVec {
    let mut ret = MyVec::new(length);
    if length <= 1 {
        ret.add_to_entry(x, MyInt::from(1));
    } else {
        let mask = (1 << (length-1)) -1;
        let msk = 1 << (length-1);
        let xx = msk & x;
        let pre = lie_to_ass(x & mask, length -1);
        for (y, v) in pre.into_iter() {
            ret.add_to_entry(xx | y, v);
            ret.add_to_entry((y<< 1) | (xx>>(length-1)), -v); 
        }
    }
    ret
}

pub fn lyndon_to_ass(x : CycInt, length : usize, lyndon_cache : &FxHashMap<usize, FxHashSet<CycInt>>) -> MyVec {
    if length == 1 {
        let mut ret = MyVec::new(1);
        ret.add_to_entry(x, MyInt::from(1));
        return ret;
    }
    // find standard factorization
    let mut ulen=0;
    let mut u = 0;
    let mut v = 0;
    for u_length in 1..length {
        let v_length = length - u_length;
        let mask = (1 << v_length)-1;
        v = mask & x;
        u = x >> v_length;
        if lyndon_cache[&u_length].contains(&u) && lyndon_cache[&v_length].contains(&v) {
            ulen = u_length;
            break;
        }
    }

    if ulen == 0 {
        //error
        panic!("Could not find lyndon uv");
    }

    let vec1 = lyndon_to_ass(u, ulen, lyndon_cache);
    let vec2 = lyndon_to_ass(v, length-ulen, lyndon_cache);
    vec1.lie_bracket(&vec2)
}


pub fn word_to_bitstr(x : CycInt, length : usize) -> String {
    // convert uint to bit string
    let mut s = String::with_capacity(length);
    for i in (0..length).rev() {
        if (x & (1 << i)) != 0 {
            s.push('1');
        } else {
            s.push('0');
        }
    }
    s
}


pub fn make_lyndonwords_cache_depth2(maxlen : usize, max_depth : usize) -> FxHashMap<usize, FxHashSet<CycInt>> {
    let mut ret = FxHashMap::default();
    for l in 1..=maxlen {
        ret.insert(l, all_lyndon_words_depth_ex(l, 0, max_depth));
    }
    ret
}



// a function that rotates the last N bits of a given number by one
pub fn cyc_rotate_left(x : CycInt, length : usize) -> CycInt {
    let mask = (1 << (length-1)) - 1; // Create a mask for the last N bits
    let msk2 = 1<< (length-1); 
    let rotated_bits = (x & mask) << 1; // Shift the last N bits left by one
    let overflow_bits = (x & msk2) >> (length - 1); // Get the bit that overflows
    rotated_bits | overflow_bits // Combine the unchanged bits with the rotated bits and overflow bit
}

pub fn canon_word(x : CycInt, length: usize) -> CycInt {
    let mut min = x;
    let mut rotated = x;
    for _ in 1..length {
        rotated = cyc_rotate_left(rotated, length);
        if rotated < min {
            min = rotated;
        }
    }
    min
}


pub fn canon_word_lyndon(x : CycInt, length: usize) -> Option<CycInt> {
    let mut min = x;
    let mut rotated = x;
    let mut lyndonflag = true;
    for _ in 1..length {
        rotated = cyc_rotate_left(rotated, length);
        if rotated < min {
            min = rotated;
            lyndonflag = true; // if we find a rotation that is smaller than the minimum, it is (possibly) a Lyndon word
        } else if rotated == min {
            lyndonflag = false; // if we find a rotation that is equal to the minimum, it is not a Lyndon word
        }
    }
    if lyndonflag {
        return Some(min);
    }
    None // if no Lyndon word was found, return None
}


pub fn all_cyclic_words_ones2(length: usize, min_ones_count: usize, max_ones_count: usize) -> FxHashSet<CycInt> {
    // println!("Generating cyclic words with length {} and ones count from {} to {}", length, min_ones_count, max_ones_count);
    let max_ones_count = max_ones_count.min(length);
    let mut result = FxHashSet::default();
    for ones in min_ones_count..=max_ones_count {
        if ones == 0 {
            result.insert(0);
            continue; // skip zero ones count
        }
        // Generate all combinations of `ones` positions out of `length`
        let mut indices = (0..ones).collect::<Vec<_>>();
        while indices[0] < length - ones + 1 {
            // Build the word with ones at the selected indices
            let mut word = 0;
            for &i in &indices {
                word |= 1 << i;
            }
            result.insert(canon_word(word, length));
            // Generate next combination
            let mut t = ones;
            while t > 0 && indices[t - 1] == length - ones + t - 1 {
                t -= 1;
            }
            if t == 0 {
                break;
            }
            indices[t - 1] += 1;
            for j in t..ones {
                indices[j] = indices[j - 1] + 1;
            }
        }
    }
    result
}


pub fn all_lyndon_words_depth_ex(length: usize, mindepth : usize, maxdepth : usize) -> FxHashSet<CycInt> {
    // println!("Generating Lyndon words with length {} and depth from {} to {}", length, mindepth, maxdepth);
    let mut result = FxHashSet::default();
    let maxdepth = maxdepth.min(length);
    for ones in mindepth..=maxdepth {
        if ones == 0 {
            result.insert(0);
            continue; // skip zero ones count
        }
        // Generate all combinations of `ones` positions out of `length`
        let mut indices = (0..ones).collect::<Vec<_>>();
        while indices[0] < length - ones + 1 {
            // Build the word with ones at the selected indices
            let mut word = 0;
            for &i in &indices {
                word |= 1 << i;
            }
            if let Some(canon_word) = canon_word_lyndon(word, length) {
                // Only insert if the word is a Lyndon word
                result.insert(canon_word);
            }
            // Generate next combination
            let mut t = ones;
            while t > 0 && indices[t - 1] == length - ones + t - 1 {
                t -= 1;
            }
            if t == 0 {
                break;
            }
            indices[t - 1] += 1;
            for j in t..ones {
                indices[j] = indices[j - 1] + 1;
            }
        }
    }
    result
}





pub fn nr_liewords_xy(nx:usize, ny:usize) -> usize {
    let n = nx + ny;
    let g = gcd(nx, ny);
    (divisors(g).iter().map(|&d| {
        (i128::from(moebiusmu(d))) * (binomial_coefficient((n/d) as i128, (nx/d) as i128))
    }).sum::<i128>() as usize )/n
}

pub fn gcd(a: usize, b: usize) -> usize {
    // returns the greatest common divisor of a and b
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let mut a = a;
    let mut b = b;
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

pub fn divisors(n: usize) -> Vec<usize> {
    // returns the divisors of n
    let mut divs = Vec::new();
    for i in 1..=((n as f64).sqrt() as usize) {
        if n % i == 0 {
            divs.push(i);
            if i != n / i {
                divs.push(n / i);
            }
        }
    }
    divs.sort();
    divs
} 

pub fn moebiusmu(n:usize) -> i32 {
    let precomp = vec![0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 0, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, 0];
    
    if n > precomp.len() {
        panic!("Not precomputed for length {}", n);
    }
    return precomp[n];
}
pub fn binomial_coefficient(n:i128, m:i128) ->i128 {
    // returns the binomial coefficient n choose m
    let nn = i256::from(n);

    if m > n {
        return 0;
    }
    if m == 0 || m == n {
        return 1;
    }
    let mut res=i256::from(1); // = 1i128;
    for i in 0..m {
        let ii = i256::from(i);
        res *= nn - ii;
        res /= ii + i256::from(1);
    }
    if res > i256::from(i128::MAX) {
        panic!("Binomial coefficient overflow for n={}, m={}", n, m);
    }
    res.as_i128()
}


// warning: nx,ny are the valences of vertices 1 and 2. The degree of the components is nx+ny-1(!)
pub fn sder2_dim_xy(nx: usize, ny: usize) -> usize {
    // the dimension of the second derived algebra of the cyclic words of length nx+ny
    // is given by the number of cyclic words of length nx+ny-1
    if nx == 0 || ny == 0 {
        return 0; // no words of this type
    }
    nr_liewords_xy(nx-1, ny)+nr_liewords_xy(nx, ny-1) - nr_liewords_xy(nx, ny)
}


#[inline(always)]
pub fn normalize_mod(x: MyInt, p: MyUInt) -> MyUInt {
    // Normalize x to be in the range [0, mod_val)
    let pp = MyInt::from(p);
    // !!! THIS LINE NEEDS TO BE CHANGED IF USING OTHER MyInt or MyUInt TYPE !!!
    (((x % pp) + pp) % pp).as_u32()
    // (((x % pp) + pp) % pp) as MyUInt
}



pub fn get_halfdiv_rank_depth(word_length: usize, depths_list : &[usize], p: MyUInt) -> (usize, FxHashMap<usize, usize>) {
    let mut rank =0;
    let min_depth = depths_list.iter().min().cloned().unwrap_or(1);
    let max_depth = depths_list.iter().max().cloned().unwrap_or(word_length+1);
    // filter depths_list to only contain values in the range [0, word_length+
    println!("computing out words...");
    let out_words = all_cyclic_words_ones2(word_length, min_depth, max_depth).into_iter().collect::<Vec<CycInt>>();
    println!("counting ones...");
    let out_onebits = out_words.iter().map(|&w| w.count_ones() as usize).collect::<Vec<usize>>();
    println!("lyndon cache ...");
    let start = Instant::now();
    let lyndon_cache = make_lyndonwords_cache_depth2(word_length, max_depth);
    let duration = start.elapsed();
    println!("** make_lyndonwords_cache for length {} took {:?}", word_length, duration);  

    // the cyclic words have one longer length
    let length1 = (word_length+1)/2;
    let length2 = word_length + 1 - length1;

    let in_words1 = lyndon_cache[&length1].iter().cloned().collect::<Vec<CycInt>>();
    let in_words2 = lyndon_cache[&length2].iter().cloned().collect::<Vec<CycInt>>();
    let in_words = in_words1
        .iter()
        .flat_map(|&w1| in_words2.iter().map(move |&w2| (w1, w2)))
        .filter(|&(w1, w2)| length1 != length2 || w2 >= w1)
        .collect::<Vec<(CycInt, CycInt)>>();

    let in_onebits = in_words
        .iter()
        .map(|&(w1, w2)| w1.count_ones() as usize + w2.count_ones() as usize)
        .collect::<Vec<usize>>();
    let mut res_ranks : FxHashMap<usize, usize> = FxHashMap::default();

    for cur_depth in depths_list.iter().cloned() {
        let cur_ones = cur_depth + 1; // the number of ones is one higher than the depth
        if cur_ones > word_length+1 {
            continue; // skip if the number of ones is greater than the word length + 1
        }
        println!("+ processing {}-ones words", cur_ones);

        let target_size = sder2_dim_xy(word_length+1-cur_ones, cur_ones);// number of legs is one higher than word length 
        // println!("  sder dimension ({}-ones): {} ", cur_ones, target_size);
        if target_size == 0 {
            continue; // no words of this type
        }

        // create only the dense block matrix
        let cur_inwords = in_words.iter().zip(in_onebits.iter())
            .filter(|(_, oc)| **oc == cur_ones)
            .map(|(w,_)| *w).collect::<Vec<(CycInt, CycInt)>>();
        // println!("  found {} input words", cur_inwords.len());
        let cur_outwords = out_words.iter().zip(out_onebits.iter())
            .filter(|(_, oc)| **oc == cur_ones-1)
            .map(|(w,_)| *w).collect::<Vec<CycInt>>();
        // println!("  found {} output words", cur_outwords.len());
        let col_compressor = VectorCompressor::new_random(&cur_outwords, target_size, 5, p, word_length);

        let pb = ProgressBar::new(cur_inwords.len() as u64);
        pb.set_style(
            ProgressStyle::default_bar()
                .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("#>-"),
        );

        let row_results: Vec<Vec<MyUInt>> = cur_inwords
            .par_iter()
            .map_init(
                || pb.clone(),
                |pb, &(word1, word2)| {
                    let cycv = CycVector::from_lyndon_words(word1, length1, word2, length2, &lyndon_cache);
                    let lapv = cycv.half_laplace();
                    pb.inc(1);
                    col_compressor.compress(&lapv.data)
                },
            )
            .collect::<Vec<_>>();

        pb.finish_with_message("Done");
        let mut matrix = DenseMatrix{ 
            rows: cur_inwords.len(),
            cols: target_size,
            data: row_results,
        };

        // compress the matrix
        matrix.compress_rows(target_size, 5, p as MyUInt);
        // matrix.compress(target_size, target_size, 5, p);

        let cur_rank = matrix.rank(p as MyUInt);
        res_ranks.insert(cur_depth, cur_rank);
        rank += cur_rank;
        println!("  rank for {}-ones words: {} (corank {})", cur_ones, cur_rank, target_size - cur_rank);
    }

    (rank, res_ranks)
}




pub fn get_lkv_dim(word_length : usize, depths_list : &[usize], p: MyUInt) -> (usize, FxHashMap<usize, usize>) {
    // Calculate the rank of the division matrix for cyclic words of a given length
    let mut res_dims = FxHashMap::default();
    let mut dim = 0;
    let (_, res_ranks) = get_halfdiv_rank_depth(word_length, depths_list, p);
    for (depth, cur_rank) in res_ranks {
        let cur_ones = depth+1; // the number of ones is one higher than the depth
        let sderdim = sder2_dim_xy(word_length+1-cur_ones, cur_ones);
        let kvdim = sderdim - cur_rank;
        res_dims.insert(depth, kvdim);
        dim += kvdim;
        
        //println!("cur_ones: {}, kvdim: {}", cur_ones, kvdim);
        // println!("cur_ones: {}, cur_rank: {}, sderdim: {}, kvdim: {}", cur_ones, cur_rank, sderdim, kvdim);
    }
    (dim, res_dims)
}

