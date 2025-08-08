
use std::ops::{Add, Sub, AddAssign, SubAssign, Mul, Rem, MulAssign, Div, DivAssign};
use num_traits::{One, Zero};
use rand::distr::uniform::SampleUniform;
use std::fmt::Display;
use std::iter::Sum;

pub trait GoodInteger: MulAssign+ Div<Output = Self> + DivAssign + Sub<Output = Self> + SubAssign + Display+Sum+PartialOrd + SampleUniform + Copy + Zero + One + Rem<Output = Self> + Add<Output = Self> + AddAssign + Mul<Output = Self> + Send + Sync +'static { 

}


impl<T: MulAssign + Div<Output = Self> + DivAssign + Sub<Output = Self> + SubAssign + Display + Sum + PartialOrd + SampleUniform + Copy + Zero + One + Rem<Output = Self> + Add<Output = Self> + AddAssign + Mul<Output = Self> + Send + Sync + 'static> GoodInteger for T {

}

#[inline(always)]
pub fn mod_add<T : GoodInteger>(mut a: T, b: T, p:T) -> T {
    a += b;
    if a >= p { a -= p; }
    a
}

#[inline(always)]
pub fn mod_sub<T : GoodInteger>(mut a: T, b: T, p:T) -> T {
    if a < b { a += p; }
    a - b
}

#[inline(always)]
pub fn mod_mul<T : GoodInteger>(a: T, b: T, p : T) -> T 
{
    (a * b) % p
}

#[inline(always)]
pub fn mod_pow<T : GoodInteger>(mut base: T, mut exp: T, p: T) -> T {
    let mut result = T::one();
    let two = T::one() + T::one(); // Used for squaring
    while exp > T::zero() {
        if exp % two == T::one() {
            result = mod_mul(result, base, p);
        }
        base = mod_mul(base, base, p);
        exp /= two;
    }
    result
}

// modular inverse with mod_pow
#[inline(always)]
pub fn mod_inv<T: GoodInteger>(a: T, p: T) -> T {
    let two = T::one() + T::one();
    mod_pow(a, p - two, p)
}