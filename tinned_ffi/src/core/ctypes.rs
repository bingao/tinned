use num_complex::Complex64;
use num_rational::Rational64;
use safer_ffi::prelude::*;

#[derive_ReprC]
#[repr(C)]
pub struct CComplex64 {
    re: f64,
    im: f64,
}

impl CComplex64 {
    #[inline]
    pub fn new(re: f64, im: f64) -> Self {
        Self {
            re,
            im,
        }
    }

    #[inline]
    pub fn re(&self) -> f64 {
        self.re
    }

    #[inline]
    pub fn im(&self) -> f64 {
        self.im
    }

    #[inline]
    pub fn to_complex64(&self) -> Complex64 {
        Complex64::new(self.re, self.im)
    }
}

#[derive_ReprC]
#[repr(C)]
pub struct CRational64 {
    num: i64,
    den: i64,
}

impl CRational64 {
    #[inline]
    pub fn new(num: i64, den: i64) -> Self {
        Self {
            num,
            den,
        }
    }

    #[inline]
    pub fn num(&self) -> i64 {
        self.num
    }

    #[inline]
    pub fn den(&self) -> i64 {
        self.den
    }

    #[inline]
    pub fn to_rational64(&self) -> Option<Rational64> {
        if self.den == 0 {
            None
        } else {
            Some(Rational64::new(self.num, self.den))
        }
    }
}
