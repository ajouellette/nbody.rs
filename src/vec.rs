pub type NdVec<const NDIM: usize> = [f64; NDIM];

pub fn norm2<const NDIM: usize>(vec: NdVec<NDIM>) -> f64 {
    vec.iter().copied().map(|x| x * x).sum()
}

pub fn norm<const NDIM: usize>(vec: NdVec<NDIM>) -> f64 {
    vec.iter().copied().map(|x| x * x).sum::<f64>().sqrt()
}
