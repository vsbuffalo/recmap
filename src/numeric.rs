use ndarray::{Array1, Array2};
use num_traits::{cast::ToPrimitive, Float, NumCast};
use std::{
    cmp::Ordering,
    fmt::{Debug, Display},
};

/// Assert to float iterables are the same up to `eps`.
#[allow(dead_code)]
pub fn assert_float_eq<T>(left: T, right: T, eps: T)
where
    T: Float + Display,
{
    if left.is_nan() {
        assert!(right.is_nan(), "left is NaN, but right is not");
    } else {
        let diff = (left - right).abs();
        assert!(
            diff < eps,
            "values |{} - {}| ≥ {} (diff: {})",
            left,
            right,
            eps,
            diff
        );
    }
}

/// Assert to float values are the same up to `eps`.
#[allow(dead_code)]
pub fn assert_floats_eq<T>(left: &[T], right: &[T], eps: T)
where
    T: Float + Display,
{
    assert_eq!(left.len(), right.len());
    for (l, r) in left.iter().zip(right.iter()) {
        assert_float_eq(*l, *r, eps)
    }
}

#[allow(dead_code)]
pub fn absolute_differences<T>(a: &Array1<T>, b: &Array1<T>) -> Array1<T>
where
    T: Float,
{
    a.iter()
        .zip(b.iter())
        .map(|(&a_i, &b_i)| (a_i - b_i).abs())
        .collect::<Array1<T>>()
}

// Generic Euclidean distance function
#[allow(dead_code)]
pub fn euclidean_distance<T: Float>(x: T, y: T) -> T {
    (x - y).abs()
}

// Inverse Haldane's mapping function
//
// Brings map distances to recombination fractions
#[allow(dead_code)]
fn haldane_inverse<T: Float>(map_distance: T) -> T {
    T::from(0.5).unwrap() * (T::one() - (T::from(-2.0).unwrap() * map_distance).exp())
}

/// Build the pairwise recombination distance matrix, using some float-like type `T`.
///
/// Creates a `positions_x.len() x positions_y.len()` matrix of recombination
/// *distances* (in Morgans), for the supplied set of positions on the physical
/// map.
///
/// # Arguments
///  * `positions_x`: an [`ArrayView1`] of the first set of marker positions.
///  * `positions_y`: an [`ArrayView1`] of the second set of marker positions.
///  * `haldane`: whether to convert the recombination distances in *Morgans* to a
///      unit-less recombination *fraction*.
///  * `rec_floor`: an optional *floor* value; all elements in the matrix less than
///      this value will be set to this value. This is sometimes useful in downstream
///      processing when zero values create problems.
///
pub fn recomb_dist_matrix<T: Float>(
    positions_x: &[T],
    positions_y: &[T],
    haldane: bool,
    rec_floor: Option<T>,
) -> Array2<T> {
    let mut dist_matrix = Array2::<T>::zeros((positions_x.len(), positions_y.len()));
    for (i, &pos_x) in positions_x.iter().enumerate() {
        for (j, &pos_y) in positions_y.iter().enumerate() {
            let dist = euclidean_distance(pos_x, pos_y);
            let mut rf = if !haldane {
                dist
            } else {
                haldane_inverse(dist)
            };
            if let Some(min_rec_rate) = rec_floor {
                rf = rf.max(min_rec_rate);
            }
            dist_matrix[[i, j]] = rf;
        }
    }
    dist_matrix
}

#[derive(Debug, PartialEq)]
pub enum SearchResult {
    Exact(usize),
    LowerBound(usize),
    UpperBound(usize),
    LeftOf(usize),
}

impl SearchResult {
    #[allow(dead_code)]
    pub fn get_index(&self) -> usize {
        match self {
            SearchResult::Exact(idx) => *idx,
            SearchResult::LeftOf(idx) => *idx,
            SearchResult::LowerBound(idx) => *idx,
            SearchResult::UpperBound(idx) => *idx,
        }
    }
}

pub fn search_sorted<T: PartialOrd>(vec: &[T], new_val: T) -> SearchResult {
    let mut left = 0;
    let mut right = vec.len();
    while left < right {
        let mid = left + (right - left) / 2;

        match vec[mid].partial_cmp(&new_val).unwrap() {
            Ordering::Less => left = mid + 1,
            Ordering::Greater => right = mid,
            Ordering::Equal => return SearchResult::Exact(mid),
        }
    }

    if left == 0 {
        SearchResult::LowerBound(left)
    } else if left < vec.len() {
        SearchResult::LeftOf(left)
    } else {
        SearchResult::UpperBound(left)
    }
}

pub fn interp1d<Tx, Ty>(x: &[Tx], y: &[Ty], x0: Tx) -> Option<Ty>
where
    Tx: PartialOrd + ToPrimitive + Copy + Debug,
    Ty: ToPrimitive + NumCast + Copy + Debug,
{
    assert!(x.len() == y.len());
    let index = search_sorted(x, x0);
    match index {
        SearchResult::Exact(idx) => Some(y[idx]),
        SearchResult::LeftOf(idx) => {
            if idx == 0 || idx >= x.len() {
                return None;
            }

            let x1 = ToPrimitive::to_f64(&x[idx - 1])?;
            let x2 = ToPrimitive::to_f64(&x[idx])?;
            let y1 = ToPrimitive::to_f64(&y[idx - 1])?;
            let y2 = ToPrimitive::to_f64(&y[idx])?;
            let x0 = ToPrimitive::to_f64(&x0)?;

            // linear interpolation
            let y0 = y1 + (y2 - y1) * (x0 - x1) / (x2 - x1);

            NumCast::from(y0)
        }
        SearchResult::LowerBound(_) => Some(y[0]),
        SearchResult::UpperBound(idx) => Some(y[idx - 1]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_search_sorted_empty() {
        let vec: Vec<i32> = vec![];
        assert_eq!(search_sorted(&vec, 5), SearchResult::LowerBound(0));
    }

    #[test]
    fn test_search_sorted_exact_match() {
        let vec = vec![1, 2, 3, 4, 5];
        assert_eq!(search_sorted(&vec, 3), SearchResult::Exact(2));
    }

    #[test]
    fn test_search_sorted_no_exact_match_left_of() {
        let vec = vec![1, 3, 5, 7, 9];
        assert_eq!(search_sorted(&vec, 4), SearchResult::LeftOf(2));
    }

    #[test]
    fn test_search_sorted_no_exact_match_lower_bound() {
        let vec = vec![10, 20, 30, 40, 50];
        assert_eq!(search_sorted(&vec, 5), SearchResult::LowerBound(0));
    }

    #[test]
    fn test_search_sorted_no_exact_match_upper_bound() {
        let vec = vec![10, 20, 30, 40, 50];
        assert_eq!(search_sorted(&vec, 55), SearchResult::UpperBound(5));
    }

    #[test]
    fn test_search_sorted_with_floats() {
        let vec = vec![1.0, 2.5, 4., 4.8, 5.9];
        assert_eq!(search_sorted(&vec, 3.5), SearchResult::LeftOf(2));
    }
}
