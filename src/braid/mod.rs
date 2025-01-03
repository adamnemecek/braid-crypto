pub mod garside;
pub mod random;

// pub use crate::prelude::*;

use {
    bincode::{
        deserialize,
        serialize,
    },
    indexmap::set::IndexSet,
};

use crate::permutation::*;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BrGen {
    Sigma(usize),
    SigmaInv(usize),
}

impl From<isize> for BrGen {
    fn from(s: isize) -> Self {
        match s {
            ..0 => Self::SigmaInv(s.abs() as _),
            0 => panic!("Invalid s given to from_sigmas"),
            1.. => Self::Sigma(s as _),
        }
    }
}

impl BrGen {
    pub fn inverse(&self) -> Self {
        match self {
            Self::Sigma(a) => Self::SigmaInv(*a),
            Self::SigmaInv(a) => Self::Sigma(*a),
        }
    }

    pub fn shift(&self, n: usize) -> Self {
        match *self {
            Self::Sigma(i) => Self::Sigma(n - i),
            Self::SigmaInv(i) => Self::SigmaInv(n - i),
        }
    }

    pub fn permute(&self, v: &mut Vec<usize>) {
        let Self::Sigma(a) = self else {
            panic!("The braid given was not positive!");
        };
        v.swap_(*a + 1, *a);
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Braid {
    pub gens: Vec<BrGen>,
    // Our braid is an element of B_n
    pub n: usize,
}

impl std::ops::Mul for Braid {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        debug_assert_eq!(
            self.n, other.n,
            "Attempted to compose two different sized braids!"
        );

        Self {
            gens: self.iter().chain(other.iter()).cloned().collect(),
            n: self.n,
        }
    }
}

#[inline]
fn do_swap(fwd: &mut Vec<usize>, inv: &mut Vec<usize>, i: usize) {
    // Swap elements in inv
    inv.swap(i - 1, i);

    // Swap elements in fwd, adjust for 1-based index
    fwd.swap(inv[i - 1] - 1, inv[i] - 1);
}

impl Permutation for Braid {
    fn id(n: usize) -> Self {
        Self { n, gens: vec![] }
    }

    fn size(&self) -> usize {
        self.n as _
    }

    fn swap_(&mut self, _a: usize, _b: usize) {
        unimplemented!()
    }

    fn position(&self, x: usize) -> usize {
        let mut ret = x;
        for g in &self.gens {
            let a = *match g {
                BrGen::Sigma(val) | BrGen::SigmaInv(val) => val,
            };
            if ret == a {
                ret += 1;
            } else if ret == a + 1 {
                ret -= 1;
            }
        }
        ret
    }

    fn as_vec(&self) -> Vec<usize> {
        let mut starting: Vec<usize> = (1..=self.n).collect();
        self.permute(&mut starting);
        starting
    }

    // http://hackage.haskell.org/package/combinat-0.2.8.2/docs/src/Math-Combinat-Groups-Braid.html
    // O(n^2) where n is the length of the permutation
    // Allow many single char names since it's directly adapted from the Haskell source
    #[allow(many_single_char_names)]
    fn from_slice(perm: &[usize]) -> Self {
        // Assuming that perm is a valid permutation
        let n = perm.len();
        let mut fwd: Vec<_> = (1..=n).collect();
        let mut inv: Vec<_> = fwd.clone();
        let mut gens = Vec::with_capacity(n);

        // Pass a reference of cfwdRef each time to doswap
        // O(1)

        //O(n^2)
        for phase in 1..n {
            let target = perm[phase - 1];
            let source = fwd[target - 1];
            // Note that this is immutable
            for num in (phase..source).rev() {
                do_swap(&mut fwd, &mut inv, num);
                gens.push(num);
            }
        }

        Self::from_positive_sigmas(&gens, n)
    }
}

impl Braid {
    pub fn from_positive_sigmas(sigmas: &[usize], n: usize) -> Self {
        let gens = sigmas.iter().map(|s| BrGen::Sigma(*s)).collect();
        Self { gens, n }
    }

    pub fn from_sigmas(sigmas: &[isize], n: usize) -> Self {
        let gens = sigmas.iter().cloned().map(<_>::from).collect();
        Self { gens, n }
    }

    pub fn make_half_twist(n: usize) -> Self {
        // TODO: Use with_capacity
        let gens = (1..n + 1)
            .rev()
            .flat_map(|k| (1..k).map(BrGen::Sigma))
            .collect();

        Self { gens, n }
    }

    pub fn iter(&self) -> std::slice::Iter<'_, BrGen> {
        self.gens.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, BrGen> {
        self.gens.iter_mut()
    }

    pub fn invert(&mut self) {
        self.gens.reverse();
        for g in self.iter_mut() {
            *g = g.inverse();
        }
    }

    pub fn inverse(&self) -> Self {
        let mut ret = self.clone();
        ret.invert();
        ret
    }

    pub fn shift(&mut self) {
        for g in &mut self.gens {
            *g = g.shift(self.n);
        }
    }

    /**
     * Find the starting set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     * NOTE: See Proposition 2.4 (2) of Elrifai Morton for info
     */
    pub fn starting_set(&self) -> IndexSet<usize> {
        let n = self.n;
        let mut res = IndexSet::with_capacity(n);
        let mut string_pos = VecPermutation::id(n);
        // Iterate through each of our generators
        for g in &self.gens {
            let BrGen::Sigma(a) = g else {
                panic!("The braid given was not positive!");
            };
            let sa = string_pos[*a - 1];
            let sb = string_pos[*a];
            if sa == sb + 1 {
                res.insert(sb);
            } else if sb == sa + 1 {
                res.insert(sa);
            }
            // swap the strings
            string_pos.swap_(*a + 1, *a);
        }
        res
    }

    pub fn permute(&self, v: &mut Vec<usize>) {
        for g in &self.gens {
            g.permute(v);
        }
    }

    /**
     * Find the finishing set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     */
    pub fn finishing_set(&self) -> IndexSet<usize> {
        let mut p = VecPermutation::id(self.n);
        self.permute(&mut p);
        // Iterate through each of our generators
        (1..self.n).filter(|i| p[*i] < p[*i - 1]).collect()
    }

    pub fn as_vec_ser(&self) -> Vec<u8> {
        serialize(&self).unwrap()
    }

    pub fn from_vec_ser(vec: &[u8]) -> Self {
        deserialize(vec).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use crate::braid::*;

    #[test]
    fn permut_tests() {
        // Based off of the Haskell permutationBraid
        let p = vec![3, 4, 1, 2];
        let b = Braid::from_slice(&p[..]);
        let b2 = Braid::from_sigmas(&[2, 1, 3, 2], 4);
        assert_eq!(b.gens, b2.gens);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b = Braid::from_slice(&p[..]);
        let b2 = Braid::from_sigmas(&[2, 6, 5, 4, 3, 5], 7);
        assert_eq!(b.gens, b2.gens);
    }

    #[test]
    fn starting_finishing_set_tests() {
        let b = Braid::make_half_twist(6);
        println!("{:?}", b.starting_set());
        println!("{:?}", b.finishing_set());
        println!();
        let b1 = Braid::from_sigmas(&[2, 1, 3, 2, 1], 4);
        let b2 = Braid::from_sigmas(&[1, 2], 4);

        println!("{:?}, {:?}", b1.starting_set(), b2.starting_set());
        println!("{:?}, {:?}", b1.finishing_set(), b2.finishing_set());
        println!();

        let b1 = Braid::from_sigmas(&vec![1, 2], 3);
        let b2 = Braid::from_sigmas(&vec![2, 1, 2], 3);
        println!("{:?}, {:?}", b1.starting_set(), b2.starting_set());
        println!("{:?}, {:?}", b1.finishing_set(), b2.finishing_set());
    }

    #[test]
    fn test_serialization() {
        //TODO expand test cases, improve memory usage
        let b = Braid::random_positive(60, 30, 3, 0.3);
        // let b = Braid::from_sigmas(&[1,2], 3);
        let as_vec = b.as_vec_ser();
        let from_vec = Braid::from_vec_ser(&as_vec);
        let as_string = format!("{:?}", b);
        println!("{:?}", as_vec);
        println!("{:?}", b.gens);
        println!("as_string: {}, as_vec: {}", as_string.len(), as_vec.len());
        assert_eq!(b, from_vec);
    }
}
