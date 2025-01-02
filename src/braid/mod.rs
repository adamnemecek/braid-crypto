pub mod garside;
pub mod random;

// pub use crate::prelude::*;

use {
    bincode::{
        deserialize,
        serialize,
    },
    indexmap::set::IndexSet,
    std::ops::Mul,
};

use crate::permutation::*;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BrGen {
    Sigma(usize),
    SigmaInv(usize),
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Braid {
    pub contents: Vec<BrGen>,
    pub n: usize, // Our braid is an element of B_n
}

impl Mul for Braid {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        debug_assert_eq!(
            self.n, other.n,
            "Attempted to compose two different sized braids!"
        );

        let mut contents = self.contents.clone();
        let mut new_other_contents = other.contents.clone();
        contents.append(&mut new_other_contents);

        Self {
            contents,
            n: self.n,
        }
    }
}

// O(L)
pub fn braid_to_permutation_with_starting(b: &Braid, starting: &mut Vec<usize>) {
    let string_pos = starting;
    // Iterate through each of our generators
    for g in &b.contents {
        if let BrGen::Sigma(a) = g {
            string_pos.swap((*a + 1) as usize, *a as usize);
        } else {
            panic!("The braid given was not positive!");
        }
    }
}

impl Permutation for Braid {
    fn id(n: usize) -> Self {
        Self {
            n: n as usize,
            contents: vec![],
        }
    }

    fn size(&self) -> usize {
        self.n as usize
    }

    fn swap(&mut self, _a: usize, _b: usize) {
        panic!("unimplimented")
    }

    fn follow_starting(&self, x: usize) -> usize {
        let mut ret = x;
        for g in &self.contents {
            let a = match g {
                BrGen::Sigma(val) => val,
                BrGen::SigmaInv(val) => val,
            };
            if ret == *a {
                ret += 1;
            } else if ret == *a + 1 {
                ret -= 1;
            }
        }
        ret
    }

    fn as_vec(&self) -> Vec<usize> {
        let mut starting: Vec<usize> = (1..=self.n).collect();
        braid_to_permutation_with_starting(self, &mut starting);
        starting
    }

    // http://hackage.haskell.org/package/combinat-0.2.8.2/docs/src/Math-Combinat-Groups-Braid.html
    // O(n^2) where n is the length of the permutation
    // Allow many single char names since it's directly adapted from the Haskell source
    #[allow(many_single_char_names)]
    fn from_slice(perm: &[usize]) -> Self {
        // Assuming that perm is a valid permutation
        let n = perm.len();
        let mut cfwd: Vec<usize> = (1..=n).collect();
        let mut cinv: Vec<usize> = (1..=n).collect();
        let mut contents: Vec<usize> = Vec::with_capacity(n);

        // Pass a reference of cfwdRef each time to doswap
        // O(1)
        let mut do_swap = |i: usize, cfwd_ref: &mut Vec<usize>| {
            let a = cinv[i - 1];
            let b = cinv[i];
            cinv[i - 1] = b;
            cinv[i] = a;

            let u = cfwd_ref[a - 1];
            let v = cfwd_ref[b - 1];
            cfwd_ref[a - 1] = v;
            cfwd_ref[b - 1] = u;
        };

        let mut phase = 1;

        //O(n^2)
        while phase < n {
            let target = perm[(phase - 1) as usize];
            let source = cfwd[(target - 1) as usize];
            // Note that this is immutable
            let this: Vec<usize> = (phase..source).rev().collect();
            for num in this {
                do_swap(num, &mut cfwd);
                contents.push(num as usize);
            }
            phase += 1;
        }

        Self::from_positive_sigmas(&contents, n as usize)
    }
}

fn invert_gens(gens: &mut Vec<BrGen>) {
    gens.reverse();
    for g in gens {
        let new_gen = match g {
            BrGen::Sigma(a) => BrGen::SigmaInv(*a),
            BrGen::SigmaInv(a) => BrGen::Sigma(*a),
        };
        *g = new_gen;
    }
}

impl Braid {
    pub fn from_positive_sigmas(sigmas: &[usize], n: usize) -> Self {
        let contents = sigmas.iter().map(|s| BrGen::Sigma(*s)).collect();
        Self { contents, n }
    }

    pub fn from_sigmas(sigmas: &[isize], n: usize) -> Self {
        let contents = sigmas
            .iter()
            .map(|s| {
                if *s < 0 {
                    BrGen::SigmaInv((*s).abs() as usize)
                } else if *s > 0 {
                    BrGen::Sigma(*s as usize)
                } else {
                    panic!("Invalid s given to from_sigmas")
                }
            })
            .collect();
        Self { contents, n }
    }

    pub fn make_half_twist(n: usize) -> Self {
        // TODO: Use with_capacity
        let mut contents: Vec<BrGen> = Vec::new();

        for k in (1..n + 1).rev() {
            for j in 1..k {
                contents.push(BrGen::Sigma(j));
            }
        }

        Self { contents, n }
    }

    pub fn invert(&mut self) {
        invert_gens(&mut self.contents);
    }

    pub fn inverse(&self) -> Self {
        let mut contents = self.contents.clone();
        invert_gens(&mut contents);
        Self {
            contents,
            n: self.n,
        }
    }

    pub fn shift(&mut self) {
        for g in &mut self.contents {
            *g = match g {
                BrGen::Sigma(i) => BrGen::Sigma(self.n - *i),
                BrGen::SigmaInv(i) => BrGen::SigmaInv(self.n - *i),
            }
        }
    }

    /**
     * Find the starting set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     * NOTE: See Proposition 2.4 (2) of Elrifai Morton for info
     */
    pub fn starting_set(&self) -> IndexSet<usize> {
        let n = self.n as usize;
        let mut res = IndexSet::with_capacity(n);
        let mut string_pos = VecPermutation::id(n);
        // Iterate through each of our generators
        for g in &self.contents {
            if let BrGen::Sigma(a) = g {
                let sa = string_pos[(*a - 1) as usize];
                let sb = string_pos[(*a) as usize];
                if sa == sb + 1 {
                    res.insert(sb as usize);
                } else if sb == sa + 1 {
                    res.insert(sa as usize);
                }
                // swap the strings
                string_pos.swap((*a + 1) as usize, (*a) as usize);
            } else {
                panic!("The braid given was not positive!");
            }
        }
        res
    }

    /**
     * Find the finishing set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     */
    pub fn finishing_set(&self) -> IndexSet<usize> {
        let n = self.n as usize;
        let mut res = IndexSet::with_capacity(n);
        let mut string_pos = VecPermutation::id(n);
        // Iterate through each of our generators
        for g in &self.contents {
            if let BrGen::Sigma(a) = g {
                // swap the strings
                string_pos.swap((*a + 1) as usize, (*a) as usize);
            } else {
                panic!("The braid given was not positive!");
            }
        }

        for i in 1..n {
            if string_pos[i] < string_pos[i - 1] {
                res.insert(i as usize);
            }
        }

        res
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
    use crate::{
        braid::*,
        permutation::*,
    };
    #[test]
    fn permut_tests() {
        // Based off of the Haskell permutationBraid
        let p = vec![3, 4, 1, 2];
        let b: Braid = Permutation::from_slice(&p[..]);
        let b2 = Braid::from_sigmas(&[2, 1, 3, 2], 4);
        assert_eq!(b.contents, b2.contents);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b: Braid = Permutation::from_slice(&p[..]);
        let b2 = Braid::from_sigmas(&[2, 6, 5, 4, 3, 5], 7);
        assert_eq!(b.contents, b2.contents);
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
        println!("{:?}", b.contents);
        println!("as_string: {}, as_vec: {}", as_string.len(), as_vec.len());
        assert_eq!(b, from_vec);
    }
}
