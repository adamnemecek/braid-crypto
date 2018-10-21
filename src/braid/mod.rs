pub mod garside;
pub mod random;

use std::ops::Mul;
use bincode::{serialize, deserialize};
use indexmap::set::IndexSet;
use self::BrGen::*;

type BSize = u16;

type Permutation = Vec<usize>;

#[derive(Debug, Copy, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum BrGen {
    Sigma(BSize),
    SigmaInv(BSize)
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Braid {
    pub contents: Vec<BrGen>,
    pub n: BSize // Our braid is an element of B_n
}

impl Mul for Braid {
    type Output = Braid;

    fn mul(self, other: Braid) -> Braid {
        debug_assert_eq!(self.n, other.n, "Attempted to compose two different sized vectors!");

        let mut new_contents = self.contents.clone();
        let mut new_other_contents = other.contents.clone();
        new_contents.append(&mut new_other_contents);
        
        let mut ret = Braid { contents: new_contents, n: self.n};
        ret.free_reduce();
        ret
    }
}

fn invert_gens(gens: &mut Vec<BrGen>) {
    gens.reverse();
    for gen in gens {
        let new_gen = match gen {
            BrGen::Sigma(a) => BrGen::SigmaInv(*a),
            BrGen::SigmaInv(a) => BrGen::Sigma(*a),
        };
        *gen = new_gen;
    }
}

impl Braid {
    // http://hackage.haskell.org/package/combinat-0.2.8.2/docs/src/Math-Combinat-Groups-Braid.html
    // O(n^2) where n is the length of the permutation
    pub fn from_permutation(perm: Permutation) -> Braid {
        // Assuming that perm is a valid permutation
        let n = perm.len();
        let mut cfwd: Permutation = (1..=n).collect();
        let mut cinv: Permutation = (1..=n).collect();
        let mut contents: Vec<BSize> = Vec::with_capacity(n);

        // Pass a reference of cfwdRef each time to doswap
        // O(1)
        let mut do_swap = |i: usize, cfwd_ref: &mut Permutation| {
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
            let this: Permutation = (phase..source).rev().collect();
            for num in this {
                do_swap(num, &mut cfwd);
                contents.push(num as BSize);
            }
            phase += 1;
        }

        Braid::from_positive_sigmas(contents, n as BSize)
    }

    pub fn from_positive_sigmas(sigmas: Vec<BSize>, n:BSize) -> Braid {
        let contents = sigmas.iter().map(|s| {
            Sigma(*s)
        }).collect();
        Braid {contents: contents, n: n}
    }

    pub fn from_sigmas(sigmas: Vec<isize>, n:BSize) -> Braid {
        let contents = sigmas.iter().map(|s| {
            if *s < 0 {
                SigmaInv((*s).abs() as BSize)
            } else if *s > 0 {
                Sigma(*s as BSize)
            } else {
                panic!("Invalid s given to from_sigmas")
            }
        }).collect();
        Braid {contents: contents, n: n}
    }

    pub fn make_half_twist(n: BSize) -> Braid {
        // TODO: Use with_capacity
        let mut contents: Vec<BrGen> = Vec::new();

        for k in (1..n+1).rev() {
            for j in 1..k {
                contents.push(BrGen::Sigma(j));
            }
        }

        Braid {contents: contents, n: n}
    }

    pub fn invert(&mut self) {
        invert_gens(&mut self.contents);
    }

    pub fn inverse(&self) -> Braid {
        let mut new_contents = self.contents.clone();
        invert_gens(&mut new_contents);
        Braid {contents: new_contents, n:self.n}
    }

    pub fn shift(&mut self) {
        for gen in &mut self.contents {
            *gen = match gen {
                BrGen::Sigma(i) => BrGen::Sigma(self.n - *i),
                BrGen::SigmaInv(i) => BrGen::SigmaInv(self.n - *i)
            };
        }
    }

    pub fn free_reduce_once(&mut self) -> bool {
        if self.contents.len() == 0 {
            return false; 
        }
        let mut changed = false;
        let mut new_contents = Vec::with_capacity(self.contents.len());
        let mut i = 0;

        while i < self.contents.len() - 1 {
            if let (BrGen::Sigma(a), BrGen::SigmaInv(b)) = (self.contents[i], self.contents[i + 1]) {
                if a == b {
                    changed = true;
                    i += 2;
                    continue;
                }
            }
            new_contents.push(self.contents[i]);
            i += 1;
        }

        self.contents = new_contents;
        changed
    }

    pub fn free_reduce(&mut self) {
        while self.free_reduce_once() {}
    }

    /**
     * Find the starting set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     * ASSUMPTION: i is in the starting set iff strand i and strand i + 1
     * cross once.
     */
    pub fn starting_set(&self) -> IndexSet<BSize> {
        let n = self.n as usize;
        let mut res = IndexSet::with_capacity(n);
        let mut string_pos: Permutation = (1..=n).collect();
        // Iterate through each of our generators
        for gen in self.contents.iter() {
            if let BrGen::Sigma(a) = gen {
                let sa = string_pos[(*a - 1) as usize];
                let sb = string_pos[(*a) as usize];
                if sa == sb + 1 {
                    res.insert(sb as BSize);
                } else if sb == sa + 1{
                    res.insert(sa as BSize);
                }
                // swap the strings
                string_pos[(*a) as usize] = sa;
                string_pos[(*a - 1) as usize] = sb;
            } else {
                panic!("The braid given was not positive!");
            }
        }
        res
    }

    /**
     * Find the finishing set of a braid. Original algorithm. Assumes
     * this is a positive permutation braid.
     * ASSUMPTION: i is in the finishing set iff the strands which end up at
     * indices i and i + 1 cross once
     */
    pub fn finishing_set(&self) -> IndexSet<BSize> {
        let n = self.n as usize;
        let mut res = IndexSet::with_capacity(n);
        let mut string_pos: Permutation = (1..=n).collect();
        // Iterate through each of our generators
        for gen in self.contents.iter() {
            if let BrGen::Sigma(a) = gen {
                let sa = string_pos[(*a - 1) as usize];
                let sb = string_pos[(*a) as usize];
                // swap the strings
                string_pos[(*a) as usize] = sa;
                string_pos[(*a - 1) as usize] = sb;
            } else {
                panic!("The braid given was not positive!");
            }
        }

        for i in 1..n {
            if string_pos[i] < string_pos[i - 1] {
                res.insert(i as BSize);
            }
        }

        res
    }

    pub fn as_vec(&self) -> Vec<u8> {
        serialize(&self).unwrap()
    }

    pub fn from_vec(vec: &Vec<u8>) -> Braid {
        deserialize(vec).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reduction_tests() {
        // TODO: Do more tests for reduction
        // I don't believe it's doing it correctly
        let b = Braid::from_sigmas(vec![1, 2, 3], 3);
        let b2 = b.inverse();
        let mut b3 = b * b2;
        b3.free_reduce();
        assert!(b3.contents.len() == 0);
    }

    #[test]
    fn permut_tests() {
        // Based off of the Haskell permutationBraid
        let p = vec![3, 4, 1, 2];
        let b = Braid::from_permutation(p);
        let b2 = Braid::from_sigmas(vec![2,1,3,2], 4);
        assert_eq!(b.contents, b2.contents);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b = Braid::from_permutation(p);
        let b2 = Braid::from_sigmas(vec![2, 6, 5, 4, 3, 5], 7);
        assert_eq!(b.contents, b2.contents);
    }

    #[test]
    fn starting_finishing_set_tests() {
        let b = Braid::make_half_twist(6);
        println!("{:?}", b.starting_set());
        println!("{:?}", b.finishing_set());
        println!("");
        let b1 = Braid::from_sigmas(vec![2,1,3,2,1], 4);
        let b2 = Braid::from_sigmas(vec![1,2], 4);

        println!("{:?}, {:?}", b1.starting_set(), b2.starting_set());
        println!("{:?}, {:?}", b1.finishing_set(), b2.finishing_set());
        println!("");

        let b1 = Braid::from_sigmas(vec![1,2], 3);
        let b2 = Braid::from_sigmas(vec![2,1,2], 3);
        println!("{:?}, {:?}", b1.starting_set(), b2.starting_set());
        println!("{:?}, {:?}", b1.finishing_set(), b2.finishing_set());
    }

    #[test]
    fn test_serialization() {
        //TODO expand test cases, improve memory usage
        let b = Braid::random_positive(60, 30, 3, 0.3);
        // let b = Braid::from_sigmas(vec![1,2], 3);
        let as_vec = b.as_vec();
        let from_vec = Braid::from_vec(&as_vec);
        let as_string = format!("{:?}", b);
        println!("{:?}", as_vec);
        println!("{:?}", b.contents);
        println!("as_string: {}, as_vec: {}", as_string.len(), as_vec.len());
        assert_eq!(b, from_vec);
    }
}