#![allow(dead_code)] // Temporary while implimenting the full normal form code

// pub use braid_crypto::prelude::*;
use std::{
    collections::HashSet,
    fmt,
};

// use BrGen::*;
use crate::{
    braid::*,
    permutation::*,
};

pub struct GarsideForm {
    delta_exp: isize,
    permutations: Vec<VecPermutation>,
}

impl fmt::Display for GarsideForm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let vec_string: String = format!("{:?}", self.permutations);
        let vec_string: String = vec_string[1..vec_string.len() - 1].to_string();
        let vec_string = vec_string.replace("[", "(").replace("]", ")");
        write!(f, "[{};{}]", self.delta_exp, vec_string)
    }
}

impl Braid {
    /**
     * Calculate what B_i should be given that sigma_i^-1 = delta_n^-1 B_i
     * and that B_i is a permutation Braid
     * NOTE: This is based on step 1 of Garside Normal Form
     * O(n^2)
     */
    fn neg_pow_to_permute(i: usize, n: usize) -> Self {
        let mut p: VecPermutation = (1..=n).rev().collect();
        // Swap i + 1 with i
        p.swap(i + 1, i);
        Self::from_slice(&p[..])
    }
}

/**
 * Decompose a braid into a power of delta_n and a positive braid on the right.
 * NOTE: This is based on step 1 and 2 of Garside Normal Form
 * O(L*n^2 + L^2) where L is the length of b
 */
fn left_slide_delta_form(b: &Braid) -> (isize, Braid) {
    let n = b.n;
    let mut final_vec: Vec<BrGen> = b.contents.clone();
    let mut counter = 0;
    let mut acting_index = 0;
    // O(Ln^2 + L^2)
    while acting_index < final_vec.len() {
        if let BrGen::SigmaInv(i) = final_vec[acting_index] {
            // This generator needs to be replaced
            let mut loc_of_delta = acting_index as isize - 1;
            // remove the "bad" inverse generator
            final_vec.remove(acting_index);
            // replacement = the replacement (minus the delta)
            // O(n^2)
            let replacement = Braid::neg_pow_to_permute(i, n);
            for symb in &replacement.contents {
                final_vec.insert(acting_index, *symb);
                acting_index += 1;
            }
            // Now go backwards and replace sigma_a with sigma_{n - a}
            // O(L) where L is the length of b in terms of generators
            while loc_of_delta != -1 {
                let BrGen::Sigma(a) = final_vec[loc_of_delta as usize] else {
                    panic!("There was a negative sigma?");
                };
                final_vec[loc_of_delta as usize] = BrGen::Sigma(n - a);
                loc_of_delta -= 1;
            }
            counter -= 1;
        } else {
            // Skip this generator and move on the check the next
            acting_index += 1;
        }
    }

    (
        counter,
        Braid {
            contents: final_vec,
            n: b.n,
        },
    )
}

// A helpful function for filtering our original braid into what we need
fn sym_to_i(sym: &BrGen) -> usize {
    let BrGen::Sigma(a) = *sym else {
        panic!("The given braid was not positive");
    };
    a
}
/**
 * Break a braid into a series of permutation braids Q_i
 * Choose the longest permutation braids possible.
 * Original algorithm by me
 * O(L)
 */
pub fn break_into_permutations(b: &Braid) -> Vec<Braid> {
    let n = b.n;
    // String at position i is string # string_pos[i - 1]
    let mut string_pos = VecPermutation::id(n);
    // String i has crossed String j if has_crossed.contains((i, j)) is true (and i < j)
    let mut has_crossed = HashSet::<(usize, usize)>::new();
    // Our current index in breaking the original braid
    let mut working_index = 0;
    // The result we'll return
    let mut res: Vec<Braid> = vec![];
    let mut contents: Vec<BrGen> = vec![];

    let symbols: Vec<usize> = b.contents.iter().map(sym_to_i).collect();

    // The algorithm
    // O(L)
    while working_index < symbols.len() {
        let swap = symbols[working_index];
        let string1_name = string_pos[swap - 1];
        let string2_name = string_pos[swap];
        // Have these strings crossed before?
        if has_crossed.contains(&(string1_name, string2_name))
            || has_crossed.contains(&(string2_name, string1_name))
        {
            // They have. Let's seperate into a new Q
            res.push(Braid {
                contents: contents.clone(),
                n,
            });
            contents.clear();
            // Clear the has_crossed set
            has_crossed.clear();
            // New: these strings have still crossed now
            has_crossed.insert((string1_name, string2_name));
        } else {
            // They have not. Now they have
            has_crossed.insert((string1_name, string2_name));
        }

        contents.push(BrGen::Sigma(swap));

        // Update the string_pos with the swap
        string_pos.swap(swap, swap + 1);

        working_index += 1;
    }

    if !contents.is_empty() {
        res.push(Braid { contents, n });
    }

    res
}

impl Braid {
    // TODO: Re-evaluate O() of as_garside form with workingindex -= 1 change
    // O(Ln^2 + L^2 + p(L^2 + Ln^2)) where p is the number of permutations which make up self
    pub fn as_garside_form(&self) -> GarsideForm {
        let n = self.n;
        // O(L*n^2 + L^2)
        let (exponent, braid) = left_slide_delta_form(&self);
        // O(L)
        let mut bs = break_into_permutations(&braid);
        let mut working_index = 0;
        while working_index < bs.len() - 1 {
            let mut changed = false;
            {
                // Scope for slices of bs
                // Split bs to take multiple mutable references
                let (head, tail) = bs.split_at_mut(working_index + 1);
                let bi = &mut head[working_index];
                let bi1 = &mut tail[0];
                // O(L)
                let mut next_starting = bi1.starting_set();
                // O(L)
                let mut prev_finishing = bi.finishing_set();
                // O(kL + kn^2) where k is the average runtime until a superset is found (probably small)
                // If we assume k = L for an upper bound, then this is O(L^2 + Ln^2)
                while !prev_finishing.is_superset(&next_starting) {
                    {
                        changed = true;
                        let j = next_starting.difference(&prev_finishing).next().unwrap();
                        // j is in S(B_i+1) but not F(B_i)
                        // bi is easy, just push a sigma on the end
                        bi.contents.push(BrGen::Sigma(*j));
                        // bi1 is harder.
                        // We want to put a sigma -j on the beginning, but we want it to stay positive
                        // Instead, let's consider bi1 as a permutation with j and j + 1 switched
                        let mut perm: Vec<usize> = (1..=n).collect();
                        perm[*j - 1] = *j + 1;
                        perm[*j] = *j;
                        // TODO: O(?)
                        braid_to_permutation_with_starting(bi1, &mut perm);

                        // Now, we turn it back into a permutation braid
                        // O(n^2)
                        let pb = Braid::from_slice(&perm[..]);
                        // and replace bi1 with it
                        bi1.contents = pb.contents.clone();
                    } // For j to go out of scope (j borrows bi1 and bi)
                      // O(L)
                    next_starting = bi1.starting_set();
                    // O(L)
                    prev_finishing = bi.finishing_set();
                }
            } // Slices of bs go out of scope
            if bs[working_index + 1].is_identity() {
                // bi1 is the identity permutation
                // don't add it
                bs.remove(working_index + 1);
            }
            if changed && working_index != 0 {
                // Go back to check if the starting set of this was
                working_index -= 1;
            } else {
                working_index += 1;
            }
        }

        // Recombine bs
        // New: remove intial factors of Delta
        let mut delta_exp = exponent;
        let mut result: Vec<VecPermutation> = Vec::new();
        for braid in &mut bs {
            let perm = braid.as_vec();
            if perm.is_twist() {
                delta_exp += 1;
            } else {
                result.push(perm);
            }
        }

        GarsideForm {
            delta_exp,
            permutations: result,
        }
    }

    pub fn is_left_weighted(&self) -> bool {
        let ps = break_into_permutations(self);
        for i in 0..ps.len() - 1 {
            if !ps[i].finishing_set().is_superset(&ps[i + 1].starting_set()) {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn neg_pow_to_perm_tests() {
        // Based on example 1.2 from https://arxiv.org/pdf/0711.3941.pdf
        let expected = Braid::from_sigmas(&[3, 2, 1, 3, 2], 4);
        let actual = Braid::neg_pow_to_permute(3, 4);
        assert_eq!(expected.contents, actual.contents);
    }

    #[test]
    fn left_slide_delta_form_tests() {
        // Based on example 1.2 from https://arxiv.org/pdf/0711.3941.pdf
        let w = Braid::from_sigmas(&[1, -3, 2], 4);
        let lsdf = left_slide_delta_form(&w);
        assert_eq!(-1, lsdf.0);
        let expected = Braid::from_sigmas(&vec![3, 3, 2, 1, 3, 2, 2], 4);
        assert_eq!(expected.contents, lsdf.1.contents);
    }

    #[test]
    fn break_into_permutations_tests() {
        let b = Braid::from_sigmas(&[1, 2, 2, 1, 2], 3);
        let ps = break_into_permutations(&b);
        assert_eq!(ps[0].contents, Braid::from_sigmas(&[1, 2], 3).contents);
        assert_eq!(ps[1].contents, Braid::from_sigmas(&[2, 1, 2], 3).contents);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b = Braid::from_slice(&p[..]);
        let old_contents = b.contents.clone();
        let ps = break_into_permutations(&b);
        assert_eq!(ps.len(), 1);
        assert_eq!(ps[0].contents, old_contents);
    }

    #[test]
    fn new_break_into_permutations_tests() {
        let a = Braid::from_sigmas(&[2, 2, 2, 2], 3);
        println!("{:?}", break_into_permutations(&a));
    }

    #[test]
    fn is_left_weighted_tests() {
        // From Page 12/13 of Garber
        assert!(!Braid::from_sigmas(&[1, 2, 2, 1, 2], 3).is_left_weighted());
        assert!(Braid::from_sigmas(&[1, 2, 2, 1], 3).is_left_weighted());

        assert!(!Braid::from_sigmas(&[3, 3, 2, 1, 3, 2, 2], 4).is_left_weighted());
        assert!(Braid::from_sigmas(&[2, 1, 3, 2, 1, 1, 2], 4).is_left_weighted());
    }

    #[test]
    fn garside_form_tests() {
        // From page 13/14 of Garber
        let b = Braid::from_sigmas(&[1, -3, 2], 4);
        let gform = b.as_garside_form();

        let expected_exp = -1;
        let expected_perm1 = Braid::from_sigmas(&[2, 1, 3, 2, 1], 4).as_vec();
        let expected_perm2 = Braid::from_sigmas(&[1, 2], 4).as_vec();

        assert_eq!(expected_exp, gform.delta_exp);
        assert_eq!(expected_perm1, gform.permutations[0]);
        assert_eq!(expected_perm2, gform.permutations[1]);
        println!("{}", gform);
    }

    #[test]
    fn more_garside_tests() {
        // To test breaking into permutations
        let a1 = Braid::from_sigmas(&[2, 1, 2, 1, 2], 3);
        let a3 = Braid::from_sigmas(&[2, 2, 1, 2, 2], 3);

        assert_eq!(
            format!("{}", a1.as_garside_form()),
            format!("{}", a3.as_garside_form())
        );

        let a1 = Braid::from_sigmas(&[2, 1, 2, 1, 2, 2, 2], 3);
        let a3 = Braid::from_sigmas(&[2, 2, 1, 2, 2, 2, 2], 3);

        assert_eq!(
            format!("{}", a1.as_garside_form()),
            format!("{}", a3.as_garside_form())
        );

        // To test removing final twists
        let a1 = Braid::from_sigmas(&[1, 3, -3, 2, 1], 4);
        let a3 = Braid::from_sigmas(&[2, 1, 2], 4);

        assert_eq!(
            format!("{}", a1.as_garside_form()),
            format!("{}", a3.as_garside_form())
        );
    }
}
