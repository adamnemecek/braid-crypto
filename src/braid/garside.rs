#![allow(dead_code)] // Temporary while implimenting the full normal form code

use braid::internals::*;
use std::collections::HashSet;
use std::fmt;

type Permutation = Vec<usize>;

pub struct GarsideForm {
    delta_exp: usize,
    permutations: Vec<Permutation>,
}

impl fmt::Display for GarsideForm {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let vec_string: String = format!("{:?}", self.permutations);
        let vec_string: String = vec_string[1..vec_string.len() - 1].to_string();
        let vec_string = vec_string.replace("[", "(").replace("]", ")");
        write!(f, "[{};{}]", self.delta_exp, vec_string)
    }
}

/**
 * Calculate what B_i should be given that sigma_i^-1 = delta_n^-1 B_i
 * and that B_i is a permutation Braid
 * NOTE: This is based on step 1 of Garside Normal Form
 */
fn neg_pow_to_permute(i:usize, n:usize) -> Braid {
    let mut my_permutation: Vec<usize> = (1..=n).rev().collect();
    // Swap i - 1 with i
    let a = my_permutation[i];
    let b = my_permutation[i - 1];
    my_permutation[i] = b;
    my_permutation[i - 1] = a;
    Braid::from_permutation(my_permutation)
}

// Decompose a braid into a power of delta_n and a
// positive braid on the right.
/**
 * Decompose a braid into a power of delta_n and a positive braid on the right.
 * NOTE: This is based on step 1 and 2 of Garside Normal Form
 */
fn left_slide_delta_form(b: &Braid) -> (usize, Braid) {
    let n = b.n;
    let mut final_vec: Vec<BrGen> = b.contents.clone();
    let mut counter = 0;
    let mut acting_index = 0;
    while acting_index < final_vec.len() {
        if let BrGen::SigmaInv(i) = final_vec[acting_index] {
            // This generator needs to be replaced
            let mut loc_of_delta = (acting_index - 1) as isize;
            // remove the "bad" inverse generator
            final_vec.remove(acting_index);
            // replacement = the replacement (minus the delta)
            let replacement = neg_pow_to_permute(i, n);
            for symb in replacement.contents.iter() {
                final_vec.insert(acting_index, *symb);
                acting_index += 1;
            }
            // Now go backwards and replace sigma_a with sigma_{n - a}
            while loc_of_delta != -1 {
                if let BrGen::Sigma(a) = final_vec[loc_of_delta as usize] {
                    final_vec[loc_of_delta as usize] = BrGen::Sigma(n - a);
                } else {
                    panic!("There was a negative sigma?");
                }
                loc_of_delta -= 1;
            }
            counter += 1;
        } else {
            // Skip this generator and move on the check the next
            acting_index += 1;
        }
    }

    (counter, Braid {contents: final_vec, n: b.n})
}

/**
 * Break a braid into a series of permutation braids Q_i
 * Choose the longest permutation braids possible.
 * Original algorithm by me
 */
pub fn break_into_permutations(b: &Braid) -> Vec<Braid> {
    let n = b.n;
    // String at position i is string # string_pos[i - 1]
    let mut string_pos: Vec<usize> = (1..=n).collect();
    // String i has crossed String j if has_crossed.contains((i, j)) is true (and i < j)
    let mut has_crossed: HashSet<(usize, usize)> = HashSet::new();
    // Our current index in breaking the original braid
    let mut working_index = 0;
    // The result we'll return
    let mut res: Vec<Braid> = vec![];
    let mut tmp_to_add: Vec<BrGen> = vec![];

    // A helpful function for filtering our original braid into what we need
    let sym_to_i = |sym: &BrGen| {
        if let BrGen::Sigma(a) = *sym {
            return a;
        }
        panic!("The given braid was not positive");
    };

    let symbols: Vec<usize> = b.contents.iter().map(sym_to_i).collect();

    // The algorithm
    while working_index < symbols.len() {
        let swap = symbols[working_index];
        let string1_name = string_pos[swap - 1];
        let string2_name = string_pos[swap];
        // Have these strings crossed before?
        if has_crossed.contains(&(string1_name, string2_name))
        || has_crossed.contains(&(string2_name, string1_name)) {
            // They have. Let's seperate into a new Q
            res.push(Braid {contents: tmp_to_add.clone(), n:n});
            tmp_to_add.clear();
            // Clear the has_crossed set
            has_crossed.clear();
        } else {
            // They have not. Now they have
            has_crossed.insert((string1_name, string2_name));
        }

        tmp_to_add.push(BrGen::Sigma(swap));

        // Update the string_pos with the swap
        string_pos[swap] = string1_name;
        string_pos[swap - 1] = string2_name;

        working_index += 1;
    }

    if tmp_to_add.len() != 0 {
        res.push(Braid {contents: tmp_to_add, n:n});
    }

    return res;
}

pub fn is_left_weighted(b: &Braid) -> bool {
    let ps = break_into_permutations(b);
    for i in 0..ps.len() - 1 {
        if !ps[i].finishing_set().is_superset(&ps[i + 1].starting_set()) {
            return false;
        }
    }
    true
}

fn braid_to_permutation(b: &Braid) -> Vec<usize> {
    let mut starting: Vec<usize> = (1..=b.n).collect();
    braid_to_permutation_with_starting(b, &mut starting);
    starting
}

fn braid_to_permutation_with_starting(b: &Braid, starting: &mut Vec<usize>){
    // let mut string_pos: Vec<usize> = (1..=n).collect();
    let string_pos = starting;
    // Iterate through each of our generators
    for gen in b.contents.iter() {
        if let BrGen::Sigma(a) = gen {
            let sa = string_pos[*a - 1];
            let sb = string_pos[*a];
            // swap the strings
            string_pos[*a] = sa;
            string_pos[*a - 1] = sb;
        } else {
            panic!("The braid given was not positive!");
        }
    }
}

pub fn garside_form(b: &Braid) -> GarsideForm {
    let n = b.n;
    let (exponent, braid) = left_slide_delta_form(b);
    let mut bs = break_into_permutations(&braid);
    let mut working_index = 0;
    while working_index < bs.len() - 1 {
        // Split bs to take multiple mutable references
        let (head, tail) = bs.split_at_mut(working_index + 1); 
        let bi = &mut head[working_index];
        let bi1 = &mut tail[0];
        let mut next_starting = bi1.starting_set();
        let mut prev_finishing = bi.finishing_set();
        while !prev_finishing.is_superset(&next_starting) {
            {
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
                braid_to_permutation_with_starting(bi1, &mut perm);
                // Now, we turn it back into a permutation braid
                let pb = Braid::from_permutation(perm);
                // and replace bi1 with it
                bi1.contents = pb.contents.clone();
            } // For j to go out of scope (j borrows bi1 and bi)
            next_starting = bi1.starting_set();
            prev_finishing = bi.finishing_set();
        }
        working_index += 1;
    }
    // Recombine bs
    let mut result: Vec<Permutation> = Vec::new();
    for braid in &mut bs {
        result.push(braid_to_permutation(braid));
    }
    GarsideForm {delta_exp: exponent, permutations: result}
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn neg_pow_to_perm_tests() {
        // Based on example 1.2 from https://arxiv.org/pdf/0711.3941.pdf
        let expected = Braid::make_positive(vec![3, 2, 1, 3, 2], 4);
        let actual = neg_pow_to_permute(3, 4);
        assert_eq!(expected.contents, actual.contents);
    }

    #[test]
    fn left_slide_delta_form_tests() {
        // Based on example 1.2 from https://arxiv.org/pdf/0711.3941.pdf
        let contents = vec![BrGen::Sigma(1), BrGen::SigmaInv(3), BrGen::Sigma(2)];
        let w = Braid {contents: contents, n: 4};
        let lsdf = left_slide_delta_form(&w);
        assert_eq!(1, lsdf.0);
        let expected = Braid::make_positive(vec![3, 3, 2, 1, 3, 2, 2], 4);
        assert_eq!(expected.contents, lsdf.1.contents);
    }

    #[test]
    fn break_into_permutations_tests() {
        let b = Braid::make_positive(vec![1, 2, 2, 1, 2], 3);
        let ps = break_into_permutations(&b);
        assert_eq!(ps[0].contents, Braid::make_positive(vec![1, 2], 3).contents);
        assert_eq!(ps[1].contents, Braid::make_positive(vec![2, 1, 2], 3).contents);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b = Braid::from_permutation(p);
        let old_contents = b.contents.clone();
        let ps = break_into_permutations(&b);
        assert_eq!(ps.len(), 1);
        assert_eq!(ps[0].contents, old_contents);
    }

    #[test]
    fn is_left_weighted_tests() {
        // From Page 12/13 of Garber
        assert!(!is_left_weighted(&Braid::make_positive(vec![1, 2, 2, 1, 2], 3)));
        assert!(is_left_weighted(&Braid::make_positive(vec![1, 2, 2, 1], 3)));

        assert!(!is_left_weighted(&Braid::make_positive(vec![3, 3, 2, 1, 3, 2, 2], 4)));
        assert!(is_left_weighted(&Braid::make_positive(vec![2, 1, 3, 2, 1, 1, 2], 4)));
    }

    #[test]
    fn garside_form_tests() {
        // From page 13/14 of Garber
        let b = Braid {contents: vec![BrGen::Sigma(1), BrGen::SigmaInv(3), BrGen::Sigma(2)], n:4};
        let gform = garside_form(&b);
        
        let expected_exp = 1;
        let expected_perm1 = braid_to_permutation(&Braid::make_positive(vec![2, 1, 3, 2, 1], 4));
        let expected_perm2 = braid_to_permutation(&Braid::make_positive(vec![1, 2], 4));

        assert_eq!(expected_exp, gform.delta_exp);
        assert_eq!(expected_perm1, gform.permutations[0]);
        assert_eq!(expected_perm2, gform.permutations[1]);
        println!("{}", gform);
    }
}