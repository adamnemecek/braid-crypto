#![allow(dead_code)] // Temporary while implimenting the full normal form code

use braid::internals::*;
use std::collections::HashSet;

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
fn left_slide_delta_form(b: Braid) -> (usize, Braid) {
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
pub fn break_into_permutations(b: Braid) -> Vec<Braid> {
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
        let lsdf = left_slide_delta_form(w);
        assert_eq!(1, lsdf.0);
        let expected = Braid::make_positive(vec![3, 3, 2, 1, 3, 2, 2], 4);
        assert_eq!(expected.contents, lsdf.1.contents);
    }

    #[test]
    fn break_into_permutations_tests() {
        let b = Braid::make_positive(vec![1, 2, 2, 1, 2], 3);
        let ps = break_into_permutations(b);
        assert_eq!(ps[0].contents, Braid::make_positive(vec![1, 2], 3).contents);
        assert_eq!(ps[1].contents, Braid::make_positive(vec![2, 1, 2], 3).contents);

        let p = vec![1, 3, 7, 2, 5, 4, 6];
        let b = Braid::from_permutation(p);
        let old_contents = b.contents.clone();
        let ps = break_into_permutations(b);
        assert_eq!(ps.len(), 1);
        assert_eq!(ps[0].contents, old_contents);
    }
}