
use rand::prng::hc128::*;
use rand::SeedableRng;
use rand::RngCore;
use rand::CryptoRng;
use rand::os::OsRng;
use rand::distributions::{Distribution, Standard};
// Import Braid and members
use super::*;

fn make_rng() -> Hc128Rng {
    let mut seed: [u8; 32] = [0; 32];
    let mut osrng = OsRng::new().unwrap();
    osrng.fill_bytes(&mut seed);
    Hc128Rng::from_seed(seed)
}

/**
 * Generate a random permutation of
 * size n (using [1..n])
 * using Rng "rng"
 * complexity is how many "twists" there should be
 * miss rate is the probability of skipping a given twist
 */
fn random_permutation<CR: CryptoRng + RngCore>(n: usize, complexity: usize, miss_rate: f32, rng: &mut CR) -> VecPermutation {
    let mut perm: VecPermutation = Permutation::id(n);
    for _ in 0..complexity {
        let drawing: f32 = Standard.sample(rng);
        if drawing < miss_rate {
            continue;
        }
        // TODO: Remove the below %n's, as they give a very slightly
        // non-uniform distribution
        let to_swap_a = (rng.next_u32() as usize) % n;
        let to_swap_b = (rng.next_u32() as usize) % n;
        perm.swap(to_swap_a + 1, to_swap_b + 1);
    }
    perm
}

impl Braid {
    pub fn random_positive(n: usize, num_perms: usize, complexity: usize, miss_rate: f32) -> Braid {
        let mut rng = make_rng();
        let mut result = Braid::from_sigmas(&[], n as BSize);
        
        for _ in 0..num_perms {
            let this_permutation = random_permutation(n, complexity, miss_rate, &mut rng);
            let to_add: Braid = Permutation::from_slice(&this_permutation[..]);
            result = result * to_add;
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn random_tests() {
        let b = Braid::random_positive(60, 30, 3, 0.3);
        println!("{:?}", b);
        println!("{}", b.as_garside_form());
    }
}