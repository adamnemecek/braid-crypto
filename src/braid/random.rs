use rand::{
    distributions::{
        Distribution,
        Standard,
    },
    os::OsRng,
    prng::hc128::*,
    CryptoRng,
    Rng,
    RngCore,
    SeedableRng,
};
// Import Braid and members
use crate::{
    braid::*,
    permutation::*,
};

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
fn random_permutation<CR: CryptoRng + RngCore>(
    n: usize,
    complexity: usize,
    miss_rate: f32,
    rng: &mut CR,
) -> VecPermutation {
    let mut perm = VecPermutation::id(n);
    for _ in 0..complexity {
        let drawing: f32 = Standard.sample(rng);
        if drawing < miss_rate {
            continue;
        }
        let to_swap_a = rng.gen_range(1, n + 1);
        let to_swap_b = rng.gen_range(1, n + 1);
        perm.swap_(to_swap_a, to_swap_b);
    }
    perm
}

impl Braid {
    pub fn random_positive(n: usize, num_perms: usize, complexity: usize, miss_rate: f32) -> Self {
        let mut rng = make_rng();
        let mut result = Self::from_sigmas(&[], n);

        for _ in 0..num_perms {
            let this_permutation = random_permutation(n, complexity, miss_rate, &mut rng);
            result = result * Self::from_slice(&this_permutation[..]);
        }

        result
    }

    pub fn mutate(&mut self, n: usize) {
        let mut rng = make_rng();
        for _ in 0..n {
            match rng.gen_range(1, 4) {
                1 => {
                    let idx = rng.gen_range(0, self.gens.len());
                    let v = rng.gen_range(1, self.n);
                    self.insert_mutation(idx, v);
                }
                2 => {
                    self.swap_mutation();
                }
                3 => {
                    self.exchange_mutation();
                }
                _ => {}
            }
        }
    }

    pub fn insert_mutation(&mut self, idx: usize, v: usize) {
        self.gens.insert(idx, BrGen::Sigma(v));
        self.gens.insert(idx + 1, BrGen::SigmaInv(v));
    }

    pub fn swap_mutation(&mut self) {
        let mut last = self.gens[0];
        for indx in 1..self.gens.len() {
            let curr = self.gens[indx];
            match (last, curr) {
                (BrGen::Sigma(a), BrGen::Sigma(b)) => {
                    if a.abs_diff(b) > 1 {
                        let tmp = self.gens[indx].clone();
                        self.gens[indx] = self.gens[indx - 1];
                        self.gens[indx - 1] = tmp;
                    }
                }
                (BrGen::SigmaInv(a), BrGen::SigmaInv(b)) => {
                    if a.abs_diff(b) > 1 {
                        let tmp = self.gens[indx].clone();
                        self.gens[indx] = self.gens[indx - 1];
                        self.gens[indx - 1] = tmp;
                    }
                }
                _ => {}
            }
            last = self.gens[indx];
        }
    }

    // reidemeister 3
    pub fn exchange_mutation(&mut self) {
        if self.gens.len() < 3 {
            return;
        }

        for idx in 0..self.gens.len() - 2 {
            let kern1 = self.gens[idx];
            let kern2 = self.gens[idx + 1];
            let kern3 = self.gens[idx + 2];

            match (kern1, kern2, kern3) {
                (BrGen::Sigma(a), BrGen::Sigma(b), BrGen::Sigma(c)) => {
                    if a == c && b == a + 1 {
                        self.gens[idx] = BrGen::Sigma(a + 1);
                        self.gens[idx + 1] = BrGen::Sigma(a);
                        self.gens[idx + 2] = BrGen::Sigma(a + 1);
                    }
                }
                (BrGen::SigmaInv(a), BrGen::SigmaInv(b), BrGen::SigmaInv(c)) => {
                    if a == c && b == a + 1 {
                        self.gens[idx] = BrGen::SigmaInv(a + 1);
                        self.gens[idx + 1] = BrGen::SigmaInv(a);
                        self.gens[idx + 2] = BrGen::SigmaInv(a + 1);
                    }
                }
                _ => {}
            }
        }
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

    #[test]
    fn mutation_tests() {
        let mut b = Braid::from_sigmas(&[3, 1, 1, 4, 1, 3, 2, 4], 5);
        b.swap_mutation();
        println!("{:?}", b);

        let mut c = Braid::from_sigmas(&[1, 2, 1, 2, 3, 4, 3], 5);
        c.exchange_mutation();
        println!("{:?}", c);
    }
}
