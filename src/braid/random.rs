
use rand::Rng;
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
        let to_swap_a = rng.gen_range(1, n+1);
        let to_swap_b = rng.gen_range(1, n+1);
        perm.swap(to_swap_a, to_swap_b);
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

    pub fn mutate(&mut self, n: usize) {
        let mut rng = make_rng();
        for _ in 0..n {
            match rng.gen_range(1, 4) {
                1 => {
                    self.insert_mutation(&mut rng);
                },
                2 => {
                    self.swap_mutation();
                },
                3 => {
                    self.exchange_mutation();
                },
                _ => {}
            }
        }
    }

    pub fn insert_mutation<CR: CryptoRng + RngCore>(&mut self, rng: &mut CR) {
        let insertion_point = rng.gen_range(0, self.contents.len());
        let to_insert = rng.gen_range(1, self.n);
        self.contents.insert(insertion_point, Sigma(to_insert));
        self.contents.insert(insertion_point + 1, SigmaInv(to_insert));
    }

    pub fn swap_mutation(&mut self) {
        let mut last = self.contents[0];
        for indx in 1..self.contents.len() {
            let curr = self.contents[indx];
            match (last, curr) {
                (Sigma(a), Sigma(b)) => {
                    if (a as isize - b as isize).abs() > 1 {
                        let tmp = self.contents[indx].clone();
                        self.contents[indx] = self.contents[indx - 1];
                        self.contents[indx - 1] = tmp;
                    }
                },
                (SigmaInv(a), SigmaInv(b)) => {
                    if (a as isize - b as isize).abs() > 1 {
                        let tmp = self.contents[indx].clone();
                        self.contents[indx] = self.contents[indx - 1];
                        self.contents[indx - 1] = tmp;
                    }
                }
                _ => {}
            }
            last = self.contents[indx];
        }
    }

    pub fn exchange_mutation(&mut self) {
        if self.contents.len() < 3 {
            return;
        }
        let mut indx = 0;
        while indx < self.contents.len() - 2 {
            let kern1 = self.contents[indx];
            let kern2 = self.contents[indx + 1];
            let kern3 = self.contents[indx + 2];

            match (kern1, kern2, kern3) {
                (Sigma(a), Sigma(b), Sigma(c)) => {
                    if a == c && b == a + 1 {
                        self.contents[indx] = Sigma(a + 1);
                        self.contents[indx + 1] = Sigma(a);
                        self.contents[indx + 2] = Sigma(a + 1);
                    }
                },
                (SigmaInv(a), SigmaInv(b), SigmaInv(c)) => {
                    if a == c && b == a + 1 {
                        self.contents[indx] = SigmaInv(a + 1);
                        self.contents[indx + 1] = SigmaInv(a);
                        self.contents[indx + 2] = SigmaInv(a + 1);
                    }
                },
                _ => {},
            }
            indx += 1;
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