#![allow(unknown_lints)]

#[macro_use]
extern crate serde_derive;

pub mod braid;
pub mod permutation;

use braid_crypto::prelude::*;

fn main() {
    println!("running Diffie-Hellman-type integration test");
    println!("computing random braids...");
    let public = Braid::from_sigmas(&[1, 2, 3, 4, 5, 6, 7], 8);

    let mut s_alice = Braid::random_positive(3, 3, 3, 0.1);
    let mut r_bob = Braid::random_positive(3, 3, 3, 0.1);

    s_alice.n = 8;
    r_bob.n = 8;

    r_bob.shift();

    println!("computing public messages...");
    let p_prime = s_alice.clone() * public.clone() * s_alice.inverse();
    let p_prime_prime = r_bob.clone() * public.clone() * r_bob.inverse();

    println!("computing shared keys...");
    let k_alice = s_alice.clone() * p_prime_prime.clone() * s_alice.inverse();
    let k_bob = r_bob.clone() * p_prime.clone() * r_bob.inverse();

    println!("k_alice is length: {}", k_alice.gens.len());
    println!("s_alice is length: {}", s_alice.gens.len());

    println!("computing reduced form k1...");
    let k1 = k_alice.as_garside_form();
    println!("computing reduced form k2...");
    let k2 = k_bob.as_garside_form();

    println!("results:");
    println!("{}", k1);
    println!();
    println!("{}", k2);

    println!("All Data:");
    println!("{:?}\n\n{:?}\n\n{:?}", public, k_alice, k_bob);
}
