#![allow(unknown_lints)]
pub mod braid;
pub mod permutation;

#[macro_use]
extern crate serde_derive;

extern crate serde;
extern crate bincode;
extern crate indexmap;
extern crate rand;

use braid::*;

fn main() {
    let braid1 = Braid::from_sigmas(&[1, -3, 2], 4);
    let braid2 = Braid::from_sigmas(&[2, 3, 1], 4);
    let sum = braid1 * braid2;

    let gform = sum.as_garside_form();
    println!("{}", gform);
    // outputs: [1;(4, 3, 1, 2), (2, 3, 1, 4), (1, 3, 4, 2)]
}
