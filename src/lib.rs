#![allow(unknown_lints)]

#[macro_use]
extern crate serde_derive;

pub mod braid;
mod permutation;

pub mod prelude {
    pub use crate::{braid::*, permutation::*};
}
