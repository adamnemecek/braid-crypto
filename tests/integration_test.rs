extern crate braid_crypto;

use braid_crypto::braid::*;

#[test]
fn test_garsides_mutations() {
    let b = Braid::random_positive(20, 10, 3, 0.2);
    let mut c = b.clone();
    c.mutate(50);

    let s1 = format!("{}", b.as_garside_form());
    let s2 = format!("{}", c.as_garside_form());
    assert_eq!(s1, s2);
}