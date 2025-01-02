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

#[test]
fn test_key_exchange() {
    for i in 0..2 {
        println!("trial {}", i);
        println!("computing random braids...");
        let public = Braid::random_positive(15, 5, 2, 0.0);

        let mut s_alice = Braid::random_positive(7, 3, 2, 0.1);
        let mut r_bob = Braid::random_positive(7, 3, 2, 0.1);

        s_alice.n = 15;
        r_bob.n = 15;

        r_bob.shift();

        println!("computing public messages...");
        let p_prime = s_alice.clone() * public.clone() * s_alice.inverse();
        let p_prime_prime = r_bob.clone() * public.clone() * r_bob.inverse();

        println!("computing shared keys...");
        let k_alice = s_alice.clone() * p_prime_prime * s_alice.inverse();
        let k_bob = r_bob.clone() * p_prime * r_bob.inverse();

        println!("computing reduced form k1...");
        let k1 = k_alice.as_garside_form();
        println!("computing reduced form k2...");
        let k2 = k_bob.as_garside_form();

        println!("results:\n{}\n{}", k1, k2);

        assert_eq!(k1.to_string(), k2.to_string());
    }
}

#[test]
fn test_slides() {
    let braid1 = Braid::from_sigmas(&[1, -3, 2], 4);
    let braid2 = Braid::from_sigmas(&[2, 3, 1], 4);
    let sum = braid1 * braid2;

    let gform = sum.as_garside_form();
    println!("{}", gform);
    // outputs: [-1;(4, 3, 1, 2), (2, 3, 1, 4), (3, 1, 4, 2)]
}
