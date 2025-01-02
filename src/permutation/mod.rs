pub type VecPermutation = Vec<usize>;

pub trait Permutation {
    fn id(n: usize) -> Self;
    fn size(&self) -> usize;
    fn swap(&mut self, a: usize, b: usize);
    // Where does the "strand" starting at position x end up?
    fn position(&self, x: usize) -> usize;
    // Where did the "strand" that ended up at position x start from?
    fn follow_ending(&self, x: usize) -> usize {
        self.as_vec()[x - 1]
    }

    fn from_slice(v: &[usize]) -> Self;

    fn as_vec(&self) -> Vec<usize> {
        let n = self.size();
        (1..=n).map(|i| self.follow_ending(i)).collect()
    }

    fn is_identity(&self) -> bool {
        (1..=self.size()).all(|i| self.position(i) == i)
    }

    fn is_twist(&self) -> bool {
        let n = self.size();
        (1..=n).all(|i| self.position(i) == n - i + 1)
    }

    fn compose<B: Permutation, C: Permutation>(&self, second: &B) -> C {
        debug_assert_eq!(self.size(), second.size());
        let n = self.size();
        let mut res = VecPermutation::id(n);
        for strand_number in 1..=n {
            let after_first = self.position(strand_number);
            println!("{:?}", after_first);
            let result_place = second.position(after_first);
            res[result_place - 1] = strand_number;
        }
        C::from_slice(&res[..])
    }
}

// TODO: More tests on this function's correctness
pub fn from_slice_slow<P: Permutation>(v: &[usize]) -> P {
    let n = v.len();
    // Default implimentation using swap
    let mut res = P::id(n);
    let mut reference = VecPermutation::id(n);
    for (i, target) in v.iter().enumerate().take(n) {
        let place_to = reference
            .iter()
            .position(|&pos| pos == *target)
            .expect("A non-permutation slice was passed to from_slice");
        res.swap(i, place_to);
        reference.swap(i, place_to);
    }
    res
}

impl Permutation for VecPermutation {
    fn id(n: usize) -> Self {
        (1..=n).collect()
    }

    fn size(&self) -> usize {
        self.len()
    }

    fn swap(&mut self, a: usize, b: usize) {
        let at_a = self[a - 1];
        let at_b = self[b - 1];
        self[b - 1] = at_a;
        self[a - 1] = at_b;
    }

    fn position(&self, x: usize) -> usize {
        self.iter()
            .position(|&strand_number| strand_number == x)
            .expect("Invalid x given to position. X has to be in the range (1..=n)")
            + 1
    }

    fn follow_ending(&self, x: usize) -> usize {
        self[x - 1]
    }

    fn from_slice(v: &[usize]) -> Self {
        v.to_vec()
    }

    fn as_vec(&self) -> Vec<usize> {
        self.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_module() {
        let perm1: VecPermutation = vec![4, 2, 3, 1];
        let mut perm2 = VecPermutation::id(4);
        perm2.swap(2, 3);
        let composed: VecPermutation = perm1.compose(&perm2);
        assert_eq!(composed, vec![4, 3, 2, 1]);
    }
}
