#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BrGen {
    Sigma(usize),
    SigmaInv(usize)
}

#[derive(Debug, Clone)]
pub struct Braid {
    pub contents: Vec<BrGen>,
    pub n: usize // Our braid is an element of B_n
}