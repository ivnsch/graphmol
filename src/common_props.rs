#[derive(PartialEq, Eq, Clone, Debug)]
pub enum PropsValue {
    U32(u32),
    String(String),
}

#[derive(PartialEq, Eq, Hash, Debug)]
pub enum PropsKey {
    ChiralPermutation,
}
