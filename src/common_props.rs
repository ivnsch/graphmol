use anyhow::{anyhow, Result};
use once_cell::sync::OnceCell;
use std::{
    collections::HashMap,
    sync::{Mutex, MutexGuard},
};

#[derive(Clone, Debug)]
pub enum CommonPropsValue {
    U32(u32),
    String(String),
}

#[derive(PartialEq, Eq, Hash, Debug)]
pub enum CommonPropsKey {
    ChiralPermutation,
}

static COMMON_PROPS: OnceCell<Mutex<HashMap<CommonPropsKey, CommonPropsValue>>> = OnceCell::new();

pub fn common_props<'a>() -> Result<MutexGuard<'a, HashMap<CommonPropsKey, CommonPropsValue>>> {
    COMMON_PROPS
        .get_or_init(|| Mutex::new(HashMap::new()))
        .lock()
        .map_err(|err| anyhow!("Error accessing mutex: {:?}", err))
}
