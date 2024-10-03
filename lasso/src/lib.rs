#![allow(non_snake_case)]
#![allow(clippy::needless_range_loop)]
#![feature(generic_arg_infer)]

pub mod lasso;
pub mod memory_checking;
pub mod table;

pub use lasso::{LassoNode, LassoPreprocessing};
