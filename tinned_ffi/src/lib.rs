#![deny(unsafe_op_in_unsafe_fn)]

#[macro_use]
mod macros;

mod c_support;
mod core;
mod expressions;
mod perturbations;
mod public;

// Re-export `ExprTag`` for convenience
pub use tinned::ExprTag;

pub use crate::c_support::*;
pub use crate::core::*;
pub use crate::perturbations::*;
pub use crate::public::*;

#[cfg(all(test, feature = "c-headers"))]
mod header_gen {
    use safer_ffi::headers::{Language, NamingConvention, builder};
    use std::{fs, io, path::PathBuf};

    // Generates include/tinned.h
    #[test]
    #[safer_ffi::cfg_headers]
    fn generate_c_header() -> io::Result<()> {
        let out_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("include");
        fs::create_dir_all(&out_dir)?;
        let out = out_dir.join("tinned.h");

        builder()
            .with_language(Language::C)
            .with_stable_header(true)
            .with_guard("__TINNED_H__")
            .with_naming_convention(NamingConvention::Prefix("tinned_".into()))
            .to_file(&out)?
            .generate()?;

        eprintln!("Generated header at {}", out.display());
        Ok(())
    }
}
