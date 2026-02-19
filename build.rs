fn main() {
    cc::Build::new()
        .file("primer3/src/thal.c")
        .include("primer3/src")
        .compile("thal");

    println!("cargo:rerun-if-changed=primer3/src/thal.c");
    println!("cargo:rerun-if-changed=primer3/src/thal.h");
    println!("cargo:rerun-if-changed=primer3/src/thal_default_params.h");
}
