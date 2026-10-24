use std::path::Path;

fn main() {
    let root = Path::new("EPANET");
    let src = root.join("src");
    if !src.join("epanet.c").exists() {
        panic!(
            "EPANET sources not found in {}. Run `git submodule update --init` first.",
            root.display()
        );
    }

    let mut build = cc::Build::new();
    for dir in [src.clone(), src.join("util")] {
        for entry in std::fs::read_dir(&dir).expect("read EPANET source dir") {
            let path = entry.unwrap().path();
            if path.extension().is_some_and(|e| e == "c") && path.file_name().unwrap() != "epanet_py.c" {
                println!("cargo:rerun-if-changed={}", path.display());
                build.file(path);
            }
        }
    }
    build
        .include(root.join("include"))
        .include(&src)
        // Static linking: no dllimport/dllexport decoration, default calling convention.
        .define("DLLEXPORT", Some(""))
        .define("_CRT_SECURE_NO_WARNINGS", None)
        .warnings(false)
        .compile("epanet2");
}
