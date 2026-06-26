use std::process::Command;

#[test]
fn creates_output_file() {
    let status = Command::new(env!("CARGO_BIN_EXE_rustyomestats"))
        .args([
            "genome",
            "--fasta", "tests/data/JAPAEP01-1.fna",
            "--outdir", "tests/output",
        ])
        .status()
        .unwrap();

    assert!(status.success());
    assert!(std::path::Path::new("tests/output/report.html").exists());
}
