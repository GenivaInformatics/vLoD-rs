//! CLI binary for VCF integration - equivalent to merge_vcf_lod.py

use clap::Parser;
use env_logger::Env;
use std::path::PathBuf;
use vlod_rs::{
    merge::merge_detectability_into_vcf,
    utils::{validate_file_readable, Timer},
    VlodError, VlodResult,
};

#[derive(Parser)]
#[command(name = "merge_vcf_lod")]
#[command(about = "Merge detectability results into VCF files")]
#[command(long_about = "
This tool merges detectability analysis results back into the original VCF file.
It reads detectability scores from a TSV file (produced by lod_edit) and adds
them as INFO fields to the corresponding variants in the VCF file.

Two new INFO fields are added:
- DET: Detectability status (Yes/No)
- DETS: Detectability score (float)

The tool supports both compressed and uncompressed VCF files.
")]
struct Args {
    /// Path to the input VCF file
    #[arg(value_name = "VCF_FILE")]
    vcf_file: PathBuf,

    /// Path to the detectability results TSV file
    #[arg(value_name = "DETECTABILITY_FILE")]
    detectability_file: PathBuf,

    /// Path to the output VCF file
    #[arg(value_name = "OUTPUT_FILE")]
    output_file: PathBuf,

    /// Enable verbose logging
    #[arg(short, long)]
    verbose: bool,

    /// Enable debug logging
    #[arg(short, long)]
    debug: bool,

    /// Force overwrite of output file if it exists
    #[arg(short, long)]
    force: bool,
}

fn run() -> VlodResult<()> {
    let args = Args::parse();

    // Initialize logging
    let log_level = if args.debug {
        "debug"
    } else if args.verbose {
        "info"
    } else {
        "warn"
    };

    env_logger::Builder::from_env(Env::default().default_filter_or(log_level))
        .format_timestamp_secs()
        .init();

    log::info!("Starting VCF merge operation");
    log::info!("VCF file: {:?}", args.vcf_file);
    log::info!("Detectability file: {:?}", args.detectability_file);
    log::info!("Output file: {:?}", args.output_file);

    // Validate input files
    validate_file_readable(&args.vcf_file)?;
    validate_file_readable(&args.detectability_file)?;

    // Check if output file exists and handle accordingly
    if args.output_file.exists() && !args.force {
        return Err(VlodError::Io(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("Output file {:?} already exists. Use --force to overwrite.", args.output_file),
        )));
    }

    // Create output directory if it doesn't exist
    if let Some(parent) = args.output_file.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Perform the merge operation
    let _timer = Timer::new("Merging detectability results into VCF");
    merge_detectability_into_vcf(&args.vcf_file, &args.detectability_file, &args.output_file)?;

    log::info!("Merge operation completed successfully");
    log::info!("Output written to: {:?}", args.output_file);

    // Log file sizes for reference
    if let Ok(input_size) = std::fs::metadata(&args.vcf_file).map(|m| m.len()) {
        if let Ok(output_size) = std::fs::metadata(&args.output_file).map(|m| m.len()) {
            log::info!("Input VCF size: {} bytes", input_size);
            log::info!("Output VCF size: {} bytes", output_size);
            
            if output_size > input_size {
                let size_increase = output_size - input_size;
                log::info!("Size increase: {} bytes ({:.1}%)", 
                          size_increase, 
                          (size_increase as f64 / input_size as f64) * 100.0);
            }
        }
    }

    Ok(())
}

/// Handle application errors and provide user-friendly messages
fn handle_error(error: VlodError) -> ! {
    match error {
        VlodError::FileNotFound(path) => {
            eprintln!("Error: File not found: {}", path);
            eprintln!("Please check that the file exists and is readable.");
        }
        VlodError::InvalidVariant(msg) => {
            eprintln!("Error: Invalid variant data: {}", msg);
            eprintln!("Please check that your VCF or detectability file is properly formatted.");
        }
        VlodError::InvalidConfig(msg) => {
            eprintln!("Error: Invalid configuration: {}", msg);
        }
        VlodError::ThreadPool(msg) => {
            eprintln!("Error: Failed to build thread pool: {}", msg);
        }
        VlodError::Htslib(ref e) => {
            eprintln!("Error: VCF processing error: {}", e);
            eprintln!("Please check that your VCF file is valid.");
        }
        VlodError::Io(ref e) => {
            eprintln!("Error: I/O error: {}", e);
            eprintln!("Please check file permissions and disk space.");
        }
        VlodError::Csv(ref e) => {
            eprintln!("Error: CSV processing error: {}", e);
            eprintln!("Please check the detectability file format.");
        }
    }
    std::process::exit(1);
}

fn main() {
    if let Err(e) = run() {
        handle_error(e);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_merge_vcf_integration() {
        // Create test detectability file
        let mut detectability_file = NamedTempFile::new().unwrap();
        writeln!(detectability_file, "Chrom\tPos\tRef\tAlt\tDetectability_Score\tDetectability_Condition\tCoverage\tVariant_Reads").unwrap();
        writeln!(detectability_file, "chr1\t100\tA\tT\t3.5\tDetectable\t30\t15").unwrap();
        
        // Create test VCF file
        let mut vcf_file = NamedTempFile::new().unwrap();
        writeln!(vcf_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf_file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(vcf_file, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30").unwrap();
        
        let output_file = NamedTempFile::new().unwrap();
        
        // Test the merge operation
        let result = merge_detectability_into_vcf(
            vcf_file.path(),
            detectability_file.path(),
            output_file.path(),
        );
        
        assert!(result.is_ok());
        
        // Verify the output contains the expected modifications
        let output_content = std::fs::read_to_string(output_file.path()).unwrap();
        assert!(output_content.contains("DET=Yes"));
        assert!(output_content.contains("DETS=3.5"));
        assert!(output_content.contains("##INFO=<ID=DET,Number=1,Type=String"));
        assert!(output_content.contains("##INFO=<ID=DETS,Number=1,Type=Float"));
    }
}