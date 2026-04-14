//! Combined CLI binary for vLoD - performs detectability analysis and VCF annotation in one step

use clap::Parser;
use env_logger::Env;
use std::path::PathBuf;
use vlod_rs::{
    lod::{calculate_detectability_scores, validate_lod_config},
    merge::merge_detectability_results_into_vcf,
    utils::{get_num_cpus, validate_file_readable, Timer},
    vcf::read_vcf_variants,
    LodConfig, VlodError, VlodResult,
};

#[derive(Parser)]
#[command(name = "vlod")]
#[command(about = "vLoD - Variant Limit of Detection analysis and VCF annotation tool")]
#[command(long_about = "
vLoD (Variant Limit of Detection) analyzes the detectability of variants in a VCF file
by examining corresponding BAM alignment data and directly annotates the VCF with results.

This tool combines detectability analysis and VCF annotation in a single step:
1. Reads variants from the input VCF file
2. Analyzes BAM alignment data to calculate detectability scores
3. Annotates the VCF with detectability information and writes the output

The BAM index file (.bai) must be present next to the BAM file. The tool will
automatically look for files with .bam.bai or .bai extensions.

Two new INFO fields are added to the output VCF:
- DET: Detectability status (Yes if detectable, No if non-detectable)
- DETS: Detectability score (float)

For advanced use cases requiring separate analysis and annotation steps,
use the individual tools: lod_edit and merge_vcf_lod.
")]
struct Args {
    /// Path to the input VCF file
    #[arg(long, value_name = "FILE")]
    input_vcf: PathBuf,

    /// Path to the input BAM file
    #[arg(long, value_name = "FILE")]
    input_bam: PathBuf,

    /// Path to the output annotated VCF file
    #[arg(long, value_name = "FILE")]
    output: PathBuf,

    /// Probability of true positive result
    #[arg(long = "TP", default_value = "0.999")]
    tp: f64,

    /// Probability of false positive result
    #[arg(long = "FP", default_value = "0.001")]
    fp: f64,

    /// Probability of sequencing error
    #[arg(long = "SE", default_value = "0.0001")]
    se: f64,

    /// Number of processes to use for parallel processing
    #[arg(long, default_value_t = get_num_cpus())]
    num_processes: usize,

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

    log::info!("Starting vLoD combined analysis");
    log::info!("Input VCF: {:?}", args.input_vcf);
    log::info!("Input BAM: {:?}", args.input_bam);
    log::info!("Output VCF: {:?}", args.output);
    log::info!("Number of processes: {}", args.num_processes);

    // Validate input files
    validate_file_readable(&args.input_vcf)?;
    validate_file_readable(&args.input_bam)?;

    // Check if output file exists and handle accordingly
    if args.output.exists() && !args.force {
        return Err(VlodError::Io(std::io::Error::new(
            std::io::ErrorKind::AlreadyExists,
            format!("Output file {:?} already exists. Use --force to overwrite.", args.output),
        )));
    }

    // Create output directory if it doesn't exist
    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Create LOD configuration
    let config = LodConfig {
        p_tp: args.tp,
        p_fp: args.fp,
        p_se: args.se,
    };

    // Validate configuration
    validate_lod_config(&config)?;
    log::info!("Configuration: TP={}, FP={}, SE={}", config.p_tp, config.p_fp, config.p_se);

    // Step 1: Read VCF variants
    let _timer = Timer::new("Reading VCF variants");
    let variants = read_vcf_variants(&args.input_vcf)?;
    log::info!("Read {} variants from VCF file", variants.len());

    if variants.is_empty() {
        log::warn!("No variants found in the input VCF file");
        // Copy input VCF to output with detectability headers but no annotations
        std::fs::copy(&args.input_vcf, &args.output)?;
        log::info!("Copied input VCF to output (no variants to analyze)");
        return Ok(());
    }

    // Step 2: Calculate detectability scores
    let _timer = Timer::new("Calculating detectability scores");
    let results = calculate_detectability_scores(
        variants,
        &args.input_bam,
        &config,
        args.num_processes,
    )?;

    log::info!("Calculated detectability scores for {} variants", results.len());

    // Log statistics
    let detectable_count = results.iter().filter(|r| r.detectability_condition == "Detectable").count();
    let non_detectable_count = results.len() - detectable_count;
    
    log::info!("Detectability summary:");
    log::info!("  Detectable: {} ({:.1}%)", detectable_count, (detectable_count as f64 / results.len() as f64) * 100.0);
    log::info!("  Non-detectable: {} ({:.1}%)", non_detectable_count, (non_detectable_count as f64 / results.len() as f64) * 100.0);

    if !results.is_empty() {
        let scores: Vec<f64> = results.iter().map(|r| r.detectability_score).collect();
        let min_score = scores.iter().copied().fold(f64::INFINITY, f64::min);
        let max_score = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let avg_score = scores.iter().sum::<f64>() / scores.len() as f64;
        
        log::info!("  Score range: {:.3} to {:.3}", min_score, max_score);
        log::info!("  Average score: {:.3}", avg_score);
    }

    // Step 3: Merge results directly into VCF
    let _timer = Timer::new("Merging results into VCF");
    merge_detectability_results_into_vcf(&args.input_vcf, &results, &args.output)?;

    log::info!("Analysis completed successfully");
    log::info!("Annotated VCF written to: {:?}", args.output);

    // Log file sizes for reference
    if let Ok(input_size) = std::fs::metadata(&args.input_vcf).map(|m| m.len()) {
        if let Ok(output_size) = std::fs::metadata(&args.output).map(|m| m.len()) {
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
            eprintln!("For BAM files, ensure the index file (.bai) is present.");
        }
        VlodError::InvalidVariant(msg) => {
            eprintln!("Error: Invalid variant data: {}", msg);
            eprintln!("Please check that your VCF file is properly formatted.");
        }
        VlodError::InvalidConfig(msg) => {
            eprintln!("Error: Invalid configuration: {}", msg);
            eprintln!("Please check your probability parameters (TP, FP, SE).");
        }
        VlodError::ThreadPool(msg) => {
            eprintln!("Error: Failed to build thread pool: {}", msg);
        }
        VlodError::Htslib(ref e) => {
            eprintln!("Error: BAM/VCF processing error: {}", e);
            eprintln!("Please check that your BAM file is valid and has an index (.bai) file.");
            eprintln!("Also verify that your VCF file is properly formatted.");
        }
        VlodError::Io(ref e) => {
            eprintln!("Error: I/O error: {}", e);
            eprintln!("Please check file permissions and disk space.");
        }
        VlodError::Csv(ref e) => {
            eprintln!("Error: Data processing error: {}", e);
            eprintln!("This is unexpected in the combined workflow. Please report this issue.");
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
    fn test_combined_workflow_empty_vcf() {
        // Create empty VCF file
        let mut vcf_file = NamedTempFile::new().unwrap();
        writeln!(vcf_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        
        // Create empty BAM file (will fail but that's expected)
        let _bam_file = NamedTempFile::new().unwrap();
        
        let _output_file = NamedTempFile::new().unwrap();
        
        // This should handle empty VCF gracefully
        let _config = LodConfig::default();
        let variants = read_vcf_variants(vcf_file.path()).unwrap();
        assert!(variants.is_empty());
    }

    #[test]
    fn test_config_validation() {
        let config = LodConfig::default();
        assert!(validate_lod_config(&config).is_ok());
        
        let invalid_config = LodConfig {
            p_tp: 0.0,
            p_fp: 0.001,
            p_se: 0.0001,
        };
        assert!(validate_lod_config(&invalid_config).is_err());
    }

    #[test]
    fn test_combined_workflow_integration() {
        use vlod_rs::merge::merge_detectability_results_into_vcf;
        use vlod_rs::DetectabilityResult;
        use vlod_rs::Variant;

        // Create test VCF with variants
        let mut vcf_file = NamedTempFile::new().unwrap();
        writeln!(vcf_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(vcf_file, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
        writeln!(vcf_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(vcf_file, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30").unwrap();
        writeln!(vcf_file, "chr2\t200\t.\tG\tC\t.\tPASS\tDP=40").unwrap();
        
        // Read variants from VCF
        let variants = read_vcf_variants(vcf_file.path()).unwrap();
        assert_eq!(variants.len(), 2);
        
        // Create mock detectability results
        let results = vec![
            DetectabilityResult::new(
                Variant::new("chr1".to_string(), 100, "A".to_string(), "T".to_string()),
                3.5,
                "Detectable".to_string(),
                30,
                15,
            ),
            DetectabilityResult::new(
                Variant::new("chr2".to_string(), 200, "G".to_string(), "C".to_string()),
                1.2,
                "Non-detectable".to_string(),
                40,
                8,
            ),
        ];
        
        // Test the merge functionality (core of the combined workflow)
        let output_file = NamedTempFile::new().unwrap();
        let merge_result = merge_detectability_results_into_vcf(
            vcf_file.path(),
            &results,
            output_file.path(),
        );
        
        assert!(merge_result.is_ok());
        
        // Verify the output contains the expected annotations
        let output_content = std::fs::read_to_string(output_file.path()).unwrap();
        assert!(output_content.contains("DET=Yes"));
        assert!(output_content.contains("DETS=3.5"));
        assert!(output_content.contains("DET=No"));
        assert!(output_content.contains("DETS=1.2"));
        assert!(output_content.contains("##INFO=<ID=DET,Number=1,Type=String"));
        assert!(output_content.contains("##INFO=<ID=DETS,Number=1,Type=Float"));
    }
}