//! CLI binary for LOD analysis - equivalent to LOD_edit.py

use clap::Parser;
use env_logger::Env;
use std::path::PathBuf;
use vlod_rs::{
    lod::{calculate_detectability_scores, validate_lod_config, write_detectability_results},
    utils::{get_num_cpus, validate_file_readable, Timer},
    vcf::read_vcf_variants,
    LodConfig, VlodError, VlodResult,
};

#[derive(Parser)]
#[command(name = "lod_edit")]
#[command(about = "Detectability analysis tool for VCF variants using BAM alignment data")]
#[command(long_about = "
This tool analyzes the detectability of variants in a VCF file by examining
the corresponding BAM alignment data. It calculates a detectability score
for each variant based on variant allele frequency (VAF) and statistical
parameters including true positive rate, false positive rate, and sequencing
error rate.

The BAM index file (.bai) must be present next to the BAM file. The tool will
automatically look for files with .bam.bai or .bai extensions.

The output is a TSV file containing detectability scores and classifications
for each variant, along with coverage and read count information.
")]
struct Args {
    /// Path to the input VCF file
    #[arg(long, value_name = "FILE")]
    input_vcf: PathBuf,

    /// Path to the input BAM file
    #[arg(long, value_name = "FILE")]
    input_bam: PathBuf,

    /// Path to the output TSV file
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

    log::info!("Starting vLoD analysis");
    log::info!("VCF file: {:?}", args.input_vcf);
    log::info!("BAM file: {:?}", args.input_bam);
    log::info!("Output file: {:?}", args.output);
    log::info!("Number of processes: {}", args.num_processes);

    // Validate input files
    validate_file_readable(&args.input_vcf)?;
    validate_file_readable(&args.input_bam)?;

    // Create LOD configuration
    let config = LodConfig {
        p_tp: args.tp,
        p_fp: args.fp,
        p_se: args.se,
    };

    // Validate configuration
    validate_lod_config(&config)?;

    log::info!(
        "Configuration: TP={}, FP={}, SE={}",
        config.p_tp,
        config.p_fp,
        config.p_se
    );

    // Create output directory if it doesn't exist
    if let Some(parent) = args.output.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Read VCF variants
    let _timer = Timer::new("Reading VCF variants");
    let variants = read_vcf_variants(&args.input_vcf)?;
    log::info!("Read {} variants from VCF file", variants.len());

    if variants.is_empty() {
        log::warn!("No variants found in the input VCF file");
        // Create empty output file with header
        write_detectability_results(&[], &args.output)?;
        return Ok(());
    }

    // Calculate detectability scores
    let _timer = Timer::new("Calculating detectability scores");
    let results =
        calculate_detectability_scores(variants, &args.input_bam, &config, args.num_processes)?;

    log::info!(
        "Calculated detectability scores for {} variants",
        results.len()
    );

    // Log statistics
    let detectable_count = results
        .iter()
        .filter(|r| r.detectability_condition == "Detectable")
        .count();
    let non_detectable_count = results.len() - detectable_count;

    log::info!("Results summary:");
    log::info!(
        "  Detectable: {} ({:.1}%)",
        detectable_count,
        (detectable_count as f64 / results.len() as f64) * 100.0
    );
    log::info!(
        "  Non-detectable: {} ({:.1}%)",
        non_detectable_count,
        (non_detectable_count as f64 / results.len() as f64) * 100.0
    );

    if !results.is_empty() {
        let scores: Vec<f64> = results.iter().map(|r| r.detectability_score).collect();
        let min_score = scores.iter().copied().fold(f64::INFINITY, f64::min);
        let max_score = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let avg_score = scores.iter().sum::<f64>() / scores.len() as f64;

        log::info!("  Score range: {:.3} to {:.3}", min_score, max_score);
        log::info!("  Average score: {:.3}", avg_score);
    }

    // Write results
    let _timer = Timer::new("Writing results");
    write_detectability_results(&results, &args.output)?;

    log::info!("Results written to: {:?}", args.output);
    log::info!("Analysis completed successfully");

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
        }
        VlodError::Io(ref e) => {
            eprintln!("Error: I/O error: {}", e);
            eprintln!("Please check file permissions and disk space.");
        }
        VlodError::Csv(ref e) => {
            eprintln!("Error: CSV processing error: {}", e);
            eprintln!("Please check the output file format.");
        }
    }
    std::process::exit(1);
}

fn main() {
    if let Err(e) = run() {
        handle_error(e);
    }
}
