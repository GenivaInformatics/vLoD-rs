//! # vLoD - Variant Limit of Detection Tool
//!
//! A Rust implementation of the vLoD tool for assessing the detectability status
//! of alleles from variant call files (VCF) using matched sequencing data.

pub mod bam;
pub mod lod;
pub mod merge;
pub mod utils;
pub mod vcf;

use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Represents a genomic variant with its position and alleles
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Variant {
    pub chrom: String,
    pub pos: u32,
    pub ref_allele: String,
    pub alt_allele: String,
}

impl Variant {
    pub fn new(chrom: String, pos: u32, ref_allele: String, alt_allele: String) -> Self {
        Self {
            chrom,
            pos,
            ref_allele,
            alt_allele,
        }
    }
}

/// Represents the detectability analysis result for a variant
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectabilityResult {
    pub variant: Variant,
    pub detectability_score: f64,
    pub detectability_condition: String,
    pub coverage: u32,
    pub variant_reads: u32,
}

impl DetectabilityResult {
    pub fn new(
        variant: Variant,
        detectability_score: f64,
        detectability_condition: String,
        coverage: u32,
        variant_reads: u32,
    ) -> Self {
        Self {
            variant,
            detectability_score,
            detectability_condition,
            coverage,
            variant_reads,
        }
    }

    /// Determine detectability condition based on score
    pub fn condition_from_score(score: f64) -> String {
        if score >= 2.50 {
            "Detectable".to_string()
        } else {
            "Non-detectable".to_string()
        }
    }
}

/// Configuration parameters for LOD calculation
#[derive(Debug, Clone)]
pub struct LodConfig {
    pub p_tp: f64,  // Probability of true positive
    pub p_fp: f64,  // Probability of false positive
    pub p_se: f64,  // Probability of sequencing error
}

impl Default for LodConfig {
    fn default() -> Self {
        Self {
            p_tp: 0.999,
            p_fp: 0.001,
            p_se: 0.0001,
        }
    }
}

/// Error types for the vLoD library
#[derive(Debug, thiserror::Error)]
pub enum VlodError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("HTSlib error: {0}")]
    Htslib(#[from] rust_htslib::errors::Error),
    
    #[error("CSV error: {0}")]
    Csv(#[from] csv::Error),
    
    #[error("Invalid variant format: {0}")]
    InvalidVariant(String),
    
    #[error("File not found: {0}")]
    FileNotFound(String),
    
    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),

    #[error("Thread pool error: {0}")]
    ThreadPool(String),
}

pub type VlodResult<T> = Result<T, VlodError>;