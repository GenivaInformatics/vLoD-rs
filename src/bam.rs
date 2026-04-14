//! BAM file processing and pileup analysis

use crate::{LodConfig, Variant, VlodError, VlodResult};
use rust_htslib::bam::{pileup::Alignment, IndexedReader, Read};
use std::collections::HashMap;
use std::path::Path;

/// Represents allele counts at a specific position
#[derive(Debug, Clone)]
pub struct AlleleCounts {
    pub ref_count: u32,
    pub alt_counts: HashMap<String, u32>,
    pub total_count: u32,
}

impl AlleleCounts {
    pub fn new() -> Self {
        Self {
            ref_count: 0,
            alt_counts: HashMap::new(),
            total_count: 0,
        }
    }

    pub fn add_ref(&mut self) {
        self.ref_count += 1;
        self.total_count += 1;
    }

    pub fn add_alt(&mut self, allele: String) {
        *self.alt_counts.entry(allele).or_insert(0) += 1;
        self.total_count += 1;
    }

    pub fn get_alt_count(&self, allele: &str) -> u32 {
        self.alt_counts.get(allele).copied().unwrap_or(0)
    }

    pub fn get_vaf(&self, allele: &str) -> f64 {
        if self.total_count == 0 {
            0.0
        } else {
            self.get_alt_count(allele) as f64 / self.total_count as f64
        }
    }
}

/// BAM analyzer for processing variants
pub struct BamAnalyzer {
    bam_reader: IndexedReader,
}

impl BamAnalyzer {
    pub fn new<P: AsRef<Path>>(bam_path: P) -> VlodResult<Self> {
        let bam_path = bam_path.as_ref();
        
        // Check for BAI index file next to the BAM file
        let bai_path = bam_path.with_extension("bam.bai");
        let alt_bai_path = bam_path.with_extension("bai");
        
        let bam_reader = if bai_path.exists() {
            IndexedReader::from_path_and_index(bam_path, &bai_path)?
        } else if alt_bai_path.exists() {
            IndexedReader::from_path_and_index(bam_path, &alt_bai_path)?
        } else {
            return Err(VlodError::FileNotFound(format!(
                "BAM index file not found. Expected {} or {}",
                bai_path.display(),
                alt_bai_path.display()
            )));
        };
        
        Ok(BamAnalyzer { bam_reader })
    }

    /// Analyze a single variant and return allele counts
    pub fn analyze_variant(&mut self, variant: &Variant) -> VlodResult<AlleleCounts> {
        let tid = self.bam_reader.header().tid(variant.chrom.as_bytes())
            .ok_or_else(|| VlodError::InvalidVariant(format!("Unknown chromosome: {}", variant.chrom)))?;

        // Fetch only the specific region around the variant
        // For indels, we need a slightly larger window
        let ref_len = variant.ref_allele.len();
        let alt_lens: Vec<usize> = variant.alt_allele.split(',').map(|a| a.len()).collect();
        let max_len = (*alt_lens.iter().max().unwrap_or(&1)).max(ref_len) as u32;
        
        // Fetch region with some padding for indels
        let start = variant.pos.saturating_sub(1); // Convert to 0-based
        let end = variant.pos.saturating_add(max_len); // Inclusive end
        
        self.bam_reader.fetch((tid, start, end))?;

        let mut pileup = self.bam_reader.pileup();
        pileup.set_max_depth(1_000_000);

        let mut allele_counts = AlleleCounts::new();
        let alt_alleles: Vec<&str> = variant.alt_allele.split(',').collect();

        for p in pileup {
            let p = p?;
            
            // Check if this is the position we're interested in
            if p.pos() as u32 != variant.pos - 1 {
                continue;
            }

            for alignment in p.alignments() {
                if alignment.is_refskip() {
                    continue;
                }

                let ref_len = variant.ref_allele.len();
                let alt_len = alt_alleles.iter().map(|a| a.len()).max().unwrap_or(0);

                if ref_len == alt_len {
                    // SNV or MNV
                    Self::process_snv_mnv(&alignment, variant, &alt_alleles, &mut allele_counts)?;
                } else {
                    // Indel
                    Self::process_indel(&alignment, variant, &alt_alleles, &mut allele_counts)?;
                }
            }
            
            // Since we fetched a specific region and found our position, we can break
            break;
        }

        Ok(allele_counts)
    }

    fn process_snv_mnv(
        alignment: &Alignment,
        variant: &Variant,
        alt_alleles: &[&str],
        allele_counts: &mut AlleleCounts,
    ) -> VlodResult<()> {
        if alignment.is_del() {
            return Ok(());
        }

        let qpos = alignment.qpos();
        if qpos.is_none() {
            return Ok(());
        }

        let qpos = qpos.unwrap();
        let record = alignment.record();
        let seq = record.seq();
        let ref_len = variant.ref_allele.len();

        if ref_len == 1 {
            // SNV
            if qpos < seq.len() {
                let base = seq[qpos] as char;
                let base_str = base.to_string();
                
                if base_str == variant.ref_allele {
                    allele_counts.add_ref();
                } else if alt_alleles.contains(&base_str.as_str()) {
                    allele_counts.add_alt(base_str);
                }
            }
        } else {
            // MNV
            if qpos + ref_len <= seq.len() {
                let read_seq: String = (qpos..qpos + ref_len)
                    .map(|i| seq[i] as char)
                    .collect();
                
                if read_seq == variant.ref_allele {
                    allele_counts.add_ref();
                } else if alt_alleles.contains(&read_seq.as_str()) {
                    allele_counts.add_alt(read_seq);
                }
            }
        }

        Ok(())
    }

    fn process_indel(
        alignment: &Alignment,
        variant: &Variant,
        alt_alleles: &[&str],
        allele_counts: &mut AlleleCounts,
    ) -> VlodResult<()> {
        use rust_htslib::bam::pileup::Indel;
        
        let indel = alignment.indel();
        
        for &alt_allele in alt_alleles {
            let expected_indel = alt_allele.len() as i32 - variant.ref_allele.len() as i32;
            
            match indel {
                Indel::Ins(n) if expected_indel > 0 && n == expected_indel as u32 => {
                    allele_counts.add_alt(alt_allele.to_string());
                }
                Indel::Del(n) if expected_indel < 0 && n == expected_indel.abs() as u32 => {
                    allele_counts.add_alt(alt_allele.to_string());
                }
                Indel::None => {
                    allele_counts.add_ref();
                }
                _ => {}
            }
        }

        Ok(())
    }
}

/// Process a single variant using an existing analyzer
pub fn process_single_variant(
    analyzer: &mut BamAnalyzer,
    variant: &Variant,
    config: &LodConfig,
) -> VlodResult<Vec<(Variant, f64, u32, u32)>> {
    let allele_counts = analyzer.analyze_variant(variant)?;
    let mut results = Vec::new();

    for alt_allele in variant.alt_allele.split(',') {
        let alt_count = allele_counts.get_alt_count(alt_allele);
        let vaf = allele_counts.get_vaf(alt_allele);

        let lod_value = (config.p_tp * vaf) / ((1.0 - vaf) * config.p_se + vaf * config.p_fp);
        let lod = if lod_value > 0.0 { lod_value.log10() } else { f64::NEG_INFINITY };

        results.push((
            Variant::new(
                variant.chrom.clone(),
                variant.pos,
                variant.ref_allele.clone(),
                alt_allele.to_string(),
            ),
            lod,
            allele_counts.total_count,
            alt_count,
        ));
    }

    Ok(results)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;
    use std::fs::File;

    #[test]
    fn test_allele_counts() {
        let mut counts = AlleleCounts::new();
        
        counts.add_ref();
        counts.add_ref();
        counts.add_alt("T".to_string());
        
        assert_eq!(counts.ref_count, 2);
        assert_eq!(counts.get_alt_count("T"), 1);
        assert_eq!(counts.total_count, 3);
        assert_eq!(counts.get_vaf("T"), 1.0 / 3.0);
    }

    #[test]
    fn test_vaf_calculation() {
        let mut counts = AlleleCounts::new();
        counts.add_ref();
        counts.add_alt("G".to_string());
        counts.add_alt("G".to_string());
        
        assert_eq!(counts.get_vaf("G"), 2.0 / 3.0);
        assert_eq!(counts.get_vaf("T"), 0.0);
    }

    #[test]
    fn test_empty_allele_counts() {
        let counts = AlleleCounts::new();
        assert_eq!(counts.get_vaf("A"), 0.0);
        assert_eq!(counts.total_count, 0);
    }

    #[test]
    fn test_bam_analyzer_index_detection() {
        // Test with missing BAM file (should fail early)
        let temp_bam = NamedTempFile::new().unwrap();
        let bam_path = temp_bam.path();
        
        // No index file exists, should return error
        let result = BamAnalyzer::new(bam_path);
        assert!(result.is_err());
        
        if let Err(VlodError::FileNotFound(msg)) = result {
            assert!(msg.contains("BAM index file not found"));
            assert!(msg.contains(".bam.bai"));
            assert!(msg.contains(".bai"));
        } else {
            panic!("Expected FileNotFound error");
        }
    }

    #[test]
    fn test_bam_analyzer_with_bai_extension() {
        // Create a temporary BAM file
        let temp_bam = NamedTempFile::new().unwrap();
        let bam_path = temp_bam.path();
        
        // Create a BAI index file with .bam.bai extension
        let bai_path = bam_path.with_extension("bam.bai");
        let _temp_bai = File::create(&bai_path).unwrap();
        
        // The BamAnalyzer should find the index but still fail because it's not a real BAM file
        let result = BamAnalyzer::new(bam_path);
        assert!(result.is_err());
        
        // Clean up
        std::fs::remove_file(bai_path).ok();
    }

    #[test]
    fn test_bam_analyzer_with_bai_only_extension() {
        // Create a temporary BAM file
        let temp_bam = NamedTempFile::new().unwrap();
        let bam_path = temp_bam.path();
        
        // Create a BAI index file with .bai extension
        let bai_path = bam_path.with_extension("bai");
        let _temp_bai = File::create(&bai_path).unwrap();
        
        // The BamAnalyzer should find the index but still fail because it's not a real BAM file
        let result = BamAnalyzer::new(bam_path);
        assert!(result.is_err());
        
        // Clean up
        std::fs::remove_file(bai_path).ok();
    }
}