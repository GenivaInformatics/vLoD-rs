//! VCF file processing functionality

use crate::{Variant, VlodError, VlodResult};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Column indices for VCF parsing
#[derive(Debug, Clone)]
pub struct VcfColumnIndices {
    pub chrom: usize,
    pub pos: usize,
    pub id: usize,
    pub ref_allele: usize,
    pub alt: usize,
    pub qual: usize,
    pub filter: usize,
    pub info: usize,
    pub format: Option<usize>,
    pub samples_start: usize,
}

impl VcfColumnIndices {
    pub fn from_header(header_line: &str) -> VlodResult<Self> {
        let fields: Vec<&str> = header_line.trim().split('\t').collect();

        let chrom = fields
            .iter()
            .position(|&col| col == "CHROM" || col == "#CHROM")
            .ok_or_else(|| {
                VlodError::InvalidVariant("CHROM column not found in VCF header".to_string())
            })?;
        let pos = fields.iter().position(|&col| col == "POS").ok_or_else(|| {
            VlodError::InvalidVariant("POS column not found in VCF header".to_string())
        })?;
        let id = fields.iter().position(|&col| col == "ID").ok_or_else(|| {
            VlodError::InvalidVariant("ID column not found in VCF header".to_string())
        })?;
        let ref_allele = fields.iter().position(|&col| col == "REF").ok_or_else(|| {
            VlodError::InvalidVariant("REF column not found in VCF header".to_string())
        })?;
        let alt = fields.iter().position(|&col| col == "ALT").ok_or_else(|| {
            VlodError::InvalidVariant("ALT column not found in VCF header".to_string())
        })?;
        let qual = fields
            .iter()
            .position(|&col| col == "QUAL")
            .ok_or_else(|| {
                VlodError::InvalidVariant("QUAL column not found in VCF header".to_string())
            })?;
        let filter = fields
            .iter()
            .position(|&col| col == "FILTER")
            .ok_or_else(|| {
                VlodError::InvalidVariant("FILTER column not found in VCF header".to_string())
            })?;
        let info = fields
            .iter()
            .position(|&col| col == "INFO")
            .ok_or_else(|| {
                VlodError::InvalidVariant("INFO column not found in VCF header".to_string())
            })?;
        let format = fields.iter().position(|&col| col == "FORMAT");
        let samples_start = format.map(|f| f + 1).unwrap_or(fields.len());

        Ok(VcfColumnIndices {
            chrom,
            pos,
            id,
            ref_allele,
            alt,
            qual,
            filter,
            info,
            format,
            samples_start,
        })
    }
}

/// Represents a VCF record with essential information
#[derive(Debug, Clone)]
pub struct VcfRecord {
    pub variant: Variant,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

impl VcfRecord {
    pub fn from_line_with_indices(line: &str, indices: &VcfColumnIndices) -> VlodResult<Self> {
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() <= indices.info {
            return Err(VlodError::InvalidVariant(format!(
                "Invalid VCF line format - not enough columns: {}",
                line
            )));
        }

        let chrom = fields[indices.chrom].to_string();
        let pos = fields[indices.pos].parse::<u32>().map_err(|_| {
            VlodError::InvalidVariant(format!("Invalid position: {}", fields[indices.pos]))
        })?;
        let ref_allele = fields[indices.ref_allele].to_string();
        let alt_allele = fields[indices.alt].to_string();

        let variant = Variant::new(chrom, pos, ref_allele, alt_allele);
        let info = fields[indices.info].to_string();
        let format = indices.format.and_then(|f| {
            if f < fields.len() {
                Some(fields[f].to_string())
            } else {
                None
            }
        });
        let samples = if indices.samples_start < fields.len() {
            fields[indices.samples_start..]
                .iter()
                .map(|s| s.to_string())
                .collect()
        } else {
            Vec::new()
        };

        Ok(VcfRecord {
            variant,
            info,
            format,
            samples,
        })
    }

    // Legacy method for backward compatibility
    pub fn from_line(line: &str) -> VlodResult<Self> {
        // Assume standard VCF column order for backward compatibility
        let fields: Vec<&str> = line.split('\t').collect();

        if fields.len() < 8 {
            return Err(VlodError::InvalidVariant(format!(
                "Invalid VCF line format: {}",
                line
            )));
        }

        let chrom = fields[0].to_string();
        let pos = fields[1]
            .parse::<u32>()
            .map_err(|_| VlodError::InvalidVariant(format!("Invalid position: {}", fields[1])))?;
        let ref_allele = fields[3].to_string();
        let alt_allele = fields[4].to_string();

        let variant = Variant::new(chrom, pos, ref_allele, alt_allele);
        let info = fields[7].to_string();
        let format = if fields.len() > 8 {
            Some(fields[8].to_string())
        } else {
            None
        };
        let samples = if fields.len() > 9 {
            fields[9..].iter().map(|s| s.to_string()).collect()
        } else {
            Vec::new()
        };

        Ok(VcfRecord {
            variant,
            info,
            format,
            samples,
        })
    }

    pub fn to_line(&self) -> String {
        let mut line = format!(
            "{}\t{}\t.\t{}\t{}\t.\tPASS\t{}",
            self.variant.chrom,
            self.variant.pos,
            self.variant.ref_allele,
            self.variant.alt_allele,
            self.info
        );

        if let Some(format) = &self.format {
            line.push('\t');
            line.push_str(format);

            for sample in &self.samples {
                line.push('\t');
                line.push_str(sample);
            }
        }

        line
    }
}

/// VCF file reader that handles both compressed and uncompressed files
pub struct VcfReader {
    reader: Box<dyn BufRead>,
}

impl VcfReader {
    pub fn new<P: AsRef<Path>>(path: P) -> VlodResult<Self> {
        let file = File::open(&path)
            .map_err(|_| VlodError::FileNotFound(path.as_ref().to_string_lossy().to_string()))?;

        let reader: Box<dyn BufRead> = if is_gzipped(&path)? {
            let gz_decoder = MultiGzDecoder::new(file);
            Box::new(BufReader::new(gz_decoder))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(VcfReader { reader })
    }

    pub fn records(&mut self) -> VcfRecordIterator<'_> {
        VcfRecordIterator {
            reader: &mut self.reader,
        }
    }

    pub fn header_lines(&mut self) -> VlodResult<Vec<String>> {
        let mut header_lines = Vec::new();
        let mut line = String::new();

        loop {
            line.clear();
            match self.reader.read_line(&mut line)? {
                0 => break, // EOF
                _ => {
                    if line.starts_with('#') {
                        header_lines.push(line.trim_end().to_string());
                    } else {
                        // We've reached the first data line, stop reading header
                        break;
                    }
                }
            }
        }

        Ok(header_lines)
    }
}

/// Iterator over VCF records
pub struct VcfRecordIterator<'a> {
    reader: &'a mut Box<dyn BufRead>,
}

impl<'a> Iterator for VcfRecordIterator<'a> {
    type Item = VlodResult<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut line = String::new();

        loop {
            line.clear();
            match self.reader.read_line(&mut line) {
                Ok(0) => return None, // EOF
                Ok(_) => {
                    let line = line.trim_end();
                    if line.starts_with('#') {
                        continue; // Skip header lines
                    }
                    if line.is_empty() {
                        continue; // Skip empty lines
                    }

                    return Some(VcfRecord::from_line(line));
                }
                Err(e) => return Some(Err(VlodError::Io(e))),
            }
        }
    }
}

/// Check if a file is gzipped
pub fn is_gzipped<P: AsRef<Path>>(path: P) -> VlodResult<bool> {
    let mut file = File::open(path)?;
    let mut buffer = [0; 2];

    match file.read_exact(&mut buffer) {
        Ok(()) => Ok(buffer == [0x1f, 0x8b]),
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(VlodError::Io(e)),
    }
}

/// Read VCF variants from a file and return them as a vector
pub fn read_vcf_variants<P: AsRef<Path>>(path: P) -> VlodResult<Vec<Variant>> {
    let file = File::open(&path)
        .map_err(|_| VlodError::FileNotFound(path.as_ref().to_string_lossy().to_string()))?;

    let reader: Box<dyn BufRead> = if is_gzipped(&path)? {
        let gz_decoder = MultiGzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut variants = Vec::new();
    let mut column_indices: Option<VcfColumnIndices> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        if line.starts_with("##") {
            continue; // Skip metadata lines
        }

        if line.starts_with("#CHROM") || line.starts_with("#") {
            // Parse header to get column indices
            column_indices = Some(VcfColumnIndices::from_header(&line)?);
            continue;
        }

        if line.is_empty() {
            continue;
        }

        // Parse variant line
        let record = if let Some(ref indices) = column_indices {
            // Use header-based parsing if we found a header
            VcfRecord::from_line_with_indices(&line, indices)
        } else {
            // Fall back to standard VCF column order if no header found
            VcfRecord::from_line(&line)
        };

        match record {
            Ok(record) => {
                // Handle multiple alternative alleles
                let alt_alleles: Vec<&str> = record.variant.alt_allele.split(',').collect();
                for alt_allele in alt_alleles {
                    let variant = Variant::new(
                        record.variant.chrom.clone(),
                        record.variant.pos,
                        record.variant.ref_allele.clone(),
                        alt_allele.to_string(),
                    );
                    variants.push(variant);
                }
            }
            Err(e) => {
                log::warn!("Skipping invalid VCF record: {}", e);
                continue;
            }
        }
    }

    Ok(variants)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_vcf_record_from_line() {
        let line = "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30";
        let record = VcfRecord::from_line(line).unwrap();

        assert_eq!(record.variant.chrom, "chr1");
        assert_eq!(record.variant.pos, 100);
        assert_eq!(record.variant.ref_allele, "A");
        assert_eq!(record.variant.alt_allele, "T");
        assert_eq!(record.info, "DP=30");
    }

    #[test]
    fn test_vcf_record_to_line() {
        let variant = Variant::new("chr1".to_string(), 100, "A".to_string(), "T".to_string());
        let record = VcfRecord {
            variant,
            info: "DP=30".to_string(),
            format: None,
            samples: Vec::new(),
        };

        let line = record.to_line();
        assert_eq!(line, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30");
    }

    #[test]
    fn test_read_vcf_variants() {
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "##fileformat=VCFv4.2").unwrap();
        writeln!(temp_file, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();
        writeln!(temp_file, "chr1\t100\t.\tA\tT\t.\tPASS\tDP=30").unwrap();
        writeln!(temp_file, "chr2\t200\t.\tG\tC,A\t.\tPASS\tDP=40").unwrap();

        let variants = read_vcf_variants(temp_file.path()).unwrap();
        assert_eq!(variants.len(), 3); // One variant with single alt + one with two alts

        assert_eq!(variants[0].chrom, "chr1");
        assert_eq!(variants[0].alt_allele, "T");

        assert_eq!(variants[1].chrom, "chr2");
        assert_eq!(variants[1].alt_allele, "C");

        assert_eq!(variants[2].chrom, "chr2");
        assert_eq!(variants[2].alt_allele, "A");
    }
}
