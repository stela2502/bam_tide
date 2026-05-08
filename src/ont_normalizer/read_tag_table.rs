use anyhow::{Context, Result};
use clap::Args;
use flate2::read::MultiGzDecoder;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, Read},
    path::{Path, PathBuf},
};
use mapping_info::MappingInfo;

use csv::WriterBuilder;
use std::io::Write;
use std::fmt;

#[derive(Debug, Clone, Args)]
pub struct ReadTagTableCli {
    /// Optional external read-tag table, TSV or TSV.GZ.
    ///
    /// Maps BAM query names / read IDs to observed cell barcode and UMI
    /// information. This is useful when read names are preserved after
    /// preprocessing and alignment.
    #[arg(long = "read-tag-table")]
    pub read_tag_table: Option<PathBuf>,

    /// Column containing the read id / BAM query name.
    #[arg(long = "rt-read-id-column", default_value = "output_read_id")]
    pub rt_read_id_column: String,

    /// Column containing the observed cell barcode.
    #[arg(long = "rt-cell-column", default_value = "raw_cb")]
    pub rt_cell_column: String,

    /// Column containing the observed cell barcode quality string.
    #[arg(long = "rt-cell-qual-column", default_value = "quality_cb")]
    pub rt_cell_qual_column: String,

    /// Column containing the observed UMI.
    #[arg(long = "rt-umi-column", default_value = "raw_umi")]
    pub rt_umi_column: String,

    /// Column containing the observed UMI quality string.
    #[arg(long = "rt-umi-qual-column", default_value = "quality_umi")]
    pub rt_umi_qual_column: String,

    /// Optional provenance column containing the original read id.
    #[arg(long = "rt-original-read-id-column", default_value = "original_read_id")]
    pub rt_original_read_id_column: String,
}

impl ReadTagTableCli {
    pub fn to_config(&self) -> Option<ReadTagTableConfig> {
        Some(ReadTagTableConfig {
            path: self.read_tag_table.clone()?,
            read_id_column: self.rt_read_id_column.clone(),
            original_read_id_column: self.rt_original_read_id_column.clone(),
            cell_column: self.rt_cell_column.clone(),
            cell_qual_column: self.rt_cell_qual_column.clone(),
            umi_column: self.rt_umi_column.clone(),
            umi_qual_column: self.rt_umi_qual_column.clone(),
        })
    }
}

#[derive(Debug, Clone)]
pub struct ReadTagTableConfig {
    pub path: PathBuf,
    pub read_id_column: String,
    pub original_read_id_column: String,
    pub cell_column: String,
    pub cell_qual_column: String,
    pub umi_column: String,
    pub umi_qual_column: String,
}

#[derive(Debug, Clone)]
pub struct ReadTagRecord {
    pub read_id: String,
    pub original_read_id: Option<String>,
    pub cell: String,
    pub cell_qual: Option<String>,
    pub umi: String,
    pub umi_qual: Option<String>,
}

#[derive(Debug)]
pub struct ReadTagTable {
    records: HashMap<String, ReadTagRecord>,
    //mapping_info: MappingInfo,
}

impl ReadTagTable {
    pub fn from_config(config: &ReadTagTableConfig) -> Result<Self> {
        let mut mapping_info = mapping_info::MappingInfo::new( None, 0.0, 0);        

        let reader = open_maybe_gz(&config.path)?;

        let mut rdr = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .flexible(true)
            .from_reader(reader);

        let headers = rdr.headers()?.clone();

        let read_id_ix = required_column_ix(&headers, &config.read_id_column)?;
        let cell_ix = required_column_ix(&headers, &config.cell_column)?;
        let umi_ix = required_column_ix(&headers, &config.umi_column)?;

        let original_read_id_ix = optional_column_ix(&headers, &config.original_read_id_column);
        let cell_qual_ix = optional_column_ix(&headers, &config.cell_qual_column);
        let umi_qual_ix = optional_column_ix(&headers, &config.umi_qual_column);

        let mut records = HashMap::new();

        for rec in rdr.records() {
            let rec = rec?;

            let read_id = rec.get(read_id_ix).unwrap_or("").trim();
            let cell = rec.get(cell_ix).unwrap_or("").trim();
            let umi = rec.get(umi_ix).unwrap_or("").trim();

            if read_id.is_empty() || cell.is_empty() || umi.is_empty() {
                continue;
            }

            let record = ReadTagRecord {
                read_id: read_id.to_string(),
                original_read_id: get_optional(&rec, original_read_id_ix),
                cell: cell.to_string(),
                cell_qual: get_optional(&rec, cell_qual_ix),
                umi: umi.to_string(),
                umi_qual: get_optional(&rec, umi_qual_ix),
            };

            records.insert(read_id.to_string(), record);
        }

        mapping_info.stop_file_io_time();

        
        let (h, m, s, ms) =
            MappingInfo::split_duration(mapping_info.file_io_time);
        println!(
            "Read-tag table loaded: {} entries in {}:{:02}:{:02}.{:03}",
            records.len(),
            h,
            m,
            s,
            ms,
        );

        Ok(Self { 
            records,
           //mapping_info,
        })
    }

    pub fn get(&self, read_id: &str) -> Option<&ReadTagRecord> {
        self.records.get(read_id)
    }

    pub fn cell_umi_for_read(&self, read_id: &str) -> Option<(&str, &str)> {
        let rec = self.records.get(read_id)?;
        Some((rec.cell.as_str(), rec.umi.as_str()))
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn pair_counts(&self) -> HashMap<(String, String), u64> {
        let mut counts = HashMap::new();

        for rec in self.records.values() {
            if rec.cell.is_empty() || rec.umi.is_empty() {
                continue;
            }

            *counts
                .entry((rec.cell.clone(), rec.umi.clone()))
                .or_insert(0) += 1;
        }

        counts
    }

    pub fn summarize_pairs(&self, min_pair_count: u64, min_cell_umis: u64) -> PairStats {
        let counts = self.pair_counts();
        PairStats::from_counts(&counts, min_pair_count, min_cell_umis)
    }
}

#[derive(Debug, Clone)]
pub struct PairStats {
    pub cell_entries: usize,
    pub unique_cell_umi_combos: usize,
    pub total_pair_observations: u64,
    pub umis_per_cell: Summary,
    pub detections_per_cell_umi: Summary,
}

impl PairStats {
    pub fn from_counts(
        cb_umi_counts: &HashMap<(String, String), u64>,
        min_pair_count: u64,
        min_cell_umis: u64,
    ) -> Self {
        let mut cell_to_umis: HashMap<&str, HashSet<&str>> = HashMap::new();

        for ((cell, umi), count) in cb_umi_counts {
            if *count < min_pair_count {
                continue;
            }

            cell_to_umis
                .entry(cell.as_str())
                .or_default()
                .insert(umi.as_str());
        }

        let valid_cells: HashSet<&str> = cell_to_umis
            .iter()
            .filter_map(|(cell, umis)| {
                if umis.len() as u64 >= min_cell_umis {
                    Some(*cell)
                } else {
                    None
                }
            })
            .collect();

        let mut final_cell_to_umi_count: HashMap<&str, u64> = HashMap::new();
        let mut final_pair_counts = Vec::new();
        let mut total_pair_observations = 0u64;

        for ((cell, _umi), count) in cb_umi_counts {
            if *count < min_pair_count {
                continue;
            }

            if !valid_cells.contains(cell.as_str()) {
                continue;
            }

            *final_cell_to_umi_count.entry(cell.as_str()).or_insert(0) += 1;
            final_pair_counts.push(*count);
            total_pair_observations += *count;
        }

        let umis_per_cell: Vec<u64> = final_cell_to_umi_count.values().copied().collect();

        Self {
            cell_entries: final_cell_to_umi_count.len(),
            unique_cell_umi_combos: final_pair_counts.len(),
            total_pair_observations,
            umis_per_cell: summarize(umis_per_cell),
            detections_per_cell_umi: summarize(final_pair_counts),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Summary {
    pub n: usize,
    pub mean: f64,
    pub median: f64,
    pub min: u64,
    pub max: u64,
}

fn summarize(mut values: Vec<u64>) -> Summary {
    if values.is_empty() {
        return Summary {
            n: 0,
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
        };
    }

    values.sort_unstable();

    let n = values.len();
    let min = values[0];
    let max = values[n - 1];

    let sum: u128 = values.iter().map(|&x| x as u128).sum();
    let mean = sum as f64 / n as f64;

    let median = if n % 2 == 0 {
        let a = values[n / 2 - 1];
        let b = values[n / 2];
        (a as f64 + b as f64) / 2.0
    } else {
        values[n / 2] as f64
    };

    Summary {
        n,
        mean,
        median,
        min,
        max,
    }
}

fn required_column_ix(headers: &csv::StringRecord, name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|h| h == name)
        .with_context(|| format!("Could not find required column '{name}'"))
}

fn optional_column_ix(headers: &csv::StringRecord, name: &str) -> Option<usize> {
    headers.iter().position(|h| h == name)
}

fn get_optional(rec: &csv::StringRecord, ix: Option<usize>) -> Option<String> {
    let value = rec.get(ix?)?.trim();

    if value.is_empty() {
        None
    } else {
        Some(value.to_string())
    }
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader = BufReader::new(file);

    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(MultiGzDecoder::new(reader)))
    } else {
        Ok(Box::new(reader))
    }
}


pub const READ_TAG_TABLE_COLUMNS: [&str; 15] = [
    "output_read_id",
    "original_read_id",
    "molecule_index",
    "orientation",
    "adapter_start",
    "adapter_end",
    "segment_start",
    "segment_end",
    "raw_cb",
    "quality_cb",
    "raw_umi",
    "quality_umi",
    "poly_t_start",
    "poly_t_len",
    "status",
];

#[derive(Debug, Clone)]
pub struct ReadTagWriteRecord<'a> {
    pub output_read_id: &'a str,
    pub original_read_id: Option<&'a str>,
    pub molecule_index: Option<usize>,
    pub orientation: Option<&'a str>,
    pub adapter_start: Option<usize>,
    pub adapter_end: Option<usize>,
    pub segment_start: Option<usize>,
    pub segment_end: Option<usize>,
    pub raw_cb: Option<String>,
    pub quality_cb: Option<String>,
    pub raw_umi: Option<String>,
    pub quality_umi: Option<String>,
    pub poly_t_start: Option<usize>,
    pub poly_t_len: Option<usize>,
    pub status: &'a str,
}

pub struct ReadTagTableWriter<W: Write> {
    writer: csv::Writer<W>,
}

impl<W: Write> ReadTagTableWriter<W> {
    pub fn new(inner: W) -> Result<Self> {
        let mut writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(inner);

        writer
            .write_record(READ_TAG_TABLE_COLUMNS)
            .context("writing read-tag table header")?;

        Ok(Self { writer })
    }

    pub fn write_record(&mut self, rec: &ReadTagWriteRecord<'_>) -> Result<()> {
        let molecule_index = opt_usize(rec.molecule_index);
        let adapter_start = opt_usize(rec.adapter_start);
        let adapter_end = opt_usize(rec.adapter_end);
        let segment_start = opt_usize(rec.segment_start);
        let segment_end = opt_usize(rec.segment_end);
        let poly_t_start = opt_usize(rec.poly_t_start);
        let poly_t_len = opt_usize(rec.poly_t_len);

        self.writer
            .write_record([
                rec.output_read_id,
                rec.original_read_id.unwrap_or(""),
                molecule_index.as_str(),
                rec.orientation.unwrap_or(""),
                adapter_start.as_str(),
                adapter_end.as_str(),
                segment_start.as_str(),
                segment_end.as_str(),
                rec.raw_cb.as_deref().unwrap_or(""),
                rec.quality_cb.as_deref().unwrap_or(""),
                rec.raw_umi.as_deref().unwrap_or(""),
                rec.quality_umi.as_deref().unwrap_or(""),
                poly_t_start.as_str(),
                poly_t_len.as_str(),
                rec.status,
            ])
            .with_context(|| {
                format!(
                    "writing read-tag table row for output_read_id '{}'",
                    rec.output_read_id
                )
            })?;

        Ok(())
    }

    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("flushing read-tag table writer")?;
        Ok(())
    }
}


impl fmt::Display for ReadTagRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "read_id={}, original_read_id={}, cell={}, cell_qual={}, umi={}, umi_qual={}",
            self.read_id,
            self.original_read_id.as_deref().unwrap_or("-"),
            self.cell,
            self.cell_qual.as_deref().unwrap_or("-"),
            self.umi,
            self.umi_qual.as_deref().unwrap_or("-"),
        )
    }
}

impl fmt::Display for ReadTagTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "ReadTagTable: {} records", self.records.len())?;

        for (i, record) in self.records.values().take(3).enumerate() {
            writeln!(f, "  [{}] {}", i, record)?;
        }

        Ok(())
    }
}


fn opt_usize(value: Option<usize>) -> String {
    value.map(|v| v.to_string()).unwrap_or_default()
}


#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::fs;
    use std::io::Write;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp_path(name: &str) -> PathBuf {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();

        std::env::temp_dir().join(format!("read_tag_table_test_{stamp}_{name}"))
    }

    fn write_text(path: &Path, text: &str) {
        fs::write(path, text).unwrap();
    }

    fn default_config(path: PathBuf) -> ReadTagTableConfig {
        ReadTagTableConfig {
            path,
            read_id_column: "output_read_id".to_string(),
            original_read_id_column: "original_read_id".to_string(),
            cell_column: "raw_cb".to_string(),
            cell_qual_column: "quality_cb".to_string(),
            umi_column: "raw_umi".to_string(),
            umi_qual_column: "quality_umi".to_string(),
        }
    }

    #[test]
    fn reads_plain_tsv_table() {
        let path = tmp_path("table.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\toriginal_read_id\traw_cb\tquality_cb\traw_umi\tquality_umi\n",
                "read1\torig1\tCELL1\tIIII\tUMI1\tJJJJ\n",
                "read2\torig2\tCELL1\tIIII\tUMI2\tJJJJ\n",
                "read3\t\tCELL2\t\tUMI3\t\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 3);
        assert!(!table.is_empty());

        let rec = table.get("read1").unwrap();
        assert_eq!(rec.read_id, "read1");
        assert_eq!(rec.original_read_id.as_deref(), Some("orig1"));
        assert_eq!(rec.cell, "CELL1");
        assert_eq!(rec.cell_qual.as_deref(), Some("IIII"));
        assert_eq!(rec.umi, "UMI1");
        assert_eq!(rec.umi_qual.as_deref(), Some("JJJJ"));

        let rec = table.get("read3").unwrap();
        assert_eq!(rec.original_read_id, None);
        assert_eq!(rec.cell_qual, None);
        assert_eq!(rec.umi_qual, None);

        fs::remove_file(path).ok();
    }

    #[test]
    fn skips_rows_with_missing_required_values() {
        let path = tmp_path("missing.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "\tCELL2\tUMI2\n",
                "read3\t\tUMI3\n",
                "read4\tCELL4\t\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 1);
        assert!(table.get("read1").is_some());
        assert!(table.get("read3").is_none());

        fs::remove_file(path).ok();
    }

    #[test]
    fn missing_required_column_returns_error() {
        let path = tmp_path("bad_columns.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\traw_cb\n",
                "read1\tCELL1\n",
            ),
        );

        let err = ReadTagTable::from_config(&default_config(path.clone()))
            .unwrap_err()
            .to_string();

        assert!(
            err.contains("Could not find required column 'raw_umi'"),
            "unexpected error: {err}"
        );

        fs::remove_file(path).ok();
    }

    #[test]
    fn reads_gzipped_tsv_table() {
        let path = tmp_path("table.tsv.gz");

        let file = File::create(&path).unwrap();
        let mut gz = GzEncoder::new(file, Compression::default());

        gz.write_all(
            concat!(
                "output_read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL2\tUMI2\n",
            )
            .as_bytes(),
        )
        .unwrap();

        gz.finish().unwrap();

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        assert_eq!(table.len(), 2);
        assert_eq!(table.cell_umi_for_read("read1"), Some(("CELL1", "UMI1")));
        assert_eq!(table.cell_umi_for_read("missing"), None);

        fs::remove_file(path).ok();
    }

    #[test]
    fn pair_counts_counts_cell_umi_observations() {
        let path = tmp_path("pairs.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL1\tUMI1\n",
                "read3\tCELL1\tUMI2\n",
                "read4\tCELL2\tUMI3\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();
        let counts = table.pair_counts();

        assert_eq!(counts.get(&("CELL1".to_string(), "UMI1".to_string())), Some(&2));
        assert_eq!(counts.get(&("CELL1".to_string(), "UMI2".to_string())), Some(&1));
        assert_eq!(counts.get(&("CELL2".to_string(), "UMI3".to_string())), Some(&1));

        fs::remove_file(path).ok();
    }

    #[test]
    fn summarize_pairs_applies_pair_and_cell_thresholds() {
        let path = tmp_path("summary.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\traw_cb\traw_umi\n",
                "r1\tCELL1\tUMI1\n",
                "r2\tCELL1\tUMI1\n",
                "r3\tCELL1\tUMI2\n",
                "r4\tCELL1\tUMI2\n",
                "r5\tCELL2\tUMI3\n",
                "r6\tCELL2\tUMI3\n",
                "r7\tCELL3\tUMI4\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();

        let stats = table.summarize_pairs(2, 2);

        assert_eq!(stats.cell_entries, 1);
        assert_eq!(stats.unique_cell_umi_combos, 2);
        assert_eq!(stats.total_pair_observations, 4);

        assert_eq!(stats.umis_per_cell.n, 1);
        assert_eq!(stats.umis_per_cell.min, 2);
        assert_eq!(stats.umis_per_cell.max, 2);

        assert_eq!(stats.detections_per_cell_umi.n, 2);
        assert_eq!(stats.detections_per_cell_umi.min, 2);
        assert_eq!(stats.detections_per_cell_umi.max, 2);

        fs::remove_file(path).ok();
    }

    #[test]
    fn display_for_record_is_informative() {
        let rec = ReadTagRecord {
            read_id: "read1".to_string(),
            original_read_id: Some("orig1".to_string()),
            cell: "CELL1".to_string(),
            cell_qual: None,
            umi: "UMI1".to_string(),
            umi_qual: Some("JJJJ".to_string()),
        };

        let text = rec.to_string();

        assert!(text.contains("read_id=read1"));
        assert!(text.contains("original_read_id=orig1"));
        assert!(text.contains("cell=CELL1"));
        assert!(text.contains("cell_qual=-"));
        assert!(text.contains("umi=UMI1"));
        assert!(text.contains("umi_qual=JJJJ"));
    }

    #[test]
    fn display_for_table_shows_count_and_at_most_three_records() {
        let path = tmp_path("display.tsv");

        write_text(
            &path,
            concat!(
                "output_read_id\traw_cb\traw_umi\n",
                "read1\tCELL1\tUMI1\n",
                "read2\tCELL2\tUMI2\n",
                "read3\tCELL3\tUMI3\n",
                "read4\tCELL4\tUMI4\n",
            ),
        );

        let table = ReadTagTable::from_config(&default_config(path.clone())).unwrap();
        let text = table.to_string();

        assert!(text.contains("ReadTagTable: 4 records"));

        let shown_records = text.lines().filter(|line| line.trim_start().starts_with('[')).count();

        assert_eq!(shown_records, 3);

        fs::remove_file(path).ok();
    }

    #[test]
    fn writer_writes_header_and_rows() {
        let mut out = Vec::new();

        {
            let mut writer = ReadTagTableWriter::new(&mut out).unwrap();

            writer
                .write_record(&ReadTagWriteRecord {
                    output_read_id: "read1",
                    original_read_id: Some("orig1"),
                    molecule_index: Some(7),
                    orientation: Some("Forward"),
                    adapter_start: Some(10),
                    adapter_end: Some(30),
                    segment_start: Some(31),
                    segment_end: Some(100),
                    raw_cb: Some("CELL1".to_string()),
                    quality_cb: Some("IIII".to_string()),
                    raw_umi: Some("UMI1".to_string()),
                    quality_umi: Some("JJJJ".to_string()),
                    poly_t_start: Some(101),
                    poly_t_len: Some(20),
                    status: "ok",
                })
                .unwrap();

            writer.flush().unwrap();
        }

        let text = String::from_utf8(out).unwrap();

        assert!(text.starts_with(&READ_TAG_TABLE_COLUMNS.join("\t")));
        assert!(text.contains("read1\torig1\t7\tForward\t10\t30\t31\t100\tCELL1\tIIII\tUMI1\tJJJJ\t101\t20\tok"));
    }
}