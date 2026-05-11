use crate::ont_normalizer::primer::model::{
    BarcodeMatcherSpec, PrimerPart, PrimerStructure,
};
use anyhow::{bail, Result};
use std::fmt;

#[derive(Debug, Clone)]
pub struct PrimerStructureCli {
    pub name: String,
    pub adapter: String,
    pub adapter_suffix_min_len: usize,
    pub max_adapter_mismatches: usize,
    pub structure: String,
}

impl PrimerStructureCli {
    pub fn tenx_v3() -> Self {
        PrimerStructure::tenx_v3().to_cli()
    }

    pub fn to_structure(&self) -> Result<PrimerStructure> {
        Ok(PrimerStructure {
            name: self.name.clone(),
            adapter: self.adapter.as_bytes().to_vec(),
            adapter_suffix_min_len: self.adapter_suffix_min_len,
            max_adapter_mismatches: self.max_adapter_mismatches,
            parts: PrimerStructure::parse_parts(&self.structure)?,
        })
    }

    pub fn flattened_cli(&self) -> String {
        format!(
            "--primer-name {} --primer-adapter {} --primer-adapter-suffix-min-len {} --primer-max-adapter-mismatches {} --primer-structure '{}'",
            self.name,
            self.adapter,
            self.adapter_suffix_min_len,
            self.max_adapter_mismatches,
            self.structure,
        )
    }
}

impl PrimerStructure {
    pub fn tenx_v3() -> Self {
        Self {
            name: "10x-v3".to_string(),
            adapter: b"CTACACGACGCTCTTCCGATCT".to_vec(),
            adapter_suffix_min_len: 13,
            max_adapter_mismatches: 4,
            parts: vec![
                PrimerPart::CellId {
                    name: "CB".to_string(),
                    len: 16,
                    matcher: BarcodeMatcherSpec::Any,
                },
                PrimerPart::Umi {
                    name: "UB".to_string(),
                    len: 12,
                },
                PrimerPart::PolyT {
                    min_len: 10,
                    max_non_t: 2,
                },
                PrimerPart::Insert,
            ],
        }
    }

    pub fn to_cli(&self) -> PrimerStructureCli {
        PrimerStructureCli {
            name: self.name.clone(),
            adapter: String::from_utf8_lossy(&self.adapter).to_string(),
            adapter_suffix_min_len: self.adapter_suffix_min_len,
            max_adapter_mismatches: self.max_adapter_mismatches,
            structure: self.parts_as_cli_string(),
        }
    }

    pub fn parts_as_cli_string(&self) -> String {
        self.parts
            .iter()
            .map(PrimerPart::as_cli_token)
            .collect::<Vec<_>>()
            .join(",")
    }

    pub fn flattened_cli(&self) -> String {
        self.to_cli().flattened_cli()
    }

    pub fn parse_parts(s: &str) -> Result<Vec<PrimerPart>> {
        let mut parts = Vec::new();

        for raw in s.split(',') {
            let token = raw.trim();

            if token.is_empty() {
                continue;
            }

            if let Some(rest) = token.strip_prefix("FIXED:") {
                let mut pieces = rest.split(":mm=");
                let seq = pieces.next().unwrap_or_default();
                let max_mismatches = pieces
                    .next()
                    .map(|x| x.parse::<usize>())
                    .transpose()?
                    .unwrap_or(0);

                parts.push(PrimerPart::Fixed {
                    name: "FIXED".to_string(),
                    seq: seq.as_bytes().to_vec(),
                    max_mismatches,
                });
            } else if let Some(rest) = token.strip_prefix("RANDOM:") {
                let (min_len, max_len) = Self::parse_range(rest)?;

                parts.push(PrimerPart::Random {
                    name: "RANDOM".to_string(),
                    min_len,
                    max_len,
                });
            } else if let Some(rest) = token.strip_prefix("CELL:") {
                let len = rest
                    .split(':')
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("bad CELL token: {token}"))?
                    .parse::<usize>()?;

                parts.push(PrimerPart::CellId {
                    name: "CB".to_string(),
                    len,
                    matcher: BarcodeMatcherSpec::Any,
                });
            } else if let Some(rest) = token.strip_prefix("UMI:") {
                parts.push(PrimerPart::Umi {
                    name: "UMI".to_string(),
                    len: rest.parse::<usize>()?,
                });
            } else if let Some(rest) = token.strip_prefix("POLYT:") {
                let mut pieces = rest.split(":non_t=");
                let min_len = pieces
                    .next()
                    .ok_or_else(|| anyhow::anyhow!("bad POLYT token: {token}"))?
                    .parse::<usize>()?;
                let max_non_t = pieces
                    .next()
                    .map(|x| x.parse::<usize>())
                    .transpose()?
                    .unwrap_or(2);

                parts.push(PrimerPart::PolyT { min_len, max_non_t });
            } else if token == "INSERT" {
                parts.push(PrimerPart::Insert);
            } else {
                bail!("unknown primer grammar token: {token}");
            }
        }

        Ok(parts)
    }

    fn parse_range(s: &str) -> Result<(usize, usize)> {
        let Some((a, b)) = s.split_once('-') else {
            let n = s.parse::<usize>()?;
            return Ok((n, n));
        };

        let min_len = a.parse::<usize>()?;
        let max_len = b.parse::<usize>()?;

        if min_len > max_len {
            bail!("invalid range {s}: min > max");
        }

        Ok((min_len, max_len))
    }
}

impl PrimerPart {
    pub fn as_cli_token(&self) -> String {
        match self {
            PrimerPart::Fixed {
                seq,
                max_mismatches,
                ..
            } => format!(
                "FIXED:{}:mm={}",
                String::from_utf8_lossy(seq),
                max_mismatches
            ),

            PrimerPart::Random { min_len, max_len, .. } => {
                format!("RANDOM:{min_len}-{max_len}")
            }

            PrimerPart::CellId { len, matcher, .. } => match matcher {
                BarcodeMatcherSpec::Any => format!("CELL:{len}"),
                BarcodeMatcherSpec::Whitelist {
                    path,
                    max_mismatches,
                } => format!(
                    "CELL:{len}:whitelist={}:mm={}",
                    path.display(),
                    max_mismatches
                ),
            },

            PrimerPart::Umi { len, .. } => format!("UMI:{len}"),

            PrimerPart::PolyT { min_len, max_non_t } => {
                format!("POLYT:{min_len}:non_t={max_non_t}")
            }

            PrimerPart::Insert => "INSERT".to_string(),
        }
    }
}

impl fmt::Display for PrimerStructure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "PrimerStructure: {}", self.name)?;
        writeln!(
            f,
            "  adapter                : {}",
            String::from_utf8_lossy(&self.adapter)
        )?;
        writeln!(
            f,
            "  adapter_suffix_min_len : {}",
            self.adapter_suffix_min_len
        )?;
        writeln!(
            f,
            "  max_adapter_mismatches : {}",
            self.max_adapter_mismatches
        )?;
        writeln!(f, "  parts                  : {}", self.parts_as_cli_string())?;
        writeln!(f, "  flattened_cli          : {}", self.flattened_cli())
    }
}