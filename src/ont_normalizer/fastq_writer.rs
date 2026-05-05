use anyhow::Result;
use flate2::Compression;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const BUFFER_SIZE: usize = 8 * 1024 * 1024;

pub enum FastqWriter {
    Gzip(GzEncoder<BufWriter<File>>),
    Plain(BufWriter<File>),
}

impl FastqWriter {
    pub fn new<P: AsRef<Path>>(path: P, gzip: bool, gzip_level: u32) -> Result<Self> {
        let file = File::create(path)?;
        let buffered = BufWriter::with_capacity(BUFFER_SIZE, file);

        if gzip {
            let level = Compression::new(gzip_level.min(9));
            Ok(Self::Gzip(GzEncoder::new(buffered, level)))
        } else {
            Ok(Self::Plain(buffered))
        }
    }

    pub fn write_record(&mut self, id: &str, seq: &[u8], qual: &[u8]) -> Result<()> {
        self.write_all(b"@")?;
        self.write_all(id.as_bytes())?;
        self.write_all(b"\n")?;
        self.write_all(seq)?;
        self.write_all(b"\n+\n")?;

        for q in qual {
            self.write_all(&[q.saturating_add(33)])?;
        }

        self.write_all(b"\n")?;
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        match self {
            Self::Gzip(writer) => {
                writer.finish()?;
            }
            Self::Plain(mut writer) => {
                writer.flush()?;
            }
        }

        Ok(())
    }

    fn write_all(&mut self, bytes: &[u8]) -> Result<()> {
        match self {
            Self::Gzip(writer) => writer.write_all(bytes)?,
            Self::Plain(writer) => writer.write_all(bytes)?,
        }

        Ok(())
    }
}
