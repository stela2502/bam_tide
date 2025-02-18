

use bigtools;
use crate::bed_data::BedData;



pub struct DataIter<'a> {
    /// all data
    pub data: &'a BedData,
    /// the pointer to the BedData::coverage_data id
    pub current_bin: usize,
    /// the current chromosome name and offset for that chromosome
    pub current_chr: Option<( String, usize, usize )>, 
}

impl<'a> DataIter<'a> {
  pub fn new ( data: &'a BedData) ->Self {
    let ret = Some( (
        data.genome_info[0].0.clone(), 
        data.genome_info[0].1,
        data.genome_info[0].2 
      )
    );
    Self{
      data,
      current_bin: 0,
      current_chr: ret,
    }
  }
}



impl<'a> Iterator for DataIter<'a> {
  type Item = (String, bigtools::Value);

  fn next(&mut self) -> Option<Self::Item> {
    match &self.current_chr {
      Some((chr, size, offset)) => {

        let mut rel_bin = self.current_bin - offset;
        let start:u32 = (rel_bin * self.data.bin_width).try_into().unwrap();
        let mut skipped = false;
        let val = self.data.coverage_data[self.current_bin];
        while self.data.coverage_data[self.current_bin] == val{
          self.current_bin +=1;
          rel_bin +=1;
          skipped = true;
          if (rel_bin * self.data.bin_width) >= *size {
            // shit - we overreached out chr!
            // chould not be an issue as we flip back one in the next step
            break;
          }
        }
        if skipped {
          // we need to report the last zero!
          self.current_bin -=1;
        }
       
        rel_bin = self.current_bin - offset;
        // Create the return value based on current values
        let ret = (
          chr.to_string(),
          bigtools::Value {
            start: start,
            end: (rel_bin * self.data.bin_width + self.data.bin_width)
              .min(*size)
              .try_into()
              .unwrap(),
            value: val as f32,
          },
        );

        // Increment the bin index
        self.current_bin += 1;

        // Check if the current bin exceeds the size of the chromosome
        let rel_bin = self.current_bin - offset;
        if (rel_bin * self.data.bin_width) >= *size {
          // Move to the next chromosome
          self.current_chr = self.data.current_chr_for_id(self.current_bin);
        }
        Some(ret)
        
      }
      None => None, // No more chromosomes to process
    }
  }
}