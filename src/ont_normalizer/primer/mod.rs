pub mod detector;
pub mod grammer;
pub mod model;
//pub mod primer_table; //no time for this!

pub use detector::PrimerDetector;
pub use grammer::PrimerStructureCli;
pub use model::{
    BarcodeMatcherSpec, Orientation, PrimerHit, PrimerPart, PrimerSplit, PrimerStructure,
};
