use scdata::FeatureIndex;
use std::collections::HashMap;

use crate::tags::fast_tag_mapper::FastTagMapper;

pub struct FastTagFeatureIndex<'a> {
    mapper: &'a FastTagMapper,
    name_to_id: HashMap<String, u64>,
}

impl<'a> FastTagFeatureIndex<'a> {
    pub fn new(mapper: &'a FastTagMapper) -> Self {
        let mut name_to_id = HashMap::new();

        for tag in mapper.tags() {
            name_to_id.insert(tag.name.clone(), tag.id as u64);
        }

        Self {
            mapper,
            name_to_id,
        }
    }
}

impl FeatureIndex for FastTagFeatureIndex<'_> {
	fn feature_name(&self, feature_id: u64) -> &str {
	    self.mapper
	        .tag(feature_id as usize)
	        .map(|t| t.name.as_str())
	        .unwrap_or("NA")
	}

    fn feature_id(&self, name: &str) -> Option<u64> {
        self.name_to_id.get(name).copied()
    }

    fn ordered_feature_ids(&self) -> Vec<u64> {
        (0..self.mapper.tag_count() as u64).collect()
    }

    fn to_10x_feature_line(&self, feature_id: u64) -> String {
        let name = &self.mapper.tag(feature_id as usize).expect("Feature id {feature_id} was not found in the feature index!").name;
        format!("{name}\t{name}\tAntibody Capture")
    }
}