/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

// a library function to load zgy files, preprocess it, and extract geobodies
#pragma once

void extract_and_preprocess_geobodies(const char *zgy_mask_filepath, const char *zgy_field_filepath)
{
	// 1. read zgy files
	float *mask =
	float *field = 

	// 2. check for noninteger labels
	auto last = std::unique();
	for (auto id : mask) {
		// print?
	}
	
	// 3. extract regions
	extract_labeled_components();

	// 4. construct new field by combining mask and seismic/probability
	if (high_values) {
		for (int64_t i = 0; i < voxel_count; i++) {
			float mask_normalized = mask[];
			float field_normalized = []
			combined[i] = mask_normalized*(amplitude_normalized + 0.001f);
		}
	}

	// 5. dump combined volume and seismic/probability
}