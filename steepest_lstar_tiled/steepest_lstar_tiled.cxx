/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#include <stdint.h>
#include <stdio.h>

#include "localized_morse_complex.h"


int main(int argc, char** argv) {
	if (argc < 8) {
		printf("Usage: X Y Z filename per_x per_y per_z [numthreads=max_on_machine]\n");
		return EXIT_FAILURE;
	}

	int64_t dims[3] = {
		strtoll(argv[1], NULL, 10),
		strtoll(argv[2], NULL, 10),
		strtoll(argv[3], NULL, 10),
	};
	const char* filepath = argv[4];
	int64_t periodic[3] = {
		strtoll(argv[5], NULL, 10),
		strtoll(argv[6], NULL, 10),
		strtoll(argv[7], NULL, 10),
	};
	if (argc >= 9) {
		int64_t thread_count = strtoll(argv[8], NULL, 10);
		omp_set_num_threads((int)thread_count);
	}

	if (periodic[0] != 0 || periodic[1] != 0 || periodic[2] != 0) {
		printf("Periodic grids are not supported\n");
		return EXIT_FAILURE;
	}

	// read dataset
	data_t* data = NULL;
	{
		data = new data_t[dims[0]*dims[1]*dims[2]];
		assert(data != NULL);

		FILE* fp = fopen(filepath, "rb");
		assert(fp != NULL);
		int64_t read = fread(data, dims[0]*dims[1]*dims[2]*sizeof *data, 1, fp);
		assert(read == 1);
		fclose(fp);
	}

	char gradient_filepath[4*4096] = {0};
	snprintf(gradient_filepath, sizeof gradient_filepath, "%s.grad", filepath);

	struct localized_morse_complex1 localized_morse_complex;
	compute_localized_morse_complex_from_grid1(data, dims, &localized_morse_complex);
	save_discrete_gradient_to_file1(&localized_morse_complex, gradient_filepath);
	localized_morse_complex_free1(&localized_morse_complex);

	delete[] data;

	return EXIT_SUCCESS;
}