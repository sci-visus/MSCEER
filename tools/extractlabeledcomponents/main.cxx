// extract all labeled components in the 'filename' file as subregion volumes, and do likewise for the optional additional fields;
//	values without label are set to zero in the mask file, but the additional fields are kept intact;
//	all volumes must have the same dimensions
//	assumes label 0 is the background, and other labels are integers

#include <inttypes.h>

#include "extractlabeledcomponents.h"

int main(int argc, char** argv) {
	float background_value = 0.0f;

	// READ IN THE COMMAND LINE ARGUMENTS
	int64_t dims[3] = {0};
	if (argc < 5) { printf("Usage: X Y Z filename [additional_field]\n"); return 0; }
	sscanf(argv[1], "%" SCNi64, &dims[0]);
	sscanf(argv[2], "%" SCNi64, &dims[1]);
	sscanf(argv[3], "%" SCNi64, &dims[2]);
	char const* mask_filename = argv[4];

	float* mask_field = new float[dims[0]*dims[1]*dims[2]];
	FILE* fin = fopen(mask_filename, "rb");
	assert(fin != NULL);
	fread(mask_field, sizeof mask_field[0], dims[0]*dims[1]*dims[2], fin);
	fclose(fin);

	float *additional_field = NULL;
	if (argc == 6) {
		additional_field = new float[dims[0]*dims[1]*dims[2]];
		FILE *fp = fopen(argv[5], "rb");
		assert(fp != NULL);
		fread(additional_field, sizeof additional_field[0], dims[0]*dims[1]*dims[2], fp);
		fclose(fp);
	}

	auto labeled_components = extract_labeled_components(mask_field, additional_field, dims, background_value);
	for (auto component : labeled_components) {
		// save mask field
		{
			char filename[1024];
			auto ret = snprintf(filename, sizeof filename, "label%" PRIi64 "_mask_position%" PRIi64 "x%" PRIi64 "x%" PRIi64 "_%" PRIi64 "x%" PRIi64 "x%" PRIi64 "_float32.raw",
				component.label, component.offset[0], component.offset[1], component.offset[2], component.dims[0], component.dims[1], component.dims[2]);
			assert(ret < sizeof filename);
			FILE *fp = fopen(filename, "wb");
			assert(fp != NULL);
			fwrite(component.mask_field, sizeof component.mask_field[0], component.dims[0]*component.dims[1]*component.dims[2], fp);
			fclose(fp);
		}

		// save additional field
		{
			char filename[1024];
			auto ret = snprintf(filename, sizeof filename, "label%" PRIi64 "_field0_position%" PRIi64 "x%" PRIi64 "x%" PRIi64 "_%" PRIi64 "x%" PRIi64 "x%" PRIi64 "_float32.raw",
				component.label, component.offset[0], component.offset[1], component.offset[2], component.dims[0], component.dims[1], component.dims[2]);
			assert(ret < sizeof filename);
			FILE *fp = fopen(filename, "wb");
			assert(fp != NULL);
			fwrite(component.additional_field, sizeof component.additional_field[0], component.dims[0]*component.dims[1]*component.dims[2], fp);
			fclose(fp);
		}

		delete[] component.mask_field;
		delete[] component.additional_field;
	}

	return 0;
}
