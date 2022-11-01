/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#include "extractandpreprocessgeobodies.h"

int main(int argc, char **argv)
{
	extract_and_preprocess_geobodies(zgy_mask_filepath, zgy_field_filepath);

	return 0;
}