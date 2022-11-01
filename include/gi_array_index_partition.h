/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef ARRAY_INDEX_PARTITION_H
#define ARRAY_INDEX_PARTITION_H


#include <math.h>
#include <stdio.h>

namespace GInt {

	class ArrayIndexPartitioner {
	protected:

	public:

		// fills the partition array with numbers such that the 0th element is 0, and last element
		// is the number of elements total. A task with index i should do work between
		// the interval from index >= partition[i] to index < partition[i+1] 
		// equalizes the amount of work between the numthreads indices 
		static void EvenChunkSplit(INDEX_TYPE num_elements, int numThreads, std::vector<INDEX_TYPE>& partition)  {
			partition.clear();
			partition.reserve(numThreads + 1);
			INDEX_TYPE chunksize = num_elements / numThreads;
			INDEX_TYPE remainder = num_elements % numThreads;
			INDEX_TYPE start = 0;
			partition.push_back(start);
			for (int i = 0; i < numThreads; i++) {
				start += chunksize;
				if (i < remainder) start += 1;
				partition.push_back(start);
			}
		}
	};


}

#endif