/*
*
* Copyright (C) 2018 Attila Gyulassy <jediati@sci.utah.edu>
* All rights reserved.
*
* This software may be modified and distributed under the terms
* of the BSD license.  See the LICENSE file for details.
*/

#ifndef INDEX_COMPARERS_H
#define INDEX_COMPARERS_H


#include "gi_basic_types.h"

namespace GInt {

	template<class GridFuncType>
	class IndexCompareLessThan {
	protected:
		GridFuncType* m_func;
	public:
		IndexCompareLessThan(GridFuncType* func) : m_func(func) {}
		inline bool Compare(INDEX_TYPE a, INDEX_TYPE b) const {
			return m_func->IsGreater(b, a);
		}
		inline bool CompareF(FLOATTYPE a, FLOATTYPE b) const {
			return a < b;
		}
		bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
			return m_func->IsGreater(a, b);
		}

		void PrintRule() {
			printf("Comparer::IndexCompareLessThan: Using regulartrilineargrid's IsGreater function\n");
			printf("\t\t returns true if a is smaller than b (func then index)\n");
		}

	};
	template<class GridFuncType>
	class IndexCompareGreaterThan {
	protected:
		GridFuncType* m_func;
	public:
		IndexCompareGreaterThan(GridFuncType* func) : m_func(func) {}
		inline bool Compare(INDEX_TYPE a, INDEX_TYPE b) const {
			return m_func->IsGreater(a, b);
		}
		inline bool CompareF(FLOATTYPE a, FLOATTYPE b) const {
			return a > b;
		}
		bool operator()(INDEX_TYPE a, INDEX_TYPE b) const {
			return m_func->IsGreater(b, a);
		}

		void PrintRule() {
			printf("Comparer::IndexCompareGreaterThan: Using regulartrilineargrid's IsGreater function\n");
			printf("\t\t returns true if a is greater than b (func then index)\n");
		}
	};

}
 
#endif
