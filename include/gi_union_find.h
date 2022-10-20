#ifndef UNION_FIND_H
#define UNION_FIND_H


#include "gi_basic_types.h"
#include "gi_labeling.h"
#include <unordered_map>

namespace GInt{
	//// for now just need a FIND algorithm;
	//INDEX_TYPE Find(INDEX_TYPE id, DenseLabeling<INDEX_TYPE>* id_labeling) {
	//	INDEX_TYPE tlabel = id_labeling->GetLabel(id);
	//	if (id_labeling->GetLabel(id) == id) {
	//		return id;
	//	}
	//	// hopefully we don't have 1-cycles!!!!
	//	//else if (a[a[s]] == s) {
	//	//	return s;
	//	//}
	//	tlabel = Find(tlabel, id_labeling);
	//	id_labeling->SetLabel(id, tlabel);
	//	return tlabel;
	//}

	//should probably use encapsulation, but i'm really lazy here
	class SparseUnionFind : protected std::unordered_map<INDEX_TYPE, INDEX_TYPE> {
	protected:
		// to enable range-based iteration over keys
		typedef std::unordered_map<INDEX_TYPE, INDEX_TYPE>::iterator SUIter;

	public:
		

		int Size() { return this->size(); }


		void ClearAll() {
			clear();
		}
		bool Exists(INDEX_TYPE id) {
			return this->count(id) != 0;
		}

		void MakeSet(INDEX_TYPE id) {
			this->insert_or_assign(id, id);
		}

		
		INDEX_TYPE Find(INDEX_TYPE id) {
			INDEX_TYPE mval = at(id);
			if (mval == id) return id;
			INDEX_TYPE tval = Find(mval);
			at(id) = tval;
			return tval;
		}

		void Union(INDEX_TYPE id, INDEX_TYPE other) {
			INDEX_TYPE tvalid = Find(id);
			INDEX_TYPE tvalother = Find(other);
			if (tvalid < tvalother) at(tvalother) = at(tvalid);
			else at(tvalid) = at(tvalother);
		}


		// to enable range-based iteration over keys
		class SparseUnionFindKeyIterator : public SUIter {

		public:
			SparseUnionFindKeyIterator() : SUIter() {};
			SparseUnionFindKeyIterator(SUIter s) : SUIter(s) {};
			INDEX_TYPE operator->() { return SUIter::operator->()->first; }
			INDEX_TYPE operator*() { return SUIter::operator*().first; }
		};
		SparseUnionFindKeyIterator begin() {
			return std::unordered_map<INDEX_TYPE, INDEX_TYPE>::begin();
		}
		SparseUnionFindKeyIterator end() {
			return std::unordered_map<INDEX_TYPE, INDEX_TYPE>::end();
		}
	};


	//// have block # and id #
	//class SparseBlockedUnionFind : protected SparseBlockedLabeling<INDEX_TYPE> {
	//protected:
	//	// to enable range-based iteration over keys
	//	typedef std::unordered_map<INDEX_TYPE, INDEX_TYPE>::iterator SUIter;

	//public:


	//	int Size() { return this->size(); }


	//	void ClearAll() {
	//		clear();
	//	}
	//	bool Exists(INDEX_TYPE id) {
	//		return this->count(id) != 0;
	//	}

	//	void MakeSet(INDEX_TYPE id) {
	//		this->insert_or_assign(id, id);
	//	}


	//	INDEX_TYPE Find(INDEX_TYPE id) {
	//		INDEX_TYPE mval = at(id);
	//		if (mval == id) return id;
	//		INDEX_TYPE tval = Find(mval);
	//		at(id) = tval;
	//		return tval;
	//	}

	//	void Union(INDEX_TYPE id, INDEX_TYPE other) {
	//		INDEX_TYPE tvalid = Find(id);
	//		INDEX_TYPE tvalother = Find(other);
	//		if (tvalid < tvalother) at(tvalother) = at(tvalid);
	//		else at(tvalid) = at(tvalother);
	//	}


	//	// to enable range-based iteration over keys
	//	class SparseUnionFindKeyIterator : public SUIter {

	//	public:
	//		SparseUnionFindKeyIterator() : SUIter() {};
	//		SparseUnionFindKeyIterator(SUIter s) : SUIter(s) {};
	//		INDEX_TYPE operator->() { return SUIter::operator->()->first; }
	//		INDEX_TYPE operator*() { return SUIter::operator*().first; }
	//	};
	//	SparseUnionFindKeyIterator begin() {
	//		return std::unordered_map<INDEX_TYPE, INDEX_TYPE>::begin();
	//	}
	//	SparseUnionFindKeyIterator end() {
	//		return std::unordered_map<INDEX_TYPE, INDEX_TYPE>::end();
	//	}
	//};


};

#endif