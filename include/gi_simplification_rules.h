#ifndef GI_SIMPLIFICATION_RULES_H
#define GI_SIMPLIFICATION_RULES_H

#include <cstdint>
#include <memory>
#include <vector>

#include "gi_basic_types.h"
#include "gi_extremum_merge_forest.h"

namespace GInt {

	// Custom simplification rules: predicates on per-region statistics that
	// can VETO a cancellation / merge the persistence order would otherwise
	// perform. Used by MorseSmaleComplexWithRules (MSC hierarchy) and
	// ForestStatisticsGate (ExtremumMergeForest).
	//
	// Statistics live on the manifolds of the restricted extrema: ascending
	// manifolds of minima and/or descending manifolds of maxima. A merge of
	// extrema A and B combines their regions, so a rule is asked about the
	// two regions about to be joined.
	//
	// REQUIREMENT - rules must be growth-monotone: if a rule refuses joining A
	// and B, it must also refuse joining any region containing A with any
	// region containing B. Both engines drop a refused arc for good (it is
	// never re-examined), and both rely on the resulting hierarchy being a
	// persistence-ordered prefix structure (MSC select-by-persistence,
	// forest BuildRemap). MaxAreaRule satisfies this; a "min area" or
	// "mean value in range" rule would not.

	// Which extrema carry the restriction.
	enum RestrictedManifolds : unsigned {
		RESTRICT_NONE = 0,
		RESTRICT_MINIMA_ASC = 1,  // ascending manifolds of minima
		RESTRICT_MAXIMA_DSC = 2   // descending manifolds of maxima
	};

	// Statistic concept (any type with these members works):
	//   default-constructible as the empty region
	//   void AddCell(INDEX_TYPE cell)          -- accumulate one base cell
	//   static Stat Merge(const Stat&, const Stat&)
	// The cells fed to AddCell are the lattice points of the base region:
	// vertices for ascending manifolds of minima, top-dimensional cells for
	// descending manifolds of maxima.
	struct CellCountStat {
		long long count = 0;
		void AddCell(INDEX_TYPE) { ++count; }
		void AddCells(long long n) { count += n; }
		static CellCountStat Merge(const CellCountStat& a, const CellCountStat& b) {
			CellCountStat r;
			r.count = a.count + b.count;
			return r;
		}
	};

	template<class Stat>
	struct MergeRule {
		virtual ~MergeRule() {}
		// The two regions about to be joined. The MSC passes (dying,
		// survivor); the merge forest passes its two representatives in arc
		// order, so rules should be symmetric in their arguments.
		virtual bool Allows(const Stat& dying, const Stat& survivor) const = 0;
	};

	// Conjunction of rules; an empty set allows everything.
	template<class Stat>
	class MergeRuleSet {
	public:
		void Add(std::shared_ptr<MergeRule<Stat> > rule) { if (rule) mRules.push_back(rule); }
		bool Empty() const { return mRules.empty(); }
		bool Allows(const Stat& dying, const Stat& survivor) const {
			for (size_t i = 0; i < mRules.size(); ++i)
				if (!mRules[i]->Allows(dying, survivor)) return false;
			return true;
		}
	private:
		std::vector<std::shared_ptr<MergeRule<Stat> > > mRules;
	};

	// Refuse any merge whose combined region would exceed maxCells cells.
	struct MaxAreaRule : public MergeRule<CellCountStat> {
		explicit MaxAreaRule(long long maxCells) : maxCells(maxCells) {}
		bool Allows(const CellCountStat& a, const CellCountStat& b) const override {
			return a.count + b.count <= maxCells;
		}
		long long maxCells;
	};

	// ExtremumMergeForest gate holding one statistic per forest node (region
	// id == node id). Prime it before Build, then SetMergeGate(&gate).
	template<class Stat>
	class ForestStatisticsGate : public ExtremumMergeForest::MergeGate {
	public:
		explicit ForestStatisticsGate(const MergeRuleSet<Stat>& rules) : mRules(rules) {}

		// Existing base labeling: latToNode[i] = region/node id of lattice
		// point i (or < 0), cellOf(i) = the cell id handed to AddCell.
		template<class CellOf>
		void PrimeFromLabeling(int32_t numNodes, const int32_t* latToNode, size_t numPoints, CellOf cellOf) {
			mStats.assign((size_t)numNodes, Stat());
			for (size_t i = 0; i < numPoints; ++i) {
				const int32_t r = latToNode[i];
				if (r >= 0) mStats[(size_t)r].AddCell(cellOf(i));
			}
		}

		// Direct access, e.g. to seed CellCountStat from precomputed weights.
		std::vector<Stat>& Stats() { return mStats; }
		const std::vector<Stat>& Stats() const { return mStats; }

		bool AllowMerge(int32_t r1, int32_t r2) override {
			return mRules.Allows(mStats[(size_t)r1], mStats[(size_t)r2]);
		}
		void OnMerge(int32_t survivor, int32_t dying) override {
			mStats[(size_t)survivor] = Stat::Merge(mStats[(size_t)survivor], mStats[(size_t)dying]);
		}

	private:
		MergeRuleSet<Stat> mRules;
		std::vector<Stat> mStats;
	};

} // namespace GInt

#endif
