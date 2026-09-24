#ifndef GI_MORSE_SMALE_COMPLEX_RULES_H
#define GI_MORSE_SMALE_COMPLEX_RULES_H

#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

#include "gi_labeling.h"
#include "gi_morse_smale_complex_basic.h"
#include "gi_simplification_rules.h"

namespace GInt {

	// Mixin over any MorseSmaleComplexBasic-derived class whose cancellation
	// hierarchy should obey custom simplification rules
	// (gi_simplification_rules.h) evaluated on per-manifold statistics of the
	// restricted extrema. BASE_MSC must be constructible from
	// (grad, mesh, func); MorseSmaleComplexBasic, PartitionLocalMsc and
	// ReconciledGlobalMsc all are. BASE_MSC::isValid runs first, so a base's
	// own restrictions (e.g. PartitionLocalMsc's frozen/interior checks) keep
	// applying.
	//
	// Usage (serial):
	//   MorseSmaleComplexWithRules<...> msc(grad, mesh, func,
	//       RESTRICT_MINIMA_ASC, rules);
	//   msc.ComputeFromGrad();
	//   msc.PrimeStatistics();          // ends the compute phase
	//   msc.ComputeHierarchy(persLimit);
	//
	// PrimeStatistics() builds the base manifold labeling itself by tracing
	// each restricted extremum's unsimplified manifold through the GLOBAL
	// gradient - exact even for a partition-local MSC whose manifolds cross
	// the partition boundary. baseSeedCells() says which extremum cells make
	// up a node's base region (a reconciled node imports several); alternatively
	// PrimeStatisticsFromLabeling() takes an existing per-cell labeling.
	// ComputeHierarchy primes automatically if neither was called.
	//
	// Statistics are kept per manifold record (the mans array), so merged
	// manifolds produced by cancel() inherit Stat::Merge of their two parts
	// through the onManifoldCreated hook - cancel() itself is untouched.
	template<class BASE_MSC, class STAT_TYPE>
	class SimplificationRulesMixin : public BASE_MSC {
	public:
		typedef BASE_MSC Base;
		typedef typename BASE_MSC::ScalarType SCALAR_TYPE;
		typedef typename BASE_MSC::MeshType MESH_TYPE;
		typedef typename BASE_MSC::FuncType FUNC_TYPE;
		typedef typename BASE_MSC::GradType GRAD_TYPE;
		typedef STAT_TYPE StatType;

		SimplificationRulesMixin(GRAD_TYPE* grad, MESH_TYPE* mesh, FUNC_TYPE* func,
			unsigned restrictMask, const MergeRuleSet<STAT_TYPE>& rules) :
			Base(grad, mesh, func), mRestrictMask(restrictMask), mRules(rules) {}

		// Build base statistics by tracing every restricted extremum's
		// unsimplified manifold. Call after ComputeFromGrad, before
		// ComputeHierarchy.
		void PrimeStatistics() {
			const DIM_TYPE maxdim = this->mMesh->maxDim();
			std::vector<INT_TYPE> todo;
			for (INT_TYPE i = 0; i < (INT_TYPE)this->nodes.size(); ++i)
				if (isRestricted(this->nodes[i].dim, maxdim)) todo.push_back(i);
			const long long n = (long long)todo.size();
			// A partition-local MSC primes from inside a worker thread: stay
			// serial there rather than spawning nested teams.
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic, 16) if(!omp_in_parallel())
#endif
			for (long long k = 0; k < n; ++k) {
				const INT_TYPE nid = todo[(size_t)k];
				const node<SCALAR_TYPE>& nd = this->nodes[nid];
				const bool ascending = (nd.dim == 0);
				const INT_TYPE man = ascending ? nd.amanifoldid : nd.dmanifoldid;
				std::vector<INDEX_TYPE> seeds;
				baseSeedCells(nid, ascending, seeds);
				STAT_TYPE s;
				for (size_t j = 0; j < seeds.size(); ++j) traceBaseManifold(seeds[j], ascending, s);
				mManStat[(size_t)man] = s;
			}
			mPrimed = true;
		}

		// Existing labeling: labels[cell] = region id (< 0 = none) over mesh
		// cell ids, regionToExtremumCell[region] = that region's extremum
		// cell. Only cells of the extremum's dimension are accumulated.
		// Call once per restricted direction (ascending = minima).
		void PrimeStatisticsFromLabeling(bool ascending, const DenseLabeling<int>& labels,
			const std::vector<INDEX_TYPE>& regionToExtremumCell) {
			const DIM_TYPE want = ascending ? 0 : this->mMesh->maxDim();
			std::vector<INT_TYPE> regionMan(regionToExtremumCell.size(), -1);
			for (size_t r = 0; r < regionToExtremumCell.size(); ++r) {
				typename std::unordered_map<INDEX_TYPE, INT_TYPE>::const_iterator it =
					this->nodemap.find(regionToExtremumCell[r]);
				if (it == this->nodemap.end()) continue;
				const node<SCALAR_TYPE>& nd = this->nodes[it->second];
				regionMan[r] = ascending ? nd.amanifoldid : nd.dmanifoldid;
				mManStat[(size_t)regionMan[r]] = STAT_TYPE();
			}
			for (INDEX_TYPE c = 0; c < labels.GetNumLabels(); ++c) {
				const int r = labels.GetLabel(c);
				if (r < 0 || (size_t)r >= regionMan.size() || regionMan[(size_t)r] < 0) continue;
				if (this->mMesh->dimension(c) != want) continue;
				mManStat[(size_t)regionMan[(size_t)r]].AddCell(c);
			}
			mPrimed = true;
		}

		bool IsPrimed() const { return mPrimed; }

		// Statistic of a node's CURRENT manifold (at the last cancellation).
		const STAT_TYPE& NodeStatistic(INT_TYPE nodeId, bool ascending) const {
			const node<SCALAR_TYPE>& nd = this->nodes[nodeId];
			return mManStat[(size_t)(ascending ? nd.amanifoldid : nd.dmanifoldid)];
		}

		long long NumVetoed() const { return mNumVetoed; }

		void ComputeHierarchy(SCALAR_TYPE pers_limit) override {
			if (!mPrimed) PrimeStatistics();
			mNumVetoed = 0;
			Base::ComputeHierarchy(pers_limit);
			printf("  -- Simplification rules vetoed %lld cancellations\n", mNumVetoed);
		}

	protected:
		// Extremum cells whose unsimplified manifolds together form node's
		// base region. Default: the node's own cell.
		virtual void baseSeedCells(INT_TYPE nodeId, bool ascending, std::vector<INDEX_TYPE>& out) const {
			(void)ascending;
			out.assign(1, this->nodes[nodeId].cellindex);
		}

		void onManifoldCreated(INT_TYPE manId) override {
			const merged_manifold& m = this->mans[(size_t)manId];
			STAT_TYPE s;
			if (m.merged[0] >= 0 && m.merged[1] >= 0)
				s = STAT_TYPE::Merge(mManStat[(size_t)m.merged[0]], mManStat[(size_t)m.merged[1]]);
			else if (m.merged[0] >= 0)
				s = mManStat[(size_t)m.merged[0]];
			mManStat.push_back(s);
		}

		bool isValid(INT_TYPE a, arc<SCALAR_TYPE>& ap) const override {
			if (!Base::isValid(a, ap)) return false;
			if (!mPrimed || mRules.Empty()) return true;
			const DIM_TYPE maxdim = this->mMesh->maxDim();
			const node<SCALAR_TYPE>& lo = this->nodes[ap.lower];
			const node<SCALAR_TYPE>& up = this->nodes[ap.upper];

			// minimum dies: its ascending manifold joins that of every other
			// minimum below the saddle (cancel()'s ascending update)
			if ((mRestrictMask & RESTRICT_MINIMA_ASC) && lo.dim == 0) {
				const STAT_TYPE& dying = mManStat[(size_t)lo.amanifoldid];
				INT_TYPE b = up.firstarc;
				while (b != NULLID) {
					const arc<SCALAR_TYPE>& bp = this->arcs[b];
					if (bp.upper == ap.upper && bp.lower != ap.lower &&
						this->isAlive(bp, this->num_cancelled) &&
						!mRules.Allows(dying, mManStat[(size_t)this->nodes[bp.lower].amanifoldid])) {
						++mNumVetoed;
						return false;
					}
					b = this->nextArc(bp, ap.upper);
				}
			}
			// maximum dies: its descending manifold joins that of every other
			// maximum above the saddle
			if ((mRestrictMask & RESTRICT_MAXIMA_DSC) && up.dim == maxdim) {
				const STAT_TYPE& dying = mManStat[(size_t)up.dmanifoldid];
				INT_TYPE b = lo.firstarc;
				while (b != NULLID) {
					const arc<SCALAR_TYPE>& bp = this->arcs[b];
					if (bp.lower == ap.lower && bp.upper != ap.upper &&
						this->isAlive(bp, this->num_cancelled) &&
						!mRules.Allows(dying, mManStat[(size_t)this->nodes[bp.upper].dmanifoldid])) {
						++mNumVetoed;
						return false;
					}
					b = this->nextArc(bp, ap.lower);
				}
			}
			return true;
		}

	private:
		bool isRestricted(DIM_TYPE dim, DIM_TYPE maxdim) const {
			return ((mRestrictMask & RESTRICT_MINIMA_ASC) && dim == 0) ||
				((mRestrictMask & RESTRICT_MAXIMA_DSC) && dim == maxdim);
		}

		// Iterative form of rec_man_trace_up/down (base manifolds can be far
		// deeper than an OpenMP worker's stack), feeding AddCell only the
		// cells of the extremum's own dimension.
		void traceBaseManifold(INDEX_TYPE start, bool ascending, STAT_TYPE& s) const {
			std::vector<INDEX_TYPE> stack(1, start);
			while (!stack.empty()) {
				const INDEX_TYPE cur = stack.back();
				stack.pop_back();
				const DIM_TYPE cdim = this->mMesh->dimension(cur);
				s.AddCell(cur);
				if (ascending) {
					typename MESH_TYPE::CofacetsIterator it(this->mMesh);
					for (it.begin(cur); it.valid(); it.advance()) {
						const INDEX_TYPE t = it.value();
						if (this->mGrad->getCritical(t)) continue;
						const INDEX_TYPE p = this->mGrad->getPair(t);
						if (p == cur || this->mMesh->dimension(p) != cdim) continue;
						stack.push_back(p);
					}
				}
				else {
					typename MESH_TYPE::FacetsIterator it(this->mMesh);
					for (it.begin(cur); it.valid(); it.advance()) {
						const INDEX_TYPE t = it.value();
						if (this->mGrad->getCritical(t)) continue;
						const INDEX_TYPE p = this->mGrad->getPair(t);
						if (p == cur || this->mMesh->dimension(p) != cdim) continue;
						stack.push_back(p);
					}
				}
			}
		}

		unsigned mRestrictMask;
		MergeRuleSet<STAT_TYPE> mRules;
		std::vector<STAT_TYPE> mManStat;  // parallel to mans
		bool mPrimed = false;
		mutable long long mNumVetoed = 0;
	};

	template<typename SCALAR_TYPE, class MESH_TYPE, class FUNC_TYPE, class GRAD_TYPE, class STAT_TYPE>
	using MorseSmaleComplexWithRules = SimplificationRulesMixin<
		MorseSmaleComplexBasic<SCALAR_TYPE, MESH_TYPE, FUNC_TYPE, GRAD_TYPE>, STAT_TYPE>;

} // namespace GInt

#endif
