#ifndef GI_MORSE_SMALE_COMPLEX_PARTITIONED_RULES_H
#define GI_MORSE_SMALE_COMPLEX_PARTITIONED_RULES_H

#include <vector>

#include "gi_morse_smale_complex_partitioned.h"
#include "gi_morse_smale_complex_rules.h"

namespace GInt {

	// Custom simplification rules on the partitioned MSC pipeline.
	//
	// Local phase: SimplificationRulesMixin<PartitionLocalMsc, STAT>. Each
	// local MSC primes by tracing its own extrema through the GLOBAL gradient,
	// so region statistics are exact even where a manifold leaves the
	// partition, and every local cancellation is checked against them.
	//
	// Reconciled phase: RulesReconciledMsc below. A reconciled base node
	// stands for a locally merged region; its base statistic is the sum over
	// the imported lineage cells, so the global hierarchy continues from the
	// true region sizes.
	//
	// The rule therefore holds for every applied cancellation. What differs
	// from a serial build is the ORDER cancellations are considered in (local
	// phases run ahead of the global one), and so possibly which merges win -
	// the same way partitioned output already differs from serial.
	template<class PARTITIONED, class STAT_TYPE>
	class RulesReconciledMsc
		: public SimplificationRulesMixin<typename PARTITIONED::ReconciledGlobalMsc, STAT_TYPE> {
	public:
		typedef SimplificationRulesMixin<typename PARTITIONED::ReconciledGlobalMsc, STAT_TYPE> Mixin;

		RulesReconciledMsc(typename Mixin::GRAD_TYPE* grad, typename Mixin::MESH_TYPE* mesh,
			typename Mixin::FUNC_TYPE* func, unsigned restrictMask, const MergeRuleSet<STAT_TYPE>& rules) :
			Mixin(grad, mesh, func, restrictMask, rules) {}

	protected:
		void baseSeedCells(INT_TYPE nodeId, bool ascending, std::vector<INDEX_TYPE>& out) const override {
			const INDEX_TYPE rep = this->nodes[nodeId].cellindex;
			const bool found = ascending ? this->GetAscendingLineageCells(rep, out)
				: this->GetDescendingLineageCells(rep, out);
			if (!found || out.empty()) out.assign(1, rep);
		}
	};

	// Install rule-aware local and reconciled MSCs on a partitioned pipeline.
	template<class PARTITIONED, class STAT_TYPE>
	void InstallSimplificationRules(PARTITIONED& pipeline, unsigned restrictMask,
		const MergeRuleSet<STAT_TYPE>& rules) {
		typedef typename PARTITIONED::PartitionLocalMsc LocalType;
		typedef typename PARTITIONED::ReconciledGlobalMsc ReconciledType;
		typedef SimplificationRulesMixin<LocalType, STAT_TYPE> LocalRules;
		typedef RulesReconciledMsc<PARTITIONED, STAT_TYPE> ReconciledRules;
		pipeline.SetMscFactories(
			[restrictMask, rules](typename LocalRules::GRAD_TYPE* g, typename LocalRules::MESH_TYPE* m,
				typename LocalRules::FUNC_TYPE* f) -> LocalType* {
				return new LocalRules(g, m, f, restrictMask, rules);
			},
			[restrictMask, rules](typename LocalRules::GRAD_TYPE* g, typename LocalRules::MESH_TYPE* m,
				typename LocalRules::FUNC_TYPE* f) -> ReconciledType* {
				return new ReconciledRules(g, m, f, restrictMask, rules);
			});
	}

} // namespace GInt

#endif
