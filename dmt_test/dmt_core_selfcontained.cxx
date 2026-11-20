// Compiles every include/dmt_core header with ONLY include/dmt_core on the
// include path (see CMakeLists.txt), so a core header that reaches for a gi_*
// header, or for anything other than <cstdint>, fails to build here.
//
// Exit code: 0 on success.

#include "dmt_hd.h"
#include "dmt_cubical3d.h"
#include "dmt_star27_tables.h"
#include "dmt_lowerstar_pq.h"
#include "dmt_dmf_cell.h"
#include "dmt_dmf_keycache.h"

#include <cstdio>

int main() {
	int bad = 0;
	if (dmt::GradPairCode(dmt::MakeGradByte(5)) != 5) bad++;
	if (!dmt::GradAssigned(dmt::MakeGradByte(0))) bad++;
	if (!dmt::GradCritical(dmt::MakeGradByte(dmt::kPairCritical))) bad++;
	if (dmt::MakeGradByte(dmt::kPairCritical) != (1u | (7u << 2))) bad++;

	// One star through the core alone.
	dmt::Star27Tables T;
	dmt::BuildStar27TablesCanonical(T);
	const dmt::Cubical3D G = dmt::Cubical3D::Make(3, 3, 3);
	uint32_t ranks[27], key[27], vvalid;
	for (int i = 0; i < 27; i++) ranks[i] = (uint32_t)((i * 7) % 27);
	dmt::GatherStarKeys(G, ranks, 1, 1, 1, key, vvalid);
	dmt::StarResult R;
	dmt::PLSStar(T, key, vvalid, R);
	if (R.pair[dmt::kStarCentre] == -2) bad++;

	printf("dmt_core_selfcontained: %s\n", bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
