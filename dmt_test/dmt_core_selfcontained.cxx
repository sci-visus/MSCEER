// Compiles every include/dmt_core header with ONLY include/dmt_core on the
// include path (see CMakeLists.txt), so a core header that reaches for a gi_*
// header, or for anything other than <cstdint>, fails to build here.
//
// Exit code: 0 on success.

#include "dmt_hd.h"

#include <cstdio>

int main() {
	int bad = 0;
	if (dmt::GradPairCode(dmt::MakeGradByte(5)) != 5) bad++;
	if (!dmt::GradAssigned(dmt::MakeGradByte(0))) bad++;
	if (!dmt::GradCritical(dmt::MakeGradByte(dmt::kPairCritical))) bad++;
	if (dmt::MakeGradByte(dmt::kPairCritical) != (1u | (7u << 2))) bad++;
	printf("dmt_core_selfcontained: %s\n", bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
