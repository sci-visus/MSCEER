// The method IDs of gradient/gi_gradient_methods.h: names, parsing and the
// aligned-first helpers.
//
// Exit code: 0 on success.

#include "gradient/gi_gradient_methods.h"

#include <cstdio>
#include <cstring>

int main() {
	using GInt::GradientMethod;
	int bad = 0;
	for (int i = 0; i < GInt::kNumGradientMethods; i++) {
		const GradientMethod m = (GradientMethod)i;
		const char* name = GInt::MethodName(m);
		GradientMethod back = GradientMethod::PLS;
		if (!name[0] || !GInt::ParseMethod(name, back) || back != m) { printf("round trip failed for %d\n", i); bad++; }
		// an aligned-first ID is its base ID plus "_AF"
		const GradientMethod base = GInt::BaseMethod(m);
		if (GInt::IsAlignedFirst(m)) {
			char want[32];
			std::snprintf(want, sizeof(want), "%s_AF", GInt::MethodName(base));
			if (std::strcmp(want, name) != 0) { printf("%s is not %s\n", name, want); bad++; }
			GradientMethod af;
			if (!GInt::WithAlignedFirst(base, true, af) || af != m) { printf("WithAlignedFirst failed for %s\n", name); bad++; }
		}
		else if (base != m) { printf("BaseMethod changed %s\n", name); bad++; }
		GradientMethod same;
		if (!GInt::WithAlignedFirst(m, false, same) || same != m) bad++;
	}
	GradientMethod x = GradientMethod::CLS;
	if (GInt::ParseMethod("pls", x) || GInt::ParseMethod("PLS_", x) || GInt::ParseMethod("", x) ||
		GInt::ParseMethod(nullptr, x) || x != GradientMethod::CLS) { printf("ParseMethod accepted a non-ID\n"); bad++; }
	if (GInt::WithAlignedFirst(GradientMethod::DMF, true, x) || GInt::WithAlignedFirst(GradientMethod::PLS_BC, true, x)) {
		printf("WithAlignedFirst accepted a method with no aligned-first form\n"); bad++;
	}
	printf("gradient_methods_smoke: %d failures -> %s\n", bad, bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
