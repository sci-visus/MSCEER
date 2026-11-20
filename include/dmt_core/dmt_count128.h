#ifndef DMT_COUNT128_H
#define DMT_COUNT128_H

// 128-bit path counts with a sticky overflow flag.
//
// Gradient-path counts grow exponentially on some fields, so 64 bits are not
// enough; MSVC has no __int128, so these are two 64-bit limbs with explicit
// carries (portable, also on device). U128 is unsigned; S128 is two's
// complement. Every operation that cannot represent its exact result sets
// `overflow` on the result and keeps it set through later operations, so a
// count is either exact or flagged -- never silently wrapped.

#include "dmt_hd.h"

namespace dmt {

	struct U128 {
		uint64_t lo = 0, hi = 0;
		bool overflow = false;

		DMT_HD U128() {}
		DMT_HD explicit U128(uint64_t v) : lo(v), hi(0), overflow(false) {}
		DMT_HD bool IsZero() const { return lo == 0 && hi == 0; }
		DMT_HD bool FitsU64() const { return hi == 0 && !overflow; }
	};

	DMT_HD U128 Add(const U128& a, const U128& b) {
		U128 r;
		r.lo = a.lo + b.lo;
		const uint64_t c0 = r.lo < a.lo ? 1u : 0u;
		const uint64_t h1 = a.hi + b.hi;
		const bool c1 = h1 < a.hi;
		r.hi = h1 + c0;
		const bool c2 = r.hi < h1;
		r.overflow = a.overflow || b.overflow || c1 || c2;
		return r;
	}

	DMT_HD bool Less(const U128& a, const U128& b) {
		return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo;
	}

	struct S128 {
		uint64_t lo = 0, hi = 0;   // two's complement
		bool overflow = false;

		DMT_HD S128() {}
		DMT_HD static S128 FromInt(int64_t v) {
			S128 r;
			r.lo = (uint64_t)v;
			r.hi = v < 0 ? ~(uint64_t)0 : 0;
			return r;
		}
		DMT_HD bool Negative() const { return (hi >> 63) != 0; }
		DMT_HD bool IsZero() const { return lo == 0 && hi == 0; }
		DMT_HD bool FitsI64() const {
			if (overflow) return false;
			return Negative() ? (hi == ~(uint64_t)0 && (lo >> 63) == 1) : (hi == 0 && (lo >> 63) == 0);
		}
		DMT_HD int64_t ToI64() const { return (int64_t)lo; }
	};

	DMT_HD S128 Add(const S128& a, const S128& b) {
		S128 r;
		r.lo = a.lo + b.lo;
		const uint64_t c0 = r.lo < a.lo ? 1u : 0u;
		r.hi = a.hi + b.hi + c0;
		const bool sa = a.Negative(), sb = b.Negative(), sr = r.Negative();
		r.overflow = a.overflow || b.overflow || (sa == sb && sr != sa);
		return r;
	}

	DMT_HD S128 Neg(const S128& a) {
		S128 r;
		r.lo = ~a.lo + 1;
		r.hi = ~a.hi + (r.lo == 0 ? 1u : 0u);
		// -INT128_MIN is not representable
		r.overflow = a.overflow || (a.hi == ((uint64_t)1 << 63) && a.lo == 0);
		return r;
	}

	DMT_HD S128 Sub(const S128& a, const S128& b) { return Add(a, Neg(b)); }

	// a times a sign (+1 or -1).
	DMT_HD S128 MulSign(const S128& a, int sign) { return sign >= 0 ? a : Neg(a); }

	// Unsigned value as signed; overflows if it exceeds 2^127 - 1.
	DMT_HD S128 ToSigned(const U128& u) {
		S128 r;
		r.lo = u.lo;
		r.hi = u.hi;
		r.overflow = u.overflow || (u.hi >> 63) != 0;
		return r;
	}

	// |a| as unsigned (exact for every S128 value).
	DMT_HD U128 Abs(const S128& a) {
		U128 r;
		if (a.Negative()) {
			r.lo = ~a.lo + 1;
			r.hi = ~a.hi + (r.lo == 0 ? 1u : 0u);
		}
		else {
			r.lo = a.lo;
			r.hi = a.hi;
		}
		r.overflow = a.overflow;
		return r;
	}

	// 64 x 64 -> 128 unsigned, portable.
	DMT_HD void MulU64(uint64_t a, uint64_t b, uint64_t& lo, uint64_t& hi) {
		const uint64_t a0 = a & 0xFFFFFFFFull, a1 = a >> 32;
		const uint64_t b0 = b & 0xFFFFFFFFull, b1 = b >> 32;
		const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
		const uint64_t mid = (p00 >> 32) + (p01 & 0xFFFFFFFFull) + (p10 & 0xFFFFFFFFull);
		lo = (mid << 32) | (p00 & 0xFFFFFFFFull);
		hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
	}

	// Signed product, exact or flagged.
	DMT_HD S128 Mul(const S128& a, const S128& b) {
		const U128 ma = Abs(a), mb = Abs(b);
		S128 r;
		r.overflow = a.overflow || b.overflow;
		if (ma.hi != 0 && mb.hi != 0) { r.overflow = true; return r; }
		uint64_t lo, hi;
		// (ma.hi*2^64 + ma.lo) * mb.lo, or the symmetric case
		const U128& big = ma.hi != 0 ? ma : mb;
		const U128& small = ma.hi != 0 ? mb : ma;
		MulU64(big.lo, small.lo, lo, hi);
		uint64_t cross_lo, cross_hi;
		MulU64(big.hi, small.lo, cross_lo, cross_hi);
		if (cross_hi != 0) r.overflow = true;
		const uint64_t h = hi + cross_lo;
		if (h < hi) r.overflow = true;
		S128 m;
		m.lo = lo;
		m.hi = h;
		if ((m.hi >> 63) != 0) r.overflow = true;   // magnitude >= 2^127
		const bool neg = a.Negative() != b.Negative();
		S128 out = neg ? Neg(m) : m;
		out.overflow = out.overflow || r.overflow;
		return out;
	}

} // namespace dmt

#endif
