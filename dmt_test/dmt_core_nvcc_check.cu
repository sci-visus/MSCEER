// The device-safe core compiled by nvcc and run on the GPU.
//
// One thread per vertex runs PLSStar / PLSBRStar from dmt_core and writes the
// GradBitfield bytes of its lower star; the same core functions, compiled for
// the host, produce the reference. Only include/dmt_core is on the include
// path (see CMakeLists.txt), so this also proves the core needs nothing else.
//
// Exit code: 0 iff the GPU bytes equal the host bytes on every field, or 2 if
// no CUDA device is available (the compile itself is then the check).

#include "dmt_hd.h"
#include "dmt_cubical3d.h"
#include "dmt_star27_tables.h"
#include "dmt_lowerstar_pq.h"

#include <cuda_runtime.h>

#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
#include <algorithm>

__device__ void WriteStarDevice(const dmt::Star27Tables& T, const dmt::Cubical3D& G,
	int64_t x, int64_t y, int64_t z, const dmt::StarResult& R, uint8_t* out) {
	for (int s = 0; s < 27; s++) {
		if (R.pair[s] == -2) continue;
		const int64_t id = G.Id(2 * x + T.off[s][0], 2 * y + T.off[s][1], 2 * z + T.off[s][2]);
		const int code = R.pair[s] < 0 ? dmt::kPairCritical : dmt::StarDir(T, s, R.pair[s]);
		out[id] = dmt::MakeGradByte((uint8_t)code);
	}
}

__global__ void LowerStarKernel(dmt::Star27Tables T, dmt::Cubical3D G, const uint32_t* ranks,
	uint8_t* out, int boundary_restricted) {
	const int64_t v = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
	if (v >= G.NumVertices()) return;
	const int64_t x = v % G.n[0], y = (v / G.n[0]) % G.n[1], z = v / (G.n[0] * G.n[1]);
	uint32_t key[27], vvalid;
	dmt::GatherStarKeys(G, ranks, x, y, z, key, vvalid);
	dmt::StarResult R;
	if (boundary_restricted)
		dmt::PLSBRStar(T, key, vvalid, dmt::StarBoundaryMask(T, G, x, y, z), R);
	else
		dmt::PLSStar(T, key, vvalid, R);
	WriteStarDevice(T, G, x, y, z, R, out);
}

static void HostReference(const dmt::Star27Tables& T, const dmt::Cubical3D& G, const uint32_t* ranks,
	uint8_t* out, int boundary_restricted) {
	for (int64_t v = 0; v < G.NumVertices(); v++) {
		const int64_t x = v % G.n[0], y = (v / G.n[0]) % G.n[1], z = v / (G.n[0] * G.n[1]);
		uint32_t key[27], vvalid;
		dmt::GatherStarKeys(G, ranks, x, y, z, key, vvalid);
		dmt::StarResult R;
		if (boundary_restricted)
			dmt::PLSBRStar(T, key, vvalid, dmt::StarBoundaryMask(T, G, x, y, z), R);
		else
			dmt::PLSStar(T, key, vvalid, R);
		for (int s = 0; s < 27; s++) {
			if (R.pair[s] == -2) continue;
			const int64_t id = G.Id(2 * x + T.off[s][0], 2 * y + T.off[s][1], 2 * z + T.off[s][2]);
			const int code = R.pair[s] < 0 ? dmt::kPairCritical : dmt::StarDir(T, s, R.pair[s]);
			out[id] = dmt::MakeGradByte((uint8_t)code);
		}
	}
}

int main() {
	int ndev = 0;
	const cudaError_t derr = cudaGetDeviceCount(&ndev);
	if (derr != cudaSuccess || ndev == 0) {
		printf("dmt_core_nvcc_check: compiled; no CUDA device to run on (%s)\n",
			derr != cudaSuccess ? cudaGetErrorString(derr) : "0 devices");
		return 2;
	}
	dmt::Star27Tables T;
	dmt::BuildStar27TablesCanonical(T);
	int bad = 0, fields = 0;
	std::mt19937 pick(7u);
	const int sizes[][3] = { {3, 3, 4}, {5, 6, 7}, {9, 8, 7}, {16, 17, 18}, {64, 48, 40}, {129, 97, 65} };
	for (int k = 0; k < 6 + 20; k++) {
		int X, Y, Z;
		if (k < 6) { X = sizes[k][0]; Y = sizes[k][1]; Z = sizes[k][2]; }
		else { X = 2 + (int)(pick() % 12); Y = 2 + (int)(pick() % 12); Z = 2 + (int)(pick() % 12); }
		const dmt::Cubical3D G = dmt::Cubical3D::Make(X, Y, Z);
		const int64_t nv = G.NumVertices(), nc = G.NumCells();
		std::vector<uint32_t> ranks((size_t)nv);
		for (int64_t i = 0; i < nv; i++) ranks[(size_t)i] = (uint32_t)i;
		std::mt19937 rng(100u + (unsigned)k);
		std::shuffle(ranks.begin(), ranks.end(), rng);

		uint32_t* d_ranks = nullptr;
		uint8_t* d_out = nullptr;
		cudaMalloc(&d_ranks, sizeof(uint32_t) * (size_t)nv);
		cudaMalloc(&d_out, (size_t)nc);
		cudaMemcpy(d_ranks, ranks.data(), sizeof(uint32_t) * (size_t)nv, cudaMemcpyHostToDevice);
		for (int rule = 0; rule < 2; rule++) {
			cudaMemset(d_out, 0, (size_t)nc);
			const int block = 128;
			const int grid = (int)((nv + block - 1) / block);
			LowerStarKernel<<<grid, block>>>(T, G, d_ranks, d_out, rule);
			cudaError_t err = cudaDeviceSynchronize();
			std::vector<uint8_t> gpu((size_t)nc), cpu((size_t)nc, 0);
			cudaMemcpy(gpu.data(), d_out, (size_t)nc, cudaMemcpyDeviceToHost);
			HostReference(T, G, ranks.data(), cpu.data(), rule);
			const bool same = err == cudaSuccess && std::memcmp(gpu.data(), cpu.data(), (size_t)nc) == 0;
			if (!same) {
				bad++;
				printf("  %s %dx%dx%d: GPU differs from host core (%s)\n", rule ? "PLS-BR" : "PLS", X, Y, Z,
					cudaGetErrorString(err));
			}
			fields++;
		}
		cudaFree(d_ranks);
		cudaFree(d_out);
	}
	printf("dmt_core_nvcc_check: %d runs, %d mismatches -> %s\n", fields, bad, bad ? "FAIL" : "ok");
	return bad ? 1 : 0;
}
