# Fixing AdaptiveCpp Build on Windows

## Problem

AdaptiveCpp build fails with:
```
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.
hipSYCL_OpenMP_libomp_LIBRARY
```

This occurs because AdaptiveCpp's OpenMP backend requires the LLVM OpenMP library (`libomp`), which is not found.

## Solutions

### Option 1: Disable OpenMP Backend (Recommended if you only need SYCL)

If you only need SYCL support and don't need the OpenMP backend, disable it:

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp\build
cmake -DCMAKE_INSTALL_PREFIX=C:\Users\jediati\bin ^
      -DWITH_CPU_BACKEND=OFF ^
      -DDEFAULT_TARGETS="cuda:sm_75" ..
```

Note: Fix the install prefix typo (remove the double `C:`).

### Option 2: Build LLVM OpenMP Library

If you need the OpenMP backend, you need to build LLVM's OpenMP library:

1. **Download and build LLVM OpenMP:**
   ```powershell
   git clone https://github.com/llvm/llvm-project.git
   cd llvm-project\openmp
   mkdir build
   cd build
   cmake -G "Visual Studio 18 2026" -A x64 -DCMAKE_INSTALL_PREFIX=C:\llvm-openmp ..
   cmake --build . --config Release
   cmake --install . --config Release
   ```

2. **Configure AdaptiveCpp with OpenMP path:**
   ```powershell
   cd C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp\build
   cmake -DCMAKE_INSTALL_PREFIX=C:\Users\jediati\bin ^
         -DLLVM_OpenMP_DIR=C:\llvm-openmp\lib\cmake\openmp ^
         -DACPP_TARGETS="cuda:sm_75;cuda:sm_86;cuda:sm_89" ..
   ```

### Option 3: Use Pre-built LLVM (Easiest)

Download a pre-built LLVM distribution that includes OpenMP:

1. Download LLVM from: https://github.com/llvm/llvm-project/releases
   - Look for Windows x64 builds
   - Extract to `C:\llvm` or similar

2. Configure AdaptiveCpp:
   ```powershell
   cd C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp\build
   cmake -DCMAKE_INSTALL_PREFIX=C:\Users\jediati\bin ^
         -DLLVM_DIR=C:\llvm\lib\cmake\llvm ^
         -DLLVM_OpenMP_DIR=C:\llvm\lib\cmake\openmp ^
         -DACPP_TARGETS="cuda:sm_75;cuda:sm_86;cuda:sm_89" ..
   ```

### Option 4: Disable OpenMP Backend and Use CUDA Only

Since you have CUDA 12.1 installed, you can build AdaptiveCpp with CUDA support only:

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp\build
cmake -DCMAKE_INSTALL_PREFIX=C:\Users\jediati\bin ^
      -DENABLE_OPENMP_BACKEND=OFF ^
      -DACPP_TARGETS="cuda:sm_75;cuda:sm_86;cuda:sm_89" ..
```

## Recommended Build Command (CUDA Only)

Based on your setup (CUDA 12.1, Visual Studio 2026), here's the recommended command:

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp\build
cmake -DCMAKE_INSTALL_PREFIX=C:\Users\jediati\bin ^
      -DWITH_CPU_BACKEND=OFF ^
      -DDEFAULT_TARGETS="cuda:sm_75" ..
cmake --build . --config Release
cmake --install . --config Release
```

**Note:** 
- Fix the install prefix: use `C:\Users\jediati\bin` not `C:C:\Users\jediati\bin`
- Replace `sm_75`, `sm_86`, `sm_89` with your actual GPU compute capability
- To find your GPU compute capability, run: `nvidia-smi --query-gpu=compute_cap --format=csv`

## After Building AdaptiveCpp

Once AdaptiveCpp is built and installed, configure your test_convolution project:

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\code\MSCEER\test_convolution
mkdir build_acpp
cd build_acpp
cmake -DENABLE_SYCL=ON ^
      -DPREFER_ADAPTIVECPP=ON ^
      -DCMAKE_PREFIX_PATH=C:\Users\jediati\bin ..
cmake --build . --config Release
```

## Troubleshooting

### If AdaptiveCpp still can't find OpenMP:
- Set `-DWITH_CPU_BACKEND=OFF` to disable the CPU/OpenMP backend (recommended if only using CUDA)
- Or provide explicit path: `-DLLVM_OpenMP_DIR=<path-to-openmp>` and `-DOpenMP_ROOT=<path-to-openmp>`
- Or manually set: `-DhipSYCL_OpenMP_libomp_LIBRARY=<path-to-libomp.lib>`

### If CUDA targets fail:
- Verify CUDA is in PATH: `nvcc --version`
- Check GPU compute capability matches your targets
- Try with a single target first: `-DACPP_TARGETS="cuda:sm_75"`

### If CMake can't find AdaptiveCpp after installation:
- Set `CMAKE_PREFIX_PATH` to the install directory
- Or set `ACPP_ROOT` environment variable
