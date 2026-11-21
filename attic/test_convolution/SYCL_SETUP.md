# SYCL GPU Backend Configuration for test_convolution

This document describes the SYCL configuration fixes and how to set up the SYCL + GPU backend for building test_convolution.

## What Was Fixed

### 1. CMakeLists.txt Improvements
- **Enhanced SYCL Detection**: Now automatically searches for Intel DPC++ compiler (icpx/dpcpp) in multiple standard installation locations
- **Windows Environment Variable Handling**: Fixed issues with Windows environment variables containing parentheses (`ProgramFiles(x86)`)
- **Multiple SYCL Implementation Support**: Supports Intel oneAPI DPC++, hipSYCL, and other SYCL implementations
- **Better Error Messages**: Provides clear instructions when SYCL is not found

### 2. SYCL Header Compatibility
- **Dual Header Support**: Updated `execution_strategy_sycl.h` to support both:
  - Old SYCL 1.2.1 format: `<CL/sycl.hpp>`
  - Modern DPC++ format: `<sycl/sycl.hpp>`
- Uses `__has_include` for automatic detection

### 3. GPU Device Selection
- **Improved Device Selection**: Enhanced GPU selector with fallback to default device if GPU is unavailable
- **Better Device Information**: Prints detailed device information including vendor, name, and type
- **Error Handling**: Gracefully handles cases where GPU is not available

### 4. Compiler Flags for GPU Offloading
- **Intel GPU Targeting**: Configures `-fsycl-targets=spir64_gen` for Intel GPU support
- **Device Backend**: Sets `-device gen12` for Intel integrated/discrete GPUs

## Setup Instructions

### Option 1: Install Intel oneAPI Base Toolkit (Recommended)

1. Download and install Intel oneAPI Base Toolkit:
   ```
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html
   ```

2. After installation, configure your build:
   ```powershell
   cd test_convolution
   mkdir build
   cd build
   cmake -DENABLE_SYCL=ON -DENABLE_OPENMP=OFF ..
   cmake --build .
   ```

   The CMake script will automatically detect the Intel DPC++ compiler.

### Option 2: Set Environment Variable

If Intel oneAPI is installed but not detected:

```powershell
$env:ONEAPI_ROOT = "C:\Program Files (x86)\Intel\oneAPI"
cd test_convolution\build
cmake -DENABLE_SYCL=ON ..
```

### Option 3: Specify Compiler Explicitly

```powershell
cd test_convolution\build
cmake -DENABLE_SYCL=ON -DCMAKE_CXX_COMPILER="C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin\icpx.exe" ..
```

### Option 4: Add Compiler to PATH

```powershell
$env:Path += ";C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin"
cd test_convolution\build
cmake -DENABLE_SYCL=ON ..
```

## Diagnostic Script

A PowerShell diagnostic script is available to check your SYCL setup:

```powershell
cd test_convolution
.\check_sycl_setup.ps1
```

This script will:
- Check for Intel oneAPI installation
- Look for SYCL compilers in PATH
- Test CMake configuration
- Provide setup recommendations

## Build Configuration

### Enable SYCL
```powershell
cmake -DENABLE_SYCL=ON ..
```

### Build with SYCL
```powershell
cmake --build . --config Release
```

### Run with SYCL Execution
```powershell
.\test_convolution.exe input.raw output.raw 512 512 512 sycl
```

## Troubleshooting

### SYCL Not Found
If CMake reports "SYCL requested but not found":
1. Verify Intel oneAPI is installed
2. Run the diagnostic script: `.\check_sycl_setup.ps1`
3. Check that `icpx` or `dpcpp` is accessible
4. Set `ONEAPI_ROOT` environment variable if needed

### GPU Device Not Available
If the program fails to find a GPU:
- The code will automatically fall back to the default device
- Check that your system has a compatible GPU (Intel, NVIDIA, AMD)
- Verify GPU drivers are installed
- For Intel GPUs, ensure you have Intel GPU drivers installed

### Compilation Errors
If you encounter compilation errors:
- Ensure you're using Intel DPC++ (icpx) or another SYCL-capable compiler
- Check that `GINT_USE_SYCL` is defined (CMake should handle this)
- Verify SYCL headers are accessible

## Supported SYCL Implementations

The configuration supports:
- **Intel oneAPI DPC++** (recommended for Intel GPUs)
- **hipSYCL** (alternative implementation)
- Other SYCL implementations via find_package(SYCL)

## GPU Targets

Currently configured for:
- **Intel GPU**: Uses `spir64_gen` target with `gen12` device backend
- Can be modified in CMakeLists.txt for other GPU vendors

To target NVIDIA GPUs (with hipSYCL), modify the CMakeLists.txt to use appropriate flags.

## Notes

- The SYCL backend requires a GPU-capable device or will fall back to CPU
- OpenMP and SYCL can be enabled simultaneously, but the execution target must be specified at runtime
- SYCL execution is only available if compiled with `ENABLE_SYCL=ON`
