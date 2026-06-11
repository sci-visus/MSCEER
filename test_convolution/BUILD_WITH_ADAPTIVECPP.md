# Building test_convolution with AdaptiveCpp

## Recommended Setup: Intel oneAPI DPC++ with Ninja

AdaptiveCpp works best with clang-based compilers. Intel oneAPI DPC++ (icpx) is recommended for Windows.

### Option 1: Using the Build Script (Easiest)

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\code\MSCEER\test_convolution
.\build_with_adaptivecpp.ps1
```

The script will automatically:
- Detect Intel oneAPI if installed
- Use Ninja if available, otherwise Visual Studio
- Configure AdaptiveCpp correctly

### Option 2: Manual Build with Intel oneAPI and Ninja

#### Step 1: Set up oneAPI Environment

```powershell
# Source oneAPI environment (adjust path if needed)
& "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
```

Or if you have it in a different location:
```powershell
& "C:\Program Files\Intel\oneAPI\setvars.bat"
```

#### Step 2: Build with Ninja

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\code\MSCEER\test_convolution
mkdir build_acpp_ninja
cd build_acpp_ninja

# Configure with Intel DPC++ and Ninja
# Note: Ensure ninja.exe is in PATH or specify CMAKE_MAKE_PROGRAM
cmake -G "Ninja" `
      -DCMAKE_CXX_COMPILER="C:\Program Files (x86)\Intel\oneAPI\2025.3\bin\icpx.exe" `
      -DCMAKE_MAKE_PROGRAM="C:\Users\jediati\bin\ninja.exe" `
      -DENABLE_SYCL=ON `
      -DPREFER_ADAPTIVECPP=ON `
      -DCMAKE_PREFIX_PATH="C:\Users\jediati\bin\AdaptiveCpp-LLVM20-Win" `
      -DCMAKE_BUILD_TYPE=Release `
      ..

# Build
cmake --build . --config Release
```

### Option 3: Manual Build with Visual Studio Generator (MSVC)

If you don't have Intel oneAPI or prefer MSVC:

```powershell
cd C:\Users\jediati\Desktop\JEDIATI\code\MSCEER\test_convolution
mkdir build_acpp
cd build_acpp

cmake -G "Visual Studio 18 2026" -A x64 `
      -DENABLE_SYCL=ON `
      -DPREFER_ADAPTIVECPP=ON `
      -DCMAKE_PREFIX_PATH="C:\Users\jediati\bin\AdaptiveCpp-LLVM20-Win" `
      ..

cmake --build . --config Release
```

**Note:** MSVC may have limitations with SYCL. Intel DPC++ is recommended.

## Installing Required Tools

### Install Ninja (if not available)

**Option A: Using Chocolatey**
```powershell
choco install ninja
```

**Option B: Manual Download**
1. Download from: https://github.com/ninja-build/ninja/releases
2. Extract `ninja.exe` to a directory in your PATH (e.g., `C:\Program Files\CMake\bin\`)

### Install Intel oneAPI (if not available)

1. Download Intel oneAPI Base Toolkit from:
   https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html

2. Install and note the installation path (typically `C:\Program Files (x86)\Intel\oneAPI\`)

3. Source the environment before building:
   ```powershell
   & "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
   ```

## Verifying the Setup

Run the diagnostic script:
```powershell
cd C:\Users\jediati\Desktop\JEDIATI\code\MSCEER\test_convolution
.\check_adaptivecpp_setup.ps1
```

## Troubleshooting

### CMake can't find AdaptiveCpp
- Ensure `CMAKE_PREFIX_PATH` points to the AdaptiveCpp installation directory
- Check that `lib\cmake\AdaptiveCpp\adaptivecpp-config.cmake` exists

### Compiler not found
- For Intel DPC++: Source oneAPI setvars.bat before running CMake
- Or provide full path: `-DCMAKE_CXX_COMPILER="C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin\icpx.exe"`

### Build errors with MSVC
- Try using Intel DPC++ instead: `-DCMAKE_CXX_COMPILER=icpx`
- Or ensure you're using a recent Visual Studio version

### Ninja not found
- Install Ninja (see above)
- Or use Visual Studio generator: `-G "Visual Studio 18 2026"`
