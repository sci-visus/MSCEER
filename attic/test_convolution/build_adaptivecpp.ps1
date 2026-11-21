# PowerShell script to build AdaptiveCpp with proper configuration
# Usage: .\build_adaptivecpp.ps1 [options]

param(
    [string]$AdaptiveCppSource = "C:\Users\jediati\Desktop\JEDIATI\libraries\AdaptiveCpp\AdaptiveCpp",
    [string]$InstallPrefix = "C:\Users\jediati\bin",
    [switch]$DisableOpenMP = $true,
    [string]$CudaTargets = "cuda:sm_75;cuda:sm_86;cuda:sm_89",
    [string]$LLVMOpenMPDir = ""
)

$ErrorActionPreference = "Stop"

Write-Host "=== Building AdaptiveCpp ===" -ForegroundColor Cyan
Write-Host ""

# Check if AdaptiveCpp source exists
if (-not (Test-Path $AdaptiveCppSource)) {
    Write-Host "Error: AdaptiveCpp source not found at: $AdaptiveCppSource" -ForegroundColor Red
    exit 1
}

# Check for CUDA
Write-Host "Checking CUDA installation..." -ForegroundColor Yellow
$cudaPath = Get-Command nvcc -ErrorAction SilentlyContinue
if ($cudaPath) {
    $cudaVersion = & nvcc --version 2>&1 | Select-String "release" | Select-Object -First 1
    Write-Host "  Found CUDA: $cudaVersion" -ForegroundColor Green
} else {
    Write-Host "  Warning: CUDA not found in PATH" -ForegroundColor Yellow
    Write-Host "    AdaptiveCpp CUDA backend may not work" -ForegroundColor Yellow
}

Write-Host ""

# Setup build directory
$buildDir = Join-Path $AdaptiveCppSource "build"
if (Test-Path $buildDir) {
    Write-Host "Cleaning existing build directory..." -ForegroundColor Yellow
    Remove-Item $buildDir -Recurse -Force
}
New-Item -ItemType Directory -Path $buildDir | Out-Null

Push-Location $buildDir

try {
    Write-Host "Configuring AdaptiveCpp..." -ForegroundColor Yellow
    Write-Host "  Source: $AdaptiveCppSource" -ForegroundColor Gray
    Write-Host "  Install: $InstallPrefix" -ForegroundColor Gray
    Write-Host "  OpenMP Backend: $(if ($DisableOpenMP) { 'DISABLED' } else { 'ENABLED' })" -ForegroundColor Gray
    Write-Host "  CUDA Targets: $CudaTargets" -ForegroundColor Gray
    Write-Host ""
    
    # Build CMake arguments
    $cmakeArgs = @(
        "-G", "Visual Studio 18 2026",
        "-A", "x64",
        "-DCMAKE_INSTALL_PREFIX=$InstallPrefix",
        "-DACPP_TARGETS=$CudaTargets",
        "-DDEFAULT_TARGETS=cuda:sm_75"
    )
    
    if ($DisableOpenMP) {
        # Disable CPU backend (which includes OpenMP backend) to avoid libomp requirement
        $cmakeArgs += "-DWITH_CPU_BACKEND=OFF"
    } elseif ($LLVMOpenMPDir) {
        $cmakeArgs += "-DLLVM_OpenMP_DIR=$LLVMOpenMPDir"
        $cmakeArgs += "-DOpenMP_ROOT=$LLVMOpenMPDir"
    }
    
    $cmakeArgs += ".."
    
    # Run CMake configuration
    & cmake $cmakeArgs
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "CMake configuration failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "Building AdaptiveCpp (this may take a while)..." -ForegroundColor Yellow
    cmake --build . --config Release
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Build failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "Installing AdaptiveCpp..." -ForegroundColor Yellow
    cmake --install . --config Release
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Installation failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "=== Build Complete ===" -ForegroundColor Green
    Write-Host "AdaptiveCpp installed to: $InstallPrefix" -ForegroundColor Cyan
    Write-Host ""
    Write-Host "Next steps:" -ForegroundColor Yellow
    Write-Host "1. Set environment variable: `$env:CMAKE_PREFIX_PATH = '$InstallPrefix'" -ForegroundColor White
    Write-Host "2. Or use: cmake -DCMAKE_PREFIX_PATH=$InstallPrefix .." -ForegroundColor White
    Write-Host "3. Configure test_convolution with:" -ForegroundColor White
    Write-Host "   cmake -DENABLE_SYCL=ON -DPREFER_ADAPTIVECPP=ON -DCMAKE_PREFIX_PATH=$InstallPrefix .." -ForegroundColor White
    
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
    exit 1
} finally {
    Pop-Location
}
