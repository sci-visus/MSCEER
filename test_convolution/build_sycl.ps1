# Build script for SYCL-enabled test_convolution
# This script sources oneAPI environment and builds with SYCL support

$ErrorActionPreference = "Stop"

Write-Host "=== Building test_convolution with SYCL Support ===" -ForegroundColor Cyan
Write-Host ""

# Set oneAPI root
$oneApiRoot = "C:\Program Files (x86)\Intel\oneAPI"
$setvarsPath = Join-Path $oneApiRoot "setvars.bat"

if (-not (Test-Path $setvarsPath)) {
    Write-Host "Error: oneAPI setvars.bat not found at: $setvarsPath" -ForegroundColor Red
    exit 1
}

# Find icpx compiler
$icpxPath = Get-ChildItem "$oneApiRoot\compiler" -Recurse -Filter "icpx.exe" -ErrorAction SilentlyContinue | Select-Object -First 1

if (-not $icpxPath) {
    Write-Host "Error: icpx.exe not found in oneAPI installation" -ForegroundColor Red
    exit 1
}

Write-Host "Found Intel DPC++ compiler: $($icpxPath.FullName)" -ForegroundColor Green
Write-Host ""

# Setup build directory
$buildDir = Join-Path $PSScriptRoot "build"
if (Test-Path $buildDir) {
    Write-Host "Cleaning build directory..." -ForegroundColor Yellow
    Remove-Item $buildDir -Recurse -Force
}
New-Item -ItemType Directory -Path $buildDir | Out-Null

Push-Location $buildDir

try {
    Write-Host "Configuring CMake with SYCL support..." -ForegroundColor Yellow
    
    # Use cmd to source setvars and run cmake with Ninja generator
    $cmakeCmd = @"
@echo off
call "$setvarsPath" >nul
set ONEAPI_ROOT=$oneApiRoot
cmake -G "Ninja" -DENABLE_SYCL=ON -DENABLE_OPENMP=OFF -DCMAKE_CXX_COMPILER="$($icpxPath.FullName)" ..
"@
    
    $cmakeCmd | cmd /c
    if ($LASTEXITCODE -ne 0) {
        Write-Host "CMake configuration failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "Building project with Ninja..." -ForegroundColor Yellow
    cmake --build . --config Release
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Build failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "=== Build Complete ===" -ForegroundColor Green
    # With Ninja, the executable is directly in build directory (no Release subdirectory)
    $exePath = Join-Path $buildDir "test_convolution.exe"
    if (Test-Path $exePath) {
        Write-Host "Executable: $exePath" -ForegroundColor Cyan
    } else {
        Write-Host "Executable location unknown (check build directory)" -ForegroundColor Yellow
    }
    
} finally {
    Pop-Location
}
