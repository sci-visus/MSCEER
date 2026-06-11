# PowerShell script to build test_convolution with AdaptiveCpp
# Usage: .\build_with_adaptivecpp.ps1 [options]

param(
    [string]$AdaptiveCppPath = "C:\Users\jediati\bin\AdaptiveCpp-LLVM20-Win",
    [string]$Compiler = "auto",  # auto, msvc, icpx, clang++
    [string]$Generator = "auto",  # auto, Ninja, Visual Studio 18 2026
    [switch]$Clean = $false
)

$ErrorActionPreference = "Stop"

Write-Host "=== Building test_convolution with AdaptiveCpp ===" -ForegroundColor Cyan
Write-Host ""

# Check AdaptiveCpp installation
if (-not (Test-Path $AdaptiveCppPath)) {
    Write-Host "Error: AdaptiveCpp not found at: $AdaptiveCppPath" -ForegroundColor Red
    exit 1
}

$configFile = "$AdaptiveCppPath\lib\cmake\AdaptiveCpp\adaptivecpp-config.cmake"
if (-not (Test-Path $configFile)) {
    $configFile = "$AdaptiveCppPath\lib\cmake\AdaptiveCpp\AdaptiveCppConfig.cmake"
    if (-not (Test-Path $configFile)) {
        Write-Host "Error: AdaptiveCpp CMake config not found" -ForegroundColor Red
        exit 1
    }
}

Write-Host "AdaptiveCpp found: $AdaptiveCppPath" -ForegroundColor Green
Write-Host ""

# Detect compiler
$useCompiler = $null
$useGenerator = $null

if ($Compiler -eq "auto" -or $Generator -eq "auto") {
    # Check for Intel oneAPI DPC++
    $icpxPaths = @(
        "C:\Program Files (x86)\Intel\oneAPI\2025.3\bin\icpx.exe",
        "C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin\icpx.exe",
        "C:\Program Files\Intel\oneAPI\compiler\latest\windows\bin\icpx.exe"
    )
    
    # Also check version-specific directories
    $oneApiBase = "C:\Program Files (x86)\Intel\oneAPI"
    if (Test-Path $oneApiBase) {
        $versions = Get-ChildItem -Path $oneApiBase -Directory -Filter "20*" -ErrorAction SilentlyContinue
        foreach ($version in $versions) {
            $icpxPath = Join-Path $version.FullName "bin\icpx.exe"
            if (Test-Path $icpxPath) {
                $icpxPaths = @($icpxPath) + $icpxPaths
            }
        }
    }
    
    $icpxFound = $false
    foreach ($path in $icpxPaths) {
        if (Test-Path $path) {
            $useCompiler = $path
            $icpxFound = $true
            Write-Host "Found Intel DPC++ compiler: $path" -ForegroundColor Green
            break
        }
    }
    
    # Check for Ninja (check multiple locations)
    $ninjaFound = $false
    $ninjaPaths = @(
        "C:\Users\jediati\bin\ninja.exe",
        "C:\Users\jediari\bin\ninja.exe",  # Check for typo variant
        "$env:USERPROFILE\bin\ninja.exe"
    )
    
    foreach ($ninjaPath in $ninjaPaths) {
        if (Test-Path $ninjaPath) {
            $env:Path = "$(Split-Path $ninjaPath -Parent);$env:Path"
            $ninjaFound = $true
            Write-Host "Found Ninja at: $ninjaPath" -ForegroundColor Green
            break
        }
    }
    
    # Also check if ninja is in PATH
    if (-not $ninjaFound) {
        $ninja = Get-Command ninja -ErrorAction SilentlyContinue
        if ($ninja) {
            $ninjaFound = $true
            Write-Host "Found Ninja in PATH: $($ninja.Source)" -ForegroundColor Green
        }
    }
    
    if ($ninjaFound -and $icpxFound) {
        $useGenerator = "Ninja"
        Write-Host "Using Ninja generator with Intel DPC++" -ForegroundColor Green
    } elseif ($icpxFound) {
        $useGenerator = "Visual Studio 18 2026"
        Write-Host "Using Visual Studio generator with Intel DPC++" -ForegroundColor Yellow
    } else {
        $useGenerator = "Visual Studio 18 2026"
        Write-Host "Using Visual Studio generator with MSVC" -ForegroundColor Yellow
    }
} else {
    $useGenerator = $Generator
    if ($Compiler -ne "auto") {
        $useCompiler = $Compiler
    }
}

Write-Host ""

# Setup build directory
$buildDir = Join-Path $PSScriptRoot "build_acpp"
if ($Clean -and (Test-Path $buildDir)) {
    Write-Host "Cleaning build directory..." -ForegroundColor Yellow
    Remove-Item $buildDir -Recurse -Force
}

if (-not (Test-Path $buildDir)) {
    New-Item -ItemType Directory -Path $buildDir | Out-Null
}

Push-Location $buildDir

try {
    Write-Host "Configuring CMake..." -ForegroundColor Yellow
    Write-Host "  Generator: $useGenerator" -ForegroundColor Gray
    if ($useCompiler) {
        Write-Host "  Compiler: $useCompiler" -ForegroundColor Gray
    }
    Write-Host "  AdaptiveCpp: $AdaptiveCppPath" -ForegroundColor Gray
    Write-Host ""
    
    # Build CMake arguments
    $cmakeArgs = @(
        "-G", $useGenerator
    )
    
    if ($useGenerator -eq "Visual Studio 18 2026") {
        $cmakeArgs += "-A", "x64"
    }
    
    # Set Ninja path if using Ninja generator
    if ($useGenerator -eq "Ninja") {
        $ninjaExe = "C:\Users\jediati\bin\ninja.exe"
        if (Test-Path $ninjaExe) {
            $cmakeArgs += "-DCMAKE_MAKE_PROGRAM=$ninjaExe"
        }
    }
    
    if ($useCompiler) {
        $cmakeArgs += "-DCMAKE_CXX_COMPILER=$useCompiler"
    }
    
    $cmakeArgs += @(
        "-DENABLE_SYCL=ON",
        "-DPREFER_ADAPTIVECPP=ON",
        "-DCMAKE_PREFIX_PATH=$AdaptiveCppPath",
        "-DCMAKE_BUILD_TYPE=Release",
        ".."
    )
    
    # Run CMake configuration
    & cmake $cmakeArgs
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "CMake configuration failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "Building project..." -ForegroundColor Yellow
    
    if ($useGenerator -eq "Ninja") {
        cmake --build . --config Release
    } else {
        cmake --build . --config Release
    }
    
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Build failed!" -ForegroundColor Red
        exit 1
    }
    
    Write-Host ""
    Write-Host "=== Build Complete ===" -ForegroundColor Green
    $exePath = Join-Path $buildDir "test_convolution.exe"
    if (Test-Path $exePath) {
        Write-Host "Executable: $exePath" -ForegroundColor Cyan
    } else {
        $exePath = Join-Path $buildDir "Release\test_convolution.exe"
        if (Test-Path $exePath) {
            Write-Host "Executable: $exePath" -ForegroundColor Cyan
        }
    }
    
} catch {
    Write-Host "Error: $_" -ForegroundColor Red
    exit 1
} finally {
    Pop-Location
}
