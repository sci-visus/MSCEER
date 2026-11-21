# PowerShell script to diagnose AdaptiveCpp setup for test_convolution
# Usage: .\check_adaptivecpp_setup.ps1

Write-Host "=== AdaptiveCpp Setup Diagnostic ===" -ForegroundColor Cyan
Write-Host ""

# Check CMake version
Write-Host "Checking CMake version..." -ForegroundColor Yellow
$cmakePath = Get-Command cmake -ErrorAction SilentlyContinue
if ($cmakePath) {
    $cmakeVersion = & cmake --version | Select-Object -First 1
    Write-Host "  Found: $($cmakePath.Source)" -ForegroundColor Green
    Write-Host "  Version: $cmakeVersion" -ForegroundColor Gray
    
    # Check if version is >= 3.22 (required for AdaptiveCpp)
    $versionMatch = $cmakeVersion -match "version (\d+)\.(\d+)"
    if ($versionMatch) {
        $major = [int]$matches[1]
        $minor = [int]$matches[2]
        if ($major -gt 3 -or ($major -eq 3 -and $minor -ge 22)) {
            Write-Host "  CMake version is sufficient for AdaptiveCpp (>= 3.22)" -ForegroundColor Green
        } else {
            Write-Host "  CMake version is too old for AdaptiveCpp (requires >= 3.22)" -ForegroundColor Red
        }
    }
} else {
    Write-Host "  CMake not found in PATH" -ForegroundColor Red
}

Write-Host ""

# Check for AdaptiveCpp installation
Write-Host "Checking for AdaptiveCpp installation..." -ForegroundColor Yellow
$acppPaths = @(
    $env:ACPP_ROOT,
    $env:CMAKE_PREFIX_PATH,
    "C:\Users\jediati\bin\AdaptiveCpp-LLVM20-Win",
    "$env:ProgramFiles\AdaptiveCpp",
    "$env:ProgramFiles(x86)\AdaptiveCpp",
    "C:\AdaptiveCpp",
    "$env:USERPROFILE\AdaptiveCpp",
    "$env:USERPROFILE\.local\AdaptiveCpp"
)

$acppFound = $false
foreach ($path in $acppPaths) {
    if ($path -and (Test-Path $path)) {
        # Check for AdaptiveCppConfig.cmake (try both naming conventions)
        $configFiles = @(
            "$path\lib\cmake\AdaptiveCpp\AdaptiveCppConfig.cmake",
            "$path\lib\cmake\AdaptiveCpp\adaptivecpp-config.cmake",
            "$path\share\AdaptiveCpp\AdaptiveCppConfig.cmake",
            "$path\AdaptiveCppConfig.cmake"
        )
        
        foreach ($configFile in $configFiles) {
            if (Test-Path $configFile) {
                Write-Host "  Found AdaptiveCpp at: $path" -ForegroundColor Green
                Write-Host "    Config file: $configFile" -ForegroundColor Gray
                $acppFound = $true
                break
            }
        }
        
        if ($acppFound) {
            break
        }
    }
}

if (-not $acppFound) {
    Write-Host "  AdaptiveCpp not found in standard locations" -ForegroundColor Yellow
    Write-Host "    Checked paths:" -ForegroundColor Gray
    foreach ($path in $acppPaths) {
        if ($path) {
            Write-Host "      - $path" -ForegroundColor Gray
        }
    }
}

Write-Host ""

# Check for clang++ compiler (recommended for AdaptiveCpp)
Write-Host "Checking for clang++ compiler..." -ForegroundColor Yellow
$clangPath = Get-Command clang++ -ErrorAction SilentlyContinue
if ($clangPath) {
    Write-Host "  Found: $($clangPath.Source)" -ForegroundColor Green
    
    try {
        $version = & clang++ --version 2>&1 | Select-Object -First 1
        Write-Host "  Version: $version" -ForegroundColor Gray
        
        # Check LLVM version (AdaptiveCpp requires LLVM 14+)
        $llvmMatch = $version -match "clang version (\d+)"
        if ($llvmMatch) {
            $llvmVersion = [int]$matches[1]
            if ($llvmVersion -ge 14) {
                Write-Host "  LLVM version is sufficient for AdaptiveCpp (>= 14)" -ForegroundColor Green
            } else {
                Write-Host "  LLVM version may be too old for AdaptiveCpp (recommends >= 14)" -ForegroundColor Yellow
            }
        }
    } catch {
        Write-Host "  Could not get version info" -ForegroundColor Yellow
    }
} else {
    Write-Host "  clang++ not found in PATH" -ForegroundColor Yellow
    Write-Host "    Note: AdaptiveCpp works best with clang++" -ForegroundColor Gray
}

Write-Host ""

# Check environment variables
Write-Host "Checking environment variables..." -ForegroundColor Yellow
if ($env:ACPP_ROOT) {
    Write-Host "  ACPP_ROOT: $env:ACPP_ROOT" -ForegroundColor Green
} else {
    Write-Host "  ACPP_ROOT: Not set" -ForegroundColor Yellow
}

if ($env:CMAKE_PREFIX_PATH) {
    Write-Host "  CMAKE_PREFIX_PATH: $env:CMAKE_PREFIX_PATH" -ForegroundColor Green
} else {
    Write-Host "  CMAKE_PREFIX_PATH: Not set" -ForegroundColor Yellow
}

Write-Host ""

# Test CMake configuration
Write-Host "Testing CMake AdaptiveCpp configuration..." -ForegroundColor Yellow
$buildDir = Join-Path $PSScriptRoot "build_test_acpp"
if (Test-Path $buildDir) {
    Remove-Item -Recurse -Force $buildDir
}
New-Item -ItemType Directory -Path $buildDir | Out-Null

try {
    Push-Location $buildDir
    
    # Try to configure with AdaptiveCpp
    $cmakeArgs = @(
        "-DENABLE_SYCL=ON",
        "-DENABLE_OPENMP=OFF",
        "-DPREFER_ADAPTIVECPP=ON"
    )
    
    # Add CMAKE_PREFIX_PATH if ACPP_ROOT is set
    if ($env:ACPP_ROOT) {
        $cmakeArgs += "-DCMAKE_PREFIX_PATH=$env:ACPP_ROOT"
    }
    
    $cmakeOutput = & cmake $cmakeArgs .. 2>&1 | Out-String
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  CMake configuration: SUCCESS" -ForegroundColor Green
        
        # Check for AdaptiveCpp-related messages
        if ($cmakeOutput -match "AdaptiveCpp") {
            Write-Host "  AdaptiveCpp detected in configuration" -ForegroundColor Green
        }
        if ($cmakeOutput -match "SYCL enabled via AdaptiveCpp") {
            Write-Host "  SYCL enabled via AdaptiveCpp" -ForegroundColor Green
        } elseif ($cmakeOutput -match "SYCL.*not found") {
            Write-Host "  SYCL not found by CMake" -ForegroundColor Red
        }
    } else {
        Write-Host "  CMake configuration: FAILED" -ForegroundColor Red
        Write-Host "  Error output:" -ForegroundColor Red
        $cmakeOutput -split "`n" | Select-Object -Last 10 | ForEach-Object { Write-Host "    $_" -ForegroundColor Red }
    }
} catch {
    Write-Host "  CMake test failed: $_" -ForegroundColor Red
} finally {
    Pop-Location
    if (Test-Path $buildDir) {
        Remove-Item -Recurse -Force $buildDir
    }
}

Write-Host ""
Write-Host "=== Recommendations ===" -ForegroundColor Cyan
if (-not $acppFound) {
    Write-Host "1. Install AdaptiveCpp:" -ForegroundColor Yellow
    Write-Host "   - Download from: https://github.com/AdaptiveCpp/AdaptiveCpp" -ForegroundColor White
    Write-Host "   - Or build from source following the AdaptiveCpp documentation" -ForegroundColor White
    Write-Host ""
    Write-Host "2. After installation, set one of:" -ForegroundColor Yellow
    Write-Host "   - ACPP_ROOT environment variable pointing to AdaptiveCpp install directory" -ForegroundColor White
    Write-Host "   - CMAKE_PREFIX_PATH environment variable" -ForegroundColor White
    Write-Host "   - Use: cmake -DCMAKE_PREFIX_PATH=<path-to-acpp> .." -ForegroundColor White
} else {
    Write-Host "Your system appears to have AdaptiveCpp support. Try building with:" -ForegroundColor Green
    Write-Host "  cmake -DENABLE_SYCL=ON -DPREFER_ADAPTIVECPP=ON -DCMAKE_PREFIX_PATH=<acpp-path> .." -ForegroundColor White
    Write-Host "  cmake --build ." -ForegroundColor White
}

Write-Host ""
