# PowerShell script to diagnose SYCL setup for test_convolution
# Usage: .\check_sycl_setup.ps1

Write-Host "=== SYCL Setup Diagnostic ===" -ForegroundColor Cyan
Write-Host ""

# Check for Intel oneAPI installation
Write-Host "Checking for Intel oneAPI installation..." -ForegroundColor Yellow
$oneApiPaths = @(
    "$env:ONEAPI_ROOT",
    "C:\Program Files (x86)\Intel\oneAPI",
    "C:\Program Files\Intel\oneAPI",
    "$env:ProgramFiles(x86)\Intel\oneAPI",
    "$env:ProgramFiles\Intel\oneAPI"
)

$oneApiFound = $false
foreach ($path in $oneApiPaths) {
    if ($path -and (Test-Path $path)) {
        Write-Host "  Found oneAPI at: $path" -ForegroundColor Green
        $oneApiFound = $true
        
        # Check for compiler
        $compilerPaths = @(
            "$path\compiler\latest\windows\bin\icpx.exe",
            "$path\compiler\latest\windows\bin\dpcpp.exe",
            "$path\compiler\latest\linux\bin\icpx",
            "$path\compiler\latest\linux\bin\dpcpp"
        )
        
        foreach ($compilerPath in $compilerPaths) {
            if (Test-Path $compilerPath) {
                Write-Host "    Found compiler: $compilerPath" -ForegroundColor Green
                
                # Check compiler version
                try {
                    $version = & $compilerPath --version 2>&1
                    Write-Host "    Version info:" -ForegroundColor Gray
                    $version | Select-Object -First 3 | ForEach-Object { Write-Host "      $_" -ForegroundColor Gray }
                } catch {
                    Write-Host "    Could not get version info" -ForegroundColor Yellow
                }
                break
            }
        }
        break
    }
}

if (-not $oneApiFound) {
    Write-Host "  Intel oneAPI not found in standard locations" -ForegroundColor Red
}

Write-Host ""

# Check for SYCL compilers in PATH
Write-Host "Checking for SYCL compilers in PATH..." -ForegroundColor Yellow
$compilers = @("icpx", "dpcpp", "clang++")
$foundCompiler = $false
foreach ($compiler in $compilers) {
    $compilerPath = Get-Command $compiler -ErrorAction SilentlyContinue
    if ($compilerPath) {
        Write-Host "  Found: $($compilerPath.Source)" -ForegroundColor Green
        
        # Check if it supports SYCL
        try {
            $version = & $compiler --version 2>&1 | Out-String
            if ($version -match "DPC|oneAPI|SYCL") {
                Write-Host "    Supports SYCL: Yes" -ForegroundColor Green
                $foundCompiler = $true
            } else {
                Write-Host "    Supports SYCL: Unknown" -ForegroundColor Yellow
            }
        } catch {
            Write-Host "    Could not check SYCL support" -ForegroundColor Yellow
        }
    }
}

if (-not $foundCompiler) {
    Write-Host "  No SYCL-capable compilers found in PATH" -ForegroundColor Red
}

Write-Host ""

# Check environment variables
Write-Host "Checking environment variables..." -ForegroundColor Yellow
if ($env:ONEAPI_ROOT) {
    Write-Host "  ONEAPI_ROOT: $env:ONEAPI_ROOT" -ForegroundColor Green
} else {
    Write-Host "  ONEAPI_ROOT: Not set" -ForegroundColor Yellow
}

Write-Host ""

# Check for CMake
Write-Host "Checking CMake..." -ForegroundColor Yellow
$cmakePath = Get-Command cmake -ErrorAction SilentlyContinue
if ($cmakePath) {
    Write-Host "  Found: $($cmakePath.Source)" -ForegroundColor Green
    try {
        $cmakeVersion = & cmake --version | Select-Object -First 1
        Write-Host "  Version: $cmakeVersion" -ForegroundColor Gray
    } catch {
        Write-Host "  Could not get version" -ForegroundColor Yellow
    }
} else {
    Write-Host "  CMake not found in PATH" -ForegroundColor Red
}

Write-Host ""

# Test CMake configuration
Write-Host "Testing CMake SYCL configuration..." -ForegroundColor Yellow
$buildDir = Join-Path $PSScriptRoot "build_test"
if (Test-Path $buildDir) {
    Remove-Item -Recurse -Force $buildDir
}
New-Item -ItemType Directory -Path $buildDir | Out-Null

try {
    Push-Location $buildDir
    $cmakeOutput = & cmake -DENABLE_SYCL=ON -DENABLE_OPENMP=OFF .. 2>&1 | Out-String
    
    if ($LASTEXITCODE -eq 0) {
        Write-Host "  CMake configuration: SUCCESS" -ForegroundColor Green
        
        # Check for SYCL-related messages
        if ($cmakeOutput -match "SYCL enabled") {
            Write-Host "  SYCL enabled in configuration" -ForegroundColor Green
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
if (-not $oneApiFound -and -not $foundCompiler) {
    Write-Host "1. Install Intel oneAPI Base Toolkit from:" -ForegroundColor Yellow
    Write-Host "   https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html" -ForegroundColor White
    Write-Host ""
    Write-Host "2. After installation, either:" -ForegroundColor Yellow
    Write-Host "   a) Set ONEAPI_ROOT environment variable" -ForegroundColor White
    Write-Host "   b) Add icpx/dpcpp to PATH" -ForegroundColor White
    Write-Host "   c) Use: cmake -DCMAKE_CXX_COMPILER=icpx .." -ForegroundColor White
} else {
    Write-Host "Your system appears to have SYCL support. Try building with:" -ForegroundColor Green
    Write-Host "  cmake -DENABLE_SYCL=ON -DENABLE_OPENMP=OFF .." -ForegroundColor White
    Write-Host "  cmake --build ." -ForegroundColor White
}

Write-Host ""
