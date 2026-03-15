$ErrorActionPreference = "Stop"

$Mode = if ( $args.Count -ge 1 ) { $args[0] } else { "release" }
$BuildDir = "build"
$Jobs = [Environment]::ProcessorCount

switch ( $Mode.ToLower() ) {
    { $_ -in "release", "r" } { $BuildType = "Release" }
    { $_ -in "debug",   "d" } { $BuildType = "Debug"   }
    "clean" {
        Write-Host "Cleaning build directory..."
        if ( Test-Path $BuildDir ) { Remove-Item -Recurse -Force $BuildDir }
        Write-Host "Done."
        exit 0
    }
    default {
        Write-Host "Usage: .\build.ps1 [release|debug|clean]"
        exit 1
    }
}

Write-Host "=== $BuildType Build ==="

# Configure (only if needed or build type changed):
$NeedsConfigure = $true
$CacheFile = "$BuildDir\CMakeCache.txt"

if ( Test-Path $CacheFile ) {
    $CacheContent = Get-Content $CacheFile -Raw
    if ( $CacheContent -match "CMAKE_BUILD_TYPE:STRING=$BuildType" ) {
        $NeedsConfigure = $false
    }
}

if ( $NeedsConfigure ) {
    Write-Host "Configuring ($BuildType)..."
    if ( Test-Path $BuildDir ) { Remove-Item -Recurse -Force $BuildDir }
    cmake -B $BuildDir -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=$BuildType
    if ( $LASTEXITCODE -ne 0 ) { exit $LASTEXITCODE }
}

# Build:
Write-Host "Building with $Jobs threads..."
cmake --build $BuildDir -j $Jobs
if ( $LASTEXITCODE -ne 0 ) { exit $LASTEXITCODE }

Write-Host ""
Write-Host "=== Build Complete ==="
Write-Host "Binary: $BuildDir\main.exe"

# Run:
Write-Host ""
& ".\$BuildDir\main.exe"