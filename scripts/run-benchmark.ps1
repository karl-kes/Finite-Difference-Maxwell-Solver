if (-not (Test-Path build_release)) {
    New-Item -ItemType Directory build_release | Out-Null
}
Push-Location build_release
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j benchmark
Pop-Location

Write-Host ""
Write-Host "Running benchmark..."
Write-Host ""

# Default: PML on, 100 steps, max grid 200^3
& .\build_release\benchmark.exe @args
