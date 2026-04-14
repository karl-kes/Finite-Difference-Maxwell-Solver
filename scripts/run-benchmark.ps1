if (-not (Test-Path build_release)) {
    New-Item -ItemType Directory build_release | Out-Null
}
Push-Location build_release
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j benchmark benchmark_serial
Pop-Location

if (-not (Test-Path benchmark\results)) {
    New-Item -ItemType Directory benchmark\results | Out-Null
}
$Stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$HostTag = $env:COMPUTERNAME
$Tag = "${HostTag}_${Stamp}"

# Thread pinning to cut variance. Override by setting these env vars
# in the shell before running.
if (-not $env:OMP_PLACES)    { $env:OMP_PLACES    = "cores" }
if (-not $env:OMP_PROC_BIND) { $env:OMP_PROC_BIND = "close" }

Write-Host ""
Write-Host "=== OpenMP benchmark (--threads = max) ==="
& .\build_release\benchmark.exe --output "benchmark\results\omp_${Tag}.csv" @args

Write-Host ""
Write-Host "=== Serial benchmark (true single-thread, no -fopenmp) ==="
& .\build_release\benchmark_serial.exe --output "benchmark\results\serial_${Tag}.csv" @args

Write-Host ""
Write-Host "Results written to benchmark\results\"
Write-Host "  benchmark\results\omp_${Tag}.csv"
Write-Host "  benchmark\results\serial_${Tag}.csv"