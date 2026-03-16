$ErrorActionPreference = "Stop"

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectDir = Split-Path -Parent $ScriptDir
Set-Location $ProjectDir

$Jobs = (Get-CimInstance Win32_Processor).NumberOfLogicalProcessors
if (-not $Jobs) { $Jobs = 4 }

foreach ($Mode in @("Debug", "Test", "Release")) {
    $Dir = "build_$($Mode.ToLower())"

    Write-Host ""
    Write-Host "================================"
    Write-Host "  Building: $Mode -> $Dir"
    Write-Host "================================"

    if (Test-Path $Dir) { Remove-Item -Recurse -Force $Dir }

    cmake -B $Dir -G "MinGW Makefiles" "-DCMAKE_BUILD_TYPE=$Mode"
    if ($LASTEXITCODE -ne 0) { throw "cmake configure failed for $Mode" }

    cmake --build $Dir -j $Jobs
    if ($LASTEXITCODE -ne 0) { throw "cmake build failed for $Mode" }

    if ($Mode -ne "Release") {
        Write-Host ""
        Write-Host "--- Running tests ($Mode) ---"
        ctest --test-dir $Dir --output-on-failure
        if ($LASTEXITCODE -ne 0) { throw "Tests failed for $Mode" }
    }
}

Write-Host ""
Write-Host "================================"
Write-Host "  All builds complete."
Write-Host "  build_debug\    -> Debug"
Write-Host "  build_test\     -> Test"
Write-Host "  build_release\  -> Release"
Write-Host "================================"