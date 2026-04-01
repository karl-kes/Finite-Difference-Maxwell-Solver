if (-not (Test-Path build_debug)) {
    New-Item -ItemType Directory build_debug | Out-Null
}
Push-Location build_debug
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Debug
mingw32-make -j
Pop-Location