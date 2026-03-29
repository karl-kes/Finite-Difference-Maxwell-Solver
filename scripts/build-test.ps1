Remove-Item -Recurse -Force build_test -ErrorAction SilentlyContinue
New-Item -ItemType Directory build_test | Out-Null
Push-Location build_test
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Test
mingw32-make -j run_tests
ctest --output-on-failure
Pop-Location
