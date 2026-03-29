Remove-Item -Recurse -Force build_release -ErrorAction SilentlyContinue
New-Item -ItemType Directory build_release | Out-Null
Push-Location build_release
cmake .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
mingw32-make -j
Pop-Location
