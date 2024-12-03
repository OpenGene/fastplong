git submodule init
./vcpkg/bootstrap-vcpkg.sh
make -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake -S . -B builds/ninja-multi-vcpkg -G "Ninja Multi-Config"
