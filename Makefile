ifeq ($(OS),Windows_NT)
  SHELL := powershell.exe
  .SHELLFLAGS := -NoProfile -NoLogo
  MKDIR := @$$null = new-item -itemtype directory -force
  TOUCH := @$$null = new-item -force
  RM := remove-item -force
  CMAKE := cmake
  define rmdir
    if (Test-Path $1) { remove-item -recurse $1 }
  endef
else
  MKDIR := mkdir -p
  TOUCH := touch
  RM := rm -rf
  CMAKE := $(shell (command -v cmake3 || command -v cmake || echo cmake))
  define rmdir
    rm -rf $1
  endef
endif

CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE)
# Extra CMake flags which extend the default set
CMAKE_EXTRA_FLAGS := -DCMAKE_TOOLCHAIN_FILE=vcpkg/scripts/buildsystems/vcpkg.cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

all: build

build: builds/.ran-cmake
	cmake --build --preset ninja-vcpkg-release

test: build
	./builds/ninja-multi-vcpkg/Release/fastplong_tests

benchmark: build
	./builds/ninja-multi-vcpkg/Release/fastplong_benchmarks

format:
	find src -name "*.cpp" | xargs clang-format -i
	find src -name "*.h" | xargs clang-format -i

lint: builds/.ran-cmake
	python3 scripts/run-clang-tidy.py -p builds/ninja-multi-vcpkg

vcpkg/.vcpkg-root:
	git submodule update --init --recursive

builds/.ran-cmake: vcpkg/.vcpkg-root
	$(CMAKE) -B builds/ninja-multi-vcpkg -G "Ninja Multi-Config" $(CMAKE_FLAGS) $(CMAKE_EXTRA_FLAGS) -S .
	$(TOUCH) $@

.PHONY: lint build test format benchmark