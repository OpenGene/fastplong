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

build: builds/ninja-multi-vcpkg/build.ninja
	cmake --build --preset ninja-vcpkg-release

clean:
	cmake --build --preset ninja-vcpkg-release --target clean

install: build
	cmake --install builds/ninja-multi-vcpkg
	@echo "Installed."

test: build
	./builds/ninja-multi-vcpkg/Release/fastplong_tests

benchmark: build
	./builds/ninja-multi-vcpkg/Release/fastplong_benchmarks

format:
	find src -name "*.cpp" | xargs clang-format -i
	find src -name "*.h" | xargs clang-format -i

lint: builds/.ran-cmake
	run-clang-tidy -p builds/ninja-multi-vcpkg

vcpkg/.vcpkg-root:
	git submodule update --init --recursive

builds/ninja-multi-vcpkg/build.ninja: vcpkg/.vcpkg-root CMakeLists.txt CMakePresets.json
	$(CMAKE) --preset "ninja-multi-vcpkg"

.PHONY: lint build test format benchmark