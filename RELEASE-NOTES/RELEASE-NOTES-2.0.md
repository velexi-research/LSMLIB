Release Notes: LSMLIB v2.0
==========================

--------------------------------------------------------------------------------------------
v2.0.0 (2022-05-08)
===================

LSMLIB v2.0.0 is intended to be a stable version of LSMLIB that can serve as a foundation
for specialized applications and packages. The library has been refactored and only
includes support for serial computations implemented in C/C++ (the parallel LSMLIB and
Python interfaces have been moved to their own separate codebases and repositories). There
are currently no plans to further develop LSMLIB parallel LSMLIB, or the Python interface.

## Enhancements

* Refactor codebase to support only serial computations implemented in C/C++. Support for
  parallel computations and Python implementations have been moved to their own separate
  repositories.
* Migrate to CMake (from Autotools) for build system.
* Added basic unit tests and CI (continuous integration).
* Update license.

## Bug Fixes

* Fix various bugs (e.g., fast marching method declaration and implementation, grid
  management example).
* Fix include warnings identified by include-what-you-use (IWYU).

--------------------------------------------------------------------------------------------
