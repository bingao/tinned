# Tinned

Tinned contains a set of nonnumerical routines for computational chemistry, in
which symbolic routines are built on top of
[SymEngine library](https://github.com/symengine/symengine).

## License

Tinned is licensed under Mozilla Public License Version 2.0, see the
[LICENSE](LICENSE) for more information.

## Installation

Both SymEngine and Tinned require CMake and C++ compiler which supports C++11
standard.

Clone and build [forked SymEngine library](https://github.com/bingao/symengine)
first, which has implmented derivatives for different matrix expressions.

Then clone Tinned library and build it by setting `SymEngine_DIR` to the
SymEngine installation or build directory.

## Try Tinned

Tinned currently provides C++ interface and the following classes and functions
for computational chemistry:

* `class Perturbation`,
* `class ElectronicState`
