add_library(tinned
            ${LIB_TINNED_PATH}/src/Perturbation.cpp
            ${LIB_TINNED_PATH}/src/NonElecFunction.cpp
            ${LIB_TINNED_PATH}/src/ElectronicState.cpp
            ${LIB_TINNED_PATH}/src/OneElecDensity.cpp
            ${LIB_TINNED_PATH}/src/OneElecOperator.cpp
            ${LIB_TINNED_PATH}/src/TwoElecEnergy.cpp
            ${LIB_TINNED_PATH}/src/TwoElecOperator.cpp
            ${LIB_TINNED_PATH}/src/CompositeFunction.cpp
            ${LIB_TINNED_PATH}/src/ExchCorrContraction.cpp
            ${LIB_TINNED_PATH}/src/ExchCorrEnergy.cpp
            ${LIB_TINNED_PATH}/src/ExchCorrPotential.cpp
            ${LIB_TINNED_PATH}/src/TemporumOperator.cpp
            ${LIB_TINNED_PATH}/src/TemporumOverlap.cpp
            ${LIB_TINNED_PATH}/src/KeepVisitor.cpp
            ${LIB_TINNED_PATH}/src/RemoveVisitor.cpp
            ${LIB_TINNED_PATH}/src/ReplaceVisitor.cpp
            ${LIB_TINNED_PATH}/src/FindAllVisitor.cpp
            ${LIB_TINNED_PATH}/src/StringifyVisitor.cpp)

target_include_directories(tinned INTERFACE
    $<BUILD_INTERFACE:${TINNED_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
