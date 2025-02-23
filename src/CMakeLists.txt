add_library(tinned
            ${LIB_TINNED_PATH}/src/perturbations/Perturbation.cpp
            ${LIB_TINNED_PATH}/src/operators/PerturbedParameter.cpp
            ${LIB_TINNED_PATH}/src/operators/ZeroOperator.cpp
            ${LIB_TINNED_PATH}/src/operators/ConjugateTranspose.cpp
            ${LIB_TINNED_PATH}/src/operators/NonElecFunction.cpp
            ${LIB_TINNED_PATH}/src/operators/ElectronicState.cpp
            ${LIB_TINNED_PATH}/src/operators/OneElecDensity.cpp
            ${LIB_TINNED_PATH}/src/operators/OneElecOperator.cpp
            ${LIB_TINNED_PATH}/src/operators/TwoElecEnergy.cpp
            ${LIB_TINNED_PATH}/src/operators/TwoElecOperator.cpp
            ${LIB_TINNED_PATH}/src/operators/CompositeFunction.cpp
            ${LIB_TINNED_PATH}/src/operators/ExchCorrContraction.cpp
            ${LIB_TINNED_PATH}/src/operators/ExchCorrEnergy.cpp
            ${LIB_TINNED_PATH}/src/operators/ExchCorrPotential.cpp
            ${LIB_TINNED_PATH}/src/operators/TemporumOperator.cpp
            ${LIB_TINNED_PATH}/src/operators/TemporumOverlap.cpp
            ${LIB_TINNED_PATH}/src/operators/AdjointMap.cpp
            ${LIB_TINNED_PATH}/src/operators/ClusterConjHamiltonian.cpp
            ${LIB_TINNED_PATH}/src/visitors/ZerosRemover.cpp
            ${LIB_TINNED_PATH}/src/visitors/RemoveVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/KeepVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/ReplaceVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/FindAllVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/EliminationVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/TemporumCleaner.cpp
            ${LIB_TINNED_PATH}/src/visitors/ExistAnyVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/LaTeXifyVisitor.cpp
            ${LIB_TINNED_PATH}/src/visitors/StringifyVisitor.cpp)

target_include_directories(tinned INTERFACE
    $<BUILD_INTERFACE:${TINNED_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

add_library(tinned_models
            ${LIB_TINNED_PATH}/src/models/TwoLevelAtom.cpp)

target_include_directories(tinned_models INTERFACE
    $<BUILD_INTERFACE:${TINNED_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
