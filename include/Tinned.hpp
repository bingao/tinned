/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of Tinned library.

   2024-05-08, Bin Gao:
   * add more visitors for response theory

   2024-05-04, Bin Gao:
   * add latexify visitor

   2024-04-28, Bin Gao:
   * add elimination visitor

   2024-04-17, Bin Gao:
   * support coupled-cluster response theory

   2023-10-22, Bin Gao:
   * first version
*/

#pragma once

#include "Tinned/Perturbation.hpp"
#include "Tinned/PertDependency.hpp"
#include "Tinned/PertTuple.hpp"
#include "Tinned/PerturbedParameter.hpp"
#include "Tinned/ZeroOperator.hpp"
#include "Tinned/ConjugateTranspose.hpp"

#include "Tinned/OneElecDensity.hpp"
#include "Tinned/OneElecOperator.hpp"
#include "Tinned/TwoElecEnergy.hpp"
#include "Tinned/TwoElecOperator.hpp"
#include "Tinned/ExchCorrContraction.hpp"
#include "Tinned/ExchCorrEnergy.hpp"
#include "Tinned/ExchCorrPotential.hpp"
#include "Tinned/NonElecFunction.hpp"
#include "Tinned/TemporumOperator.hpp"
#include "Tinned/TemporumOverlap.hpp"

#include "Tinned/AdjointMap.hpp"
#include "Tinned/ClusterConjHamiltonian.hpp"

#include "Tinned/ZerosRemover.hpp"
#include "Tinned/RemoveVisitor.hpp"
#include "Tinned/KeepVisitor.hpp"
#include "Tinned/ReplaceVisitor.hpp"
#include "Tinned/FindAllVisitor.hpp"
#include "Tinned/EliminationVisitor.hpp"
#include "Tinned/TemporumCleaner.hpp"
#include "Tinned/ExistAnyVisitor.hpp"
#include "Tinned/LaTeXifyVisitor.hpp"
#include "Tinned/StringifyVisitor.hpp"
