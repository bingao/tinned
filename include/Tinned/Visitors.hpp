/* Tinned: a set of nonnumerical routines for computational chemistry
   Copyright 2023-2024 Bin Gao

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   This file is the header file of different visitors.

   2025-02-23, Bin Gao:
   * first version
*/

#pragma once

#include "Tinned/visitors/ZerosRemover.hpp"
#include "Tinned/visitors/RemoveVisitor.hpp"
#include "Tinned/visitors/KeepVisitor.hpp"
#include "Tinned/visitors/ReplaceVisitor.hpp"
#include "Tinned/visitors/FindAllVisitor.hpp"
#include "Tinned/visitors/EliminationVisitor.hpp"
#include "Tinned/visitors/TemporumCleaner.hpp"
#include "Tinned/visitors/ExistAnyVisitor.hpp"
#include "Tinned/visitors/LaTeXifyVisitor.hpp"
#include "Tinned/visitors/StringifyVisitor.hpp"
