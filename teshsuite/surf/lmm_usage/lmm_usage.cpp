/* A few tests for the maxmin library                                       */

/* Copyright (c) 2007-2019. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/msg.h"
#include "src/kernel/lmm/maxmin.hpp"
#include "src/surf/surf_interface.hpp"
#include "xbt/log.h"
#include "xbt/module.h"
#include "xbt/sysdep.h"
#include <algorithm>
#include <cmath>

XBT_LOG_NEW_DEFAULT_CATEGORY(surf_test, "Messages specific for surf example");

namespace lmm = simgrid::kernel::lmm;

#define PRINT_VAR(var) XBT_DEBUG(_XBT_STRINGIFY(var) " = %g", (var)->get_value())
#define SHOW_EXPR(expr) XBT_DEBUG(_XBT_STRINGIFY(expr) " = %g", (expr))

/*        ______                 */
/*  ==l1==  L2  ==L3==           */
/*        ------                 */

static void test1()
{
  double a = 1.0;
  double b = 10.0;

  lmm::System* Sys    = lmm::make_new_maxmin_system(false);
  lmm::Constraint* L1 = Sys->constraint_new(nullptr, a);
  lmm::Constraint* L2 = Sys->constraint_new(nullptr, b);
  lmm::Constraint* L3 = Sys->constraint_new(nullptr, a);

  lmm::Variable* R_1_2_3 = Sys->variable_new(nullptr, 1.0, -1.0, 3);
  lmm::Variable* R_1     = Sys->variable_new(nullptr, 1.0, -1.0, 1);
  lmm::Variable* R_2     = Sys->variable_new(nullptr, 1.0, -1.0, 1);
  lmm::Variable* R_3     = Sys->variable_new(nullptr, 1.0, -1.0, 1);

  Sys->update_variable_penalty(R_1_2_3, 1.0);
  Sys->update_variable_penalty(R_1, 1.0);
  Sys->update_variable_penalty(R_2, 1.0);
  Sys->update_variable_penalty(R_3, 1.0);

  Sys->expand(L1, R_1_2_3, 1.0);
  Sys->expand(L2, R_1_2_3, 1.0);
  Sys->expand(L3, R_1_2_3, 1.0);

  Sys->expand(L1, R_1, 1.0);
  Sys->expand(L2, R_2, 1.0);
  Sys->expand(L3, R_3, 1.0);

  Sys->solve();

  PRINT_VAR(R_1_2_3);
  PRINT_VAR(R_1);
  PRINT_VAR(R_2);
  PRINT_VAR(R_3);

  Sys->variable_free(R_1_2_3);
  Sys->variable_free(R_1);
  Sys->variable_free(R_2);
  Sys->variable_free(R_3);
  delete Sys;
}

static void test2()
{
  lmm::System* Sys = lmm::make_new_maxmin_system(false);

  lmm::Constraint* CPU1 = Sys->constraint_new(nullptr, 200.0);
  lmm::Constraint* CPU2 = Sys->constraint_new(nullptr, 100.0);

  lmm::Variable* T1 = Sys->variable_new(nullptr, 1.0, -1.0, 1);
  lmm::Variable* T2 = Sys->variable_new(nullptr, 1.0, -1.0, 1);

  Sys->update_variable_penalty(T1, 1.0);
  Sys->update_variable_penalty(T2, 1.0);

  Sys->expand(CPU1, T1, 1.0);
  Sys->expand(CPU2, T2, 1.0);

  Sys->solve();

  PRINT_VAR(T1);
  PRINT_VAR(T2);

  Sys->variable_free(T1);
  Sys->variable_free(T2);
  delete Sys;
}

static void test3()
{
  int flows = 11;
  int links = 10;

  double** A = new double*[links + 5];
  /* array to add the constraints of fictitious variables */
  double B[15] = { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 1, 1, 1, 1, 1 };

  for (int i = 0; i < links + 5; i++) {
    A[i] = new double[flows + 5];
    for (int j = 0; j < flows + 5; j++) {
      A[i][j] = 0.0;

      if (i >= links || j >= flows) {
        A[i][j] = 0.0;
      }
    }
  }

  /*matrix that store the constraints/topology */
  A[0][1] = A[0][7] =                                1.0;
  A[1][1] = A[1][7] = A[1][8] =                      1.0;
  A[2][1] = A[2][8] =                                1.0;
  A[3][8] =                                          1.0;
  A[4][0] = A[4][3] = A[4][9] =                      1.0;
  A[5][0] = A[5][3] = A[5][4] = A[5][9] =            1.0;
  A[6][0] = A[6][4] = A[6][9] = A[6][10] =           1.0;
  A[7][2] = A[7][4] = A[7][6] = A[7][9] = A[7][10] = 1.0;
  A[8][2] = A[8][10] =                               1.0;
  A[9][5] = A[9][6] = A[9][9] =                      1.0;
  A[10][11] =                                        1.0;
  A[11][12] =                                        1.0;
  A[12][13] =                                        1.0;
  A[13][14] =                                        1.0;
  A[14][15] =                                        1.0;

  lmm::System* Sys = lmm::make_new_maxmin_system(false);

  /* Creates the constraints */
  lmm::Constraint** tmp_cnst = new lmm::Constraint*[15];
  for (int i = 0; i < 15; i++)
    tmp_cnst[i] = Sys->constraint_new(nullptr, B[i]);

  /* Creates the variables */
  lmm::Variable** tmp_var = new lmm::Variable*[16];
  for (int j = 0; j < 16; j++) {
    tmp_var[j] = Sys->variable_new(nullptr, 1.0, -1.0, 15);
    Sys->update_variable_penalty(tmp_var[j], 1.0);
  }

  /* Link constraints and variables */
  for (int i = 0; i < 15; i++)
    for (int j = 0; j < 16; j++)
      if (A[i][j])
        Sys->expand(tmp_cnst[i], tmp_var[j], 1.0);

  Sys->solve();

  for (int j = 0; j < 16; j++)
    PRINT_VAR(tmp_var[j]);

  for (int j = 0; j < 16; j++)
    Sys->variable_free(tmp_var[j]);
  delete[] tmp_var;
  delete[] tmp_cnst;
  delete Sys;
  for (int i = 0; i < links + 5; i++)
    delete[] A[i];
  delete[] A;
}

int main(int argc, char** argv)
{
  MSG_init(&argc, argv);
  XBT_INFO("***** Test 1");
  test1();

  XBT_INFO("***** Test 2");
  test2();

  XBT_INFO("***** Test 3");
  test3();

  return 0;
}
