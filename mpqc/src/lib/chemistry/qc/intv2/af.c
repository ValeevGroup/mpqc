/*
 * af.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#include <stdio.h>
#include <stdarg.h>
#include <strings.h>
#include <util/sgen/sgen.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv2/atomsfree.h>
#include <chemistry/qc/intv2/atomsinit.h>
#define NAMES_LENGTH 256

/* Generated free function for shell_type: */
void
free_shell_type(_shell_type)
shell_type_t *_shell_type;
{
  }

/* Generated free function for shell: */
void
free_shell(_shell)
shell_t *_shell;
{
  int i;
  if (_shell->nprim) {
    free(_shell->exp);
    }
  if (_shell->ncon) {
    int ii;
    for (ii=0; ii<_shell->ncon; ii++) {
      free_shell_type(&(_shell->type[ii]));
      }
    free(_shell->type);
    }
  if (_shell->ncon && _shell->nprim) {
    for (i=0; i<_shell->ncon; i++)  {
      free(_shell->coef[i]);
      }
    free(_shell->coef);
    }
  }

/* Generated free function for basis: */
void
free_basis(_basis)
basis_t *_basis;
{
  if (_basis->n) {
    int ii;
    for (ii=0; ii<_basis->n; ii++) {
      free_shell(&(_basis->shell[ii]));
      }
    free(_basis->shell);
    }
  }

/* Generated free function for center: */
void
free_center(_center)
center_t *_center;
{
  if (3) {
    free(_center->r);
    }
  free_basis(&_center->basis);
  }

/* Generated free function for centers: */
void
free_centers(_centers)
centers_t *_centers;
{
  if (_centers->n) {
    int ii;
    for (ii=0; ii<_centers->n; ii++) {
      free_center(&(_centers->center[ii]));
      }
    free(_centers->center);
    }
  if (_centers->nshell) {
    free(_centers->center_num);
    }
  if (_centers->nshell) {
    free(_centers->shell_num);
    }
  if (_centers->nshell) {
    free(_centers->func_num);
    }
  }
