
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
