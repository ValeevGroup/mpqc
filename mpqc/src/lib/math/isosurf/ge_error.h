/****************************************************************************
 * ge_error.h
 * Author Joel Welling
 * Copyright 1988, Pittsburgh Supercomputing Center, Carnegie Mellon University
 *
 * Permission use, copy, and modify this software and its documentation
 * without fee for personal use or use within your organization is hereby
 * granted, provided that the above copyright notice is preserved in all
 * copies and that that copyright and this permission notice appear in
 * supporting documentation.  Permission to redistribute this software to
 * other organizations or individuals is not granted;  that must be
 * negotiated with the PSC.  Neither the PSC nor Carnegie Mellon
 * University make any representations about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *****************************************************************************/

/*
        This file provides definitions for the generic error handler.
*/

extern void ger_init();
extern void ger_toggledebug();
extern void ger_debug();
extern void ger_error();
extern void ger_fatal();

