//------------------------------------------------------------------------------
// UMFPACK/Source/umfpack_col_to_triplet: convert CSC sparse to triplet
//------------------------------------------------------------------------------

// UMFPACK, Copyright (c) 2005-2022, Timothy A. Davis, All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0+

//------------------------------------------------------------------------------

/*
    User callable.  Converts a column-oriented input matrix to triplet form by
    constructing the column indices Tj from the column pointers Ap.  The matrix
    may be singular.  See umfpack.h for details.

*/

#include "umf_internal.h"

GLOBAL int UMFPACK_col_to_triplet
(
    Int n_col,
    const Int Ap [ ],
    Int Tj [ ]
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Int nz, j, p, p1, p2, length ;

    /* ---------------------------------------------------------------------- */
    /* construct the column indices */
    /* ---------------------------------------------------------------------- */

    if (!Ap || !Tj)
    {
	return (UMFPACK_ERROR_argument_missing) ;
    }
    if (n_col <= 0)
    {
	return (UMFPACK_ERROR_n_nonpositive) ;
    }
    if (Ap [0] != 0)
    {
	return (UMFPACK_ERROR_invalid_matrix) ;
    }
    nz = Ap [n_col] ;
    if (nz < 0)
    {
	return (UMFPACK_ERROR_invalid_matrix) ;
    }

    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	length = p2 - p1 ;
	if (length < 0 || p2 > nz)
	{
	    return (UMFPACK_ERROR_invalid_matrix) ;
	}
	for (p = p1 ; p < p2 ; p++)
	{
	    Tj [p] = j ;
	}
    }

    return (UMFPACK_OK) ;
}
