/**
 * @file matrix.c
 * @author Benjamin Ming Yang @ Department of Geoscience, National Taiwan University (b98204032@gmail.com)
 * @brief
 * @date 2018-03-05
 *
 * @copyright Copyright (c) 2018
 *
 */

/**
 * @name
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#ifdef __USE_AVX_INTRIN
/* */
#include <immintrin.h>
/* */
#define STEP_AVX_PD  4
#endif

/**
 * @name
 *
 */
#include "../include/matrix.h"

/**
 * @name Internal macros
 *
 */
#define ARE_SAME_SIZE( _MATRIX_A, _MATRIX_B ) \
		((_MATRIX_A)->i == (_MATRIX_B)->i && (_MATRIX_A)->j == (_MATRIX_B)->j)

#define IS_SQUARE( _MATRIX ) \
		((_MATRIX)->i == (_MATRIX)->j)

#define IS_EXIST( _MATRIX ) \
		((_MATRIX) ? true : false)

/**
 * @name
 *
 */
static matrix_t *_matrix_new( const size_t, const size_t );
static matrix_t *_matrix_idt( const size_t );
static matrix_t *_matrix_dup( const matrix_t * );
static matrix_t *_matrix_tps( const matrix_t * );
static matrix_t *_matrix_inv( const matrix_t * );

static matrix_t *_matrix_add( const matrix_t *, const matrix_t * );
static matrix_t *_matrix_sub( const matrix_t *, const matrix_t * );
static matrix_t *_matrix_mul( const matrix_t *, const matrix_t * );
static matrix_t *_matrix_wls( const matrix_t *, const matrix_t *, const matrix_t *, const matrix_t * );

static matrix_t *_matrix_scalar_add( matrix_t *, const double );
static matrix_t *_matrix_scalar_mul( matrix_t *, const double );

static matrix_t *_matrix_lup_dec( const matrix_t *, matrix_t **, matrix_t **, matrix_t ** );

static matrix_t *_matrix_set( matrix_t *, const double, const size_t, const size_t );
static matrix_t *_matrix_set_all( matrix_t *, const double *, const size_t );
static matrix_t *_matrix_set_row( matrix_t *, const double *, const size_t, const size_t );
static matrix_t *_matrix_set_col( matrix_t *, const double *, const size_t, const size_t );
static matrix_t *_matrix_set_diag( matrix_t *, const double *, const size_t );
static matrix_t *_matrix_set_va( matrix_t *, const size_t, va_list ap );

static matrix_t *_matrix_get( const matrix_t *, double *, const size_t, const size_t );
static matrix_t *_matrix_get_all( const matrix_t *, double *, const size_t );
static matrix_t *_matrix_get_row( const matrix_t *, double *, const size_t, const size_t );
static matrix_t *_matrix_get_col( const matrix_t *, double *, const size_t, const size_t );
static matrix_t *_matrix_get_diag( const matrix_t *, double *, const size_t );

static matrix_t *_matrix_apply( matrix_t *, double (*)( double, void * ), void *, const size_t, const size_t );
static matrix_t *_matrix_apply_all( matrix_t *, double (*)( double, void * ), void * );
static matrix_t *_matrix_apply_row( matrix_t *, double (*)( double, void * ), void *, const size_t );
static matrix_t *_matrix_apply_col( matrix_t *, double (*)( double, void * ), void *, const size_t );
static matrix_t *_matrix_apply_diag( matrix_t *, double (*)( double, void * ), void * );

static size_t _matrix_rk( const matrix_t * );
static double _matrix_det( const matrix_t * );
static double _matrix_tr( const matrix_t * );

static int  _matrix_cmp( const matrix_t *, const matrix_t *, const double );
static void _matrix_free( matrix_t * );

/**
 * @name
 *
 */
static void *(* MatrixAllocFunc)( size_t ) = malloc;
static void  (* MaxtrixFreeFunc)( void * ) = free;

/**
 * @brief
 *
 * @param alloc_func
 * @param free_func
 * @return int
 */
int matrix_lib_init( void *(*alloc_func)( size_t ), void (*free_func)( void * ) )
{
	matrix_t *test_matrix = NULL;

/* */
	MatrixAllocFunc = alloc_func;
	MaxtrixFreeFunc = free_func;
/* */
	if ( (test_matrix = matrix_new( 1, 1 )) ) {
		matrix_free( test_matrix );
		return 0;
	}
/* Reset to default functions */
	MatrixAllocFunc = malloc;
	MaxtrixFreeFunc = free;

	return -1;
}

/**
 * @brief
 *
 * @param a
 */
void matrix_print( const matrix_t *a, FILE *stream )
{
/* */
	if ( a ) {
		for ( register size_t i = 0; i < a->i; i++ ) {
			fprintf(stream, "| ");
		/* */
			for ( register size_t j = 0; j < a->j; j++ )
				fprintf(stream, "%lf ", a->element[i * a->j + j]);
		/* */
			fprintf(stream, "|\n");
		}
	}
	else {
		fprintf(stream, "| Null Matrix |\n");
	}

	return;
}

/**
 * @name
 *
 */

/**
 * @brief
 *
 * @param row
 * @param column
 * @return matrix_t*
 */
matrix_t *matrix_new( const size_t row, const size_t column )
{
	return ((row * column) > 1) ? _matrix_new( row, column ) : NULL;
}

/**
 * @brief Given a rank, return (create) the rank * rank identity matrix.
 *
 * @param rank
 * @return matrix_t*
 */
matrix_t *matrix_idt( const size_t rank )
{
	return rank > 1 ? _matrix_idt( rank ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @return matrix_t*
 */
matrix_t *matrix_dup( const matrix_t *a )
{
	return IS_EXIST( a ) ? _matrix_dup( a ) : NULL;
}

/**
 * @brief Given a matrix a, return (create) the transpose of matrix a.
 *
 * @param a
 * @return matrix_t*
 */
matrix_t *matrix_tps( const matrix_t *a )
{
	return IS_EXIST( a ) ? _matrix_tps( a ) : NULL;
}

/**
 * @brief Given a matrix a, return (create) the inverse of matrix a.
 *
 * @param a
 * @return matrix_t*
 */
matrix_t *matrix_inv( const matrix_t *a )
{
	return (IS_EXIST( a ) && IS_SQUARE( a )) ? _matrix_inv( a ) : NULL;
}

/**
 * @name
 *
 */

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
matrix_t *matrix_add( const matrix_t *a, const matrix_t *b )
{
	return (IS_EXIST( a ) && IS_EXIST( b ) && ARE_SAME_SIZE( a, b )) ?
		_matrix_add( a, b ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
matrix_t *matrix_sub( const matrix_t *a, const matrix_t *b )
{
	return (IS_EXIST( a ) && IS_EXIST( b ) && ARE_SAME_SIZE( a, b )) ?
		_matrix_sub( a, b ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
matrix_t *matrix_mul( const matrix_t *a, const matrix_t *b )
{
	return (IS_EXIST( a ) && IS_EXIST( b ) && (a->j == b->i)) ?
		_matrix_mul( a, b ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
matrix_t *matrix_div( const matrix_t *a, const matrix_t *b )
{
	return matrix_wls( b, a, NULL, NULL );
}

/**
 * @brief Return (create) the result Weighted Least Squares(WLS) matrix from 'W * A * x = W * B',
 *        where x should be '(A^T * W * A + Opt)^-1 * A^T * W * B'.
 *
 * @param a
 * @param b
 * @param w
 * @param opt
 * @return matrix_t*
 */
matrix_t *matrix_wls( const matrix_t *a, const matrix_t *b, const matrix_t *w, const matrix_t *opt )
{
	return (IS_EXIST( a ) && IS_EXIST( b ) && (a->i == b->i)) ?
		_matrix_wls( a, b, w, opt ) : NULL;
}

/**
 * @name scalar's functions
 *
 */

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
matrix_t *matrix_scalar_add( matrix_t *a, const double scalar )
{
	return IS_EXIST( a ) ? _matrix_scalar_add( a, scalar ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
matrix_t *matrix_scalar_sub( matrix_t *a, const double scalar )
{
	return matrix_scalar_add( a, -scalar );
}

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
matrix_t *matrix_scalar_mul( matrix_t *a, const double scalar )
{
	return IS_EXIST( a ) ? _matrix_scalar_mul( a, scalar ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
matrix_t *matrix_scalar_div( matrix_t *a, const double scalar )
{
	return matrix_scalar_mul( a, 1.0 / scalar );
}

/**
 * @name Decompose functions
 *
 */

/**
 * @brief
 *
 * @param a
 * @param l
 * @param u
 * @param p
 * @return matrix_t*
 */
matrix_t *matrix_lup_dec( const matrix_t *a, matrix_t **l, matrix_t **u, matrix_t **p )
{
	return (IS_EXIST( a ) && IS_SQUARE( a )) ?
		_matrix_lup_dec( a, l, u, p ) : NULL;
}

/**
 * @name Assignment functions
 *
 */

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_set( matrix_t *dest, const double src, const size_t row_index, const size_t col_index )
{
	return (IS_EXIST( dest ) && (row_index <= dest->i) && (col_index <= dest->j) ) ?
		_matrix_set( dest, src, row_index, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @return matrix_t*
 */
matrix_t *matrix_set_all( matrix_t *dest, const double *src, const size_t data_size )
{
	return (IS_EXIST( dest ) && IS_EXIST( src ) && (data_size <= dest->total)) ?
		_matrix_set_all( dest, src, data_size ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @param row_index
 * @return matrix_t*
 */
matrix_t *matrix_set_row( matrix_t *dest, const double *src, const size_t data_size, size_t row_index )
{
	return (IS_EXIST( dest ) && IS_EXIST( src ) && (data_size <= dest->j) && (row_index <= dest->i)) ?
		_matrix_set_row( dest, src, data_size, row_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_set_col( matrix_t *dest, const double *src, const size_t data_size, const size_t col_index )
{
	return (IS_EXIST( dest ) && IS_EXIST( src ) && (data_size <= dest->i) && (col_index <= dest->j)) ?
		_matrix_set_col( dest, src, data_size, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @return matrix_t*
 */
matrix_t *matrix_set_diag( matrix_t *dest, const double *src, const size_t data_size )
{
	return (IS_EXIST( dest ) && IS_EXIST( src ) && IS_SQUARE( dest ) && (data_size <= dest->i)) ?
		_matrix_set_diag( dest, src, data_size ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param data_size
 * @param ...
 * @return matrix_t*
 */
matrix_t *matrix_set_va( matrix_t *dest, size_t data_size, ... )
{
	matrix_t *result = NULL;
	va_list   ap;

	if ( IS_EXIST( dest ) && (data_size <= dest->total) ) {
		va_start(ap, data_size);
		result = _matrix_set_va( dest, data_size, ap );
		va_end(ap);
	}

	return result;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_get( const matrix_t *src, double *dest, const size_t row_index, const size_t col_index )
{
	return (IS_EXIST( src ) && IS_EXIST( dest ) && (row_index <= src->i) && (col_index <= src->j)) ?
		_matrix_get( src, dest, row_index, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @return matrix_t*
 */
matrix_t *matrix_get_all( const matrix_t *src, double *dest, const size_t dest_size )
{
	return (IS_EXIST( src ) && IS_EXIST( dest )) ?
		_matrix_get_all( src, dest, dest_size ) : NULL;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @param row_index
 * @return matrix_t*
 */
matrix_t *matrix_get_row( const matrix_t *src, double *dest, const size_t dest_size, const size_t row_index )
{
	return (IS_EXIST( src ) && IS_EXIST( dest ) && (row_index <= src->i)) ?
		_matrix_get_row( src, dest, dest_size, row_index ) : NULL;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_get_col( const matrix_t *src, double *dest, const size_t dest_size, const size_t col_index )
{
	return (IS_EXIST( src ) && IS_EXIST( dest ) && (col_index <= src->j)) ?
		_matrix_get_col( src, dest, dest_size, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @return matrix_t*
 */
matrix_t *matrix_get_diag( const matrix_t *src, double *dest, const size_t dest_size )
{
	return (IS_EXIST( src ) && IS_EXIST( dest ) && IS_SQUARE( src )) ?
		_matrix_get_diag( src, dest, dest_size ) : NULL;
}

/**
 * @name Appling functions
 *
 */

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_apply( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t row_index, const size_t col_index )
{
	return (IS_EXIST( dest ) && IS_EXIST( func ) && (row_index <= dest->i) && (col_index <= dest->j)) ?
		_matrix_apply( dest, func, arg, row_index, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @return matrix_t*
 */
matrix_t *matrix_apply_all( matrix_t *dest, double (*func)( double, void * ), void *arg )
{
	return (IS_EXIST( dest ) && IS_EXIST( func )) ? _matrix_apply_all( dest, func, arg ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param row_index
 * @return matrix_t*
 */
matrix_t *matrix_apply_row( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t row_index )
{
	return (IS_EXIST( dest ) && IS_EXIST( func ) && (row_index <= dest->i)) ?
		_matrix_apply_row( dest, func, arg, row_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param col_index
 * @return matrix_t*
 */
matrix_t *matrix_apply_col( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t col_index )
{
	return (IS_EXIST( dest ) && IS_EXIST( func ) && (col_index <= dest->j)) ?
		_matrix_apply_col( dest, func, arg, col_index ) : NULL;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @return matrix_t*
 */
matrix_t *matrix_apply_diag( matrix_t *dest, double (*func)( double, void * ), void *arg )
{
	return (IS_EXIST( dest ) && IS_EXIST( func ) && IS_SQUARE( dest )) ?
		_matrix_apply_diag( dest, func, arg ) : NULL;
}

/**
 * @brief
 *
 * @param a
 * @return size_t
 */
size_t matrix_rk( const matrix_t *a )
{
	return IS_EXIST( a ) ? _matrix_rk( a ) : 0;
}

/**
 * @brief
 *
 * @param a
 * @return double
 */
double matrix_det( const matrix_t *a )
{
	return (IS_EXIST( a ) && IS_SQUARE( a )) ? _matrix_det( a ) : NAN;
}

/**
 * @brief
 *
 * @param a
 * @return double
 */
double matrix_tr( const matrix_t *a )
{
	return (IS_EXIST( a ) && IS_SQUARE( a )) ? _matrix_tr( a ) : NAN;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @param tol
 * @return int
 */
int matrix_cmp( const matrix_t *a, const matrix_t *b, const double tol )
{
	return (IS_EXIST( a ) && IS_EXIST( b ) && ARE_SAME_SIZE( a, b )) ?
		_matrix_cmp( a, b, tol ) : -1;
}

/**
 * @brief
 *
 * @param a
 */
void matrix_free( matrix_t *a )
{
	if ( IS_EXIST( a ) )
		_matrix_free( a );

	return;
}

/**
 * @name
 *
 */

/**
 * @brief
 *
 * @param row
 * @param column
 * @return matrix_t*
 */
static matrix_t *_matrix_new( const size_t row, const size_t column )
{
	matrix_t *result = (matrix_t *)MatrixAllocFunc( sizeof(matrix_t) );
	size_t    elem_size;

/* */
	if ( result ) {
		result->i     = row;
		result->j     = column;
		result->total = row * column;
		elem_size     = result->total * sizeof(double);
	/* Then use the memory allocating function to allocate the elements' space */
		if ( !(result->element = (double *)MatrixAllocFunc( elem_size )) ) {
			MaxtrixFreeFunc( result );
			result = NULL;
		}
		else {
		/* Initialize the all the elements to zero */
			memset(result->element, 0, elem_size);
		}
	}

	return result;
}

/**
 * @brief
 *
 * @param rank
 * @return matrix_t*
 */
static matrix_t *_matrix_idt( const size_t rank )
{
	matrix_t    *result = _matrix_new( rank, rank );
	double      *elem_res;
	const size_t step = rank + 1;

	if ( result ) {
		elem_res = result->element;
		for ( register size_t i = 0; i < rank; i++, elem_res += step )
			*elem_res = 1.0;
	}

	return result;
}

/**
 * @brief
 *
 * @param a
 * @return matrix_t*
 */
static matrix_t *_matrix_dup( const matrix_t *a )
{
	matrix_t *result = _matrix_new( a->i, a->j );

	if ( result )
		memcpy(result->element, a->element, a->total * sizeof(double));

	return result;
}

/**
 * @brief Given a matrix a, return (create) the transpose of matrix a.
 *
 * @param a
 * @return matrix_t*
 */
static matrix_t *_matrix_tps( const matrix_t *a )
{
	matrix_t     *result = _matrix_new( a->j, a->i );
	double       *elem_res;
	const double *elem_a = a->element;

	if ( result ) {
		elem_res = result->element;
		for ( register size_t i = 0; i < a->i; i++, elem_res = result->element + i )
			for ( register size_t j = 0; j < a->j; j++, elem_res += a->i, elem_a++ )
				*elem_res = *elem_a;
	}

	return result;
}

/**
 * @brief Given a matrix a, return (create) the inverse of matrix a.
 *
 * @param a
 * @return matrix_t*
 */
static matrix_t *_matrix_inv( const matrix_t *a )
{
	matrix_t    *result = NULL;
	matrix_t    *tmp    = NULL;
	const size_t size_r = sizeof(double) * a->j;
	size_t       prow;
	double       pivot;
	double       swap_r[a->j];

/* */
	if ( (tmp = _matrix_dup( a )) && (result = _matrix_idt( a->i ))) {
	/* */
		for ( register size_t i = 0; i < a->i; i++ ) {
		/* */
			pivot = fabs(tmp->element[i * a->j + i]);
			prow  = i;
			for ( register size_t j = i + 1; j < a->i; j++ ) {
				if ( (swap_r[0] = fabs(tmp->element[j * a->j + i])) > pivot ) {
					pivot = swap_r[0];
					prow  = j;
				}
			}
		/* Failure, matrix is degenerate */
			if ( pivot <= DBL_EPSILON ) {
				_matrix_free( result );
				_matrix_free( tmp );
				return NULL;
			}
		/* */
			if ( prow != i ) {
				prow *= a->j;
			/* Swap the row for tmp matrix */
				memcpy(swap_r, tmp->element + (i * a->j), size_r);
				memcpy(tmp->element + (i * a->j), tmp->element + prow, size_r);
				memcpy(tmp->element + prow, swap_r, size_r);
			/* Swap the row for result matrix */
				memcpy(swap_r, result->element + (i * a->j), size_r);
				memcpy(result->element + (i * a->j), result->element + prow, size_r);
				memcpy(result->element + prow, swap_r, size_r);
			}
		/* */
			pivot = tmp->element[i * a->j + i];
			for ( register size_t j = 0, p_idx = i * a->j; j < a->j; j++, p_idx++ ) {
			/* */
				if ( fabs(tmp->element[p_idx]) > DBL_EPSILON )
					tmp->element[p_idx] /= pivot;
			/* */
				if ( fabs(result->element[p_idx]) > DBL_EPSILON )
					result->element[p_idx] /= pivot;
			}
		/* Eliminate the lower triangular */
			for ( register size_t j = i + 1; j < a->i; j++ ) {
				pivot = tmp->element[j * a->j + i];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( register size_t k = 0, t_idx = j * a->j, p_idx = i * a->j; k < a->j; k++, t_idx++, p_idx++ ) {
					/* */
						if ( fabs(tmp->element[p_idx]) > DBL_EPSILON )
							tmp->element[t_idx] -= pivot * tmp->element[p_idx];
					/* */
						if ( fabs(result->element[p_idx]) > DBL_EPSILON )
							result->element[t_idx] -= pivot * result->element[p_idx];
					}
				}
			}
		}
	/* Eliminate the upper triangular */
		for ( register size_t i = a->i - 1; i > 0; i-- ) {
			for ( register size_t j = 0; j < i; j++ ) {
				pivot = tmp->element[j * a->j + i];
				if ( fabs(pivot) > DBL_EPSILON ) {
					for ( register size_t k = 0, t_idx = j * a->j, p_idx = i * a->j; k < a->j; k++, t_idx++, p_idx++ ) {
					/* */
						if ( fabs(tmp->element[p_idx]) > DBL_EPSILON )
							tmp->element[t_idx] -= pivot * tmp->element[p_idx];
					/* */
						if ( fabs(result->element[p_idx]) > DBL_EPSILON )
							result->element[t_idx] -= pivot * result->element[p_idx];
					}
				}
			}
		}
	}
/* */
	if ( tmp )
		_matrix_free( tmp );

	return result;
}

/**
 * @name
 *
 */

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
static matrix_t *_matrix_add( const matrix_t *a, const matrix_t *b )
{
	matrix_t     *result = NULL;
	double       *elem_res;
	const double *elem_b = b->element;

	if ( (result = _matrix_dup( a )) ) {
		elem_res = result->element;
		for ( register size_t i = 0; i < result->total; i++, elem_res++, elem_b++  )
			*elem_res += *elem_b;
	}

	return result;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
static matrix_t *_matrix_sub( const matrix_t *a, const matrix_t *b )
{
	matrix_t     *result = NULL;
	double       *elem_res;
	const double *elem_b = b->element;

	if ( (result = _matrix_dup( a )) ) {
		elem_res = result->element;
		for ( register size_t i = 0; i < result->total; i++, elem_res++, elem_b++ )
			*elem_res -= *elem_b;
	}

	return result;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return matrix_t*
 */
static matrix_t *_matrix_mul( const matrix_t *a, const matrix_t *b )
{
/* */
#ifdef __USE_AVX_INTRIN
	matrix_t *result = NULL;
	__m256d   elem_a_vec, row_vec[STEP_AVX_PD], sum_vec[STEP_AVX_PD];
	__m256i   mask;

	if ( (result = _matrix_new( a->i, b->j )) ) {
		for ( register size_t i = 0; i < a->i; i += STEP_AVX_PD ) {
			for ( register size_t j = 0; j < b->j; j += STEP_AVX_PD) {
			/* Set to zero */
				for ( register int k = 0; k < STEP_AVX_PD; k++ )
					sum_vec[k] = _mm256_setzero_pd();
			/* */
				for ( register int k = 0; k < a->j; k += STEP_AVX_PD ) {
				/* Load the matrix b elements */
					for ( register int l = 0; l < STEP_AVX_PD; l++ ) {
						if ( (j + STEP_AVX_PD) <= b->j ) {
							row_vec[l] = _mm256_loadu_pd(b->element + (k + l) * b->j + j);
						}
						else {
							int64_t _mask[STEP_AVX_PD];
							for ( register int ll = 0; ll < STEP_AVX_PD; ll++ )
								_mask[ll] = (j + ll < b->j) ? -1 : 0;
							mask       = _mm256_setr_epi64x(_mask[0], _mask[1], _mask[2], _mask[3]);
							row_vec[l] = _mm256_maskload_pd(b->element + (k + l) * b->j + j, mask);
						}
					}

					for ( register int l = 0; l < STEP_AVX_PD && (i + l) < a->i; l++ ) {
						for ( register int ll = 0; ll < STEP_AVX_PD && (k + ll) < a->j; ll++ ) {
							elem_a_vec = _mm256_set1_pd(a->element[(i + l) * a->j + k + ll]);
							sum_vec[l] = _mm256_fmadd_pd(row_vec[ll], elem_a_vec, sum_vec[l]);
						}
					}
				}

				for ( register int k = 0; k < STEP_AVX_PD && i + k < a->i; k++ ) {
					if ( (j + STEP_AVX_PD) <= b->j )
						_mm256_storeu_pd(result->element + (i + k) * b->j + j, sum_vec[k]);
					else
						_mm256_maskstore_pd(result->element + (i + k) * b->j + j, mask, sum_vec[k]);
				}
			}
		}
	}
#else
	matrix_t     *result = NULL;
	double       *elem_res;
	const double *elem_a = a->element;
	const double *elem_b;

	if ( (result = _matrix_new( a->i, b->j )) ) {
		for ( register size_t i = 0; i < a->i; i++ ) {
			for ( register size_t j = 0; j < a->j; j++, elem_a++ ) {
				elem_b   = b->element + j * b->j;
				elem_res = result->element + i * b->j;
				if ( fabs(*elem_a) > DBL_EPSILON ) {
					for ( register size_t k = 0; k < b->j; k++, elem_b++, elem_res++ ) {
						if ( fabs(*elem_b) > DBL_EPSILON )
						#ifdef __USE_FMA_INTRIN
							*elem_res = fma(*elem_a, *elem_b, *elem_res);
						#else
							*elem_res += (*elem_a) * (*elem_b);
						#endif
					}
				}
			}
		}
	}
#endif

	return result;
}

/**
 * @brief Return (create) the result Weighted Least Squares(WLS) matrix from 'W * A * x = W * B',
 *        where x should be '(A^T * W * A + Opt)^-1 * A^T * W * B'.
 *
 * @param a
 * @param b
 * @param w
 * @param opt
 * @return matrix_t*
 */
static matrix_t *_matrix_wls( const matrix_t *a, const matrix_t *b, const matrix_t *w, const matrix_t *opt )
{
	matrix_t *matrix_1 = NULL;
	matrix_t *matrix_2 = NULL;
	matrix_t *result = NULL;
	double    tmp;

/* Find the transpose matrix of A */
	if ( w && IS_SQUARE( w ) && w->i == a->i ) {
	/* Here, if we got weighting matrix, just multiply with A transpose */
		if ( !(matrix_1 = _matrix_dup( a )) )
			return NULL;
		for ( register size_t i = 0; i < matrix_1->i; i++ ) {
			tmp = w->element[i * matrix_1->j + i];
			for ( register size_t j = 0; j < matrix_1->j; j++ )
				matrix_1->element[i * matrix_1->j + j] *= tmp;
		}
		result = _matrix_tps( matrix_1 );
		_matrix_free( matrix_1 );
	}
	else {
		result = _matrix_tps( a );
	}
/* */
	if ( !result )
		return NULL;
/* Multiply A transpose with A and B */
	if ( (matrix_1 = _matrix_mul( result, a )) ) {
		if ( (matrix_2 = _matrix_mul( result, b )) ) {
		/* */
			_matrix_free( result );
		/* If we got optimaztion matrix, add it to the A^T*A */
			if ( opt && ARE_SAME_SIZE( matrix_1, opt ) ) {
				result = _matrix_add( matrix_1, opt );
				_matrix_free( matrix_1 );
			}
			else {
			/* Or we should directly use the A^T*A matrix */
				result = matrix_1;
				matrix_1 = NULL;
			}
		/* Then get the inverse of (A^T * A + Opt) */
			if ( result && (matrix_1 = _matrix_inv( result )) ) {
				_matrix_free( result );
				result = _matrix_mul( matrix_1, matrix_2 );
				_matrix_free( matrix_1 );
				_matrix_free( matrix_2 );

				return result;
			}
		/* */
			if ( result )
				_matrix_free( result );
		/* */
			_matrix_free( matrix_2 );

			return NULL;
		}
	/* */
		_matrix_free( matrix_1 );
	}
/* */
	_matrix_free( result );

	return NULL;
}

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
static matrix_t *_matrix_scalar_add( matrix_t *a, const double scalar )
{
	double *elem_a = a->element;

	for ( register size_t i = 0; i < a->total; i++, elem_a++ )
		*elem_a += scalar;

	return a;
}

/**
 * @brief
 *
 * @param a
 * @param scalar
 * @return matrix_t*
 */
static matrix_t *_matrix_scalar_mul( matrix_t *a, const double scalar )
{
	double *elem_a = a->element;

	for ( register size_t i = 0; i < a->total; i++, elem_a++ )
		if ( fabs(*elem_a) > DBL_EPSILON )
			*elem_a *= scalar;

	return a;
}

/**
 * @brief
 *
 * @param a
 * @param l
 * @param u
 * @param p
 * @return matrix_t*
 */
static matrix_t *_matrix_lup_dec( const matrix_t *a, matrix_t **l, matrix_t **u, matrix_t **p )
{
	matrix_t    *result = NULL;
	matrix_t    *_p     = NULL;
	matrix_t    *_l     = NULL;
	matrix_t    *_u     = NULL;
	const size_t size_r = sizeof(double) * a->j;
	size_t       prow;
	double       pivot;
	double       swap_r[a->j];

/* */
	if ( (_p = _matrix_new( a->i + 1, 1 )) && (result = _matrix_dup( a )) ) {
	/* */
		for ( register size_t i = 0; i < a->i; i++ )
			_p->element[i] = i;
	/* */
		for ( register size_t i = 0; i < a->i; i++ ) {
		/* Find the pivoting row */
			pivot = fabs(result->element[i * a->j + i]);
			prow  = i;
			for ( register size_t j = i + 1; j < a->i; j++ ) {
				if ( (swap_r[0] = fabs(result->element[j * a->j + i])) > pivot ) {
					pivot = swap_r[0];
					prow  = j;
				}
			}
		/* Failure, matrix is degenerate */
			if ( pivot <= DBL_EPSILON ) {
				_matrix_free( result );
				_matrix_free( _p );
				return NULL;
			}
		/* Do the permutation */
			if ( prow != i ) {
			/* Swap the row for permutation matrix */
				pivot             = _p->element[i];
				_p->element[i]    = _p->element[prow];
				_p->element[prow] = pivot;
			/* Swap the row for result matrix */
				memcpy(swap_r, result->element + (i * a->j), size_r);
				memcpy(result->element + (i * a->j), result->element + (prow * a->j), size_r);
				memcpy(result->element + (prow * a->j), swap_r, size_r);
			/* */
				_p->element[a->i] += 1.0;
			}
		/* Do the elimination */
			for ( register size_t j = i + 1; j < a->i; j++ ) {
			/* Pivoting, store the L matrix elements */
				pivot = (result->element[j * a->j + i] /= result->element[i * a->j + i]);
			/* Eliminate the rest elements, store the U matrix elements */
				if ( fabs(pivot) > DBL_EPSILON )
					for ( register size_t k = i + 1; k < a->j; k++ )
						if ( fabs(result->element[i * a->j + k]) > DBL_EPSILON )
							result->element[j * a->j + k] -= pivot * result->element[i * a->j + k];
			}
		}
	/* */
		if ( p )
			*p = _p;
		else
			_matrix_free( _p );
	/* Extract the L matrix elements */
		if ( l && (_l = _matrix_idt( a->i )) ) {
			for ( register size_t i = 1; i < a->i; i++ )
				for ( register size_t j = 0; j < i; j++ )
					_l->element[i * a->j + j] = result->element[i * a->j + j];
		/* */
			*l = _l;
		}
	/* Extract the U matrix elements */
		if ( u && (_u = _matrix_new( a->i, a->j )) ) {
			for ( register size_t i = 0; i < a->i; i++ )
				for ( register size_t j = i; j < a->j; j++ )
					_u->element[i * a->j + j] = result->element[i * a->j + j];
		/* */
			*u = _u;
		}
	}
	else if ( _p ) {
		_matrix_free( _p );
	}

	return result;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_set( matrix_t *dest, const double src, const size_t row_index, const size_t col_index )
{
	dest->element[(row_index - 1) * dest->j + (col_index - 1)] = src;

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @return matrix_t*
 */
static matrix_t *_matrix_set_all( matrix_t *dest, const double *src, const size_t data_size )
{
	double *elem_dest = dest->element;

/* */
	memcpy(elem_dest, src, data_size * sizeof(double));
/* */
	elem_dest += data_size;
	for ( register size_t i = data_size; i < dest->total; i++, elem_dest++ )
		*elem_dest = 0.0;

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @param row_index
 * @return matrix_t*
 */
static matrix_t *_matrix_set_row( matrix_t *dest, const double *src, const size_t data_size, const size_t row_index )
{
	double *elem_dest = dest->element + (row_index - 1) * dest->j;

/* */
	memcpy(elem_dest, src, data_size * sizeof(double));
/* */
	elem_dest += data_size;
	for ( register size_t i = data_size; i < dest->j; i++, elem_dest++ )
		*elem_dest = 0.0;

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_set_col( matrix_t *dest, const double *src, const size_t data_size, const size_t col_index )
{
	double *elem_dest = dest->element + (col_index - 1);

/* */
	for ( register size_t i = 0; i < dest->i; i++, elem_dest += dest->j )
		*elem_dest = i < data_size ? *src++ : 0.0;

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param src
 * @param data_size
 * @return matrix_t*
 */
static matrix_t *_matrix_set_diag( matrix_t *dest, const double *src, const size_t data_size )
{
	double      *elem_dest = dest->element;
	const size_t step      = dest->j + 1;

/* */
	for ( register size_t i = 0; i < dest->i; i++, elem_dest += step )
		*elem_dest = i < data_size ? *src++ : 0.0;

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param data_size
 * @param ap
 * @return matrix_t*
 */
static matrix_t *_matrix_set_va( matrix_t *dest, size_t data_size, va_list ap )
{
	double *elem_dest = dest->element;

/* */
	for ( ; data_size > 0; data_size--, elem_dest++ )
		*elem_dest = va_arg(ap, double);

	return dest;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_get( const matrix_t *src, double *dest, const size_t row_index, const size_t col_index )
{
/* */
	*dest = src->element[(row_index - 1) * src->j + (col_index - 1)];

	return (matrix_t *)src;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @return matrix_t*
 */
static matrix_t *_matrix_get_all( const matrix_t *src, double *dest, const size_t dest_size )
{
/* */
	memcpy(dest, src->element, (dest_size < src->total ? dest_size : src->total) * sizeof(double));
/* Fill zero to the rest space of dest */
	for ( register size_t i = src->total; i < dest_size; i++ )
		dest[i] = 0.0;
/* */
	return (matrix_t *)src;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @param row_index
 * @return matrix_t*
 */
static matrix_t *_matrix_get_row( const matrix_t *src, double *dest, const size_t dest_size, const size_t row_index )
{
/* */
	memcpy(dest, src->element + (row_index - 1) * src->j, (dest_size < src->j ? dest_size : src->j) * sizeof(double));
/* Fill zero to the rest space of dest */
	for ( register size_t i = src->j; i < dest_size; i++ )
		dest[i] = 0.0;
/* */
	return (matrix_t *)src;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_get_col( const matrix_t *src, double *dest, const size_t dest_size, const size_t col_index )
{
	const double *elem_src = src->element + (col_index - 1);

/* */
	for ( register size_t i = 0; i < src->i && i < dest_size; i++, dest++, elem_src += src->j )
		*dest = *elem_src;
/* Fill zero to the rest space of dest */
	for ( register size_t i = src->i; i < dest_size; i++, dest++ )
		*dest = 0.0;
/* */
	return (matrix_t *)src;
}

/**
 * @brief
 *
 * @param src
 * @param dest
 * @param dest_size
 * @return matrix_t*
 */
static matrix_t *_matrix_get_diag( const matrix_t *src, double *dest, const size_t dest_size )
{
	const double *elem_src = src->element;
	const size_t  step = src->j + 1;

/* */
	for ( register size_t i = 0; i < src->i && i < dest_size; i++, dest++, elem_src += step )
		*dest = *elem_src;
/* Fill zero to the rest space of dest */
	for ( register size_t i = src->i; i < dest_size; i++, dest++ )
		*dest = 0.0;
/* */
	return (matrix_t *)src;
}

/**
 * @name Appling functions
 *
 */

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param row_index
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_apply( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t row_index, const size_t col_index )
{
	double *elem_dest = dest->element + (row_index - 1) * dest->j + (col_index - 1);

/* */
	*elem_dest = func( *elem_dest, arg );

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @return matrix_t*
 */
static matrix_t *_matrix_apply_all( matrix_t *dest, double (*func)( double, void * ), void *arg )
{
	double *elem_dest = dest->element;

/* */
	for ( register size_t i = 0; i < dest->total; i++, elem_dest++ )
		*elem_dest = func( *elem_dest, arg );

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param row_index
 * @return matrix_t*
 */
static matrix_t *_matrix_apply_row( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t row_index )
{
	double *elem_dest = dest->element + (row_index - 1) * dest->j;

/* */
	for ( register size_t i = 0; i < dest->j; i++, elem_dest++ )
		*elem_dest = func( *elem_dest, arg );

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @param col_index
 * @return matrix_t*
 */
static matrix_t *_matrix_apply_col( matrix_t *dest, double (*func)( double, void * ), void *arg, const size_t col_index )
{
	double *elem_dest = dest->element + (col_index - 1);

/* */
	for ( register size_t i = 0; i < dest->i; i++, elem_dest += dest->j )
		*elem_dest = func( *elem_dest, arg );

	return dest;
}

/**
 * @brief
 *
 * @param dest
 * @param func
 * @param arg
 * @return matrix_t*
 */
static matrix_t *_matrix_apply_diag( matrix_t *dest, double (*func)( double, void * ), void *arg )
{
	double      *elem_dest = dest->element;
	const size_t step = dest->j + 1;

	for ( register size_t i = 0; i < dest->i; i++, elem_dest += step )
		*elem_dest = func( *elem_dest, arg );

	return dest;
}

/**
 * @brief
 *
 * @param a
 * @return size_t
 */
static size_t _matrix_rk( const matrix_t *a )
{
	matrix_t    *u      = NULL;
	const size_t size_r = sizeof(double) * a->j;
	size_t       result = 0;
	size_t       r_idx  = 0;
	size_t       c_idx  = 0;
	size_t       prow;
	double       pivot;
	double       swap_r[a->j];

/* */
	if ( (u = matrix_dup( a )) ) {
	/* Upper triangularize */
		while ( r_idx < a->i && c_idx < a->j ) {
		/* Find the pivoting row */
			pivot = fabs(u->element[r_idx * a->j + c_idx]);
			prow  = r_idx;
			for ( register size_t i = r_idx + 1; i < a->i; i++ ) {
				if ( (swap_r[0] = fabs(u->element[i * a->j + c_idx])) > pivot ) {
					pivot = swap_r[0];
					prow  = i;
				}
			}
		/* */
			if ( pivot <= DBL_EPSILON ) {
			/* */
				c_idx++;
				continue;
			}
		/* Do the permutation */
			if ( prow != r_idx ) {
			/* */
				prow *= a->j;
			/* Swap the row for result matrix */
				memcpy(swap_r, u->element + (r_idx * a->j), size_r);
				memcpy(u->element + (r_idx * a->j), u->element + prow, size_r);
				memcpy(u->element + prow, swap_r, size_r);
			}
		/* Do the elimination */
			for ( register size_t i = r_idx + 1; i < a->i; i++ ) {
			/* Pivoting */
				if ( fabs(pivot = u->element[i * a->j + c_idx] / u->element[r_idx * a->j + c_idx]) > DBL_EPSILON )
				/* Eliminate the rest elements */
					for ( register size_t k = c_idx + 1; k < a->j; k++ )
						if ( fabs(u->element[r_idx * a->j + k]) > DBL_EPSILON )
							u->element[i * a->j + k] -= pivot * u->element[r_idx * a->j + k];
			/* Set the lower triangular part to zero */
				u->element[i * a->j + c_idx] = 0.0;
			}
		/* */
			r_idx++;
			c_idx++;
		}

	/* */
		for ( register size_t i = 0; i < a->i; i++ ) {
			for ( register size_t j = 0; j < a->j; j++ ) {
				if ( fabs(u->element[i * a->j + j]) > DBL_EPSILON ) {
					result++;
					break;
				}
			}
		}
	/* */
		_matrix_free( u );
	}

	return result;
}

/**
 * @brief
 *
 * @param a
 * @return double
 */
static double _matrix_det( const matrix_t *a )
{
	matrix_t     *p      = NULL;
	matrix_t     *lu     = NULL;
	double        result = 0.0;
	const double *elem_lu;
	const size_t  step   = a->j + 1;

/* */
	if ( (lu = matrix_lup_dec( a, NULL, NULL, &p )) ) {
		elem_lu = lu->element;
	/* */
		p->element[p->total - 1] += p->element[p->total - 1] - (int)p->element[p->total - 1];
		result = ((int)p->element[p->total - 1] % 2) ? -*elem_lu : *elem_lu;
	/* */
		for ( register size_t i = 1; i < a->i; i++ )
			result *= *(elem_lu += step);
	/* */
		_matrix_free( lu );
		_matrix_free( p );
	}

	return result;
}

/**
 * @brief
 *
 * @param a
 * @return double
 */
static double _matrix_tr( const matrix_t *a )
{
	const double *elem_a = a->element;
	const size_t  step   = a->j + 1;
	double        result = *elem_a;

	for ( register size_t i = 1; i < a->i; i++ )
		result += *(elem_a += step);

	return result;
}


/**
 * @brief
 *
 * @param a
 * @param b
 * @param tol
 * @return int
 */
static int _matrix_cmp( const matrix_t *a, const matrix_t *b, const double tol )
{
	const double *elem_a;
	const double *elem_b;

/* */
	elem_a = a->element;
	elem_b = b->element;
/* */
	if ( tol > DBL_EPSILON ) {
		for ( register size_t i = 0; i < a->total; i++ ) {
			if ( fabs(*elem_a++ - *elem_b++) > tol )
				return 1;
		}
	}
	else if ( memcmp(elem_a, elem_b, a->total * sizeof(double)) ) {
		return 1;
	}

/* */
	return 0;
}

/**
 * @brief
 *
 * @param a
 */
static void _matrix_free( matrix_t *a )
{
/* */
	if ( a->element )
		MaxtrixFreeFunc( a->element );
/* */
	MaxtrixFreeFunc( a );

	return;
}
