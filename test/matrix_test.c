/**
 * @file matrix_test.c
 * @author Benjamin Ming Yang @ Department of Geoscience, National Taiwan University (b98204032@gmail.com)
 * @brief
 * @date 2024-09-06
 *
 * @copyright Copyright (c) 2024
 *
 */

/**
 * @name
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

/**
 * @name
 *
 */
#include "./munit/munit.h"
#include "../include/matrix.h"

/**
 * @brief
 *
 */
typedef struct {
	int        num;
	matrix_t **input;
	void      *answer;
} TestData;

/**
 * @name
 *
 */
static void test_matrix_matrix_tear_down( void * );
static void test_matrix_scalar_tear_down( void * );

static MunitResult test_matrix_dup( const MunitParameter [], void * );
static void *test_matrix_dup_setup( const MunitParameter [], void * );
static MunitResult test_matrix_tps( const MunitParameter [], void * );
static void *test_matrix_tps_setup( const MunitParameter [], void * );
static MunitResult test_matrix_inv( const MunitParameter [], void * );
static void *test_matrix_inv_setup( const MunitParameter [], void * );

static MunitResult test_matrix_det( const MunitParameter [], void * );
static void *test_matrix_det_setup( const MunitParameter [], void * );
static MunitResult test_matrix_rk( const MunitParameter [], void * );
static void *test_matrix_rk_setup( const MunitParameter [], void * );
static MunitResult test_matrix_tr( const MunitParameter [], void * );
static void *test_matrix_tr_setup( const MunitParameter [], void * );

/**
 * @brief
 *
 */
static MunitTest test_suite_tests[] = {
	{
		(char *)"/matrix_dup",
		test_matrix_dup,
		test_matrix_dup_setup,
		test_matrix_matrix_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{
		(char *)"/matrix_tps",
		test_matrix_tps,
		test_matrix_tps_setup,
		test_matrix_matrix_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{
		(char *)"/matrix_inv",
		test_matrix_inv,
		test_matrix_inv_setup,
		test_matrix_matrix_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{
		(char *)"/matrix_det",
		test_matrix_det,
		test_matrix_det_setup,
		test_matrix_scalar_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{
		(char *)"/matrix_rk",
		test_matrix_rk,
		test_matrix_rk_setup,
		test_matrix_scalar_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{
		(char *)"/matrix_tr",
		test_matrix_tr,
		test_matrix_tr_setup,
		test_matrix_scalar_tear_down,
		MUNIT_TEST_OPTION_NONE, NULL
	},
	{ NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

/**
 * @brief
 *
 */
static const MunitSuite test_suite = {
	(char*) "libmatrix",
	test_suite_tests,
	NULL,
	1,
	MUNIT_SUITE_OPTION_NONE
};

/**
 * @brief
 *
 * @return int
 */
int main( int argc, char* argv[MUNIT_ARRAY_PARAM(argc + 1)] )
{
	// matrix_t *a = matrix_new( 1024, 1024 );

	// struct timespec t1;
	// struct timespec t2;

	// for ( int i = 0; i < a->total; i++ )
	// 	a->element[i] = -100.0 + 200.0 * (rand() / (double)RAND_MAX);
	// timespec_get(&t1, TIME_UTC);
	// timespec_get(&t2, TIME_UTC);
	// printf("Timer1: %lf\n", (double)(t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec) * 1e-9);

	return munit_suite_main(&test_suite, NULL, argc, argv);
}

/**
 * @brief
 *
 * @param fixture
 */
static void test_matrix_scalar_tear_down( void * fixture )
{
	TestData *data = (TestData *)fixture;

/* */
	for ( int i = 0; i < data->num; i++ )
		matrix_free( data->input[i] );
/* */
	free(data->input);
/* */
	free(data->answer);

	return;
}

/**
 * @brief
 *
 * @param fixture
 */
static void test_matrix_matrix_tear_down( void * fixture )
{
	TestData *data = (TestData *)fixture;

/* */
	for ( int i = 0; i < data->num; i++ )
		matrix_free( data->input[i] );
/* */
	free(data->input);
/* */
	for ( int i = 0; i < data->num; i++ )
		matrix_free( ((matrix_t **)data->answer)[i] );
/* */
	free(data->answer);

	return;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_det( const MunitParameter params[], void *data )
{
	TestData *_data = (TestData *)data;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_double_equal(matrix_det( _data->input[i] ), ((double *)_data->answer)[i], 12);
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_det_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 8;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(double) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	((double *)result->answer)[0] = 46.0;
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	((double *)result->answer)[1] = -106.0;
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	((double *)result->answer)[2] = 1270.0;
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	((double *)result->answer)[3] = -56964.0;
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 0.0,
		-6.0, 7.0, 0.0, -8.0, 0.0,
		-2.0, 11.0, 9.0, 4.0, 0.0,
		0.0, 2.0, 0.0, 1.0, 0.0,
		-2.0, 10.0, -9.0, 4.0, 0.0
	));
	((double *)result->answer)[4] = 0.0;
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	((double *)result->answer)[5] = NAN;
/* */
	munit_assert_not_null(result->input[6] = matrix_idt( 16 ));
	((double *)result->answer)[6] = 1.0;
/* */
	munit_assert_null(result->input[7] = NULL);
	((double *)result->answer)[7] = NAN;

	return result;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_rk( const MunitParameter params[], void *data )
{
	TestData *_data = (TestData *)data;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_ulong(matrix_rk( _data->input[i] ), ==, ((size_t *)_data->answer)[i]);
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_rk_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 10;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(size_t) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	((size_t *)result->answer)[0] = 2;
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	((size_t *)result->answer)[1] = 3;
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	((size_t *)result->answer)[2] = 4;
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	((size_t *)result->answer)[3] = 5;
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 0.0,
		-6.0, 7.0, 0.0, -8.0, 0.0,
		-2.0, 11.0, 9.0, 4.0, 0.0,
		0.0, 2.0, 0.0, 1.0, 0.0,
		-2.0, 10.0, -9.0, 4.0, 0.0
	));
	((size_t *)result->answer)[4] = 4;
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	((size_t *)result->answer)[5] = 4;
/* */
	munit_assert_not_null(result->input[6] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[6], result->input[6]->total,
		4.0, 3.0, 0.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 0.0, 4.0, 3.0,
		8.0, 6.0, 0.0, 0.0, 10.0
	));
	((size_t *)result->answer)[6] = 3;
/* */
	munit_assert_not_null(result->input[7] = matrix_new( 4, 5 ));
	((size_t *)result->answer)[7] = 0;
/* */
	munit_assert_not_null(result->input[8] = matrix_idt( 16 ));
	((size_t *)result->answer)[8] = 16;
/* */
	munit_assert_null(result->input[9] = NULL);
	((size_t *)result->answer)[9] = 0;

	return result;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_tr( const MunitParameter params[], void *data )
{
	TestData *_data = (TestData *)data;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_double_equal(matrix_tr( _data->input[i] ), ((double *)_data->answer)[i], 12);
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_tr_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 10;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(double) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	((double *)result->answer)[0] = 11.0;
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	((double *)result->answer)[1] = 20.0;
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	((double *)result->answer)[2] = 21.0;
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	((double *)result->answer)[3] = 17.0;
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 0.0,
		-6.0, 7.0, 0.0, -8.0, 0.0,
		-2.0, 11.0, 9.0, 4.0, 0.0,
		0.0, 2.0, 0.0, 1.0, 0.0,
		-2.0, 10.0, -9.0, 4.0, 0.0
	));
	((double *)result->answer)[4] = 21.0;
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	((double *)result->answer)[5] = NAN;
/* */
	munit_assert_not_null(result->input[6] = matrix_new( 5, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[6], result->input[6]->total,
		4.0, 3.0, 0.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 0.0, 4.0,
		8.0, 6.0, 0.0, 0.0,
		8.0, 6.0, 0.0, 0.0
	));
	((double *)result->answer)[6] = NAN;
/* */
	munit_assert_not_null(result->input[7] = matrix_new( 4, 4 ));
	((double *)result->answer)[7] = 0.0;
/* */
	munit_assert_not_null(result->input[8] = matrix_idt( 16 ));
	((double *)result->answer)[8] = 16.0;
/* */
	munit_assert_null(result->input[9] = NULL);
	((double *)result->answer)[9] = 0;

	return result;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_dup( const MunitParameter params[], void *data )
{
	TestData *_data   = (TestData *)data;
	matrix_t *mtx_dup = NULL;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_true(matrix_cmp( (mtx_dup = matrix_dup( _data->input[i] )), ((matrix_t **)_data->answer)[i], 0.0 ) == 0);
		matrix_free( mtx_dup );
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_dup_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 9;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(matrix_t *) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[0], ((matrix_t **)result->answer)[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[1], ((matrix_t **)result->answer)[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[2], ((matrix_t **)result->answer)[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[3], ((matrix_t **)result->answer)[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[4] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[4], ((matrix_t **)result->answer)[4]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 5, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, -6.0, -2.0, 0.0,
		3.0, 7.0, 11.0, 2.0,
		10.0, 0.0, 9.0, 0.0,
		0.0, -8.0, 4.0, 1.0,
		5.0, 28.0, 3.0, 5.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[5] = matrix_new( 5, 4 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[5], ((matrix_t **)result->answer)[5]->total,
		4.0, -6.0, -2.0, 0.0,
		3.0, 7.0, 11.0, 2.0,
		10.0, 0.0, 9.0, 0.0,
		0.0, -8.0, 4.0, 1.0,
		5.0, 28.0, 3.0, 5.0
	));
/* */
	munit_assert_not_null(result->input[6] = matrix_new( 4, 4 ));
	munit_assert_not_null(((matrix_t **)result->answer)[6] = matrix_new( 4, 4 ));
/* */
	munit_assert_not_null(result->input[7] = matrix_idt( 16 ));
	munit_assert_not_null(((matrix_t **)result->answer)[7] = matrix_idt( 16 ));
/* */
	munit_assert_null(result->input[8] = NULL);
	((matrix_t **)result->answer)[8] = NULL;

	return result;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_tps( const MunitParameter params[], void *data )
{
	TestData *_data   = (TestData *)data;
	matrix_t *mtx_tps = NULL;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_true(matrix_cmp( (mtx_tps = matrix_tps( _data->input[i] )), ((matrix_t **)_data->answer)[i], 0.0 ) == 0);
		matrix_free( mtx_tps );
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_tps_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 9;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(matrix_t *) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[0], ((matrix_t **)result->answer)[0]->total,
		4.0, -6.0,
		3.0, 7.0
	));
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[1], ((matrix_t **)result->answer)[1]->total,
		4.0, -6.0, -2.0,
		3.0, 7.0, 11.0,
		10.0, 0.0, 9.0
	));
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[2], ((matrix_t **)result->answer)[2]->total,
		4.0, -6.0, -2.0, 0.0,
		3.0, 7.0, 11.0, 2.0,
		10.0, 0.0, 9.0, 0.0,
		0.0, -8.0, 4.0, 1.0
	));
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[3], ((matrix_t **)result->answer)[3]->total,
		4.0, -6.0, -2.0, 0.0, -2.0,
		3.0, 7.0, 11.0, 2.0, 10.0,
		10.0, 0.0, 9.0, 0.0, -9.0,
		0.0, -8.0, 4.0, 1.0, 4.0,
		5.0, 28.0, 3.0, 5.0, -4.0
	));
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[4] = matrix_new( 5, 4 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[4], ((matrix_t **)result->answer)[4]->total,
		4.0, -6.0, -2.0, 0.0,
		3.0, 7.0, 11.0, 2.0,
		10.0, 0.0, 9.0, 0.0,
		0.0, -8.0, 4.0, 1.0,
		5.0, 28.0, 3.0, 5.0
	));
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 5, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, -6.0, -2.0, 0.0,
		3.0, 7.0, 11.0, 2.0,
		10.0, 0.0, 9.0, 0.0,
		0.0, -8.0, 4.0, 1.0,
		5.0, 28.0, 3.0, 5.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[5] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[5], ((matrix_t **)result->answer)[5]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
/* */
	munit_assert_not_null(result->input[6] = matrix_new( 4, 4 ));
	munit_assert_not_null(((matrix_t **)result->answer)[6] = matrix_new( 4, 4 ));
/* */
	munit_assert_not_null(result->input[7] = matrix_idt( 16 ));
	munit_assert_not_null(((matrix_t **)result->answer)[7] = matrix_idt( 16 ));
/* */
	munit_assert_null(result->input[8] = NULL);
	((matrix_t **)result->answer)[8] = NULL;

	return result;
}

/**
 * @brief
 *
 * @param params
 * @param data
 * @return MunitResult
 */
static MunitResult test_matrix_inv( const MunitParameter params[], void *data )
{
	TestData *_data   = (TestData *)data;
	matrix_t *mtx_inv = NULL;

/*
 * These are just to silence compiler warnings about the parameters
 * being unused.
 */
	(void)params;
/* */
	//munit_logf(MUNIT_LOG_INFO, "Total testing data is %d\n", _data->num);
/* Real testing */
	for ( int i = 0; i < _data->num; i++ ) {
		munit_assert_true(matrix_cmp( (mtx_inv = matrix_inv( _data->input[i] )), ((matrix_t **)_data->answer)[i], 1e-12 ) == 0);
		matrix_free( mtx_inv );
	}

	return MUNIT_OK;
}

/**
 * @brief
 *
 * @param params
 * @param user_data
 * @return void*
 */
static void *test_matrix_inv_setup( const MunitParameter params[], void *user_data )
{
	TestData *result = munit_malloc(sizeof(TestData));

	(void)params;

/* */
	result->num    = 9;
	result->input  = munit_malloc(sizeof(matrix_t *) * result->num);
	result->answer = munit_malloc(sizeof(matrix_t *) * result->num);
/* */
	munit_assert_not_null(result->input[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( result->input[0], result->input[0]->total,
		4.0, 3.0,
		-6.0, 7.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[0] = matrix_new( 2, 2 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[0], ((matrix_t **)result->answer)[0]->total,
		7.0 / 46.0, -3.0 / 46.0,
		3.0 / 23.0, 2.0 / 23.0
	));
/* */
	munit_assert_not_null(result->input[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( result->input[1], result->input[1]->total,
		4.0, 3.0, 10.0,
		-6.0, 7.0, 0.0,
		-2.0, 11.0, 9.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[1] = matrix_new( 3, 3 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[1], ((matrix_t **)result->answer)[1]->total,
		-63.0 / 106.0, -83.0 / 106.0, 35.0 / 53.0,
		-27.0 / 53.0, -28.0 / 53.0, 30.0 / 53.0,
		26.0 / 53.0, 25.0 / 53.0, -23.0 / 53.0
	));
/* */
	munit_assert_not_null(result->input[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( result->input[2], result->input[2]->total,
		4.0, 3.0, 10.0, 0.0,
		-6.0, 7.0, 0.0, -8.0,
		-2.0, 11.0, 9.0, 4.0,
		0.0, 2.0, 0.0, 1.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[2] = matrix_new( 4, 4 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[2], ((matrix_t **)result->answer)[2]->total,
		207.0 / 1270.0, 3.0 / 1270.0, -23.0 / 127.0, 472.0 / 635.0,
		27.0 / 635.0, 28.0 / 635.0, -6.0 / 127.0, 344.0 / 635.0,
		14.0 / 635.0, -9.0 / 635.0, 11.0 / 127.0, -292.0 / 635.0,
		-54.0 / 635.0, -56.0 / 635.0, 12.0 / 127.0, -53.0 / 635.0
	));
/* */
	munit_assert_not_null(result->input[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[3], result->input[3]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, -9.0, 4.0, -4.0
	));
	munit_assert_not_null(((matrix_t **)result->answer)[3] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( ((matrix_t **)result->answer)[3], ((matrix_t **)result->answer)[3]->total,
		3849.0 / 18988.0, -238.0 / 14241.0, -7681.0 / 56964.0, 628.0 / 14241.0, 5149.0 / 56964.0,
		777.0 / 9494.0, 359.0 / 14241.0, -41.0 / 28482.0, -2144.0 / 14241.0, 2549.0 / 28482.0,
		-7.0 / 9494.0, -46.0 / 14241.0, 1711.0 / 28482.0, -836.0 / 14241.0,-1477.0 / 28482.0,
		-1089.0 / 9494.0, -351.0 / 4747.0, 569.0 / 9494.0, 2083.0 / 4747.0, -641.0 / 9494.0,
		-93.0 / 9494.0, 67.0 / 14241.0, -325.0 / 28482.0, 2456.0 / 14241.0, -635.0 / 28482.0
	));
/* */
	munit_assert_not_null(result->input[4] = matrix_new( 4, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[4], result->input[4]->total,
		4.0, 3.0, 10.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 9.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0
	));
	munit_assert_null(((matrix_t **)result->answer)[4] = NULL);
/* */
	munit_assert_not_null(result->input[5] = matrix_new( 5, 5 ));
	munit_assert_not_null(matrix_set_va( result->input[5], result->input[5]->total,
		4.0, 3.0, 0.0, 0.0, 5.0,
		-6.0, 7.0, 0.0, -8.0, 28.0,
		-2.0, 11.0, 0.0, 4.0, 3.0,
		0.0, 2.0, 0.0, 1.0, 5.0,
		-2.0, 10.0, 0.0, 4.0, -4.0
	));
	munit_assert_null(((matrix_t **)result->answer)[5] = NULL);
/* */
	munit_assert_not_null(result->input[6] = matrix_new( 4, 4 ));
	munit_assert_null(((matrix_t **)result->answer)[6] = NULL);
/* */
	munit_assert_not_null(result->input[7] = matrix_idt( 16 ));
	munit_assert_not_null(((matrix_t **)result->answer)[7] = matrix_idt( 16 ));
/* */
	munit_assert_null(result->input[8] = NULL);
	munit_assert_null(((matrix_t **)result->answer)[8] = NULL);

	return result;
}
