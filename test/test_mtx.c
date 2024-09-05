#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../include/matrix.h"

int main()
{
	matrix_t *a = matrix_new( 1024, 1024 );
	matrix_t *b = matrix_new( 4, 3 );
	matrix_t *c = NULL;

	struct timespec t1;
	struct timespec t2;

	for ( int i = 0; i < a->total; i++ )
		a->element[i] = -100.0 + 200.0 * (rand() / (double)RAND_MAX);
	// matrix_set_va( a, 25,
	// 	4.0, 3.0, 10.0, 0.0, 5.0,
	// 	-6.0, 7.0, 0.0, -8.0, 28.0,
	// 	-2.0, 11.0, 9.0, 4.0, 3.0,
	// 	0.0, 2.0, 0.0, 1.0, 5.0,
	// 	-2.0, 10.0, -9.0, 4.0, -4.0
	// );
	// matrix_set_va( b, 12,
	// 	4.0, 0.0, 0.0,
	// 	0.0, 2.0, 0.0,
	// 	1.0, 0.0, 3.0,
	// 	0.0, -7.0, 0.0
	// );
	// matrix_print( a, stdout );
	timespec_get(&t1, TIME_UTC);
	//c = matrix_tps( a );
	//printf("%lf\n", matrix_det(a));
	//printf("%ld\n", matrix_rk(a));
	//printf("%lf\n", matrix_tr(a));
	//b = matrix_lup_dec(a, NULL, NULL, NULL);
	//for ( int i = 0; i < 100; i++ )
	//printf("%ld\n", matrix_rk(b));
	c = matrix_mul(a, a);
	timespec_get(&t2, TIME_UTC);
	printf("Timer1: %lf\n", (double)(t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec) * 1e-9);
	// timespec_get(&t1, TIME_UTC);
	// at = matrix_mul(a, a);
	// timespec_get(&t2, TIME_UTC);
	// printf("Timer2: %lf\n", (double)(t2.tv_sec - t1.tv_sec) + (t2.tv_nsec - t1.tv_nsec) * 1e-9);
	//b = matrix_tps(a);
	//matrix_print( c, stdout );
	//matrix_t *c = matrix_idt( 12 );
	//matrix_print( c, stdout );


	return 0;
}
