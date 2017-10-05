/**
 * Calculate pi value.
 *
 * An exercise to print the pi value upto
 * 1 million decimal places using nonic iterations.
 *
 *
 * Author: Abhineet Gupta
 *      
 */
#include <iostream>
#include <iomanip>
#include <mpir.h>
#include <stdlib.h>
#include <string.h>


using namespace std;

const int MAX_ITERATIONS = 100;
const int PLACES = 1000000;        // desired decimal places
const int PRECISION = PLACES + 1;  // +1 for the digit 3 before the decimal

const int BASE = 10;  // base 10 numbers
const int BIT_COUNT = 8;   // bits per machine word

const int BLOCK_SIZE = 10;                // print digits in blocks
const int LINE_SIZE = 100;               // digits to print per line
const int LINE_COUNT = PLACES / LINE_SIZE;  // lines to print
const int GROUP_SIZE = 5;                 // line grouping size

/**
 * Instead of using newton's method now halley's method is used to
 * compute cuber root as done by assignment solution provided by Ron_Mak
 * @param x where to store the result.
 * @param a the number whose cube root to compute.
 */
void cube_root(mpf_t& x, const mpf_t a);

/**
 * This method converts the mpir floating variable pi to
 * char pointer variable containing the whole pi values
 * Using this pointer , expected output is printed.
 *
 */
void print_big_pi(mpf_t& pi);

/**
 * The main.
 * Applies nonic iterations to calculate pi value.
 */
int main() {



	bool first_time = true;
	int equalityFlag = 1;

	mpf_set_default_prec(BIT_COUNT * PRECISION);  // precision in bits

	mpf_t pi, cons1, cons2, cons3, cons4, cons5, cons6, a, r, s, t, u, v, w,
			exponential, interMedPi, localVarOne, localVarTwo;

	mpf_init_set_str(pi, "0", BASE);
	mpf_init_set_str(cons1, "1", BASE);
	mpf_init_set_str(cons2, "2", BASE);
	mpf_init_set_str(cons3, "3", BASE);
	mpf_init_set_str(cons4, "0.5", BASE);
	mpf_init_set_str(cons5, "9", BASE);
	mpf_init_set_str(cons6, "27", BASE);
	mpf_init_set_str(a, "1", BASE);
	mpf_init_set_str(r, "3", BASE);
	mpf_init_set_str(s, "5", BASE);
	mpf_init_set_str(t, "0", BASE);
	mpf_init_set_str(u, "0", BASE);
	mpf_init_set_str(v, "0", BASE);
	mpf_init_set_str(w, "0", BASE);
	mpf_init_set_str(exponential, "0", BASE);
	mpf_init_set_str(interMedPi, "0", BASE);
	mpf_init_set_str(localVarOne, "0", BASE);
	mpf_init_set_str(localVarTwo, "0", BASE);

	mpf_div(a, a, cons3);

	mpf_sqrt(localVarOne, r);
	mpf_sub(localVarOne, localVarOne, cons1);
	mpf_mul(r, cons4, localVarOne);

	mpf_pow_ui(localVarOne, r, mpf_get_ui(cons3));
	mpf_sub(localVarOne, cons1, localVarOne);
	cube_root(s, localVarOne);
	for (int i = 0; i < MAX_ITERATIONS; i++) {

		mpf_mul(localVarOne, cons2, r);
		mpf_add(t, cons1, localVarOne);

		mpf_mul(localVarOne, r, r);
		mpf_add(localVarOne, localVarOne, r);
		mpf_add(u, cons1, localVarOne);

		mpf_mul(localVarOne, u, cons5);
		mpf_mul(localVarOne, localVarOne, r);

		cube_root(u, localVarOne);

		mpf_add(localVarOne, u, t);
		mpf_mul(localVarOne, u, localVarOne);
		mpf_mul(localVarTwo, t, t);
		mpf_add(v, localVarTwo, localVarOne);

		mpf_mul(localVarOne, s, s);
		mpf_add(localVarOne, s, localVarOne);
		mpf_add(localVarOne, cons1, localVarOne);
		mpf_mul(w, cons6, localVarOne);

		mpf_div(w, w, v);

		if (first_time == true) {

			mpf_set(exponential, a);
			first_time = false;

		} else {

			mpf_mul(exponential, exponential, cons5);
		}

		mpf_sub(localVarOne, cons1, w);
		mpf_mul(localVarOne, exponential, localVarOne);
		mpf_mul(localVarTwo, w, a);
		mpf_add(interMedPi, localVarTwo, localVarOne);

		mpf_sub(s, cons1, r);

		mpf_pow_ui(s, s, mpf_get_ui(cons3));

		mpf_mul(localVarOne, cons2, u);
		mpf_add(localVarOne, t, localVarOne);
		mpf_mul(localVarOne, localVarOne, v);
		mpf_div(s, s, localVarOne);

		mpf_pow_ui(localVarOne, s, mpf_get_ui(cons3));
		mpf_sub(localVarOne, cons1, localVarOne);

		cube_root(r, localVarOne);

		equalityFlag = mpf_cmp(interMedPi, a);

		if (i > 0 && (equalityFlag == 0))
			break;

		mpf_set(a, interMedPi);

	}
	mpf_div(pi, cons1, interMedPi);
	print_big_pi(pi);

	// clearing big size variables
	mpf_clears(pi, a, interMedPi);

	return 0;
}


 // Method changed from newton's method to halley's method
void cube_root(mpf_t& x, const mpf_t a) {

	 mpf_t x_prev; mpf_init(x_prev);

	 mpf_t temp1;  mpf_init(temp1);
	  mpf_t temp2; mpf_init(temp2);
	  mpf_t two_a;  mpf_init(two_a);
	   mpf_t x_cubed; mpf_init(x_cubed);


	   mpf_t three; mpf_init(three); mpf_set_str(three, "3", BASE);

	     // Set an initial estimate for x.
	     mpf_div(x, a, three);  // x = a/3

	     int n = 0; // iteration counter
	   // Loop until two consecutive values are equal
	     // or up to MAX_ITERATIONS times.
	     do
	     {
	       mpf_set(x_prev, x);

	       mpf_mul(x_cubed, x, x);
	       mpf_mul(x_cubed, x_cubed, x);    // x_cubed = x^3
	       mpf_add(two_a, a, a);        // two_a = 2a
	       mpf_add(temp1, x_cubed, two_a);   // temp1 = x^3 + 2a
	       mpf_add(temp2, x_cubed, x_cubed);  // temp2 = 2x^3
	       mpf_add(temp2, temp2, a);      // temp2 = 2x^3 + a
	       mpf_div(temp1, temp1, temp2);    // temp1 = (x^3 + 2a)/(2x^3 + a)
	       mpf_mul(x, x, temp1);        // x = x((x^3 + 2a)/(2x^3 + a))

	       n++;
	     } while ((mpf_cmp(x, x_prev) != 0) && (n < MAX_ITERATIONS));
	   }

/*  Commenting out newton's implementation to calculate cube root
 * int compareFlag = 1;
	mpf_set_default_prec(BIT_COUNT * PRECISION);  // precision in bits
	mpf_t interCubeVal, cubeRootVal, localPrecision, localVarOne, localVarTwo,
			localVarThree, cons2;
	mpf_init_set_str(interCubeVal, "0", BASE);
	mpf_init_set_str(cubeRootVal, "0", BASE);
	mpf_init_set_str(localPrecision,
			"0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001",
			BASE);
	mpf_init_set_str(localVarOne, "0", BASE);
	mpf_init_set_str(localVarTwo, "0", BASE);
	mpf_init_set_str(localVarThree, "0", BASE);
	mpf_init_set_str(cons2, "2", BASE);

	mpf_set(interCubeVal, a);
	mpf_set(cubeRootVal, a);

	while (compareFlag >= 0) {

		mpf_set(interCubeVal, cubeRootVal);

		mpf_mul(localVarTwo, interCubeVal, interCubeVal);
		mpf_div(localVarOne, a, interCubeVal);
		mpf_sub(localVarOne, localVarTwo, localVarOne);

		mpf_mul(localVarThree, cons2, interCubeVal);
		mpf_div(localVarTwo, a, localVarTwo);
		mpf_add(localVarThree, localVarThree, localVarTwo);

		mpf_div(localVarOne, localVarOne, localVarThree);
		mpf_sub(cubeRootVal, interCubeVal, localVarOne);

		mpf_sub(localVarOne, cubeRootVal, interCubeVal);
		mpf_abs(localVarTwo, localVarOne);
		mpf_div(localVarThree, localVarTwo, interCubeVal);
		compareFlag = mpf_cmp(localVarThree, localPrecision);

	}

	mpf_set(q, cubeRootVal);*/


void print_big_pi(mpf_t& pi) {

	char *str2 = 0;
	long int exponent = 01;
	str2 = mpf_get_str(NULL, &exponent, BASE, 0, pi);


	cout << *str2 << "."; // string returned removed . operator so need to append that
	str2++;

	for (int i = 0; i < LINE_COUNT; i++) {

		if (i == GROUP_SIZE) {
			cout << endl;
		}
		for (int j = 0; j < BLOCK_SIZE; j++) {
			int count = 0;
			while (count != 10) {
				if (*str2 == '\0') {
					break;
				}
				cout << *str2;
				str2++;
				count++;
			}
			cout << " ";
		}
		cout << setw(3) << endl;
	}

	// clearing up big pi string
	delete str2;

}

