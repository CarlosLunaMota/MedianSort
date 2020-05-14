/** *********************************************************************** **/
/**                                                                         **/
/**  MedianSort.c                                                           **/
/**  ------------                                                           **/
/**                                                                         **/
/** Content: A generic implementation of the MedianSort algorithm           **/
/**                                                                         **/
/** Author:  Carlos Luna Mota                                               **/
/**                                                                         **/
/** Source:  <https://github.com/CarlosLunaMota/MedianSort>                 **/
/**                                                                         **/
/** License: The Unlicense                                                  **/
/**                                                                         **/
/** This is free and unencumbered software released into the public domain. **/
/**                                                                         **/
/** Anyone is free to copy, modify, publish, use, compile, sell, or         **/
/** distribute this software, either in source code form or as a compiled   **/
/** binary, for any purpose, commercial or non-commercial, and by any       **/
/** means.                                                                  **/
/**                                                                         **/
/** In jurisdictions that recognize copyright laws, the author or authors   **/
/** of this software dedicate any and all copyright interest in the         **/
/** software to the public domain. We make this dedication for the benefit  **/
/** of the public at large and to the detriment of our heirs and            **/
/** successors. We intend this dedication to be an overt act of             **/
/** relinquishment in perpetuity of all present and future rights to this   **/
/** software under copyright law.                                           **/
/**                                                                         **/
/** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,         **/
/** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF      **/
/** MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  **/
/** IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR       **/
/** OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   **/
/** ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR   **/
/** OTHER DEALINGS IN THE SOFTWARE.                                         **/
/**                                                                         **/
/** For more information, please refer to <http://unlicense.org/>           **/
/**                                                                         **/
/** *********************************************************************** **/


/** LIBRARIES ************************************************************* **/

#include <assert.h>     /* Include assertions unless NDEBUG is defined.      */
#include <stdio.h>      /* Input and output functions.                       */
#include <stdlib.h>     /* Memory allocation and random functions.           */
#include <string.h>     /* The size_t type.                                  */
#include <time.h>       /* Time functions.                                   */
#include <math.h>       /* The pow function.                                 */ 


/** GENERIC SORTING FUNCTION TEMPLATES ************************************ **/

#define LESS_THAN(i, j) ((i) < (j))

#define IMPORT_MEDIAN_SORT(type_t, less_than, power, suffix)                    \
                                                                                \
    static void quick_select##suffix(type_t *A, const size_t length,            \
                                     const size_t rank) {                       \
                                                                                \
        size_t l, low  = 0;                                                     \
        size_t h, high = length-1;                                              \
        type_t temp, pivot;                                                     \
                                                                                \
        while (low < high) {                                                    \
            pivot = A[rank];                                                    \
            l     = low;                                                        \
            h     = high;                                                       \
            do {while (less_than(A[l], pivot)) { ++l; }                         \
                while (less_than(pivot, A[h])) { --h; }                         \
                if (l <= h) {                                                   \
                    temp = A[l]; A[l] = A[h]; A[h] = temp;                      \
                    ++l;                                                        \
                    --h;                                                        \
                }                                                               \
            } while (l <= h);                                                   \
            if (h < rank) { low  = l; }                                         \
            if (rank < l) { high = h; }                                         \
        }                                                                       \
    }                                                                           \
                                                                                \
    static void median_sort##suffix(type_t *A, const size_t length) {           \
                                                                                \
        const size_t MIN_SIZE = 1 << power;                                     \
                                                                                \
        size_t l, r, rank, step;                                                \
        type_t t;                                                               \
                                                                                \
        /* MEDIAN SORT (down to 2^(power+1) intervals) */                       \
        for (step   = 1; step <  length;   step <<= 1);                         \
        for (step >>= 1; step >= MIN_SIZE; step >>= 1) {                        \
            for (rank = step; rank < length; rank += (step << 1)) {             \
                l = (rank-step);                                                \
                r = (rank+step) > length ? length : (rank+step);                \
                                                                                \
                /* QUICK SELECT rank in the interval [l, r) */                  \
                quick_select##suffix(A+l, r-l, step);                           \
            }                                                                   \
        }                                                                       \
                                                                                \
        /* INSERTION SORT (up to 2^(power+1) distances) */                      \
        if (MIN_SIZE > 1) {                                                     \
            for (r = 1; r < length; ++r) {                                      \
                t = A[r];                                                       \
                for (l=r; l>0 && less_than(t, A[l-1]); --l) { A[l] = A[l-1]; }  \
                A[l] = t;                                                       \
            }                                                                   \
        }                                                                       \
    }                                                                           \
                                                                                \
                                                                                
#define IMPORT_QUICK_SORT(type_t, less_than, power, suffix)                     \
                                                                                \
    static void quick_sort##suffix(type_t *A, const size_t length) {            \
                                                                                \
        size_t l, r;                                                            \
        type_t t;                                                               \
                                                                                \
        /* QUICK SORT (down to 2^power intervals) */                            \
        if (length > (1 << power)) {                                            \
                                                                                \
            const type_t pivot = A[length/2];                                   \
                                                                                \
            for (l = 0, r = length-1; ; ++l, --r) {                             \
                while (A[l] < pivot) { ++l; }                                   \
                while (A[r] > pivot) { --r; }                                   \
                if (l >= r) { break; }                                          \
                t = A[l]; A[l] = A[r]; A[r] = t;                                \
            }                                                                   \
                                                                                \
            quick_sort##suffix(A,   l);                                         \
            quick_sort##suffix(A+l, length-l);                                  \
        }                                                                       \
                                                                                \
        /* INSERTION SORT (up to 2^power distances) */                          \
        else {                                                                  \
            for (r = 1; r < length; ++r) {                                      \
                t = A[r];                                                       \
                for (l=r; l>0 && less_than(t, A[l-1]); --l) { A[l] = A[l-1]; }  \
                A[l] = t;                                                       \
            }                                                                   \
        }                                                                       \
    }                                                                           \
                                                                                \

#define IMPORT_HEAP_SORT(type_t, less_than, suffix)                             \
                                                                                \
    static void heap_sort##suffix(type_t *A, const size_t length) {             \
                                                                                \
        size_t i, j, k;                                                         \
        type_t temp;                                                            \
                                                                                \
        /* Max-Heapify the array A[0..length-1] in O(length) */                 \
        for (i = (length>>1); i > 0; i--) {                                     \
            j    = i-1;                                                         \
            temp = A[j];                                                        \
            for (;;) {                                                          \
                k = (j<<1) + 1;                                                 \
                if (k+1 >  length) { break; }                                   \
                if (k+2 <= length && less_than(A[k], A[k+1])) { k++; }          \
                if (less_than(A[k], temp)) { break; }                           \
                else { A[j] = A[k]; j = k; }                                    \
            }                                                                   \
            A[j] = temp;                                                        \
        }                                                                       \
                                                                                \
        /* Selection Sort phase in O(length * log(length)) */                   \
        for (i = length-1; i > 0; i--) {                                        \
            temp = A[i];                                                        \
            A[i] = A[0];                                                        \
            j    = 0;                                                           \
            for (;;) {                                                          \
                k = (j<<1) + 1;                                                 \
                if (k+1 >  i) { break; }                                        \
                if (k+2 <= i && less_than(A[k], A[k+1])) { k++; }               \
                if (less_than(A[k], temp)) { break; }                           \
                else { A[j] = A[k]; j = k; }                                    \
            }                                                                   \
            A[j] = temp;                                                        \
        }                                                                       \
    }                                                                           \
                                                                                \

/** AUXILIARY FUNCTIONS *************************************************** **/

int comp_int(const void *i, const void *j) {
    int ii = *(int *) i;
    int jj = *(int *) j;
    return (((ii) > (jj)) - ((ii) < (jj)));
}

int rand_int(const int n) {

        /* Preconditions */
        assert(n > 0);
        assert(RAND_MAX >= n);

        /* Monte-Carlo uniformly random generator */
        int r, range = RAND_MAX - (RAND_MAX % n);
        do { r = rand(); } while (r >= range);
        return r % n;
}

IMPORT_MEDIAN_SORT(int, LESS_THAN, 0, _0)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 1, _1)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 2, _2)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 3, _3)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 4, _4)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 5, _5)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 6, _6)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 7, _7)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 8, _8)
IMPORT_MEDIAN_SORT(int, LESS_THAN, 9, _9)

IMPORT_QUICK_SORT(int, LESS_THAN, 0, _0)
IMPORT_QUICK_SORT(int, LESS_THAN, 1, _1)
IMPORT_QUICK_SORT(int, LESS_THAN, 2, _2)
IMPORT_QUICK_SORT(int, LESS_THAN, 3, _3)
IMPORT_QUICK_SORT(int, LESS_THAN, 4, _4)
IMPORT_QUICK_SORT(int, LESS_THAN, 5, _5)
IMPORT_QUICK_SORT(int, LESS_THAN, 6, _6)
IMPORT_QUICK_SORT(int, LESS_THAN, 7, _7)
IMPORT_QUICK_SORT(int, LESS_THAN, 8, _8)
IMPORT_QUICK_SORT(int, LESS_THAN, 9, _9)

IMPORT_HEAP_SORT(int, LESS_THAN, )

/** MAIN ****************************************************************** **/

int main (void) {

    clock_t crono;
    size_t i, j, step, size;

    const size_t repeat = 10;
    const size_t steps  = 5;

    double T[22][steps];
    size_t S[steps];
    for (step = 1, S[0] = 1000; step < steps; S[step] = S[step-1]*10, step++);

    int *buffer = (int *) malloc(S[steps-1] * sizeof(int));  assert(buffer);
    int *random = (int *) malloc(S[steps-1] * sizeof(int));  assert(random);
    int *sorted = (int *) malloc(S[steps-1] * sizeof(int));  assert(sorted);
    int *array  = (int *) malloc(S[steps-1] * sizeof(int));  assert(array);
    
    for (step = 0, size = 10; step < steps; step++) {

        size = S[step];
        for (i = 0; i < 22; i++) { T[i][step] = 0.0; }

        fprintf(stderr, "\nSORTING %zu RANDOM INTS IN THE RANGE [0,%zu)\n", size, size);

        for (j = 0; j < repeat; j++) {

            /*** GENERATE INSTANCE *******************************************/

            for (i = 0; i < size; i++) { random[i] = rand_int(size); }

            /*** TEST QSORT **************************************************/

            for (i = 0; i < size; i++) { sorted[i] = random[i]; }
            crono = clock();
            qsort(sorted, size, sizeof(int), &comp_int);
            T[20][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 1; i < size; i++) { assert(sorted[i-1] <= sorted[i]); }

            /*** TEST HEAP SORT **********************************************/

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            heap_sort(array, size);
            T[21][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            /*** TEST MEDIAN SORT ********************************************/

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_0(array, size);
            T[0][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_1(array, size);
            T[1][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_2(array, size);
            T[2][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_3(array, size);
            T[3][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_4(array, size);
            T[4][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_5(array, size);
            T[5][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_6(array, size);
            T[6][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_7(array, size);
            T[7][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_8(array, size);
            T[8][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            median_sort_9(array, size);
            T[9][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            /*** TEST QUICK SORT *********************************************/

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_0(array, size);
            T[10][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_1(array, size);
            T[11][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_2(array, size);
            T[12][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_3(array, size);
            T[13][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_4(array, size);
            T[14][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_5(array, size);
            T[15][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_6(array, size);
            T[16][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_7(array, size);
            T[17][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_8(array, size);
            T[18][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            for (i = 0; i < size; i++) { array[i] = random[i]; }
            crono = clock();
            quick_sort_9(array, size);
            T[19][step] += ((double) (clock() - crono)) / CLOCKS_PER_SEC;
            for (i = 0; i < size; i++) { assert(array[i] == sorted[i]); }

            /*****************************************************************/

        }

        fprintf(stderr, "\n    heap_sort   vs qsort = %+.2f %%\n", 100.0 * (T[21][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_0 vs qsort = %+.2f %%\n", 100.0 * (T[ 0][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_0 vs qsort = %+.2f %%\n", 100.0 * (T[10][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_1 vs qsort = %+.2f %%\n", 100.0 * (T[ 1][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_1 vs qsort = %+.2f %%\n", 100.0 * (T[11][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_2 vs qsort = %+.2f %%\n", 100.0 * (T[ 2][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_2 vs qsort = %+.2f %%\n", 100.0 * (T[12][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_3 vs qsort = %+.2f %%\n", 100.0 * (T[ 3][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_3 vs qsort = %+.2f %%\n", 100.0 * (T[13][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_4 vs qsort = %+.2f %%\n", 100.0 * (T[ 4][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_4 vs qsort = %+.2f %%\n", 100.0 * (T[14][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_5 vs qsort = %+.2f %%\n", 100.0 * (T[ 5][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_5 vs qsort = %+.2f %%\n", 100.0 * (T[15][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_6 vs qsort = %+.2f %%\n", 100.0 * (T[ 6][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_6 vs qsort = %+.2f %%\n", 100.0 * (T[16][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_7 vs qsort = %+.2f %%\n", 100.0 * (T[ 7][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_7 vs qsort = %+.2f %%\n", 100.0 * (T[17][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_8 vs qsort = %+.2f %%\n", 100.0 * (T[ 8][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_8 vs qsort = %+.2f %%\n", 100.0 * (T[18][step]-T[20][step]) / T[20][step]);
        fprintf(stderr, "\n  median_sort_9 vs qsort = %+.2f %%\n", 100.0 * (T[ 9][step]-T[20][step]) / T[20][step]);
        fprintf(stderr,   "   quick_sort_9 vs qsort = %+.2f %%\n", 100.0 * (T[19][step]-T[20][step]) / T[20][step]);
    }

    printf("Size qsort HeapSort");
    for (i = 0; i < 10; i++) { printf(" MedianSort(%zu)", i); }
    for (i = 0; i < 10; i++) { printf(" QuickSort(%zu)", i);  }
    for (step = 0; step < steps; step++) {
        printf("\n%zu %.2f %.2f", S[step], 100.0, 100.0 * T[21][step] / T[20][step]);
        for (i =  0; i < 10; i++) { printf(" %.2f", 100.0 * T[i][step] / T[20][step]); }
        for (i = 10; i < 20; i++) { printf(" %.2f", 100.0 * T[i][step] / T[20][step]); }
    }
    
    free(buffer);
    free(random);
    free(sorted);
    free(array);

    return 0;
}
