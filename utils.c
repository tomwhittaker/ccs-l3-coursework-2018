#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include "utils.h"
#include <gmodule.h>


#ifdef _MSC_VER
double drand48()
{
  return (double)rand()/(RAND_MAX + 1);
}
#endif  /* _MSC_VER */

/*
 * Allocate a dense matrix
 * m - number of rows
 * n - number of columns
 * dense - newly allocated matrix.
 */
void alloc_dense(int m, int n, double **dense)
{
  *dense = malloc(m*n*sizeof(**dense));
}

/*
 * Free a dense matrix
 * dense - dense matrix, may be NULL
 */
void free_dense(double **dense)
{
    if (!*dense) {
        return;
    }
    free(*dense);
    *dense = NULL;
}

/*
 * Zero a dense matrix
 * m - number of rows
 * n - number of columns
 * dense - matrix to zero.
 */
void zero_dense(int m, int n, double *dense)
{
    int i, j;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            dense[j*m + i] = 0;
        }
    }
}

/*
 * Allocate a sparse matrix in coordinate format.
 * m - number of rows
 * n - number of columns
 * NZ - number of nonzeros
 * sparse - newly allocated matrix.
 */
void alloc_sparse(int m, int n, int NZ, COO *sparse)
{
    COO sp = calloc(1, sizeof(struct _p_COO));
    sp->m = m;
    sp->n = n;
    sp->NZ = NZ;
    sp->coords = calloc(NZ, sizeof(struct coord));
    sp->data = calloc(NZ, sizeof(double));
    *sparse = sp;
}

void get_values_and_coords(const CSR sparse, struct coord *coords, double *data){
    int NZ;
    NZ = sparse->NZ;
    int n;
    int counter;
    int row = 0;
    counter = 0;

    for (n=0; n<NZ;n++){
        while (counter >= sparse->IA[row+1]){
            row++;
        }
        data[n] = sparse->A[n];
        coords[n].i = sparse->JA[n];
        coords[n].j = row;
        counter++;
    }
}

void alloc_sparseCSR(int m, int n, int NZ, CSR *sparse)
{
    CSR sp = calloc(1, sizeof(struct _p_CSR));
    sp->m = m;
    sp->n = n;
    sp->NZ = NZ;
    sp->A = calloc(NZ, sizeof(double));
    sp->IA = calloc(m+1, sizeof(int));
    sp->JA = calloc(NZ, sizeof(int));
    *sparse = sp;
}

//TODO Combine all into single function so fewer for loops
void get_IA_from_sorted_COO(const COO sparse, int *IA)
{
    int size;
    size = sparse->NZ;
    size = size + 1;
    int n;
    int acc;
    IA[0]=0;
    for (n = 0; n < sparse->NZ; n++) {
        int j = sparse->coords[n].j;
        IA[j+1] = IA[j+1] + 1;
    }
    for (n = 1; n < size; n++) {
        IA[n] = IA[n] + IA[n-1];
    }
}

void get_A_from_sorted_COO(const COO sparse, double *A)
{
    int n;
    int i, j;
    int counter;
    counter = 0;
    for (n = 0; n < sparse->NZ; n++) {
        double data = sparse->data[n];
        A[counter] = data;
        counter = counter + 1;
    }
}

void get_JA_from_sorted_C00(const COO sparse, int *JA){
    int n = 0;
    for (n=0; n<sparse->NZ; n++){
        JA[n] = sparse->coords[n].i;
    }
}

/*
 * Free a sparse matrix.
 * sparse - sparse matrix, may be NULL
 */
void free_sparse(COO *sparse)
{
    COO sp = *sparse;
    if (!sp) {
        return;
    }
    free(sp->coords);
    free(sp->data);
    free(sp);
    *sparse = NULL;
}

void free_sparseCSR(CSR *sparse)
{
    CSR sp = *sparse;
    if (!sp) {
        return;
    }
    free(sp->A);
    free(sp->IA);
    free(sp->JA);
    free(sp);
    *sparse = NULL;
}

// void sort_COO(const COO sparse, COO *sparseSorted){
//     COO sp;
//     alloc_sparse(sparse->m, sparse->n, sparse->NZ, &sp)
//     int i;
//     for (i=0;i<sp->NZ;i++){
//         sp->coords.i = sparse->coords.j
//         sp->coords.j = sparse->coords.i
//     }

// }

void multiply_CSRVector(const CSR A, const double * B, double *result){
    int counter;
    int check;
    counter = 1;
    check = A->IA[counter];
    for (int x = 0; x<3;x++){
        result[x]=0;
    }
    for (int x = 0; x<A->NZ; x++){
        if (x == check){
            counter++;
            check = A->IA[counter];
        }
        result[counter-1] = result[counter-1] + A->A[x] * B[A->JA[x]];
    }

}

void multiply_CSR_CSR_to_COO(const CSR A, const CSR B){
    // COO sp;
    // alloc_sparse(A->m, B->n, A->m*B->n , &sp);
    

    // for (int i = 0; i< A->m*B->n ; i++){
    //     sp->values[i] = 0
    // }
    double *dense;
    alloc_dense(A->m,B->n,&dense);
    zero_dense(A->m,B->n,dense);
    double *row;

    int counter = 0;
    int counter2 = 0;
    int rA = 0;
    int numberInRowA;
    for (int cA = 0; cA<A->NZ; cA++){
        // if (cA >= A->IA[rA+1]){

        // }
        //Go through add to correct place in AI and add to array of array for both IA and JA and dealloc mem here so free for other threads

        //Maybe chance outer loop to go through IA and then do a for loop for each in that row. 
        // Means easy parallel also asymtopically the same??? as O(numRows * lengthRows) == O(length A)
        while (cA >= A->IA[rA+1]){
            rA++;
        }
        int colA = A->JA[cA];
        int rB = 0;
        for (int cB = 0; cB<B->NZ; cB++){
            while (cB >= B->IA[rB+1]){
                rB++;
            }
            int colB = B->JA[cB];
            counter2++;
            if (colA ==  rB){
                counter ++;
                dense[rA*A->m + colB] += A->A[cA] * B->A[cB];
            }
        }
        
    }
    for (int x=0; x<A->m; x++){
        for (int y=0; y<B->n; y++){
            printf("%.2f ", dense[x*A->m+y]);
        }
        printf("\n");
    }
    int expect = A->m * A->m * B->n;
    printf("Expected: %d -> Improved: %d\n", expect,counter );
    printf("Hit: %d \n", counter2 );

}

void multiply_CSR_CSR_to_COO2(const CSR A, const CSR B){
    // COO sp;
    // alloc_sparse(A->m, B->n, A->m*B->n , &sp);
    

    // for (int i = 0; i< A->m*B->n ; i++){
    //     sp->values[i] = 0
    // }
    double *dense;
    alloc_dense(A->m,B->n,&dense);
    zero_dense(A->m,B->n,dense);
    double *row;

    int counter = 0;
    int counter2 = 0;
    int rA = 0;
    int numberInRowA;
    for (int rA = 0; rA<A->m; rA++){
        
        // if (cA >= A->IA[rA+1]){
        for (int cA = A->IA[rA]; cA<A->IA[rA+1]; cA++){
            
            int colA = A->JA[cA];
            for (int cB = B->IA[colA]; cB<B->IA[colA+1]; cB++){
                counter2++;
                int colB = B->JA[cB];
                counter ++;
                dense[rA*A->m + colB] += A->A[cA] * B->A[cB];
                //Maybe need to add if statement back 
            }
        }
        // }/for (int cA = 0; cA<A->NZ; cA++){
        //Go through add to correct place in AI and add to array of array for both IA and JA and dealloc mem here so free for other threads

        //Maybe chance outer loop to go through IA and then do a for loop for each in that row. 
        // Means easy parallel also asymtopically the same??? as O(numRows * lengthRows) == O(length A)
    }
    for (int x=0; x<A->m; x++){
        for (int y=0; y<B->n; y++){
            printf("%.2f ", dense[x*A->m+y]);
        }
        printf("\n");
    }
    int expect = A->m * A->m * B->n;
    printf("Expected: %d -> Improved: %d\n", expect,counter );
    printf("Hit: %d \n", counter2 );

}

void multiply_CSR_CSR_to_COO3(const CSR A, const CSR B){
    // COO sp;
    // alloc_sparse(A->m, B->n, A->m*B->n , &sp);
    

    // for (int i = 0; i< A->m*B->n ; i++){
    //     sp->values[i] = 0
    // }
    int counter = 0;
    int counter2 = 0;
    int rA = 0;
    int numberInRowA;
    int IARes[A->m+1];

    IARes[0]=0;
    for (int rA = 0; rA<A->m; rA++){
        IARes[rA+1]=IARes[rA];
        // if (cA >= A->IA[rA+1]){
        for (int cA = A->IA[rA]; cA<A->IA[rA+1]; cA++){
            
            int colA = A->JA[cA];
            for (int cB = B->IA[colA]; cB<B->IA[colA+1]; cB++){
                counter2++;
                int colB = B->JA[cB];
                counter ++;
                IARes[rA+1] += 1;
                //Maybe need to add if statement back 
            }
        }
        // }/for (int cA = 0; cA<A->NZ; cA++){
        //Go through add to correct place in AI and add to array of array for both IA and JA and dealloc mem here so free for other threads

        //Maybe chance outer loop to go through IA and then do a for loop for each in that row. 
        // Means easy parallel also asymtopically the same??? as O(numRows * lengthRows) == O(length A)
    }
    int JARes[IARes[A->m]];
    int ARes[IARes[A->m]];
    for (int rA = 0; rA<A->m; rA++){
        int row[B->n];
        for (int i=0; i<B->n;i++){
            row[i] = 0;
        }
        // if (cA >= A->IA[rA+1]){
        for (int cA = A->IA[rA]; cA<A->IA[rA+1]; cA++){
            int colA = A->JA[cA];
            for (int cB = B->IA[colA]; cB<B->IA[colA+1]; cB++){
                counter2++;
                int colB = B->JA[cB];
                counter ++;
                row[colB] += A->A[cA] * B->A[cB];
                //Maybe need to add if statement back 
            }
        }
        int count = 0;
        for (int i=0; i<B->n;i++){
            if (row[i] !=0 ){
                JARes[IARes[rA]+count]=i;
                ARes[IARes[rA]+count]=row[i];
                count++;
            }
        }
    }
    printf("A:\n");
    for (int i = 0; i<IARes[A->m];i++){
        printf("%d ", ARes[i]);
    }
    printf("\n");
    printf("IA:\n");
    for (int i = 0; i<A->m+1;i++){
        printf("%d ", IARes[i]);
    }
    printf("\n");
    printf("JA:\n");
    for (int i = 0; i<IARes[A->m];i++){
        printf("%d ", JARes[i]);
    }
    printf("\n");
    int expect = A->m * A->m * B->n;
    printf("Expected: %d -> Improved: %d\n", expect,counter );
    printf("Hit: %d \n", counter2 );

}




/*
 * Convert a sparse matrix to dense format in column major format.
 *
 * sparse - The sparse matrix to convert
 * dense - pointer to output dense matrix (will be allocated)
 */
void convert_sparse_to_dense(const COO sparse, double **dense)
{
    int n;
    int i, j;
    alloc_dense(sparse->m, sparse->n, dense);
    zero_dense(sparse->m, sparse->n, *dense);
    for (n = 0; n < sparse->NZ; n++) {
        i = sparse->coords[n].i;
        j = sparse->coords[n].j;
        (*dense)[j * sparse->m + i] = sparse->data[n];
    }
}

void convert_sparse_to_CSR(const COO I, CSR *sparse)
{
    CSR sp;
    alloc_sparseCSR(I->m, I->n, I->NZ, &sp);
    get_A_from_sorted_COO(I, sp->A);
    get_IA_from_sorted_COO(I, sp->IA);
    get_JA_from_sorted_C00(I, sp->JA);
    *sparse = sp;
}

void convert_CSR_to_sparse(const CSR I, COO *sparse)
{
    COO sp;
    alloc_sparse(I->m, I->n, I->NZ, &sp);
    get_values_and_coords(I,sp->coords,sp->data);
    *sparse = sp;
}

void print_CSR(const CSR I){
    int NZ;
    NZ = I->NZ;
    double *A;
    int *IA;
    int *JA;
    A = I->A;
    IA = I->IA;
    JA = I->JA;

    printf("A\n");
    int n;
    for (n=0; n< NZ; n++){
        printf("%.2f,", A[n]);
    }
    printf("\n");

    printf("IA\n");
    for (n=0; n<= NZ; n++){
        printf("%d,", IA[n]);
    }
    printf("\n");

    printf("JA\n");
    for (n=0; n< NZ; n++){
        printf("%d,", JA[n]);
    }
    printf("\n");
}
/*
 * Convert a dense matrix in column major format to sparse.
 * Entries with absolute value < 1e-15 are flushed to zero and not
 * stored in the sparse format.
 *
 * dense - the dense array
 * m - number of rows
 * n - number of columns
 * sparse - output sparse matrix (allocated by this routine)
 */
void convert_dense_to_sparse(const double *dense, int m, int n,
                             COO *sparse)
{
    int i, j, NZ;
    COO sp;
    NZ = 0;
    /* Figure out how many nonzeros we're going to have. */
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double val = dense[j*m + i];
            if (fabs(val) > 1e-15) {
                NZ++;
            }
        }
    }
    alloc_sparse(m, n, NZ, &sp);

    NZ = 0;
    /* Fill up the sparse matrix */
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            double val = dense[j*m + i];
            if (fabs(val) > 1e-15) {
                sp->coords[NZ].i = i;
                sp->coords[NZ].j = j;
                sp->data[NZ] = val;
                NZ++;
            }
        }
    }
    *sparse = sp;
}

/*
 * Create a random sparse matrix
 *
 * m - number of rows
 * n - number of columns
 * frac - fraction of entries that should be nonzero
 * sparse - newly allocated random matrix.
 */
void random_matrix(int m, int n, double frac, COO *sparse)
{
    int i, j;
    double *d;
    alloc_dense(m, n, &d);
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            if (drand48() < frac) {
                d[j*m + i] = drand48();
            } else {
              d[j*m + i] = 0.0;
            }
        }
    }
    convert_dense_to_sparse(d, m, n, sparse);
    free_dense(&d);
}

/*
 * Read a sparse matrix from a file.
 *
 * file - The filename to read
 * sparse - The newly read sparse matrix (allocated here)
 */
void read_sparse(const char *file, COO *sparse)
{
    COO sp;
    int i, j, k, m, n, NZ;
    double val;
    int c;
    FILE *f = fopen(file, "r");
    if (!f) {
        fprintf(stderr, "Unable to open %s for reading.\n", file);
        exit(1);
    }
    c = fscanf(f, "%d %d %d\n", &m, &n, &NZ);
    if (c != 3) {
        fprintf(stderr, "File format incorrect on line 1, expecting 3 integers, got %d\n", c);
        fclose(f);
        exit(1);
    }
    if (NZ > (uint64_t)m*n) {
        fprintf(stderr, "More nonzeros (%d) than matrix entries (%d x %d)!\n", NZ, m, n);
        fclose(f);
        exit(1);
    }
    alloc_sparse(m, n, NZ, &sp);
    k = 0;
    while ((c = fscanf(f, "%d %d %lg\n", &i, &j, &val)) == 3) {
        if (k >= NZ) {
            fprintf(stderr, "File has nonzero lines than expected (%d)\n", NZ);
            fclose(f);
            free_sparse(&sp);
            exit(1);
        }
        if (i >= m || j >= n) {
            fprintf(stderr, "Entry on line %d incorrect, index (%d, %d) out of bounds for %d x %d matrix\n", k + 2, i, j, m, n);
            fclose(f);
            free_sparse(&sp);
            exit(1);
        }
        sp->coords[k].i = i;
        sp->coords[k].j = j;
        sp->data[k] = val;
        k++;
    }

    if (k != NZ) {
        fprintf(stderr, "File has fewer lines (%d) than expected (%d)\n",
                k, NZ);
        fclose(f);
        free_sparse(&sp);
        exit(1);
    }
    *sparse = sp;
    fclose(f);
}

/*
 * Write a sparse matrix to a file.
 *
 * f - The file handle.
 * sp - The sparse matrix to write.
 */
void write_sparse(FILE *f, COO sp)
{
    int i;
    fprintf(f, "%d %d %d\n", sp->m, sp->n, sp->NZ);
    for (i = 0; i < sp->NZ; i++) {
        fprintf(f, "%d %d %g\n", sp->coords[i].i, sp->coords[i].j, sp->data[i]);
    }
}

/*
 * Print a sparse matrix to stdout
 *
 * sp - The sparse matrix to print.
 */
void print_sparse(COO sp)
{
    write_sparse(stdout, sp);
}

void read_sparse_binary(const char *file, COO *sparse)
{
    COO sp;
    int m, n, NZ;
    size_t nread;
    FILE *f = fopen(file, "r");

    nread = fread(&m, sizeof(m), 1, f);
    if (nread != 1) {
      fprintf(stderr, "Did not read rows from file\n");
      exit(1);
    }
    nread = fread(&n, sizeof(n), 1, f);
    if (nread != 1) {
      fprintf(stderr, "Did not read columns from file\n");
      exit(1);
    }
    nread = fread(&NZ, sizeof(NZ), 1, f);
    if (nread != 1) {
      fprintf(stderr, "Did not read number of nonzeros from file\n");
      exit(1);
    }
    alloc_sparse(m, n, NZ, &sp);
    nread = fread(sp->coords, sizeof(*sp->coords), NZ, f);
    if (nread != NZ) {
      fprintf(stderr, "Did not read nonzero locations from file\n");
      exit(1);
    }
    nread = fread(sp->data, sizeof(*sp->data), NZ, f);
    if (nread != NZ) {
      fprintf(stderr, "Did not read nonzero values from file\n");
      exit(1);
    }
    *sparse = sp;
    fclose(f);
}

void write_sparse_binary(FILE *f, COO sp)
{
  size_t nwrite;
  nwrite = fwrite(&(sp->m), sizeof(sp->m), 1, f);
  if (nwrite != 1) {
    fprintf(stderr, "Could not write rows to output file\n");
    exit(1);
  }

  nwrite = fwrite(&(sp->n), sizeof(sp->n), 1, f);
  if (nwrite != 1) {
    fprintf(stderr, "Could not write columns to output file\n");
    exit(1);
  }

  nwrite = fwrite(&(sp->NZ), sizeof(sp->NZ), 1, f);
  if (nwrite != 1) {
    fprintf(stderr, "Could not write number of nonzeros to output file\n");
    exit(1);
  }
  nwrite = fwrite(sp->coords, sizeof(*sp->coords), sp->NZ, f);
  if (nwrite != sp->NZ) {
    fprintf(stderr, "Could not write nonzero locations to output file\n");
    exit(1);
  }
  nwrite = fwrite(sp->data, sizeof(*sp->data), sp->NZ, f);
  if (nwrite != sp->NZ) {
    fprintf(stderr, "Could not write nonzero values to output file\n");
    exit(1);
  }
}
