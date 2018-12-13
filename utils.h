#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>

struct coord {
    int i, j;
};

struct _p_COO {
    int m, n, NZ;
    struct coord *coords;
    double *data;
};

struct _p_CSR {
	int m,n,NZ;
	double *A;
	int *IA;
	int *JA;
};

typedef struct _p_COO *COO;
typedef struct _p_CSR *CSR;


void alloc_sparse(int, int, int, COO*);
void free_sparse(COO*);
void alloc_dense(int, int, double **);
void free_dense(double **);
void zero_dense(int, int, double *);

void alloc_sparseCSR(int,int,int, CSR*);
void free_sparseCSR(CSR*);

void convert_sparse_to_dense(const COO, double **);
void convert_dense_to_sparse(const double *, int, int, COO *);
void convert_CSR_to_sparse(const CSR, COO *);
void convert_sparse_to_CSR(const COO, CSR *);

void read_sparse(const char *, COO *);
void write_sparse(FILE *, COO);
void read_sparse_binary(const char *, COO *);
void write_sparse_binary(FILE *, COO);
void print_sparse(COO);
void random_matrix(int, int, double, COO *);
void print_CSR(const CSR);

void get_A_from_sorted_COO(const COO, double *);
void get_IA_from_sorted_COO(const COO, int *);
void get_JA_from_sorted_C00(const COO, int *);
void get_values_and_coords(const CSR, struct coord *, double *);

void multiply_CSRVector(const CSR, const double * , double *);
void multiply_CSR_CSR_to_COO2(const CSR, const CSR);
void multiply_CSR_CSR_to_COO(const CSR, const CSR);
void multiply_CSR_CSR_to_COO3(const CSR A, const CSR B, CSR *);

#endif
