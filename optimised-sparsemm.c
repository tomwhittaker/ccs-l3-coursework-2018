#include "utils.h"

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);

/* Computes C = A*B.
 * C should be allocated by this routine.
 */
void optimised_sparsemm(const COO A, const COO B, COO *C)
{
	CSR sp1;
	CSR sp2;
	CSR out;
	convert_sparse_to_CSR(A, &sp1);
	convert_sparse_to_CSR(B, &sp2);
	// printf("A\n");
	// for (int i = 0; i<sp1->NZ;i++){
 //        printf("%.2f ", sp1->A[i]);
 //    }
 //    printf("\n");
 //    printf("IA:\n");
 //    for (int i = 0; i<sp1->m+1;i++){
 //        printf("%d ", sp1->IA[i]);
 //    }
 //    printf("\n");
 //    printf("JA:\n");
 //    for (int i = 0; i<sp1->NZ;i++){
 //        printf("%d ", sp1->JA[i]);
 //    }
 //    printf("\n");
 //    printf("A\n");
 //    for (int i = 0; i<sp2->NZ;i++){
 //        printf("%.2f ", sp2->A[i]);
 //    }
 //    printf("\n");
 //    printf("IA:\n");
 //    for (int i = 0; i<sp2->m+1;i++){
 //        printf("%d ", sp2->IA[i]);
 //    }
 //    printf("\n");
 //    printf("JA:\n");
 //    for (int i = 0; i<sp2->NZ;i++){
 //        printf("%d ", sp2->JA[i]);
 //    }
 //    printf("\n");


	// printf("Attempt 1\n");
 //    multiply_CSR_CSR_to_COO(sp1,sp2);
 //    printf("Attempt 2\n");
 //    multiply_CSR_CSR_to_COO2(sp1,sp2);
 //    printf("Attempt 3\n");
    multiply_CSR_CSR_to_COO3(sp1,sp2, &out);
    convert_CSR_to_sparse(out, C);
    free_sparseCSR(&sp1);
    free_sparseCSR(&sp2);
    free_sparseCSR(&out);
    // return basic_sparsemm(A, B, C);
}

/* Computes O = (A + B + C) (D + E + F).
 * O should be allocated by this routine.
 */
void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
                            const COO D, const COO E, const COO F,
                            COO *O)
{
    return basic_sparsemm_sum(A, B, C, D, E, F, O);
}
