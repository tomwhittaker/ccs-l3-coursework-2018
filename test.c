#include "utils.h"
#include <stdlib.h>


int main(int argc, char **argv)
{
    COO I;
    FILE *f;
    CSR sp;
    CSR out;
    if (argc == 2) {
        read_sparse(argv[1], &I);
        
        convert_sparse_to_CSR(I, &sp);
        // print_CSR(sp);

        free_sparse(&I);
        convert_CSR_to_sparse(sp,&I);
        // print_sparse(I);
        printf("Attempt 1\n");
        multiply_CSR_CSR_to_COO(sp,sp);
        printf("Attempt 2\n");
        multiply_CSR_CSR_to_COO2(sp,sp);
        printf("Attempt 3\n");
        multiply_CSR_CSR_to_COO3(sp,sp, &out);

    }
    free_sparseCSR(&sp);
    free_sparse(&I);
    return 0;
}