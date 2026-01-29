#include "test_Constraints.h"
#include "test_CoreTransformations.h"
#include "test_Matrix.h"
//  #include "test_Parallel_cols.h"
#include "test_Parallel_rows.h"
#include "test_Presolver.h"
#include "test_SimpleReductions.h"
#include "test_domain_propagation.h"
#include "test_dton.h"
#include "test_iVec.h"
#include "test_pathological.h"
//  #include "test_postsolve.h"
#include "test_ston.h"

const char *run_all_tests()
{
    mu_assert("matrix error", test_matrix());
    mu_assert("constraints error", test_constraints());
    mu_assert("iVec error", test_iVec());
    mu_assert("dton error", test_dton());
    mu_assert("core error", test_core());
    mu_assert("ston error", test_ston());
    mu_assert("parallel_rows error", test_parallel_rows());
    mu_assert("simple reductions error", test_simple());
    mu_assert("domain propagation error", test_domain());
    //     mu_assert("parallel_cols error", test_parallel_cols());
    //     mu_assert("postsolve error", test_postsolve());
    mu_assert("presolver error", test_presolver());
    mu_assert("pathological error", test_pathological());

    return NULL;
}

int main()
{
    const char *result = run_all_tests();
    if (result != NULL)
    {
        printf("Test failed: %s\n", result);
        return -1;
    }
    printf("All tests passed!\n");
    return 0;
}
