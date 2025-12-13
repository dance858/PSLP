#include "PSLP_API.h"
#include "PSLP_stats.h"
#include <cassert>
#include <iostream>
#include <limits>

// run affiro presolver test in C++
int main()
{
    double INF = std::numeric_limits<double>::infinity();

    double Ax[] = {
        -1.,   1.,    1.,    -1.06, 1.,    1.,    -1.,   1.4,   -1.,   -1.,   -1.,
        -1.,   1.,    1.,    -1.06, -1.06, -0.96, -0.86, 1.,    1.,    -1.,   1.,
        -1.,   1.,    -1.,   1.,    -1.,   -1.,   1.,    1.,    1.,    -0.43, 1.,
        1.,    -1.,   1.4,   -0.43, -0.43, -0.39, -0.37, 1.,    1.,    1.,    1.,
        1.,    -1.,   1.,    1.,    1.,    -1.,   1.,    -1.,   1.,    -1.,   1.,
        -1.,   2.364, 2.386, 2.408, 2.429, -1.,   2.191, 2.219, 2.249, 2.279, -1.,
        0.109, -1.,   0.109, 0.108, 0.108, 0.107, 0.301, -1.,   0.301, 0.313, 0.313,
        0.326, -1.,   1.,    1.,    1.,    1.};
    int Ai[] = {0,  1,  2,  0,  3,  0,  1,  12, 4,  5,  6,  7,  12, 13, 4,  5,  6,
                7,  14, 4,  8,  5,  9,  6,  10, 7,  11, 15, 16, 17, 18, 15, 19, 15,
                16, 28, 20, 21, 22, 23, 30, 20, 21, 22, 23, 28, 29, 31, 20, 24, 21,
                25, 22, 26, 23, 27, 8,  9,  10, 11, 18, 24, 25, 26, 27, 2,  15, 13,
                20, 21, 22, 23, 0,  17, 4,  5,  6,  7,  29, 3,  19, 14, 30};
    int Ap[] = {0,  3,  5,  6,  8,  14, 19, 21, 23, 25, 27, 31, 33, 34,
                36, 41, 48, 50, 52, 54, 56, 65, 67, 72, 74, 79, 81, 83};
    int nnz = 83;
    int n_rows = 27;
    int n_cols = 32;
    double lhs[] = {0.,   0.,   -INF, -INF, 0.,   0.,   -INF, -INF, -INF,
                    -INF, 0.,   0.,   -INF, -INF, 0.,   44.,  -INF, -INF,
                    -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF, -INF};
    double rhs[] = {0., 0.,  80.,  0., 0., 0., 80., 0., 0., 0., 0., 0.,   500., 0.,
                    0., 44., 500., 0., 0., 0., 0.,  0., 0., 0., 0., 310., 300.};
    double lbs[] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    double ubs[] = {INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                    INF, INF, INF, INF, INF, INF, INF, INF, INF, INF, INF,
                    INF, INF, INF, INF, INF, INF, INF, INF, INF, INF};
    double c[] = {0., -0.4,  0., 0., 0., 0.,   0.,    0., 0., 0., 0.,
                  0., -0.32, 0., 0., 0., -0.6, 0.,    0., 0., 0., 0.,
                  0., 0.,    0., 0., 0., 0.,   -0.48, 0., 0., 10.};

    Settings *stgs = default_settings();
    Presolver *presolver =
        new_presolver(Ax, Ai, Ap, n_rows, n_cols, nnz, lhs, rhs, lbs, ubs, c, stgs);
    assert(presolver != nullptr && "Presolver initialization failed");

    PresolveStatus status = run_presolver(presolver);

    std::cout << "removed "
              << presolver->stats->nnz_original - presolver->stats->nnz_reduced
              << " non-zeros during presolving." << std::endl;

    free_settings(stgs);
    free_presolver(presolver);

    std::cout << "Presolver ran succesfully in C++ build" << std::endl;
    return 0;
}
