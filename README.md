# PSLP — A Lightweight C Presolver for Linear Programs

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build and Test Status](https://img.shields.io/github/actions/workflow/status/dance858/PSLP/cmake.yml?branch=main)](https://github.com/dance858/PSLP/actions)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Coverage Status](https://coveralls.io/repos/github/dance858/PSLP/badge.svg?branch=main)](https://coveralls.io/github/dance858/PSLP?branch=main)


**PSLP** is a fast, dependency-free presolver for linear programs (LPs) and convex
quadratic programs (QPs) of the form

$$
\begin{array}{ll}
\text{minimize} & \tfrac{1}{2} x^T P x + c^T x \\
\text{subject to} & \underline{b} \leq Ax \leq \overline{b} \\
& \underline{x} \leq x \leq \overline{x},
\end{array}
$$

with $P$ symmetric positive semidefinite ($P = 0$ for an LP).

**PSLP** simplifies LPs and QPs before they are passed to a solver by detecting redundant
constraints and fixing variables etc. It is written in C99 and is lightweight with no external
dependencies. Happy presolving!

---
## Installation
**PSLP** is built using CMake. If you wish to use **PSLP** in your own CMake-based project, a simple approach is to do as follows.
First, create a folder `third_party` in your project root and add PSLP as a submodule by running the following command
in your terminal:
```code
git submodule add https://github.com/dance858/PSLP.git third_party/PSLP
```
Next, add **PSLP** as a subdirectory by adding the following line to your CMakeLists.txt:
```code
add_subdirectory(third_party/PSLP EXCLUDE_FROM_ALL)
```
Finally, link your project to PSLP by adding the following line to your CMakeLists.txt
(replace YOUR_PROJECT with your actual target name):
```code
target_link_libraries(YOUR_PROJECT PRIVATE PSLP)
```

After this, you should be able to import the public API using `#include "PSLP_API.h"`.

---
## API
The public C API is defined in the header file `PSLP_API.h`. The other files 
in that folder contain some trivial public data structures. The API consists of three main operations:

1. **Initialization** — performed using `new_presolver_from_problem()`, which
   creates and initializes internal data structures used by the presolver. The
   problem is given in the struct `ProblemData`, with the quadratic
   objective matrix $P$ as the full symmetric matrix in CSR form (left unset for
   an LP) and an optional constant objective term `obj_offset`. For an LP,
   `new_presolver()` takes the problem data as arguments.

2. **Presolve** — performed using `run_presolver()`, which executes the presolve
   routines to reduce the problem. The reduced problem is available as a
   `ProblemData` in the field `reduced_prob` of the presolver struct (including
   the reduced $P$ in `Px`, `Pi`, `Pp`, `Pnnz` and the objective offset
   `obj_offset`, which accounts for the input offset and the variables fixed by
   the presolver).

3. **Postsolve** — performed using `postsolve()`, which recovers a primal-dual
   solution to the original problem from a primal-dual solution to the reduced
   problem, where the bound multipliers are $z = Px + c - A^T y$.

In addition, `map_original_sol_to_reduced()` maps a primal-dual point of the
original LP to the reduced LP, e.g. to warm start a solver on the reduced
problem from a solution of the original one. See `PSLP_API.h` for details.

Detailed descriptions of each function and all associated data structures are
documented in the `PSLP` folder. A video tutorial is available [here](https://www.youtube.com/watch?v=ASYi21eCB-8).
Note that since the video was recorded, `PSLP_inf` has been removed. Infinite variable bounds or constraint limits
should now be represented using the IEEE-754 encoding. (For example, in C, use `INFINITY` defined in `math.h`).


---
## Citation
If you find this repository useful, please consider giving it a star. If you wish to cite this work, you may use the following bibtex. 
```bibtex
@misc{Cederberg26
  title = {{Presolving for {GPU}-Accelerated First-Order LP Solvers}},
  author = {Daniel Cederberg and Stephen Boyd},
  year = {2026},
  eprint = {2604.23951},
  archivePrefix = {arXiv},
  primaryClass = {math.OC},
  note = {arXiv:2604.23951}
}
