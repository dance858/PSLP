# PSLP — A Lightweight C Presolver for Linear Programs

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build and Test Status](https://img.shields.io/github/actions/workflow/status/dance858/PSLP/cmake.yml?branch=main)](https://github.com/dance858/PSLP/actions)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Coverage Status](https://coveralls.io/repos/github/dance858/PSLP/badge.svg?branch=main)](https://coveralls.io/github/dance858/PSLP?branch=main)


**PSLP** is a fast, dependency-free presolver for linear programs (LPs) of the form

$$
\begin{array}{ll}
\text{minimize} & c^T x \\
\text{subject to} & \underline{b} \leq Ax \leq \overline{b} \\
& \underline{x} \leq x \leq \overline{x}.
\end{array}
$$

**PSLP** simplifies LPs before they are passed to a solver by detecting redundant
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

1. **Initialization** — performed using `new_presolver()`, which creates and
   initializes internal data structures used by the presolver.  

2. **Presolve** — performed using `run_presolver()`, which executes the presolve
   routines to reduce the LP. The reduced LP is available in the struct
   `reduced_prob` of the presolver struct.

4. **Postsolve** — performed using `postsolve()`, which recovers a primal-dual
   solution to the original LP from a primal-dual solution to the reduced problem.

Detailed descriptions of each function and all associated data structures are
documented in the `PSLP` folder. A video tutorial is available [here](https://www.youtube.com/watch?v=ASYi21eCB-8).


---
## Citation
A paper is under preparation. If you wish to cite this software, you may use the following bibtex. 
```bibtex
@software{Cederberg2025,
  author = {Cederberg, Daniel and Boyd, Stephen},
  title = {PSLP — A Lightweight C Presolver for Linear Programs},
  year = {2025},
  url = {https://github.com/dance858/PSLP},
}
