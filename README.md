# PSLP — A Lightweight C Presolver for Linear Programs

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build Status](https://img.shields.io/github/actions/workflow/status/dance858/PSLP/ci.yml?branch=main)](https://github.com/dance858/PSLP/actions)
[![Language: C](https://img.shields.io/badge/Language-C-blue.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Coverage Status](https://coveralls.io/repos/github/XXX/badge.svg?branch=master)](https://coveralls.io/github/XXXbranch=master)

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
---
## API
The public C API is defined in the header file `PSLP/API.h"
`, which contains all public data structures and function declarations.
After installing PSLP, you can include it in your project with:
```c
#include <PSLP/API.h"
>
```
The API consists of three main operations:

1. **Initialization** — performed using `new_presolver()`, which creates and
   initializes internal data structures used by the presolver.  

2. **Presolve** — performed using `run_presolver()`, which executes the presolve
   routines to reduce the LP. The reduced LP is available in the struct
   `reduced_problem` of the presolver object.

4. **Postsolve** — performed using `postsolve()`, which recovers a primal-dual
   solution to the original LP from a primal-dual solution to the reduced problem.
   (The exact convention we use for the dual variables can be found in our XXX.)

Detailed descriptions of each function and all associated data structures are
documented in `PSLP_"API.h"
`. A video tutorial is available at XXX.

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
