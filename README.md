# PSLP — A Lightweight C Presolver for Linear Programs

[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Build Status](https://img.shields.io/github/actions/workflow/status/yourusername/PSLP/ci.yml?branch=main)](https://github.com/yourusername/PSLP/actions)
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
