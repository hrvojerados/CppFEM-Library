# C-FEM-Library
I'm designing a complete FEM library for C++. My goal is to make it efficient and more importantly simple and flexible to use.

## Current plan 
I will probably start working on this project in March (academic reasons) but in the meantime I plan on using this README file to write down my ideas.
Currently I want to rewrite my python code in https://github.com/hrvojerados/FiniteElementMethods to C++ or at least parts of it. From my current experience with FEM the most computationally intensive part is solving the linear system. So maybe the best first course of action would be to look at ways to parallelise that process. I will definitely implement a SparseMatrix class that is much more memory efficient since in FEM sparse matricies are really common. 
As for generating meshes, I plan to reasearch what C++ has already implemented before implementation. Generating meshes isn't really that computationally intensive so it's definitely more important to generate a good mesh (maybe add some scaling inside it) rather than generate a poor mesh quicky. Adding mesh adaptivity also adds complexity to this project.
I should also look into graphing and evaluating the solution.
