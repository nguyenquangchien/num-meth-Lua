# Numerical methods with Lua programming

The computer programs accompanying the book [_Numerical Methods for Nonlinear Engineering Models_](https://link.springer.com/book/10.1007/978-1-4020-9920-5) by John R. Hauser (2009, Springer). The programming language used is Lua, which according to the book author, is "a modern scripting language that is very easy to program and understand by the reader". However, the code need auxiliary modules. For example, the generation of triangular grids is managed by the [EasyMesh](http://www-dinma.univ.trieste.it/nirftc/research/easymesh/Default.htm) software. Plotting can be performed with [`gnuplot`](http://www.gnuplot.info) and code editing with [SciTE](http://scintilla.sourceforge.net/SciTE.html). However, such dependence on seprate tools may raise problems of compatibility when collaborating on code development.

In this code repository, I try to use features the of web notebook based library [`iTorch`](https://github.com/facebookarchive/iTorch) to help delivering a visual experience for user. Currently, iTorch is archived but still being used by many users.

## Contents
1. Introduction to Nonlinear Engineering Problems and Models
2. Numerical Fundamentals and Computer Programming
3. Roots of Nonlinear Equations
4. Solving Sets of Equations: Linear and Nonlinear
5. Numerical Derivatives and Numerical Integration
6. Interpolation
7. Curve Fitting and Data Plotting
8. Statistical Methods and Basic Statistical Functions
9. Data Models and Parameter Estimation
10. Differential Equations: Initial Value Problems
11. Differential Equations: Boundary Value Problems
12. Partial Differential Equations: Finite Difference Approaches
13. Partial Differential Equations: The Finite Element Method
