Final project for the sage class
================================

This repository contains my final project for the special topics course about Sage. It contains a PDF that explains the problem
and some methods that are beeing used. 

Briefly speaking, the goal of the project is (the precise definitions will be available soon in the pfd file): 

1. Write a function that returns a cyclic polytope of dimension d with n>d vertices (This is already implemented in the polyhedal class in sage)
2. Write a function that verifies wether a finite integer sequence is a Macaulay sequence. (Done in the worksheet Macualay sequence check)
3. Write a function that returns a line shelling of a polytope, given the polytope, an interior point and an (approximate) direction
4. Produce code that constructs a Billera-Lee polytope given a suitable g-vector.
5. Write a function that returns the r-stacked triangulation of a given polytope with g_r = 0
6. Given the combinatorial data of an r-stacked triangulation, write a function that writes the linear program that verifies stackedness of the triangulation

What has been done so far: 
==========================
- The program takes a finite sequence of integer numbers and determines wether it is an M-sequence or not, that is, wethere it 
is the g-vector of a simplicial polytope or not. 
- It also produces line shellings of a polytope in a direction that is approximately normal to the first facet of the polytope, 
that is the first facet returned by the facet method in the polyhedral class
- It constructs Billera-Lee spheres with a prescribed g-vector, following the original combinatorial construction of Billera and Lee
- The program verifies that if g_r = 0, the triangulation induced by the cyclic polytope on the Billera-Lee sphere is (r-1)-stacked

The experiments done using this program encouraged a proof that stacked triangulations of Billera-Lee polytopes are
indeed regular. To my knowledge a new result and the first general result for the conjecture for values bigger than r. 

Goals: 
======
- Use polynomial optimization methods to determine the appropiate entries of the cyclic polytope to realize the Billera-Lee sphere as a polytope
- Study the geometry of the feasible set of solutions for regular triangulations of a Billera-Lee polytope.
