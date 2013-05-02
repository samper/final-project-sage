Final project for the sage class
================================
Currently I am deciding what I want to do. I have four different research projects for which some coding
would be useful and I plan to implement one of them. 
Possible titles: 

1. Billera-Lee polytopes and regularity of stacked triangulations
2. The combinatorial coalgebra of PS-ear decomposable complexes
3. Random polytopes approximating spheres: computational aspects. 
4. Unimodality of h-vectors of simplicial polytopes: dimensions 16-19.


1. Billera Lee polytopes and regularity of stacked triangulations
=================================================================
The goal of this projects is to produce examples to support a recent conjecture by S. Murai and E. Nevo, 
that states that if P is a simplicial polytope with g_r= 0, where (g_0,g_1, \dots g_[d/2]) is the simplicial 
g-vector, then the stacked triangulation of P (the they showed exists) is regular, that is, it is the projection 
of the lower pointing faces of a polytope Q with dim Q = dim P + 1 that projects to P. 

Regular triangulations are highly studied and very interesting for many reasons. I will briefly mention : 
- They are good for computing invariants in computational aspects of discrete geometry(via shellability)
- The realization space of a polytope and slight deformations of the vertices are better understood in terms of regular triangulations.

The stacked triangulation of a polytope with g_2 = 0 is easily shown to be regular. The key behind the proof is that 
this polytopes are easily constructed inductively. On the other hand, polytopes with those characteristics do not arise
naturally in the study of geometry, but are interesting for theoretical reasons. 

The only big family of polytopes I know that has g_r=0 is the family of Billera-Lee polytopes. These are polytopes that 
were constructed to classify the family of all face numbers of simplicial polytopes. To construct them, one choose a 
suitable g-vector (the entries of this vector satisfy some prescribed inequalites) and constructs this polytope as the
boundary of a line shelling subcomplex of a special cyclic polytope of one dimension higher than the desired dimension. 

