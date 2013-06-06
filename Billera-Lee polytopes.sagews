︠b6c127e3-5019-4ae0-a0a5-8a8b1741a104︠
#Computes the canonical binomial decomposition of n given k.
def decomp(n, k):
    list=[]
    if k==0 or n <=0:
        return []
    elif k == 1:
        list.append([n,1])
    else:
        stop = false
        a = k-1
        while not stop:
            if binomial(a+1, k) > n:
                list.append([a,k])

                list.extend(decomp(n-binomial(a, k), k-1))
                stop = true
            else:
                a+=1
    return list


#This computes the k-th shadow operator of n.
def shadow_op(n,k):
    if n < k:
        return 0
    sum = 0
    list = decomp(n,k)
    for i in range(len(list)):
        sum += binomial(list[i][0]-1,list[i][1]-1)
    return sum

# This function verifies if a list of integers is an M-sequence, or equivalently, if it is the g-vector of a simplicial polytope.
def is_M_sequence(list):
    if not list[0] == 1:
        return False
    for i in range(len(list)-1):
        if i == 0:
            if list[i+1]<0:
                return False

        else:
            if shadow_op(list[i+1], i+1) > list[i]:
                return False
    return True

#is_M_sequence([1,23,45,34])

# Now we start with line shellings.
# This returns the scalar distance in the randomDirection from the affine plane to the point.
# affinePlane is a list of length n+1 where the 0th coordinate defines the intercept.
# startingPoint and randomDirection are lists of length n
def find_intersection(startingPoint,randomDirection,affinePlane):
    return -(affinePlane[0]+vector(affinePlane[1:]).dot_product(vector(startingPoint)))/vector(affinePlane[1:]).dot_product(vector(randomDirection))

# This function finds the vertices of a facet of a polytope given the inequality that defines it.
def find_facet_vertices(inequality,homoVertices):
    facetVertices=[]
    for i in range(len(homoVertices)):
        if vector(homoVertices[i]).dot_product(vector(inequality))==0:
            facetVertices.append(homoVertices[i][1:])
    return facetVertices

# This function generates a starting point
# and direction for a line shelling that is
# roughly the outward pointing normal from the
# 0th facet
def generate_starting_point_direction(P):
    dim=P.dim()
    homoVerts=P.vertices_list()[:]
    for i in range(len(homoVerts)):
        homoVerts[i].insert(0,1)

    startingFacet=find_facet_vertices(P.inequalities_list()[0],homoVerts)
    startingPoint=reduce(lambda a, b: [a[i]+b[i] for i in range(len(a))],startingFacet)
    startingPoint=[startingPoint[i]/len(startingFacet) for i in range(len(startingPoint))]
    perturbationVector=[ZZ.random_element(90,110)/1001 for i in range(dim)]
    directionVector=[perturbationVector[i]-P.inequalities_list()[0][i+1] for i in range(dim)]

    return [startingPoint, directionVector]


# This function calculates the line shelling
# with direction directionVector and starting point startingPoint
def line_shelling(P, k, startingPoint, directionVector):
    genericity_tolerance=.00001

    for j in P.inequalities_list():
        if abs(vector(j[1:]).dot_product(vector(directionVector))) < genericity_tolerance:
            print "Genericity of directionVector failed!"
            return

    dim=P.dim()
    homoVerts=P.vertices_list()[:]
    for i in range(len(homoVerts)):
        homoVerts[i].insert(0,1)

    ShellingOrder=[]
    for i in range(len(P.inequalities_list())):
        ShellingOrder.append([P.inequalities_list()[i],find_facet_vertices(P.inequalities_list()[i],homoVerts),find_intersection(startingPoint,directionVector,P.inequalities_list()[i])])
    ShellingOrder.sort(key=lambda x: x[2])

    for i in range(len(ShellingOrder)-1):
        if (ShellingOrder[i+1][2]-ShellingOrder[i][2]) < genericity_tolerance:
            print "Genericity of Shelling Order failed!"
            return

    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]

    if k <= len(ShellingOrder) and k > 0:
        return [ShellingOrder[i][1] for i in range(k)]

    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]


<<<<<<< HEAD
=======
#Q=polytopes.cyclic_polytope(4,7)
#D=generate_starting_point_direction(Q)
#line_shelling(Q,4,D[0],D[1])

#R=Polyhedron(vertices=[[0,0],[1,0],[0,1],[1,1]])
#line_shelling(R,4,[0,.5],[-1,.5])




>>>>>>> f6cf5781015ac51b70399241e7f8ba595cae093b


#This method computes the g-vector of the boudary of a shelling.Right now it only works for cyclic polytopes (or polytopes whose projection on the first coordinate is injective at the level of vertices).
def shelling_facets(P,k):
    facets=[]
    for i in line_shelling(P,k):
        facet = []
        for j in i:
            facet.append(j[0])
        facets.append(facet)
    return facets


#This method returns the g-vector of the boundary complex of a simplicial ball obtained by truncating a line shelling of a cyclic polytope.
def shelling_g_vector(P, k):
    Delta = SimplicialComplex(boundary(shelling_facets(P, k)))
    return Delta.g_vector()


#This method computes the boundary complex of a pseudomanifold (we will only use it for d-balls). Since we intent to work with complexes that are guaranteed to be pseudomanifolds this method does NOT check wether the imput facets are the facets of a pseudomanifold.
def boundary(facets):
    ridges = []
    for i in range(len(facets)):
        for j in range(len(facets[i])):
            ridge = []
            for s in range(len(facets[i])):
                if s != j:
                    ridge.append(facets[i][s])
            ridges.append(ridge)

    bad_ridges = []
    for i in range(len(facets)):
        for j in range(len(facets)):
            if i < j:
                vertices = []
                for s in range(len(facets[i])):
                    if (s>0 and facets[i][s] ==facets[j][s-1]) or (facets[i][s] == facets[j][s]) or (s<len(facets[j])-1 and facets[i][s]==facets[j][s+1]):
                        vertices.append(facets[i][s])
                if len(vertices) + 1 == len(facets[i]):
                    bad_ridges.append(vertices)

    for i in range(len(bad_ridges)):
          ridges.remove(bad_ridges[i])
          ridges.remove(bad_ridges[i])
    return ridges

<<<<<<< HEAD
=======
# boundary([[1,2,3,4],[2,3,4,5],[3,4,5,6]])
>>>>>>> f6cf5781015ac51b70399241e7f8ba595cae093b


import itertools
#This method returns the list of facets of a simplicial ball, whose boundary has the g-vector that we imput and has the dimension passed as a parameter:
#Warning: For this to properly work we need and not print an exeption the dimension and the g-vector have to be compatible: i.e it must have
def ball_constructor(g_vector, dimension):
    if not is_M_sequence(g_vector) or (dimension != 2*len(g_vector)-1 and dimension != 2*len(g_vector)-2):
        print "Invalid input"
        return


    if dimension%2 == 0:
        num_link_vertices = g_vector[1]+2*dimension+1
        mu = dimension/2
    if dimension%2 != 0:
        num_link_vertices = g_vector[1]+2*dimension
        mu = (dimension-1)/2

    link_facets=Combinations(range(1,num_link_vertices),mu).list()
    for i in range(mu-1):
        bad_facets=[l for l in link_facets if l[i+1]-l[i] <= 1]
        for j in range(len(bad_facets)):
            link_facets.remove(bad_facets[j])

    ordered_link_facets=[]
    for i in range(mu):
        if i==mu-1:
            same_degree_facets=[l for l in link_facets if l[i]==2*i+1]
        else:
            same_degree_facets=[l for l in link_facets if (l[i]==2*i+1 and l[i+1]>2*i+3)]
        ordered_link_facets.insert(0,same_degree_facets)
    same_degree_facets=[l for l in link_facets if l[0]!=1]
    ordered_link_facets.append(same_degree_facets)


    for i in range(len(ordered_link_facets)):
        for j in range(mu):
            ordered_link_facets[i].sort(key=lambda x: x[j])


    final_list = []
    for i in range(len(g_vector)):
        for j in range(g_vector[i]):
            final_list.append(ordered_link_facets[i][j])

    for i in range(len(final_list)):
        a = len(final_list[i])
        for j in range(a):
            final_list[i].insert(a-j, final_list[i][a-j-1]+1)

    if dimension%2 == 0:
        for i in range(len(final_list)):
            for j in range(len(final_list[i])):
                final_list[i][j]+=1
            final_list[i].insert(0,1)




    return final_list

#This generates rational cyclic polytopes. It takes as imput the dimension, and a set of base points on the moment curve.
def rational_cyclic_polytope(dimension, base_points):
    vertex_set = [[i^j for j in range(1,dimension+1)] for i in base_points]
    return Polyhedron(vertices = vertex_set)


︡a25e4dbf-b828-4aad-9c02-a9d17528d9d5︡{"stdout":"A 4-dimensional polyhedron in RDF^4 defined as the convex hull of 6 vertices\n"}︡
︠ea6a4b92-0f6a-4624-965e-2f1caec5fe54︠
Q = polytopes.cyclic_polytope(5, 11)
show(line_shelling(Q,10))
facets = shelling_facets(Q,10)
show(facets)
show(SimplicialComplex(facets))
show(shelling_g_vector(Q,10))
︡b7dfcfd1-05a8-4495-bf7f-3018d3eae515︡{"tex":{"tex":"\\left[\\left[\\left[0, 0, 0, 0, 0\\right], \\left[1, 1, 1, 1, 1\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[4, 16, 64, 256, 1024\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[1, 1, 1, 1, 1\\right], \\left[2, 4, 8, 16, 32\\right], \\left[4, 16, 64, 256, 1024\\right], \\left[5, 25, 125, 625, 3125\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[4, 16, 64, 256, 1024\\right], \\left[5, 25, 125, 625, 3125\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[5, 25, 125, 625, 3125\\right], \\left[6, 36, 216, 1296, 7776\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[6, 36, 216, 1296, 7776\\right], \\left[7, 49, 343, 2401, 16807\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[3, 9, 27, 81, 243\\right], \\left[4, 16, 64, 256, 1024\\right], \\left[5, 25, 125, 625, 3125\\right], \\left[6, 36, 216, 1296, 7776\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[7, 49, 343, 2401, 16807\\right], \\left[8, 64, 512, 4096, 32768\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[3, 9, 27, 81, 243\\right], \\left[4, 16, 64, 256, 1024\\right], \\left[6, 36, 216, 1296, 7776\\right], \\left[7, 49, 343, 2401, 16807\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[2, 4, 8, 16, 32\\right], \\left[3, 9, 27, 81, 243\\right], \\left[8, 64, 512, 4096, 32768\\right], \\left[9, 81, 729, 6561, 59049\\right]\\right], \\left[\\left[0, 0, 0, 0, 0\\right], \\left[3, 9, 27, 81, 243\\right], \\left[4, 16, 64, 256, 1024\\right], \\left[7, 49, 343, 2401, 16807\\right], \\left[8, 64, 512, 4096, 32768\\right]\\right]\\right]","display":true}}︡{"tex":{"tex":"\\left[\\left[0, 1, 2, 3, 4\\right], \\left[0, 1, 2, 4, 5\\right], \\left[0, 2, 3, 4, 5\\right], \\left[0, 2, 3, 5, 6\\right], \\left[0, 2, 3, 6, 7\\right], \\left[0, 3, 4, 5, 6\\right], \\left[0, 2, 3, 7, 8\\right], \\left[0, 3, 4, 6, 7\\right], \\left[0, 2, 3, 8, 9\\right], \\left[0, 3, 4, 7, 8\\right]\\right]","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|10|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|10|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\left[1, 5, 4\\right]","display":true}}︡
︠852e3e67-6657-4cf3-a259-b9d803e641b8i︠
Q=polytopes.cyclic_polytope(4,12)
show(line_shelling(Q,5))
complex_1 = shelling_facets(Q,5)
show(complex_1)
show(SimplicialComplex(complex_1))

<<<<<<< HEAD
ball = ball_constructor([1,12,20,30,10], 8)
delta = SimplicialComplex(ball)
show(delta)
boundary_complex =SimplicialComplex(boundary(ball))
show(boundary_complex)
show(boundary_complex.g_vector())

︡cef6f078-8fae-4e23-b9db-5368aa38f2ea︡{"tex":{"tex":"\\left[\\left[\\left[0, 0, 0, 0\\right], \\left[1, 1, 1, 1\\right], \\left[2, 4, 8, 16\\right], \\left[3, 9, 27, 81\\right]\\right], \\left[\\left[0, 0, 0, 0\\right], \\left[1, 1, 1, 1\\right], \\left[3, 9, 27, 81\\right], \\left[4, 16, 64, 256\\right]\\right], \\left[\\left[1, 1, 1, 1\\right], \\left[2, 4, 8, 16\\right], \\left[3, 9, 27, 81\\right], \\left[4, 16, 64, 256\\right]\\right], \\left[\\left[1, 1, 1, 1\\right], \\left[2, 4, 8, 16\\right], \\left[4, 16, 64, 256\\right], \\left[5, 25, 125, 625\\right]\\right], \\left[\\left[0, 0, 0, 0\\right], \\left[1, 1, 1, 1\\right], \\left[4, 16, 64, 256\\right], \\left[5, 25, 125, 625\\right]\\right]\\right]","display":true}}︡{"tex":{"tex":"\\left[\\left[0, 1, 2, 3\\right], \\left[0, 1, 3, 4\\right], \\left[1, 2, 3, 4\\right], \\left[1, 2, 4, 5\\right], \\left[0, 1, 4, 5\\right]\\right]","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|vertex|\\phantom{\\verb!x!}\\verb|set|\\phantom{\\verb!x!}\\verb|(0,|\\phantom{\\verb!x!}\\verb|1,|\\phantom{\\verb!x!}\\verb|2,|\\phantom{\\verb!x!}\\verb|3,|\\phantom{\\verb!x!}\\verb|4,|\\phantom{\\verb!x!}\\verb|5)|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|5|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|21|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|73|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|21|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|293|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\left[1, 12, 20, 30, 10\\right]","display":true}}︡
︠953a70d9-2e21-4a42-89fc-53ef910328d5︠
ball = ball_constructor([1,12,20,0,0], 8)
beta = SimplicialComplex(ball)
show(beta)
beta_skel = beta.n_skeleton(7-3)
show(beta_skel)
gamma= SimplicialComplex(boundary(ball))
gamma_skel = gamma.n_skeleton(7-3)
show(gamma_skel)
︡0ad5b3d4-e680-4613-9777-c3e1b3510c6e︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|21|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|33|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|21|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|1666|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡{"tex":{"tex":"\\verb|Simplicial|\\phantom{\\verb!x!}\\verb|complex|\\phantom{\\verb!x!}\\verb|with|\\phantom{\\verb!x!}\\verb|21|\\phantom{\\verb!x!}\\verb|vertices|\\phantom{\\verb!x!}\\verb|and|\\phantom{\\verb!x!}\\verb|1666|\\phantom{\\verb!x!}\\verb|facets|","display":true}}︡
︠05736660-bd6f-4070-a8c0-b8be1102b23b︠
=======
#revlex_list_sort([5,1,3,4],[6,0,3,4])
#Combinations(range(1,10),2).list()
#subsets_of_list(range(1,10),1)
#ball = ball_constructor([1,9,13,0,0], 8)
#delta = SimplicialComplex(ball)
#delta_skel = delta.n_skeleton(7-3)
#delta_skel
# boundary(ball_constructor([1,9,20,10], 8))
#gamma= SimplicialComplex(boundary(ball))
#gamma.g_vector()
#gamma_skel = gamma.n_skeleton(7-3)
#gamma_skel
︡7976ab98-a14e-4772-b0bf-21cea75f1ccd︡{"stdout":"[[[0, 0, 0, 0], [1, 1, 1, 1], [2, 4, 8, 16], [3, 9, 27, 81]], [[0, 0, 0, 0], [1, 1, 1, 1], [3, 9, 27, 81], [4, 16, 64, 256]], [[1, 1, 1, 1], [2, 4, 8, 16], [3, 9, 27, 81], [4, 16, 64, 256]], [[1, 1, 1, 1], [2, 4, 8, 16], [4, 16, 64, 256], [5, 25, 125, 625]]]\n"}︡
︠512ba0b8-bd4d-4d76-ace5-a9ebd2696364︠

︠f74f6ec5-e0b8-4723-b061-506aba0f53df︠

︠852e3e67-6657-4cf3-a259-b9d803e641b8i︠
︠953a70d9-2e21-4a42-89fc-53ef910328d5︠
>>>>>>> f6cf5781015ac51b70399241e7f8ba595cae093b

