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
def find_intersection(startingPoint,randomDirection,affinePlane):
    return -vector(affinePlane).dot_product(vector(startingPoint))/vector(affinePlane).dot_product(vector(randomDirection))

# This function finds the vertices of a facet of a polytope given the inequality that defines it.
def find_facet_vertices(inequality,homoVertices):
    facetVertices=[]
    for i in range(len(homoVertices)):
        if vector(homoVertices[i]).dot_product(vector(inequality))==0:
            facetVertices.append(homoVertices[i][1:])
    return facetVertices

# This function calculates the line shelling
# with direction of the outward pointing normal vector from the 0th facet
# plus a small random vector.
def line_shelling(P, k):
    dim=P.dim()
    homoVerts=P.vertices_list()[:]
    for i in range(len(homoVerts)):
        homoVerts[i].insert(0,1)

    startingFacet=find_facet_vertices(P.inequalities_list()[0],homoVerts)
    startingPoint=reduce(lambda a, b: [a[i]+b[i] for i in range(len(a))],startingFacet)
    startingPoint=[startingPoint[i]/len(startingFacet) for i in range(len(startingPoint))]
    startingPoint.insert(0,1)
    perturbationVector=[ZZ.random_element(90,110)/1001 for i in range(dim+1)]
    directionVector=[perturbationVector[i]-P.inequalities_list()[0][i] for i in range(dim+1)]

    ShellingOrder=[]
    for i in range(len(P.inequalities_list())):
        ShellingOrder.append([P.inequalities_list()[i],find_facet_vertices(P.inequalities_list()[i],homoVerts),find_intersection(startingPoint,directionVector,P.inequalities_list()[i])])
    ShellingOrder.sort(key=lambda x: x[2])
    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]
    if k <= ShellingOrder and k > 0:
        return [ShellingOrder[i][1] for i in range(k)]

    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]

#Q=polytopes.cyclic_polytope(4,12)
#line_shelling(Q,5)


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

# boundary([[1,2,3],[2,3,4],[1,3,4]])

import itertools
def sphere_constructor(g_vector, dimension):
    if not is_M_sequence(g_vector) or (dimension != 2*len(g_vector) and dimension != 2*len(g_vector)+1):
        print "Invalid input"
        return

    if dimension%2 == 0:
        num_link_vertices = g_vector[1]+2*dimension-1
        mu = dimension/2
    if dimension%2 != 0:
        num_link_vertices = g_vector[1]+2*dimension-2
        mu = (dimension-1)/2

    link_facets=Combinations(range(1,num_link_vertices),mu).list()
    for i in range(mu-1):
        bad_facets=[l for l in link_facets if l[i+1]-l[i] <= 1]
        for j in range(len(bad_facets)):
            link_facets.remove(bad_facets[j])



    return num_link_vertices, mu, link_facets




#Combinations(range(1,10),2).list()
#subsets_of_list(range(1,10),1)
sphere_constructor([1,2,3,3], 8)

︡faa5ad75-10ed-4872-a0c9-e5f4789e3a45︡{"stdout":"(17, 4, [[1, 3, 5, 7], [1, 3, 5, 8], [1, 3, 5, 9], [1, 3, 5, 10], [1, 3, 5, 11], [1, 3, 5, 12], [1, 3, 5, 13], [1, 3, 5, 14], [1, 3, 5, 15], [1, 3, 5, 16], [1, 3, 6, 8], [1, 3, 6, 9], [1, 3, 6, 10], [1, 3, 6, 11], [1, 3, 6, 12], [1, 3, 6, 13], [1, 3, 6, 14], [1, 3, 6, 15], [1, 3, 6, 16], [1, 3, 7, 9], [1, 3, 7, 10], [1, 3, 7, 11], [1, 3, 7, 12], [1, 3, 7, 13], [1, 3, 7, 14], [1, 3, 7, 15], [1, 3, 7, 16], [1, 3, 8, 10], [1, 3, 8, 11], [1, 3, 8, 12], [1, 3, 8, 13], [1, 3, 8, 14], [1, 3, 8, 15], [1, 3, 8, 16], [1, 3, 9, 11], [1, 3, 9, 12], [1, 3, 9, 13], [1, 3, 9, 14], [1, 3, 9, 15], [1, 3, 9, 16], [1, 3, 10, 12], [1, 3, 10, 13], [1, 3, 10, 14], [1, 3, 10, 15], [1, 3, 10, 16], [1, 3, 11, 13], [1, 3, 11, 14], [1, 3, 11, 15], [1, 3, 11, 16], [1, 3, 12, 14], [1, 3, 12, 15], [1, 3, 12, 16], [1, 3, 13, 15], [1, 3, 13, 16], [1, 3, 14, 16], [1, 4, 6, 8], [1, 4, 6, 9], [1, 4, 6, 10], [1, 4, 6, 11], [1, 4, 6, 12], [1, 4, 6, 13], [1, 4, 6, 14], [1, 4, 6, 15], [1, 4, 6, 16], [1, 4, 7, 9], [1, 4, 7, 10], [1, 4, 7, 11], [1, 4, 7, 12], [1, 4, 7, 13], [1, 4, 7, 14], [1, 4, 7, 15], [1, 4, 7, 16], [1, 4, 8, 10], [1, 4, 8, 11], [1, 4, 8, 12], [1, 4, 8, 13], [1, 4, 8, 14], [1, 4, 8, 15], [1, 4, 8, 16], [1, 4, 9, 11], [1, 4, 9, 12], [1, 4, 9, 13], [1, 4, 9, 14], [1, 4, 9, 15], [1, 4, 9, 16], [1, 4, 10, 12], [1, 4, 10, 13], [1, 4, 10, 14], [1, 4, 10, 15], [1, 4, 10, 16], [1, 4, 11, 13], [1, 4, 11, 14], [1, 4, 11, 15], [1, 4, 11, 16], [1, 4, 12, 14], [1, 4, 12, 15], [1, 4, 12, 16], [1, 4, 13, 15], [1, 4, 13, 16], [1, 4, 14, 16], [1, 5, 7, 9], [1, 5, 7, 10], [1, 5, 7, 11], [1, 5, 7, 12], [1, 5, 7, 13], [1, 5, 7, 14], [1, 5, 7, 15], [1, 5, 7, 16], [1, 5, 8, 10], [1, 5, 8, 11], [1, 5, 8, 12], [1, 5, 8, 13], [1, 5, 8, 14], [1, 5, 8, 15], [1, 5, 8, 16], [1, 5, 9, 11], [1, 5, 9, 12], [1, 5, 9, 13], [1, 5, 9, 14], [1, 5, 9, 15], [1, 5, 9, 16], [1, 5, 10, 12], [1, 5, 10, 13], [1, 5, 10, 14], [1, 5, 10, 15], [1, 5, 10, 16], [1, 5, 11, 13], [1, 5, 11, 14], [1, 5, 11, 15], [1, 5, 11, 16], [1, 5, 12, 14], [1, 5, 12, 15], [1, 5, 12, 16], [1, 5, 13, 15], [1, 5, 13, 16], [1, 5, 14, 16], [1, 6, 8, 10], [1, 6, 8, 11], [1, 6, 8, 12], [1, 6, 8, 13], [1, 6, 8, 14], [1, 6, 8, 15], [1, 6, 8, 16], [1, 6, 9, 11], [1, 6, 9, 12], [1, 6, 9, 13], [1, 6, 9, 14], [1, 6, 9, 15], [1, 6, 9, 16], [1, 6, 10, 12], [1, 6, 10, 13], [1, 6, 10, 14], [1, 6, 10, 15], [1, 6, 10, 16], [1, 6, 11, 13], [1, 6, 11, 14], [1, 6, 11, 15], [1, 6, 11, 16], [1, 6, 12, 14], [1, 6, 12, 15], [1, 6, 12, 16], [1, 6, 13, 15], [1, 6, 13, 16], [1, 6, 14, 16], [1, 7, 9, 11], [1, 7, 9, 12], [1, 7, 9, 13], [1, 7, 9, 14], [1, 7, 9, 15], [1, 7, 9, 16], [1, 7, 10, 12], [1, 7, 10, 13], [1, 7, 10, 14], [1, 7, 10, 15], [1, 7, 10, 16], [1, 7, 11, 13], [1, 7, 11, 14], [1, 7, 11, 15], [1, 7, 11, 16], [1, 7, 12, 14], [1, 7, 12, 15], [1, 7, 12, 16], [1, 7, 13, 15], [1, 7, 13, 16], [1, 7, 14, 16], [1, 8, 10, 12], [1, 8, 10, 13], [1, 8, 10, 14], [1, 8, 10, 15], [1, 8, 10, 16], [1, 8, 11, 13], [1, 8, 11, 14], [1, 8, 11, 15], [1, 8, 11, 16], [1, 8, 12, 14], [1, 8, 12, 15], [1, 8, 12, 16], [1, 8, 13, 15], [1, 8, 13, 16], [1, 8, 14, 16], [1, 9, 11, 13], [1, 9, 11, 14], [1, 9, 11, 15], [1, 9, 11, 16], [1, 9, 12, 14], [1, 9, 12, 15], [1, 9, 12, 16], [1, 9, 13, 15], [1, 9, 13, 16], [1, 9, 14, 16], [1, 10, 12, 14], [1, 10, 12, 15], [1, 10, 12, 16], [1, 10, 13, 15], [1, 10, 13, 16], [1, 10, 14, 16], [1, 11, 13, 15], [1, 11, 13, 16], [1, 11, 14, 16], [1, 12, 14, 16], [2, 4, 6, 8], [2, 4, 6, 9], [2, 4, 6, 10], [2, 4, 6, 11], [2, 4, 6, 12], [2, 4, 6, 13], [2, 4, 6, 14], [2, 4, 6, 15], [2, 4, 6, 16], [2, 4, 7, 9], [2, 4, 7, 10], [2, 4, 7, 11], [2, 4, 7, 12], [2, 4, 7, 13], [2, 4, 7, 14], [2, 4, 7, 15], [2, 4, 7, 16], [2, 4, 8, 10], [2, 4, 8, 11], [2, 4, 8, 12], [2, 4, 8, 13], [2, 4, 8, 14], [2, 4, 8, 15], [2, 4, 8, 16], [2, 4, 9, 11], [2, 4, 9, 12], [2, 4, 9, 13], [2, 4, 9, 14], [2, 4, 9, 15], [2, 4, 9, 16], [2, 4, 10, 12], [2, 4, 10, 13], [2, 4, 10, 14], [2, 4, 10, 15], [2, 4, 10, 16], [2, 4, 11, 13], [2, 4, 11, 14], [2, 4, 11, 15], [2, 4, 11, 16], [2, 4, 12, 14], [2, 4, 12, 15], [2, 4, 12, 16], [2, 4, 13, 15], [2, 4, 13, 16], [2, 4, 14, 16], [2, 5, 7, 9], [2, 5, 7, 10], [2, 5, 7, 11], [2, 5, 7, 12], [2, 5, 7, 13], [2, 5, 7, 14], [2, 5, 7, 15], [2, 5, 7, 16], [2, 5, 8, 10], [2, 5, 8, 11], [2, 5, 8, 12], [2, 5, 8, 13], [2, 5, 8, 14], [2, 5, 8, 15], [2, 5, 8, 16], [2, 5, 9, 11], [2, 5, 9, 12], [2, 5, 9, 13], [2, 5, 9, 14], [2, 5, 9, 15], [2, 5, 9, 16], [2, 5, 10, 12], [2, 5, 10, 13], [2, 5, 10, 14], [2, 5, 10, 15], [2, 5, 10, 16], [2, 5, 11, 13], [2, 5, 11, 14], [2, 5, 11, 15], [2, 5, 11, 16], [2, 5, 12, 14], [2, 5, 12, 15], [2, 5, 12, 16], [2, 5, 13, 15], [2, 5, 13, 16], [2, 5, 14, 16], [2, 6, 8, 10], [2, 6, 8, 11], [2, 6, 8, 12], [2, 6, 8, 13], [2, 6, 8, 14], [2, 6, 8, 15], [2, 6, 8, 16], [2, 6, 9, 11], [2, 6, 9, 12], [2, 6, 9, 13], [2, 6, 9, 14], [2, 6, 9, 15], [2, 6, 9, 16], [2, 6, 10, 12], [2, 6, 10, 13], [2, 6, 10, 14], [2, 6, 10, 15], [2, 6, 10, 16], [2, 6, 11, 13], [2, 6, 11, 14], [2, 6, 11, 15], [2, 6, 11, 16], [2, 6, 12, 14], [2, 6, 12, 15], [2, 6, 12, 16], [2, 6, 13, 15], [2, 6, 13, 16], [2, 6, 14, 16], [2, 7, 9, 11], [2, 7, 9, 12], [2, 7, 9, 13], [2, 7, 9, 14], [2, 7, 9, 15], [2, 7, 9, 16], [2, 7, 10, 12], [2, 7, 10, 13], [2, 7, 10, 14], [2, 7, 10, 15], [2, 7, 10, 16], [2, 7, 11, 13], [2, 7, 11, 14], [2, 7, 11, 15], [2, 7, 11, 16], [2, 7, 12, 14], [2, 7, 12, 15], [2, 7, 12, 16], [2, 7, 13, 15], [2, 7, 13, 16], [2, 7, 14, 16], [2, 8, 10, 12], [2, 8, 10, 13], [2, 8, 10, 14], [2, 8, 10, 15], [2, 8, 10, 16], [2, 8, 11, 13], [2, 8, 11, 14], [2, 8, 11, 15], [2, 8, 11, 16], [2, 8, 12, 14], [2, 8, 12, 15], [2, 8, 12, 16], [2, 8, 13, 15], [2, 8, 13, 16], [2, 8, 14, 16], [2, 9, 11, 13], [2, 9, 11, 14], [2, 9, 11, 15], [2, 9, 11, 16], [2, 9, 12, 14], [2, 9, 12, 15], [2, 9, 12, 16], [2, 9, 13, 15], [2, 9, 13, 16], [2, 9, 14, 16], [2, 10, 12, 14], [2, 10, 12, 15], [2, 10, 12, 16], [2, 10, 13, 15], [2, 10, 13, 16], [2, 10, 14, 16], [2, 11, 13, 15], [2, 11, 13, 16], [2, 11, 14, 16], [2, 12, 14, 16], [3, 5, 7, 9], [3, 5, 7, 10], [3, 5, 7, 11], [3, 5, 7, 12], [3, 5, 7, 13], [3, 5, 7, 14], [3, 5, 7, 15], [3, 5, 7, 16], [3, 5, 8, 10], [3, 5, 8, 11], [3, 5, 8, 12], [3, 5, 8, 13], [3, 5, 8, 14], [3, 5, 8, 15], [3, 5, 8, 16], [3, 5, 9, 11], [3, 5, 9, 12], [3, 5, 9, 13], [3, 5, 9, 14], [3, 5, 9, 15], [3, 5, 9, 16], [3, 5, 10, 12], [3, 5, 10, 13], [3, 5, 10, 14], [3, 5, 10, 15], [3, 5, 10, 16], [3, 5, 11, 13], [3, 5, 11, 14], [3, 5, 11, 15], [3, 5, 11, 16], [3, 5, 12, 14], [3, 5, 12, 15], [3, 5, 12, 16], [3, 5, 13, 15], [3, 5, 13, 16], [3, 5, 14, 16], [3, 6, 8, 10], [3, 6, 8, 11], [3, 6, 8, 12], [3, 6, 8, 13], [3, 6, 8, 14], [3, 6, 8, 15], [3, 6, 8, 16], [3, 6, 9, 11], [3, 6, 9, 12], [3, 6, 9, 13], [3, 6, 9, 14], [3, 6, 9, 15], [3, 6, 9, 16], [3, 6, 10, 12], [3, 6, 10, 13], [3, 6, 10, 14], [3, 6, 10, 15], [3, 6, 10, 16], [3, 6, 11, 13], [3, 6, 11, 14], [3, 6, 11, 15], [3, 6, 11, 16], [3, 6, 12, 14], [3, 6, 12, 15], [3, 6, 12, 16], [3, 6, 13, 15], [3, 6, 13, 16], [3, 6, 14, 16], [3, 7, 9, 11], [3, 7, 9, 12], [3, 7, 9, 13], [3, 7, 9, 14], [3, 7, 9, 15], [3, 7, 9, 16], [3, 7, 10, 12], [3, 7, 10, 13], [3, 7, 10, 14], [3, 7, 10, 15], [3, 7, 10, 16], [3, 7, 11, 13], [3, 7, 11, 14], [3, 7, 11, 15], [3, 7, 11, 16], [3, 7, 12, 14], [3, 7, 12, 15], [3, 7, 12, 16], [3, 7, 13, 15], [3, 7, 13, 16], [3, 7, 14, 16], [3, 8, 10, 12], [3, 8, 10, 13], [3, 8, 10, 14], [3, 8, 10, 15], [3, 8, 10, 16], [3, 8, 11, 13], [3, 8, 11, 14], [3, 8, 11, 15], [3, 8, 11, 16], [3, 8, 12, 14], [3, 8, 12, 15], [3, 8, 12, 16], [3, 8, 13, 15], [3, 8, 13, 16], [3, 8, 14, 16], [3, 9, 11, 13], [3, 9, 11, 14], [3, 9, 11, 15], [3, 9, 11, 16], [3, 9, 12, 14], [3, 9, 12, 15], [3, 9, 12, 16], [3, 9, 13, 15], [3, 9, 13, 16], [3, 9, 14, 16], [3, 10, 12, 14], [3, 10, 12, 15], [3, 10, 12, 16], [3, 10, 13, 15], [3, 10, 13, 16], [3, 10, 14, 16], [3, 11, 13, 15], [3, 11, 13, 16], [3, 11, 14, 16], [3, 12, 14, 16], [4, 6, 8, 10], [4, 6, 8, 11], [4, 6, 8, 12], [4, 6, 8, 13], [4, 6, 8, 14], [4, 6, 8, 15], [4, 6, 8, 16], [4, 6, 9, 11], [4, 6, 9, 12], [4, 6, 9, 13], [4, 6, 9, 14], [4, 6, 9, 15], [4, 6, 9, 16], [4, 6, 10, 12], [4, 6, 10, 13], [4, 6, 10, 14], [4, 6, 10, 15], [4, 6, 10, 16], [4, 6, 11, 13], [4, 6, 11, 14], [4, 6, 11, 15], [4, 6, 11, 16], [4, 6, 12, 14], [4, 6, 12, 15], [4, 6, 12, 16], [4, 6, 13, 15], [4, 6, 13, 16], [4, 6, 14, 16], [4, 7, 9, 11], [4, 7, 9, 12], [4, 7, 9, 13], [4, 7, 9, 14], [4, 7, 9, 15], [4, 7, 9, 16], [4, 7, 10, 12], [4, 7, 10, 13], [4, 7, 10, 14], [4, 7, 10, 15], [4, 7, 10, 16], [4, 7, 11, 13], [4, 7, 11, 14], [4, 7, 11, 15], [4, 7, 11, 16], [4, 7, 12, 14], [4, 7, 12, 15], [4, 7, 12, 16], [4, 7, 13, 15], [4, 7, 13, 16], [4, 7, 14, 16], [4, 8, 10, 12], [4, 8, 10, 13], [4, 8, 10, 14], [4, 8, 10, 15], [4, 8, 10, 16], [4, 8, 11, 13], [4, 8, 11, 14], [4, 8, 11, 15], [4, 8, 11, 16], [4, 8, 12, 14], [4, 8, 12, 15], [4, 8, 12, 16], [4, 8, 13, 15], [4, 8, 13, 16], [4, 8, 14, 16], [4, 9, 11, 13], [4, 9, 11, 14], [4, 9, 11, 15], [4, 9, 11, 16], [4, 9, 12, 14], [4, 9, 12, 15], [4, 9, 12, 16], [4, 9, 13, 15], [4, 9, 13, 16], [4, 9, 14, 16], [4, 10, 12, 14], [4, 10, 12, 15], [4, 10, 12, 16], [4, 10, 13, 15], [4, 10, 13, 16], [4, 10, 14, 16], [4, 11, 13, 15], [4, 11, 13, 16], [4, 11, 14, 16], [4, 12, 14, 16], [5, 7, 9, 11], [5, 7, 9, 12], [5, 7, 9, 13], [5, 7, 9, 14], [5, 7, 9, 15], [5, 7, 9, 16], [5, 7, 10, 12], [5, 7, 10, 13], [5, 7, 10, 14], [5, 7, 10, 15], [5, 7, 10, 16], [5, 7, 11, 13], [5, 7, 11, 14], [5, 7, 11, 15], [5, 7, 11, 16], [5, 7, 12, 14], [5, 7, 12, 15], [5, 7, 12, 16], [5, 7, 13, 15], [5, 7, 13, 16], [5, 7, 14, 16], [5, 8, 10, 12], [5, 8, 10, 13], [5, 8, 10, 14], [5, 8, 10, 15], [5, 8, 10, 16], [5, 8, 11, 13], [5, 8, 11, 14], [5, 8, 11, 15], [5, 8, 11, 16], [5, 8, 12, 14], [5, 8, 12, 15], [5, 8, 12, 16], [5, 8, 13, 15], [5, 8, 13, 16], [5, 8, 14, 16], [5, 9, 11, 13], [5, 9, 11, 14], [5, 9, 11, 15], [5, 9, 11, 16], [5, 9, 12, 14], [5, 9, 12, 15], [5, 9, 12, 16], [5, 9, 13, 15], [5, 9, 13, 16], [5, 9, 14, 16], [5, 10, 12, 14], [5, 10, 12, 15], [5, 10, 12, 16], [5, 10, 13, 15], [5, 10, 13, 16], [5, 10, 14, 16], [5, 11, 13, 15], [5, 11, 13, 16], [5, 11, 14, 16], [5, 12, 14, 16], [6, 8, 10, 12], [6, 8, 10, 13], [6, 8, 10, 14], [6, 8, 10, 15], [6, 8, 10, 16], [6, 8, 11, 13], [6, 8, 11, 14], [6, 8, 11, 15], [6, 8, 11, 16], [6, 8, 12, 14], [6, 8, 12, 15], [6, 8, 12, 16], [6, 8, 13, 15], [6, 8, 13, 16], [6, 8, 14, 16], [6, 9, 11, 13], [6, 9, 11, 14], [6, 9, 11, 15], [6, 9, 11, 16], [6, 9, 12, 14], [6, 9, 12, 15], [6, 9, 12, 16], [6, 9, 13, 15], [6, 9, 13, 16], [6, 9, 14, 16], [6, 10, 12, 14], [6, 10, 12, 15], [6, 10, 12, 16], [6, 10, 13, 15], [6, 10, 13, 16], [6, 10, 14, 16], [6, 11, 13, 15], [6, 11, 13, 16], [6, 11, 14, 16], [6, 12, 14, 16], [7, 9, 11, 13], [7, 9, 11, 14], [7, 9, 11, 15], [7, 9, 11, 16], [7, 9, 12, 14], [7, 9, 12, 15], [7, 9, 12, 16], [7, 9, 13, 15], [7, 9, 13, 16], [7, 9, 14, 16], [7, 10, 12, 14], [7, 10, 12, 15], [7, 10, 12, 16], [7, 10, 13, 15], [7, 10, 13, 16], [7, 10, 14, 16], [7, 11, 13, 15], [7, 11, 13, 16], [7, 11, 14, 16], [7, 12, 14, 16], [8, 10, 12, 14], [8, 10, 12, 15], [8, 10, 12, 16], [8, 10, 13, 15], [8, 10, 13, 16], [8, 10, 14, 16], [8, 11, 13, 15], [8, 11, 13, 16], [8, 11, 14, 16], [8, 12, 14, 16], [9, 11, 13, 15], [9, 11, 13, 16], [9, 11, 14, 16], [9, 12, 14, 16], [10, 12, 14, 16]])"}︡
︠852e3e67-6657-4cf3-a259-b9d803e641b8i︠
︡87b0fc8e-0e2d-4289-a47d-4db3c5f7106f︡
︠174e42e3-97db-4e96-8f63-bd85444a4613︠

