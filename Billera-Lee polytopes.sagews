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

is_M_sequence([1,23,45,34])

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

# This function re
def line_shelling(P, k):
    homoVerts=P.vertices_list()[:]
    for i in range(len(homoVerts)):
        homoVerts[i].insert(0,1)

    startingFacet=find_facet_vertices(P.inequalities()[0],homoVerts)
    startingPoint=reduce(lambda a, b: [a[i]+b[i] for i in range(len(a))],startingFacet)
    startingPoint=[startingPoint[i]/len(startingFacet) for i in range(len(startingPoint))]
    startingPoint.insert(0,1)
    perturbationVector=[ZZ.random_element(90,110)/1001 for i in range(dim+1)]
    directionVector=[perturbationVector[i]+P.inequalities()[0][i] for i in range(dim+1)]
    ShellingOrder=[]
    for i in range(len(P.inequalities())):
        ShellingOrder.append([P.inequalities()[i],find_facet_vertices(P.inequalities()[i],homoVerts),find_intersection(startingPoint,directionVector,P.inequalities()[i])])
    ShellingOrder.sort(key=lambda x: x[2])
    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]
    if k <= ShellingOrder and k > 0:
        return [ShellingOrder[i][1] for i in range(k)]

    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]


#This method computes the g-vector of the boudary of a shelling.Right now it only works for cyclic polytopes (or polytopes whose projection on the first coordinate is injective at the level of vertices).
def shelling_facets(P,k):
    facets=[]
    for i in line_shelling(P,k):
        facet = []
        for j in i:
            facet.append(j[0])
        facets.append(facet)


    return facets


#This method computes the boundary complex of a pseudomanifold (we will only use it for d-balls). Since we intent to work with complexes that are guaranteed to be pseudomanifolds this method does NOT check wether the imput facets are the facets of a pseudomanifold.
def boundary(facets):
    ridges = []
    for i in range(len(facets)):
        for j in range(len(facets[i])):
            ridge = []
            for s in range(len(facets[i])):
                if not s == j:
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

boundary([[1,2,3],[2,3,4],[1,3,4]])

#This method returns the g-vector of the boundary complex of a simplicial ball obtained by truncating a line shelling of a cyclic polytope.
def shelling_g_vector(P, k):
    Delta = SimplicialComplex(boundary(shelling_facets(P, k)))
    return Delta.g_vector()

dim=8
numVerts=30
Q = polytopes.cyclic_polytope(dim,numVerts)
type(Q.vertices_list()[0])
3+4
Q
show(shelling_g_vector(Q,15))#Computes the canonical binomial decomposition of n given k.
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

is_M_sequence([1,23,45,34])

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

# This function ca
def line_shelling(P, k):
    homoVerts=P.vertices_list()[:]
    for i in range(len(homoVerts)):
        homoVerts[i].insert(0,1)

    startingFacet=find_facet_vertices(P.inequalities()[0],homoVerts)
    startingPoint=reduce(lambda a, b: [a[i]+b[i] for i in range(len(a))],startingFacet)
    startingPoint=[startingPoint[i]/len(startingFacet) for i in range(len(startingPoint))]
    startingPoint.insert(0,1)
    perturbationVector=[ZZ.random_element(90,110)/1001 for i in range(dim+1)]
    directionVector=[perturbationVector[i]+P.inequalities()[0][i] for i in range(dim+1)]
    ShellingOrder=[]
    for i in range(len(P.inequalities())):
        ShellingOrder.append([P.inequalities()[i],find_facet_vertices(P.inequalities()[i],homoVerts),find_intersection(startingPoint,directionVector,P.inequalities()[i])])
    ShellingOrder.sort(key=lambda x: x[2])
    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]
    if k <= ShellingOrder and k > 0:
        return [ShellingOrder[i][1] for i in range(k)]

    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]


#This method computes the g-vector of the boudary of a shelling.Right now it only works for cyclic polytopes (or polytopes whose projection on the first coordinate is injective at the level of vertices).
def shelling_facets(P,k):
    facets=[]
    for i in line_shelling(P,k):
        facet = []
        for j in i:
            facet.append(j[0])
        facets.append(facet)


    return facets


#This method computes the boundary complex of a pseudomanifold (we will only use it for d-balls). Since we intent to work with complexes that are guaranteed to be pseudomanifolds this method does NOT check wether the imput facets are the facets of a pseudomanifold.
def boundary(facets):
    ridges = []
    for i in range(len(facets)):
        for j in range(len(facets[i])):
            ridge = []
            for s in range(len(facets[i])):
                if not s == j:
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

boundary([[1,2,3],[2,3,4],[1,3,4]])

#This method returns the g-vector of the boundary complex of a simplicial ball obtained by truncating a line shelling of a cyclic polytope.
def shelling_g_vector(P, k):
    Delta = SimplicialComplex(boundary(shelling_facets(P, k)))
    return Delta.g_vector()

dim=4
numVerts=10
Q = polytopes.cyclic_polytope(dim,numVerts)



def sphere_constructor(g_vector, dimension):
    if is_M_sequence(g_vector) and (dimension == 2*len(g_vector) or dimension == 2*len(g_vector)+1) :
        return True

    return false

sphere_constructor([1,2,3,3], 8)

︡cfebf2ba-dd36-4cc4-8d43-ae7ee094634b︡{"stdout":"True\n[[1, 2], [2, 4], [1, 4]]\n<type 'list'>\n7\nA 8-dimensional polyhedron in QQ^8 defined as the convex hull of 30 vertices\n"}︡{"tex":{"tex":"\\left[1, 5, 7, 2\\right]","display":true}}︡{"stdout":"True\n[[1, 2], [2, 4], [1, 4]]\nTrue\n"}︡
︠852e3e67-6657-4cf3-a259-b9d803e641b8︠
︡87b0fc8e-0e2d-4289-a47d-4db3c5f7106f︡
︠174e42e3-97db-4e96-8f63-bd85444a4613︠

