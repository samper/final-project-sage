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
︡8ef4cb73-238d-43c7-a0ed-3b30bfeb3110︡
︠6a089b0c-7368-41c0-9a7e-785dd3bab053︠
is_M_sequence([1,23,45,34])
︡bb62279d-2a5e-476e-8d28-ec51b081e680︡{"stdout":"True\n"}︡
︠35fbf10c-5b62-489a-b006-9793a1a1ce35︠
dim=8
numVerts=30
Q = polytopes.cyclic_polytope(dim,numVerts)
type(Q.vertices_list()[0])

︡e0a851c2-1f4c-4945-ac56-cf89b397d0c8︡{"stdout":"<type 'list'>\n"}︡
︠44bf6803-9bd4-468a-bcb8-7e5366176d90︠
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

︡aca36ef0-23de-4651-b428-45ac46dc899b︡
︠ac7879e3-c0a8-473f-ac63-396439106d62︠

︡8a0d6630-c3ad-4c9d-8e59-8d2ab48fb33b︡
︠aeb4b701-3a74-451b-ba76-0f47d1994c12︠
#This method computes the g-vector of the boudary of a shelling.Right now it only works for cyclic polytopes (or polytopes whose projection on the first coordinate is injective at the level of vertices).
def shelling_facets(P,k):
    facets=[]
    for i in line_shelling(P,k):
        facet = []
        for j in i:
            facet.append(j[0])
        facets.append(facet)


    return facets


︡61378f7b-3509-4a66-9a8e-ed5af3fe689c︡
︠79133403-81ee-4681-96c5-e18d8cb0b6ba︠
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

boundary(shelling_facets(Q,5))

︡863d78bc-ad74-4fb5-acbc-0cc4700354fc︡{"stdout":"[[1, 2, 3, 4, 5, 6, 7], [0, 1, 3, 4, 5, 6, 7], [0, 1, 2, 3, 5, 6, 7], [0, 1, 2, 3, 4, 5, 7], [0, 1, 2, 3, 4, 5, 6], [1, 2, 4, 5, 6, 7, 19], [0, 1, 4, 5, 6, 7, 19], [0, 1, 2, 5, 6, 7, 19], [0, 1, 2, 4, 5, 6, 19], [1, 2, 3, 4, 6, 7, 19], [0, 1, 3, 4, 6, 7, 19], [0, 1, 2, 3, 6, 7, 19], [0, 1, 2, 3, 4, 7, 19], [0, 1, 2, 3, 4, 6, 19], [2, 3, 4, 5, 6, 7, 19], [0, 3, 4, 5, 6, 7, 19], [0, 2, 3, 5, 6, 7, 19], [0, 2, 3, 4, 5, 7, 19], [0, 2, 3, 4, 5, 6, 19], [1, 2, 4, 5, 7, 8, 19], [0, 2, 4, 5, 7, 8, 19], [0, 1, 4, 5, 7, 8, 19], [0, 1, 2, 5, 7, 8, 19], [0, 1, 2, 4, 7, 8, 19], [0, 1, 2, 4, 5, 8, 19], [0, 1, 2, 4, 5, 7, 8]]"}︡
︠8d0f903b-2af3-4f16-b1f0-fa3db1f0788c︠
boundary([[1,2,3],[2,3,4],[1,3,4]])
︡30e25666-11fc-44cc-a697-c72681af3469︡{"stdout":"[[1, 2], [2, 4], [1, 4]]\n"}︡
︠95c28218-91af-45e6-977a-6ddbf1d2ae32r︠
#This method returns the g-vector of the boundary complex of a simplicial ball obtained by truncating a line shelling of a cyclic polytope.
def shelling_g_vector(P, k):
    Delta = SimplicialComplex(boundary(shelling_facets(P, k)))
    return Delta.g_vector()

for i in range(30):
    shelling_g_vector(Q, i)
︡e6a9b292-cc08-4bdb-8fd1-dd51ce9a489c︡
︠2e6eba68-c791-4094-b8e3-7020d53b2a80i︠
︡55bfa0c1-eb77-43ab-b385-94c5c1c7fddc︡
︠9b857ef4-9aac-47f2-a41a-6b7120d53200︠

