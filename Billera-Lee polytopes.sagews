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

#Old Code: if abs(vector(j[1:]).dot_product(vector(directionVector)))/(vector(j[1:]).norm()*vector(directionVector).norm()) < genericity_tolerance:
def line_shelling(P, k, startingPoint, directionVector):
    genericity_tolerance=.00000001

    for j in P.inequalities_list():
        if vector(j[1:]).dot_product(vector(directionVector)) == 0:
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

    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]

    if k <= len(ShellingOrder) and k > 0:
        return [ShellingOrder[i][1] for i in range(k)]

    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]


#Q=polytopes.cyclic_polytope(4,7)
#D=generate_starting_point_direction(Q)
#line_shelling(Q,4,D[0],D[1])

#R=Polyhedron(vertices=[[0,0],[1,0],[0,1],[1,1]])
#line_shelling(R,4,[0,.5],[-1,0])




import itertools

#This method computes the boundary complex of a pseudomanifold (we will only use it for d-balls).
#Since we intent to work with complexes that are guaranteed to be pseudomanifolds
#this method does NOT check wether the imput facets are the facets of a pseudomanifold.
def boundary(facets):
    ridges=[]
    dimension=len(facets[0])
    for i in range(len(facets)):
        ridges.append(Combinations(facets[i],dimension-1).list())

    badridges=[]
    for ridgeList in ridges:
        for current in ridgeList:
            count=0
            for i in range(len(ridges)):
                if current in ridges[i]:
                    count+=1
            if count >= 2:
                badridges.append(current)

    return [i for j in ridges for i in j if i not in badridges]

#def boundary(facets):
#    ridges = []
#    for i in range(len(facets)):
#        for j in range(len(facets[i])):
#            ridge = []
#            for s in range(len(facets[i])):
#                if s != j:
#                    ridge.append(facets[i][s])
#            ridges.append(ridge)
#
#    print ridges
#    bad_ridges = []
#    for i in range(len(facets)):
#        for j in range(len(facets)):
#            if i < j:
#                vertices = []
#                for s in range(len(facets[i])):
#                    if (s>0 and facets[i][s] ==facets[j][s-1]) or (facets[i][s] == facets[j][s]) or (s<len(facets[j])-1 and facets[i][s]==facets[j][s+1]):
#                        vertices.append(facets[i][s])
#                if len(vertices) + 1 == len(facets[i]):
#                    bad_ridges.append(vertices)
#    print bad_ridges
#
#    for i in range(len(bad_ridges)):
#          ridges.remove(bad_ridges[i])
#          ridges.remove(bad_ridges[i])
#    return ridges

#We say a facet of a simplicial complex is an exterior facet if one of its ridges is a boundary ridge
#and it is an interior facet otherwise.  For these functions to work, the simplicial complex
#MUST BE SORTED e.g. [[1,2,3,4],[5,3,4,2],[3,4,5,6]] is NOT OKAY
def interior(facets):
    Interior=[]
    Boundary=boundary(facets)
    for face in facets:
        isInterior=True
        faceRidges=Combinations(face,len(face)-1).list()
        for b in Boundary:
            if b in faceRidges:
                isInterior=False
                break
        if isInterior:
            Interior.append(face)
    return Interior

def exterior(facets):
    Exterior=[]
    Boundary=boundary(facets)
    for face in facets:
        isInterior=True
        faceRidges=Combinations(face,len(face)-1).list()
        for b in Boundary:
            if b in faceRidges:
                isInterior=False
                break
        if not isInterior:
            Exterior.append(face)
    return Exterior

def exterior_of_degree(facets,degree):
    Exterior=[]
    Boundary=boundary(facets)
    for face in facets:
        exteriorDegree=0
        faceRidges=Combinations(face,len(face)-1).list()
        for b in Boundary:
            if b in faceRidges:
                exteriorDegree+=1
        if exteriorDegree>= degree:
            Exterior.append(face)
    return Exterior

def exteriorDegree(face,facets):
    exteriorDegree=0
    Boundary=boundary(facets)
    faceRidges=Combinations(face,len(face)-1).list()
    for b in Boundary:
        if b in faceRidges:
            exteriorDegree+=1
    return exteriorDegree

interior([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
exterior([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
boundary([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
exterior_of_degree([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]],2)
exteriorDegree([4,5,6],[[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
︡b0e0fdcf-cc01-496d-9d84-b359e8bbb389︡{"stdout":"[[2, 4, 6], [2, 7, 8]]\n[[1, 2, 6], [2, 3, 4], [4, 5, 6], [7, 8, 9], [1, 2, 7], [2, 3, 8]]\n[[1, 6], [3, 4], [4, 5], [5, 6], [7, 9], [8, 9], [1, 7], [3, 8]]\n[[4, 5, 6], [7, 8, 9]]\n2\n"}︡
︠05736660-bd6f-4070-a8c0-b8be1102b23b︠


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

    if dimension%2 == 1:
        for i in range(len(final_list)):
            for j in range(len(final_list[i])):
                final_list[i][j]+=2
            final_list[i].insert(0,2)
            final_list[i].insert(0,1)


    return final_list

#b=SimplicialComplex(ball_constructor([1,4,5,5,0],8))
#s=SimplicialComplex(boundary(ball_constructor([1,4,5,6,6],8)))
#b.h_vector()
#s.h_vector()
#s.g_vector()
#print "hello"

desiredFacets = ball_constructor([1,4,5,5,0], 8)
#print desiredFacets
exteriorDesiredFacets=exterior_of_degree(desiredFacets,6)
exteriorFacetIndices=[i for i in range(len(desiredFacets)) if desiredFacets[i] in exteriorDesiredFacets]
print desiredFacets
print exteriorFacetIndices
print len(desiredFacets)
print len(exteriorDesiredFacets)
︡3aa4416c-61ab-4f55-ac31-76bf265d4551︡{"stdout":"[[1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 4, 5, 6, 7, 9, 10], [1, 2, 3, 4, 5, 6, 7, 10, 11], [1, 2, 3, 4, 5, 6, 7, 11, 12], [1, 2, 3, 4, 5, 6, 7, 12, 13], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 10, 11], [1, 2, 3, 4, 5, 8, 9, 10, 11], [1, 2, 3, 4, 5, 7, 8, 11, 12], [1, 2, 3, 4, 5, 8, 9, 11, 12], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 5, 6, 7, 8, 10, 11], [1, 2, 3, 5, 6, 8, 9, 10, 11], [1, 2, 3, 6, 7, 8, 9, 10, 11], [1, 2, 3, 5, 6, 7, 8, 11, 12]]\n[0, 4, 9, 13, 14]\n15\n5\n"}︡
︠554e8b08-eb31-4786-b638-d83f137fc0e0︠
︠f686a673-2a08-43eb-9ef6-b56e56bd7c1b︠



#This generates rational cyclic polytopes. It takes as imput the dimension, and a set of base points on the moment curve.
def rational_cyclic_polytope(dimension, base_points):
    vertex_set = [[i^j for j in range(1,dimension+1)] for i in base_points]
    return Polyhedron(vertices = vertex_set)




#TESTING THE LINE SHELLING CODE


#dimension = 6
#g_vector = [1,8,0,0]

#This is an attempt to build Billera Lee polytopes using a non-standard method,
#i.e something that slightly differs from the original proof from Billera-Lee to prove that these
#polytopes indeed exist.  It takes a long time because we are picking large numbers to
#do computations with, guided by the intuition from the paper by Billera and Lee.

def BilleraLeePolytope(g_vector, dimension):
    numVertices = g_vector[1] + dimension + 1
    stopPoint = sum(g_vector)
    #print stopPoint

    #This is the simplicial complex that we want to find via a line shelling
    desiredFacets = ball_constructor(g_vector, dimension)

    #Here we note the "exterior degree" of the desired Facets
    desiredFacetDegrees=[exteriorDegree(face,desiredFacets) for face in desiredFacets]
    print desiredFacetDegrees

    #Here we construct the cyclic polytope that we will shell
    vertex_base = [i for i in range(numVertices)]
    #print vertex_base
    Q = rational_cyclic_polytope(dimension+1, vertex_base)
    #print Q.vertices_list()

    #Here we find the facets of Q that correspond to the ones we want.
    #desiredFacetCenters is a list of the centers of these facets
    #desiredFacetInequalities is a list of the facet normal directions of these facets
    #desiredFacetCenters = []
    desiredFacetInequalities = []
    ieqList=Q.inequalities_list()[:]
    for i in range(len(desiredFacets)):
        polytopeFacet = [vertex_base[j-1] for j in desiredFacets[i]]
        P = rational_cyclic_polytope(dimension+1, polytopeFacet)
        targetVertices=P.vertices_list()[:]
        inequality=[]
        for k in range(len(ieqList)):
            totalDotProduct=0
            for j in range(len(targetVertices)):
                totalDotProduct += abs(ieqList[k][0]+vector(ieqList[k][1:]).dot_product(vector(targetVertices[j])))
            if totalDotProduct==0:
                inequality = ieqList[k][1:]
        coordinates = list(P.center())
        #desiredFacetCenters.append(coordinates)
        desiredFacetInequalities.append(inequality)
        #desiredFacetInequalities.append(list(vector(inequality).normalized()))
    #print desiredFacetCenters
    desiredLengthCenters = []
    for i in range(len(desiredFacets)):
        polytopeFacet = [vertex_base[j-1] for j in desiredFacets[i] if not j==1]
        P = rational_cyclic_polytope(dimension+1, polytopeFacet)
        coordinates = list(P.center())
        desiredLengthCenters.append(coordinates)
    #We want a shelling that only captures the facets with centers in desiredFacetCenters.
    #We need to choose a starting point and direction that accomplishes this...
    #This has proven elusive so far.
    sumDegrees=sum(desiredFacetDegrees)
    weightedAverage=[]
    for j in range(len(desiredLengthCenters[0])):
        tempSum=0
        for i in range(len(desiredLengthCenters)):
            tempSum += desiredLengthCenters[i][j]*desiredFacetDegrees[i]/sumDegrees
        weightedAverage.append(tempSum)
    #print weightedAverage

    centerDirection=[vertex_base[0]-weightedAverage[k]*10000 for k in range(dimension+1)]
    outwardDirections=[-sum(i)/len(desiredFacetInequalities) for i in zip(*desiredFacetInequalities)]
    print centerDirection
    print outwardDirections
    weightFactor = dimension
    combinedDirection=[1/(weightFactor+1)*outwardDirections[k]+weightFactor/(weightFactor+1)*centerDirection[k] for k in range(dimension+1)]

    perturbationVector1=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    perturbationVector2=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]

    direction=[combinedDirection[i]-perturbationVector1[i] for i in range(dimension+1)]
    startingPoint=[weightedAverage[i]-perturbationVector2[i] for i in range(dimension+1)]


    #polytopeCenter=list(Q.center())
    #averageCenter=[sum(i)/len(desiredFacetCenters) for i in zip(*desiredFacetCenters)]
    #directionList=[[d[j]/d[0]-averageCenter[j]/d[0] for j in range(dimension+1)] for d in desiredFacetCenters if desiredFacetDegrees[desiredFacetCenters.index(d)] >= 8]
    #directionList.append([(1-averageCenter[k]) for k in range(dimension+1)])
    #print directionList
    #print averageCenter
    #specialVertex=[1 for i in range(dimension+1)]
    #preStartingPoint=[((averageCenter[0]-(2**17+116800))/(averageCenter[0]+1))*specialVertex[i]+((1+(2**17+116800))/(averageCenter[0]+1))*averageCenter[i] for i in range(dimension+1)]
    #perturbationVector1=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    #perturbationVector2=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    #preDirection = [sum(i) for i in zip(*directionList)]
    #preDirection = [-sum(i) for i in zip(*desiredFacetInequalities)]
    #direction = [preDirection[i]-perturbationVector1[i] for i in range(dimension+1)]
    #startingPoint = [specialVertex[i]-perturbationVector2[i] for i in range(dimension+1)]
    #print startingPoint
    #print direction
    #directions=[d for d in desiredFacetInequalities if desiredFacetDegrees[desiredFacetInequalities.index(d)] >= 5]
    #preDirection = [sum(i) for i in zip(*directions)]
    #direction = [preDirection[i]-perturbationVector1[i] for i in range(dimension+1)]


    #weightList = []
    #for i in range(len(initFacets)+1):
    #    if i == 0:
    #        weightList.append(len(initFacets))
    #    else:
    #        weightList.append(i/2)

    #totalWeight=0
    #for i in range(len(weightList)):
    #    totalWeight+=weightList[i]
    #startingPoint = list(Q.center())
    #startingPoint
    #for i in range(len(startingPoint)):
    #    startingPoint[i]=startingPoint[i]*weightList[0]
    #startingPoint
    #preDirection = Q.vertices_list()[0]
    #for i in range(len(initFacets)):
    #    for j in range(len(startingPoint)):
    #        startingPoint[j] += weightList[i+1]*initFacets[i][j]
    #for i in range(len(startingPoint)):
    #    startingPoint[i]=startingPoint[i]/totalWeight

    #perturbationVector=[ZZ.random_element(90,110)/100001  for i in range(len(startingPoint))]
    #direction = [Q.vertices_list()[0][i]-startingPoint[i]-perturbationVector[i] for i in range(len(startingPoint))]
    #perturbationVector2=[ZZ.random_element(90,110)/100001 for i in range(len(startingPoint))]
    #for i in range(len(startingPoint)):
    #    startingPoint[i]+= perturbationVector2[i]


    #vector([5,1,1]).norm()
    #vector([5,1,1]).norm()*vector([0,-1,1]).norm()
    #direction
    #[n(direction[i]/vector(direction).norm()) for i in range(len(direction))]
    #for j in Q.inequalities_list():
    #    v=vector(j[1:])
    #    l=[n(j[i]/v.norm()) for i in range(1,len(j))]
    #    print l
    #for j in Q.inequalities_list():
    #    print n(abs(vector(j[1:]).dot_product(vector(direction)))/(vector(j[1:]).norm()*vector(direction).norm()))
    #for i in range(len(temporal)):
    #    starting_point.append(temporal[i] + perturbationVector[i] )
    #direction = [ZZ.random_element(90,110)/1001 for i in range(dimension)]
    #for i in range(1,30):


    #Now we construct the shelling and calculate the g-vector of the corresponding sphere
    shelling_order = line_shelling(Q, stopPoint, startingPoint, direction )
    if not shelling_order == None:
        abstract_complex = []
        for i in shelling_order:
            facet = []
            for j in i:
                facet.append(j[0])
            abstract_complex.append(facet)
        #show(abstract_complex)

        #ball = SimplicialComplex(abstract_complex)
        #show(ball.h_vector())
        sphere = SimplicialComplex(boundary(abstract_complex))
        #show(sphere.h_vector())
        show(sphere.h_vector())
        show(sphere.g_vector())
    return



BilleraLeePolytope([1,4,5,5,0],8)
︡6ce0abd8-0ec0-4ef8-ba5c-9b60272fd9d9︡{"stdout":"[6, 5, 5, 5, 8, 4, 4, 5, 5, 7, 3, 3, 5, 6, 6]"}︡{"stdout":"\n"}︡{"stderr":"Error in lines 150-150\n"}︡{"stderr":"Traceback (most recent call last):\n  File \"/mnt/home/xWnefKLQ/.sagemathcloud/sage_server.py\", line 412, in execute\n    exec compile(block, '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\n  File \"\", line 29, in BilleraLeePolytope\n  File \"free_module_element.pyx\", line 502, in sage.modules.free_module_element.vector (build/cythonized/sage/modules/free_module_element.c:4929)\n  File \"/usr/local/sage/sage-5.10.rc1/local/lib/python2.7/site-packages/sage/modules/free_module.py\", line 951, in __call__\n    return self._element_class(self, x, coerce, copy)\n  File \"vector_integer_dense.pyx\", line 108, in sage.modules.vector_integer_dense.Vector_integer_dense.__init__ (build/cythonized/sage/modules/vector_integer_dense.c:3126)\n  File \"/usr/local/sage/sage-5.10.rc1/local/lib/python2.7/site-packages/sage/structure/sequence.py\", line 571, in __getitem__\n    if isinstance(n, slice):\n  File \"c_lib.pyx\", line 70, in sage.ext.c_lib.sage_python_check_interrupt (build/cythonized/sage/ext/c_lib.c:925)\nKeyboardInterrupt\n"}︡
︠5d14a65a-a635-4b40-ba59-38c54316f186︠
︠eeee17d2-c96a-41f4-9df1-cf141111c4a7︠

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
︡30b4222b-7ddd-4f15-b22e-e64306776e4d︡
︠512ba0b8-bd4d-4d76-ace5-a9ebd2696364︠
#DEMO OF HOW TO USE:
#Decompositions:
'The 4-binomial decompostion of 10 is given by:'
decomp(10, 4)
'The 4-binomial decompostion of 100 is given by:'
decomp(100, 4)
'The 4th shadow of 10 is:'
shadow_op(10, 4)
'The 4th shadow of 100 is:'
shadow_op(100, 4)
#M sequences:
#This is easy to check by hand!
'The list [1,4,9,15,13,10,10] is an M-sequence?'
is_M_sequence([1,4,9,15,13,10,10])
#The first entry is not 1
'The list [2,3,4,5,6,7] is an M-sequence?'
is_M_sequence([2,3,4,5,6,7])
#First check the conditions without the 7 disturbing, then add the 7 condition
'The set of values i such that [1,5,i,7] is an M-sequnece is:'
i = 0
stop = False
indicesTemp = []
while not stop:
    if not is_M_sequence([1,5,i]):
        stop = True
    else:
        indicesTemp.append(i)
        i+=1
indices = []
for j in indicesTemp:
    if is_M_sequence([1,5,j,7]):
        indices.append(j)
indices
#Now we construct a rational cyclic polytope
Q = rational_cyclic_polytope(8, [1,3,5,7,9,11,13,15,17])
Q
# Now we generate a direction and a starting point for a Line shelling of Q
startingPoint = list(Q.center)
#Next we generate a direction for the line shelling
facets = ball_constructor([1,5,10,0],6)

SimplicialComplex(facets)

︡e5befac3-ab46-42bd-a0fa-a21b661985a5︡{"stdout":"'The 4-binomial decompostion of 10 is given by:'\n[[5, 4], [4, 3], [2, 2]]\n'The 4-binomial decompostion of 100 is given by:'\n[[8, 4], [6, 3], [5, 2]]\n'The 4th shadow of 10 is:'\n8\n'The 4th shadow of 100 is:'\n49\n'The list [1,4,9,15,13,10,10] is an M-sequence?'\nTrue\n'The list [2,3,4,5,6,7] is an M-sequence?'\nFalse\n'The set of values i such that [1,5,i,7] is an M-sequnece is:'\n[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]\nA 8-dimensional polyhedron in ZZ^8 defined as the convex hull of 9 vertices\n"}︡{"stderr":"Error in lines 38-38\nTraceback (most recent call last):\n  File \"/mnt/home/xWnefKLQ/.sagemathcloud/sage_server.py\", line 412, in execute\n    exec compile(block, '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\nTypeError: 'sage.misc.cachefunc.CachedMethodCallerNoArgs' object is not iterable\n"}︡
︠f74f6ec5-e0b8-4723-b061-506aba0f53df︠

#Returns the product of all the elements in a list
def list_product(l):
    p=1
    for i in l:
        p *= i
    return p

#Implements Gales Evenness Criterion
def is_gales_even(facet, numTotalVertices):
    complement=range(1,numTotalVertices+1)
    for j in facet:
        complement.remove(j)
    result=True
    for i in range(len(complement)-1):
        if (complement[i+1]-complement[i])%2 == 0:
            result=False
            break
    return result


#This function takes in the first coordinates of the verticies of a facet of the cyclic polytope
#and returns the equality defining that facet
#input: for the facet ((1,1,1), (2,4,8), (5,25,125)), we input [1,2,5]
#output: the output is a d+1 length vector such that the homogenized versions of the vertices
#        in facet "dot" to zero with this vector
import itertools
def cyclic_facet_vector(facet):
    combVector=[]
    for i in range(len(facet)+1):
        combVector.append(Combinations(facet,i).list())
    combVector.reverse()
    #print combVector
    prodVector = [[list_product(j) for j in i] for i in combVector]
    #print prodVector
    facetVector=[]
    for i in range(len(prodVector)):
        facetVector.append(sum(prodVector[i])*(-1)**(i+1))
    #print facetVector
    return facetVector

#This functions gives the facets (in an abstract simplicial complex sense)
# of C(n,d)
def cyclic_polytope_facets(numVertices,dimension):
    potentialFacets=Combinations(range(1,numVertices+1),dimension).list()
    return [facet for facet in potentialFacets if is_gales_even(facet,numVertices)]

g_vector=[1,4,5,5,0]
dimension=8
numVertices = g_vector[1] + dimension + 1
stopPoint = sum(g_vector)
vertex_base = [i+1 for i in range(numVertices)]
Q = rational_cyclic_polytope(dimension+1, vertex_base)
print Q.vertices()
center=list(Q.center())
center.insert(0,1)

desiredFacets = ball_constructor(g_vector, dimension)
print desiredFacets

allCyclicFacets=cyclic_polytope_facets(numVertices,dimension+1)
print len(allCyclicFacets)

#For each facet we want to construct an inequality reflecting whether we want to be above or below that facet
ieqSystem=[]
for facet in allCyclicFacets:
    realizedFacet=[vertex_base[j-1] for j in facet]
    facetVector=cyclic_facet_vector(realizedFacet)
    #need to decide if we require >=0 or <=0
    #orientation=1 corresponds to the polytope facet being >=0
    orientation=1
    if vector(center).dot_product(vector(facetVector)) < 0:
        orientation=-1
    #now we reverse orientation for those facets in desired facets
    if facet in desiredFacets:
        orientation *= -1

    if orientation == 1:
        ieqSystem.append(facetVector)
    if orientation == -1:
        ieqSystem.append([-j for j in facetVector])

print ieqSystem

#Now we just need to find a strictly feasible point for ieqsystem >= 0
#Alternatively, we could adjust the constants by a small value and then we would only need to find a feasible point

#This is not the complete linear program since we only enter the first three values of our inequalities.  Just trying
#to get a feel for what this is doing.
p=MixedIntegerLinearProgram( maximization=True )
p.set_objective(p[1]+p[2]+p[3])
for ieq in ieqSystem:
    p.add_constraint( ieq[1]*p[1]+ieq[2]*p[2]+ieq[3]*p[3], min = -ieq[0] )
p.solve()
p.get_values(p[1])

#HOW CAN WE DO LP FEASIBILITY IN SAGE?

#cyclic_polytope_facets(6,4)

#print vector(cyclic_facet_vector([1,3,5,7])).dot_product(vector([1,5,25,125,625]))

︡88e31da1-94bf-4aad-b1dd-2a9444b296f7︡{"stdout":"(A vertex at (1, 1, 1, 1, 1, 1, 1, 1, 1), A vertex at (2, 4, 8, 16, 32, 64, 128, 256, 512), A vertex at (3, 9, 27, 81, 243, 729, 2187, 6561, 19683), A vertex at (4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144), A vertex at (5, 25, 125, 625, 3125, 15625, 78125, 390625, 1953125), A vertex at (6, 36, 216, 1296, 7776, 46656, 279936, 1679616, 10077696), A vertex at (7, 49, 343, 2401, 16807, 117649, 823543, 5764801, 40353607), A vertex at (8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728), A vertex at (9, 81, 729, 6561, 59049, 531441, 4782969, 43046721, 387420489), A vertex at (10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000), A vertex at (11, 121, 1331, 14641, 161051, 1771561, 19487171, 214358881, 2357947691), A vertex at (12, 144, 1728, 20736, 248832, 2985984, 35831808, 429981696, 5159780352), A vertex at (13, 169, 2197, 28561, 371293, 4826809, 62748517, 815730721, 10604499373))\n[[1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 4, 5, 6, 7, 9, 10], [1, 2, 3, 4, 5, 6, 7, 10, 11], [1, 2, 3, 4, 5, 6, 7, 11, 12], [1, 2, 3, 4, 5, 6, 7, 12, 13], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 10, 11], [1, 2, 3, 4, 5, 8, 9, 10, 11], [1, 2, 3, 4, 5, 7, 8, 11, 12], [1, 2, 3, 4, 5, 8, 9, 11, 12], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 5, 6, 7, 8, 10, 11], [1, 2, 3, 5, 6, 8, 9, 10, 11], [1, 2, 3, 6, 7, 8, 9, 10, 11], [1, 2, 3, 5, 6, 7, 8, 11, 12]]\n140\n[[362880, -1026576, 1172700, -723680, 269325, -63273, 9450, -870, 45, -1], [524160, -1464912, 1645196, -992816, 359121, -81417, 11634, -1014, 49, -1], [453600, -1271880, 1435212, -871786, 318143, -72989, 10598, -944, 47, -1], [554400, -1543320, 1723988, -1033430, 370881, -83349, 11802, -1020, 49, -1], [665280, -1840896, 2039028, -1208612, 427539, -94353, 13062, -1098, 51, -1], [786240, -2164608, 2380332, -1397332, 488117, -106001, 14378, -1178, 53, -1], [673920, -1862064, 2058516, -1217432, 429639, -94605, 13074, -1098, 51, -1], [842400, -2306520, 2517588, -1464406, 506177, -108689, 14582, -1184, 53, -1], [1029600, -2798280, 3022412, -1733738, 588735, -123669, 16158, -1272, 55, -1], [1235520, -3337344, 3572988, -2025428, 677313, -139545, 17802, -1362, 57, -1], [604800, -1670640, 1847156, -1093724, 387201, -85809, 11994, -1026, 49, -1], [739200, -2026960, 2218044, -1295564, 450819, -97809, 13326, -1106, 51, -1], [887040, -2417568, 2622592, -1514222, 519117, -110541, 14718, -1188, 53, -1], [-1048320, 2842464, -3060800, 1749698, -592095, 124005, -16170, 1272, -55, 1], [950400, -2575920, 2773348, -1586396, 538077, -113289, 14922, -1194, 53, -1], [1140480, -3072096, 3278400, -1853170, 619015, -127853, 16450, -1280, 55, -1], [-1347840, 3611808, -3825408, 2140382, -705453, 143241, -18042, 1368, -57, 1], [-1425600, 3804480, -4006452, 2225456, -727233, 146289, -18258, 1374, -57, 1], [-1684800, 4472640, -4674156, 2569384, -828191, 163709, -19994, 1466, -59, 1], [-2059200, 5424960, -5607044, 3036704, -960309, 185409, -22026, 1566, -61, 1], [943488, -2552976, 2743740, -1567360, 531509, -112073, 14810, -1190, 53, -1], [1179360, -3161736, 3353580, -1882970, 624935, -128413, 16470, -1280, 55, -1], [1441440, -3835224, 4023988, -2226918, 725593, -145765, 18202, -1372, 57, -1], [1729728, -4573440, 4754964, -2599204, 833483, -164129, 20006, -1466, 59, -1], [1572480, -4150128, 4306708, -2351772, 755113, -149521, 18442, -1378, 57, -1], [1921920, -5033552, 5165596, -2778956, 875451, -169377, 20334, -1474, 59, -1], [2306304, -6001824, 6101888, -3241118, 1004325, -190365, 22302, -1572, 61, -1], [2471040, -6393264, 6447236, -3389852, 1038309, -194505, 22554, -1578, 61, -1], [2965248, -7622496, 7613760, -3951170, 1189839, -218253, 24690, -1680, 63, -1], [3706560, -9435456, 9290580, -4729360, 1390025, -247793, 27170, -1790, 65, -1], [907200, -2430360, 2580804, -1457174, 489939, -103029, 13686, -1116, 51, -1], [1108800, -2948040, 3096796, -1723658, 569205, -117117, 15162, -1200, 53, -1], [1330560, -3515472, 3659412, -2012128, 654199, -132041, 16702, -1286, 55, -1], [-1572480, 4132656, -4268652, 2322584, -744921, 147801, -18306, 1374, -57, 1], [1425600, -3745080, 3867732, -2105998, 677215, -135149, 16918, -1292, 55, -1], [-1710720, 4465584, -4569228, 2457048, -777553, 152145, -18602, 1382, -57, 1], [-2021760, 5249232, -5328756, 2834736, -884591, 170077, -20354, 1474, -59, 1], [-2138400, 5528520, -5578668, 2945130, -910931, 173509, -20582, 1480, -59, 1], [-2527200, 6498360, -6504804, 3396414, -1035529, 193729, -22486, 1576, -61, 1], [-3088800, 7880040, -7796796, 4007618, -1197735, 218757, -24702, 1680, -63, 1], [1995840, -5129064, 5138028, -2694338, 830433, -158565, 19002, -1392, 57, -1], [-2395008, 6114960, -6067044, 3140304, -951923, 178129, -20846, 1486, -59, 1], [-2830464, 7187184, -7072668, 3619848, -1081405, 198745, -22762, 1582, -61, 1], [-2993760, 7568856, -7402068, 3758502, -1112629, 202573, -23002, 1588, -61, 1], [-3538080, 8895528, -8627292, 4330482, -1262927, 225741, -25078, 1688, -63, 1], [-4324320, 10784952, -10334532, 5103118, -1457729, 254261, -27482, 1796, -65, 1], [-3991680, 9925488, -9476652, 4665572, -1332327, 233481, -25518, 1698, -63, 1], [-4717440, 11664144, -11041620, 5371588, -1510385, 259745, -27770, 1802, -65, 1], [-5765760, 14139696, -13220252, 6323204, -1740291, 291921, -30366, 1914, -67, 1], [-7413120, 17944272, -16453956, 7667876, -2045421, 331737, -33354, 2034, -69, 1], [1572480, -4045296, 4075460, -2168360, 685515, -135453, 16890, -1290, 55, -1], [1965600, -5007480, 4974052, -2598126, 803173, -154609, 18718, -1384, 57, -1], [2402400, -6071720, 5961148, -3065810, 929691, -174909, 20622, -1480, 59, -1], [2882880, -7238016, 7036748, -3571412, 1065069, -196353, 22602, -1578, 61, -1], [2620800, -6567440, 6372076, -3231164, 965391, -179109, 20874, -1486, 59, -1], [3203200, -7962160, 7633124, -3809084, 1115709, -202209, 22946, -1586, 61, -1], [3843840, -9490528, 9006912, -4433522, 1276407, -226581, 25098, -1688, 63, -1], [4118400, -10106320, 9507708, -4629596, 1317267, -231189, 25362, -1694, 63, -1], [4942080, -12045216, 11215360, -5384750, 1505205, -258633, 27690, -1800, 65, -1], [6177600, -14902080, 13662092, -6425096, 1751463, -292509, 30378, -1914, 67, -1], [3931200, -9523560, 8819084, -4234994, 1197129, -210729, 23406, -1596, 61, -1], [4804800, -11542840, 10554516, -4983278, 1380015, -237237, 25662, -1700, 63, -1], [5765760, -13755312, 12444172, -5790968, 1575249, -265161, 28002, -1806, 65, -1], [6177600, -14644680, 13126972, -6039578, 1623405, -270249, 28278, -1812, 65, -1], [7413120, -17450064, 15471828, -7013008, 1850743, -301565, 30802, -1922, 67, -1], [9266400, -21580920, 18823428, -8347310, 2146641, -339969, 33702, -2040, 69, -1], [8648640, -20008344, 17305028, -7619798, 1956723, -311685, 31302, -1932, 67, -1], [10378368, -23837040, 20383164, -8835944, 2226453, -347049, 34026, -2046, 69, -1], [12972960, -29471976, 24774588, -10496002, 2575559, -390173, 37142, -2168, 71, -1], [17297280, -38575248, 31515572, -12871172, 3038217, -443121, 40698, -2298, 73, -1], [-1814400, 4407120, -4173228, 2118136, -649397, 126329, -15722, 1214, -53, 1], [-2217600, 5341680, -4996772, 2497328, -751575, 143049, -17358, 1302, -55, 1], [-2661120, 6365664, -5893728, 2907098, -860913, 160725, -19062, 1392, -57, 1], [-3144960, 7479072, -6864096, 3347446, -977411, 179357, -20834, 1484, -59, 1], [-2851200, 6777360, -6219324, 3036200, -889329, 164241, -19290, 1398, -57, 1], [-3421440, 8075808, -7333344, 3532038, -1017611, 184261, -21146, 1492, -59, 1], [-4043520, 9487584, -8538336, 4064682, -1154209, 205345, -23074, 1588, -61, 1], [-4276800, 9987840, -8927676, 4215756, -1186549, 209209, -23314, 1594, -61, 1], [-5054400, 11733120, -10392228, 4849116, -1344707, 232869, -25402, 1694, -63, 1], [-6177600, 14215680, -12425772, 5700748, -1548905, 261929, -27818, 1802, -65, 1], [-3991680, 9260208, -8210484, 3852448, -1081811, 191441, -21566, 1502, -59, 1], [-4790016, 11032416, -9675360, 4476450, -1235773, 214333, -23590, 1600, -61, 1], [-5660928, 12959136, -11259360, 5146350, -1399559, 238413, -25690, 1700, -63, 1], [-5987520, 13640832, -11768148, 5333964, -1437527, 242697, -25942, 1706, -63, 1], [-7076160, 16022016, -13691340, 6128940, -1626625, 269633, -28210, 1810, -65, 1], [-8648640, 19407744, -16357668, 7194668, -1869683, 302545, -30822, 1922, -67, 1], [-7983360, 17855136, -14988480, 6575230, -1709825, 278213, -28670, 1820, -65, 1], [-9434880, 20969568, -17430528, 7548722, -1932203, 308581, -31122, 1928, -67, 1], [-11531520, 25396512, -20812096, 8850486, -2216941, 345517, -33934, 2044, -69, 1], [-14826240, 32181984, -25789056, 10668202, -2590679, 390845, -37154, 2168, -71, 1], [-6652800, 14546640, -11921996, 5183576, -1359309, 226569, -24234, 1614, -61, 1], [-7983360, 17322912, -14028768, 6009614, -1548567, 252957, -26442, 1716, -63, 1], [-9434880, 20340576, -16305120, 6895330, -1749605, 280673, -28730, 1820, -65, 1], [-9979200, 21404160, -17025804, 7137388, -1794737, 285425, -28994, 1826, -65, 1], [-11793600, 25130880, -19782612, 8184476, -2025863, 316309, -31458, 1934, -67, 1], [-14414400, 30424320, -23590588, 9580140, -2320981, 353809, -34282, 2050, -69, 1], [-13305600, 27984480, -21604352, 8753426, -2123643, 325749, -31938, 1944, -67, 1], [-15724800, 32852640, -25089856, 10028094, -2393881, 360409, -34594, 2056, -69, 1], [-19219200, 39764960, -29897344, 11722394, -2737599, 402309, -37626, 2176, -71, 1], [-24710400, 50341920, -36928448, 14067206, -3184797, 453369, -41082, 2304, -73, 1], [-19958400, 40313520, -29324268, 11200160, -2573781, 376089, -35322, 2070, -69, 1], [-23587200, 47313360, -34019604, 12809704, -2895419, 415229, -38186, 2186, -71, 1], [-28828800, 57245040, -40475996, 14938928, -3302265, 462297, -41442, 2310, -73, 1], [-37065600, 72424080, -49872132, 17864888, -3827775, 519309, -45138, 2442, -75, 1], [-51891840, 98428464, -64620108, 22061248, -4519613, 588665, -49322, 2582, -77, 1], [4717440, -8990928, 7280748, -3299792, 925113, -166761, 19362, -1398, 57, -1], [5896800, -11091240, 8838396, -3930034, 1077611, -189389, 21374, -1496, 59, -1], [7207200, -13410360, 10544804, -4613774, 1241109, -213309, 23466, -1596, 61, -1], [8648640, -15948288, 12399972, -5351012, 1415607, -238521, 25638, -1698, 63, -1], [7862400, -14460720, 11222948, -4842620, 1284717, -218001, 23730, -1602, 61, -1], [9609600, -17480080, 13381452, -5678924, 1477287, -245049, 25998, -1706, 63, -1], [11531520, -20783904, 15727360, -6580118, 1682625, -273525, 28350, -1812, 65, -1], [12355200, -22082160, 16547284, -6849212, 1732185, -278649, 28626, -1818, 65, -1], [14826240, -26251488, 19439808, -7929802, 1970563, -310541, 31162, -1928, 67, -1], [18532800, -32351040, 23537316, -9400064, 2279421, -349569, 34074, -2046, 69, -1], [11793600, -20708280, 15272532, -6251534, 1574847, -254469, 26478, -1716, 63, -1], [14414400, -25018920, 18187468, -7316882, 1806441, -285285, 28938, -1824, 65, -1], [17297280, -29734416, 21353412, -8463664, 2053051, -317681, 31486, -1934, 67, -1], [18532800, -31578840, 22446756, -8798950, 2110843, -323309, 31774, -1940, 67, -1], [22239360, -37523952, 26341596, -10169256, 2395981, -359457, 34514, -2054, 69, -1], [27799200, -46209960, 31841244, -12024114, 2763119, -403429, 37646, -2176, 71, -1], [25945920, -42727752, 29195676, -10968746, 2521221, -370557, 35034, -2064, 69, -1], [31135104, -50754384, 34232148, -12658848, 2856455, -411145, 37982, -2182, 71, -1], [38918880, -62470008, 41325732, -14936862, 3285817, -460261, 41338, -2308, 73, -1], [51891840, -81131184, 51990780, -18138308, 3847515, -520065, 45150, -2442, 75, -1], [23587200, -35519760, 23139324, -8569672, 1990625, -301049, 29666, -1838, 65, -1], [28828800, -42830640, 27469076, -9992960, 2274843, -336369, 32334, -1950, 67, -1], [34594560, -50820192, 32163936, -11522066, 2576901, -373437, 35094, -2064, 69, -1], [37065600, -53891280, 33737292, -11952632, 2644845, -379617, 35394, -2070, 69, -1], [44478720, -63928224, 39481056, -13768782, 2992199, -420805, 38354, -2188, 71, -1], [55598400, -78520320, 47527308, -16205196, 3435697, -470569, 41722, -2314, 73, -1], [51891840, -72482544, 43513956, -14778352, 3137639, -432905, 38894, -2198, 71, -1], [62270208, -85941216, 50870880, -16998330, 3543169, -478933, 42070, -2320, 73, -1], [77837760, -105480576, 61146180, -19963500, 4058315, -534273, 45670, -2450, 75, -1], [103783680, -136316448, 76388928, -24077542, 4725413, -601181, 49742, -2588, 77, -1], [86486400, -109272720, 60259868, -19048616, 3805977, -498729, 42882, -2334, 73, -1], [103783680, -129397536, 70299360, -21854390, 4286595, -550413, 46290, -2460, 75, -1], [129729600, -158503680, 84235932, -25575916, 4893053, -612209, 50138, -2594, 77, -1], [172972800, -204131040, 104710016, -30688874, 5671911, -686469, 54474, -2736, 79, -1], [259459200, -284574960, 136954044, -37972304, 6687009, -775929, 59346, -2886, 81, -1]]"}︡{"stdout":"\n"}︡{"stderr":"Error in lines 70-70\nTraceback (most recent call last):\n  File \"/mnt/home/xWnefKLQ/.sagemathcloud/sage_server.py\", line 412, in execute\n    exec compile(block, '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\n  File \"mip.pyx\", line 1821, in sage.numerical.mip.MixedIntegerLinearProgram.solve (build/cythonized/sage/numerical/mip.c:10099)\n  File \"glpk_backend.pyx\", line 851, in sage.numerical.backends.glpk_backend.GLPKBackend.solve (build/cythonized/sage/numerical/backends/glpk_backend.cpp:7246)\nMIPSolverException: 'GLPK : Solution is undefined'\n"}︡
︠535aeeb0-a08f-4366-91c7-a81b65c53565︠
︡e9fbb8ea-cc6d-4104-9095-fc41fb2cfcfb︡
︠49f13941-8a64-4070-a820-414789f75201︠


