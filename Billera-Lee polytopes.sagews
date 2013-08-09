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
︡b58e4132-60e3-4b52-80a0-98c1497b78aa︡{"stdout":"[[2, 4, 6], [2, 7, 8]]\n"}︡{"stdout":"[[1, 2, 6], [2, 3, 4], [4, 5, 6], [7, 8, 9], [1, 2, 7], [2, 3, 8]]\n"}︡{"stdout":"[[1, 6], [3, 4], [4, 5], [5, 6], [7, 9], [8, 9], [1, 7], [3, 8]]\n"}︡{"stdout":"[[4, 5, 6], [7, 8, 9]]\n"}︡{"stdout":"2\n"}︡
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
︡9ad61cce-0d79-4811-a1c2-47cc3cc36f3f︡{"stdout":"[[1, 2, 3, 4, 5, 6, 7, 8, 9], [1, 2, 3, 4, 5, 6, 7, 9, 10], [1, 2, 3, 4, 5, 6, 7, 10, 11], [1, 2, 3, 4, 5, 6, 7, 11, 12], [1, 2, 3, 4, 5, 6, 7, 12, 13], [1, 2, 3, 4, 5, 7, 8, 9, 10], [1, 2, 3, 4, 5, 7, 8, 10, 11], [1, 2, 3, 4, 5, 8, 9, 10, 11], [1, 2, 3, 4, 5, 7, 8, 11, 12], [1, 2, 3, 4, 5, 8, 9, 11, 12], [1, 2, 3, 5, 6, 7, 8, 9, 10], [1, 2, 3, 5, 6, 7, 8, 10, 11], [1, 2, 3, 5, 6, 8, 9, 10, 11], [1, 2, 3, 6, 7, 8, 9, 10, 11], [1, 2, 3, 5, 6, 7, 8, 11, 12]]\n"}︡{"stdout":"[0, 4, 9, 13, 14]\n"}︡{"stdout":"15\n"}︡{"stdout":"5\n"}︡
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
︡dbd03671-4781-47ed-b21e-796931fa47e7︡{"stdout":"[6, 5, 5, 5, 8, 4, 4, 5, 5, 7, 3, 3, 5, 6, 6]"}︡{"stdout":"\n[-4137500/77, -2725000/7, -35877500/11, -2279885000/77, -3109587500/11, -215139875000/77, -2180850042500/77, -22540161785000/77, -33799893087500/11]"}︡{"stdout":"\n[-611712/5, 1535464/5, -892396/3, 447358/3, -129503/3, 112336/15, -3838/5, 214/5, -1]\n"}︡{"tex":{"tex":"\\left[1, 5, 12, 15, 15, 15, 12, 5, 1\\right]","display":true}}︡{"tex":{"tex":"\\left[1, 4, 7, 3, 0\\right]","display":true}}︡
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



def BilleraLeePolytopeWithLP(g_vector, dimension, epsilon=0.01):
    #epsilon governs "how feasible" we require the LP solution to be
    numVertices = g_vector[1] + dimension + 1
    stopPoint = sum(g_vector)

    #This is the polytope we will eventually shell
    vertex_base = [i+1 for i in range(numVertices)]
    Q = rational_cyclic_polytope(dimension+1, vertex_base)

    #The polytope center that will serve as the starting point of the shelling
    center=list(Q.center())
    center.insert(0,1)

    #These are the facets that we want our shelling to pick out.
    #These are not written using the points of "vertex_base" yet.
    desiredFacets = ball_constructor(g_vector, dimension)

    #These are the facets (combinatorially) of our polytope Q
    allCyclicFacets=cyclic_polytope_facets(numVertices,dimension+1)

    #For each facet we want to construct an inequality reflecting whether we want to be above or below that facet
    ieqSystem=[]
    for facet in allCyclicFacets:
        realizedFacet=[vertex_base[j-1] for j in facet] #here we construct the actual facet from the combinatorial one
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

    #Now we just need to find a strictly feasible point for ieqsystem >= 0
    #We use an LP class of SAGE to do this
    successCheck=-1
    count=0
    while successCheck <= 0 and count <= 10:
        p=MixedIntegerLinearProgram()
        p.set_objective(None)
        for ieq in ieqSystem:
            p.add_constraint( sum(ieq[i]*p[i] for i in range(1,dimension+2)), min = -ieq[0]+epsilon*(10**count))
        _ = [ p.set_min(p[i], None) for i in range(1,dimension+2) ]

        #p.show()
        s = p.solve()
        solutionList=[1] #This list will be our feasible point
        for i in range(1,dimension+2):
            solutionList.append(p.get_values(p[i]))

        #This is a check.  If the LP was successful, all these values should be strictly positive.
        lpTest=[]
        for ieq in ieqSystem:
            lpTest.append(vector(ieq).dot_product(vector(solutionList)))
        successCheck = min(lpTest)
        count = count+1
    #print count

    #Now our direction becomes the vector from the center to our feasible point
    direction=[solutionList[i]-center[i] for i in range(1,dimension+2)]

    #Now we are ready to calculate the line shelling
    shelling_order = line_shelling(Q, stopPoint, center[1:], direction )

    #Here we check the g vector of the calculated shelling
    if not shelling_order == None:
        abstract_complex = []
        for i in shelling_order:
            facet = []
            for j in i:
                facet.append(j[0])
            abstract_complex.append(facet)

        sphere = SimplicialComplex(boundary(abstract_complex))
        show(sphere.h_vector())
        show(sphere.g_vector())

    return shelling_order


shelling = BilleraLeePolytopeWithLP([1,4,5,5,0],9)


︡1a66b1b2-e9d8-4732-b01d-edf9a515b1ce︡{"stdout":"3"}︡{"stdout":"\n"}︡{"tex":{"tex":"\\left[1, 5, 10, 15, 15, 15, 15, 10, 5, 1\\right]","display":true}}︡{"tex":{"tex":"\\left[1, 4, 5, 5, 0\\right]","display":true}}︡
︠535aeeb0-a08f-4366-91c7-a81b65c53565︠
︡e9fbb8ea-cc6d-4104-9095-fc41fb2cfcfb︡
︠49f13941-8a64-4070-a820-414789f75201︠


try:
    p = MixedIntegerLinearProgram( maximization=True )
    x,y= p[0], p[1]
    p.set_objective(None) #maximize this
    p.add_constraint(-x + 4*y <= 6)
    p.add_constraint(-x -2*y <= -5)
    #p.add_constraint(x <= -5)
    s = p.solve(log=0)
    p.get_values([x,y])
    val123 = "good"
except sage.numerical.mip.MIPSolverException:
    val123 = "oops"


print val123
︡a3f51044-abd0-4a4f-8334-db3487bbd31a︡{"stdout":"[1.3333333333333333, 1.8333333333333333]\n"}︡{"stdout":"good\n"}︡
︠085a69f8-a12c-4762-9cfd-6135879299b6︠

p = MixedIntegerLinearProgram( maximization=True )
x,y= p[0], p[1]
p.set_objective(None) #maximize this
p.add_constraint(-x + 4*y <= 6)
p.add_constraint(-x -2*y <= -5)
#p.add_constraint(x <= -5)
s = p.solve()
s
p.get_values([x,y])
︡80277841-37a2-4091-9e3c-8b90634d9514︡{"stdout":"0.0\n"}︡{"stdout":"[1.3333333333333333, 1.8333333333333333]\n"}︡
︠a34b9c58-ccee-4d9e-8062-b2fd01e88cb3︠









