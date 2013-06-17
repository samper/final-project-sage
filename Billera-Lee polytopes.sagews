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

interior([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
exterior([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
boundary([[1,2,6],[2,3,4],[2,4,6],[4,5,6],[7,8,9],[1,2,7],[2,7,8],[2,3,8]])
︡d475a598-9edc-4f07-8ffb-1261fac16831︡{"stdout":"[[2, 4, 6], [2, 7, 8]]\n[[1, 2, 6], [2, 3, 4], [4, 5, 6], [7, 8, 9], [1, 2, 7], [2, 3, 8]]\n[[1, 6], [3, 4], [4, 5], [5, 6], [7, 9], [8, 9], [1, 7], [3, 8]]\n"}︡
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
s=SimplicialComplex(boundary(ball_constructor([1,4,5,6,6],8)))
#b.h_vector()
#s.h_vector()
s.g_vector()
#print "hello"
︡6a5c38b3-4311-41be-b6d0-2719874002fd︡{"stdout":"[1, 4, 5, 6, 6]"}︡
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
    vertices = g_vector[1] + dimension + 1
    stopPoint = sum(g_vector)
    #print stopPoint

    #This is the simplicial complex that we want to find via a line shelling
    facets = ball_constructor(g_vector, dimension)

    #Here we construct the cyclic polytope that we will shell
    vertex_base = [stopPoint^(i) for i in range(vertices)]
    #print vertex_base
    Q = rational_cyclic_polytope(dimension+1, vertex_base)
    #print Q.vertices_list()

    #Here we find the facets of Q that correspond to the ones we want.
    #desiredFacetCenters is a list of the centers of these facets
    #desiredFacetInequalities is a list of the facet normal directions of these facets
    desiredFacetCenters = []
    desiredFacetInequalities = []
    ieqList=Q.inequalities_list()[:]
    for i in range(len(facets)):
        polytopeFacet = [vertex_base[j-1] for j in facets[i]]
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
        desiredFacetCenters.append(coordinates)
        desiredFacetInequalities.append(inequality)
        #desiredFacetInequalities.append(list(vector(inequality).normalized()))
    #print desiredFacetCenters

    #We want a shelling that only captures the facets with centers in desiredFacetCenters.
    #We need to choose a starting point and direction that accomplishes this...
    polytopeCenter=list(Q.center())
    averageCenter=[sum(i)/len(desiredFacetCenters) for i in zip(*desiredFacetCenters)]
    directionList=[[d[j]/d[0]-averageCenter[j]/d[0] for j in range(dimension+1)] for d in desiredFacetCenters]
    #directionList.append([(1-averageCenter[k]) for k in range(dimension+1)])
    #print directionList
    print averageCenter
    specialVertex=[1 for i in range(dimension+1)]
    preStartingPoint=[((averageCenter[0]-(2**17+116800))/(averageCenter[0]+1))*specialVertex[i]+((1+(2**17+116800))/(averageCenter[0]+1))*averageCenter[i] for i in range(dimension+1)]
    perturbationVector1=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    perturbationVector2=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    preDirection = [sum(i) for i in zip(*directionList)]
    #preDirection = [-sum(i) for i in zip(*desiredFacetInequalities)]
    direction = [preDirection[i]-perturbationVector1[i] for i in range(dimension+1)]
    startingPoint = [preStartingPoint[i]-perturbationVector2[i] for i in range(dimension+1)]
    #print startingPoint
    #print direction



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
︡9ff4f6b8-bec7-4b59-95dc-a856c3e75ed9︡{"stdout":"[11942009262241/9, 1147436145919276337223894601/9, 145826808414745895549872188951322417644001/9, 18894354893777116826184036213817600311300087964965285001/9, 2451247387632267087864267798139939800485876868204951063986416494900001/9, 318038417329841963265357999658460271236128317726673660234944172858926239060314125001/9, 41264303051351914414003391843513762822216947808825585467497813757774320556213158423052789477500001/9, 5353892060291482213479368418903808853555409623212915042448805289693174238500646974351159638629209971773803125001/9, 694647887019656087867617578876720457809539331847226137806397421164223267553989028590323640914540949600753235513645210562500001/9]"}︡{"stdout":"\n"}︡{"tex":{"tex":"\\left[1, 4, 10, 15, 15, 15, 10, 4, 1\\right]","display":true}}︡{"tex":{"tex":"\\left[1, 3, 6, 5, 0\\right]","display":true}}︡
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
︡aed8d970-0872-4aa8-9c8f-3b51e91908e8︡
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

︡6fc4321f-fe9d-4f5f-9cb6-ef8972c0e935︡{"stdout":"'The 4-binomial decompostion of 10 is given by:'\n[[5, 4], [4, 3], [2, 2]]\n'The 4-binomial decompostion of 100 is given by:'\n[[8, 4], [6, 3], [5, 2]]\n'The 4th shadow of 10 is:'\n8\n'The 4th shadow of 100 is:'\n49\n'The list [1,4,9,15,13,10,10] is an M-sequence?'\nTrue\n'The list [2,3,4,5,6,7] is an M-sequence?'\nFalse\n'The set of values i such that [1,5,i,7] is an M-sequnece is:'\n[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]\nA 8-dimensional polyhedron in ZZ^8 defined as the convex hull of 9 vertices\nSimplicial complex with 12 vertices and 16 facets\n"}︡
︠f74f6ec5-e0b8-4723-b061-506aba0f53df︠

︡2d1ad23c-c7d1-4f96-bf4a-b7f77191b3cb︡
︠773e2e17-565b-49a8-9edf-be6935a0f518i︠
︠535aeeb0-a08f-4366-91c7-a81b65c53565︠


