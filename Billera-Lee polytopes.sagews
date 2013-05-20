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
︡97aa8aee-6a7a-45eb-ab8b-3cd61934d473︡
︠6a089b0c-7368-41c0-9a7e-785dd3bab053︠
is_M_sequence([1,23,45,34])
︡4139841c-0c93-450a-b751-8ecf15595676︡{"stdout":"True\n"}︡
︠35fbf10c-5b62-489a-b006-9793a1a1ce35︠
dim=4
numVerts=8
Q = polytopes.cyclic_polytope(dim,numVerts)
type(Q.vertices_list()[0])

︡62308e0e-75d3-4174-9135-3600139a76c7︡{"stdout":"<type 'list'>\n"}︡
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
def line_shelling(P):
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
    for i in range(len(P.ieqs())):
        ShellingOrder.append([P.inequalities()[i],find_facet_vertices(P.ieqs()[i],homoVerts),find_intersection(startingPoint,directionVector,P.ieqs()[i])])
    ShellingOrder.sort(key=lambda x: x[2])
    iter=ShellingOrder[0][2]
    while iter<0:
        item=ShellingOrder.pop(0)
        ShellingOrder.append(item)
        iter=ShellingOrder[0][2]
    return [ShellingOrder[i][1] for i in range(len(ShellingOrder))]


line_shelling(Q)


︡0f6181e1-b887-441f-96f5-975ecc98dfe9︡{"stdout":"[[[0, 0, 0, 0], [1, 1, 1, 1], [2, 4, 8, 16], [3, 9, 27, 81]], [[0, 0, 0, 0], [2, 4, 8, 16], [3, 9, 27, 81], [7, 49, 343, 2401]], [[0, 0, 0, 0], [3, 9, 27, 81], [4, 16, 64, 256], [7, 49, 343, 2401]], [[0, 0, 0, 0], [1, 1, 1, 1], [2, 4, 8, 16], [7, 49, 343, 2401]], [[0, 0, 0, 0], [4, 16, 64, 256], [5, 25, 125, 625], [7, 49, 343, 2401]], [[0, 0, 0, 0], [5, 25, 125, 625], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[4, 16, 64, 256], [5, 25, 125, 625], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[3, 9, 27, 81], [4, 16, 64, 256], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[3, 9, 27, 81], [4, 16, 64, 256], [5, 25, 125, 625], [6, 36, 216, 1296]], [[0, 0, 0, 0], [1, 1, 1, 1], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[2, 4, 8, 16], [3, 9, 27, 81], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[2, 4, 8, 16], [3, 9, 27, 81], [5, 25, 125, 625], [6, 36, 216, 1296]], [[0, 0, 0, 0], [1, 1, 1, 1], [5, 25, 125, 625], [6, 36, 216, 1296]], [[2, 4, 8, 16], [3, 9, 27, 81], [4, 16, 64, 256], [5, 25, 125, 625]], [[1, 1, 1, 1], [2, 4, 8, 16], [6, 36, 216, 1296], [7, 49, 343, 2401]], [[1, 1, 1, 1], [2, 4, 8, 16], [5, 25, 125, 625], [6, 36, 216, 1296]], [[0, 0, 0, 0], [1, 1, 1, 1], [4, 16, 64, 256], [5, 25, 125, 625]], [[1, 1, 1, 1], [2, 4, 8, 16], [4, 16, 64, 256], [5, 25, 125, 625]], [[1, 1, 1, 1], [2, 4, 8, 16], [3, 9, 27, 81], [4, 16, 64, 256]], [[0, 0, 0, 0], [1, 1, 1, 1], [3, 9, 27, 81], [4, 16, 64, 256]]]\n"}︡
︠eb2dbb54-ef9f-49cd-b6eb-ae90c8107af2︠

