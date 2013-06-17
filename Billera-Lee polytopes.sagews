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





#This method computes the boundary complex of a pseudomanifold (we will only use it for d-balls).
#Since we intent to work with complexes that are guaranteed to be pseudomanifolds
#this method does NOT check wether the imput facets are the facets of a pseudomanifold.
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
s=SimplicialComplex(boundary(ball_constructor([1,4,5,5,0],8)))
#b.h_vector()
#s.h_vector()
#s.g_vector()
print "hello"
︡caea3fff-aa12-47eb-8665-6e579b054731︡{"stdout":"hello\n"}︡
︠554e8b08-eb31-4786-b638-d83f137fc0e0︠
︠05736660-bd6f-4070-a8c0-b8be1102b23b︠



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
    print desiredFacetCenters

    #We want a shelling that only captures the facets with centers in desiredFacetCenters.
    #We need to choose a starting point and direction that accomplishes this...
    polytopeCenter=list(Q.center())
    averageCenter=[sum(i)/len(desiredFacetCenters) for i in zip(*desiredFacetCenters)]
    directionList=[[d[j]/d[0]-averageCenter[j]/d[0] for j in range(dimension+1)] for d in desiredFacetCenters]
    perturbationVector=[ZZ.random_element(90,110)/100001  for i in range(dimension+1)]
    preDirection = [sum(i) for i in zip(*directionList)]
    #preDirection = [-sum(i) for i in zip(*desiredFacetInequalities)]
    direction = [preDirection[i]-perturbationVector[i] for i in range(dimension+1)]
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
    shelling_order = line_shelling(Q, stopPoint, averageCenter, direction )
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



BilleraLeePolytope([1,4,6,6,0],8)
︡a2ad119b-8f5d-48d4-931d-140ee443a66d︡{"stdout":"[[7411742281/9, 5425572705273042049, 339517777314472371236468999881/9, 263104438465821872126630483794358196481, 16517988560338364532505825903765048219418595395401/9, 12802822803438995170325754166326162573979473264164811869569, 803784446590736159342877641324229568279192559228641051416219539634121/9, 623000591267204855323161712806746338446216178427441280215829282434667478479361, 39113109087239253520113099166094486365063175811912461098426842157477571739815353280206921/9], [125589280105/9, 1567971803176513663169, 1668050770854068846876039315588137/9, 21974745801953787195569671205018578998727041, 23453181683304715700365066274663297975440290768849593065/9, 309029018812781652803812217048264690193709392205457332287717630529, 329823843192082047759041155569012066103824863475766835442386482735615466074537/9, 4345901010279603888702564395021087504660603547164263662915985163205631787657229417767681, 4638340549851216995133314885538338341289635211981673161433088501560677140615277717132448265575608425/9], [2134607423113/9, 453143832409365053166849, 8195133437136948310788972424603689865/9, 1835352744124979110239788190080924766540751628801, 33300164185311983708539692240486430694682630691448643182207113/9, 7459209264595815226485530441453038671343748379369540242674670573387136769, 135339478139199031584567569569743163352689517011830426123250600089800046352781381082505/9, 30315951310307364317149449414260705129811449112927013358076555733937133316132359596731117902156801, 550050956276783192744016995696431828751612707927061944358642643275665901685775670348745473094076954602955659913/9], [36287915854249/9, 130958567547597852969730369, 40262690576653757958972308513345047479529/9, 153290496542062380263187228037430283994280104859525761, 47281471219664517252456030171950676968684417133573957424925566558249/9, 180047178309620747140545118531642770596801430755135013168088313541359991451067329, 55535021864151439903296543776083588642133039618303164665217280897714594666668987937981754628969/9, 211476722933870296632453155660582205323099777693756584461753837770867630430079160758061164510512734081598721, 65229374870007912209572834021317438211576101167406312353373796528195491475070226826147213096006320951853936315013400737449/9], [616894159183561/9, 37847026021237070860856587649, 197810598803099912783339017813055485386098761/9, 12802975561689592061961657322792828430852836667960521521921, 67132927881539202472520461631843738799060292282446086373243802323217465801/9, 4345901189703774147936460496170706058101055121724175014125199328449919171447961754275969, 22788167176761888120957952098564548900490741407599139327083490609014550718666192651654625658142096231241/9, 1475210323604241072662712342683237589174803083432579003919087435479681632660867192350940217027610637023260698748925441, 7735413051061013722521641622046272227866815240957157290850273361573356718190707609261445314499121392983636169505297208398654157808841/9], [125998198921/9, 1567990511599909829889, 1668050839945999897461720686599561/9, 21974745805103908581436723369583305092467201, 23453181683316349250030118679614717721276649662757886601/9, 309029018812782183214404643665824840554001055435594040426097157889, 329823843192082049717872963342354312161497103207178893799362057832302638376841/9, 4345901010279603888791873777757965765413861554190164910611377891726040972248779113190401, 4638340549851216995133644709380726638892762325332507472787770407095071358112454227746792998420153481/9], [2135016341929/9, 453143851117788449333569, 8195133437206040241839558105974701289/9, 1835352744124982260361174057133089331266845368961, 33300164185311983720173241905538835646102376527807537090500649/9, 7459209264595815226486060852045465288903898739661203472811378711766664129, 135339478139199031584567571528574971126031763069502665854662658446775621449468553384809/9, 30315951310307364317149449414350014512548327373680271365102456981632526044652768781322667597579521, 550050956276783192744016995696432158575455096224665071471993477587020583591310064566242649604691299335800204969/9], [2141967961801/9, 453149257852149941515649, 8195133776654697493367010681753827401/9, 1835352744388083548630176120969699819342117272321, 33300164185328501697100014218075943632166351557835485025645001/9, 7459209264595828029308333880404260029875813621003391144934376101606051969, 135339478139199032388352016160479313148176471389828166297062576225362082428118331587401/9, 30315951310307364317772450005527909985109005181184032273444776909026063334432014114888990278039041, 550050956276783192744056108805519068005132815463705930618832540163253945841576189053185972379842667936880153801/9], [36288324773065/9, 130958567566306276365897089, 40262690576653827050903359099026418490953/9, 153290496542062380266337349423297336158844830953265921, 47281471219664517252456041805500342021089368553319793783819474851785/9, 180047178309620747140545118532173181189228048315285373459751543678068129830594689, 55535021864151439903296543776083590600964847391645410722889520629126653023644563034668926931273/9, 211476722933870296632453155660582205323189087076493462722507095777893531677774553486581573695104283777021441, 65229374870007912209572834021317438211576101497230154741671399655308842309381581508052747490223818128364550659746245282505/9], [36295276392937/9, 130958572973040637858079169, 40262690576993275708154886551602197617065/9, 153290496542062643367625618425361172769332906225169281, 47281471219664517268974018732272654558197354617294823811767409996137/9, 180047178309620747140557921354446209548022789287200254801939215801065519669982529, 55535021864151439903296544579868035232869189413790119043215021071526570802231024013318705133865/9, 211476722933870296632453155661205205914366982549054140530010856686235851605168090776360819028670606457480961, 65229374870007912209572834021317477324685188406659832460910440514455681371957814870303013614710761451139702028347325231337/9], [36413453930761/9, 130960135519271109098700289, 40262692244704528784909261056405044205321/9, 153290496564037126065140939799058717343867690865699841, 47281471219687970434139334887634655091966192011505215883117664193801/9, 180047178309621056169563931313295574365069511797724122348757441785133642575743489, 55535021864151440233120386968165636401174185639924543822810316268100847015976455332714631574281/9, 211476722933870296636799056670861809211802342088752066321864713894953899089656634659350621533893168396769281, 65229374870007912209577472361867289428571234482291850686152564055181297981943436046343475865009219251285094661259620632841/9], [125999613865/9, 1567990511823906470209, 1668050839946002759884653608538857/9, 21974745805103908581888304636923855124152961, 23453181683316349250030124450242130066816605550818061865/9, 309029018812782183214404643666735225634647138934868434755104627649, 329823843192082049717872963342354323795046768265354403539522383561403235834217/9, 4345901010279603888791873777757965765415696884959807189056099486410287302484816137995521, 4638340549851216995133644709380726638892762348785672638115559318758995591206511844090088646425220265/9], [2135017756873/9, 453143851118012445973889, 8195133437206040244701981038896640585/9, 1835352744124982260361174508714356671816877054721, 33300164185311983720173241905544606273514722067763425150675913/9, 7459209264595815226486060852045465289814283820307286972085773040774133889, 135339478139199031584567571528574971126031774703052330912838168186935947178569150842185/9, 30315951310307364317149449414350014512548327373682106695872099260077247639337015111558704622384641, 550050956276783192744016995696432158575455096224665071495446642752348372502973988799336707221034594983805271753/9], [2141969376745/9, 453149257852373938155969, 8195133776654697496229433614675766697/9, 1835352744388083548630176572550967159892148958081, 33300164185328501697100014218081714259578697097791373085820265/9, 7459209264595828029308333880404260030786198701649474644208770430613521729, 135339478139199032388352016160479313148176483023377831355238085965522408157218929044777/9, 30315951310307364317772450005527909985109005181185867604214419187470784929116260445125027302844161, 550050956276783192744056108805519068005132815463705930642285705328581734753240113286280029996185963584885220585/9], [2141993430793/9, 453149257917108967208449, 8195133776654711559313303060163527945/9, 1835352744388083548667893091580517239088575319041, 33300164185328501697100022411547440070280022245144829355637513/9, 7459209264595828029308333880426234513486864106692971787258396933460935169, 135339478139199032388352016160479317921871814862944037823114048832449304801385770946825/9, 30315951310307364317772450005527909985121808003458895984983379758762802312737239007260844717742081, 550050956276783192744056108805519068005132818244966984731920123307069224477763658077800983055768407796114196233/9], [36288326188009/9, 130958567566306500362537409, 40262690576653827050906221521959340430249/9, 153290496542062380266337349423748917426185380984951681, 47281471219664517252456041805500342026859995965665333739707535027049/9, 180047178309620747140545118532173181189228049225670454105835042952462458838064449, 55535021864151439903296543776083590600964847391645422356439185687302162763804888763769524388649/9, 211476722933870296632453155660582205323189087076493462722508931108663173956219275081265820025340320801826561, 65229374870007912209572834021317438211576101497230154741671399655332295474546909296964411414456912185980893955394250349289/9], [36295277807881/9, 130958572973040861854719489, 40262690576993275708157748974535119556361/9, 153290496542062643367625618425812754036673456256855041, 47281471219664517268974018732272654563967982029640363767655470171401/9, 180047178309620747140557921354446209548022790197585335448022715075459848677452289, 55535021864151439903296544579868035232869189413790130676764686129702080542391349742419302591241/9, 211476722933870296632453155661205205914366982549054140530012692017005493883612812371045065358906643482286081, 65229374870007912209572834021317477324685188406659832460910440514479134537123142659214677538943855508756045323995330298121/9]]"}︡{"stdout":"\n"}︡{"tex":{"tex":"\\left[1, 5, 11, 17, 17, 17, 11, 5, 1\\right]","display":true}}︡{"tex":{"tex":"\\left[1, 4, 6, 6, 0\\right]","display":true}}︡
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


