import numpy as np
from IPython.display import Image
from fatiando import mesher, gridder, utils
from fatiando.gravmag import prism, sphere, polyprism
from fatiando.utils import ang2vec, vec2ang
from fatiando.vis import mpl, myv
from fatiando.constants import CM, T2NT

def magnetic_data(x, y, z, model, alpha, eff_area = None, grains = None):
    '''
    Calculates the magnetic indution on the plane alpha
    located around the sample.
    
    input
    
    x, y, z: numpy arrays - Cartesian coordinates (in m) of 
             the points on which the magnetic field is calculated.
    model: list - geometrical elements of the Fatiando a Terra 
           class mesher - interpretation model.
    grains: None or list - if not None, is a list of geometrical elements 
            of the Fatiando a Terra class mesher - randomly magnetized grains.
    alpha: int - index of the plane on which the data are calculated.
    eff_area: None or tuple of floats - effective area of the simulated sensor.
              If None, calculates the field at the points x, y, z.
              If not None, eff_area = (lx, ly), where lx and ly
              are the side lengths (in microns) of the effective area of the sensor
              along the x and y axes.
    ns: None or tuple of ints - number of points on which the field is
        averaged within the effective area of the sensor.
        Is None if eff_area is none.
        If effe_area is not None, ns = (nsx,nsy), where nsx and nsy are
        the number of points on which the filed is averaged within the
        the effective area of the sensor along the x and y axes.
    
    output
    
    B: numpy array - magnetic data
    '''
    
    assert (alpha == 0) or (alpha == 1) or (alpha == 2) or (alpha == 3), \
           'alpha must be equal to 0, 1, 2 or 3'
           
    ns = 7
    ns2 = ns/2
    
    if eff_area is not None:
        
        dx = 0.000001*eff_area[0]/ns
        dy = 0.000001*eff_area[1]/ns
            
        xmin = x - ns2*dx
        
        B = np.zeros_like(x)
            
        if (alpha == 0) or (alpha == 2):
            ymin = y - ns2*dy
            if grains is None:
                for i in range(ns):
                    for j in range(ns):
                        B += prism.bz(xmin + i*dx, ymin + j*dy, z, model)
            else:
                for i in range(ns):
                    for j in range(ns):
                        B += prism.bz(xmin + i*dx, ymin + j*dy, z, model)
                        B += sphere.bz(xmin + i*dx, ymin + j*dy, z, grains)
        if (alpha == 1) or (alpha == 3):
            #zmin = z - ns2*dz
            zmin = z - ns2*dy
            if grains is None:
                for i in range(ns):
                    for j in range(ns):
                        #B += prism.by(xmin + i*dx, y, zmin + j*dz, model)
                        B += prism.by(xmin + i*dx, y, zmin + j*dy, model)
            else:
                for i in range(ns):
                    for j in range(ns):
                        #B += prism.by(xmin + i*dx, y, zmin + j*dz, model)
                        B += prism.by(xmin + i*dx, y, zmin + j*dy, model)
                        #B += sphere.by(xmin + i*dx, y, zmin + j*dz, grains)
                        B += sphere.by(xmin + i*dx, y, zmin + j*dy, grains)
        B /= ns*ns
        
    else:
        if (alpha == 0) or (alpha == 2):
            B = prism.bz(x, y, z, model)
            if grains is not None:
                B += sphere.bz(x, y, z, grains)
        if (alpha == 1) or (alpha == 3):
            B = prism.by(x, y, z, model)
            if grains is not None:
                B += sphere.by(x, y, z, grains)

    return B

def sample(Lx,Ly,Lz,N,m=None,inc=None,dec=None):
    '''
    Define o modelo interpretativo formado por prismas retangulares e justapostos ao longo da direcao x.
    
    input
    
    Lx: float - dimensao x (m) de todos os prismas.
    
    Ly: float - dimensao y (m) de todos os prismas.
    
    Lz: float - dimensao z (m) de todos os prismas.
    
    N: int - numero de prismas.
    
    m: list - intensidade de magnetizacao (A/m) de cada prisma.
    inc: list - inclinacao da magnetizacao (graus) de cada prisma.
    dec: list - declinacao da magnetizacao (graus) de cada prisma.
    
    output
    
    model: list de elementos geometricos da classe mesher da biblioteca Fatiando a Terra - modelo interpretativo.
    '''
    
    sizex = Lx
    sizey = Ly
    sizez = Lz

    L = N*sizex
    a = -0.5*L
    
    model = []

    if ((m is not None) & (inc is not None) & (dec is not None)):
        intensity = np.array(m)
        inclinacao = np.array(inc)
        declinacao = np.array(dec)
        mag = []
        for i in range(N):
            mag.append(ang2vec(intensity[i],inclinacao[i],declinacao[i]))
    
        for i in range(N):
            model.append(mesher.Prism(a+i*sizex, a+(i+1)*sizex, -0.5*sizey, 0.5*sizey, -0.5*sizez, 0.5*sizez, {'magnetization': mag[i]}))
    else:
        for i in range(N):
            model.append(mesher.Prism(a+i*sizex, a+(i+1)*sizex, -0.5*sizey, 0.5*sizey, -0.5*sizez, 0.5*sizez))
        
    return model



def sensitivity(N,x,y,z,model,alpha):
    '''
    Funcao que calcula a matriz de sensibilidade para a inversao
    
    input
    
    N: int - numero de prismas
    
    x: list - coordenadas x (in meters)
    
    y: list - coordenadas y (in meters)
    
    z: list - coordenadas z (in meters)
    model: list - Elementos geometricos da classe mesher
        da biblioteca Fatiando a Terra - modelo interpretativo.
    
    G: Matriz de sensibilidade com dimensao (n x m) em que n e o numero de observacos e m o numero de parametros 
    
    '''
    g = []
    for i in range(N):
        if (alpha == 0 or alpha == 2):
            g.append(prism.kernelxz(x,y,z,model[i]))
            g.append(prism.kernelyz(x,y,z,model[i]))
            g.append(prism.kernelzz(x,y,z,model[i]))
        if (alpha == 1 or alpha == 3):
            g.append(prism.kernelxy(x,y,z,model[i]))
            g.append(prism.kernelyy(x,y,z,model[i]))
            g.append(prism.kernelyz(x,y,z,model[i]))
    
    G = CM*np.transpose(g)*T2NT
    return G
    
    


def parameters_sph(N,coord):
    '''
    Funcao que transforma os parametros estimados em coordenadas Cartesianas 
    para intensidade, declinacao e inclinacao
    
    input
    
    N: int - numero de prismas
    
    coord: list - lista com os parametros estimados em coordenadas cartesianas
    
    return
    
    mag_sph: array - Matriz do tipo (N x 3) com cada linha contendo os valores 
        de intensidade, declinacao e inclinacao.
    '''
    
    mag_cart = []
    for i in range(N):
        mag_cart.append(np.array([coord[3*i],coord[3*i+1],coord[3*i+2]]))
        
        
    mag_spheric = []
    for i in range(N):
        aux = vec2ang(mag_cart[i])
        if (aux[2] > 180.):
            aux[2] -= 360.
        if (aux[2] <= -180.):
            aux[2] += 360.
        mag_spheric.append(aux)
    mag_sph = np.array(mag_spheric)
    return mag_sph



def proj_polar(N,mag,color='k',simbolo='o',size=10):
    
    '''
    Funcao que gera uma figura com a projecao polar de vetores
    
    input
    
    N: int - numero de prismas
    
    mag: array - Matriz do tipo (N x 3) com os valores de intensidade, declinacao e inclinacao
    
    proj: array - lista com valores das projecoes horizontais
    
    return
    
    Figura com as projecoes polares para os vetores
    
    '''
       
    projections = []
    for i in range(N):
        projections.append(np.sqrt((np.cos(np.deg2rad(mag[i,1]))*np.cos(np.deg2rad(mag[i,2])))**2 +(np.cos(np.deg2rad(mag[i,1]))*np.sin(np.deg2rad(mag[i,2])))**2))
        proj = np.array(projections)
    for i in range(N):
        if mag[i,1] >= 0.0:
            mpl.plot(np.deg2rad(mag[i,2]), proj[i], marker=simbolo, ms=size,mec= color, mew=3, mfc= color, fillstyle='full')
        else:
            mpl.plot(np.deg2rad(mag[i,2]), proj[i], marker=simbolo, ms=size, mec= color, mew=3, fillstyle='none') 



def dipolesrand(N, seed,n,deg_dec,deg_inc,std,mag,raio,Lx,Ly,Lz):
    '''
    Funcao que gera um conjunto de esferas 
    distribuidas aleatoriamente dentro da amostra.
    
    input
    
    N: int - numero de prismas
        
    s: int - valor da semente para a geracao de uma sequencia de numeros aleatorios
    
    n: int - numero de esferas magnetizadas
    
    deg_dec: float - valor da declinacao (graus) media do conjunto de esferas
    
    deg_inc: float - valor da inclinacao (graus) media do conjunto de esferas
    
    std: float - desvio padrao (graus) da declinacao e inclinacao do conjunto de esferas
    
    mag: float - intensidade de magnetizacao (A/m) em cada esfera
    
    raio: float - raio (m) de cada esfera 
    
    Lx,Ly,Lz: int - dimensoes do prisma (m)
    
    return
    
    modelrand: list de elementos geometricos da classe mesher da biblioteca Fatiando a Terra.
    '''
    np.random.seed(seed=seed)
    sizex = Lx
    sizey = Ly
    sizez = Lz
    L = N*sizex
    R = raio
    Coordx = np.random.uniform(-0.5*L+R,+0.5*L-R,n) 
    Coordy = np.random.uniform(-0.5*sizey+R,+0.5*sizey-R,n) 
    Coordz = np.random.uniform(-0.5*sizez+R,+0.5*sizez-R,n)
    Dec_rand = np.random.normal(deg_dec, std,n)
    Inc_rand = np.random.normal(deg_inc, std,n)
    
    magrand = []
    for i in range(n):
        magrand.append(ang2vec(mag,Inc_rand[i],Dec_rand[i]))
    
    modelrand = []
    for i in range(n):
        modelrand.append(mesher.Sphere(Coordx[i], Coordy[i], Coordz[i], R , {'magnetization': magrand[i]})) 

    return modelrand,Coordx,Coordy,Coordz

def L1_norm(A,d,n,std):
    At = A.T
    AtA = np.dot(At,A)
    Atd = np.dot(At,d)
    m0 = np.linalg.solve(AtA,Atd)
    for k in range(n):
        r = d - np.dot(A,m0)
        R = np.diag(1/np.abs(r))
        AtR = np.dot(At,R)
        AtRA= np.dot(AtR,A)
        AtRd = np.dot(AtR,d)
        m = np.linalg.solve(AtRA,AtRd)
        a = np.linalg.norm(m - m0)
        b = 1+np.linalg.norm(m)
        c = a/b
        tau = 2*std
        if c < tau:
            break
        else:
            m0=m
            continue
    return m

def L2_norm(G,d):
    Gt = G.T
    GtG = np.dot(Gt,G)
    Gtd = np.dot(Gt,d)
    m = np.linalg.solve(GtG,Gtd)
    return m


def residual(do,dp):
    r = do - dp
    r_mean = np.mean(r)
    r_std = np.std(r)
    r_norm = (r - r_mean)/r_std
    return r_norm, r_mean, r_std

def coordplane(h,L,Nx,Ny,area,alpha,theta):
    '''
    Funcao que calcula as coordenadas das medidas em um plano 
    a uma distancia h acima da amostra, sobre um grid regular
    de Nx x Ny pontos.
    
    input
    
    h: int - distancia (10**-6 m) da amostra ao sensor
    
    L: float - dimensao em (mm) da aresta do prisma paralela ao plano
    
    Nx: int - numero de observacoes ao longo do eixo x 
    
    Ny: int - numero de observacoes ao longo do eixo y
    
    area: list - lista com os valores da area que sera feita a medicao
    
    alpha: int - numero sobre qual o plano vao ser feitas as medidas
    
    theta: int - numero (em graus) com o valor de quanto o plano correspondente foi rotacionado 
    
    return
    
    x,y,z: list com as coordenadas de observacao acima
    da amostra
    
    '''
    shape = (Nx, Ny)
    areaxy = area
    dist = h*0.000001
    size = L
    voo = dist + 0.5*size
    ang = np.deg2rad(theta)
    cos = np.cos(ang)
    sin = np.sin(ang)
    
    if (alpha == 0):
        x, y, z = gridder.regular(areaxy, shape, -voo)
    if (alpha == 1):
        x, z, y = gridder.regular(areaxy, shape, voo)
    if (alpha == 2):
        x, y, z = gridder.regular(areaxy, shape, voo)
    if (alpha == 3):
        x, z, y = gridder.regular(areaxy, shape, -voo)
    
    x_rot = []
    y_rot = []
    z_rot = []
    if (alpha==0 or alpha ==2):
        x_rot.append(cos*x - sin*y)
        y_rot.append(sin*x + cos*y)
        z_rot.append(z)
    if (alpha == 1 or alpha == 3):
        x_rot.append(cos*x - sin*z)
        z_rot.append(sin*x + cos*z)
        y_rot.append(y)
    
    return x_rot[0], y_rot[0], z_rot[0]
