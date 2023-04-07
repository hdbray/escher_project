
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import time
from math import sqrt
from math import pi 
from math import sin 
from math import cos
start_time = time.time()

### generators of PSL(2,Z)

A=[[1,1],[0,1]]
B=[[0,-1],[1,0]]
Ainv=np.linalg.inv(A)

# this is the standard conformal map sending the disk to the upper half
# plane
J=np.array([[0-1j,1+0j],[1+0j,0-1j]])
Jinv=np.linalg.inv(J)


# pick two endpoints of a geodesic in the upper half plane 
# define inversion through that geodesic
def mobius_reflection(a,b):
    return np.matrix([[(a+b)/2,-a*b],[1,-(a+b)/2]])

# conjugate an element of PSL(2,R) into the isometry group of the Poincare
# disk

def conj_to_poincare(A):
    A=np.dot(A,J)
    A=np.dot(Jinv,A)
    return A 




        

### function for matrix product

def mp(list_of_matrices):
    X = np.identity(len(list_of_matrices[0]))
    for matrix in list_of_matrices:
        X = np.dot(X,matrix)
    return X.tolist()

### action of PSL(2,R) on upper half plane by mobius transformations
### this function is more general than needed for this picture but maybe it
### will be useful later?

#def mobius_transformation(M,v):
#    # M is a matrix which represents a mobius transformation
#    # v is a vector in R^2 which will be our model for the complex plane
#    # via the bijection [x,y] <-> x+iy
#    # output will be a vector
#    a=M[0][0]
#    b=M[0][1]
#    c=M[1][0]
#    d=M[1][1]
#    x=v[0]
#    y=v[1]
#    numerator_real_part = (a*x+b)*(c*x+d)+a*c*y*y
#    numerator_imaginary_part = (a*d-b*c-a*c*x)*y
#    denominator = (c*x+d)**2+(c*y)**2
##    if denominator!=0:
#    return [numerator_real_part/denominator,numerator_imaginary_part/denominator]
##    else:

def mobius_transformation(M,z):
    # M is a matrix which represents a mobius transformation
    # z is a complex number
    a=M[0][0]
    b=M[0][1]
    c=M[1][0]
    d=M[1][1]

    #output is a complex number
    w=(a*z+b)/(c*z+d)

    return round(w.real, 5) + round(w.imag, 5) * 1j




# this function seems unnecessary

def mobius_of_geodesic(M,v,w,height=4):
    # M is a matrix which represents a mobius transformation
    # v,w are distinct vectors that represent the endpoints of the initial
    # geodesic in the boundary of the Poincare disk
    # output will be a tuple of distinct vectors representing the image
    # geodesic


#    r_hor=r[0]
#    r_vert=r[1]
#    t_hor=t[0]
#    t_vert=t[1]
#    a=M[0][0]
#    b=M[0][1]
#    c=M[1][0]
#    d=M[1][1]

    endpoint1=mobius_transformation(M,v)
    endpoint2=mobius_transformation(M,w)
#    if r_hor==t_hor:
#        if c==0:
#            endpoint1=mobius_transformation(M,[r_hor,0])
#            endpoint2=mobius_transformation(M,[r_hor,height])
#        else:
#            if c*r_hor+d==0:
#                endpoint1=[a/c,height]#mobius_transformation(M,[r_hor,0])
#                endpoint2=[a/c,0.]
#            else:
#                endpoint1=mobius_transformation(M,[r_hor,0])
#                endpoint2=[a/c,0.]
#    else:
#        endpoint1=mobius_transformation(M,r)
#        endpoint2=mobius_transformation(M,t)
    return (endpoint1,endpoint2)


#print(mobius_of_geodesic(mp([A,B]),[1,3],[1,0]))

### generate the group

def group_elts(n,generators=[A,B]):
    list_these_matrices=[ [] for i in range(n+1)]
    list_these_matrices[0]=[(np.identity(len(generators[0]))).tolist()]
    with_duplicate_matrices=[ [] for i in range(n+1)]
    with_duplicate_matrices[0]=[(np.identity(len(generators[0]))).tolist()]
    flat_list = []
    for j in range(n):
        for i in list_these_matrices[j]:
#        range(len(list_these_matrices[j])):
            for M in generators:
                with_duplicate_matrices[j+1].append(mp([M,i]))
                with_duplicate_matrices[j+1].append(mp([i,M]))
        for k in with_duplicate_matrices[j+1]:
            if k not in list_these_matrices[j+1]:
                list_these_matrices[j+1].append(k)
        flat_list.extend(list_these_matrices[j+1])
    return flat_list



#### draw some lines
#
def draw_lines(list_of_vectors, default_color='black',default_linewidth=1):
    for i in range(len(list_of_vectors)-1):
        plt.plot([list_of_vectors[i][0],list_of_vectors[i+1][0]],[list_of_vectors[i][1],list_of_vectors[i+1][1]], color=default_color, linewidth=default_linewidth)
#
#
#draw_lines([[-3,0],[3,0]])
#for i in range(-2,3):
#    draw_lines([[i,0],[i,3]])


##### draw a geodesic orthocircle
def upper_half_orthocircle(a,b,default_color='black',default_linewidth=1):
    # a and b are real numbers
    center=(a+b)/2
    radius=abs(a-b)
    return patches.Arc((center,0), radius, radius, 0.0, 0.0, 180, color=default_color, linewidth=default_linewidth)


# find intersection of lines ax+b with cx+d assuming a\neq c
# return a list of two entries for the (x,y) coordinates of the
# intersection

def intersect_lines(a,b,c,d):
    return ((d-b)/(a-c),a*(d-b)/(a-c)+b)

# find the slope of the unit circle at a point (a,b) in the circle 
# assume both a, b are nonzero

def slope_of_circ(v):
    a=v[0]
    b=v[1]
    if abs(a)==1:
        pos_slope=1
    else: 
        pos_slope=a/sqrt(1-a**2)

    if b<0:
        return pos_slope
    else: 
        return -pos_slope

def eucl_dist(v,w):
    #v, w are lists with two coordinates
    v_hor=v[0]
    v_vert=v[1]
    w_hor=w[0]
    w_vert=w[1]

    return sqrt((v_hor-w_hor)**2+(v_vert-w_vert)**2)


def poincare_orthocircle(v,w,default_color='black',default_linewidth=1):
    # z,w are complex number

#    v_hor=v[0]
#    v_vert=v[1]
#    w_hor=w[0]
#    w_vert=w[1]
    v_hor=round(v.real,5)
    v_vert=round(v.imag,5)
    w_hor=round(w.real,5)
    w_vert=round(w.imag,5)

    v_vec=[v_hor,v_vert]
    w_vec=[w_hor,w_vert]

    if v_vert!=0:
        v_slope=slope_of_circ(v_vec)
        v_height=-v_hor*v_slope+v_vert
    if w_vert!=0:
        w_slope=slope_of_circ(w_vec)
        w_height=-w_hor*w_slope+w_vert

    if v_vert==0: 
        if w_vert==0:
            draw_circle=draw_lines([v_vec,w_vec],default_color,default_linewidth)
        elif w_hor==0: 
            center=(v_hor,w_vert)
            radius=eucl_dist(v_vec,center)
#            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
            draw_circle=ax.add_patch(draw_circle)
        else: 
            center=(v_hor, w_slope*v_hor+w_height)
            radius=eucl_dist(v_vec,center)
#            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
            draw_circle=ax.add_patch(draw_circle)
            
    elif w_vert==0:
        if v_hor==0:
            center=(w_hor,v_vert)
            radius = eucl_dist(v_vec,center)
#            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
            draw_circle=ax.add_patch(draw_circle)
        else: 
            center=(w_hor, v_slope*w_hor+v_height)
            radius=eucl_dist(v_vec,center)
#            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
            draw_circle=ax.add_patch(draw_circle)
    else: 
        if v_slope==w_slope:
            draw_circle=draw_lines([v_vec,w_vec],default_color,default_linewidth)
        else: 
            center=intersect_lines(v_slope,v_height,w_slope,w_height)
            radius=eucl_dist(center,w_vec)
#            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
            draw_circle=ax.add_patch(draw_circle)

#            print(v_slope,v_height,w_slope,w_height)
#            print(center)
#            print(radius)

    
    return draw_circle



#geodesic= patches.Arc((0,0), 1, 1, 0.0, 0.0, 180)
#patches.Arc((xcenter,ycenter), width, height, angle=0.0, theta1=0.0,
#theta2=360.0)

def draw_geodesic(vectors, default_linewidth=1):
    v=vectors[0]
    w=vectors[1]
    v_hor=v[0]
    v_vert=abs(v[1])
    w_hor=w[0]
    w_vert=abs(w[1])
    if v_hor==w_hor:
        draw_lines([[v_hor,v_vert],[w_hor,w_vert]],'black',default_linewidth)
    else:
        ax.add_patch(orthocircle(v_hor,w_hor,'black',default_linewidth))

# starting to draw

fig = plt.figure()
ax = fig.add_subplot()


# the boundary of the disk
ax.add_patch(patches.Circle( (0,0),1,fill=False))

# init_geodesics is a list of generating geodesics in the tiling
# list of vectors with distinct complex entries in the boundary of the
# unit circle 
# poincare_orthocircle will draw a geodesic between those endpoints for
# each vector in the list

init_geodesics=[[complex(sqrt(2)/2,sqrt(2)/2),complex(-sqrt(2)/2,-sqrt(2)/2)],[-1,1],[-1j,1j],[complex(-sqrt(2)/2,sqrt(2)/2),complex(sqrt(2)/2,-sqrt(2)/2)],[complex(-sqrt(3)/2,1/2),complex(-1/2,sqrt(3)/2)]]

#for c in init_geodesics:
#    poincare_orthocircle(c[0],c[1])


# this is inversion through the unit circle in the upper half plane
A=[[0,-1],[1,0]]

# this will be inversion through the real axis in the Poincare disk
# note that this function preserves any line passing through the origin
# as inversion in the upper half plane does for any geodesic through i
# which intersects the unit disk orthogonally
#A=[[-1,0],[0,1]]
##A=conj_to_poincare(A) is not necessary here
#
#A=[[1,1],[0,1]]
#B=[[0,-1],[1,0]]
#

A=mobius_reflection(-2.414,.414)
B=mobius_reflection(-.613,3.573)

A=conj_to_poincare(A)
B=conj_to_poincare(B)

#print(A,B)

n=2
#print(group_elts(n,[A,B]))
for G in group_elts(n,[A,B]):
    for c in init_geodesics:
        v=mobius_transformation(G,c[0])
        w=mobius_transformation(G,c[1])
        poincare_orthocircle(v,w)

#G=group_elts(n,[A,B])[2]
#c=init_geodesics[3]
#v=mobius_transformation(G,c[0])
#w=mobius_transformation(G,c[1])
#print(v,w)
#poincare_orthocircle(v,w)

#
#for c in init_geodesics: v=mobius_transformation(G,c[0])
#    w=mobius_transformation(G,c[1])
#    print (v,w)
##print(G)

#for c in init_geodesics:
#    print(c)
#    v=mobius_transformation(G,c[0])
#    w=mobius_transformation(G,c[1])
#    poincare_orthocircle(v,w)
#c=init_geodesics[4]
#v=mobius_transformation(G,c[0])
#w=mobius_transformation(G,c[1])
#poincare_orthocircle(v,w)




axes = plt.gca()
axes.set_xlim([-1.1,1.1])
axes.set_ylim([-1.1,1.1])

ax.set_aspect("equal")


plt.show()
