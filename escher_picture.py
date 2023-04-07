
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import time
from math import sqrt
start_time = time.time()

### generators of PSL(2,Z)

A=[[1,1],[0,1]]
B=[[0,1],[1,0]]
Ainv=np.linalg.inv(A)
# note that B is its own inverse so we do not need to introduce it


def mobius_reflection(a,b):
    return np.matrix([[(a+b)/2,-a*b],[1,-(a+b)/2]])

A=mobius_reflection(-2.414,.414)
B=mobius_reflection(-.613,3.573)

#print(A)
# this is the standard conformal map sending the disk to the upper half
# plane
#a= complex(0,-1)
#a=0-1j
#print(a)
J=np.array([[0-1j,1+0j],[1+0j,0-1j]])
Jinv=np.linalg.inv(J)

A=np.dot(A,J)
A=np.dot(Jinv,A)
B=np.dot(B,J)
B=np.dot(Jinv,B)
#print(J)
#print(np.dot(J,A))

#print(A)
        

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

    return (a*z+b)/(c*z+d)




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

#n=3
#for M in group_elts(n):

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
    v_hor=v.real
    v_vert=v.imag
    w_hor=w.real
    w_vert=w.imag

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
#            print(v,w)
        elif w_hor==0: 
            center=(v_hor,w_vert)
            radius=eucl_dist(v_vec,center)
            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
        else: 
            center=(v_hor, w_slope*v_hor+w_height)
            radius=eucl_dist(v_vec,center)
            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            
    elif w_vert==0:
        if v_hor==0:
            center=(w_hor,v_vert)
            radius = eucl_dist(v_vec,center)
            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
        else: 
            center=(w_hor, v_slope*w_hor+v_height)
            radius=eucl_dist(v_vec,center)
            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
            draw_circle=patches.Circle(center,radius,color=default_color,linewidth=default_linewidth,fill=False)
    else: 
        if v_slope==w_slope:
            draw_circle=draw_lines([v_vec,w_vec],default_color,default_linewidth)
        else: 
            center=intersect_lines(v_slope,v_height,w_slope,w_height)
            radius=eucl_dist(center,w_vec)
            draw_circle= patches.Arc(center, radius, radius, 0.0, 0.0, 360, color=default_color, linewidth=default_linewidth)
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


fig = plt.figure()
ax = fig.add_subplot()

#ax.add_patch(orthocircle(-3,1))
#ax.add_patch(orthocircle(-1,1))

lw=.5


#draw_lines([[-3,0],[3,0]],'black',lw)
#draw_lines([[0,-3],[0,3]],'black',lw)

ax.add_patch(patches.Circle( (0,0),1,fill=False))

#A=mobius_reflection(-2.414,.414)
#B=mobius_reflection(-.613,3.573)
#A=[[0,1],[1,0]]



#v=(sqrt(2)/2,sqrt(2)/2)


#z= complex(sqrt(2)/2,sqrt(2)/2)
#print(z)
#
#print(mobius_transformation(A,z))
#
#print(mobius_transformation(A,-1j))
#print(mobius_transformation(A,-1))

#n=2
#for G in group_elts(n,[A,B]):
    # apply G to v,w

#w=1

#Az=mobius_transformation(A,z)
#Aw=mobius_transformation(A,w)




A=[[0,1],[1,0]]
v=complex(sqrt(2)/2,sqrt(2)/2) 
w=1j

#print()
#print(v.real)
#

poincare_orthocircle(v,w)


poincare_orthocircle(mobius_transformation(A,v),mobius_transformation(A,w))

#w=complex(1/2,-sqrt(3)/2)
#
##poincare_orthocircle(v,w)



#poincare_orthocircle([1,0],[-1,0])
#poincare_orthocircle([0,1],[0,-1])

axes = plt.gca()
axes.set_xlim([-1.1,1.1])
axes.set_ylim([-1.1,1.1])

ax.set_aspect("equal")


plt.show()
