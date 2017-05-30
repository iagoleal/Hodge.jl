import numpy as np
from itertools import permutations

def d_0(v, edges):
	x = np.zeros([v.shape[-1]]*2)
	for edge in edges:
		for i,j in permutations(edge):
			x[i,j] = v[j] - v[i]
	return x

#Formula for adjoint 
# (d_0^* X)[i] = -2\sum_j X[i][j]
def ad_0(v, edges):
	n = v.shape[-1]
	x = np.zeros([n])

	for edge in edges:
		for i,j in permutations(edge):
			x[i] += v[i,j]
	return -2*x


def d_1(v, faces):
	x = np.zeros([v.shape[-1]]*3)
	for face in faces:
		for i,j, k in permutations(face):
			x[i, j, k] = v[j,k] - v[i,k] + v[i,j]
	return x


def ad_1(v, faces):
	n = v.shape[-1]
	x = np.zeros([n]*2)
	for face in faces:
		for i,j,k in permutations(face):
			x[i,j] += v[i,j,k]
	return 3*x

#Canonical inner product
def inner(v,w):
	assert(v.ndim == w.ndim)

	return np.tensordot(v,w, v.ndim)

#Laplace operator for 1-forms. Receives the 1 form, the edges and the faces
def laplace1(x, edges, faces):

	return d_0(ad_0(x, edges), edges) + ad_1(d_1(x, faces), faces)

# Symmetric Conjugate gradient for least squares
# Applies cg to A^2 = A b (least squares on the self-adjoint case)
def scgls(f, b, eps=1e-10):
	# f MUST be a function
	x = np.zeros(b.shape)
	r = f(b)
	p = f(b)
	#i = 1
	while (inner(r,r) > eps):
		#print(i, inner(r,r))
		#i+= 1
		#Applies twice for least squares
		y = f(f(p))
		norm_old = inner(r, r)
		a = norm_old / inner(p, y)
		x = x + a*p
		
		r = r - a*y

		beta = inner(r, r) / norm_old

		p = r + beta*p
	return x


#Discrete Hodge decomposition
#How it works:
#Receives an 1-form x and the edges and faces of the graph as input (necessary for the derivatives)
#returns 
#a 0-form
#b 1-form
#c 2-form
#where x = d(a) + d*(b) + c
#and laplace(c) = 0
def hodge(x, edges, faces):

	# Just an alias for calling the graph
	def laplace(x):
		return laplace1(x, edges, faces)



	#Uses symmetric conjugate gradient for least square problem to solve laplace(z) = x
	#The idea is to write x = laplace(z) + c = d(  d*(a) ) + d*( d(b)  ) + c
	#Which gives the unique solution to the problem
	u = scgls(laplace, x)	


	
	return ad_0(u, edges), d_1(u, faces), x - laplace(u)





def adjacency(n, edges):
	adj = np.full((n,n), False, dtype=bool)
	for i,j in edges:
			adj[i,j] = True
			adj[j,i] = True
	return adj


def find_triangles(n, edges):
	adj = adjacency(n, edges)
	faces = []
	for i in range(0,n):
		for j in range(i+1, n):
			if adj[i,j]:
				for k in range(j+1,n):
					if adj[j,k] and adj[i,k]:
						faces.append((i,j,k))

	return faces



# Return edges of a nxn grid with horizontal, vertical and principal diagonal components
def ngrid(n):
	edges = []
	for k in range(0, n*n-1):
		# horizontal
		if (k+1)%n != 0:
			edges.append((k, k+1))
		if k < n*(n-1):
		# vertical
			edges.append((k, k+n))
		# diagonal
		if k < n*(n-1) and (k+1)%n != 0:
			edges.append((k, k+n+1))
	return edges



def example_grid1(m):
	n = m*m
	edges = ngrid(m)
	faces = find_triangles(n, edges)

	v = np.random.rand(n)*10.0

	x = np.zeros((n,n))
	x = d_0(v, edges)
	return x, edges, faces



def example_harmonic():
	n = 6
	edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0), (0,2), (4,0)]
	faces = find_triangles(n, edges)

	x = np.zeros((n,n))
	x[0,1] = 1
	x[0,2] = 2
	x[0,4] =-2
	x[0,5] =-1
	x[1,2] = 1
	x[2,3] = 3
	x[3,4] = 3
	x[4,5] = 1

	x = (x - x.T)

	return x, edges, faces





x, edges, faces = example_grid1(20)
print(x)


a = hodge(x, edges, faces)

print('Gradiente:')
print(a[0])
print("Rotacional:")
print(a[1])
print("Harmonica:")
print(a[2])
print("Teste se e mesmo Harmonica:")
print(laplace1(a[2], edges, faces))

print('\n\n Revisao:')
print(d_0(a[0], edges) + ad_1(a[1], faces) + a[2] - x)