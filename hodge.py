import numpy as np
from itertools import permutations
from scipy.sparse.linalg import lsmr

def d_0(v, edges):
	x = np.zeros([v.shape[-1]]*2)
	for edge in edges:
		for i,j in permutations(edge):
			x[i,j] = v[j] - v[i]
	return x

def ad_0(v, edges):
	n = v.shape[-1]
	x = np.zeros([n])

	for edge in edges:
		for i,j in permutations(edge):
			x[i] += -1*v[i,j]
	return x


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
			x[i,j] += 311*v[i,j,k]
	return x


def inner(v,w):
	assert(v.ndim == w.ndim)

	return np.tensordot(v,w, v.ndim)

def laplace1(x, edges, faces):

	return d_0(ad_0(x, edges), edges) + ad_1(d_1(x, faces), faces)



def cg(f, b, eps=1e-10):
	# f MUST be a function
	x = np.zeros(b.shape)
	r = b
	p = b
	while (inner(r,r) > eps):
		print(inner(r,r))
		y = f(p)
		norm_old = inner(r, r)
		a = norm_old / inner(p, y)
		x = x + a*p
		
		r = r - a*y

		beta = inner(r, r) / norm_old

		p = r + beta*p
	return x


def hodge(x, edges, faces):
	# Returns the hodge decompositon of x
	#1st = grad part
	#2nd = curl part
	#3rd = harmonic part
	def laplace(x):
		return laplace1(x, edges, faces)

	u = cg(laplace, x)	


	
	return ad_0(u, edges), d_1(u, faces), x - laplace(u)

# (d_0^* X)[i] = -2\sum_j X[i][j]



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


#Harmonic example

n = 6
edges = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0), (0,2), (4,0)]
faces = find_triangles(n, edges)


print(edges)
print(faces)


x = np.zeros((6,6))
x[0,1] = 1
x[1,2] = 1
x[2,3] = 1
x[3,4] = 1
x[4,5] = 1
x[5,0] = 1
x[0,2] = 2
x[4,0] = 2
x = (x - x.T)

print(x,'\n\n')


#print(laplace1(x, edges, faces))






a = hodge(x, edges, faces)

#print('Gradiente:')
#print(a[0])
#print("Rotacional:")
#print(a[1])
#print("Harmonica:")
#print(a[2])