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


def inner(v,w):
	assert(v.ndim == w.ndim)

	return np.tensordot(v,w, v.ndim)

def laplace1(x, edges, faces):

	return d_0(ad_0(x, edges), edges) + ad_1(d_1(x, faces), faces)



def  crls(f,b,eps=1e-16):
	x = np.zeros(b.shape)
	r = b
	s = f(r)
	p = s
	norm = inner(s,s)
	normx_old = 500
	while abs(inner(x,x) - normx_old) > eps:
		normx_old = inner(x,x)
		y = f(p)
		alpha = norm / inner( f(y), f(y) )
		x = x + alpha*p
		r = r - alpha*y

		s = f(r)
		beta = inner(f(s), f(s)) / norm
		p = s + beta*p
	return x



def steep(f,b, eps=1e-10):
	x = np.zeros(b.shape)
	r = b
	s = f(r)

	while inner(s,s) > eps:
		print(inner(s,s))
		alpha = inner(s,s) / inner(f(s), f(s))
		x = x + alpha*s
		r = r - alpha*f(s)
		s = f(r)
	return x


def hodge(x, edges, faces):
	# Returns the hodge decompositon of x
	#1st = grad part
	#2nd = curl part
	#3rd = harmonic part
	def laplace(x):
		return laplace1(x, edges, faces)

	u = steep(laplace, x)	


	
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



def example_grid1(m):
	n = m*m
	edges = ngrid(m)
	faces = find_triangles(n, edges)

	v = np.random.rand(n)*10.0

	x = np.zeros((n,n))
	x = d_0(v, edges)
	return x, edges, faces




x, edges, faces = example_grid1(6)
print(x)





a = hodge(x, edges, faces)

print('Gradiente:')
print(a[0])
print("Rotacional:")
print(a[1])
print("Harmonica:")
print(a[2])



print('\n\n Revisao:')
print(d_0(a[0], edges) + ad_1(a[1], faces) + a[2] - x)