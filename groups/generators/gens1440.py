#!/usr/bin/env python3

# Generates 360x4 = 1440 subgroup of su4 (Hanany-He group III, Lift(A6))
# Generators: F1, F2, F3, F4
# Derived from gens720x4.py by stripping fpp (F''); IX = (III) + F''.
# Note: Hanany-He §3.1.1 prints III's generators as just F1,F2,F3 — that's a typo
# (those produce group I, order 240). F4 is the corrected addition.

from numpy import *
from numpy.linalg import norm

def equal(A, B):
	return norm(A-B) < 1e-5

def matrix_key(M):
    return tuple(round(M.real, 5).flatten()) + tuple(round(M.imag, 5).flatten())

def extend_unique(Ms):
    seen_keys = set(matrix_key(M) for M in Ms)
    result = Ms.copy()

    for A in Ms:
        for B in Ms:
            C = A * B
            key = matrix_key(C)
            if key not in seen_keys:
                seen_keys.add(key)
                result.append(C)

    return result

def generate(Ms):
	l = []
	for A in Ms:
		for B in Ms:
			l.append(A*B)
	return l

w=exp(2.0*pi*1j/3.0)
b=exp(2.0*pi*1j/7.0)
p=b+b*b+b*b*b*b
q=b*b*b+b*b*b*b*b+b*b*b*b*b*b
s=b*b+b*b*b*b*b
t=b*b*b+b*b*b*b
u=b+b*b*b*b*b*b

f1= matrix([[1,0,0,0],[0,1,0,0],[0,0,w,0],[0,0,0,w*w]])
f2= 1.0/sqrt(3.0)*matrix([[1,0,0,sqrt(2.0)],[0,-1.0,sqrt(2.0),0],[0,sqrt(2.0),1.0,0],[sqrt(2.0),0,0,-1.0]])
f3= matrix([[sqrt(3.0)/2.0,0.5,0,0],[0.5,-sqrt(3.0)/2.0,0,0],[0,0,0,1],[0,0,1,0]])
f4= matrix([[0,1.0,0,0],[1.0,0,0,0],[0,0,0,-1.0],[0,0,-1.0,0]])

gen = [f1,f2,f3,f4]
els = extend_unique(gen)

l = len(gen)

while True:
    l = len(els)
    print("Elements:", l)
    els = extend_unique(els)
    if len(els) == l:
        break

print(l)

trs = []
for i in range(l):
	trs.append(-els[i].trace()[0,0].real)
print(' '.join(str(x) for x in trs))

trs = []
for i in range(l):
        trs.append(-els[i].trace()[0,0].imag)
print(' '.join(str(x) for x in trs))

for A in els:
    for B in els:
        M = A*B
        for i in range(l):
            if equal(M, els[i]):
                print(i, end=' ')
                break
    print('')
