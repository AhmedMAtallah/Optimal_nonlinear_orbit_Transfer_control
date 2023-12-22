clear
clc
r = sym('r',[3 1])
assumeAlso(r,'real')
v = sym('v',[3 1])
assumeAlso(v,'real')
L = cross(r,v)
A = cross(v,L)-r/(r(1)^2+r(2)^2+r(3)^2)^0.5
E = 1/(r(1)^2+r(2)^2+r(3)^2)^0.5

L_T = sym('L_T',[3 1])
assumeAlso(L_T,'real')
dL = sym('dL',[3 1])
assumeAlso(dL,'real')
dA = sym('dA',[3 1])
assumeAlso(dA,'real')
dE = sym('dE',[1 1])
assumeAlso(dE,'real')

E_T = sym('E_T',[1 1])
assumeAlso(E_T,'real')
A_T = sym('A_T',[3 1])
assumeAlso(A_T,'real')

V = 1/2 *(L-L_T)'*(L-L_T)
V= simplify(V)

L_der = dL'*[diff(L,r(1)) diff(L,r(2)) diff(L,r(3))]

A_der = dA'*[diff(A,r(1)) diff(A,r(2)) diff(A,r(3))]
E_der = dE' * [diff(E,r(1)) diff(E,r(2)) diff(E,r(3))]

A_f = r *(v'*v) - v *(v'*r)
A_f2 = r/(r'*r)%(r(1)^2+r(2)^2+r(3)^2)^0.5
%diff(A_f,r)
diff(A_f2,r)