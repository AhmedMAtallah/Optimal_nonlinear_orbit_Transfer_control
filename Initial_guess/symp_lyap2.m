clear
clc

syms r [3 1] matrix
syms v [3 1] matrix

D = r*(r(1)^2+r(2)^2+r(3)^2)

diff(D,r)
syms A [3 4] matrix
alpha = cross(r,v)
diff(alpha,r)

