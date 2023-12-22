clear
clc

r_0 = [1.04703;0;0]
v_0 = [0;0.97728;0]

r_T = [1.0627; 0;0]
v_T = [0;0.97728;0]

L_0 = cross(r_0,v_0)
L_T = cross(r_T,v_T)
 
E_0 = v_0'*v_0/2 - 1/norm(r_0)
E_T = v_T'*v_T/2 - 1/norm(r_T)

dL = L_0-L_T
dE = E_0-E_T
lambda = [cross(-dL,v_0);cross(dL,r_0)]+dE*[r_0/norm(r_0)^3;v_0]