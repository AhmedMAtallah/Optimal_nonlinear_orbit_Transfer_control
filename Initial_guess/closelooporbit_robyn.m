function [dy,normG,dA_norm,dL_norm,dE_norm,r_norm]=closelooporbit_robyn(...
    t,y,...
   F_max,r_f_norm,v_f_norm,k) %...
   % )
%F_max =0.01;
mu =1;
epsil= 0.0001;
%k=2;

%F_max = 0.01;
dy = y;
r=y(1:3);
v=y(4:6);
dy(1:3)=v;

r_f = [r_f_norm,0,0]';
v_f = [0,v_f_norm,0]';

%r_f = [1.152320630706853, 0, 0]';
%v_f = [0, 0.550215519767527, 1.098755003610472]';
L_T = cross(r_f,v_f);
%L_T = [0, 0, 2.56612389857378]';
%r_f_conical = r_f;
%v_f_conical = v_f;
E_T = 0.5*norm(v_f)^2-1/norm(r_f);
A_T = cross(v_f,L_T)-r_f/norm(r_f);
%A_T= [0,0,0.001]';

L= cross(r,v);
dL = L-L_T;
dL_norm = norm(dL);

A= cross(v,L)-r/norm(r);
dA= A-A_T;
dA_norm = norm(dA);

E = 0.5*norm(v)^2-1/norm(r);
dE = (E-E_T);
G= - diag([1.0,1,1.0])*(cross(k*dL , r) + cross(L,dA)+ cross(cross(dA,v), r));

%G= - (cross(k*dL , r) + 2*dE*v);
dE_norm = norm(dE);

normG = norm(G);
if normG <= epsil*F_max
   F=1/epsil*G;
else
    F=F_max*G/norm(G);
end

dy(4:6)=-mu/norm(r)^3.*(r)+F;%+a_J2;
r_norm = norm(r);

end 