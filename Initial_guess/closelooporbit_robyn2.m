function [dy,normG,dA_norm,dL_norm,r_norm]=closelooporbit_robyn2(t,y)
global mu k1 k2 rd    
F_max = 1.7e-4;
epsil= 0.0000001;
J_2= 0; %1.08264*10^(-3);
dy = y;
r=y(1:3);
v=y(4:6);
norm_r= norm(r);
mu =1;
kk=0; %J_2*mu/norm_r^3*3/2*(1/norm_r)^2*(2*[0;0;r(3)]);
%k11 = k1+ mu/(norm(r))^3;
%k2=1;
dy(1:3)=v; %-[0.3;0.7];
k=1 ;
r_f = [35786+6378.137,0,0 ]';
v_f = [0, 3.07460,0]';
r_f = [2.730128898256979e+04,0,0]';
v_f = [0,3.826879302585730,0]';
L_T = cross(r_f,v_f)/6378.127/6378.137*806.8111;
r_f_conical = r_f/6378.137;
v_f_conical = v_f/6378.137*806.8111;
E_T = 0.5*norm(v_f_conical)^2-1/norm(r_f_conical);
%L_T =  [0, 0, 2.56612389857378]';
%L_T= [ 0,   -2.416,    0.8794]';
%L_T = [    0,  -0.9694,    0.3528]';

L= cross(r,v);
dL = L-L_T;
dL_norm = norm(dL);
A_T= [0,0,0]';
%A_T= [4.0421,  0,      0]';
A= cross(v,L)-r/norm(r);
dA= A-A_T;
dA_norm = norm(dA);

E = 0.5/norm(v^2)-1/norm(r);
dE = norm(E-E_T);
G= - (cross(k*dL , r) + cross(L,dA)+ cross(cross(dA,v), r));
%dG = eye(3);
% for i =1:3
%     if abs(G(i))<eps*F_max
%         dG(i,i) = 1/F_max^2/0.001/.001;
%     else
%         dG(i,i) = 1/G(i)/G(i);
%     end
% end
%dG = diag(1./G)*diag(1./G);
%dG = eye(3);
%G= G/dot(G,G)*(1/dL_norm^2+dA_norm^2)+kk;

%G=(eye(3)+diag(L)*diag(L)+100*diag(dA)*diag(dA))*G;
%u=-J*((2*k2+3/2*k1)*eye(3)+k1*q*q'+4/k1*inv(J)*Sw'*J*J*inv(J))*z;
%u=-2*((1+1/k1)+mu^2/norm(r)^6/k1^3 +k2)*z+0*t;


normG = norm(G);
if normG<=epsil*F_max
    F=1/epsil*G;
else
%    G=(eye(3)+0.1*diag(L)*diag(L)+0.001*diag(dA)*diag(dA))*G;
  % G= G*(1+dL_norm^2+dA_norm^2/10)
    F=F_max*G/norm(G);
end

%F=G;
%G;
%F=[0;0;0];    
%F=G;
%a_J2 = J_2*mu/norm_r^3*3/2*(1/norm_r)^2*((5*r(3)^2/norm_r^2-1)*r-2*[0;0;r(3)]);
%aJ2= aJ2*3.94*10^5/6378/6378
dy(4:6)=-mu/norm(r)^3.*(r)+F;%+a_J2;
r_norm = norm(r);
%dy(3:4)=u;
end 