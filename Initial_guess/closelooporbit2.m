function dy=closelooporbit2(t,y)
global mu k1 k2 rd    
dy = y;
r=y(1:3);
v=y(4:6);

%k11 = k1+ mu/(norm(r))^3;
%k2=1;
dy(1:3)=v; %-[0.3;0.7];
k=2 ;
L_T =  [0, 0, 2.56612389857378]'*6378137*6378137/806.812;
L= cross(r,v);
dL = L-L_T;
A_T= [0,0,0]';
A= cross(v,L)-mu*r/norm(r);
dA= A-A_T;
G= - (cross(k*dL , r) + cross(L,dA)+ cross(cross(dA,v), r));
F_max =9.8e-2;
epsil= 0.00001;

%u=-J*((2*k2+3/2*k1)*eye(3)+k1*q*q'+4/k1*inv(J)*Sw'*J*J*inv(J))*z;
%u=-2*((1+1/k1)+mu^2/norm(r)^6/k1^3 +k2)*z+0*t;
normG = norm(G); 
%if normG<=epsil*F_max
   % dy(4:6)=-mu/norm(r)^3*(r)+1/epsil*G;
%else
    dy(4:6)=-mu/norm(r)^3*(r)+F_max*G/norm(G);
%$end
    
%dy(3:4)=u;
end 