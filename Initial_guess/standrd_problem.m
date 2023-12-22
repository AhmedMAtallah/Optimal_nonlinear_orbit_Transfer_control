function standrd_problem

%% ---Constants------
mu = 1; %---km^3 S^-2
Re = 1; %---km
canonical_distance = Re;
canonical_time = 1;
canonical_accel = 1;
canonical_vel = canonical_distance/canonical_time;

%% ---- Inputs -----
mu_earth = 1;
r_0 = [1, 0, 0]';
v_0 = [0, sqrt(mu_earth/norm(r_0)), 0]';
r_f_norm = 1.4; %1.5138; %35780+6378.137; %
v_f_norm = sqrt(mu_earth/r_f_norm); %3.738702503601947;
span =100;%631.23; %2*pi*14
F_max = 0.8405;

x0 = [r_0;v_0];
r_f_norm = r_f_norm/canonical_distance;
v_f_norm = v_f_norm/canonical_vel;
F_max = F_max/canonical_accel;
%F_max =0.01;
span = span/canonical_time;

options = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',@events);
%options = odeset('RelTol',1e-11,'AbsTol',1e-11);
[Tk2e01,Yk2e01]=ode45(@closelooporbit_standard,[0 span],x0,...
    options,F_max,r_f_norm,v_f_norm) ;

t_guess = Tk2e01(end);
nrow= size(Tk2e01,1);
Gnorm = zeros(nrow,1);
dA_norm = zeros(nrow,1);
dL_norm = zeros(nrow,1);
r_norm = zeros(nrow,1);
for row=1:nrow
     [~,Gnorm(row),dA_norm(row),dL_norm(row),dE_norm(row),r_norm(row)]=...
         closelooporbit_standard(Tk2e01(row),Yk2e01(row,:)',F_max,r_f_norm,v_f_norm);
 end
t_final = Tk2e01(end);
figure
plot(Tk2e01*canonical_time,dA_norm,'r-',Tk2e01*canonical_time,...
     dL_norm,'b-',Tk2e01*canonical_time,dE_norm,'k-')
legend('dA','dL','dE') 
figure
plot(Tk2e01,r_norm,'b-')
xlabel('Time')
ylabel('r_norm')
figure
plot(Yk2e01(:,1),Yk2e01(:,2))
xlabel('x')
ylabel('y')
xlim([-2 2])
ylim([-2 2])


% figure(10)
% set(gca,'FontSize',18)
% plot(Yk2e01(:,1),Yk2e01(:,2),'LineWidth',2)
% 
end
function [value,isterminal,direction] = events(t,y,F_max,r_f_norm,v_f_norm)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration.

r_f = [r_f_norm,0,0]';
v_f = [0,v_f_norm,0]';
L_T = cross(r_f,v_f);
A_T= [0,0,0]';

r = y(1:3);
v = y(4:6);

L= cross(r,v);
dL = L-L_T;
dL_norm = norm(dL);

A= cross(v,L)-r/norm(r);
dA= A-A_T;
dA_norm = norm(dA);
value = dA_norm <0.001 && norm(y(1:3)) > 0.999999*r_f_norm ;
%value = norm(y(1:3)) < 0.999*r_f_norm;     % detect height = 0
norm(y(1:3));
%value = norm(y(1:3)) < 0.9999999999*r_f_norm  
%value = t<4.99
value = 0;
%value = dA_norm > 9E-3;
%value = dA_norm > 0.2;
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
end

function [dy,normG,dA_norm,dL_norm,dE_norm,r_norm]=closelooporbit_standard(...
    t,y,...
   F_max,r_f_norm,v_f_norm) %...
   % )
%F_max =0.01;
mu =1;
epsil= 0.0001;
01;
k=0.1;

%F_max = 0.01;
dy = y;
r=y(1:3);
v=y(4:6);
dy(1:3)=v;

r_f = [r_f_norm,0,0]';
v_f = [0,v_f_norm,0]';

L_T = cross(r_f,v_f);
E_T = 0.5*norm(v_f)^2-mu/norm(r_f);
A_T= cross(v_f,L_T)-mu*r_f/norm(r_f);

L= cross(r,v);
dL = L-L_T;
dL_norm = norm(dL);

A= cross(v,L)-mu*r/norm(r);
dA= A-A_T;
dA_norm = norm(dA);

E = 0.5*norm(v)^2-mu/norm(r);
dE = (E-E_T);
G= -  (k*cross(dL , r) + 1*cross(L,dA)+ cross(cross(dA,v), r));
%G= -(1.*cross(dL , r) + 1.0*dE*v );
dE_norm = norm(dE);

normG = norm(G);
if normG <= epsil*F_max
   F=1/epsil*G;
else
    F=F_max*G/norm(G);
end

dy(4:6)=-mu/norm(r)^3.*(r)+F;%+a_J2;
r_norm = norm(r);

dot(r,v);
end 