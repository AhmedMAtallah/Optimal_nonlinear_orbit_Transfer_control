clc;
clear all;
mu = 3.986004418e5;            % Earth’s gravitational parameter [km^3/s^2]
R_earth = 6378.140;         % Earth radius [km]

% Molniya Orbit 1
a     = 26554;
e     = 0.72;
i     = 63.4;
omega = 0;
theta = 0;
RAAN  = 0;
h =(a*mu*(1 - e^2))^0.5;
% Calculating initial state vector
[ R0 V0 ] = Orbital2State( h, i, RAAN, e,omega,theta);
r = norm(R0);           % Initial radius  [km]
v = norm(V0);           % Initial speed   [km/s]
T = 2*pi*a^1.5/mu^0.5;  % Orbital period [s]
dt = T/10000;           % time step [s]
t = 0;                  % initial time
% Using fourth-order Runge–Kutta method to solve fundamental equation
% of relative two-body motion
F_r = @(R) -mu/(norm(R)^3)*R;
Rd = V0; R  = R0;
i = 1;
while (t <= T)
    Rv(i,:) = R;
    tv(i) =t;
    k_1 = dt*F_r(R);
    k_2 = dt*F_r(R+0.5*k_1);
    k_3 = dt*F_r(R+0.5*k_2);
    k_4 = dt*F_r(R+k_3);
    Rd  = Rd + (1/6)*(k_1+2*k_2+2*k_3+k_4);
    R   = R + Rd*dt;
    t   = t+dt;
    i = i+1;
end

% Molniya Orbit 2
a     = 26554;
e     = 0.72;
i     = 63.4;
omega = -90;
theta = 0;
RAAN  = 250;
h =(a*mu*(1 - e^2))^0.5;
% Calculating initial state vector
[ R0 V0 ] = Orbital2State( h, i, RAAN, e,omega,theta);
r = norm(R0);           % Initial radius  [km]
v = norm(V0);           % Initial speed   [km/s]
T = 2*pi*a^1.5/mu^0.5;  % Orbital period [s]
dt = T/10000;           % time step [s]
t = 0;                  % initial time
% Using fourth-order Runge–Kutta method to solve fundamental equation
% of relative two-body motion
F_r = @(R) -mu/(norm(R)^3)*R;
Rd = V0; R  = R0;
i = 1;
while (t <= T)
    Rv2(i,:) = R;
    tv(i) =t;
    k_1 = dt*F_r(R);
    k_2 = dt*F_r(R+0.5*k_1);
    k_3 = dt*F_r(R+0.5*k_2);
    k_4 = dt*F_r(R+k_3);
    Rd  = Rd + (1/6)*(k_1+2*k_2+2*k_3+k_4);
    R   = R + Rd*dt;
    t   = t+dt;
    i = i+1;
end
% Plotting
figure('Color',[0 0 0]);
figure(1);
hold on;
load('topo.mat','topo','topomap1');
colormap(topomap1);
% Create the surface.
radius_earth=6378;
[x,y,z] = sphere(50);
x =radius_earth*x;
y =radius_earth*y;
z =radius_earth*z;

props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;
surface(x,y,z,props);

% Inertial Frame Axis
Xa  = [radius_earth:100:radius_earth+2500];
Z0  = Xa*0;
plot3(-Xa,Z0,Z0,'r')
plot3(Z0,-Xa,Z0,'y')
plot3(Z0,Z0,Xa,'g')
% Plotting Orbits
plot3(Rv(:,1),Rv(:,2),Rv(:,3));
plot3(Rv2(:,1),Rv2(:,2),Rv2(:,3),'y');
axis square off
view(3)
zoom(2)