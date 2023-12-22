clear all;
clc;
% Six  orbital elements are:
h    = 7.139652961424629e+04;      % [km^2/s] Specific angular momentum
i    = 63.4;         % [deg] Inclination
RAAN = 0;         % [deg] Right ascension (RA) of the ascending node
e    = 0.74;        % Eccentricity
omega= 0;         % [deg] Argument of perigee
theta= 0;         % [deg] True anomaly
mu = 3.986004418e5;       % Earth’s gravitational parameter [km^3/s^2]
% Components of the state vector of a body relative to its perifocal
% reference
rx = h^2/mu*(1/(1 + e*cosd(theta)))*[cosd(theta);sind(theta);0];
vx = mu/h*[-sind(theta); (e +cosd(theta));0];
% Direction cosine matrix
QXx = [cosd(omega), sind(omega),0;-sind(omega),cosd(omega),0;0,0,1]*...
    [1,0,0;0,cosd(i),sind(i);0,-sind(i),cosd(i)]*...
    [cosd(RAAN), sind(RAAN),0;-sind(RAAN),cosd(RAAN),0;0,0,1];
% Transformation Matrix
QxX = inv(QXx);
% Geocentric equatorial position vector R
R = QxX*rx;
% Geocentric equatorial velocity vector V
V = QxX*vx;

canonical_distance = 6378.137;
canonical_time = 806.812;
canonical_accel = 9.8e-3;
canonical_vel = canonical_distance/canonical_time;
R =R / canonical_distance
V = V / canonical_vel
fprintf('R = %4.2f*i +  %4.2f*j + %4.2f*k [km]\n',R);
fprintf('V = %4.4f*i +  %4.4f*j + %4.4f*k [km/s]\n',V);