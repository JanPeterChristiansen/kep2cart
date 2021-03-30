function [r,v] = kep2cart(a,e,w,omega,i,M0,t)
% Calculates Cartesian state vectors from a set of Kelperian orbit elements
% 
% Inputs: 
% a - semi major axis [m]
% e - eccentricity
% w - argument of periapsis [rad]
% omega - longitude of ascending node [rad]
% i - inclination [rad]
% M0 - mean anomaly at time 0
% t - time
% 
% Outputs:
% r - cartesian position vector at time t [m]
% v - cartesian velocity vector at time t [m/s]
G = 6.67259E-11;
M = 5.9736E24;
mu = G*M;

n = sqrt(mu/a^3);
M = mod(M0 + n*t,2*pi);

E = solveKeplerEq(M,e,1e-10);

nu = 2*atan2(sqrt(1+e)*sin(E/2),sqrt(1-e)*cos(E/2));

rc = a*(1-e*cos(E));

o = rc * [cos(nu) sin(nu) 0]';
odot = sqrt(mu*a)/rc * [-sin(E) sqrt(1-e^2)*cos(E) 0]';

r = [o(1)*(cos(w)*cos(omega)-sin(w)*cos(i)*sin(omega)) - o(2)*(sin(w)*cos(omega)+cos(w)*cos(i)*sin(omega));
    o(1)*(cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega)) + o(2)*(cos(w)*cos(i)*cos(omega)-sin(w)*sin(omega));
    o(1)*(sin(w)*sin(i)) + o(2)*(cos(w)*sin(i))];

v = [odot(1)*(cos(w)*cos(omega)-sin(w)*cos(i)*sin(omega)) - odot(2)*(sin(w)*cos(omega)+cos(w)*cos(i)*sin(omega));
    odot(1)*(cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega)) + odot(2)*(cos(w)*cos(i)*cos(omega)-sin(w)*sin(omega));
    odot(1)*(sin(w)*sin(i)) + odot(2)*(cos(w)*sin(i))];
end

function E = solveKeplerEq(M,e,tol)
% Solves Kepler's equation M(t) = E(t) - e*sin(E(t)) for E(t) via the
% Newton-Raphson method. 
% 
% Inputs: 
% M - Mean anomaly at time t [rad]
% e - eccentricity
% tol - desired numerical accuracy
%
% Outputs: 
% E - eccentric anomaly at time t [rad]
E = M;
E(end+1)= E(end) - (E(end)-e*sin(E(end))-M)/(1-e*cos(E(end)));
while( abs(E(end)-E(end-1)) > tol )
    E(end+1) = E(end) - (E(end)-e*sin(E(end))-M)/(1-e*cos(E(end)));
end
E = E(end);
end
