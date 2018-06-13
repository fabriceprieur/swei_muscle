% This function computes the integrand for the gs function. It is called
% symbolically to evaluate its one dimensional integration over beta.
% Theory can be found in Gridin: "Far field asymptotics of the Green’s
% tensor for a transversely isotropic solid" 2000 and Wang 1995 "Three-
% dimensional time-harmonic elastodynamic Green’s functions for anisotropic
% solids"
% The integration here is done on the unit circle containing the origin and
% orthogonal to the vector joining the origin and the observation point.
% This coordinates of the circle have first to be defined then the
% integration can be done.
% 
% Author: Fabrice Prieur
% Creation date: 09 Jan 2018
% Copyright: University of Oslo

function f_int = getIntegrandGs(beta)

global  rho c11 c12 c13 c33 c44 r iline icolumn v psi l

% a and b are two unit vectors perpendicular to the postion vector r and
% perpendicular to each other. Any combination of a and b is orthogonal to
% the position vector.
% b is in the (x,y) plane and a is found as a=e x b where e = r/|r|
ax=-cos(v).*cos(psi);
ay=-cos(v).*sin(psi);
az=sin(v);
bx=sin(psi);
by=-cos(psi);
bz=0*bx;

% as beta varies between 0 and 2pi it describes the circle of radius 1
% containing the origin and with normal the position vector (all unit
% vectors obtained with combining a and b).
n1=ax.*cos(beta)+bx.*sin(beta);
n2=ay.*cos(beta)+by.*sin(beta);
n3=az.*cos(beta)+bz.*sin(beta);

c66=(c11-c12)/2;
MrCoef=getMrCoef(r,n1,n2,n3,c11,c13,c33,c44,iline,icolumn);
cr=getCr(r,n1,n2,n3,c11,c13,c33,c44,c66,rho);

f_int=1/(8*pi^2*rho)*MrCoef./(cr.^2)/l; % we divide by l=|r| since the
% factor contains Dirac(r.xsi)=Dirac(|r| e.xsi)
