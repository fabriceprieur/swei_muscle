% This function computes the integrand for the gr function. It is called
% symbolically to evaluate its two dimensional integration over theta and
% phi. Theory can be found in Gridin: "Far field asymptotics of the Green’s
% tensor for a transversely isotropic solid" 2000 and Wang 1995 "Three-
% dimensional time-harmonic elastodynamic Green’s functions for anisotropic
% solids"
% 
% Author: Fabrice Prieur
% Creation date: 09 Jan 2018
% Copyright: University of Oslo

function f_int = getIntegrandGr(phi,theta)

global omega rho c11 c12 c13 c33 c44 r iline icolumn v psi l


x=l*sin(v)*cos(psi);
y=l*sin(v)*sin(psi);
z=l*cos(v);

% n1, n2, and n3 are changed to their opposite when theta becomes pi-theta
% and when phi becomes pi+phi. This corresponds to changing the normal
% vector n to -n.
% When n1, n2, and n3 are changed to their opposite, cr, Mr and the factor
% do not change so the integration can be done on the half sphere (theta
% varies between 0 and pi and phi varies between 0 and pi instead of from
% -pi to pi) and be multiplied by 2 instead of the whole sphere.
n1=sin(theta).*cos(phi);
n2=sin(theta).*sin(phi);
n3=cos(theta);

dS=sin(theta); % in fact dS=r dtheta (r sin theta) d phi, but r=1

c66=(c11-c12)/2;
MrCoef=getMrCoef(r,n1,n2,n3,c11,c13,c33,c44,iline,icolumn);
cr=getCr(r,n1,n2,n3,c11,c13,c33,c44,c66,rho);

factor=exp(1i*omega*(l*abs(n1*sin(v)*cos(psi)+n2*sin(v)*sin(psi)+n3.*cos(v))./cr))./(cr.^3).*dS;

% summing along the first dimension (where n varies)
f_int=1i*omega/(16*pi^2*rho)*MrCoef.*factor;