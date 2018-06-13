% This function computes the propagation speeds for P, SH, or SV waves in
% the direction specified. The theory can be found in Gridin 200 : "Far
% field asymptotics of the Green’s tensor for a transversely isotropic 
% solid" and Wang 1995 "Three-dimensional time-harmonic elastodynamic 
% Green’s functions for anisotropic solids"
% 
% Author: Fabrice Prieur
% Creation date: 09 Jan 2018
% Copyright: University of Oslo

function cr=getCr(r,n1,n2,n3,c11,c13,c33,c44,c66,rho)
% r=1, 2, or 3 for P, SH, and SV waves,
% n1, n2, and n3 can be arrays and contain the components of unit vectors
% (the unit sphere is scanned in the integration involving cr).
% cij are the components of the stiffness tensor
% rho is the density

np2=n1.^2+n2.^2;
n32=n3.^2;
D=((c11-c44)*np2+(c44-c33)*n32).^2+4*(c13+c44)^2*np2.*n32;

if r==1
    cr=sqrt(1/2*((c11+c44)*np2+(c33+c44)*n32+sqrt(D))/rho);
elseif r==2
    cr=sqrt(1/2*((c11+c44)*np2+(c33+c44)*n32-sqrt(D))/rho);
else
    cr=sqrt((c66*np2+c44*n32)/rho);
end

end