% This function returns the coefficients of the matrix A/tr(A) where A is
% the adjunct matrix of (L-rho c_i I3). Three matrices exist for the P, SH,
% and SV waves. Each matrix is a 3x3 matrix. The first, second, and third
% columns correspond to displacement due to a forcde applied in the x, y,
% or z direction, respectively. Each line correspond to the x, y, and z,
% component of the displacement generated.
% The theory can be found in Gridin 200 : "Far field asymptotics of the
% Green’s tensor for a transversely isotropic  solid" and Wang 1995 "Three-
% dimensional time-harmonic elastodynamic  Green’s functions for 
% anisotropic solids"
% 
% Author: Fabrice Prieur
% Creation date: 09 Jan 2018
% Copyright: University of Oslo

function MrCoef=getMrCoef(r,n1,n2,n3,c11,c13,c33,c44,i,j)
% r=1, 2, or 3 for P, SH, and SV waves,
% n1, n2, and n3 can be arrays and contain the components of unit vectors
% (the unit sphere is scanned in the integration involving Mr).
% cij are the components of the stiffness tensor
% i, and j are the line and column indices of the desired matrrix coef.

np2=n1.^2+n2.^2;
n32=n3.^2;
D=((c11-c44)*np2+(c44-c33)*n32).^2+4*(c13+c44)^2*np2.*n32;


Mr=zeros(size(n1,1),size(n1,2),3,3);
if r==1
    if i==1 && j==1
    Mr(:,:,1,1)=-n1.*n1./np2.*((c44-c11)*np2+(c33-c44)*n32-sqrt(D))./(2*sqrt(D));
    elseif i==1 && j==2
    Mr(:,:,1,2)=-n1.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32-sqrt(D))./(2*sqrt(D));
    elseif i==1 && j==3
    Mr(:,:,1,3)=n1.*n3*(c44+c13)./sqrt(D);
    elseif i==2 && j==1
    Mr(:,:,2,1)=-n1.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32-sqrt(D))./(2*sqrt(D));
    elseif i==2 && j==2
    Mr(:,:,2,2)=-n2.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32-sqrt(D))./(2*sqrt(D));
    elseif i==2 && j==3
    Mr(:,:,2,3)=n2.*n3*(c44+c13)./sqrt(D);
    elseif i==3 && j==1
    Mr(:,:,3,1)=n1.*n3*(c44+c13)./sqrt(D);
    elseif i==3 && j==2
    Mr(:,:,3,2)=n2.*n3*(c44+c13)./sqrt(D);
    else
    Mr(:,:,3,3)=((c44-c11)*np2+(c33-c44)*n32+sqrt(D))./(2*sqrt(D));
    end
elseif r==2
    if i==1 && j==1
    Mr(:,:,1,1)=n1.*n1./np2.*((c44-c11)*np2+(c33-c44)*n32+sqrt(D))./(2*sqrt(D));
    elseif i==1 && j==2
    Mr(:,:,1,2)=n1.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32+sqrt(D))./(2*sqrt(D));
    elseif i==1 && j==3
    Mr(:,:,1,3)=-n1.*n3*(c44+c13)./sqrt(D);
    elseif i==2 && j==1
    Mr(:,:,2,1)=n1.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32+sqrt(D))./(2*sqrt(D));
    elseif i==2 && j==2
    Mr(:,:,2,2)=n2.*n2./np2.*((c44-c11)*np2+(c33-c44)*n32+sqrt(D))./(2*sqrt(D));
    elseif i==2 && j==3
    Mr(:,:,2,3)=-n2.*n3*(c44+c13)./sqrt(D);
    elseif i==3 && j==1
    Mr(:,:,3,1)=-n1.*n3*(c44+c13)./sqrt(D);
    elseif i==3 && j==2
    Mr(:,:,3,2)=-n2.*n3*(c44+c13)./sqrt(D);
    else
    Mr(:,:,3,3)=-((c44-c11)*np2+(c33-c44)*n32-sqrt(D))./(2*sqrt(D));
    end
else
    if i==1 && j==1
    Mr(:,:,1,1)=1-n1.*n1./np2;
    elseif i==1 && j==2
    Mr(:,:,1,2)=-n1.*n2./np2;
    elseif i==2 && j==1
    Mr(:,:,2,1)=-n1.*n2./np2;
    else
    Mr(:,:,2,2)=1-n2.*n2./np2;
    end
end
MrCoef=squeeze(Mr(:,:,i,j));
end