% This script computes the harmonic Green's functions for a transverse
% isotropic medium. It is based on the integral formulation describes in
% Eq. 2.3 in "Far-Field asymptotics of the Green’s tensor for a transversely
% isotropic solid" D. Gridin 2000.
% It consists of two integrals one on the unit sphere (first term in Eq.
% 2.3) and one on the unit circle containing the origin and normal to the
% vector joining the origin and the observation point.
%
% This script is very slow! It is due to the difficulties in evaluating the
% integrals either in the far field or for high frequency due to the fast
% variation of the complex exponential term in the first integral.
%
% Author: Fabrice Prieur
% Creation date: 7 Jun. 2018
% Copyright: University of Oslo

sf=1;   % scale factor
tsignal=300e-6; % pulse duration
tsig=sf*1.6e-3; % duration of modelled signal
fs=100e3;   % sampling frequency
t=(0:1/fs:tsig);
freq=(0:length(t)-1)/length(t)*fs;
ix=find(t<=tsignal);
signal=zeros(1,length(t));
signal(ix)=hann(length(ix)).'; % Hann window over tisgnal
specsignal=fft(signal); % ideal frequency spectrum
specRed=specsignal(1:sf*25); % truncated freq. spectrum
specRecon=0*specsignal; % Spec. of reconstructed signal
specRecon(end-length(specRed)+2:end)=conj(fliplr(specRed(2:end)));
signalRecon=real(ifft(specRecon)); % Reconstructed signal 
% amplitude has to be adjusted since only part of the spectrum is taken
% into account
coef=max(signal)/max(signalRecon);
signalRecon=signalRecon*coef;
freq(freq==0)=eps; % Green functions diverge for f = 0;

global omega rho c11 c12 c13 c33 c44 r iline icolumn v psi l

rho=1000;
% The following values of the elasticity tensor are taken from Royer 2011
% "On the elasticity of transverse isotropic soft tissues" and Genisson
% 2010 "VISCOELASTIC AND ANISOTROPIC MECHANICAL PROPERTIES OF
% IN VIVO MUSCLE TISSUE ASSESSED BY SUPERSONIC SHEAR IMAGING"
c44=2.7^2*rho; % use speed of SH waves perp to fibers
c66=1.2^2*rho; % use speed of SH waves parallel to fibers

c11=2.6e9;
c12=c11-2*c66;
c13=3.29e9;
c33=4.17e9;

% Compute the harmonic Green's function for each frequency of the reduced
% frequency spectrum.
for ii=25:length(specRed)
    omega=2*pi*freq(ii);

    fprintf('Computing for frequency %d out of %d\n',ii,length(specRed));
    % psi: angle of (x,y) plane of projection of position point
    psi=0; % we stay in the (x,z) plane

    srcExt=3; % half source extension in mm along x direction (between +/- srcExt)
    spatStep=1e-3; % domain in which we require better resolution
    x=horzcat(linspace(eps,srcExt*1e-3,20),(srcExt*1e-3+spatStep:spatStep:30e-3));
    z=horzcat(linspace(eps,srcExt*1e-3,20),(srcExt*1e-3+spatStep:spatStep:40e-3));
    [X,Z]=meshgrid(x,z);
    vlist=atan(X./Z);
    rad=sqrt(X.^2+Z.^2);
    % will compute displacements due to force along x
    GrCoef1=zeros(size(X,1),size(X,2),3,2); % will contain double integration on unit sphere
    GsCoef1=GrCoef1; % will contain single integration on unit circle perpendicular to position vector
    % will compute displacements due to force along z
    GrCoef3=GrCoef1; % will contain double integration on unit sphere
    GsCoef3=GrCoef1; % will contain single integration on unit circle perpendicular to position vector
    for i=1:size(X,1)
        fprintf('computing for line %d of %d\n',i,size(X,1));
        for j=1:size(X,2)
        % v: angle between position point and z axis
            v=vlist(i,j);
            l=rad(i,j);
            for r=1:3 % go through contribution from P, SV, and SH waves
                icolumn=1; % due to force along x
                iline=1; % x displacement
                GrCoef1(i,j,r,1)=2*integral2(@getIntegrandGr,0,pi,0,pi,'AbsTol',1e-12);
                GsCoef1(i,j,r,1)=integral(@getIntegrandGs,-pi,pi,'AbsTol',1e-12);
                iline=3; % z displacement
                GrCoef1(i,j,r,2)=2*integral2(@getIntegrandGr,0,pi,0,pi,'AbsTol',1e-12);
                GsCoef1(i,j,r,2)=integral(@getIntegrandGs,-pi,pi,'AbsTol',1e-12);
                icolumn=3; % due to force along z
                if r<3 % No need to compute for r=3 since M^3_13=M^3_33=0
                    iline=1; % x displacement
                    GrCoef3(i,j,r,1)=2*integral2(@getIntegrandGr,0,pi,0,pi,'AbsTol',1e-12);
                    GsCoef3(i,j,r,1)=integral(@getIntegrandGs,-pi,pi,'AbsTol',1e-12);
                    iline=3; % z displacement
                    GrCoef3(i,j,r,2)=2*integral2(@getIntegrandGr,0,pi,0,pi,'AbsTol',1e-12);
                    GsCoef3(i,j,r,2)=integral(@getIntegrandGs,-pi,pi,'AbsTol',1e-12);
                end
            end
        end
    end

    fname=['Data\GreenFunctions_',int2str(ii-idx+6)];
    save(fname','GrCoef1','GsCoef1','GrCoef3','GsCoef3','x','z','freq','t'...
        ,'pulse','Fc','fs','omega','c11','c12','c13','c33','c44','c66','-v7.3');
end