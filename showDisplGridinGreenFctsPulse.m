% This script uses the time-varying displacement due to a point force
% located at the origin and directed orthogonal or parallel to the fiber
% direction computed by the script displGridinPulse.m and saved in the file
% displGridinPulse.mat. This displacement is convolved with the spatial
% distribution of the ARF generated by a LT 7-4 linear array and saved in
% the file ARF_L7_4_17p8mFoc_4MHz.mat. The displacement is computed for
% various angle between the applied ARF and the fiber direciton.
% The displacement distribution is then plotted at a given time for 3
% angles. The RMS displacements along a vertical line as a function of
% angle between the ARF and fibers is also shown.
%
% Author: Fabrice Prieur
% Creation date: 05 June 2018
% Copyright: University of Oslo

clear all
load('Data\displGridinPulse.mat');

% load results file from k-wave simulation containing ARF from L7-4 probe
load('Data\ARF_L7_4_17p8mFoc_4MHz.mat');
% plotting
xax_int=(xax(1):kgrid.dx/5:xax(end));
yax_int=(yax(1):kgrid.dy/5:yax(end));
[X,Y]=meshgrid(xax,yax);
Fint=interp2(sqrt(Fx.^2+Fy.^2),5);
figure('color','w','Position',[251,378,531,420]);
imagesc(yax_int*1e3,xax_int*1e3,Fint/max(Fint(:)));hold on;
caxis([0 1]);
quiver(yax*1e3,xax*1e3,Fy,Fx,'k');
xlabel('lateral distance [mm]','fontweight','bold','fontsize',11);
ylabel('depth [mm]','fontweight','bold','fontsize',11);
axis('image');
xlim([-3 3]);ylim([15 20]);
title('Spatial distribution of the ARF');colorbar
set(gca,'position',[0.1,0.11,0.78,0.815]);
% Careful, in this file x is still vertical but y is horizontal (not z).
% Grid on which ARF was computed (beam focus is placed at origin)

src_X=kgrid.x_vec-min(kgrid.x_vec)-transducer.focus_distance;
src_Y=kgrid.y_vec;
[src_Y src_X]=meshgrid(src_Y,src_X);
% size of src area (in muscle axes). We only take into account the ARF
% within this distance of the focal point.
srcExt=3e-3; 
focus=17.8e-3; % focus depth for ARF
xax=xax-focus;
[Z_srcExtent,X_srcExtent]=meshgrid(yax(abs(yax)<=srcExt),xax(abs(xax)<=srcExt));

% Create interpolant that will be used with rotated grids
funcFx=griddedInterpolant(src_Y.',src_X.',Fx.');
funcFy=griddedInterpolant(src_Y.',src_X.',Fy.');

% angle between ARF main axis and fiber orientation (z horizontal axis)
theta=[90 60 30]/180*pi;

%%
% convolve the displacements from a point source with the ARF distribution
% to get the displacement generated by the ARF at a given time point in the
% whole field

toff=2.1e-3; % time at which displacements will be displayed
j=find(tred>=toff,1,'first');
figure('color','w','position',[520 154 1251 403]);
for i=1:length(theta)
    fprintf('Computing disp from ARF step %d out of %d\n',i,length(theta));
    % coordinates (in probe axes) where we wish to interpolate the ARF
    % components as the angle between probe dir and muscle fibre changes
    ARF_X=-Z_srcExtent*sin(theta(i)-pi/2)+X_srcExtent*cos(theta(i)-pi/2);
    ARF_Y=Z_srcExtent*cos(theta(i)-pi/2)+X_srcExtent*sin(theta(i)-pi/2);
    % x and y components of ARF in the defined src area (in probe axes)
    thetaARFx=funcFx(ARF_Y.',ARF_X.').';
    thetaARFy=funcFy(ARF_Y.',ARF_X.').';
    
    nl=round(size(thetaARFx,1)/2);
    nc=round(size(thetaARFx,2)/2);
    zfig=zint(nc:end-(size(thetaARFx,2)-nc));
    xfig=xint(nl:end-(size(thetaARFx,1)-nl))+focus;
    if i==1
        xdisp=zeros(length(xfig),length(zfig),length(tred));
        zdisp=xdisp;
    end
    
    taper=repmat(tukeywin(size(Z_srcExtent,1)),1,size(Z_srcExtent,2)).*...
        repmat(tukeywin(size(Z_srcExtent,2)).',size(Z_srcExtent,1),1);
    % x and z components of applied force field (in muscle axes)
    srcFx=taper.*(thetaARFx*cos(theta(i)-pi/2)+thetaARFy*sin(theta(i)-pi/2));
    srcFz=taper.*(-thetaARFx*sin(theta(i)-pi/2)+thetaARFy*cos(theta(i)-pi/2));

    % x and z displacement components (in muscle axes)    
    xdisp(:,:,j)=(conv2(xdispShearFxInt(:,:,j).',srcFx,'valid')+...
        conv2(xdispShearFzInt(:,:,j).',srcFz,'valid'))/sum(sqrt(srcFx(:).^2+srcFz(:).^2));
    zdisp(:,:,j)=(conv2(zdispShearFxInt(:,:,j).',srcFx,'valid')+...
        conv2(zdispShearFzInt(:,:,j).',srcFz,'valid'))/sum(sqrt(srcFx(:).^2+srcFz(:).^2));
    totaldisp=sqrt(xdisp.^2+zdisp.^2);
 
    
    subplot(1,length(theta),i);
    if i==1
        nor=max(max(abs(totaldisp(:,:,j))));
    end
    imagesc(zfig*1e3,xfig*1e3,abs(totaldisp(:,:,j))/nor);
    [Xfig,Zfig]=meshgrid(xfig,zfig);
    rad=sqrt(Xfig.^2+Zfig.^2);
    line([0 max(rad(:))*cos(theta(i))]*1e3,[focus focus+max(rad(:))*sin(theta(i))]*1e3,'color','r','linewidth',2)    
    caxis([0 1]);daspect([1 1 1]);
    ds=4;hold on;
    quiver(zfig(1:ds:end)*1e3,xfig(1:ds:end)*1e3,squeeze(zdisp(1:ds:end,1:ds:end,j)),squeeze(xdisp(1:ds:end,1:ds:end,j)),'w');
    title(['\theta = ',int2str(round((theta(i))/pi*180)),' deg, t = ',int2str(tred(j)*1e6),' \mus']);
    if i==1
        ylabel('x [mm]','fontweight','bold','fontsize',11);
    elseif i==3
        colorbar;
    end
    xlabel('z [mm]','fontweight','bold','fontsize',11);
    set(gca,'position',[0.0415+0.305*(i-1) 0.11 0.27 0.815]);
    colormap jet
    xlim([-8 8]);ylim([-8 8]+focus*1e3);    
    
end

%%
% Convolve the displacements from a point source with the spatial
% distribution of the ARF to get the RMS displacement at a given lateral
% offset as a function of the angle between the ARF and the fibers.

% Same loop out of laziness. We only look at the displacement at a given
% lateral distance though.
theta=(90:-10:0)/180*pi;
zoff=[2 4 6]*1e-3; % lateral offsets
tsignal=480e-6;
for i=1:length(theta)
    fprintf('Computing disp from ARF step %d out of %d\n',i,length(theta));
    % coordinates (in probe axes) where we wish to interpolate the ARF
    % components as the angle between probe dir and muscle fibre changes
    ARF_X=-Z_srcExtent*sin(theta(i)-pi/2)+X_srcExtent*cos(theta(i)-pi/2);
    ARF_Y=Z_srcExtent*cos(theta(i)-pi/2)+X_srcExtent*sin(theta(i)-pi/2);
    % x and y components of ARF in the defined src area (in probe axes)
    thetaARFx=funcFx(ARF_Y.',ARF_X.').';
    thetaARFy=funcFy(ARF_Y.',ARF_X.').';
    
    nl=round(size(thetaARFx,1)/2);
    nc=round(size(thetaARFx,2)/2);
    zfig=zint(nc:end-(size(thetaARFx,2)-nc));
    xfig=xint(nl:end-(size(thetaARFx,1)-nl))+focus;
    if i==1 && j==1
        zline=zeros(length(xfig),length(theta),length(zoff));
        xdisp=zeros(length(xfig),length(zfig),length(tred));
        zdisp=xdisp;
        totaldisp=zeros(length(xfig),length(zfig),length(tred),length(theta));
    end
    
    taper=repmat(tukeywin(size(Z_srcExtent,1)),1,size(Z_srcExtent,2)).*...
        repmat(tukeywin(size(Z_srcExtent,2)).',size(Z_srcExtent,1),1);
    % x and z components of applied force field (in muscle axes)
    srcFx=taper.*(thetaARFx*cos(theta(i)-pi/2)+thetaARFy*sin(theta(i)-pi/2));
    srcFz=taper.*(-thetaARFx*sin(theta(i)-pi/2)+thetaARFy*cos(theta(i)-pi/2));
    
    % Here we have to do the convolution at every time step since we're
    % interested in the RMS of the displacement.
    for jj=1:length(tred)
        % x and z displacement components (in muscle axes)    
        xdisp(:,:,jj)=(conv2(xdispShearFxInt(:,:,jj).',srcFx,'valid')+...
            conv2(xdispShearFzInt(:,:,jj).',srcFz,'valid'))/sum(sqrt(srcFx(:).^2+srcFz(:).^2));
        zdisp(:,:,jj)=(conv2(zdispShearFxInt(:,:,jj).',srcFx,'valid')+...
            conv2(zdispShearFzInt(:,:,jj).',srcFz,'valid'))/sum(sqrt(srcFx(:).^2+srcFz(:).^2));
    end
    totaldisp(:,:,:,i)=sqrt(xdisp.^2+zdisp.^2);
        
    for j=1:length(zoff)
        zline(:,i,j)=rms(totaldisp(:,find(zfig>=zoff(j),1,'first'),find(tred>=tsignal*2,1,'first'):end),3);
    end
 
end    

%% Plot

figure('color','w','position',[520 154 1251 403]);
for j=1:length(zoff)
    subplot(1,length(zoff),j);
    imagesc(rad2deg(theta),xfig*1e3,zline(:,:,j)/max(max(zline(:,:,j))));colormap copper
    xlabel('\theta [deg]','fontweight','bold','fontsize',11);caxis([0 1]);
    title(['RMS displacement at z = ',int2str(round(zoff(j)*1e3)),' mm'],...
    'fontweight','bold','fontsize',11);
    if j==1
        ylabel('x [mm]','fontweight','bold','fontsize',11);
    elseif j==3
        colorbar;
    end
    set(gca,'position',[0.0415+0.305*(j-1) 0.11 0.27 0.815]);
    set(gca,'xdir','reverse')
    colormap jet
    h=line([rad2deg(theta(end)) rad2deg(theta(1))],[focus focus]*1e3,...
        'linestyle','--','color','k','linewidth',2);
end