% This file simulate the pressure and particle velocity fields emitted by a
% 3D transducer. The transducer is supposed to be an ATLl7-4 linea array
% described in "SHEAR WAVE SPEED MEASUREMENT USING AN UNFOCUSED ULTRASOUND 
% BEAM" Zhao 2012, UMB 38 (9)
% The results for the three particle velocity components and the pressure
% are saved into a mat file that is later treated to compute the acoustic
% radiation force.
%
% To run this file you should have installed k-Wave: http://k-wave.org/ and
% the k-Wave installation folder should be in the Matlab path.
%
% Author: Fabrice Prieur
% Creation date: 05 Feb. 2018
% Copyright: University of Oslo

clear all;

% simulation settings
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 10*2;            % [grid points]
PML_Y_SIZE = 10*2;            % [grid points]
PML_Z_SIZE = 5*2;            % [grid points]
scaleFactor=2;

% set total number of grid points not including the PML
Nx = 88*scaleFactor+1;    % [grid points]
Ny = 170*scaleFactor+1;    % [grid points]
Nz = 22*scaleFactor+1;     % [grid points]
    
dx = 0.308e-3/scaleFactor;                  % [m]

% calculate the spacing between the grid points
dy = dx;                    % [m]
dz = dx*2;                    % [m]

% set desired grid size in the x-direction not including the PML
x = Nx*dx;                  % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1540;      % [m/s]
medium.density = 980;          % [kg/m^3]
medium.alpha_coeff = 0.5;      % [dB/(MHz^y cm)]
medium.alpha_power = 1;
medium.alpha_mode ='no_dispersion'; % needed if y=1;
% medium.BonA = 6;

% create the time array
t_end = 25e-6;                  % [s] 
cfl=0.2;
kgrid.t_array = makeTime(kgrid, medium.sound_speed, cfl, t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 4e6;        % [Hz] 
tone_burst_cycles = 15;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles,'Envelope','Gaussian');

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength./(medium.sound_speed*medium.density)).*input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer

transducer.number_elements = 128;
transducer.element_width =  1*scaleFactor;% 0.283 mm width + 0.025 mm kerf
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer.element_length = floor(22.7*scaleFactor/2);% 7 mm height
transducer.elevation_focus_distance = 25e-3;% focus distance in the elevation plane [m]
transducer.sound_speed = 1540;              % sound speed [m/s]
transducer.focus_distance = 17.8e-3;          % focus distance [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = floor([1, Ny/2 - transducer_width/2+1, Nz/2 - transducer.element_length/2+1]);

            % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(1:transducer.number_elements) = 1;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE SENSOR MASK
% =========================================================================

% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
lz=2;
sensor.mask(1:Nx, 1:Ny, ceil(Nz/2)-5:ceil(Nz/2)+5) = 1;
sensor.record = {'p','u_non_staggered'};

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {'DisplayMask', transducer.all_elements_mask, ...
    'PMLInside', false, 'PlotPML', true, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-source_strength/2, source_strength/2],'UseSG',false};

% run the simulation
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});
tstart=20e-6;
tstop=t_end;
nstart=find(kgrid.t_array>tstart,1,'first');
nstop=find(kgrid.t_array<tstop,1,'last');

%*****************************************************************%
% Second method, uses p and v to compute surfacic forces
%*****************************************************************%

ux=sensor_data.ux_non_staggered;
uy=sensor_data.uy_non_staggered;
uz=sensor_data.uz_non_staggered;
pr=sensor_data.p;

presh=reshape(pr,kgrid.Nx,kgrid.Ny,11,kgrid.Nt);
uxresh=reshape(ux,kgrid.Nx,kgrid.Ny,11,kgrid.Nt);
clear ux
uyresh=reshape(uy,kgrid.Nx,kgrid.Ny,11,kgrid.Nt);
clear uy
uzresh=reshape(uz,kgrid.Nx,kgrid.Ny,11,kgrid.Nt);
clear  uz

FxV=mean(presh.*uxresh,4);
clear uxresh
FyV=mean(presh.*uyresh,4);
clear uyresh
FzV=mean(presh.*uzresh,4);
clear uzresh
F=sqrt(FxV.^2+FyV.^2+FzV.^2);
Fx=squeeze(FxV(:,:,round(size(FxV,3)/2)));
Fy=squeeze(FyV(:,:,round(size(FxV,3)/2)));
Fz=squeeze(FzV(:,:,round(size(FxV,3)/2)));
xax=kgrid.x_vec-min(kgrid.x_vec);
yax=kgrid.y_vec;

save('Data\ARF_L7_4_17p8mFoc_4MHz.mat','Fx','Fy','Fz','sensor','kgrid','medium','transducer','tone_burst_freq','xax','yax','-v7.3');

%% plots

figure('color','w');
subplot(1,2,1);
imagesc(yax*1e3,xax*1e3,Fx);xlabel('azimuth [mm]');
ylabel('depth [mm]');title('F_x');
subplot(1,2,2);
imagesc(yax*1e3,xax*1e3,Fy);xlabel('azimuth [mm]');
ylabel('depth [mm]');title('F_y');

% Contour plots
figure('color','w');
surfVal=[-20 -10 -6];
for jj=1:length(surfVal)
    subplot(1,3,jj);
    p = patch(isosurface(kgrid.y_vec,kgrid.x_vec,(-5:5)*kgrid.dz,db(F/max(F(:))),surfVal(jj)));
    isonormals(kgrid.y_vec,kgrid.x_vec,(-5:5)*kgrid.dz,db(F/max(F(:))),p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3); 
    axis tight
    camlight 
    lighting gouraud
    title(['Normalized F isosurface at ',num2str(surfVal(jj),'%.1f'),' dB']);
end