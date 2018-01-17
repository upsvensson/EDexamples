% EDexample_LspCylinder16.m
%
% This example calculates the frequency response of a loudspeaker shaped
% as a polygonal cylinder, with 16 edges. This shape is chosen
% because it maximizes the effect of the edge diffraction on the on-axis
% response, as shown in the paper "Direct radiator loudspeaker enclosures"
% by Harry F. Olson, in JAES 17, pp. 22-29. The dimensions were chosen as
% in Olson's paper: cylinder diameter and length: 2ft = 0.6096 meters.
% 
% A function is called, EDmakegeo_cylinder, which generates the corners and
% planes geometry.
%
% A single point source is modeled, at a distance of 0.01 mm from one of
% the flat surfaces of the cylinder. A receiver position is at a distance
% of 6.6 m, on-axis. This distance was not specified in the paper by Olson
% but information given on the expected dip frequencies indicate that the
% distance was 6.6 m.
%
% Some settings worth mentioning:
%  - diffraction order 50. This was chosen quite arbitrarily to a high
%  value. The more polygons, and the lower the frequency, the higher does
%  the diffraction order need to be.
%  - ngauss 24 gives at least 2.8 points per wavelength up to around 5 kHz.
%
% If you have set the path to the EDtoolbox folder, you should be able to
% run this file and the calculations will get started. A folder called
% "results" will be created inside the folder where this m-file is stored,
% and results files will end up in that "results" folder.
% 
% The last part of this script, below the "%%%%%%%%%%%%%%%%%%%%" line,
% presents the results in a diagram and plots the model in a separate window.

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

radius = 0.3048;
length = 2*radius;
[corners,planecorners,ncorners,radius] = EDmakegeo_cylinder(radius,length,16,'e',1,0,-radius);

geofiledata = struct('corners',corners,'planecorners',planecorners);
geofiledata.planecornertype = 'zero';
Sindata = struct('coordinates',[0 0 0.00001]);
Rindata = struct('coordinates',[0 0 6.6]);
controlparameters = struct('frequencies',linspace(50,4000,100));
controlparameters.ngauss = 24;
controlparameters.difforder = 50;
filehandlingparameters = struct('outputdirectory',infilepath);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;

EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'])

measdistance = norm(Sindata.coordinates-Rindata.coordinates)
tftot = (tfdirect + tfgeom + tfdiff + tfinteqdiff)*measdistance;

figure
semilogx(controlparameters.frequencies,20*log10(abs(tftot)),'-o')
xlabel('Frequency   [Hz]')
ylabel('TF magnitude re. 1m   [dB]')
title('Frequency response of a cylindrical loudspeaker')
grid
axis([100 4000 -10 15])

figure
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,1+2)
