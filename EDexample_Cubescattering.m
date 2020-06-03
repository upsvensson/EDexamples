% EDexample_Cubescattering.m
%
% This example calculates the response around a cube for an incidendent plane wave.
% The case is the one in the paper "Benchmark cases in 3D diffraction with
% different methods" by U.P. Svensson, H. Brick, and J. Forssén, in proc.
% of Forum Acusticum 2014, 7-12 Sept., Krakow, Poland; as well as in 
% "A hybrid method combining the edge source integral equation and the 
% boundary element method for scattering problems" by S.R. Martin, U. P.
%  Svensson, J. Slechta, and J. O. Smith, Proc. Mtgs. Acoust. 26, 015001
%  (2016), https://doi.org/10.1121/2.0000226.
%
% True plane waves can not be generated with version 0.1 of the EDtoolbox,
% but a point source at a huge distance works fine in most cases. The 
% response must then be scaled up accordingly. The receivers are distributed
% around the cube along a circle in steps of 0.01 radians.
%
% The number of gauss points set here is 48, even if a higher value would
% still give reasonable calculation times. The lower value of 48
% illustrates better the problem for certain receiver positions.
% 
% Speed of sound: 344 m/s
%
% If you have set the path to the EDtoolbox folder, you should be able to
% run this file and the calculations will get started. A folder called
% "results" will be created inside the folder where this m-file is stored,
% and results files will end up in that "results" folder.
% 
% The last part of this script, below the "%%%%%%%%%%%%%%%%%%%%" line,
% presents the results in a diagram. This particular example illustrates
% that the ESIE approach has a singularity problem for certain receiver
% positions. For 8 directions along this receiver circle, the higher-order
% diffraction component has a numerically challenging singularity. It shows
% up as tiny little wiggles. There is one more little glitch, for the
% receiver angle 0 degrees. The 0.1 version of the EDtoolbox has a 
% numerical problem when a direct sound or a specular reflection hits
% exactly at a corner.
% The model is also plotted in a separate window.

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     -0.5000   -0.500   -0.500
    0.5000   -0.50   -0.500
    0.5000    0.5000   -0.500
   -0.5000    0.5000   -0.500
   -0.5000   -0.500         0.5
    0.5000   -0.500         0.5
    0.5000    0.5000         0.5
   -0.5000    0.5000         0.5];

planecorners = [   1     4     3     2
     5     6     7     8
     1     2     6     5
     3     4     8     7
     2     3     7     6
     1     5     8     4];

geofiledata = struct('corners',corners,'planecorners',planecorners);
soudist = 1e6;
sources = soudist/sqrt(3)*[1 1 1];
Sindata = struct('coordinates',sources);
anglevec = [0:0.01:2*pi].';
receivers = [cos(anglevec) sin(anglevec) 0*anglevec];
Rindata = struct('coordinates',receivers);
controlparameters = struct('frequencies',500);
controlparameters.ngauss = 48;
controlparameters.Rstart = soudist;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
%filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;

EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'''])

tftot = soudist*(tfdirect + tfgeom + tfdiff + tfinteqdiff);

figure
plot(anglevec*180/pi,20*log10(abs(tftot)),'-o','MarkerSize',4)
ylim([-6 4])
xlabel('Receiver angle   [deg.]')
ylabel('TF magnitude re. 1m   [-]')
title('Total field around a 1m cube, at 500 Hz, for plane wave incidence')
grid

figure
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,2)

