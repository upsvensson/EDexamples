% EDexample_LspKessel_minimal.m
%
% This example calculates the frequency response of a simplified
% loudspeaker. The loudspeaker is called "LspKessel" because the dimensions
% are taken from the paper "A simple theory of cabinet edge diffraction" by
% John Vanderkooy, in JAES 39, pp. 923-933, 1991. That paper presented BEM
% calculations done by R.T. Kessel.
%
% A single point source is modeled, at a distance of 0.01 mm from the front
% baffle. A receiver position is at a distance of 10 m, on-axis. These
% positions were also taken from the paper by Vanderkooy.
%
% One purpose with this example is to demonstrate a minimal problem, to get
% started. A number of parameters are not set in this file, but are given default
% values, including:
%   - diffraction order 15
%   - number of edge points = number of Gauss-Legendre quadrature points = 16
%     This is the number per longest edge
%     You will see a message printed out that tells you that 16 edge points
%     could make it possible to calculate up to 3 kHz. For higher
%     frequencies, you need more points (on average ≈ 2.8 points per
%     wavelength)
%   - speed of sound: 344 m/s
%
% If you have set the path to the EDtoolbox folder, you should be able to
% run this file and the calculations will get started. A folder called
% "results" will be created inside the folder where this m-file is stored,
% and results files will end up in that "results" folder.
% 
% The last part of this script, below the "%%%%%%%%%%%%%%%%%%%%" line,
% presents the results in a diagram.

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     -0.2000   -0.4400   -0.3200
    0.2000   -0.4400   -0.3200
    0.2000    0.2000   -0.3200
   -0.2000    0.2000   -0.3200
   -0.2000   -0.4400         0
    0.2000   -0.4400         0
    0.2000    0.2000         0
   -0.2000    0.2000         0];

planecorners = [   1     4     3     2
     5     6     7     8
     1     2     6     5
     3     4     8     7
     2     3     7     6
     1     5     8     4];

geofiledata = struct('corners',corners,'planecorners',planecorners);
Sindata = struct('coordinates',[0 0 0.00001]);
Rindata = struct('coordinates',[0 0 10]);
controlparameters = struct('frequencies',logspace(log10(50),log10(3000),100));
filehandlingparameters = struct('outputdirectory',infilepath);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;

EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'])

tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
semilogx(controlparameters.frequencies,20*log10(abs(tftot)),'-o')
xlabel('Frequency   [Hz]')
ylabel('TF magnitude re. 1m   [dB]')
title('Frequency response of the Kessel loudspeaker')
grid

