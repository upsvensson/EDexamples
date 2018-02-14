% EDexample_LspKessel_minimal.m
%
% This example calculates the frequency response of a simplified
% loudspeaker. The loudspeaker is called "Kessel" because the dimensions
% are taken from the paper "A simple theory of cabinet edge diffraction" by
% John Vanderkooy, in JAES 39, pp. 923-933, 1991, which presented BEM
% calculations done by R.T. Kessel (for his MSc thesis at the University of
% Waterloo, Canada.
%
% A single point source is modeled, at a distance of 0.01 mm from the front
% baffle. A receiver position is at a distance of 1 m, on-axis. These
% positions were also taken from the paper by Vanderkooy - except that here
% we choose a 1 m distance instead of the 10 m distance in Vanderkooy's
% paper.
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
% run this script, the calculations will be carried out, and you will see
% result plots. A folder called 'results' will be created inside the folder
% where this m-file is stored, and results files will end up in that
% 'results' folder.
% 
% The last part of this script, below the "%%%%%%%%%%%%%%%%%%%%" line,
% presents the results in a diagram and plots the model.

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
Rindata = struct('coordinates',[0 0 1]);
controlparameters = struct('frequencies',linspace(50,3000,100));
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.filestem = filestem;

EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'])

tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;

figure
semilogx(controlparameters.frequencies,20*log10(abs(tftot)),'-o')
xlabel('Frequency   [Hz]')
ylabel('TF magnitude re. 1m   [dB]')
title('Frequency response of the Kessel loudspeaker, at 1m distance')
axis([50 5000 0 10])
grid

figure
eddatafile = [infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_eddata.mat'];
EDplotmodel(eddatafile,3)

