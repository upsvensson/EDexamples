% EDexample_Cube_HODcost_ng.m

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

phivec = [10:1:85].';
theta = 60;
sdist = 10;
rdist = 11;

sources = sdist*[sind(theta)*cosd(phivec) sind(theta)*sind(phivec) cosd(theta)*ones(size(phivec))]; 
receivers = rdist*[sind(theta)*cosd(phivec) sind(theta)*sind(phivec) cosd(theta)*ones(size(phivec))]; 

nsources = size(sources,1);

frequencies = linspace(500,20000,20);
frequencies = frequencies(1);

% nSRvec = [1 4 7 10 13 16 19];
% ncases = length(nSRvec);

% diffordervec = [2:10];
% ncases = length(diffordervec);

ngaussvec = [8 16 20 24 32 40];
ncases = length(ngaussvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Give values to the input structs

geofiledata = struct('firstcornertoskip',1e6);
geofiledata.corners = corners;
geofiledata.planecorners = planecorners;

envdata = struct('cair',344);
controlparameters = struct('frequencies',frequencies);
controlparameters.directsound = 1;
controlparameters.difforder = 10;
controlparameters.ngauss = 32;
filehandlingparameters = struct('showtext',0);
filehandlingparameters.outputdirectory = [infilepath,filesep,'results'];
filehandlingparameters.savecadgeofile = 0;
filehandlingparameters.saveeddatafile = 0;
filehandlingparameters.saveSRdatafiles = 0;
filehandlingparameters.savelogfile = 0;
filehandlingparameters.savesubmatrixdata = 0;
filehandlingparameters.saveinteqsousigs = 0;
filehandlingparameters.savediff2result = 0;
filehandlingparameters.showtext = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timingres = zeros(ncases,5);

Sindata = struct('coordinates',sources(1,:));
Sindata.doaddsources = 1;
Rindata = struct('coordinates',receivers(1,:));

for ii = 1:ncases

    disp(['   Case ',int2str(ii)])
 
    controlparameters.ngauss = ngaussvec(ii);
    filehandlingparameters.filestem = [filestem,'_',int2str(ii)];
    
      EDmain_convexESIE(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters);        

    eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat timingstruct'])

    timingres(ii,:) = timingstruct.integralequation;
    
end

timingres(:,4) = timingres(:,4)/(controlparameters.difforder-2);

plot(ngaussvec.^3,timingres(:,[2 4]),'-o')
grid
xlabel('n_{gauss}^3   [-]')
ylabel('Calculation time   [s]')
title('Calculation time for higher-order diffraction, one frequency, one S and one R')
h = legend('Compute H-matrix','q, per iteration');
set(h,'location','best')









