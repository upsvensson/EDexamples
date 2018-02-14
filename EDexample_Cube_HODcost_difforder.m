% EDexample_Cube_HODcost_difforder.m

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

phivec = [5:1:85].';
theta = 60;
sdist = 10;
rdist = 11;

sources = sdist*[sind(theta)*cosd(phivec) sind(theta)*sind(phivec) cosd(theta)*ones(size(phivec))]; 
receivers = rdist*[sind(theta)*cosd(phivec) sind(theta)*sind(phivec) cosd(theta)*ones(size(phivec))]; 

nsources = size(sources,1);

frequencies = linspace(500,20000,20);
frequencies = frequencies(1);

nSRvec = [1 20 30 40 50 60 70 80];
ncases = length(nSRvec);

diffordervec = [2:10];
ncases = length(diffordervec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Give values to the input structs

geofiledata = struct('firstcornertoskip',1e6);
geofiledata.corners = corners;
geofiledata.planecorners = planecorners;

envdata = struct('cair',344);
controlparameters = struct('frequencies',frequencies);
controlparameters.directsound = 1;
controlparameters.difforder = 1;
controlparameters.ngauss = 64;
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

    Sindata = struct('coordinates',sources(5,:));
    Rindata = struct('coordinates',receivers(5,:));

for ii = 1:ncases
    disp(['   Case ',int2str(ii),' of ',int2str(ncases)])
    filehandlingparameters.filestem = [filestem,'_',int2str(ii)];
    controlparameters.difforder = diffordervec(ii);
    
       EDmain_convexESIE(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters);        

    eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat timingstruct'])

    timingres(ii,:) = timingstruct.integralequation;
    
end

plot(diffordervec,timingres(:,2:end),'-o')
grid
xlabel('Diffraction order   [-]')
ylabel('Calculation time   [s]')
title('Calculation time for higher-order diffraction, one frequency, one S, one R, ngauss = 16')
h = legend('Compute H-matrix','q_0','q, via iteration','p_{receiver}');
set(h,'location','best')









