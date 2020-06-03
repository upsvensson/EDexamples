% EDexample_Cube_convergenceDC.m


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
soudist = 1e9;
sources = soudist/sqrt(3)*[1 1 1];
Sindata = struct('coordinates',sources);
receivers = [0 0 0.50001;0 0 -0.50001];
Rindata = struct('coordinates',receivers);
controlparameters = struct('frequencies',[0 100 1000]);
controlparameters.Rstart = soudist;
controlparameters.difforder = 30;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);

ngvec = [24 32 48 64].';
ncases = size(ngvec,1);

nreceivers = size(Rindata.coordinates,1);
nfrequencies = prod(size(controlparameters.frequencies));
allres = zeros(nfrequencies,nreceivers,ncases);
for ii = 1:ncases

    controlparameters.ngauss = ngvec(ii);
    filehandlingparameters.filestem = [filestem,int2str(ngvec(ii))];

%     EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and present the results
    
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'''])

    tfdirect = tfdirect*soudist;
    tfgeom = tfgeom*soudist;
    tfdiff = tfdiff*soudist;
    tfinteqdiff = tfinteqdiff*soudist;

    tftot = (tfdirect + tfgeom + tfdiff + tfinteqdiff);

    allres(:,:,ii) = tftot;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Illustrate the linear extrapolation, for 0 Hz

xvec = 1./ngvec.^2;xvecplot = [0 xvec(1)].';
pfront = polyfit(xvec,squeeze(allres(1,1,:)),1);
prear = polyfit(xvec,squeeze(allres(1,2,:)),1);
xvecplot = [0 xvec(1)].';
yvecfrontplot = polyval(pfront,xvecplot);
yvecrearplot = polyval(prear,xvecplot);

figure
plot(xvec,squeeze(allres(1,1:2,:)),'o',xvecplot,yvecfrontplot,xvecplot,yvecrearplot)
xlabel('1/ngauss squared [-]')
h = legend('Computed, front','Computed, rear','Extrapol. front', 'Extrapol, rear');
set(h,'Location','best')
grid
ylabel('Sound pressure amplitude [-]')
title('Sound pressure amplitude at cube surface for the frequency 0 ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the relative error, for 0 Hz

A = [ngvec squeeze(allres(1,1:2,:)).'];
A = [A abs(A(:,2:3)-1)];
A = [A;0 pfront(2) prear(2) abs(pfront(2)-1) abs(prear(2)-1)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do several extrapolations, for each consecutive pair of points
% but only for the frontal result here

extrapolres = zeros(ncases,1);
for ii = 2:ncases
    pfront = polyfit(xvec(ii-1:ii),squeeze(allres(1,1,ii-1:ii)),1);
    extrapolres(ii) = pfront(2);
end
relerr = abs(   squeeze(allres(1,1,:)).'   -1);
relerrcompref = abs(   squeeze(   allres(1,1,:)).'   -allres(1,1,4));
relerrextrapol = abs(extrapolres-1);

figure
loglog(xvec,relerr,'-o',xvec,relerrcompref,'-o',xvec(2:end),relerrextrapol(2:end),'-*',xvec,0.1*xvec,'--',xvec,10*xvec.^2,'--')
grid
xlim([1e-4 2e-3])
h = legend('Computed results','Computed results with best computation as ref.','Extrapolated results','1/ngauss squared','1/ngauss to the fourth');
set(h,'Location','best')
title('Error convergence for the sound pressure at cube surface for f = 0')
ylabel('Relative error [-]')
xlabel('1/ngauss squared [-]')







   



    