% EDexample_Squareplate_edgepointsperwavelength.m
%
% This script computes the sound pressure in two receiver positions, one at 
% the frontal side and one at the rear side, of a square plate, for plane
% wave incidence. Twenty frequencies, and nine different
% numbers of gausspoints per edge are run. The relative errors for all of
% these are computed and plotted. The ref. results are, for each frequency,
% either the best computed result, or the extrapolated value from the two best results.

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     -0.5000   -0.500   0
    0.5000   -0.50   0
    0.5000    0.5000   0
   -0.5000    0.5000   0];

planecorners = [   1     2 3 4; 4 3 2 1];

geofiledata = struct('corners',corners,'planecorners',planecorners);
soudist = 1e9;
sources = soudist/sqrt(3)*[1 1 1];
Sindata = struct('coordinates',sources);
receivers = [0 0 -0.001;0 0 0.001];
Rindata = struct('coordinates',receivers);
controlparameters = struct('frequencies',logspace(log10(500),log10(5000),20));
controlparameters.Rstart = soudist;
controlparameters.difforder = 15;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.showtext = 0;

ngvec = [12 14 16 20 24 32 48 64 72].';
ncases = size(ngvec,1);

nreceivers = size(Rindata.coordinates,1);
nfrequencies = prod(size(controlparameters.frequencies));

allres = zeros(nfrequencies,nreceivers,ncases);
for ii = 1:ncases
    disp(['   Case ',int2str(ii),' of ',int2str(ncases)])
    controlparameters.ngauss = ngvec(ii);
    filehandlingparameters.filestem = [filestem,int2str(ngvec(ii))];

     EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load and present the results
    
    eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
    eval(['load ',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tf.mat'])

    tfdirect = tfdirect*soudist;
    tfgeom = tfgeom*soudist;
    tfdiff = tfdiff*soudist;
    tfinteqdiff = tfinteqdiff*soudist;

    tftot = ( tfdirect + tfgeom + tfdiff + tfinteqdiff);

    allres(:,:,ii) = tftot;
end

allresR = squeeze(allres(:,1,:));
allresF = squeeze(allres(:,2,:));

% Generate n-1 extrapolated results from each consecutive pair of data
% points.

xvec = 1./ngvec.^2;

extrapolresR = zeros(size(allresR));
extrapolresF = zeros(size(allresF));

nfreq = size(allres,1);
refresR = zeros(nfreq,1);
refresF = zeros(nfreq,1);
for ii = 1:nfreq
    for jj = 2:ncases
        p = polyfit(xvec(jj-1:jj).',allresR(ii,jj-1:jj),1);
        extrapolresR(ii,jj) = polyval(p,0);
        p = polyfit(xvec(jj-1:jj).',allresF(ii,jj-1:jj),1);
        extrapolresF(ii,jj) = polyval(p,0);
    end
end

refresR = extrapolresR(:,end);
refresF = extrapolresF(:,end);

relerrextraR = abs( (extrapolresR-refresR(:,ones(1,ncases)))./refresR(:,ones(1,ncases)) );
relerrextraF = abs( (extrapolresF-refresF(:,ones(1,ncases)))./refresF(:,ones(1,ncases)) );

elsize = (1./ngvec).';
lambdavec = (envdata.cair./controlparameters.frequencies).';

elperlambda = lambdavec(:,ones(1,ncases))./elsize(ones(nfreq,1),:);

relerrextraR = reshape(relerrextraR,ncases*nfrequencies,1);
relerrextraF = reshape(relerrextraF,ncases*nfrequencies,1);
elperlambda = reshape(elperlambda,ncases*nfrequencies,1);

figure(1)
h = loglog(elperlambda,relerrextraR,'bo',elperlambda,relerrextraF,'r*');
grid
legend('Rec. without GA components','Rec. with GA components')
xlim([0.7 70])
xlabel('Edge points per wavelength   [-]')
ylabel('Rel. error   [-]')
title('Rel. error for extrapol. result, plane wave incidence onto a square plate (ref. = best result)')
ylim([1e-7 0.1])











