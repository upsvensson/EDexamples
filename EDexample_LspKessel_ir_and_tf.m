% EDexample_LspKessel_ir_and_tf.m

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
Sindata = struct('coordinates',[0 0 0.01]);
Rindata = struct('coordinates',[0 0 2]);
controlparameters = struct('fs',48000);
controlparameters.difforder = 3;
controlparameters.savealldifforders = 1;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results']);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;
filehandlingparameters.showtext = 1;

EDmain_convex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

nfft = 4096;
fvec = controlparameters.fs/nfft*[0:nfft/2-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

controlparameters.difforder = 15;
controlparameters.ngauss = 24;
fvec_tf = fvec(1:100);
controlparameters.frequencies = fvec_tf;
EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
nreceivers = size(Rindata.coordinates,1);
nsources = size(Sindata.coordinates,1);

eval(['load ''',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_ir.mat'''])
eval(['load ''',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_irhod.mat'''])

ncells = size(irhod,2);
ntotlength = size(irhod{2},1);
allirtots = zeros(ntotlength,ncells+1);
for ii = 2:ncells    
   allirtots(:,ii+1) = irhod{ii}; 
end

ndir = length(irdirect);
allirtots(1:ndir,1) = irdirect + irgeom;
allirtots(1:ndir,2) = irdiff;

cumsumirtots = cumsum(allirtots.').';

eval(['load ''',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'''])
eval(['load ''',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])

tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;

tvec = 1/controlparameters.fs*[0:ntotlength-1];

figure(1)
h = plot(tvec*1e3,cumsumirtots,'-');
% set(h(1),'MarkerSize',2);
set(h(2),'LineWidth',2);
set(h(3),'LineWidth',2);
set(h(4),'LineWidth',2);
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
axis([6 12 -0.02 0.003])
g = legend('GA','1st-order diffr.','2nd-order diffr.','3rd-order diffr.');
set(g,'FontSize',14,'Location','SouthEast')


F = fft(cumsumirtots,nfft);

figure(2)
h = semilogx(fvec,20*log10(abs(F(1:nfft/2,:))),fvec_tf,20*log10(abs(tftot)),'*');
for ii = 1:4
   set(h(ii),'LineWidth',2) 
end
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Frequency   [Hz]');
set(g,'FontSize',14)
g = ylabel('Frequency response magnitude   [dB]');
set(g,'FontSize',14)
g = legend('GA only','Incl. diffr.1','Incl. diffr.2','Incl. diffr.3','Incl. diffr.15');
set(g,'FontSize',14,'Location','SouthEast')
axis([20 20000 -8 4])








% % 
% irtot = irdirect + irgeom + irdiff;
% nold = length(irtot);
% nnew = length(hodir);
% if nnew > nold
%    irtot = [irtot;zeros(nnew-nold,nreceivers,nsources)]; 
% end
% irtot(1:nnew,:,:) = irtot(1:nnew,:,:) + hodir;
% 
% nfft = 4096;
% fvec = controlparameters.fs/nfft*[0:nfft/2-1];
% F = fft(squeeze(irtot),nfft);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% controlparameters.frequencies = fvec(1:50);
% EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
% 
% eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'''])
% eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'''])
% % 
% tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
% 
% semilogx(fvec,20*log10(abs((F(1:nfft/2,:)))),'-',fvec(1:50),20*log10(abs(squeeze(tftot))),'*')
% grid
% % 
