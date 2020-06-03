% Lsp_Kessel_ir.m

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
% Sindata = struct('coordinates',[0.1 0.1 0.0001]);
% Rindata = struct('coordinates',[-0.3 0.3 -0.2]);
Sindata = struct('coordinates',[0 0 0.0001]);
Rindata = struct('coordinates',[0 0 2]);
controlparameters = struct('fs',48000);
controlparameters.difforder = 2;
controlparameters.saveindividualfirstdiff = 1;
controlparameters.savealldifforders = 1;
filehandlingparameters = struct('outputdirectory',[infilepath,filesep,'results2']);
filehandlingparameters.filestem = filestem;
filehandlingparameters.savelogfile = 1;
filehandlingparameters.showtext = 2;

EDmain_convex_time(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
nreceivers = size(Rindata.coordinates,1);
nsources = size(Sindata.coordinates,1);

eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_ir.mat'''])
if controlparameters.difforder > 1
    eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_irhod.mat'''])

    ncells = size(irhod,2);
    ntotlength = size(irhod{2},1);
    allirtots = zeros(ntotlength,ncells+1);
    for ii = 2:ncells    
       allirtots(:,ii+1) = irhod{ii}; 
    end
end

ndir = length(irdirect);
allirtots(1:ndir,1) = irdirect + irgeom;
irdiff_matrix = irdiff{1,1}.irvectors;
allirtots(1:ndir,2:1+size(irdiff_matrix,2)) = irdiff_matrix;

cumsumirtots = cumsum(allirtots.').';

tvec = 1/controlparameters.fs*[0:ntotlength-1].';

figure(1)
h = plot(tvec*1e3,sum(allirtots,2),'-o');
set(h(1),'LineWidth',1,'MarkerSize',3);
g = get(h(1),'Parent');
set(g,'FontSize',14)
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
grid
axis([5 8 -0.1 0.1])
g = legend('Direct sound (w spec. refl.) and 1st-order diffr.');
set(g,'FontSize',14)


tvec = 1/controlparameters.fs*[0:ntotlength-1].';

figure(2)
h = plot(tvec*1e3,allirtots(:,3),'-o');
set(h(1),'LineWidth',1,'MarkerSize',3);
g = get(h(1),'Parent');
set(g,'FontSize',14)
g = xlabel('Time   [ms]');
set(g,'FontSize',14)
g = ylabel('Impulse response   [-]');
set(g,'FontSize',14)
grid
axis([5 12 -0.001 0.003])
g = legend('Second-order diffr.');
set(g,'FontSize',14)

% figure(3)
% h = plot(tvec*1e3,allirtots(:,4),'-o');
% set(h(1),'LineWidth',1,'MarkerSize',3);
% g = get(h(1),'Parent');
% set(g,'FontSize',14)
% g = xlabel('Time   [ms]');
% set(g,'FontSize',14)
% g = ylabel('Impulse response   [-]');
% set(g,'FontSize',14)
% grid
% axis([5 12 -0.0025 0.0005])
% g = legend('Third-order diffr.');
% set(g,'FontSize',14)


nfft = 4096;
fvec = controlparameters.fs/nfft*[0:nfft/2-1];
F = fft(cumsumirtots,nfft);

figure(4)
h = semilogx(fvec,20*log10(abs(F(1:nfft/2,:))));
for ii = 1:3
   set(h(ii),'LineWidth',2) 
end
g = get(h(1),'Parent');
set(g,'FontSize',14)
grid
g = xlabel('Frequency   [Hz]');
set(g,'FontSize',14)
g = ylabel('Frequency response magnitude   [dB]');
set(g,'FontSize',14)
g = legend('GA only','Incl. diffr.1','Incl. diffr.2','Incl. diffr.3');
% axis([20 20000 -8 4])








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
