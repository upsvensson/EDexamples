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

diffordervec = [3 5 8];
ncases = length(diffordervec);

ngaussvec = [16 20 24 28 32 36 40 48 56 64 72 80 88 96 104 112 120 128 136 142 148];

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
% filehandlingparameters.suppressresultrecycling = 1;
filehandlingparameters.savecadgeofile = 0;
filehandlingparameters.saveeddatafile = 0;
filehandlingparameters.saveSRdatafiles = 0;
filehandlingparameters.savelogfile = 1;
filehandlingparameters.savesubmatrixdata = 0;
filehandlingparameters.saveinteqsousigs = 0;
filehandlingparameters.savediff2result = 0;
filehandlingparameters.showtext = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngausscases = length(ngaussvec);
timingres = zeros(ngausscases,3);

Sindata = struct('coordinates',sources(5,:));
Rindata = struct('coordinates',receivers(5,:));

tcase = zeros(ncases,2);

for jj = 1:ngausscases
    disp(['   ngausscase ',int2str(jj),' of ',int2str(ngausscases)])
    controlparameters.ngauss = ngaussvec(jj);
    for ii = 1:ncases
        disp(['      Case ',int2str(ii),' of ',int2str(ncases)])
        filehandlingparameters.filestem = [filestem,'_',int2str(ii),'_',int2str(jj)];
        controlparameters.difforder = diffordervec(ii);

        EDmain_convexESIE(geofiledata,Sindata,Rindata,envdata,controlparameters,filehandlingparameters);        

        eval(['load ''',filehandlingparameters.outputdirectory,filesep,filehandlingparameters.filestem,'_tfinteq.mat''','  timingstruct'])    
        tcase(ii,:) = timingstruct.integralequation([2 4]);

    end
    timingres(jj,1) = mean(tcase(:,1));
    timingres(jj,2:3) = polyfit(diffordervec.',tcase(:,2),1);    
end    
   
iv = find(ngaussvec>=56);

% p = polyfit(ngaussvec(iv).',1e6*timingres(iv,1)./(ngaussvec(iv).'.^{3.5}),1)
% const1 = p(2);
const1 = mean( 1e6*timingres(iv,1)./(ngaussvec(iv).'.^3.5) );

% p = polyfit(ngaussvec(iv).',1e6*timingres(iv,2)./(ngaussvec(iv).'.^{3.5}),1)
% const2 = p(2);
const2 = mean( 1e6*timingres(iv,2)./(ngaussvec(iv).'.^3.5) );


figure(1)
h = loglog(ngaussvec.',1e6*timingres(:,2)./(ngaussvec.'.^3.5),'b-o',ngaussvec.',...
    1e6*timingres(:,1)./(ngaussvec.'.^3.5),'r-*',ngaussvec(iv),const2*ones(size(iv)),'b-',...
    ngaussvec(iv),const1*ones(size(iv)),'r-');
set(h(3),'LineWidth',2);
set(h(4),'LineWidth',2);
grid
xlabel('ngauss, per edge   [-]')
ylabel('Time   [s]')
h = legend('t_{iteration}*1e6/ng^{3.5}','t_{setup}*1e6/ng^{3.5}');
set(h,'FontSize',14)
title('Calculation times for two HOD stages')
ylim([0.6 2])
xlim([10 200])


figure(2)
h = loglog(ngaussvec.',timingres(:,2),'b-o',ngaussvec.',...
    timingres(:,1),'r-*');
grid
xlabel('ngauss, per edge   [-]')
ylabel('Time   [s]')
h = legend('t_{iteration}','t_{setup}');
set(h,'FontSize',14,'Location','best')
title('Calculation times for two HOD stages')
% ylim([0.6 2])
xlim([10 200])


% figure(2)
% h = loglog(ngaussvec.',timingres(:,1)./(ngaussvec.'.^3),'-o',ngaussvec.',10*timingres(:,1)./(ngaussvec.'.^4),'-o',ngaussvec.',2*timingres(:,1)./(ngaussvec.'.^3.5),'-o');%,ngaussvec,0.1e-4*ngaussvec.^3,ngaussvec,0.1e-6*ngaussvec.^4);
% grid
% xlabel('ngauss, per edge   [-]')
% ylabel('Time   [s]')
% legend('Normalized against ng^3','Normalized against ng^4')
% title('Timingres(1)')
% set(h(2),'LineWidth',2)


% plot(diffordervec,timingres(:,2:end),'-o')
% grid
% xlabel('Diffraction order   [-]')
% ylabel('Calculation time   [s]')
% title('Calculation time for higher-order diffraction, one frequency, one S, one R, ngauss = 16')
% h = legend('Compute H-matrix','q_0','q, via iteration','p_{receiver}');
% set(h,'location','best')









