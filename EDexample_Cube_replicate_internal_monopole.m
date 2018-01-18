% EDexample_Cube_replicate_internal_monopole
% Temporary version

% Cube_internalmono1_32g.m

mfile = mfilename('fullpath');
[infilepath,filestem] = fileparts(mfile);

corners = [     -0.5000   -0.500   -0.500
    0.5000   -0.500   -0.500
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
geofiledata = struct('corners',corners);
geofiledata.planecorners = planecorners;
planedata = EDreadgeomatrices(geofiledata.corners,geofiledata.planecorners);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controlparameters

controlparameters = struct('frequencies',0.1);
nfrequencies = length(controlparameters.frequencies);
controlparameters.ngauss = 32;
controlparameters.difforder = 15;

envdata = struct('cair',344);
envdata.rhoair = 1.21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sources

patchpoints = [];
patchcorners = [];
patchareas = [];
patchnvecs = [];
n1 = 6;
n2 = 6;
planenvecs = planedata.planeeqs(:,1:3);
for ii = 1:6
   [patchcornercoords,patchmidcoords,areas] = EDdivrect(corners(planecorners(ii,1),:),corners(planecorners(ii,2),:),corners(planecorners(ii,3),:),corners(planecorners(ii,4),:),n1,n2);     
   nvec = planenvecs(ii,:);
   patchmidcoords = patchmidcoords + 0.0001*nvec(ones(n1*n2,1),:);
   patchpoints = [patchpoints;patchmidcoords];
   patchcorners = [patchcorners;patchcornercoords];
   patchareas = [patchareas;areas];
   patchnvecs = [patchnvecs;nvec(ones(n1*n2,1),:)];
end

k = 2*pi*controlparameters.frequencies/envdata.cair;
if nfrequencies > 1
   error('ERROR: sorry, only one freq. in this script') 
end

internalmonopole = [0.3 0.3 0.3];
quadweights = [64 25 25 25 25 40 40 40 40]/81/4;
quadshift = sqrt(3)/sqrt(5);
umono_normal_patch = 0;

% close all
% figure(1)
% Gauss point 1 = center point
patchquadraturepoints = patchpoints;
% plot3(patchquadraturepoints(1:3,1),patchquadraturepoints(1:3,2),patchquadraturepoints(1:3,3),'*')
% hold on
[distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
pmono_patch = exp(-1i*k*distances)./distances;
umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
umono_normal_patch =  quadweights(1)*umono_radial_patch.*cosfi;

% Gauss points 2-5 = from center point, towards corners
for ii = 1:4
    patchquadraturepoints = patchpoints + quadshift*(patchcorners(:,[1:3] + 3*(ii-1)) - patchpoints);
% plot3(patchquadraturepoints(1:3,1),patchquadraturepoints(1:3,2),patchquadraturepoints(1:3,3),'o')

    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch = umono_normal_patch + quadweights(ii+1)*umono_radial_patch.*cosfi;    
end

% Gauss points 6-9 = from center point, towards midpoint
for ii = 1:4
    if ii < 4
        patchedgemidpoint = 0.5*( patchcorners(:,[1:3] + 3*(ii-1)) + patchcorners(:,[1:3] + 3*(ii)) );
    else
        patchedgemidpoint = 0.5*( patchcorners(:,10:12) + patchcorners(:,1:3) );        
    end
    patchquadraturepoints = patchpoints + quadshift*(patchedgemidpoint - patchpoints);
%     plot3(patchquadraturepoints(1:3,1),patchquadraturepoints(1:3,2),patchquadraturepoints(1:3,3),'s')

    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch = umono_normal_patch + quadweights(ii+5)*umono_radial_patch.*cosfi;    
end

% plot3(patchcorners(1:3,1),patchcorners(1:3,2),patchcorners(1:3,3),'+')
% plot3(patchcorners(1:3,4),patchcorners(1:3,5),patchcorners(1:3,6),'+')
% plot3(patchcorners(1:3,7),patchcorners(1:3,8),patchcorners(1:3,9),'+')
% 
% axis equal
% 
% pause

Sindata = struct('coordinates',patchpoints);
nsources = size(Sindata.coordinates,1);
Sindata.sourceamplitudes = patchareas.*umono_normal_patch...
     .*1i*2*pi*controlparameters.frequencies*envdata.rhoair/4/pi;
Sindata.doaddsources = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receivers

rdist = 3;
phivec = [0:0.1:2*pi].';
receivers = rdist*[cos(phivec) sin(phivec) 0*phivec];
nreceivers = size(receivers,1);
receivers = receivers + internalmonopole(ones(nreceivers,1),:);

Rindata = struct('coordinates',receivers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filehandling

filehandlingparameters = struct('outputdirectory',infilepath);
filehandlingparameters.filestem = filestem;
filehandlingparameters.showtext = 1;

% if exist(['./results/',filehandlingparameters.filestem,'_tfinteq.mat'],'file') ~= 2
%     EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and present the results
    
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'])

tftot = tfinteqdiff + tfdiff + tfdirect + tfgeom;

% tftot = tftot*factor;

figure(1)
plot(phivec*180/pi,rdist*abs([tftot.' ]),'-o')
xlabel(['Receiver angle   [deg.]'])
ylabel(['Sound pressure amp., times distance   [-]'])
% legend('Computed field','Computed field with GA and diff1')
% ylim([0.95 1.05])
grid


% tftot = tfdirect + tfgeom + tfdiff + tfinteqdiff;
% semilogx(controlparameters.frequencies,20*log10(abs(tftot)),'-o')
% grid

