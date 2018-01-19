% EDexample_Cube_replicatemonopole
%
% This example demonstrates how the surface of a rigid object (a cube) can
% be covered with a mesh of point sources, the source amplitudes of which
% will be given by the normal velocity generated by a fictive internal
% (inside the cube) free-field radiating monopole source, at the locations
% of those surface sources. When all of these sources are activated, they
% should together create a perfect spherical wave, centered at the location
% of the fictive internal source. The internal source could be anywhere
% inside the rigid object.
%
% The script also demonstrates how to run a set of calculations, with a
% range of surface mesh points, as well as a range of numbers of edge points.
%
% For each combination of number of surface mesh points, and edge points,
% results are saved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% First we specify the values that will be run in the for-loops
% ngaussvec are the numbers of edge points per edge
% nsurfelvec is the number of surface mesh divisions: "6" means that a 6*6
% surface mesh will be generated per cube face.

ngaussvec = [24 32 40];
nsurfelvec = [6 7 8 10 12 14 16];

frequency = 100;
receiverradius = 3;

% Below the position of the internal, fictive monopole is given. This
% script does not check that it is actually inside the cube.

internalmonopole = [0.3 0.3 0.3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A number of input parameters are set outside the for loops.

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
planenvecs = planedata.planeeqs(:,1:3);

% Controlparameters = calculation settings

controlparameters = struct('frequencies',frequency);
nfrequencies = length(controlparameters.frequencies);
controlparameters.difforder = 15;

envdata = struct('cair',344);
envdata.rhoair = 1.21;

% Receivers

phivec = [0:0.1:2*pi].';
receivers = receiverradius*[cos(phivec) sin(phivec) 0*phivec];
nreceivers = size(receivers,1);
receivers = receivers + internalmonopole(ones(nreceivers,1),:);

Rindata = struct('coordinates',receivers);

% Filehandling

filehandlingparameters = struct('outputdirectory',infilepath);
filehandlingparameters.showtext = 0;

% We will store the median value (across receiver positions) of the errors
medianerr_results = zeros(length(nsurfelvec),length(ngaussvec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here comes the for-loop over number of surface sources

k = 2*pi*controlparameters.frequencies/envdata.cair;

% A Gauss quadrature will be used to calculate the average normal particle
% velocity across each surface patch. A 3*3 quadrature is used.

quadweights = [64 25 25 25 25 40 40 40 40]/81/4;
quadshift = sqrt(3)/sqrt(5);

for kk = 1:length(nsurfelvec)

    % Sources
    %
    % This section defines the source positions. Each cube face is subdivided
    % into an n by n mesh. Each patch centre point is selected as the source
    % position.
    
    n1 = nsurfelvec(kk);
    n2 = n1;
    patchpoints = [];
    patchcorners = [];
    patchareas = [];
    patchnvecs = [];
    for ii = 1:6
       [patchcornercoords,patchmidcoords,areas] = EDdivrect(corners(planecorners(ii,1),:),corners(planecorners(ii,2),:),corners(planecorners(ii,3),:),corners(planecorners(ii,4),:),n1,n2);     
       nvec = planenvecs(ii,:);
       patchmidcoords = patchmidcoords + 0.0001*nvec(ones(n1*n2,1),:);
       patchpoints = [patchpoints;patchmidcoords];
       patchcorners = [patchcorners;patchcornercoords];
       patchareas = [patchareas;areas];
       patchnvecs = [patchnvecs;nvec(ones(n1*n2,1),:)];
    end

    % Now, compute the average normal particle velocity across each surface
    % patch with a 9-point quadrature.
    %
    % "umono_normal_patch" is this average value.
    
    umono_normal_patch = 0;

    % Gauss point 1 = patch center point
    patchquadraturepoints = patchpoints;
    [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
    pmono_patch = exp(-1i*k*distances)./distances;
    umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
    umono_normal_patch =  quadweights(1)*umono_radial_patch.*cosfi;

    % Gauss points 2-5 = from center point, towards corners
    for ii = 1:4
        patchquadraturepoints = patchpoints + quadshift*(patchcorners(:,[1:3] + 3*(ii-1)) - patchpoints);
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
        [distances,cosfi] = EDcalccosfi(internalmonopole,patchquadraturepoints,patchnvecs);
        pmono_patch = exp(-1i*k*distances)./distances;
        umono_radial_patch = pmono_patch/envdata.cair/envdata.rhoair.*(1 + 1./(1i*k*distances));
        umono_normal_patch = umono_normal_patch + quadweights(ii+5)*umono_radial_patch.*cosfi;    
    end
    
    % Finally, the Sindata input struct can be defined, with the source
    % positions and source amplitudes. The source amplitudes are the normal
    % particle velocity * patch area *jw *rhoair /(4*pi)
    Sindata = struct('coordinates',patchpoints);
    nsources = size(Sindata.coordinates,1);
    Sindata.sourceamplitudes = patchareas.*umono_normal_patch...
         .*1i*2*pi*controlparameters.frequencies*envdata.rhoair/4/pi;
    Sindata.doaddsources = 1;

    SS = ['S',int2str(n1)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Then the for-loop over number of edge points ("ngauss")
    
    for jj = 1:length(ngaussvec)
        GG = ['EG',int2str(ngaussvec(jj))];
        controlparameters.ngauss = ngaussvec(jj);
        filehandlingparameters.filestem = [filestem,'_',GG,'_',SS];
        disp(filehandlingparameters.filestem)

        % Here the actual calculations are done.
        % To develop further post-processing, after running through the
        % script once, you can comment out the line below and run the
        % script again.

        EDmain_convexESIE(geofiledata,Sindata,Rindata,struct,controlparameters,filehandlingparameters);

        % Load and analyze the results

        eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tfinteq.mat'])
        eval(['load ',infilepath,filesep,'results',filesep,filehandlingparameters.filestem,'_tf.mat'])

        tftot = tfinteqdiff + tfdiff + tfdirect + tfgeom;
        tftot = tftot*receiverradius;
        
        % If you prefer, the results can be plotted for each case
        % Either with a pause, or with a new figure for each case, and no pause.         
        %         figure(1)
        %         plot(phivec*180/pi,abs([tftot.' ]),'-o')
        %         xlabel(['Receiver angle   [deg.]'])
        %         ylabel(['Sound pressure amp., times distance   [-]'])
        %         grid
        %         title([filestem,GG,SS])
        %         pause 
        
        relerr = abs( (abs(tftot) - 1 ));
        medianerr_results(kk,jj) = median(relerr);

    end
end



