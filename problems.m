%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE THE PROBLEMS TO BE SOLVED  %%%%
if (probtype == 1) % DOUBLE PIPE PROBLEM
    if ( mod(nely,6) > 0 > 0 )
        error('ERROR: Number of elements in y-dir. must be divisable by 6.');
    elseif ( nely < 30 )
        error('ERROR: Number of elements in y-dir. must be above 30.');
    end
    inletLength = 1/6*nely; inlet1 = 1/6*nely+1; inlet2 = 4/6*nely+1;
    outletLength = 1/6*nely; outlet1 = 1/6*nely+1; outlet2 = 4/6*nely+1; 
    nodesInlet  = [inlet1:(inlet1+inletLength)' inlet2:(inlet2+inletLength)'];
    nodesOutlet = (nodx-1)*nody+[outlet1:(outlet1+outletLength)' outlet2:(outlet2+outletLength)'];
    nodesTopBot = [1:nely+1:nodtot nody:nody:nodtot];
    nodesLefRig = [2:nely (nelx)*nody+2:nodtot-1];
    fixedDofsTBx = 2*nodesTopBot-1; fixedDofsTBy = 2*nodesTopBot;
    fixedDofsLRx = 2*setdiff(nodesLefRig,[nodesInlet nodesOutlet])-1;
    fixedDofsLRy = 2*nodesLefRig;
    fixedDofsInX  = 2*nodesInlet-1; fixedDofsInY  = 2*nodesInlet;
    fixedDofsOutY = 2*nodesOutlet; fixedDofsOutP = 2*nodtot+nodesOutlet;
    fixedDofsU = [fixedDofsTBx fixedDofsTBy fixedDofsLRx fixedDofsLRy ...
                  fixedDofsInX fixedDofsInY fixedDofsOutY];
    fixedDofsP = [fixedDofsOutP];
    fixedDofs  = [fixedDofsU fixedDofsP];
    % DIRICHLET VECTORS
    u = @(y) -4*y.^2+4*y; Uinlet = Uin*u([0:inletLength]'/inletLength);
    DIRU=zeros(nodtot*2,1); DIRU(fixedDofsInX) = [Uinlet' Uinlet'];
    DIRP=zeros(nodtot,1); DIR = [DIRU; DIRP];     
    % INLET REYNOLDS NUMBER
    Renum = Uin*(inletLength*Ly/nely)*rho/mu;
elseif (probtype == 2) % PIPE BEND PROBLEM
    if ( mod(nelx,10) > 0 || mod(nely,10) > 0 )
        error('ERROR: Number of elements per side must be divisable by 10.');
    end
    inletLength = 2/10*nely; inlet1 = 1/10*nely+1; 
    outletLength = 2/10*nelx; outlet1 = 7/10*nelx+1; 
    nodesInlet  = [inlet1:(inlet1+inletLength)];
    nodesOutlet = nody*[outlet1:(outlet1+outletLength)];
    nodesTopBot = [1:nely+1:nodtot nody:nody:nodtot];
    nodesLefRig = [2:nely (nodx-1)*nody+2:nodtot-1];
    fixedDofsTBx = 2*nodesTopBot-1; fixedDofsTBy = 2*setdiff(nodesTopBot,nodesOutlet);
    fixedDofsLRx = 2*setdiff(nodesLefRig,nodesInlet)-1; fixedDofsLRy = 2*nodesLefRig;
    fixedDofsInX  = 2*nodesInlet-1; fixedDofsInY  = 2*nodesInlet;
    fixedDofsOutX = 2*nodesOutlet-1; fixedDofsOutP = 2*nodtot+nodesOutlet;
    fixedDofsU = [fixedDofsTBx fixedDofsTBy fixedDofsLRx fixedDofsLRy ...
                  fixedDofsInX fixedDofsInY fixedDofsOutX];
    fixedDofsP = [fixedDofsOutP];
    fixedDofs  = [fixedDofsU fixedDofsP];
    % DIRICHLET VECTORS
    DIRU=zeros(nodtot*2,1); DIRP=zeros(nodtot,1);
    u = @(y) -4*y.^2+4*y; Uinlet = Uin*u([0:inletLength]'/inletLength);
    DIRU(fixedDofsInX) = Uinlet'; DIR = [DIRU; DIRP];     
    % INLET REYNOLDS NUMBER
    Renum = Uin*(inletLength*Ly/nely)*rho/mu;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by: Joe Alexandersen                              %
%                           Department of Mechanical and                  %
%                                         Electrical Engineering          %
%                           University of Southern Denmark                %
%                           DK-5230 Odense M, Denmark.                    %
% Please send your comments and questions to: joal@sdu.dk                 %
%                                                                         %
% The code is intended for educational purposes and theoretical details   %
% are discussed in the paper: "A detailed introduction to density-based   %
% topology optimisation of fluid flow problems including implementation   %
% in MATLAB", J. Alexandersen, SMO 2022, doi:                             %                          
%                                                                         %
% A preprint version of the paper can be downloaded from the author's     %
% website: joealexandersen.com                                            %
% The code is available from GitHub: github.com/sdu-multiphysics/topflow  %
%                                                                         %
% Disclaimer:                                                             %
% The author does not guarantee that the code is free from errors.        %
% Furthermore, the author shall not be liable in any event caused by the  %
% use of the program.                                                     %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
