%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% DEFINITION OF INPUT PARAMETERS
% PROBLEM TO SOLVE (1 = DOUBLE PIPE; 2 = PIPE BEND)
probtype = 1;
% DOMAIN SIZE
Lx = 1.0; Ly = 1.0;
% DOMAIN DISCRETISATION
nely = 30; nelx = nely*Lx/Ly;
% ALLOWABLE FLUID VOLUME FRACTION
volfrac = 1/3; xinit = volfrac;
% PHYSICAL PARAMETERS
Uin = 1e0; rho = 1e0; mu = 1e0;
% BRINKMAN PENALISATION
alphamax = 2.5*mu/(0.01^2); alphamin = 2.5*mu/(100^2);
% CONTINUATION STRATEGY
ainit = 2.5*mu/(0.1^2);
qinit = (-xinit*(alphamax-alphamin) - ainit + alphamax)/(xinit*(ainit-alphamin));
qavec = qinit./[1 2 10 20]; qanum = length(qavec); conit = 50;
% OPTIMISATION PARAMETERS
maxiter = qanum*conit; mvlim = 0.2; plotdes = 1;
chlim = 1e-3; chnum = 5;
% NEWTON SOLVER PARAMETERS
nltol = 1e-6; nlmax = 25; plotres = 0;
% EXPORT FILE
filename='output'; exportdxf = 0;
%% PREPARE FINITE ELEMENT ANALYSIS
dx = Lx/nelx; dy = Ly/nely;
nodx = nelx+1; nody = nely+1; nodtot = nodx*nody;
neltot = nelx*nely; doftot = 3*nodtot;
% NODAL CONNECTIVITY
nodenrs = reshape(1:nodtot,nody,nodx);
edofVecU = reshape(2*nodenrs(1:end-1,1:end-1)+1,neltot,1);
edofMatU = repmat(edofVecU,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],neltot,1);
edofVecP = reshape(nodenrs(1:end-1,1:end-1),neltot,1);
edofMatP = repmat(edofVecP,1,4)+repmat([1 nely+[2 1] 0],neltot,1);
edofMat = [edofMatU 2*nodtot+edofMatP];
iJ = reshape(kron(edofMat,ones(12,1))',144*neltot,1);
jJ = reshape(kron(edofMat,ones(1,12))',144*neltot,1);
iR = reshape(edofMat',12*neltot,1); jR = ones(12*neltot,1); 
jE = repmat(1:neltot,12,1);
%% DEFINE BOUNDARY CONDITIONS
% DEFINE THE PROBLEMS IN SEPARATE MATLAB FILE
run('problems.m');
% NULLSPACE MATRICES FOR IMPOSING BOUNDARY CONDITIONS
EN=speye(doftot); ND=EN; ND(fixedDofs,fixedDofs)=0.0; EN=EN-ND;
% VECTORS FOR FREE DOFS
alldofs = 1:doftot; freedofs = setdiff(alldofs,fixedDofs);
%% INITIALISATION
% SOLUTION VECTOR
S = zeros(doftot,1); dS = S; L = S; 
S(fixedDofs) = DIR(fixedDofs);
% DESIGN FIELD
xPhys = xinit*ones(nely,nelx); 
% COUNTERS
loop = 0; loopcont = 0; nlittot = 0; chcnt = 0;
% CHANGE
change = Inf; objOld = Inf;
% CONTINUATION
qastep = 1; qa = qavec(1);
% VECTORISED CONSTANTS
dxv = dx*ones(1,neltot); dyv = dy*ones(1,neltot);
muv = mu*ones(1,neltot); rhov = rho*ones(1,neltot);
%% OUTPUT PROBLEM INFORMATION
fprintf('=========================================================\n');
fprintf('      Problem number: %2i - Reynolds number: %3.2e\n',probtype,Renum);
fprintf('=========================================================\n');
fprintf('      Design it.:   0\n');
%% START ITERATION
destime = tic; ittime = tic;
while (loop <= maxiter)
    if (plotdes); figure(1); imagesc(xPhys); colorbar; caxis([0 1]); axis equal; axis off; drawnow; end
    %% GREYSCALE INDICATOR
    Md = 100*full(4*sum(xPhys(:).*(1-xPhys(:)))/neltot); 
    %% MATERIAL INTERPOLATION
    alpha = alphamin + (alphamax-alphamin)*(1-xPhys(:))./(1+qa*xPhys(:));
    dalpha = (qa*(alphamax - alphamin)*(xPhys(:) - 1))./(xPhys(:)*qa + 1).^2 - (alphamax - alphamin)./(xPhys(:)*qa + 1);
    %% NON-LINEAR NEWTON SOLVER
    normR = 1; nlit = 0; fail = -1; nltime = tic;
    while (fail ~= 1)
        nlit = nlit+1; nlittot = nlittot+1;
        % BUILD RESIDUAL AND JACOBIAN
        sR = RES(dxv,dyv,muv,rhov,alpha(:)',S(edofMat'));
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0;
        if (nlit == 1); r0 = norm(R); end
        r1 = norm(R); normR = r1/r0;
        if (plotres); figure(6); semilogy(nlittot,normR,'x'); axis square; grid on; hold on; end
        if (normR < nltol); break; end
        sJ = JAC(dxv,dyv,muv,rhov,alpha(:)',S(edofMat'));
        J = sparse(iJ,jJ,sJ(:)); J = (ND'*J*ND+EN);  
        % CALCULATE NEWTON STEP
        dS = -J\R;
        % L2-NORM LINE SEARCH
        Sp = S + 0.5*dS;
        sR = RES(dxv,dyv,muv,rhov,alpha(:)',Sp(edofMat'));
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0; r2 = norm(R);
        Sp = S + 1.0*dS;
        sR = RES(dxv,dyv,muv,rhov,alpha(:)',Sp(edofMat'));
        R = sparse(iR,jR,sR(:)); R(fixedDofs) = 0; r3 = norm(R);
        % SOLUTION UPDATE WITH "OPTIMAL" DAMPING
        lambda = max(0.01,min(1.0,(3*r1 + r3 - 4*r2)/(4*r1 + 4*r3 - 8*r2)));
        S = S + lambda*dS;
        % IF FAIL, RETRY FROM ZERO SOLUTION
        if (nlit == nlmax && fail < 0); nlit = 0; S(freedofs) = 0.0; normR=1; fail = fail+1; end
        if (nlit == nlmax && fail < 1); fail = fail+1; end
    end
    nltime=toc(nltime);
    fprintf('      Newton it.: %2i - Res. norm: %3.2e - Sol. time: %6.3f sec\n',nlit,normR,nltime);
    if (fail == 1); error('ERROR: Solver did not converge after retry from zero!\n      Stopping optimisation.\n'); end
    %% OBJECTIVE EVALUATION
    obj = sum( PHI(dxv,dyv,muv,alpha(:)',S(edofMat')) ); 
    change = abs(objOld-obj)/objOld; objOld = obj;
    %% VOLUME CONSTRAINT
    V = mean(xPhys(:));
    %% PRINT RESULTS
    ittime = toc(ittime);
    fprintf('      Obj.: %3.2e - Constr.: %3.2e - Md: %3.2f\n',obj,V,Md);
    fprintf('      Change: %4.3e - It. time: %6.3f sec\n',change,ittime);
    fprintf('      Contin. step: %2i - qa: %4.3e\n',qastep,qa);
    ittime = tic;
    %% EVALUATE CURRENT ITERATE - CONTINUE UNLESS CONSIDERED CONVERGED
    if (change < chlim); chcnt = chcnt + 1; else; chcnt = 0; end
    if (qastep == qanum && ( (chcnt == chnum) || (loopcont == conit) ) ); break; end
    %% PRINT HEADER FOR ITERATION
    loop = loop + 1; loopcont = loopcont + 1;
    fprintf('---------------------------------------------------------\n');
    fprintf('      Design it.:%4i\n',loop);
    %% ADJOINT SOLVER
    sR = [dPHIds(dxv,dyv,muv,alpha(:)',S(edofMat')); zeros(4,neltot)];
    RHS = sparse(iR,jR,sR(:)); RHS(fixedDofs) = 0;
    sJ = JAC(dxv,dyv,muv,rhov,alpha(:)',S(edofMat'));
    J = sparse(iJ,jJ,sJ(:)); J = (ND'*J*ND+EN);
    L = J'\RHS;
    %% COMPUTE SENSITIVITIES
    % OBJECTIVE
    sR = dRESdg(dxv,dyv,muv,rhov,alpha(:)',dalpha(:)',S(edofMat'));
    dRdg = sparse(iR(:),jE(:),sR(:));
    dphidg = dPHIdg(dxv,dyv,muv,alpha(:)',dalpha(:)',S(edofMat'));
    sens = reshape(dphidg - L'*dRdg,nely,nelx);
    % VOLUME CONSTRAINT
    dV = ones(nely,nelx)/neltot;
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    xnew = xPhys; xlow = xPhys(:)-mvlim; xupp = xPhys(:)+mvlim;
    ocfac = xPhys(:).*max(1e-10,(-sens(:)./dV(:))).^(1/3);
    l1 = 0; l2 = ( 1/(neltot*volfrac)*sum( ocfac ) )^3;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew(:) = max(0,max(xlow,min(1,min(xupp,ocfac/(lmid^(1/3))))));
        if mean(xnew(:)) > volfrac; l1 = lmid; else; l2 = lmid; end
    end
    xPhys = xnew;
    %% CONTINUATION UPDATE
    if (qastep < qanum && (loopcont == conit || chcnt == chnum) )
        loopcont = 0; chcnt = 0;
        qastep = qastep + 1; qa = qavec(qastep);
    end
end
%% PRINT FINAL INFORMATION
destime = toc(destime);
fprintf('=========================================================\n');
fprintf('      Number of design iterations: %4i\n',loop);
fprintf('      Final objective: %4.3e\n',obj);
fprintf('      Total time taken: %6.2f min\n',destime/60);
fprintf('=========================================================\n');
%% PLOT RESULTS
run('postproc.m');
if (exportdxf); run('export.m'); end
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
% The basic structure of the code is based on the 88-line code for        %
% elastic compliance from: "Efficient topology optimization in MATLAB     %
% using 88 lines of code", E. Andreassen, A. Clausen, M. Schevenels,      %
% B. S. Lazarov and O. Sigmund, SMO 2010, doi:10.1007/s00158-010-0594-7   %
%                                                                         %
% Disclaimer:                                                             %
% The author does not guarantee that the code is free from errors.        %
% Furthermore, the author shall not be liable in any event caused by the  %
% use of the program.                                                     %      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

