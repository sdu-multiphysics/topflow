%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% Switches
buildJR = 1; exportJR = 1;
buildPhi = 1; exportPhi = 1;
%% Initialisation
syms rho mu alpha dalpha xi eta dx dy;
syms u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4;
time = tic;
%% Shape functions and matrices
xv = [-1 1 1 -1]'; yv = [-1 -1 1 1]';
Np = 1/4*(1+xi*xv').*(1+eta*yv');
Nu = sym(zeros(2,8));
Nu(1,1:2:end-1) = Np;
Nu(2,2:2:end) = Np;
%% Nodal coordinates, interpolation and coordinate transforms
xc = dx/2*xv; yc = dy/2*yv;
x = Np*xc; y = Np*yc;
J = [diff(x,xi) diff(y,xi); diff(x,eta) diff(y,eta)];
iJ = inv(J); detJ = det(J);
%% Derivatives of shape functions
dNpdx = iJ(:,1)*diff(Np,xi) + iJ(:,2)*diff(Np,eta);
dNudx = sym(zeros(2,8,2));
for i = 1:2
   dNudx(i,:,:) = transpose(iJ(:,1)*diff(Nu(i,:),xi) + iJ(:,2)*diff(Nu(i,:),eta));
end
%% Nodal DOFs
u = [u1; u2; u3; u4; u5; u6; u7; u8];
p = [p1; p2; p3; p4]; s = [u; p];
ux = Nu*u; px = Np*p;
dudx = [dNudx(:,:,1)*u dNudx(:,:,2)*u];
dpdx = dNpdx*p;
%% Stabilisation parameters
h = sqrt(dx^2 + dy^2);
u0 = subs(ux,[xi,eta],[0,0]);
ue = sqrt(transpose(u0)*u0);
tau1 = h/(2*ue);
tau3 = rho*h^2/(12*mu);
tau4 = rho/alpha;
tau  = (tau1^(-2) + tau3^(-2) + tau4^(-2))^(-1/2);
%% Loop over the tensor weak form to form residual
if (buildJR)
    Ru = sym(zeros(8,1)); Rp = sym(zeros(4,1));
    % Momentum equations
    for g = 1:8
        for i = 1:2
            Ru(g) = Ru(g) + alpha*Nu(i,g)*ux(i); % Brinkman term
            for j = 1:2
                Ru(g) = Ru(g) + mu*dNudx(i,g,j)*( dudx(i,j) + dudx(j,i) ); % Viscous term
                Ru(g) = Ru(g) + rho*Nu(i,g)*ux(j)*dudx(i,j); % Convection term
                Ru(g) = Ru(g) + tau*ux(j)*dNudx(i,g,j)*alpha*ux(i); % SUPG Brinkman term
                for k = 1:2
                    Ru(g) = Ru(g) + tau*ux(j)*dNudx(i,g,j)*rho*ux(k)*dudx(i,k); % SUPG convection term
                end
                Ru(g) = Ru(g) + tau*ux(j)*dNudx(i,g,j)*dpdx(i); % SUPG pressure term
            end
            Ru(g) = Ru(g) - dNudx(i,g,i)*px; % Pressure term
        end
    end
    % Incompressibility equations
    for g = 1:4
        for i = 1:2
            Rp(g) = Rp(g) + Np(1,g)*dudx(i,i); % Divergence term
            Rp(g) = Rp(g) + tau/rho*dNpdx(i,g)*alpha*ux(i); % PSPG Brinkman term
            for j = 1:2
                Rp(g) = Rp(g) + tau*dNpdx(i,g)*ux(j)*dudx(i,j); % PSPG convection term
            end
            Rp(g) = Rp(g) + tau/rho*dNpdx(i,g)*dpdx(i); % PSPG pressure term
        end
    end
    fprintf('Simplifying Ru... \n');
    Ru = simplify(detJ*Ru);
    fprintf('Simplifying Rp... \n');
    Rp = simplify(detJ*Rp);
    %% Integrate analytically
    fprintf('Integrating Ru... \n');
    Ru = int(int(Ru,xi,[-1 1]),eta,[-1 1]);
    fprintf('Integrating Rp... \n');
    Rp = int(int(Rp,xi,[-1 1]),eta,[-1 1]);
    fprintf('Simplifying Ru... \n');
    Ru = simplify(Ru);
    fprintf('Simplifying Rp... \n');
    Rp = simplify(Rp);
    %% Differentiate residual to form Jacobian
    Re = [Ru; Rp];
    Je = sym(zeros(12,12));
    for b = 1:12
        fprintf('Computing dR/ds%2i... \n',b);
        Je(:,b) = diff(Re,s(b));
    end
end
%% Export residual and Jacobian
if (exportJR)
    fprintf('Exporting residual... \n');
    f = matlabFunction(Re,'File','RES','Comments','Version: 0.99','Vars',{dx,dy,mu,rho,alpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
    fprintf('Exporting Jacobian... \n');
    f = matlabFunction(Je(:),'File','JAC','Comments','Version: 0.99','Vars',{dx,dy,mu,rho,alpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
end
%%
time = toc(time);
fprintf('Finished analysis part in %1.3e seconds. \n',time)
%% OPTIMISATION PART
time = tic;
if (buildPhi)
    %% Compute objective functional
    fprintf('Computing phi... \n');
    phi = 1/2*alpha*transpose(ux)*ux;
    for i = 1:2
        for j = 1:2
            phi = phi + 1/2*mu*dudx(i,j)*( dudx(i,j) + dudx(j,i) );
        end
    end
    %% Intregrate and simplify objective functional
    fprintf('Integrating phi... \n');
    phi = int(int(detJ*phi,xi,[-1 1]),eta,[-1 1]);
    phi = simplify(phi);
    %% Compute the partial derivative wrt. design field
    fprintf('Computing dphi/dgamma... \n');
    dphidg = simplify( diff(phi,alpha)*dalpha );
    %% Compute the partial derivative wrt. state field
    dphids = sym(zeros(12,1));
    for a = 1:12
        fprintf('Computing dphi/ds%2i... \n',a);
        dphids(a) = simplify(diff(phi,s(a)));
    end
end
%% Compute partial derivative of residual wrt. design field
if (buildJR)
    fprintf('Computing dr/dgamma... \n');
    drdg = simplify( diff(Re,alpha)*dalpha );
end
%% Export optimisation functions
if (exportPhi)
    fprintf('Exporting phi... \n');
    f = matlabFunction(phi,'File','PHI','Comments','Version: 0.9','Vars',{dx,dy,mu,alpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
    fprintf('Exporting dphi/dg... \n');
    f = matlabFunction(dphidg,'File','dPHIdg','Comments','Version: 0.9','Vars',{dx,dy,mu,alpha,dalpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
    fprintf('Exporting dphi/ds... \n');
    f = matlabFunction(dphids(1:8),'File','dPHIds','Comments','Version: 0.9','Vars',{dx,dy,mu,alpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
end
if (exportJR)
    fprintf('Exporting dr/ds... \n');
    f = matlabFunction(drdg,'File','dRESdg','Comments','Version: 0.9','Vars',{dx,dy,mu,rho,alpha,dalpha,transpose([u1 u2 u3 u4 u5 u6 u7 u8 p1 p2 p3 p4])});
end
%%
time = toc(time);
fprintf("Finished optimisation part in %1.3e seconds. \n",time)
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




