%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022    %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESSING
% SETTING UP NODAL COORDINATES FOR STREAMLINES
[X,Y] = meshgrid(0.5:nodx,0.5:nody); sx = 0.5*ones(1,21); sy = linspace(0,nody,21);
U = S(1:(nely+1)*(nelx+1)*2); umag=reshape(sqrt(U(1:2:end).^2+U(2:2:end).^2),nely+1,nelx+1);
% DESIGN FIELD
figure(1); imagesc(xPhys); colorbar; caxis([0 1]); axis equal; axis off;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
% BRINKMAN PENALTY FACTOR
figure(2); imagesc(reshape(log10(alpha),nely,nelx)); colorbar; caxis([0 log10(alphamax)]); axis equal; axis off; %colormap turbo;
% VELOCITY MAGNITUDE FIELD
figure(3); imagesc(umag); colorbar; axis equal; axis on; hold on; %colormap turbo;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
% PRESSURE FIELD
P = S(2*nodtot+1:3*nodtot);
figure(4); imagesc(reshape(P,nody,nodx)); colorbar; axis equal; axis off; %colormap turbo;
h = streamline(X,Y,reshape(U(1:2:end),nody,nodx),-reshape(U(2:2:end),nody,nodx),sx,sy); set(h,'Color','black');
% VELOCITY ALONG A LINE
if (probtype == 2)
    uline=flipud(diag(fliplr(umag))); xline=flipud(diag(fliplr(xPhys)));
else
    uline=umag(:,floor((end-1)/2)); xline=xPhys(:,floor(end/2));
end
figure(5);
subplot(3,1,1); plot(uline,'-x'); grid on; 
subplot(3,1,2); plot(log10(uline),'-x'); grid on;
subplot(3,1,3); plot(xline,'-x'); grid on; drawnow
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