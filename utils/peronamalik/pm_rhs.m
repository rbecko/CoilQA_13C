function rhs = pm_rhs(I,K)
% PM_RHS calculates the right-hand side of the Perona-Malik equation
% Input         I = input image
%               K = contrast parameter
% Output        rhs = right-hand side of the Perona-Malik equation

% Rie Beck Hansen, April 2014

% Step size
dx = 1; 
dy = 1;

% The directional derivatives, 
% e.g. nabla_North_I(x,y,t) = I(x,y+?y,t) - I(x,y,t)
nablaN = [-I(1,:);-diff(I)];
nablaS = [-diff(I);I(end,:)];
nablaE = [diff(I,1,2) -I(:,end)];
nablaW = [I(:,1) diff(I,1,2)];

% Calculate diffusion coefficient for different directions
gIN = g(abs(nablaN),K);
gIS = g(abs(nablaS),K);
gIE = g(abs(nablaE),K);
gIW = g(abs(nablaW),K);

% Calculate the right-hand side of the Perona-Malik equation
rhs = 1/dx^2*(gIE.*nablaE)-1/dx^2*(gIW.*nablaW) ...
    + 1/dy^2*(gIN.*nablaN)-1/dy^2*(gIS.*nablaS);
