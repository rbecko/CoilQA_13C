function pm_smoothed = perona_malik(img,K,iterations)
% PERONA_MALIK calculates anisotropic diffusion of the input image
% Input         img = input image
%               K = contrast parameter
%               iterations = number of iterations
% Output        pm_smoothed = anisotropic diffusion smoothed image

% Rie Beck Hansen, April 2014

% Delta t
dt = 1/5;

% Initialization
pm_smoothed = double(img);

for t = 1:iterations
    % Update the image
    pm_smoothed = pm_smoothed + dt*pm_rhs(pm_smoothed,K);
    
    % Uncomment if the user wants the iteration process to be printed
    % fprintf('Iteration %d\n',t);
end