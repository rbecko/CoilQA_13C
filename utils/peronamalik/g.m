function gI = g(I,K)
% G, the anisotropic diffusion coefficient as a function of the input image
% Input         I = input image
%               K = contrast parameter
% Output        gI = anisotropic diffusion coefficient image

gI = 1./(1+(I/K).^2);