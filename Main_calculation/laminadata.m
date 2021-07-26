function [Q] = laminadata_2( E_L, E_T, nu_LT, G_LT)
%LAMINADATA CREATES THE NESSESSARY MATERIAL DATA FOR A LAMINA FROM FIBRE
%AND MATRIX MATERIAL DATA
%   Input:
%     E_f           -       Young's modulus for fibre
%     E_m           -       Young's modulus for matrix
%     nu_f          -       Poisson's ratio for fibre
%     nu_m          -       Poisson's ratio for matrix
%     V_f           -       Volumefraction for fibre
%     V_m           -       Volumefraction for matrix
%     xsi_E         -       xsi term in Halpin-Tsai for Young's modulus
%     xsi_G         -       xsi term in Halpin-Tsai for Shear modulus
% 
%   Output:
%     Q             -       3x3 Q matrix of a laminate
% 

% ===== RESULTS: =====


% === nu ===

nu_TL   = (nu_LT*E_T)/E_L; %Minor

% === Q ===
Q       = [E_L/(1-nu_LT*nu_TL) nu_TL*E_L/(1-nu_LT*nu_TL) 0;...
           nu_LT*E_T/(1-nu_LT*nu_TL) E_T/(1-nu_LT*nu_TL) 0;...
           0 0 G_LT]; 

end

