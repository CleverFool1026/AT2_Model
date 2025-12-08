function  Phi = NewtonItPhaseField(GaussInfo, elem, Disp, Phi_1, Para)
% 
% *** Newton Raphson (NR) approximation for nonlinear problem ***
% ------------------ Now take a brief review --------------------------
% Given a differentiable function F(x), the root x* (that is F(x*) = 0)
% is approximated by the incremental {\Delta x} and initial guess {x0}: 
% x* = x0 + \Delta x,  where \Delta x = - F(x0) / F'(x0).
% ---------------------------------------------------------------------
% Create by P.M.Hu @ BIT(CN) 2022-11-20
% Please feel free to contact us with any questions! 
% - Email: pm_hu@outlook.com
% %  ---------------------------------------

NRresidual = 1; 
tol = 1d-4;
maxIter = 20;
Phi = Phi_1;
iter = 0;
while NRresidual > tol && iter < maxIter
    
    [Kphi, Rphi] = assembleElasKPhiRes(GaussInfo, elem, Disp, Phi, Para);
    
    f =  -Rphi;
    % Solve the system by mldivide - direct method
    dphi = Kphi \ f; % 
 
    Phi = Phi + dphi;
    NRresidual = norm(dphi) / (norm(Phi) + 1e-10);
    iter = iter + 1;
end

Phi(Phi < 0) = 0; % enforce boundary constraint [0 1]
Phi(Phi >= 1) = 1; % enforce boundary constraint [0 1]
end

