function [Kphi, Fphi] = assembleElasKPhiRes(GaussInfo, elem, Disp, Phi, Para)
% -------------------------------------------------------------------
% % ** code by P.M.H. @BIT (CN) **
% Calculate Phi stiffness matrices
% %  Please feel free to contact us with any questions! 
% %  - Email: pm_hu@outlook.com
% %  ---------------------------------------

Gc = Para.Gc;    % energy release rate
Len = Para.Len;  % length scale
lambda = Para.lambda;
mu = Para.mu; % Lame Constant
PFModel = Para.PFModel;

switch PFModel  % % =======AT1/2 model parameters========
    case 1 % AT2 --Bourdin 
        dalpha = @(x) 2*x;
        ddalpha = 2;
        Calpha = 2;
        a2 = 1;
    case 2 % AT1 -- Pham
        dalpha = @(x) 1;
        ddalpha = 0;
        Calpha = 8/3;
        a2 = 1;
    case 3 % PFCZM -- Wu - Linear soften law
        error('PFCZM模型未实现');
    otherwise
        error('不支持的相场模型');
end

numEleNd  = size(elem, 2);  % 单元结点数
numEle = size(elem, 1); % 单元数
numEDofs = numEleNd * Para.ndim;

KphiVals = zeros((numEleNd)^2, numEle); % store the stiff matrix
FphiVals = zeros(numEleNd, numEle); % store the rhs vector

for ei = 1 : numEle
    elei = elem(ei,:);
    eleDOFs = reshape([2*elei-1; 2*elei], numEDofs,1);
    
    Ke = zeros(numEleNd); % element stiff-phi
    Fe = zeros(numEleNd,1); % element rhs-phi
    
    % loading FEM information
    dRdxGaussPt = GaussInfo.SpDeriv{ei};
    RGaussPt = GaussInfo.SpVal{ei};
    JW = GaussInfo.JW{ei};
    HisyGaussPt = GaussInfo.Hisy{ei};
    
    for gpti = 1 : size(dRdxGaussPt,3)
        
        %Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dRdxGaussPt( :, :, gpti);
        R = RGaussPt(gpti, :);
        
        GPphi = R * Phi(elei);       % phi at GPt
        
        Bd = zeros(2, numEleNd);
        Bd(1, :) = dRdx(1, :);
        Bd(2, :) = dRdx(2, :);  % hat{B}
        
        dGPphi = Bd * Phi(elei);   % phi-first-dirv at GPt        
        
        % %      _                                             _
        % %     | N_{1, x} 0         ... N_{m, x} 0         ...|
        % %  B =| 0        N_{1, y}  ... 0        N_{m, y}  ...|
        % %     | N_{1, y} N_{1, x}  ... N_{m, y} N_{m, x}  ...|
        % %      -                                             -
        % % \sigma = [\sigma_{xx} \sigma_{yy} \sigma_{xy}]
        B = zeros(3, 2 * numEleNd);
        B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
        B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
        
        B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);
        
        % 使用已更新的历史变量（在updateRefEnerg中已更新）
        GPHisy = HisyGaussPt(gpti);  % 直接使用已更新的历史变量
        
        switch PFModel
            case 1 % AT2
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2;          % gb = (1-d)^2
                dgdd = -2*(1-GPphi); % gb = (1-d)^2
            case 2 % AT1
                % derivatives of degradation function (g - degradation function)
                d2gdd2 = 2;          % gb = (1-d)^2
                dgdd = -2*(1-GPphi); % gb = (1-d)^2
                GPHisy = max(GPHisy, 3*Gc/(16*Len)); 
        end
        
        % compute element stiffness at quadrature point
        Ke = Ke + Gc*2*a2*Len/Calpha * (Bd' * Bd) * JW(gpti);  %
        Ke = Ke + (Gc*ddalpha/(Calpha*Len) + d2gdd2*GPHisy) * (R' * R) * JW(gpti); % %
        
        Fe = Fe + (dgdd*GPHisy) * R' * JW(gpti); % % 
        Fe = Fe + Gc/(Calpha*Len) *( R'*dalpha(GPphi) + 2*a2*Len^2*(Bd'*dGPphi) )* JW(gpti);
    end
    KphiVals(:, ei) = Ke(:);
    FphiVals(:, ei) = Fe(:);
end

% %
J = repmat(1 : numEleNd, numEleNd, 1);
I = J';
ElConn = elem;

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

Kphi = sparse(ii(:), jj(:), KphiVals(:)); % assemble Kphi
Kphi = (Kphi + Kphi')/2;

kk = reshape(ElConn',[],1);
Fphi = sparse(kk, 1, FphiVals(:));  % assemble rhs

end

