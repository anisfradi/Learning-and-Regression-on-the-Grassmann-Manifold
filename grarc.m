function sqD = grarc(Y1,Y2)
%GRARC Arc-length on Grassmann manifold (see Edelman et al., Sec. 4.3)
%
% Author: Roland Kwitt, 2013
% E-Mail: roland.kwitt@kitware.com

    [~,S,~] = svd(Y1'*Y2);
    theta = acos(diag(S)); % diagonal elements are cos(theta_i)
    %sqD = sum(theta.^2);
    sqD = norm(theta, 2).^2;
end
