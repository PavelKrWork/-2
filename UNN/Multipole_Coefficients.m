function [ AE_lm, AM_lm] = Multipole_Coefficients( f, SingleAntennaParameters)
% In this function, we calculate the Multipole Coefficients for the given frequency f

global c Z0

c = physconst('LightSpeed');

k = 2 * pi * f / c;
lm = SingleAntennaParameters.lmax;

if size(lm, 2) > 1
    l = lm(1);
    m = lm(2);
    
    [aE, aM] = Coefficients(l, m);
    
    AE_lm = ( k / sqrt(l*(l + 1)) ) .* aE;
    AM_lm = ( -1i * k.^2 / sqrt(l*(l + 1)) ) .* aM;
    
else
    lmax = lm;
    AE_lm = zeros(lmax, 2*lmax + 1);
    AM_lm = zeros(lmax, 2*lmax + 1);
    
    for l = 1 : lmax
        for m = -l : l

            [aE, aM] = Coefficients(l, m);
            
            AE_lm(l, l + m + 1) = ( k / sqrt(l*(l + 1)) ) .* aE;
            AM_lm(l, l + m + 1) = ( -1i * k.^2 / sqrt(l*(l + 1)) ) .* aM;
            
        end
    end
    
end


%% Local function for calculating the (l,m)-order Multipole coefficient 
    function [AE_lm, AM_lm] = Coefficients(l, m)
        
        x = SingleAntennaParameters.Current_density(:,1) * 1E-3;
        y = SingleAntennaParameters.Current_density(:,2) * 1E-3;
        z = SingleAntennaParameters.Current_density(:,3) * 1E-3;
        
        Jx = SingleAntennaParameters.Current_density(:,4) + ...
            1i * SingleAntennaParameters.Current_density(:,7);
        Jy = SingleAntennaParameters.Current_density(:,5) + ...
            1i * SingleAntennaParameters.Current_density(:,8);
        Jz = SingleAntennaParameters.Current_density(:,6) + ...
            1i * SingleAntennaParameters.Current_density(:,9);
        SurfaceElementArea = SingleAntennaParameters.Current_density(:,10) * 1E-6;
        
        % we calculate the coordinates of the integration points in a spherical coordinate system
        r = sqrt(x.^2 + y.^2 + z.^2);
        rho = sqrt(x.^2 + y.^2);
        theta = acos( z ./ ( 0.01 .* (r <= 0) + (r) .* (r > 0) ) ) .* (r > 0);
        phi = ( acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) .* ( y >= 0 ) + ...
              (2*pi - acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) ) .* ( y < 0 ) ) .* (rho > 0);

        kr = k .* r;  
        % jl = sqrt(pi ./ (2 * kr)) .* besselj(l + 1/2, kr);
        jl_minus1 = sqrt(pi ./ (2 * kr)) .* besselj(l - 1 + 1/2, kr);
        jl_plus1 = sqrt(pi ./ (2 * kr)) .* besselj(l + 1 + 1/2, kr);

        % we calculate the Associated Legendre polinomials
        LegendrePolinom_lm = legendre(l, cos(theta));
        
        m_minus1 = (abs(m) - 1) .* (abs(m) > 0);
        m_plus1 = (abs(m) + 1);
        
        % we use the recurence relation to calculate the derivative of the
        % Associated Legendre polinomials with respect to the variable \theta
        Plm = reshape(LegendrePolinom_lm( abs(m) + 1,:,:,:), size(theta));
        Plm_minus1 = reshape(LegendrePolinom_lm( (m_minus1 + 1),:,:,:), ...
            size(theta)) .* (abs(m) > 0);
        Plm_plus1 = reshape(LegendrePolinom_lm( (m_plus1 .* (m_plus1 <= l) + 1),:,:,:), ...
            size(theta)) .* ( 1 .* (m_plus1 <= l) + 1 .* (m == 0));
        Plm_over_theta = 0.5 * ( Plm_plus1  - (l + abs(m)) * (l - abs(m) + 1) * Plm_minus1);
                

        Alm = sqrt( ( (2*l + 1)/(4*pi) ) * ( factorial(l - abs(m)) / factorial(l + abs(m)) ) );
        Alm = Alm * ( (-1)^abs(m) .* (m < 0) + 1 .* (m >= 0) );

        % % we calculate the components of the vectors LE and LM in a spherical coordinate system
        LE_r = Alm * k * ( (l^2 + l) / (2*l + 1) ) .* ( jl_minus1 + jl_plus1) .* ...
            Plm .* exp(1i * m * phi);
        LE_theta = Alm * k * ( 1 / (2*l + 1) ) .* ( (l + 1) * jl_minus1 - l * jl_plus1) .* ...
            Plm_over_theta .* exp(1i * m * phi);
        LE_phi = Alm * (1i * k * m) * ( 1 / (2*l + 1) ) .* ( (l + 1) * jl_minus1 - l * jl_plus1) .* ...
            (Plm ./ sin(theta)) .* exp(1i * m * phi);

        LM_r = Alm * k * ( 1 / (2*l + 1) ) .* ( l * jl_minus1 - (l + 1) * jl_plus1) .* ...
            Plm .* exp(1i * m * phi);
        LM_theta = Alm * k * ( 1 / (2*l + 1) ) .* ( jl_minus1 + jl_plus1) .* ...
            Plm_over_theta .* exp(1i * m * phi);
        LM_phi = Alm * (1i * k * m) * ( 1 / (2*l + 1) ) .* ( jl_minus1 + jl_plus1) .* ...
            (Plm ./ sin(theta)) .* exp(1i * m * phi);
            
       % % we calculate the Cartesian components of the vectors LE and LM
        LE_x = LE_r .* sin(theta) .* cos(phi) + LE_theta .* cos(theta) .* cos(phi) - LE_phi .* sin(phi);
        LE_y = LE_r .* sin(theta) .* sin(phi) + LE_theta .* cos(theta) .* sin(phi) + LE_phi .* cos(phi);
        LE_z = LE_r .* cos(theta) - LE_theta .* sin(theta);
        LM_x = LM_r .* sin(theta) .* cos(phi) + LM_theta .* cos(theta) .* cos(phi) - LM_phi .* sin(phi);
        LM_y = LM_r .* sin(theta) .* sin(phi) + LM_theta .* cos(theta) .* sin(phi) + LM_phi .* cos(phi);
        LM_z = LM_r .* cos(theta) - LM_theta .* sin(theta);
        
        % % we calculate the multipole coefficients AE_lm and AM_lm
        AE_lm = sum( ( Jx .* LE_x + Jy .* LE_y + Jz .* LE_z ) .* SurfaceElementArea );
        AM_lm = sum( ( (y.*Jz - z.*Jy) .* LM_x + ...
                (z.*Jx - x.*Jz) .* LM_y + (x.*Jy - y.*Jx) .* LM_z ) .* SurfaceElementArea );
   
    end

 
end

