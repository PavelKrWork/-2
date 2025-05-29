function [ FieldComponent_x, FieldComponent_y, FieldComponent_z] = ...
    SingleAntennaFieldEH( x, y, z, f, SingleAntennaParameters, EMfield_Components )
% In this function, we calculate the EM field of a single antenna
% The calculation technique, which is used to compute the EM field, 
% is specified by the string variable SingleAntennaParameters.CalculationTechnique


global c Z0

Z0 = SingleAntennaParameters.Z0;

c = physconst('LightSpeed');

k = 2 * pi * f / c;
 
 switch SingleAntennaParameters.CalculationTechnique
     case 'PointDipole' % % EM field is calculated using the formulas for the fields of the infinitesimal dipole
         % 'PointDipole' = infinitesimal dipole
                 
         L_infdip = SingleAntennaParameters.length; 
         I_infdip = -SingleAntennaParameters.I0;
                 
         % % coordinates in formulas for the infinitesimal dipole model
         r = sqrt(x.^2 + y.^2 + z.^2);
         rho = sqrt(x.^2 + y.^2);
         theta = acos( z ./ ( 0.01 + (r - 0.01) .* (r > 0) ) ) .* (r > 0);
         phi = ( acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) .* ( y >= 0 ) + ...
             (2*pi - acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ))  .* ( y < 0 ) ) .* (rho > 0);

         Er_infdip = Z0 * I_infdip * L_infdip * cos(theta) ./ (2 * pi * r.^2) .* (1 + 1./ (1i * k * r)) .* exp(-1i*k*r);
         Etheta_infdip = 1i * Z0 * k * I_infdip * L_infdip * sin(theta) ./ (4 * pi * r) .* ...
             (1 + 1./ (1i * k * r) - 1./ (k * r).^2) .* exp(-1i*k*r);
         
         Hphi_infdip = 1i * k * I_infdip * L_infdip * sin(theta) ./ (4 * pi * r) .* ...
             (1 + 1./ (1i * k * r) ) .* exp(-1i*k*r);
        
         Ex = Er_infdip .* sin(theta) .* cos(phi) + Etheta_infdip .* cos(theta) .* cos(phi);
         Ey = Er_infdip .* sin(theta) .* sin(phi) + Etheta_infdip .* cos(theta) .* sin(phi);
         Ez = Er_infdip .* cos(theta) - Etheta_infdip .* sin(theta);
         Hx = - Hphi_infdip .* sin(phi);
         Hy = Hphi_infdip .* cos(phi);
         Hz = 0;
         
         switch EMfield_Components
             case 'E_field'
                 FieldComponent_x = Ex;
                 FieldComponent_y = Ey;
                 FieldComponent_z = Ez;
             case 'H_field'
                 FieldComponent_x = Hx;
                 FieldComponent_y = Hy;
                 FieldComponent_z = Hz;
         end
      
     
     case 'LinearDipole' % % EM field is calculated using the formulas for the fields of Linear dipole with a sine current distribution
        
         L = SingleAntennaParameters.length; 
         I_0 = SingleAntennaParameters.I0;
         I_norm = -1 * I_0 / sin(k*L/2);
         
         
         r = sqrt(x.^2 + y.^2 + z.^2);
         rho = sqrt(x.^2 + y.^2);
         % theta = acos( z ./ ( 0.01 + (r - 0.01) .* (r > 0) ) ) .* (r > 0);
         phi = ( acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) .* ( y >= 0 ) + ...
             (2*pi - acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ))  .* ( y < 0 ) ) .* (rho > 0);
         
         R1 = sqrt(rho.^2 + (z - L/2).^2);
         R2 = sqrt(rho.^2 + (z + L/2).^2);
         % r = r;
         
         exp0 = exp(-1i * k .* r);
         exp1 = exp(-1i * k .* R1);
         exp2 = exp(-1i * k .* R2);
         
         Ez = I_norm * ( -1i * Z0 ./ ( 4 * pi)) .* (exp1 ./ R1 + exp2 ./ R2 - ...
             2 * cos(k * L/2) .* exp0 ./ r);
         E_rho = I_norm * ( 1i * Z0 ./ (4 * pi * rho) ) .* (exp1 .* (z - L/2) ./ R1 + exp2 .* (z + L/2) ./ R2 - ...
             2 * cos(k * L/2) .* exp0 .* z ./ r);
         
         H_phi = I_norm * (1i ./ (4 * pi * rho)) .* (exp1 + exp2 - 2 * cos(k * L/2) .* exp0);
         Hz = 0;
         
         Ex = E_rho .* cos(phi);
         Ey = E_rho .* sin(phi);
         Hx = -H_phi .* sin(phi);
         Hy = H_phi .* cos(phi);
         
         switch EMfield_Components
             case 'E_field'
                 FieldComponent_x = Ex;
                 FieldComponent_y = Ey;
                 FieldComponent_z = Ez;
             case 'H_field'
                 FieldComponent_x = Hx;
                 FieldComponent_y = Hy;
                 FieldComponent_z = Hz;
         end
             
     
     case 'Green' % EM field is calculated according to the Green's Function Approach
         
         [FieldComponent_x, FieldComponent_y, FieldComponent_z] = GreenFunctionApproach(x, y, z, EMfield_Components);
  
     case 'Multipole' % EM field is calculated according to the Multipole Expansion technique
                  
         [FieldComponent_x, FieldComponent_y, FieldComponent_z] = MultipoleExpansion(x, y, z, EMfield_Components);
  
     case 'Combined_Green+Multipole' % Combined approach based on the Green's function and Multipole expansion technique
         
         R0 = SingleAntennaParameters.R0;
         MultExpValidity = 1 .* (sqrt(x.^2 + y.^2 + z.^2) >= R0);
         
         [GreenFun_x, GreenFun_y, GreenFun_z] = GreenFunctionApproach(x, y, z, EMfield_Components);
         global c Z0;
         [MultEpp_x, MultEpp_y, MultEpp_z] = MultipoleExpansion(x, y, z, EMfield_Components);
                  
         FieldComponent_x = GreenFun_x .* (1 - MultExpValidity) + MultEpp_x .* MultExpValidity;
         FieldComponent_y = GreenFun_y .* (1 - MultExpValidity) + MultEpp_y .* MultExpValidity;
         FieldComponent_z = GreenFun_z .* (1 - MultExpValidity) + MultEpp_z .* MultExpValidity;
         
         
 end
    

 %% % % % % Local function for calculating the EM field according to the Green's Function Approach
    function [EHx, EHy, EHz] = GreenFunctionApproach(x, y, z, EMfield_Components)
        
        Z0 = SingleAntennaParameters.Z0;

        x1 = SingleAntennaParameters.Current_density(:,1) * 1E-3;
        y1 = SingleAntennaParameters.Current_density(:,2) * 1E-3;
        z1 = SingleAntennaParameters.Current_density(:,3) * 1E-3;
       
        Jx = SingleAntennaParameters.Current_density(:,4) + ...
            1i * SingleAntennaParameters.Current_density(:,7);
        Jy = SingleAntennaParameters.Current_density(:,5) + ...
            1i * SingleAntennaParameters.Current_density(:,8);
        Jz = SingleAntennaParameters.Current_density(:,6) + ...
            1i * SingleAntennaParameters.Current_density(:,9);
        SurfaceElementArea = SingleAntennaParameters.Current_density(:,10) * 1E-6;

        
        xyz_size = size(x .* y .* z);
        [X, X1] = meshgrid(x .* ones(xyz_size), x1);
        [Y, Y1] = meshgrid(y .* ones(xyz_size), y1);
        [Z, Z1] = meshgrid(z .* ones(xyz_size), z1);

        R = sqrt( (X - X1).^2 + (Y - Y1).^2 + (Z - Z1).^2);
        Green = exp(-1i * k .* R) ./ (4 * pi * R);
        switch EMfield_Components
            case {'E_field'}
                Gfun1 = 1 - 1i * (k.*R).^(-1) - (k.*R).^(-2);
                Gfun22 = ( -1 + 3i * (k.*R).^(-1) + 3 * (k.*R).^(-2) ) .* ...                     
                     ( ( bsxfun(@times,(X - X1),Jx) + bsxfun(@times,(Y - Y1),Jy) + bsxfun(@times,(Z - Z1), Jz) ) ./ R );
                EHx = reshape(sum( - 1i * k * Z0 .* Green .* bsxfun(@times,( bsxfun(@times,Gfun1,Jx) + bsxfun(@times,Gfun22, ( (X - X1) ./ R )) ), ...
                    SurfaceElementArea )), xyz_size);
                EHy = reshape(sum( - 1i * k * Z0 .* Green .* bsxfun(@times,( bsxfun(@times,Gfun1,Jy) + bsxfun(@times,Gfun22 ,( (Y - Y1) ./ R ) )), ...
                    SurfaceElementArea )), xyz_size);
                EHz = reshape(sum( - 1i * k * Z0 .* Green .* bsxfun(@times,( bsxfun(@times,Gfun1, Jz) + bsxfun(@times,Gfun22,( (Z - Z1) ./ R ) )), ...
                    SurfaceElementArea )), xyz_size);
            case {'H_field'}
                EHx = reshape(sum( 1i * k * (1 - 1i * ( k * R ).^(-1)) .* Green .* ...
                   bsxfun(@times,( ( bsxfun(@times,(Z - Z1),Jy) - bsxfun(@times,(Y - Y1), Jz) ) ./ R ), SurfaceElementArea) ), xyz_size);
                EHy = reshape(sum( 1i * k * (1 - 1i * ( k * R).^(-1)) .* Green .* ...
                    bsxfun(@times,( ( bsxfun(@times,(X - X1), Jz) - bsxfun(@times,(Z - Z1), Jx) ) ./ R ), SurfaceElementArea) ), xyz_size);
                EHz = reshape(sum( 1i * k * (1 - 1i * (k * R).^(-1)) .* Green .* ...
                    bsxfun(@times,( ( bsxfun(@times,(Y - Y1),Jx) - bsxfun(@times,(X - X1), Jy) ) ./ R ), SurfaceElementArea) ), xyz_size);
                
        end

    end

%% % % % Local function for calculating the EM field using to the Multipole expansion technique
    function [EHx, EHy, EHz] = MultipoleExpansion(x, y, z, EMfield_Components)
        
        r = sqrt(x.^2 + y.^2 + z.^2);
        rho = sqrt(x.^2 + y.^2);
        theta = acos( z ./ ( 0.01 .* (r <= 0) + (r) .* (r > 0) ) ) .* (r > 0);
        phi = ( acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) .* ( y >= 0 ) + ...
            (2*pi - acos( x ./ ( 0.01 .* (rho <= 0) + (rho) .* (rho > 0) ) ) ) .* ( y < 0 ) ) .* (rho > 0);
        
        kr = k .* r; 
                
        EHx = 0;
        EHy = 0;
        EHz = 0;
        
        AE_lm = SingleAntennaParameters.AE_lm;
        AM_lm = SingleAntennaParameters.AM_lm;
                
        for l = 1 : SingleAntennaParameters.lmax
            for m = -l : l
                
                % spherical Hankel functions
                hl = sqrt(pi ./ (2 * kr)) .* besselh(l + 1/2, 2, kr);
                hl_plus1 = sqrt(pi ./ (2 * kr)) .* besselh(l + 1 + 1/2, 2, kr);
                
                % we use the recurence relation to calculate the derivative of the
                % Associated Legendre polinomials with respect to the variable \theta
                m_minus1 = (abs(m) - 1) .* (abs(m) > 0);
                m_plus1 = (abs(m) + 1);
                
                % we calculate the associated Legendre polynomials
                LegendrePolinom_lm = legendre(l, cos(theta));
                Plm = reshape(LegendrePolinom_lm( abs(m) + 1,:,:,:), size(theta));
                Plm_minus1 = reshape(LegendrePolinom_lm( (m_minus1 + 1),:,:,:), ...
                    size(theta)) .* (abs(m) > 0);
                Plm_plus1 = reshape(LegendrePolinom_lm( (m_plus1 .* (m_plus1 <= l) + 1),:,:,:), ...
                    size(theta)) .* ( 1 .* (m_plus1 <= l) + 1 .* (m == 0));
                Plm_over_theta = 0.5 * ( Plm_plus1  - (l + abs(m)) * (l - abs(m) + 1) * Plm_minus1);
                % Plm_over_sintheta = (Plm ./ sin(theta));
                clear Plm_minus1 Plm_plus1
                
                Alm = sqrt( ((2*l + 1 )/(4*pi)) * ( factorial(l - abs(m)) / factorial(l + abs(m)) ) );
                % Alm = Alm * ( (-1)^m .* (m < 0) + 1 .* (m >= 0) );
                Clm = ( Alm / (sqrt(l^2 + l)) ) * ( (-1)^m .* (m < 0) + 1 .* (m >= 0) );
                
                % we calculate the vector harmonics
                Mlm_r = 0;
                Mlm_theta = - Clm  * m * hl .* (Plm ./ sin(theta)) .* exp(-1i * m * phi);
                Mlm_phi = 1i * Clm  * hl .* Plm_over_theta .* exp(-1i * m * phi);

                Nlm_r = -1i * Clm  * ( l^2 + l) * (hl ./ kr) .* Plm .* exp(-1i * m * phi);
                Nlm_theta = -1i * Clm  * ( (l + 1) * hl ./ kr - hl_plus1) .* Plm_over_theta .* exp(-1i * m * phi);
                Nlm_phi = - Clm  * ( (l + 1) * hl ./ kr - hl_plus1) * m .* (Plm ./ sin(theta)) .* exp(-1i * m * phi);
                
                % we calculate the Cartesian components of the vector harmonics using the spherical-to-rectangular 
                % transformation formulas
                Mlm_x = Mlm_r .* sin(theta) .* cos(phi) + Mlm_theta .* cos(theta) .* cos(phi) - Mlm_phi .* sin(phi);
                Mlm_y = Mlm_r .* sin(theta) .* sin(phi) + Mlm_theta .* cos(theta) .* sin(phi) + Mlm_phi .* cos(phi);
                Mlm_z = Mlm_r .* cos(theta) - Mlm_theta .* sin(theta);
                Nlm_x = Nlm_r .* sin(theta) .* cos(phi) + Nlm_theta .* cos(theta) .* cos(phi) - Nlm_phi .* sin(phi);
                Nlm_y = Nlm_r .* sin(theta) .* sin(phi) + Nlm_theta .* cos(theta) .* sin(phi) + Nlm_phi .* cos(phi);
                Nlm_z = Nlm_r .* cos(theta) - Nlm_theta .* sin(theta);
                
                % we calculate the Cartesian components of the EM field
                switch EMfield_Components
                    case {'E_field'}
                        EHx = EHx + ( -1i * AE_lm(l, l + m + 1) .* Nlm_x + AM_lm(l, l + m + 1) .* Mlm_x ) * Z0;
                        EHy = EHy + ( -1i * AE_lm(l, l + m + 1) .* Nlm_y + AM_lm(l, l + m + 1) .* Mlm_y ) * Z0;
                        EHz = EHz + ( -1i * AE_lm(l, l + m + 1) .* Nlm_z + AM_lm(l, l + m + 1) .* Mlm_z ) * Z0;
                    case {'H_field'}
                        EHx = EHx + ( AE_lm(l, l + m + 1) .* Mlm_x + 1i * AM_lm(l, l + m + 1) .* Nlm_x );
                        EHy = EHy + ( AE_lm(l, l + m + 1) .* Mlm_y + 1i * AM_lm(l, l + m + 1) .* Nlm_y );
                        EHz = EHz + ( AE_lm(l, l + m + 1) .* Mlm_z + 1i * AM_lm(l, l + m + 1) .* Nlm_z );
                end
                
            end
        end
        

   
        
        
        
    end

end