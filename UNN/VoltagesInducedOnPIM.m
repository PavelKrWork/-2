function [ Voltage_f1, Voltage_f2, MutualImpedaces_f3 ] = ...
    VoltagesInducedOnPIM( f_1, f_2, f_3, SingleAntennaParameters, AntennaArrayParameters, PIMsourceParameters )
% In this function, we calculate the voltages Voltage_f1 and Voltage_f2, induced by AA elements
% at the PIMs sources at the frequencies f_1 and f_2, 
% and the mutual impedance MutualImpedaces_f3 between the antennas and PIM-sources at f_3.

% The arrays Voltage_f1, Voltage_f2, MutualImpedaces_f3 have size [N_PIM x N_array], 
% where N_PIM - number of PIMs, N_array - number of antenna elements

global c Z0

% % normalization of the vector determing the PIM-source orientation 
PIM_orientation = bsxfun(@rdivide,PIMsourceParameters.PIM_orientation, ...
     sqrt(sum(abs(PIMsourceParameters.PIM_orientation).^2,2)));




%% f = f3

    % we calculate the normalization factor to calculate the mutual impedance
    switch SingleAntennaParameters.CalculationTechnique % 'PointDipole' 'LinearDipole' 'Green' 'Combined_Green+Multipole'
        case {'PointDipole', 'LinearDipole'}
            NormalizationFactor_f3 = -SingleAntennaParameters.I0;
            
        case {'Green', 'Multipole', 'Combined_Green+Multipole'}
            [ AE_lm, AM_lm] = Multipole_Coefficients( f_3, SingleAntennaParameters);
            SingleAntennaParameters.AE_lm = AE_lm;
            SingleAntennaParameters.AM_lm = AM_lm;
            NormalizationFactor_f3 =  sqrt( (2 / real(SingleAntennaParameters.Z0)) * ...
                ( (Z0 / (2*(2 * pi * f_3 / c)^2)) * sum(sum( (abs(AE_lm)).^2 + (abs(AM_lm)).^2 )) ) );
    end
    
    MutualImpedaces_f3 = Voltage(f_3);
    MutualImpedaces_f3 = MutualImpedaces_f3 / NormalizationFactor_f3;


%% f = f1
global c Z0
    
    % we calculate the narmalization factor for the EM field of the single antenna at the frequency f = f1
    switch SingleAntennaParameters.CalculationTechnique % 'PointDipole' 'LinearDipole' 'Green' 'Combined_Green+Multipole'
        case {'PointDipole'}
            R_rad = 80 * pi^2 * (SingleAntennaParameters.length * f_1 / c)^2;
            NormalizationFactor_f1 = sqrt(abs(SingleAntennaParameters.I0)^2 * R_rad / 2);  
            
        case {'LinearDipole'}
            R_rad = 20 * pi^2 * (SingleAntennaParameters.length * f_1 / c)^2;
            NormalizationFactor_f1 = sqrt(abs(SingleAntennaParameters.I0)^2 * R_rad / 2);
            
        case {'Green', 'Multipole', 'Combined_Green+Multipole'}
            [ AE_lm, AM_lm ] = Multipole_Coefficients( f_1, SingleAntennaParameters);
            SingleAntennaParameters.AE_lm = AE_lm;
            SingleAntennaParameters.AM_lm = AM_lm;
            P_rad = (Z0 / (2 * (2*pi*f_1/c)^2)) * sum(sum( (abs(AE_lm)).^2 + (abs(AM_lm)).^2 )); % radiated power for the single antenna
            NormalizationFactor_f1 = sqrt( P_rad );
    end

    
    Voltage_f1 = Voltage(f_1);
    Voltage_f1 = Voltage_f1 / NormalizationFactor_f1;


%% f = f2
global c Z0

    % we calculate the narmalization factor for the EM field of the single antenna at the frequency f = f1
    switch SingleAntennaParameters.CalculationTechnique % 'PointDipole' 'LinearDipole' 'Green' 'Combined_Green+Multipole'
        case {'PointDipole'}
            R_rad = 80 * pi^2 * (SingleAntennaParameters.length * f_2 / c)^2;
            NormalizationFactor_f2 = sqrt(abs(SingleAntennaParameters.I0)^2 * R_rad / 2);  
            
        case {'LinearDipole'}
            R_rad = 20 * pi^2 * (SingleAntennaParameters.length * f_2 / c)^2;
            NormalizationFactor_f2 = sqrt(abs(SingleAntennaParameters.I0)^2 * R_rad / 2);
            
        case {'Green', 'Multipole', 'Combined_Green+Multipole'}
            [ AE_lm, AM_lm ] = Multipole_Coefficients( f_2, SingleAntennaParameters);
            SingleAntennaParameters.AE_lm = AE_lm;
            SingleAntennaParameters.AM_lm = AM_lm;
            P_rad = (Z0 / (2 * (2*pi*f_2/c)^2)) * sum(sum( (abs(AE_lm)).^2 + (abs(AM_lm)).^2 )); % radiated power for the single antenna
            NormalizationFactor_f2 = sqrt( P_rad );
    end

    Voltage_f2 = Voltage(f_2);
    Voltage_f2 = Voltage_f2 / NormalizationFactor_f2;



%% Local function for calculating the voltages induced on the PIM sources at a given frequency f
    function Voltage_f = Voltage(f)
      
   
        k = 2 * pi * f / c;
        Voltage_f = zeros(PIMsourceParameters.N_PIM, AntennaArrayParameters.N_array);
                        
        for ii = 1 : AntennaArrayParameters.N_array
            [ Ex, Ey, Ez] = SingleAntennaFieldEH(  PIMsourceParameters.XYZ_PIM(:, 1) - AntennaArrayParameters.XYZ_array(ii, 1) , ...
                PIMsourceParameters.XYZ_PIM(:, 2) - AntennaArrayParameters.XYZ_array(ii, 2), ...
                PIMsourceParameters.XYZ_PIM(:, 3) - AntennaArrayParameters.XYZ_array(ii, 3), ...
                f, SingleAntennaParameters, 'E_field' );
            
            Ep = (Ex .* PIM_orientation(:, 1) + ...
                  Ey .* PIM_orientation(:, 2) + ...
                  Ez .* PIM_orientation(:, 3));
            
            Voltage_f(:, ii) = - Ep .* (2 ./ k) .* tan(k .* PIMsourceParameters.L_PIM / 4);
        end
        
    end

end

