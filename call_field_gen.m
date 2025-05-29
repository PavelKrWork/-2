%function [E,H] = call_field_gen(pim,ant,dl,v,A)
function [E,H] = call_field_gen(pim,ant,freq,v,A)

%     if isempty(ant.TB)

    T = ant.T;

    switch ant.NF_model
        
        case 1 % Free space

            E = field_gen(pim.x,pim.y,pim.z,ant,freq,v,A);

            k = find(v(:));

            if numel(k) == 1

                E = sum(E(:));

            end                   

            H = [];
                
        case 2 % Matlab toolbox
        
            k = find(v(:));

            if numel(k) == 1

                ant.TB.PhaseShift = 0;

                ant.TB.AmplitudeTaper = 1;            

                k_val = ant.convMy2TB(k);
                
                [E, H] = EHfields(ant.TB, freq,[pim.x;pim.y;pim.z],'ElementNumber', k_val, 'Termination', T);
                
            elseif numel(k) > 1

                v_TB(ant.convMy2TB) = v;

                ant.TB.PhaseShift = rad2deg(angle(v_TB(:)));

                ant.TB.AmplitudeTaper = abs(v_TB(:));

                [E, H] = EHfields(ant.TB, freq,[pim.x;pim.y;pim.z]);
                        
            else 
            
                error('wrong near-field precoder');
               
            end
                        
        case 3 % UNN model
            
            SingleAntennaParameters  = Parameters();
            
            k = find(v(:));

            if numel(k) == 1

                [ Ex, Ey, Ez] = SingleAntennaFieldEH(  pim.x - ant.x(k) , pim.y - ant.y(k), pim.z - ant.z(k), ...
                                freq, SingleAntennaParameters, 'E_field' );
                 
                E = [ Ex; Ey; Ez];
                
                H = zeros(3,1);
                            
            elseif numel(k) > 1
                
                for kk = 1 : ant.M * ant.N
                    
                    [ Ex(kk), Ey(kk), Ez(kk)] = SingleAntennaFieldEH(  pim.x - ant.x(kk) , pim.y - ant.y(kk), pim.z - ant.z(kk), ...
                                    freq, SingleAntennaParameters, 'E_field' );
                                
                end
                
                E = [Ex(:).'*v(:);  Ey(:).'*v(:); Ez(:).'*v(:)];

                H = zeros(3,1);
                              
            else 
            
                error('wrong near-field precoder');

            
            end

    
        otherwise
            
            error('wrong NF model');
               
    end
        
end

