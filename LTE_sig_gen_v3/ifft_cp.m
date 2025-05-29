function txSigOFDM = ifft_cp(symbolsIn, numFFT, gard_sub_loc, cpLen)

    numDataCarriers = length(symbolsIn);

    % Pack data into an OFDM symbol
    offset = numFFT - numDataCarriers;
    
    switch  gard_sub_loc
        
        case 'start'
                
            error('not supported');

        case 'end'
            
            error('not supported');
            
        case 'middle'
            
             symbolsInOFDM = [zeros(offset/2,1); symbolsIn; zeros(offset/2,1)];
            
        otherwise
            
            error('not supported')
            
    end
        
    symbolsInOFDM = ifftshift(symbolsInOFDM);
    
    ifftOut = ifft(symbolsInOFDM);

    % Prepend cyclic prefix
    txSigOFDM = [ifftOut(end-cpLen+1:end); ifftOut];
    %%%%TEST%%%%Set CP to zero
    %txSigOFDM( 1 : cpLen ) = zeros( cpLen, 1 );
    %%%%END TEST%%%%Prepend zeros
    
end
    
