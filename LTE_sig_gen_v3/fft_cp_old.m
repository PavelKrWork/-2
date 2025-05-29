function SymbolsOut = fft_cp_old( SigIn, numFFT, numDataCarriers, EvmCalcOffset, BandShift )

    cpLen = length(SigIn) - numFFT;

    SigInSymb = SigIn( cpLen + 1 : numFFT + cpLen );

    offset = numFFT - numDataCarriers;
       
    SigFFT = fft(SigInSymb);
    
    for_sin_vec = 0 : numFFT - 1 ;
                                    
    SigFFT = SigFFT .* exp( 1i * 2 * pi * (EvmCalcOffset * for_sin_vec.' / numFFT ) );
    
    SigFFT = fftshift(SigFFT);
        
    SymbolsOut = SigFFT( offset/2 + 1 : end - offset/2 );
                                                 
end
    
