%function txSigOFDShft = slice_cont_flt( symbolsIn, numFFTLocal, BandShift, plotEnabled, plotFFTup, UpRate, DownRate, numFFTGlobal, HQ_flt, res_mat )        
function txSigOFDShft = slice_cont_flt( symbolsIn, txBufferPnt, numFFTLocal, cpLenGlobal, BandShift, plotEnabled, plotFFTup, UpRate, DownRate, numFFTGlobal, HQ_flt, res_mat, n_append, n_order )        


        txSigOFDM = symbolsIn;

        plot_psd(txSigOFDM, numFFTLocal*plotFFTup, plotEnabled);

        % Resampling
        
        if ~isempty( HQ_flt )
            
            txSigOFDM = [ txSigOFDM; zeros(n_append,1)];
            
            txSigOFDUpLong = callResampleHQ( txSigOFDM, HQ_flt, res_mat );                        
            
        else

%           txSigOFDUpLong = resample(txSigOFDM, UpRate, DownRate, n_order);

            txSigOFDUpLong = resample(txSigOFDM, 4, 3, n_order);

            txSigOFDUpLong = resample(txSigOFDUpLong, 2, 1, n_order);

            txSigOFDUpLong = resample(txSigOFDUpLong, 2, 1, n_order);
        
        end
        
        % Extract from the buffer               
        
        txSigOFDUp = txSigOFDUpLong(( txBufferPnt - 1 ) * ( numFFTGlobal + cpLenGlobal ) + 1 : txBufferPnt * ( numFFTGlobal + cpLenGlobal ) );
        
        plot_psd(txSigOFDUp, numFFTGlobal, plotEnabled);
        
        % Frequency shift
                
        for_sin_vec = 0 : 2 * length( txSigOFDUp ) - 1;
                                
        for_sin_vec = for_sin_vec * BandShift / numFFTGlobal;

        sin_samps = exp( -2 * pi * 1i * for_sin_vec ).';
        
        txSigOFDShft = [ zeros(size(txSigOFDUp)); txSigOFDUp ] .* sin_samps;
                                     
       plot_psd( txSigOFDShft, numFFTGlobal*plotFFTup, plotEnabled );

end       