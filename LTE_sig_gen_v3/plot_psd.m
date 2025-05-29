function plot_psd(txSigOFDM, numFFT, plotEnabled)

    if plotEnabled
        
        txSigOFDM = txSigOFDM / sqrt( txSigOFDM' * txSigOFDM );
                       
        % Plot power spectral density (PSD) for OFDM signal
        [psd,f] = periodogram(txSigOFDM, rectwin(length(txSigOFDM)), numFFT, ...
                              1, 'centered');
        Fig1 = figure('Position', figposition([46 15 30 30]));
        plot(f,10*log10(psd));
        grid on
        %axis([-0.5 0.5 -100 -20]);
        xlabel('Normalized frequency');
        ylabel('PSD (dBW/Hz)')
        %title(['OFDM, ' num2str(numRBs*rbSize) ' Subcarriers'])
    
    end

end