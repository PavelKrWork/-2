%h = PSDoFreqplot(Data,Color,Fs)    
function varargout = PSDoFreqplot(Data,Color,Fs,BandWidthMHz)    
    if nargin == 3
        if length(Data)<500
            [Pxx,Freq] = freqz(Data);
            lgP = 20*log10(abs(Pxx));
        else
%            [Pxx Freq] = psd(Data);
            [Pxx Freq] = pwelch(Data);
                        lgP = 10*log10(fftshift(Pxx));
        end
        maxP = max(lgP)*(nargout>0);
        h = plot(Freq/1000000,lgP-maxP,Color);
        xlabel('Frequency(MHz)');
        ylabel('Power/frequency(dB/Hz)');
        hold on;grid on;
    elseif nargin == 4
        if length(Data)<500
            [Pxx,Freq] = freqz(Data,1,4096,'whole',Fs);
            lgP = 20*log10(abs(fftshift(Pxx)));
        else
%            [Pxx Freq] = psd(Data,4096,Fs);
%            [Pxx Freq] = periodogram(Data,hamming(length(Data)),4096,Fs);
            [Pxx Freq] = pwelch(Data,[],[],4096,Fs);
            lgP = 10*log10(fftshift(Pxx));
        end   
        maxP  = max(lgP)*(nargout>0);
        meanP = mean(10*log10(Pxx(1:200)));
        h=0;
        h = plot((Freq-Fs/2)/1000000,lgP-maxP,Color);
        xlabel('Frequency(MHz)');
        ylabel('Power/frequency(dB/Hz)');
        hold on;grid on;
    else
        error('input must be 2 or 3.');
    end
    LowFreqDot  = round(2048-BandWidthMHz/2/((max(Freq)-min(Freq))/1e6/4095));
    HighFreqDot = round(2049+BandWidthMHz/2/((max(Freq)-min(Freq))/1e6/4095));
    StopBanddB  = (lgP(LowFreqDot,1)-meanP+lgP(HighFreqDot,1))/2;
    varargout{1} = StopBanddB;
    varargout{2} = h;
    varargout{3} = lgP;
    varargout = varargout(1:nargout);   
end
