function out = rx_ul(in)

        in = fft(in);
        
        out = zeros(size(in));
        
        in = circshift(in, -1300); % shift by -78 Mhz
                
        out(1:324/2) = in(1:324/2);

        out(end-324/2+1:end) = in(end-324/2+1:end);
        
%          out = ifft(out);       
        
end