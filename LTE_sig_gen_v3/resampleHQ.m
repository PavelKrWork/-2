function out = resampleHQ( in, coeffs, upRate, downRate )

   out_tmp_up = zeros( 1, length( in ) * upRate );
   
   out_tmp_up( 1 : upRate : end ) = in * upRate;
   
   out_tmp_flt = conv( out_tmp_up, coeffs );
   
   fltOrderTmp = floor( 0.5 * (length(coeffs) + 1 ));
      
   out = out_tmp_flt( fltOrderTmp : downRate : fltOrderTmp + (length( in ) * upRate / downRate - 1) * downRate );

end