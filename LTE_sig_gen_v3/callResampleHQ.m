function out = callResampleHQ( in, coeffs, res_mat )

%     out = resampleHQ( in, coeffs{ 1,1 }, 4, 3 );
%     
%     out = resampleHQ( out, coeffs{ 1,2 }, 2, 1 );
%     
%     out = resampleHQ( out, coeffs{ 1,3 }, 2, 1 );

    out = in;

    for k = 1 : size( coeffs, 2 )
        
         out = resampleHQ( out, coeffs{ 1, k }, res_mat( k, 1 ), res_mat( k, 2 ) );
        
    end        
    
    out = out(:);

end