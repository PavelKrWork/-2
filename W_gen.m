function [W, norm_W] = W_gen( far, alpha, v, pim, av, ah )
                    
    w = far_precoder(far, av, ah);
                        
    V = 0;

    for k = 1 : numel(pim)

%        V = V + alpha{av,ah}(k)'.*v(:,:,k);
        V = V + alpha(k)'.*v(:,:,k);

    end

    W = w + V;
    
    norm_W = norm(W(:));

    W = 1/norm_W*W;

end