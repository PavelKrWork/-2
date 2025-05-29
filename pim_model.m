%function out = pim_model(e_f1, e_f2, h_f1, h_f2, PIM, s_f1, s_f2)
function out = pim_model(e_f1, e_f2, h_f1, h_f2, d, s_f1, s_f2)

       e_f1 = squeeze(e_f1);

       e_f2 = squeeze(e_f2);
       
       e_f1 = reshape(e_f1(:),[3 length(e_f1(:))/3]);

       e_f2 = reshape(e_f2(:),[3 length(e_f2(:))/3]);
       
       e = e_f1 * s_f1(:) + e_f2 * s_f2(:);
       
%        out = e'*PIM.d;
       out = e'*d;
       
       out = out*out*conj(out);
               
end