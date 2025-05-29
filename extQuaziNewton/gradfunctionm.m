function Return = gradfunctionm(functname,x,nwzfp)
%
% numerical computation of gradient
% this allows automatic gradient computation
% 
%
% first forward finite difference
% hstep = 0.001; - programmed in
%
hstep = 0.001;
n = length(x);
f = feval(functname,x,nwzfp.Test,nwzfp.LTE_sb_2_up2,nwzfp.LTE_sb_5_up2,nwzfp.grid,nwzfp.ant,nwzfp.prec,nwzfp.far,nwzfp.dl2,nwzfp.ul,nwzfp.av,nwzfp.ah,nwzfp.offset);
for i = 1:n
   xs = x;
   xs(i) = xs(i) + hstep;
   gradx(i)= (feval(functname,xs,nwzfp.Test,nwzfp.LTE_sb_2_up2,nwzfp.LTE_sb_5_up2,nwzfp.grid,nwzfp.ant,nwzfp.prec,nwzfp.far,nwzfp.dl2,nwzfp.ul,nwzfp.av,nwzfp.ah,nwzfp.offset) -f)/hstep;
end
Return = gradx;
