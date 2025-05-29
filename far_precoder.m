function w = far_precoder( far, m, n )
  
  wh_mtx = dftmtx(far.n);
  
  wh = wh_mtx(n + far.offset_n,:);

  wv_mtx = dftmtx(far.m);
  
  wv = wv_mtx(:,m + far.offset_m);

  w = 1/sqrt(far.n*far.m) * kron(wh,wv);
  
end