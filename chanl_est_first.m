function chanl_est_first(file_name, step_number)

  load([file_name mat2str(step_number-1)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul','Test');
  
 
  %% Reference H gen
  
  progressbar(['Step' mat2str(step_number)]);
  
  Href = zeros(size(ant.z,1)*size(ant.y,2), numel(Test(1).foc), numel(Test));
  
  for t = 1 : numel(Test)
      
      if numel(Test(t).foc)  == numel(Test(t).pim)
      
          Href_tmp = zeros(size(ant.z,1), size(ant.y,2), numel(Test(1).foc));

         for p = 1 : numel(Test(t).foc)

            for y = 1 : size(ant.y,2)

                for z = 1 : size(ant.z,1)

                    v = zeros(size(ant.z,1),size(ant.y,2));

                    v(z,y) = 1;

                    [e,~] = call_field_gen(Test(t).foc(p),ant,dl2.small_band.freq(1),v,1.0);

                    Href_tmp(z,y,p) = e'*Test(t).d{p};

                end

            end

         end

         Href(:,:,t) = reshape(Href_tmp, [size(ant.z,1)*size(ant.y,2)  numel(Test(1).foc)]);

      end
        
      progressbar(t/numel(Test));                

  end
  
  
  
  %% Estimated H
  
  p = 1;
  
  Hest = zeros(size(ant.z,1)*size(ant.y,2),numel(Test));
  
  for t = 1 : numel(Test)
      
    if numel(Test(t).foc)  == numel(Test(t).pim)
 
      A = reshape(Test(t).W,[16 16]);
      
      A = A';
      
     [maxval,linind]=max(Test(t).alpha.fval_nonw(:));
     
      A(linind,:) = [];
      
      A =([A;ones(1,size(A,2))]);

     Hest(:,t) = A\[zeros(size(A,1)-1,1);sum(Href(:,p,t))];
      
    end
    
  end
  
  %% Comparison
  
%  for t = 1 : numel(Test)
%      
%          if numel(Test(t).foc)  == numel(Test(t).pim)
%      
%            nrm(:,t) = vecnorm(Hest(:,t)-Href(:,p,t),2,2)./mean(vecnorm(Href(:,p,t),2,2)); % ver 2
%                
%                if 1
%
%                    Rp = floor(sqrt(numel(Test)));
%                    Cp = ceil(sqrt(numel(Test)));
%
%                    if Rp*Cp < numel(Test)
%                       Rp = Rp + 1;
%                    end
%
%                    [f,x]=ecdf(nrm(:,t));
%                    if t == 1 
%                        figure; end
%                    subplot(Rp,Cp,t); plot(10*log10(x),f); grid on;
%
%                end
%            
%          end
%           
%  end
%  
%%% plot CDF
%
%% %   nrm = nrm(:)./mean(vecnorm(Href(:),2,2)); %ver 3
%
%    [f,x]=ecdf(nrm(:));
%    
%    figure; plot(10*log10(x),f); grid on;
%    title('CDF of the Channel estimate EVM ');
%    xlabel('dB');
    
%% YZ metrics    
    
  for t = 1 : numel(Test)

      if numel(Test(t).foc)  == numel(Test(t).pim)
          
          Hest_n = Hest(:,t)/norm(Hest(:,t));
          
          Href_n = Href(:,p,t)/norm(Href(:,p,t));
          
          C(t) = abs(Hest_n'*Href_n);
          
          Cdb(t) = 10 * log10(1.0 - C(t)^2);

      end

  end
  
  cdfplot(Cdb(:));
  
%    [f,x]=ecdf(Cdb(:));
%
%    figure; plot(x,f); grid on;
%    title('CDF of the Cdb ');
%    xlabel('dB');
  
end
