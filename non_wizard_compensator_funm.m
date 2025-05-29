function svd_val_out = non_wizard_compensator_funm(al,Test,LTE_sb_2_up2,LTE_sb_5_up2,grid,ant,prec,far,dl2,ul,av,ah,offset)

    alpha = al(1) + 1i*al(2);
    
    load('resample_filter_1_4','b');
    
    [h,t] = impz(b,1);
    
    imp(t+1) = h;
   
    NOFDM = 1;
    
 % index of illuminating NF precoder
    x_set = 1;
    y_set = 1;
    z_set = 1;

    numFFTGlobal = 4096;

    [Gx,Gy,Gz] = grid_gen(grid, ant); % grid gen
            
        foc = Test.foc;

        E_foc_2 = zeros(size(Gx,1),size(Gx,2),size(Gx,3));

        x = x_set;
        
        y = y_set;
        
        z = z_set;
        
        m = 1; % num foc
        
%         av = 1;
% 
%         ah = 1;
%      

        [iy, ix, iz] = find_grid_ind(grid, ant, Test.pim);
        
        for p = 1 : numel(Test.pim)

            v(:,:,p) = squeeze(prec.v_saved_f2(iy(p),ix(p),iz(p),:,:));

        end
               
        [W, norm_W] = W_gen( far, alpha, v, Test.pim, av, ah );
        

        for m = 1 % : numel(foc)

 %         [E_f2, H_f2] = pim_by_W(foc(m).x, foc(m).y, foc(m).z, ant, dl2, f, far, W);
          
          [E_f2,H_f2] = call_field_gen(foc ,ant ,dl2.freq ,W ,1.0);

        end


        for m = 1 : numel(foc)

            E_ul(:,:,m) = field_gen(foc(m).x,foc(m).y,foc(m).z,ant,ul.freq,ones(ant.N,ant.M),1); %sqrt(ant.N*ant.M));

        end   

        for of = 1 : NOFDM

            tmp_acc = zeros(numFFTGlobal,numel(foc));

            for sp = 1 : numFFTGlobal

                for m = 1 : numel(foc)

                    tmp = pim_model(Test.E.E_foc_f1.e(:,y,x,z,m,:),E_f2,Test.E.E_foc_f1.h(:,y,x,z,m,:),Test.E.E_foc_f2.h(:,y,x,z,m,:),Test.d{m}, ...
                                        LTE_sb_2_up2(sp+(of-1)*numFFTGlobal,:), LTE_sb_5_up2(sp+(of-1)*numFFTGlobal,:)); %!!! Bug, takes the same sample for different symbols (should be sp+...)

                    tmp_acc(sp,m) = tmp;

                end

            end

          tmp_p = 0;
          
          tmp_ul_acc = zeros(numFFTGlobal/4, ant.M, ant.N);

          for ah_up = 1 : ant.N

                for av_up = 1 : ant.M

                    tmp = 0;

                    for m = 1 : numel(foc)

                        tmp = tmp + E_ul(av_up,ah_up,m)*tmp_acc(:,m);

                    end

                    tmp_ul = rx_ul(tmp); % tmp_ul is in freq domain
                    
%                  tmp_ul_acc(:,av_up,ah_up) = resample(tmp_ul,1,4);
                   tmp_ul_acc(:,av_up,ah_up) = resample(tmp_ul,1,4,imp); 

                end

          end
          
          Ruu_tmp =  reshape(tmp_ul_acc,[size(tmp_ul_acc,1) ant.N * ant.M ]);

          svd_val = sum(diag(Ruu_tmp*Ruu_tmp'));

          svd_val_out = abs(svd_val(1)-offset);
        
        end
            
end