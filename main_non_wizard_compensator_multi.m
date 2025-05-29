function main_non_wizard_compensator_multi(file_name, step_number,in_file_name_ofdm_sig_1,in_file_name_ofdm_sig_2)

    load(in_file_name_ofdm_sig_1,'LTE_sb_2_acc');

    load(in_file_name_ofdm_sig_2,'LTE_sb_5_acc');

    LTE_sb_2_up2 = resample(LTE_sb_2_acc, 2, 1, 10); % columns are resampled independently

    LTE_sb_5_up2 = resample(LTE_sb_5_acc, 2, 1, 10); % columns are resampled independently

    load([file_name mat2str(step_number-1)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul','Test');
        
    progressbar(['Step' mat2str(step_number)]);
  
    [X,Y]=meshgrid(prec.comp.multistart,prec.comp.multistart);
    
    UseParallel = false;
    
    if numel(X) > 1
        
        UseParallel = true;
        
    end    

    tpoints = CustomStartPointSet([X(:) Y(:)]);
    
    offset = 0;

    for k = 1 : numel(Test)
        
        xmin = zeros(ant.M,ant.N);
        
        fval = zeros(ant.M,ant.N);
        
        for ah = 1 : ant.N

            for av = 1 : ant.M
                
                 x_wizard = Test(k).alpha.alpha(av,ah);
                 
                 if prec.comp.known_offset
                 
                    offset = non_wizard_compensator_fun([real(cell2mat(x_wizard)) imag(cell2mat(x_wizard))],Test(k),LTE_sb_2_up2,LTE_sb_5_up2,grid,ant,prec,far,dl2,ul,av,ah,0);
                 
                 end
                    
                 f = @(x)non_wizard_compensator_fun(x,Test(k),LTE_sb_2_up2,LTE_sb_5_up2,grid,ant,prec,far,dl2,ul,av,ah,offset);
                 
                 problem = createOptimProblem('fminunc','objective',f,'x0',[0 0]);
                 
                 ms = MultiStart('UseParallel',UseParallel,'Display','iter');
                 
                 [xmin_tmp, fval(av,ah), exitflag, output, manymins] = run(ms, problem, tpoints);

                 xmin(av,ah) = xmin_tmp(1) + 1i*xmin_tmp(2);
    
            end

        end
                
        Test(k).alpha.alpha_nonw = xmin;
        
        Test(k).alpha.fval_nonw = fval;
              
        progressbar(k/numel(Test));
    
    end
    
    save([file_name mat2str(step_number)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul', 'Test');    
        
end