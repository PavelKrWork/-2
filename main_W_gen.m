function main_W_gen( file_name, step_number)

   load([file_name mat2str(step_number)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul','Test');
   
    
   for t = 1 : numel(Test)
       
       if prec.comp.wizard
           
           alpha = cell2mat(Test(t).alpha.alpha);
           
       else
           
           alpha = Test(t).alpha.alpha_nonw;
    
       end
       
        [iy, ix, iz] = find_grid_ind(grid, ant, Test(t).pim);
        
        for k = 1 : numel(Test(t).pim)

            v(:,:,k) = squeeze(prec.v_saved_f2(iy(k),ix(k),iz(k),:,:));

        end
        
        for ah = 1 : far.N

            for av = 1 : far.M
                
                [W, norm_W] = W_gen( far, alpha(av,ah), v, Test(t).pim, av, ah );

                Test(t).W(:,:,av,ah) = W;

            end

        end       
                            
    end

    save([file_name mat2str(step_number)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul', 'Test');    

end