function [iy, ix, iz] = find_grid_ind(grid, ant, pim)

    [Gx,Gy,Gz] = grid_gen(grid, ant); % grid gen
    
    iy = zeros(numel(pim),1);        
    ix = zeros(numel(pim),1);
    iz = zeros(numel(pim),1);
    
    for k = 1 : numel(pim)

%        iy(k) = find(pim(k).y == Gy(:,1,1));        
%        ix(k) = find(pim(k).x == Gx(1,:,1));
%        iz(k) = find(pim(k).z == Gz(1,1,:));
        iy(k) = find(abs(pim(k).y - Gy(:,1,1))<=eps);        
        ix(k) = find(abs(pim(k).x - Gx(1,:,1))<=eps);
        iz(k) = find(abs(pim(k).z - Gz(1,1,:))<=eps);

    end

end