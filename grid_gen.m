function [Gx,Gy,Gz] = grid_gen(grid, ant)

  Gx_min = 0; % m

  Gx_max = grid.size.x; % m

  Gy_min = -grid.size.y; % m

  Gy_max = grid.size.y; % m

  Gz_min = -grid.size.z; % m

  Gz_max = grid.size.z; % m

  Gx_step = grid.step.x * ant.d_h;
  
  Gy_step = grid.step.y * ant.d_h;
  
  Gz_step = grid.step.z * ant.d_v;
  
  if grid.limit2antsize
            
      Gy_min = max(Gy_min,ant.y(1));
      
      Gy_max = min(Gy_max,ant.y(end));
      
      Gz_min = max(Gz_min,ant.z(1));
      
      Gz_max = min(Gz_max,ant.z(end));
      
  end
  
  if grid.size.x

    Gx_v = [ fliplr( ant.x(1)-Gx_step : -Gx_step : Gx_min )    ant.x(1) : Gx_step : Gx_max ];
    
  else
    
    Gx_v = 0;
    
  end
  
  if grid.size.y

    Gy_v = [ fliplr( ant.y(1)-Gy_step : -Gy_step : Gy_min )    ant.y(1) : Gy_step : Gy_max ];
    
  else
    
    Gy_v = 0;
    
  end
  
  if grid.size.z

    Gz_v = [ fliplr( ant.z(1)-Gz_step : -Gz_step : Gz_min )    ant.z(1) : Gz_step : Gz_max ];
    
  else
    
    Gz_v = 0;
    
  end


  [Gx,Gy,Gz] = meshgrid(Gx_v,Gy_v,Gz_v);
  
  Gx = Gx + grid.offset.x;
  Gy = Gy + grid.offset.y;
  Gz = Gz + grid.offset.z;


end
