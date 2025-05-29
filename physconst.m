function out = physconst(in)
  
  switch in
    
    case 'LightSpeed'
      
      out = 299792458;
      
      otherwise
      
      error('not defined');
      
      end
  
  end
  