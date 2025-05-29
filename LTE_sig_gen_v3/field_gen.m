%function E_ant = field_gen(Gx,Gy,Gz,ant,dl,w,A)
function E_ant = field_gen(Gx,Gy,Gz,ant,freq,w,A)

  E_ant = zeros([ant.M ant.N]);

  A = A/sqrt(ant.N*ant.M); % amplitude
  
  lambda = physconst('LightSpeed')/freq;

  k0 = 2*pi/lambda; % wave number
    
  for ah = 1 : ant.N
      
    for av = 1 : ant.M

        R = sqrt((ant.x(av,ah)-Gx).^2 + (ant.y(av,ah)-Gy).^2 + (ant.z(av,ah)-Gz).^2);

        E_ant(av,ah) = A./R.*exp(-1i*(k0.*R))*w(av,ah); %*exp(-1i*degtorad(0));

    end
    
  end

end
