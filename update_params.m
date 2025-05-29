function update_params(file_name, step_number_in,step_number_out)

  load([file_name mat2str(step_number_in)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul','Test');
%   load([file_name mat2str(step_number)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul');

    ant.NF_model = 3;
    
   save([file_name mat2str(step_number_out)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul','Test');
%  save([file_name mat2str(step_number)], 'ant', 'dl', 'dl2', 'grid', 'prec', 'far', 'lte', 'PIM', 'ul');
  
end
