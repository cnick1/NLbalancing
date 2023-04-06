%  The KroneckerTools can be obtained from
%       https://github.com/jborggaard/KroneckerTools
%
% The tensor_toolbox can be obtained from
%  https://gitlab.com/tensors/tensor_toolbox.git
%  Download tensor_toolbox, then adjust the paths below to suit your location.

try
    addpath('../KroneckerTools/src')
    addpath('../KroneckerTools/util')
catch
    disp("Run 'git clone https://github.com/jborggaard/KroneckerTools' in the parent directory")
    error("KroneckerTools repo not found, please clone it.")
end

try
    addpath('../tensor_toolbox')
catch
    disp("Run 'git clone https://gitlab.com/tensors/tensor_toolbox.git' in the parent directory")
    error("tensor_toolbox repo not found, please clone it.")
end
