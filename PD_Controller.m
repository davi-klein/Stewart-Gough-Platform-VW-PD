%-------------------------------------------------------------------------%
%               UFSC - Federal University of Santa Catarina               %
%               Graduate Program in Mechanical Engineering                %
%                                                                         %
%     Programmer:                                                         %
%       Davi Klein                                                        %
%                                                                         %
%   Version: 1.0                                              08/09/2022  %
%=========================================================================%
%                          Program Descriprion                            %
%=========================================================================%
%	Function file responsible for the PD Control action calculation       %
%-------------------------------------------------------------------------%

function din_mat = PD_Controller(x,dx,X,dX,din_mat)
    Kp = 1e4*diag([1 1 1 10 10 10]);
    Kv = 1e3*diag([1 1 1 1 1 1]);
    din_mat.FJ = zeros(6,10);
    din_mat.PID = zeros(6,1);
    din_mat.FD = [0;0;0;0;0;0];
    
    din_mat.PID = Kp*(X-x)+Kv*(dX-dx);
    din_mat.FJ = pinv(din_mat.Jp')*din_mat.PID;
    din_mat.F = din_mat.FD + din_mat.PID;
end
