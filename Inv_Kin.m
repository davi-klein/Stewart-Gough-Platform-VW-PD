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
%	Function file responsible for calculating the inverse kinematics of   %
%    a SPS Stewart-Gough platform                                         %
%-------------------------------------------------------------------------%

function [l,s] = Inv_Kin(P,a,b,RBA)

    L = zeros(1,6);
    S = zeros(3,6);

    RBAb = RBA*b;

    for i = 1:6
        L(i) = norm(P + RBAb(:,i) - a(:,i)); 
        S(:,i) = (P + RBAb(:,i) - a(:,i))/L(i);
    end
    
    l = L;
    s = S;
    
end

