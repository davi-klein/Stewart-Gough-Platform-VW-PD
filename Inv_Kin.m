%-------------------------------------------------------------------------%
%               UFSM - Universidade Federal de Santa Maria                %
%              Curso de Engenharia de Controle e Automação                %
%                                                                         %
%     Programador:                                                        %
%       Davi Klein                                                        %
%                                                                         %
%   Versão: 1.0                                             30/10/2020    %
%=========================================================================%
%                        Descrição do Programa                            %
%=========================================================================%
%	Function file responsible for calculating the inverse kinematics of   %
%    a SPS Stewart-Gough platform                                          %
%                                                                         %
%   v1.0 - Versão inicial.                                                %
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

