%-------------------------------------------------------------------------%
%               UFSM - Universidade Federal de Santa Maria                %
%              Curso de Engenharia de Controle e Automa��o                %
%                                                                         %
%     Programador:                                                        %
%       Davi Klein                                                        %
%                                                                         %
%   Vers�o: 1.0                                             30/10/2020    %
%=========================================================================%
%                        Descri��o do Programa                            %
%=========================================================================%
%	Function file responsible for calculating the inverse kinematics of   %
%    a SPS Stewart-Gough platform                                          %
%                                                                         %
%   v1.0 - Vers�o inicial.                                                %
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

