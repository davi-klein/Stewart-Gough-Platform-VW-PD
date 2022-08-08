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
%	Function file responsible for the Decentralized PD Control routine    %
%-------------------------------------------------------------------------%

function cost =  Control(t,X_new,N,a,b,I1ii,I2ii,m1,m2,e1,e2,mp,Ipb)

    x = X_new(1:6);
    dx = X_new(7:12);

    w = 3.0;

    syms ti;

    P = [-1.5+(0.2*sin(w*ti)); 0.2*sin(w*ti);1.0+(0.2*sin(w*ti))];
    dP = diff(P);
    sdP = size(dP);

    if sdP(1) < 3
        dP = zeros(3,1);
    end

    ea = [0;0;0];
    dea = diff(ea);
    sdea = size(dea);

    if sdea(1) < 3
        dea = zeros(3,1);
    end

    ti = t;

    P = subs(P);
    dP = subs(dP);


    ea = subs(ea);
    dea = subs(dea);

    X = [P;ea];
    dX = [dP;dea];

    X = double(X);
    dX = double(dX);

    RBA = [cos(x(6))*cos(x(5)) cos(x(6))*sin(x(5))*sin(x(4))-sin(x(6))*cos(x(4)) cos(x(6))*sin(x(5))*cos(x(4))+sin(x(6)*sin(x(4)));...
           sin(x(6))*cos(x(5)) sin(x(6))*sin(x(5))*sin(x(4))+cos(x(6))*cos(x(4)) sin(x(6))*sin(x(5))*cos(x(4))-cos(x(6))*sin(x(4));...
           -sin(x(5)) cos(x(5))*sin(x(4)) cos(x(5))*cos(x(4))];

    [L,s] = Inv_Kin(x(1:3),a,b,RBA);
    
    [din_mat] = Closed_Form_Din(dx,b,L,s,m1,m2,e1,e2,mp,Ipb,I1ii,I2ii,RBA);
  
    [din_mat]= Controle(x,dx,X,dX,din_mat);
    
    cost = [dx; pinv(din_mat.MX)*(din_mat.F-din_mat.CX*dx-din_mat.GX)];
    
    
end