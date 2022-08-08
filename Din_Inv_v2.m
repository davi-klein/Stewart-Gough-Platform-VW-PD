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
%   Function file responsible for the calculation of the Inverse Dynamics %
%   by means of the Virtual Work Principle of a SPS Stewart-Gough platform%
%-------------------------------------------------------------------------%

function [tau, din_mat] = Din_Inv_v2(L,s,RBA,dP,ddP,dea,ddea,bb,e1,e2,m1,m2,I1ii,I2ii,mp,Ipb,fe,ne)

b = RBA*bb;
Ip = RBA*Ipb*RBA';

%-------------------------- Velocity Analysis ----------------------------%

RCA = zeros(3,3,6);

cti = zeros(1,6);
sti = zeros(1,6);
spi = zeros(1,6);
cpi = zeros(1,6);

for i = 1:6
    cti(i) = s(3,i);
    sti(i) = sqrt((s(1,i)^2)+(s(2,i)^2));
    spi(i) = s(2,i)/sti(i);
    cpi(i) = s(1,i)/sti(i);
    
    RCA(:,:,i) = [cpi(i)*cti(i) -spi(i) cpi(i)*sti(i);...
                  spi(i)*cti(i) cpi(i)  spi(i)*sti(i);...
                  -sti(i) 0 cti(i)]; 
end

vbi = zeros(3,6);
vbii = zeros(3,6);
wbii = zeros(3,6);
vp = dP; 
wp = dea;
v1i = zeros(3,6);
v2i = zeros(3,6);
v1ii = zeros(3,6);
v2ii = zeros(3,6);
wbi = zeros(3,6);

for i = 1:6
    vbi(:,i) = vp + (cross(wp,b(:,i)));
    vbii(:,i) = RCA(:,:,i)'*vbi(:,i);
    
    v1ii(:,i) = (e1/(L(i)))*[vbii(1,i);vbii(2,i);0];
    v2ii(:,i) = (1/(L(i)))*[(L(i)-e2)*vbii(1,i);(L(i)-e2)*vbii(2,i);L(i)*vbii(3,i)];
    
    v1i(:,i) = RCA(:,:,i)*v1ii(:,i);
    v2i(:,i) = RCA(:,:,i)*v2ii(:,i);
end

for i = 1:6
    wbii(:,i) = ((1/L(i))*[-vbii(2,i); vbii(1,i); 0]);
    wbi(:,i) = RCA(:,:,i)*wbii(:,i); 
end
dX = [vp;wp];

%----------------------- Acceleration Analysis ---------------------------%

dvp = ddP;
dwp = ddea;

dvbi = zeros(3,6);
dvbii = zeros(3,6);
dwbii = zeros(3,6);
dwbi = zeros(3,6);

dv1ii = zeros(3,6);
dv2ii = zeros(3,6);

dv1i = zeros(3,6);
dv2i = zeros(3,6);

for i = 1:6
    dvbi(:,i) = dvp + (cross(dwp,b(:,i))) + cross(wp,cross(wp,b(:,i)));
    dvbii(:,i) = RCA(:,:,i).'*dvbi(:,i);
    
    dwbii(:,i) = ((1/(L(i)^2))*[((-L(i))*dvbii(2,i))+(2*vbii(3,i)*vbii(2,i));...
                                (L(i)*dvbii(1,i))- (2*vbii(3,i)*vbii(1,i));...
                                0]);
    
    dwbi(:,i) = RCA(:,:,i)*dwbii(:,i);
                           
    dv1ii(:,i) = ((e1/(L(i)^2))*[(L(i)*dvbii(1,i)) - (2*vbii(3,i)*vbii(1,i));...
                                 (L(i)*dvbii(2,i)) - (2*vbii(3,i)*vbii(2,i));...
                                -((vbii(1,i)^2) + (vbii(2,i)^2))]);
                          
    dv2ii(:,i) = ((1/(L(i)^2))*[(L(i)*(L(i)- e2)*dvbii(1,i)) + (2*e2*vbii(3,i)*vbii(1,i));...
                                (L(i)*(L(i)- e2)*dvbii(2,i)) + (2*e2*vbii(3,i)*vbii(2,i));...
                                ((L(i)^2)*dvbii(3,i)) + (e2*((vbii(1,i)^2) + (vbii(2,i)^2)))]);
    
    dv1i(:,i) = RCA(:,:,i)*dv1ii(:,i);
    dv2i(:,i) = RCA(:,:,i)*dv2ii(:,i);
end

ddX = [dvp;dwp];

ddX1 = [dv1i;dwbi];
ddX2 = [dv2i;dwbi];
%-------------------------- Jacobian Matrices ----------------------------% 

Jbi = zeros(3,6,6);
Jbii = zeros(3,6,6);
Jp = zeros(6,6);

for i = 1:6
    Jbi(:,:,i) = [1 0 0 0 b(3,i) -b(2,i);... 
                  0 1 0 -b(3,i) 0 b(1,i);...
                  0 0 1 b(2,i) -b(1,i) 0];
    Jbii(:,:,i) = RCA(:,:,i).'*Jbi(:,:,i);
    
    Jp(i,:) = Jbii(3,:,i);
end

% ----------------- Closed form dynamics of the limbs ------------------- %

g = [0;0;-9.8067];

s_hat = zeros(3,6);
sx = zeros(3,3,6);
s_hatx = zeros(3,3,6);

Mp =  [mp*eye(3) zeros(3); zeros(3) Ip];
Cp = [zeros(3) zeros(3); zeros(3) skws(wp)*Ip];
Gp = [-mp*g; zeros(3,1)];

M1 = zeros(6,6,6);
M2 = zeros(6,6,6);
Melos = zeros(6,6,6);
 
C1 = zeros(6,6,6);
C2 = zeros(6,6,6);
Celos = zeros(6,6,6); 

G1 = zeros(6,1,6);
G2 = zeros(6,1,6);
Gelos = zeros(6,1,6);

I1iA = zeros(3,3,6);
I2iA = zeros(3,3,6);

J1ih = zeros(6,6,6);
J2ih = zeros(6,6,6);

dJ1ih = zeros(6,6,6);
dJ2ih = zeros(6,6,6);

sumMelos = zeros(6,6);
sumCelos = zeros(6,6);
sumGelos = zeros(6,1);

for i = 1:6
    s_hat(:,i) = s(:,i)/(norm(s(:,i))); 
    sx(:,:,i) = skws(s(:,i));
    s_hatx(:,:,i) = skws(s_hat(:,i));
   
    J1ih(:,:,i) = (1/L(i))*[-e1*(s_hatx(:,:,i)^2) e1*((s_hatx(:,:,i)^2)*skws(b(:,i)));...
                             s_hatx(:,:,i) -s_hatx(:,:,i)*skws(b(:,i))];
                  
    J2ih(:,:,i) = (1/L(i))*[-e2*(s_hatx(:,:,i)^2)+(L(i)*s_hat(:,i)*s_hat(:,i).')...
                            e2*(s_hatx(:,:,i)^2)*skws(b(:,i))-(L(i)*s_hat(:,i)*s_hat(:,i).'*skws(b(:,i)));...
                             s_hatx(:,:,i)  -s_hatx(:,:,i)*skws(b(:,i))];
    % Lemus et. al%

    I1iA(:,:,i) = RCA(:,:,i)*I1ii(:,:,i)*RCA(:,:,i).';
    I2iA(:,:,i) = RCA(:,:,i)*I2ii(:,:,i)*RCA(:,:,i).';
    
    M1(:,:,i) = [m1(i)*eye(3) zeros(3); zeros(3) I1iA(:,:,i)];
    M2(:,:,i) = [m2(i)*eye(3) zeros(3); zeros(3) I2iA(:,:,i)];
    
    Melos(:,:,i) = (J1ih(:,:,i).'*M1(:,:,i)*J1ih(:,:,i)) + (J2ih(:,:,i).'*M2(:,:,i)*J2ih(:,:,i));
    
    dJ1ih(:,:,i) = (ddX1(:,i)-(J1ih(:,:,i)*ddX))*pinv(dX);
    dJ2ih(:,:,i) = (ddX2(:,i)-(J2ih(:,:,i)*ddX))*pinv(dX);
    
    C1(:,:,i) = [zeros(3) zeros(3); zeros(3) skws(wbi)*I1iA(:,:,i)];
    C2(:,:,i) = [zeros(3) zeros(3); zeros(3) skws(wbi)*I2iA(:,:,i)];
    
    Celos(:,:,i) = (J1ih(:,:,i).'*C1(:,:,i)*J1ih(:,:,i)) + (J1ih(:,:,i).'*M1(:,:,i)*dJ1ih(:,:,i)) + ...
                   (J2ih(:,:,i).'*C2(:,:,i)*J2ih(:,:,i)) + (J2ih(:,:,i).'*M2(:,:,i)*dJ2ih(:,:,i));
               
    G1(:,:,i) = [-m1(i)'*g;zeros(3,1)];
    G2(:,:,i) = [-m2(i)'*g;zeros(3,1)]; 
    
    Gelos(:,:,i) = ((J1ih(:,:,i).'*G1(:,:,i)) + ((J2ih(:,:,i).'*G2(:,:,i))));
    
    sumMelos = sumMelos + Melos(:,:,i);
    sumCelos = sumCelos + Celos(:,:,i);
    sumGelos = sumGelos + Gelos(:,:,i);

end

din_mat.MX = Mp + sumMelos;
din_mat.CX = Cp + sumCelos;
din_mat.GX = Gp + sumGelos;

din_mat.Fd = [fe;ne];

din_mat.Fc = ((din_mat.MX*ddX)+(din_mat.CX*dX)+ din_mat.GX - din_mat.Fd); %cartesian forces

tau = inv(Jp).'*din_mat.Fc; %joint forces
end
