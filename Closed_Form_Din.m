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
%	Function file responsible for the calculation of the Closed-Form      %
%   equations of motion                                                   %
%-------------------------------------------------------------------------%

function din_mat = Closed_Form_Din(dx,bb,L,s,m1,m2,e1,e2,mp,Ipb,I1ii,I2ii,RBA)


b = RBA*bb;
Ip = RBA*Ipb*RBA';


%----------------------- Análise de velocidades --------------------------%

RCA = zeros(3,3,6);% matriz de rotação que transforma o sistema de referen
                   % cia C em A 

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

vbi = zeros(3,6);% vetor que armazena os valores de velocidade de cada junta Bi em relação a A
vbii = zeros(3,6);% vetor que armazena os valores de velocidade de cada junta Bi em relação a C
wbii = zeros(3,6);% vetor que armazena os valores de velocidade angular de cada junta Bi em relação a C
vp = dx(1:3); %velocidade linear do ponto P (velocidade linear da plataforma movel)
wp = dx(4:6); %velocidade angular do ponto P (velocidade angular da plataforma movel)
wbi = zeros(3,6);

for i = 1:6
    vbi(:,i) = vp + (cross(wp,b(:,i)));
    vbii(:,i) = RCA(:,:,i)'*vbi(:,i);
end

for i = 1:6
    wbii(:,i) = ((1/L(i))*[-vbii(2,i); vbii(1,i); 0]);
    wbi(:,i) = RCA(:,:,i)*wbii(:,i); 
end

%--------------- Definições das matrizes Jacobianas dos ------------------% 
%------------------------ elos e da plataforma ---------------------------%

Jbi = zeros(3,6,6);
Jbii = zeros(3,6,6);
din_mat.Jp = zeros(6,6);

for i = 1:6
    Jbi(:,:,i) = [1 0 0 0 b(3,i) -b(2,i);... 
                  0 1 0 -b(3,i) 0 b(1,i);...
                  0 0 1 b(2,i) -b(1,i) 0];
    Jbii(:,:,i) = RCA(:,:,i)'*Jbi(:,:,i);
  
    din_mat.Jp(i,:) = Jbii(3,:,i);
end


% ----------------- Closed form dynamics of the limbs ------------------- %
mce = zeros(1,6);
mco = zeros(1,6);
mge = zeros(1,6);
Ieq = I1ii + I2ii;
s_hat = zeros(3,6);
sx = zeros(3,3,6);
s_hatx = zeros(3,3,6);
dL = zeros(1,6);
g = [0;0;-9.807];

Ml = zeros(6,6,6);
Cl = zeros(6,6,6);
Gl = zeros(6,1,6);

Mp =  [mp*eye(3) zeros(3); zeros(3) Ip];
Cp = [zeros(3) zeros(3); zeros(3) skws(wp)*Ip];
Gp = [-mp*g; zeros(3,1)];

Ji = zeros(3,6,6);
dJi = zeros(3,6,6);

M = zeros(3,3,6);
C = zeros(3,3,6);
G = zeros(3,1,6);

sumM = zeros(6,6);
sumC = zeros(6,6);
sumG = zeros(6,1);

for i = 1:6
    mce(i) = (1/(L(i)^2))*(m1(i)*(e1^2)+(m2(i)*(e2^2)));
    mco(i) = ((1/(L(i)))*m2(i)*e2)-(i/(L(i)^2))*(Ieq(1,1,i)+(L(i)^2)*mce(i));
    mge(i) = (1/L(i))*((m1(i)*e1)+m2(i)*(L(i)-e2));
    
    s_hat(:,i) = s(:,i)/(norm(s(:,i))); 
    sx(:,:,i) = skws(s(:,i));
    s_hatx(:,:,i) = skws(s_hat(:,i));
    
    dL(i)=s_hat(:,i)'*vp;

    M(:,:,i) = (m2(i)*s_hat(:,i)*s_hat(:,i)')-(1/(L(i)^2)*Ieq(1,1,i)*(s_hatx(:,:,i))^2)-mce(i)*((s_hatx(:,:,i))^2);
    C(:,:,i) = -((2/L(i))*mco(i)*dL(i)*(s_hatx(:,:,i)^2))-((1/(L(i)^2))*m2(i)*e2*s_hat(:,i)*vp'*(s_hatx(:,:,i)^2));
    G(:,:,i) = ((mge(i)*(s_hatx(:,:,i)^2))-(m2(i)*s_hat(:,i)*s_hat(:,i)'))*g;
    
    Ji(:,:,i) = [eye(3) -skws(b(:,i))];
    dJi(:,:,i) = [zeros(3) -skws(wp)*skws(b(:,i))+skws(b(:,i))*skws(wp)];   
    
    Ml(:,:,i) = Ji(:,:,i)'*M(:,:,i)*Ji(:,:,i);
    Cl(:,:,i) = Ji(:,:,i)'*M(:,:,i)*dJi(:,:,i) + Ji(:,:,i)'*C(:,:,i)*Ji(:,:,i);
    Gl(:,:,i) = Ji(:,:,i)'*G(:,:,i);
    
    sumM = sumM + Ml(:,:,i);
    sumC = sumC + Cl(:,:,i);
    sumG = sumG + Gl(:,:,i);
end

din_mat.MX = Mp + sumM;
din_mat.CX = Cp + sumC;
din_mat.GX = Gp + sumG;

end
