function result = Roe(inarray,jingdu)
%ROE 此处显示有关此函数的摘要
%   此处显示详细说明
gama = 1.4;
FL = zeros(3,length(inarray));
FR = zeros(3,length(inarray));
if jingdu == 1
    UL = inarray;
    UR = [inarray(:,2:end),inarray(:,end)];
else 
    if jingdu == 2
        UL = inarray+0.25*([inarray(:,2:end),inarray(:,end)]-[inarray(:,1),inarray(:,1:end-1)]);
        UR = [inarray(:,2:end),inarray(:,end)] - 0.25*([inarray(:,3:end),inarray(:,end),inarray(:,end)]-inarray);
    end
end
% U = inarray;
% rup = zeros(1,length(inarray));
% rup(1,:) = U(1,:);
% rup(2,:) = U(2,:)./U(1,:);
% rup(3,:) = (U(3,:)-rup(1,:).*rup(2,:).*rup(2,:)/2)*(gama-1);
rouL = UL(1,:);
uL = UL(2,:)./UL(1,:);
pL = (UL(3,:)-rouL.*uL.*uL/2)*(gama-1);
rouR = UR(1,:);
uR = UR(2,:)./UR(1,:);
pR = (UR(3,:)-rouR.*uR.*uR/2)*(gama-1);
HL = gama*pL./rouL/(gama-1)+0.5*uL.^2;
HR = gama*pR./rouR/(gama-1)+0.5*uR.^2;
FL(1,:)=UL(2,:);
FL(2,:)=UL(2,:).*UL(2,:)./UL(1,:)+(gama-1)*(UL(3,:)-UL(2,:).*UL(2,:)./UL(1,:)/2);
FL(3,:)=UL(2,:)./UL(1,:).*(UL(3,:)+(gama-1)*(UL(3,:)-UL(2,:).*UL(2,:)./UL(1,:)/2));
FR(1,:)=UR(2,:);
FR(2,:)=UR(2,:).*UR(2,:)./UR(1,:)+(gama-1)*(UR(3,:)-UR(2,:).*UR(2,:)./UR(1,:)/2);
FR(3,:)=UR(2,:)./UR(1,:).*(UR(3,:)+(gama-1)*(UR(3,:)-UR(2,:).*UR(2,:)./UR(1,:)/2));
rous = sqrt(rouL.*rouR);
us = (sqrt(rouL).*uL+sqrt(rouR).*uR)./(sqrt(rouL)+sqrt(rouR));
Hs = (sqrt(rouL).*HL+sqrt(rouR).*HR)./(sqrt(rouL)+sqrt(rouR));
ps = (gama-1)/gama*rous.*(Hs-0.5*us.^2);
as = sqrt(gama*ps./rous);
rs1 = [ones(1,length(inarray));us-as;Hs-us.*as];
rs2 = [ones(1,length(inarray));us;0.5*us.^2];
rs3 = [ones(1,length(inarray));us+as;Hs+us.*as];
alpha1 = (pR-pL-rous.*as.*(uR-uL))./(2*as.^2);
alpha2 = rouR-rouL - (pR-pL)./(as.^2);
alpha3 = (pR-pL+rous.*as.*(uR-uL))./(2*as.^2);
abslambda1 = abs(us-as);
abslambda2 = abs(us);
abslambda3 = abs(us+as);
epsilon = 1e-1;
for index = 1:length(inarray)
    if abslambda1(index)<epsilon
        abslambda1(index) = (abslambda1(index)^2+epsilon^2)/2/epsilon;
    end
    if abslambda2(index)<epsilon
        abslambda2(index) = (abslambda2(index)^2+epsilon^2)/2/epsilon;
    end
    if abslambda3(index)<epsilon
        abslambda3(index) = (abslambda3(index)^2+epsilon^2)/2/epsilon;
    end
end
F = 0.5*(FL+FR)-0.5*((rs1.*repmat(abslambda1,3,1)).*alpha1+(rs2.*repmat(abslambda2,3,1)).*alpha2+(rs3.*repmat(abslambda3,3,1)).*alpha3);
result = F - [F(:,1),F(:,1:end-1)];
end

