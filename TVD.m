function result = TVD(inarray,jingdu)
%TVD 此处显示有关此函数的摘要
%   此处显示详细说明
gama = 1.4;
delta_x = 1/(length(inarray)-1);
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
rout = sqrt(rouL.*rouR);
ut = (sqrt(rouL).*uL+sqrt(rouR).*uR)./(sqrt(rouL)+sqrt(rouR));
Ht = (sqrt(rouL).*HL+sqrt(rouR).*HR)./(sqrt(rouL)+sqrt(rouR));
pt = (gama-1)/gama*rout.*(Ht-0.5*ut.^2);
at = sqrt(gama*pt./rout);
U = inarray;
l1 = 0.5*[(gama-1)*ut.^2./(2*at.^2)+ut./at;-((gama-1)*ut./(at.^2)+ones(1,length(inarray))./at);(gama-1)*ones(1,length(inarray))./(at.^2)];
l2 = 0.5*[2*-(gama-1)*ut.^2./(at.^2);2*(gama-1)*ut./(at.^2);-2*(gama-1)*ones(1,length(inarray))./(at.^2)];
l3 = 0.5*[(gama-1)*ut.^2./(2*at.^2)-ut./at;-((gama-1)*ut./(at.^2)-ones(1,length(inarray))./at);(gama-1)*ones(1,length(inarray))./(at.^2)];
r1 = [ones(1,length(inarray));ut-at;Ht-ut.*at];
r2 = [ones(1,length(inarray));ut;0.5*ut.^2];
r3 = [ones(1,length(inarray));ut+at;Ht+ut.*at];
mL1 = sum(l1.*([U(:,2:end),U(:,end)]-U));
mR1 = sum(l1.*(U-[U(:,1),U(:,1:end-1)]));
mM1 = sum(l1.*([U(:,3:end),U(:,end),U(:,end)]-[U(:,2:end),U(:,end)]));
mL2 = sum(l2.*([U(:,2:end),U(:,end)]-U));
mR2 = sum(l2.*(U-[U(:,1),U(:,1:end-1)]));
mM2 = sum(l2.*([U(:,3:end),U(:,end),U(:,end)]-[U(:,2:end),U(:,end)]));
mL3 = sum(l3.*([U(:,2:end),U(:,end)]-U));
mR3 = sum(l3.*(U-[U(:,1),U(:,1:end-1)]));
mM3 = sum(l3.*([U(:,3:end),U(:,end),U(:,end)]-[U(:,2:end),U(:,end)]));
WL = zeros(3,length(inarray));
WR = zeros(3,length(inarray));
U = [U,U(:,end)];
for index = 1:length(inarray)
    if mL1(index)*mR1(index)<=0
        WL(1,index)=l1(:,index)'*U(:,index);
    else 
        D = sign(mL1(index))*min(abs([mL1(index),mR1(index)]));
        WL(1,index)=l1(:,index)'*U(:,index)+0.5*D*delta_x;
    end
    if mM1(index)*mL1(index)<=0
        WR(1,index) = l1(:,index)'*U(:,index+1);
    else
        D = sign(mL1(index))*min(abs([mL1(index),mM1(index)]));
        WR(1,index) = l1(:,index)'*U(:,index+1)-0.5*D*delta_x;
    end
    if mL2(index)*mR2(index)<=0
        WL(2,index)=l2(:,index)'*U(:,index);
    else 
        D = sign(mL2(index))*min(abs([mL2(index),mR2(index)]));
        WL(2,index)=l2(:,index)'*U(:,index)+0.5*D*delta_x;
    end
    if mM2(index)*mL2(index)<=0
        WR(2,index) = l2(:,index)'*U(:,index+1);
    else
        D = sign(mL2(index))*min(abs([mL2(index),mM2(index)]));
        WR(2,index) = l2(:,index)'*U(:,index+1)-0.5*D*delta_x;
    end
    if mL3(index)*mR3(index)<=0
        WL(3,index)=l3(:,index)'*U(:,index);
    else 
        D = sign(mL3(index))*min(abs([mL3(index),mR3(index)]));
        WL(3,index)=l3(:,index)'*U(:,index)+0.5*D*delta_x;
    end
    if mM3(index)*mL3(index)<=0
        WR(3,index) = l3(:,index)'*U(:,index+1);
    else
        D = sign(mL3(index))*min(abs([mL3(index),mM3(index)]));
        WR(3,index) = l3(:,index)'*U(:,index+1)-0.5*D*delta_x;
    end
    
    
end
for index = 1:length(inarray)
    UL(:,index) = [r1(:,index),r2(:,index),r3(:,index)]*WL(:,index);
    UR(:,index) = [r1(:,index),r2(:,index),r3(:,index)]*WR(:,index);
end
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

