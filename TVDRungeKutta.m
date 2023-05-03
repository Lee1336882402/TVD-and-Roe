clear
clc
total_t = 0.001;
CFL = 0.03;
Mx = 201; %domain [0,1]
delta_x = 1/(Mx-1);
delta_t = CFL*delta_x;
total_steps = int32(total_t/delta_t);
gama = 1.4;
rup = zeros(3,Mx);
rup(1,1:int32(Mx/2)) = 1;
rup(2,1:int32(Mx/2)) = 0.75;
rup(3,1:int32(Mx/2)) = 1;
rup(1,int32(Mx/2)+1:end) = 0.125;
rup(2,int32(Mx/2)+1:end) = 0;
rup(3,int32(Mx/2)+1:end) = 0.1;
U = zeros(3,Mx);
U(1,:)=rup(1,:);
U(2,:)=rup(1,:).*rup(2,:);
U(3,:)=rup(3,:)/(gama-1)+rup(1,:).*rup(2,:).*rup(2,:)/2;
F = zeros(3,Mx);
check_point = floor(total_t/5/delta_t);
rou=rup(1,:);
u = rup(2,:);
p = rup(3,:);
data = U;
for step = 1:total_steps
    temp_data = data - delta_t/delta_x*TVD(data,1);
    temp_data = 3/4*data + 1/4*(temp_data - delta_t/delta_x*TVD(temp_data,1));
    temp_data = 1/3*data + 2/3*(temp_data - delta_t/delta_x*TVD(temp_data,1));        
    data=temp_data;
    if mod(step,check_point)==0
        rup(1,:) = data(1,:);
        rup(2,:) = data(2,:)./data(1,:);
        rup(3,:) = (data(3,:)-rup(1,:).*rup(2,:).*rup(2,:)/2)*(gama-1);
        rou = [rou;rup(1,:)]; %#ok<*AGROW>
        u = [u;rup(2,:)];
        p = [p;rup(3,:)];
    end
end

plot_x = linspace(-0.5,0.5,Mx);
figure(1);
subplot(1,3,1)
plot(plot_x,rou);
legend('t=0','t=0.06','t=0.12','t=0.18','t=0.24','t=0.30');
xlabel('x');
ylabel('rou');
ylim([0,1])
subplot(1,3,2)
plot(plot_x,u);
legend('t=0','t=0.06','t=0.12','t=0.18','t=0.24','t=0.30');
xlabel('x');
ylabel('u');
subplot(1,3,3)
plot(plot_x,p);
legend('t=0','t=0.06','t=0.12','t=0.18','t=0.24','t=0.30');
xlabel('x');
ylabel('p');



