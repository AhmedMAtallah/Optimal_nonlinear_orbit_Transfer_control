function ProblemLBack2()
clear
clc
n=6


global mu k1 k2 rd
r=[-0.70545852988580, -0.73885031681775, -0.40116299069586]'*6378140;
v=[0.73122658145185, -0.53921753373056, -0.29277123328399]'* 6378140/806.812;
r=[7.0000e+06, 0, 0]';
v=[0, 7.5460e+03, 0]';

rd=[0.0;0.0];
x0=[r;v]
span=90*60*50


mu = 3.986e14
options = odeset('RelTol',1e-7);
[Tk2e01,Yk2e01]=ode45(@closelooporbit2,[0 span],x0,options) ;
plot3(Yk2e01(:,1),Yk2e01(:,2),Yk2e01(:,3))
xlabel('x')
ylabel('y')
zlabel('z')
k1 =1; k2 =10; mu = 1;
options = odeset('RelTol',1e-10);
[Tk2e10,Yk2e10]=ode45(@closelooporbit,[0 span],x0,options) ;
for i =1 :size(Tk2e10)
    r= Yk2e10(i,1:2)';    z = Yk2e10(i,3:4)'+k1*(r-[0.8 ;0.2]);    uk2e10(:,i)=-2*((1+1/k1)+mu^2/norm(r)^6/k1^3+k2)*z+mu/norm(r)^3*[0.8;0.2];
    Rk2e10(:,i) =(((1+1/k1)+mu^2/norm(r)^6/k1^3+k2))^(-1) ;
end
figure(10)
set(gca,'FontSize',18)
plot(Yk2e1(:,1),Yk2e1(:,2),'LineWidth',2)
figure(3)
set(gca,'FontSize',18)
plot(Tk2e01,Yk2e01(:,3),Tk2e1,Yk2e1(:,3),'--',Tk2e10,Yk2e10(:,3),':','LineWidth',2)
grid on
xlabel('Time (sec)','FontSize',18)
ylabel('v_1(m/s)','FontSize',18)
lgd= legend('k_2=0.1','k_2=1','k_2=10')
lgd.FontSize = 14;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16')
saveas(gca,'figure3.png')


figure(4)
set(gca,'FontSize',20)
plot(Tk2e01,Yk2e01(:,2),Tk2e1,Yk2e1(:,2),'--',Tk2e10,Yk2e10(:,2),':','LineWidth',2)
grid on 
xlabel('Time (sec)')
ylabel('r_2(m)')
legend('k_2=0.1','k_2=1','k_2=10')
saveas(gca,'figure4.png')

figure(5)
plot(Tk2e01,uk2e01(1,:),Tk2e1,uk2e1(1,:),'--',Tk2e10,uk2e10(1,:),':','LineWidth',2)
grid on
xlabel('Time (sec)','FontSize',18)
ylabel('u_1 (m/s^2)','FontSize',18)
leg=legend('k_2=0.1','k_2=1','k_2=10')
lgd.FontSize = 14;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16')
saveas(gca,'figure5.png')

figure(6)

plot(Tk2e01,Rk2e01(1,:),Tk2e1,Rk2e1(1,:),'--',Tk2e10,Rk2e10(1,:),':','LineWidth',2)
grid on
%ax = gca;
ax.FontSize = 16; 
xlabel('My Label','FontSize',20)
%xlabel('Time (sec)')
ylabel(' ||R(x)||')
legend('k_2=0.1','k_2=1','k_2=10')

saveas(gca,'figure6.png')

% subplot(311)
% plot(T,Y(:,1),T,Y(:,2))
% xlabel('t')
% ylabel('Position')
% legend('r_x','r_y')
% %figure
% subplot(312)
% plot(T,Y(:,3),T,Y(:,4))
% xlabel('t')
% ylabel('Velocity')
% legend('v_x','v_y')
% %figure
% %plot(Y(:,1),Y(:,2))
% %axis equal
% 
% %k1 =1
% %k2= 1
% subplot(313)
% 
% %figure
% subplot(313)
% plot(T,u)
% xlabel('t')
% ylabel('Inputs')
% legend('u_x','u_y')
% saveas(gca,'k1e1k2e1.png')
% a=5