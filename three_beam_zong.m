function three_beam_zong
%三根梁，身管受到两个弹簧点支撑，身管与上下导轨之间是连续分布弹簧基础,上下导轨受电磁力作用
clear;
clc;
tic
% format long
global k_Elasticbase  k_Springsupport  m  EI  L  JieDianShu  %kz
global g TotalUnit Unitlength
k_Elasticbase(1)=8.8e3;    %N/m^2   上导轨与身管间分布弹性刚度(线载荷）
k_Elasticbase(2)=9e3;       %N/m^2   下导轨与身管间分布弹性刚度(线载荷）
%试着调换一下上限导轨和身管刚度顺序、改数
k_Springsupport=5e4;     %N/m 身管受地面支撑弹簧刚度(点载荷）
% kz=0.15E3; 考不考虑扭转刚度(实验报告中没考虑）
m(1)=50;    % 上导轨线质量 kg/m
m(2)=63.37;    %身管线质量
m(3)=50;    % 下导轨线质量
m(4)=10;    % 弹丸质量kg/m^3
EI(1)=1.5133e10*1/12*0.1^4;  %上导轨抗弯截面强度
EI(2)=2.5133e10*1/12*0.1^4;  %身管抗弯截面强度
EI(3)=1.5133e10*1/12*0.1^4;  %下导轨抗弯截面强度
L(1)=5.44; 
L(2)=3.52; 
L(3)=2.24;
L(4)=11.2;  %全长11.2m
JieDianShu(1)=1361;
JieDianShu(2)=881;   %一共2801个节点2800段 1361+881+561-2
JieDianShu(3)=561;
JieDianShu(4)=2801;  %总节点数
TotalUnit=2800;      %总单元数
Unitlength=L(4)*(1/(JieDianShu(4)-1));  %单元长度
Distance_rail=0.1;       %上下导轨间距离的一半
g=-9.8; %重力加速度，方向向下
% tspan=(0:0.001:12);
tspan=(0:0.001:3)';  %时间  /ms


omega=Cal_delta(1E-10,1,40);
% freq=omega/2/pi    %角频率与圆频率之间相互转换，Comsol求出的特征频率为圆频率，传递矩阵求得是角频率
w=omega(1:8);
[SDG_ZX,SG_ZX,XDG_ZX,~]=Cal_ZX(w);   %小数点后精度影响很大
% for j=1:length(w)      %画3-8阶振型
%     b=freq(j);
%     HuaZhenXing(SDG_ZX(:,j),SG_ZX(:,j),XDG_ZX(:,j))
%     title(['第' num2str(j) '阶 w=' num2str(b)])    %中间要有空格！
% end
%模态叠加求解静位移
% [V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,M_Mode,F_generalize]=Cal_StaticDisp(SDG_ZX,SG_ZX,XDG_ZX,w);
%验证模态质量矩阵
% Test_M_Mode(M_Mode,F_generalize,w,SDG_ZX,SG_ZX,XDG_ZX,V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,tspan)
Current=Cal_current(tspan);
F_EM=Cal_Electromagnetic_Force(Current,tspan,Distance_rail); %求电磁力
M_Mode=Cal_ModalMass(SDG_ZX,SG_ZX,XDG_ZX,w);%求模态质量
F_EM_generalize=Cal_Force_Transient(SDG_ZX,SG_ZX,XDG_ZX,w,F_EM,tspan);%求模态质量和广义力
Cal_DynamicResponse(tspan,M_Mode,F_EM_generalize,w,SDG_ZX,SG_ZX,XDG_ZX) %电磁力作用下的导轨身管动态响应

toc
end
function Cal_DynamicResponse(tspan,M_Mode,F_generalize,W,SDG_ZX,SG_ZX,XDG_ZX)
global  TotalUnit Unitlength JieDianShu
q_generalize=zeros(length(tspan),length(W));     %广义坐标   length(t)*length(w)
Step=0.001; %步长
for LL=1:length(W)
  m_mode=M_Mode(LL);
  ww=W(LL);
  f_generalize=F_generalize(:,LL);
  %自编龙格库塔
  Generalize_Coordinate(1,:)=[0 0]; %初始值
  [Q1,~]=Four_RK(tspan,Step,Generalize_Coordinate,f_generalize,m_mode,ww);
  %matlab自带龙格库塔法
%   opts = odeset('RelTol',1e-17,'AbsTol',1e-20);  %精度影响很大，默认Rel=1e-3；Abs=1e-6；
%   [t,Generalize_Coordinate]=ode45(@(t,Generalize_Coordinate) Cal_ODE(t,Generalize_Coordinate,f_generalize,m_mode,ww),[0:0.001:3],[0,0],opts); %[0 tspan(QQ)]?
  q_generalize(:,LL)=Q1(:,1);
end
%%
%模态叠加
V_SDG_Disp=SDG_ZX*q_generalize';   %JieDianShu(4)*length(w) * (length(t)*length(w))'=JieDianShu(4)*length(t)
V_SG_Disp=SG_ZX*q_generalize';     %广义坐标*模态位移 ，累加得到物理位移
V_XDG_Disp=XDG_ZX*q_generalize';
% 画图
XX=(0:TotalUnit)*Unitlength;
figure
hold on
plot(XX,V_SDG_Disp(:,length(tspan)),'--g',XX,V_SG_Disp(:,length(tspan)),'-.b',XX,V_XDG_Disp(:,length(tspan)),':r')
legend('上导轨','身管','下导轨');
title(['电磁力作用下t=' num2str(tspan(end)) 'ms时导轨身管位移'])
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('距炮尾位置/m')
ylabel('位移/mm')
hold off
figure
plot(tspan,V_SDG_Disp(JieDianShu(4),:),'--g',tspan,V_SG_Disp(JieDianShu(4),:),'-.b',tspan,V_XDG_Disp(JieDianShu(4),:),':r')
legend('上导轨','身管','下导轨');
title('电磁力作用下导轨身管炮口处位移时间历程')
set(gca,'LineWidth',3)
xlabel('t/s')
ylabel('位移/mm')
end
function [Generalize_Coordinate,n]=Four_RK(t,h,Generalize_Coordinate,f_generalize,m_mode,ww)
n=size(t,1);
for i=2:length(t)			%循环：直到求完最后一个t取值
  K1 = Cal_ODE(t(i-1),Generalize_Coordinate(i-1,:),f_generalize(i-1),m_mode,ww);  		% K1
  K2 = Cal_ODE(t(i-1)+1/2*h,Generalize_Coordinate(i-1,:)+1/2*h.*K1,f_generalize(i-1),m_mode,ww);		% K2
  K3 = Cal_ODE(t(i-1)+1/2*h,Generalize_Coordinate(i-1,:)+1/2*h.*K2,f_generalize(i-1),m_mode,ww);		% K3
  K4 = Cal_ODE(t(i-1)+h,Generalize_Coordinate(i-1,:)+h.*K3,f_generalize(i-1),m_mode,ww);		% K4
  Generalize_Coordinate(i,:) = Generalize_Coordinate(i-1,:)+h/6.*(K1 + 2*K2 + 2*K3 + K4);
end% 更新Generalize_Coordinate的第i行，Generalize_Coordinate(i,1)是原函数，Generalize_Coordinate(i,2)是原函数的一阶导
end
function dq=Cal_ODE(t,Generalize_Coordinate,f_generalize,m_mode,ww)%龙格库塔法 降阶微分方程组
%Generalize_Coordinate包含两个量，第一个是q，第二个是q′
 dq=zeros(size(Generalize_Coordinate));
 dq(1)=Generalize_Coordinate(2); %q′
 dq(2)=f_generalize*(1/m_mode)-ww^2*Generalize_Coordinate(1); %q″
end
%%
function M_Mode=Cal_ModalMass(SDG_ZX,SG_ZX,XDG_ZX,W) %求模态质量
global m JieDianShu Unitlength
M_Mode=zeros(length(W),1);  %模态质量
V_SDG_Square=SDG_ZX.^2;
V_SG_Square=SG_ZX.^2;
V_XDG_Square=XDG_ZX.^2; %JieDianShu(4)*length(W)
C_simpson=Unitlength/3;  %辛普森公式系数
for ee=1:length(W) 
  V_sdg_Square=V_SDG_Square(:,ee);
  V_sg_Square=V_SG_Square(:,ee);
  V_xdg_Square=V_XDG_Square(:,ee);
  %%
  %复化辛普森公式求模态质量
  M_Mode1=0;
  M_Mode2=0;
  M_Mode3=0;
  for hh=1:(JieDianShu(4)-1)/2  %1-1400
    M_Mode1=M_Mode1+C_simpson*m(1)*(V_sdg_Square(2*hh-1)+4*V_sdg_Square(2*hh)+V_sdg_Square(2*hh+1));
    M_Mode2=M_Mode2+C_simpson*m(2)*(V_sg_Square(2*hh-1)+4*V_sg_Square(2*hh)+V_sg_Square(2*hh+1));
    M_Mode3=M_Mode3+C_simpson*m(3)*(V_xdg_Square(2*hh-1)+4*V_xdg_Square(2*hh)+V_xdg_Square(2*hh+1));
  end
  M_Mode(ee)=M_Mode1+M_Mode2+M_Mode3;
end
end
function F_generalize=Cal_Force_Transient(SDG_ZX,SG_ZX,XDG_ZX,W,F_EM,t)%求上导轨对下导轨产生的电磁力(电枢固定在炮口，最长闭合回路)对应广义力
global  TotalUnit JieDianShu Unitlength
F_generalize=zeros(length(t),length(W));   %广义力
%求节点力矩阵
F_JieDian=zeros(length(t),JieDianShu(4));
F_JieDian(:,1)=F_EM(:,1);
F_JieDian(:,JieDianShu(4))=F_EM(:,TotalUnit);
for ii=1:TotalUnit-1
  F_JieDian(:,ii+1)=(F_EM(:,ii)+F_EM(:,ii+1))*0.5;
end
% figure
% XX=((0:TotalUnit)*Unitlength)'
% plot(XX,F_JieDian(end,:))
% title('节点力')
C_simpson=Unitlength/3;  %辛普森公式系数

%求广义力
for ee=1:length(W) 
    V_sdg=SDG_ZX(:,ee);
%     V_sg=SG_ZX(:,ee);
    V_xdg=XDG_ZX(:,ee);
    %%
    %复化辛普森公式求电磁力对应的广义力（1/3法）
    F_generalize_sdg=0;
    F_generalize_xdg=0;
    for cc=1:(JieDianShu(4)-1)/2  %1-1400
      F_generalize_sdg=F_generalize_sdg+C_simpson*(V_sdg(2*cc-1)*F_JieDian(:,2*cc-1)+4*V_sdg(2*cc)*F_JieDian(:,2*cc)+V_sdg(2*cc+1)*F_JieDian(:,2*cc+1));
      F_generalize_xdg=F_generalize_xdg+C_simpson*(V_xdg(2*cc-1)*F_JieDian(:,2*cc-1)+4*V_xdg(2*cc)*F_JieDian(:,2*cc)+V_xdg(2*cc+1)*F_JieDian(:,2*cc+1));
    end
    F_generalize(:,ee)=F_generalize_sdg-F_generalize_xdg;
%     F_generalize(:,ee)=F_generalize_sdg;
%     F_generalize(:,ee)=-F_generalize_xdg;
end
end
function F_EM=Cal_Electromagnetic_Force(Current,t,h)%求解上导轨对下导轨产生的电磁力(电枢固定在炮口，最长闭合回路)
global L Unitlength TotalUnit   %TotalUnit=2800
miu0=4*pi*1e-7;  %真空磁导率 MLT^(-2)I^(-2)
XX=((0:(TotalUnit-1))*Unitlength)';  %TotalUnit段，最后一段(11.2-Unitlength),距导轨末端(炮尾)位置
YY=L(4)-XX-Unitlength;           %TotalUnit段，最后一段0m,距炮口位置
Current_Square=Current.^2*1e6;  %I^2 I单位是KA!!!
F_EM_rr=zeros(length(t),length(XX));  %length(t)*TotalUnit
F_EM_ar=zeros(length(t),length(XX));
Constant1=miu0*(1/8)*(1/pi)*(1/h);

%每时刻导轨总的电磁力
%F_EM_rr_total=miu0/4/pi/h*Current_Square*(sqrt(4*h^2+L(4)^2)-2*h);
%F_EM_rr=F_EM_rr_total/L(4);  %单位长度导轨所受电磁力 ????

for ii=1:length(XX)  %1:TotalUnit段，每一段下导轨0-3ms内所受电磁力的变化(1-TotalUnit)
    F_EM_rr(:,ii)=Constant1*Current_Square*(sqrt(4*h^2+(Unitlength+YY(ii))^2)+sqrt(4*h^2+(Unitlength+XX(ii))^2) ...
      -sqrt(XX(ii)^2+4*h^2)-sqrt(YY(ii)^2+4*h^2));
end %F_EM_rr(:,i)第i段电磁力
% figure
% plot(t,F_EM_rr(:,1))  %0-3s下导轨最末段(炮口)所受电磁力时间历程曲线
% title('（炮尾)所受电磁力时间历程曲线');
% figure
% plot(XX,F_EM_rr(length(t),:)) %3ms下导轨受力情况
% title('下导轨受力情况')
for jj=1:length(XX)  %1-(TotalUnit-1)段，  !每一段!电枢对下导轨的电磁力，最后一段不适用
    F_EM_ar(:,jj)=miu0*(1/4)*(1/pi)*Current_Square*(log(2*h^2+h*sqrt(4*h^2+(L(4)-XX(jj)-Unitlength)^2))+log(L(4)-XX(jj)) ...
      -log(L(4)-XX(jj)-Unitlength)-log(2*h^2+h*sqrt(4*h^2+(L(4)-XX(jj))^2)));
end  %F_EM_ar(:,i)第i段电磁力 结果有问题
% figure
% plot(t,F_EM_ar(:,1))  %0-3s下导轨最末段(炮尾)所受电枢作用电磁力时间历程曲线
% title('2');
% figure
% plot(XX,F_EM_ar(length(t),:)) %3ms下导轨受力情况
% title('电枢对下导轨作用力')
% F_EM=F_EM_rr+F_EM_ar;
F_EM=F_EM_rr;
% F_EM_mean=mean(F_EM,2);
% F_EM=repmat(F_EM_mean,1,2800);
F_EM=ones(3001,2800)*353.628821955585;
% figure
% plot(t,F_EM(:,1))
% title('（炮尾)所受总电磁力时间历程曲线');
% figure
% plot(XX,F_EM(length(t),:)) %3ms下导轨受力情况
% title('下导轨总受力情况')
end
function Current=Cal_current(t) %计算输入电流  单位：千安
t11=(length(t)-1)/6+1;
t1=t(1:t11);
t22=(length(t)-1)/3;
t33=(length(t)-1)/2+2;
t3=t(t33:length(t));     
%%
% %电流时间12ms
% Current1=(300*sin(3140*t1/4000))'; %正弦上升段
% Current2=300*ones(t22,1);      %恒定段    时间单位转化在此已考虑了（/1000)
% Current3=(300*exp(-(t3-6)/0.0005/4000))'; 
%%
% %电流时间3ms
Current1=300*sin(3140*t1/1000); %正弦上升段
Current2=300*ones(t22,1);      %恒定段    时间单位转化在此已考虑了（/1000)     
Current3=300*exp(-(t3-1.5)/0.0005/1000);   %指数下降段
%%
Current=[Current1;Current2;Current3];
% figure
% plot(t,Current)
% ylim([0,400])
% title('脉冲电流时间变化曲线')
% xlabel('t/ms')
% ylabel('I/KA')
end
%%
function Test_M_Mode(M_Mode,F_generalize,W,SDG_ZX,SG_ZX,XDG_ZX,V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,t)%验证模态质量矩阵
global TotalUnit Unitlength JieDianShu
q_Mode=zeros(length(t),length(W));  %广义坐标矩阵length(t)*length(W)
 for ii=1:length(W) %求广义坐标
   q_Mode(:,ii)=F_generalize(ii)/M_Mode(ii)/W(ii)^2*(1-cos(W(ii)*t));  %q的解析解代入初始值(均为0)算的
 end
 %模态叠加
 V_SDG_StaticDisp_Test=SDG_ZX*q_Mode';
 V_SG_StaticDisp_Test=SG_ZX*q_Mode';
 V_XDG_StaticDisp_Test=XDG_ZX*q_Mode';
 
 %画图
XX=(0:TotalUnit)*Unitlength;
figure
plot(XX,V_SDG_StaticDisp_Test(:,end),'--g',XX,V_SG_StaticDisp_Test(:,end),'-.b',XX,V_XDG_StaticDisp_Test(:,end),':r');
title(['t=' num2str(t(end)) 'ms时导轨身管位移'])
legend('上导轨','身管','下导轨')
% figure
% plot(XX,V_SDG_StaticDisp,'--g',XX,V_SG_StaticDisp,'-.b',XX,V_XDG_StaticDisp,':r')
% title('模态叠加法求解重力作用下的静位移')
figure
plot(t,V_SDG_StaticDisp_Test(JieDianShu(4),:),'--g',t,V_SG_StaticDisp_Test(JieDianShu(4),:),'-.b',t,V_XDG_StaticDisp_Test(JieDianShu(4),:),':r');
legend('上导轨','身管','下导轨')
% hold on
% plot([0,3],[V_SDG_StaticDisp(JieDianShu(4)),V_SDG_StaticDisp(JieDianShu(4))],'-g')
% hold on
% plot([0,3],[V_SG_StaticDisp(JieDianShu(4)),V_SG_StaticDisp(JieDianShu(4))],'-.b')
% hold on
% plot([0,3],[V_XDG_StaticDisp(JieDianShu(4)),V_XDG_StaticDisp(JieDianShu(4))],':r')
title('炮口位移时间历程')
legend('上导轨','身管','下导轨')
xlabel('t/ms')
ylabel('位移/mm')
hold off
end
function [V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,M_Mode,F_generalize]=Cal_StaticDisp(SDG_ZX,SG_ZX,XDG_ZX,W)%用模态叠加法求动力响应(静位移）  length(W)阶叠加
global m g  Unitlength JieDianShu TotalUnit
M_Mode=zeros(length(W),1);  %模态质量
F_generalize=zeros(length(W),1);   %广义力
q_t=zeros(length(W),1);     %广义坐标
% for rr=1:length(W)   按模归一化
%     SDG_ZX(:,rr)=SDG_ZX(:,rr)/sqrt(sum(SDG_ZX(:,rr).^2));
%     SG_ZX(:,rr)=SG_ZX(:,rr)/sqrt(sum(SG_ZX(:,rr).^2));
%     XDG_ZX(:,rr)=XDG_ZX(:,rr)/sqrt(sum(XDG_ZX(:,rr).^2));
% end
V_SDG_Square=SDG_ZX.^2;
V_SG_Square=SG_ZX.^2;
V_XDG_Square=XDG_ZX.^2; %JieDianShu(4)*length(W)
C_simpson1=Unitlength/3;  %辛普森公式系数
C_simpson2=Unitlength*2/45;
C_simpson3=7*Unitlength/17280;

%求广义力、模态质量和广义坐标
for ee=1:length(W) 
    V_sdg_Square=V_SDG_Square(:,ee);
    V_sg_Square=V_SG_Square(:,ee);
    V_xdg_Square=V_XDG_Square(:,ee);
    V_sdg=SDG_ZX(:,ee);
    V_sg=SG_ZX(:,ee);
    V_xdg=XDG_ZX(:,ee);
    %%
%     %直接求模态质量和广义力，不积分
%     M_Mode1=m(1)*sum(V_sdg_Square);
%     M_Mode2=m(2)*sum(V_sg_Square);
%     M_Mode3=m(3)*sum(V_xdg_Square);
%     F_generalize1=m(1)*g*sum(V_sdg);
%     F_generalize2=m(2)*g*sum(V_sg);
%     F_generalize3=m(3)*g*sum(V_xdg);
    %%
%     %牛顿-科茨区间(n=7)
%     M_Mode1=0;F_generalize1=0;
%     M_Mode2=0;F_generalize2=0;
%     M_Mode3=0;F_generalize3=0;
%     for ff=1:(JieDianShu(4)-1)/7
%       M_Mode1=M_Mode1+C_simpson3*m(1)*(751*V_sdg_Square(7*ff-6)+751*V_sdg_Square(7*ff+1)+3577*V_sdg_Square(7*ff-5)+3577*V_sdg_Square(7*ff)+1323*V_sdg_Square(7*ff-4)+1323*V_sdg_Square(7*ff-1)+2989*V_sdg_Square(7*ff-3)+2989*V_sdg_Square(7*ff-2));
%       M_Mode2=M_Mode2+C_simpson3*m(2)*(751*V_sg_Square(7*ff-6)+751*V_sg_Square(7*ff+1)+3577*V_sg_Square(7*ff-5)+3577*V_sg_Square(7*ff)+1323*V_sg_Square(7*ff-4)+1323*V_sg_Square(7*ff-1)+2989*V_sg_Square(7*ff-3)+2989*V_sg_Square(7*ff-2));
%       M_Mode3=M_Mode3+C_simpson3*m(3)*(751*V_xdg_Square(7*ff-6)+751*V_xdg_Square(7*ff+1)+3577*V_xdg_Square(7*ff-5)+3577*V_xdg_Square(7*ff)+1323*V_xdg_Square(7*ff-4)+1323*V_xdg_Square(7*ff-1)+2989*V_xdg_Square(7*ff-3)+2989*V_xdg_Square(7*ff-2));
%       F_generalize1=F_generalize1+C_simpson3*m(1)*g*(751*V_sdg(7*ff-6)+751*V_sdg(7*ff+1)+3577*V_sdg(7*ff-5)+3577*V_sdg(7*ff)+1323*V_sdg(7*ff-4)+1323*V_sdg(7*ff-1)+2989*V_sdg(7*ff-3)+2989*V_sdg(7*ff-2));
%       F_generalize2=F_generalize2+C_simpson3*m(2)*g*(751*V_sg(7*ff-6)+751*V_sg(7*ff+1)+3577*V_sg(7*ff-5)+3577*V_sg(7*ff)+1323*V_sg(7*ff-4)+1323*V_sg(7*ff-1)+2989*V_sg(7*ff-3)+2989*V_sg(7*ff-2));
%       F_generalize3=F_generalize3+C_simpson3*m(3)*g*(751*V_xdg(7*ff-6)+751*V_xdg(7*ff+1)+3577*V_xdg(7*ff-5)+3577*V_xdg(7*ff)+1323*V_xdg(7*ff-4)+1323*V_xdg(7*ff-1)+2989*V_xdg(7*ff-3)+2989*V_xdg(7*ff-2));
%     end
    %%
%     %复华辛普森32/90法则
%     M_Mode1=0;F_generalize1=0;
%     M_Mode2=0;F_generalize2=0;
%     M_Mode3=0;F_generalize3=0;
%     for gg=1:(JieDianShu(4)-1)/4
%         M_Mode1=M_Mode1+C_simpson2*m(1)*(7*V_sdg_Square(4*gg-3)+32*V_sdg_Square(4*gg-2)+12*V_sdg_Square(4*gg-1)+32*V_sdg_Square(4*gg)+7*V_sdg_Square(4*gg+1));
%         M_Mode2=M_Mode2+C_simpson2*m(2)*(7*V_sg_Square(4*gg-3)+32*V_sg_Square(4*gg-2)+12*V_sg_Square(4*gg-1)+32*V_sg_Square(4*gg)+7*V_sg_Square(4*gg+1));
%         M_Mode3=M_Mode3+C_simpson2*m(3)*(7*V_xdg_Square(4*gg-3)+32*V_xdg_Square(4*gg-2)+12*V_xdg_Square(4*gg-1)+32*V_xdg_Square(4*gg)+7*V_xdg_Square(4*gg+1));
%         F_generalize1=F_generalize1+C_simpson2*m(1)*g*(7*V_sdg(4*gg-3)+32*V_sdg(4*gg-2)+12*V_sdg(4*gg-1)+32*V_sdg(4*gg)+7*V_sdg(4*gg+1));
%         F_generalize2=F_generalize2+C_simpson2*m(2)*g*(7*V_sg(4*gg-3)+32*V_sg(4*gg-2)+12*V_sg(4*gg-1)+32*V_sg(4*gg)+7*V_sg(4*gg+1));
%         F_generalize3=F_generalize3+C_simpson2*m(3)*g*(7*V_xdg(4*gg-3)+32*V_xdg(4*gg-2)+12*V_xdg(4*gg-1)+32*V_xdg(4*gg)+7*V_xdg(4*gg+1));
%     end
    %%
    %复化辛普森公式(1/3辛普森法则）
    M_Mode1=0;F_generalize1=0;
    M_Mode2=0;F_generalize2=0;
    M_Mode3=0;F_generalize3=0;
    for hh=1:(JieDianShu(4)-1)/2
        M_Mode1=M_Mode1+C_simpson1*m(1)*(V_sdg_Square(2*hh-1)+4*V_sdg_Square(2*hh)+V_sdg_Square(2*hh+1));
        M_Mode2=M_Mode2+C_simpson1*m(2)*(V_sg_Square(2*hh-1)+4*V_sg_Square(2*hh)+V_sg_Square(2*hh+1));
        M_Mode3=M_Mode3+C_simpson1*m(3)*(V_xdg_Square(2*hh-1)+4*V_xdg_Square(2*hh)+V_xdg_Square(2*hh+1));
        F_generalize1=F_generalize1+C_simpson1*m(1)*g*(V_sdg(2*hh-1)+4*V_sdg(2*hh)+V_sdg(2*hh+1));
        F_generalize2=F_generalize2+C_simpson1*m(2)*g*(V_sg(2*hh-1)+4*V_sg(2*hh)+V_sg(2*hh+1));
        F_generalize3=F_generalize3+C_simpson1*m(3)*g*(V_xdg(2*hh-1)+4*V_xdg(2*hh)+V_xdg(2*hh+1));
    end
    %%
%     %复化梯形公式
%     M_Mode1=0.004*m(1)*(V_sdg_Square(1)+V_sdg_Square(JieDianShu(4))+2*sum(V_sdg_Square(2:JieDianShu(4)-1)))/2;
%     M_Mode2=0.004*m(2)*(V_sg_Square(1)+V_sg_Square(JieDianShu(4))+2*sum(V_sg_Square(2:JieDianShu(4)-1)))/2;
%     M_Mode3=0.004*m(3)*(V_xdg_Square(1)+V_xdg_Square(JieDianShu(4))+2*sum(V_xdg_Square(2:JieDianShu(4)-1)))/2;
%     F_generalize1=0.004*m(1)*g*(SDG_ZX(1,ee)+SDG_ZX(JieDianShu(4),ee)+2*sum(SDG_ZX(2:JieDianShu(4)-1,ee)))/2;
%     F_generalize2=0.004*m(2)*g*(SG_ZX(1,ee)+SG_ZX(JieDianShu(4),ee)+2*sum(SG_ZX(2:JieDianShu(4)-1,ee)))/2;
%     F_generalize3=0.004*m(3)*(g)*(XDG_ZX(1,ee)+XDG_ZX(JieDianShu(4),ee)+2*sum(XDG_ZX(2:JieDianShu(4)-1,ee)))/2;
   %%
   M_Mode(ee)=M_Mode1+M_Mode2+M_Mode3;
%    F_generalize(ee)=-F_generalize1+F_generalize3; %上导轨向上重力，下导轨向下重力
%    F_generalize(ee)=F_generalize2; %只加身管重力
%    F_generalize(ee)=-F_generalize1; %只加上导轨负重力
%    F_generalize(ee)=F_generalize3; %只加下导轨重力
   F_generalize(ee)=F_generalize1+F_generalize2+F_generalize3;
   q_t(ee)=F_generalize(ee)/M_Mode(ee)/W(ee)^2;
end

% % %按模态质量归一化  V'MV=1
% for rr=1:length(W)
%     SDG_ZX(:,rr)=SDG_ZX(:,rr)/sqrt(sum(SDG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
%     SG_ZX(:,rr)=SG_ZX(:,rr)/sqrt(sum(SG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
%     XDG_ZX(:,rr)=XDG_ZX(:,rr)/sqrt(sum(XDG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
% end

% 模态叠加
V_SDG_StaticDisp=SDG_ZX*q_t;   %JieDianShu(4)*length(w) * length(w)*1
V_SG_StaticDisp=SG_ZX*q_t;
V_XDG_StaticDisp=XDG_ZX*q_t;

% 画图
XX=(0:TotalUnit)*Unitlength;
figure
hold on
plot(XX,V_SDG_StaticDisp,'--g',XX,V_SG_StaticDisp,'-.b',XX,V_XDG_StaticDisp,':r')
%%%% legend('上导轨静位移','身管静位移','下导轨静位移');
%%%%  plot(,SDG_Coordinate,sdg_Height,'y:',SG_Coordinates,sg_Height,'-.',XDG_Coordinates,xdg_Height,'-')
% legend('上导轨静位移','身管静位移','下导轨静位移')%,'comsol上导轨静位移','comsol身管静位移','comsol下导轨静位移');
title('模态叠加法求解重力作用下的静位移')
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('距炮尾位置/m')
ylabel('静位移/mm')
hold off
end
%%
function HuaZhenXing(V1,V2,V3)  %画阵型 （Y)
global TotalUnit Unitlength
X=(0:TotalUnit)*Unitlength;
figure
hold on
plot(X,V1,'--g',X,V2,'-.b',X,V3,':r')
legend('上导轨振型','身管振型','下导轨振型');
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('X/m')
ylabel('挠度/mm')
hold off
end
function [SDG_ZX,SG_ZX,XDG_ZX,Z1_I_all]=Cal_ZX(W) %求振型 (Y)
global  JieDianShu Unitlength
Z1_I_all=zeros(12,length(W));
SDG_ZX=zeros(JieDianShu(4),length(W));   
SG_ZX=zeros(JieDianShu(4),length(W));
XDG_ZX=zeros(JieDianShu(4),length(W));  %前length(W)阶振型，JieDianShu(4)=2801个节点 1361+881+561-2=2801， TotalUnit段
B=zeros((JieDianShu(1)-1)*12,12);  %存储传递矩阵，节约计算量
for aa=1:length(W)  %第aa阶固有振型
    w=W(aa);
    [U_all,U1,U2,U3,U4,~,A_exp]=Cal_U(w);
    U_guiyi=[U_all([3,4,7,8,11,12],[1,2,5,6,9,10]);
                   1 1 1 1 1 1 ];   %7*6
    %[1 2 5 6 9 10]:y1,θ1，y2,θ2,y3,θ3
    %[U_bar;[1...1]]代替奇异值分解来求Z_bar，因为奇异值分解得到的奇异值矩阵S是满秩的
    Z_bar=(U_guiyi'*U_guiyi)\(U_guiyi')*[0 0 0 0 0 0 1]';% (6*7)*(7*6)*（6*7）*（7*1）
    %在求模态质量矩阵时归一化即可
    Z1_I=[Z_bar(1) Z_bar(2) 0 0 Z_bar(3) Z_bar(4) 0 0 Z_bar(5) Z_bar(6) 0 0]';
  
    Z1_I_all(:,aa)=Z1_I;
    SDG_ZX(1,aa)=Z1_I(1);
    SG_ZX(1,aa)=Z1_I(5);
    XDG_ZX(1,aa)=Z1_I(9);
    
%   向量方法（按照全部的Y全部的θ全部的M全部的Q来排列，用的时候把每一列摘出来）
%   X1=((1:(JieDianShu(1)-1))*0.08)';   %单元长度0.08   67*1  最后一个X1的传递矩阵就是U1
%   U1_X1=Cal_U_elasticbeam(m,EI,w,k_Elasticbase,X1);       %(67*12)*12
%   第一段传递矩阵  A是12*12的，怎么把X1（67*1）代入相乘????
%   %Z1_X1=U1_X1(1:67,1:4)*Z1_I;      %  67*1  第一段输出端状态矢量
%   SG_ZX(2:68,i)=U1_X1(1:67,1:4)*Z1_I;     %第一段  节点号（2-68)

%   循环
    for jj=2:JieDianShu(1)  %2-1361 JieDianShu(1)=1361
        X1=Unitlength*(jj-1);    %单元长度0.004
        U1_X1=expm(A_exp*X1);    %  最后一个X1的传递矩阵就是U1
        B(12*(jj-2)+1:12*(jj-1),1:12)=U1_X1;  %存储传递矩阵
        SDG_ZX(jj,aa)=U1_X1(1,:)*Z1_I;
        SG_ZX(jj,aa)=U1_X1(5,:)*Z1_I;   %2-1361
        XDG_ZX(jj,aa)=U1_X1(9,:)*Z1_I;
    end   %循环结束得到三梁第一段从第2个点到第1361个点的纵向位移值
     
    Z3_I=U2*U1*Z1_I;  %减少运算; 12*1
    for dd=1:(JieDianShu(2)-1)    %JieDianShu(2)=881
        U3_X2=B(12*(dd-1)+1:12*dd,1:12);   %提取传递矩阵  最后一个U3_X2的传递矩阵就是U3
        SDG_ZX(JieDianShu(1)+dd,aa)=U3_X2(1,:)*Z3_I;
        SG_ZX(JieDianShu(1)+dd,aa)=U3_X2(5,:)*Z3_I;
        XDG_ZX(JieDianShu(1)+dd,aa)=U3_X2(9,:)*Z3_I;  %1362-2241 
    end %循环结束得到三梁第二段从第1362个点到第2241个点的纵向位移值
    
    Z5_I=U4*U3*Z3_I;
    for bb=2:JieDianShu(3)      %JieDianShu(3)=561
        U5_X3=B(12*(bb-2)+1:12*(bb-1),1:12); %提取传递矩阵  最后一个U5_X3的传递矩阵就是U5
        SDG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(1,:)*Z5_I;
        SG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(5,:)*Z5_I;
        XDG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(9,:)*Z5_I;  %2241-2801
    end  %循环结束得到三梁第三段从第2241个点到第2801个点的纵向位移值
end
end  
function w=Cal_delta(Pre,Step,n)   %计算固有频率(精度，步长，阶数）  (Y)
        w=zeros(n,1);  %提前给n阶固有频率分配内存地址 
        a=0 ;    %算到第几阶固有频率
        x1=0+eps;
        x2=x1+Step;
        fun1=Cal(x1);
        fun2=Cal(x2);
        while a<n
            while fun1*fun2>0
                  x1=x2;
                  fun1=fun2;
                  x2=x2+Step;
                  fun2=Cal(x2);
            end
            while abs(x2-x1)>Pre
          %fun1*fun2<0   判断条件没有必要了，当不满足第一个判断条件，自然就是<0
                  x0=(x1+x2)/2;
                  fun0=Cal(x0);
                  if fun0*fun1>0
                     x1=x0;
                     fun1=fun0;
                  else
                     x2=x0;
                     fun2=fun0;
                  end
         end   %此处循环结束时的状态为abs(x2-x1)<Pre;fun1*fun2<0
         a=a+1;
         w(a)=(x1+x2)/2;  %得到第a阶固有频率
         x1=x2;fun1=fun2;
         x2=x1+Step;fun2=Cal(x2);   %在此基础上继续往下算，找下一阶固有频率
        end
end  
function delta=Cal(w)  %分块矩阵(Y)
[U_all,~]=Cal_U(w);
U_bar=U_all([3,4,7,8,11,12],[1,2,5,6,9,10]); %6*6 代入边界条件，两端自由
delta=det(U_bar);
end
function [U_all,U1,U2,U3,U4,U5,A_exp]=Cal_U(w)  %拼总传递矩阵(Y)
global EI m k_Elasticbase k_Springsupport L
A_exp=Pin_exp(m,EI,w,k_Elasticbase(1),k_Elasticbase(2));
U1=expm(A_exp*L(1));
U2=Cal_Springsupport(k_Springsupport);
U3=expm(A_exp*L(2));
U4=Cal_Springsupport(k_Springsupport);
U5=expm(A_exp*L(3));      
U_all=U5*U4*U3*U2*U1; %总传递矩阵
end  
function TM_Springsupport=Cal_Springsupport(k)  %虚拟元件联接弹簧支撑的传递矩阵  (Y)
TM_Springsupport=eye(12,12);
TM_Springsupport(8,5)=-k;
end
function A_exp=Pin_exp(m,EI,w,k1,k2)         %微分方程组的系数矩阵A  (Y)
%k1 上导轨与身管间分布弹性刚度    k2  下导轨与身管间分布弹性刚度
      A11=[ 0          1    0       0;
            0          0 1/EI(1)    0;
            0          0    0       1;
      (m(1)*w^2-k1)     0    0       0 ];
       A12=[zeros(3,4);
           k1 0 0 0];
       A13=zeros(4,4);
       A21=A12;
       A22=[ 0            1     0    0;
             0            0  1/EI(2) 0;
             0            0     0    1;
     m(2)*w^2-k1+k2       0     0    0 ];
       A23=[zeros(3,4);
           -k2 0 0 0];
       A31=A13;
       A32=A23;
       A33=[ 0          1     0     0;
             0          0  1/EI(3)  0;
             0          0     0     1 ;
       (m(3)*w^2+k2)     0     0     0 ];
     A_exp=[A11 A12 A13;
            A21 A22 A23;
            A31 A32 A33];
end