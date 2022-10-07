function three_beam_zong
%������������ܵ��������ɵ�֧�ţ���������µ���֮���������ֲ����ɻ���,���µ����ܵ��������
clear;
clc;
tic
% format long
global k_Elasticbase  k_Springsupport  m  EI  L  JieDianShu  %kz
global g TotalUnit Unitlength
k_Elasticbase(1)=8.8e3;    %N/m^2   �ϵ�������ܼ�ֲ����Ըն�(���غɣ�
k_Elasticbase(2)=9e3;       %N/m^2   �µ�������ܼ�ֲ����Ըն�(���غɣ�
%���ŵ���һ�����޵������ܸն�˳�򡢸���
k_Springsupport=5e4;     %N/m ����ܵ���֧�ŵ��ɸն�(���غɣ�
% kz=0.15E3; ��������Ťת�ն�(ʵ�鱨����û���ǣ�
m(1)=50;    % �ϵ��������� kg/m
m(2)=63.37;    %���������
m(3)=50;    % �µ���������
m(4)=10;    % ��������kg/m^3
EI(1)=1.5133e10*1/12*0.1^4;  %�ϵ��쿹�����ǿ��
EI(2)=2.5133e10*1/12*0.1^4;  %��ܿ������ǿ��
EI(3)=1.5133e10*1/12*0.1^4;  %�µ��쿹�����ǿ��
L(1)=5.44; 
L(2)=3.52; 
L(3)=2.24;
L(4)=11.2;  %ȫ��11.2m
JieDianShu(1)=1361;
JieDianShu(2)=881;   %һ��2801���ڵ�2800�� 1361+881+561-2
JieDianShu(3)=561;
JieDianShu(4)=2801;  %�ܽڵ���
TotalUnit=2800;      %�ܵ�Ԫ��
Unitlength=L(4)*(1/(JieDianShu(4)-1));  %��Ԫ����
Distance_rail=0.1;       %���µ��������һ��
g=-9.8; %�������ٶȣ���������
% tspan=(0:0.001:12);
tspan=(0:0.001:3)';  %ʱ��  /ms


omega=Cal_delta(1E-10,1,40);
% freq=omega/2/pi    %��Ƶ����ԲƵ��֮���໥ת����Comsol���������Ƶ��ΪԲƵ�ʣ����ݾ�������ǽ�Ƶ��
w=omega(1:8);
[SDG_ZX,SG_ZX,XDG_ZX,~]=Cal_ZX(w);   %С����󾫶�Ӱ��ܴ�
% for j=1:length(w)      %��3-8������
%     b=freq(j);
%     HuaZhenXing(SDG_ZX(:,j),SG_ZX(:,j),XDG_ZX(:,j))
%     title(['��' num2str(j) '�� w=' num2str(b)])    %�м�Ҫ�пո�
% end
%ģ̬������⾲λ��
% [V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,M_Mode,F_generalize]=Cal_StaticDisp(SDG_ZX,SG_ZX,XDG_ZX,w);
%��֤ģ̬��������
% Test_M_Mode(M_Mode,F_generalize,w,SDG_ZX,SG_ZX,XDG_ZX,V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,tspan)
Current=Cal_current(tspan);
F_EM=Cal_Electromagnetic_Force(Current,tspan,Distance_rail); %������
M_Mode=Cal_ModalMass(SDG_ZX,SG_ZX,XDG_ZX,w);%��ģ̬����
F_EM_generalize=Cal_Force_Transient(SDG_ZX,SG_ZX,XDG_ZX,w,F_EM,tspan);%��ģ̬�����͹�����
Cal_DynamicResponse(tspan,M_Mode,F_EM_generalize,w,SDG_ZX,SG_ZX,XDG_ZX) %����������µĵ�����ܶ�̬��Ӧ

toc
end
function Cal_DynamicResponse(tspan,M_Mode,F_generalize,W,SDG_ZX,SG_ZX,XDG_ZX)
global  TotalUnit Unitlength JieDianShu
q_generalize=zeros(length(tspan),length(W));     %��������   length(t)*length(w)
Step=0.001; %����
for LL=1:length(W)
  m_mode=M_Mode(LL);
  ww=W(LL);
  f_generalize=F_generalize(:,LL);
  %�Ա��������
  Generalize_Coordinate(1,:)=[0 0]; %��ʼֵ
  [Q1,~]=Four_RK(tspan,Step,Generalize_Coordinate,f_generalize,m_mode,ww);
  %matlab�Դ����������
%   opts = odeset('RelTol',1e-17,'AbsTol',1e-20);  %����Ӱ��ܴ�Ĭ��Rel=1e-3��Abs=1e-6��
%   [t,Generalize_Coordinate]=ode45(@(t,Generalize_Coordinate) Cal_ODE(t,Generalize_Coordinate,f_generalize,m_mode,ww),[0:0.001:3],[0,0],opts); %[0 tspan(QQ)]?
  q_generalize(:,LL)=Q1(:,1);
end
%%
%ģ̬����
V_SDG_Disp=SDG_ZX*q_generalize';   %JieDianShu(4)*length(w) * (length(t)*length(w))'=JieDianShu(4)*length(t)
V_SG_Disp=SG_ZX*q_generalize';     %��������*ģ̬λ�� ���ۼӵõ�����λ��
V_XDG_Disp=XDG_ZX*q_generalize';
% ��ͼ
XX=(0:TotalUnit)*Unitlength;
figure
hold on
plot(XX,V_SDG_Disp(:,length(tspan)),'--g',XX,V_SG_Disp(:,length(tspan)),'-.b',XX,V_XDG_Disp(:,length(tspan)),':r')
legend('�ϵ���','���','�µ���');
title(['�����������t=' num2str(tspan(end)) 'msʱ�������λ��'])
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('����βλ��/m')
ylabel('λ��/mm')
hold off
figure
plot(tspan,V_SDG_Disp(JieDianShu(4),:),'--g',tspan,V_SG_Disp(JieDianShu(4),:),'-.b',tspan,V_XDG_Disp(JieDianShu(4),:),':r')
legend('�ϵ���','���','�µ���');
title('����������µ�������ڿڴ�λ��ʱ������')
set(gca,'LineWidth',3)
xlabel('t/s')
ylabel('λ��/mm')
end
function [Generalize_Coordinate,n]=Four_RK(t,h,Generalize_Coordinate,f_generalize,m_mode,ww)
n=size(t,1);
for i=2:length(t)			%ѭ����ֱ���������һ��tȡֵ
  K1 = Cal_ODE(t(i-1),Generalize_Coordinate(i-1,:),f_generalize(i-1),m_mode,ww);  		% K1
  K2 = Cal_ODE(t(i-1)+1/2*h,Generalize_Coordinate(i-1,:)+1/2*h.*K1,f_generalize(i-1),m_mode,ww);		% K2
  K3 = Cal_ODE(t(i-1)+1/2*h,Generalize_Coordinate(i-1,:)+1/2*h.*K2,f_generalize(i-1),m_mode,ww);		% K3
  K4 = Cal_ODE(t(i-1)+h,Generalize_Coordinate(i-1,:)+h.*K3,f_generalize(i-1),m_mode,ww);		% K4
  Generalize_Coordinate(i,:) = Generalize_Coordinate(i-1,:)+h/6.*(K1 + 2*K2 + 2*K3 + K4);
end% ����Generalize_Coordinate�ĵ�i�У�Generalize_Coordinate(i,1)��ԭ������Generalize_Coordinate(i,2)��ԭ������һ�׵�
end
function dq=Cal_ODE(t,Generalize_Coordinate,f_generalize,m_mode,ww)%��������� ����΢�ַ�����
%Generalize_Coordinate��������������һ����q���ڶ�����q��
 dq=zeros(size(Generalize_Coordinate));
 dq(1)=Generalize_Coordinate(2); %q��
 dq(2)=f_generalize*(1/m_mode)-ww^2*Generalize_Coordinate(1); %q��
end
%%
function M_Mode=Cal_ModalMass(SDG_ZX,SG_ZX,XDG_ZX,W) %��ģ̬����
global m JieDianShu Unitlength
M_Mode=zeros(length(W),1);  %ģ̬����
V_SDG_Square=SDG_ZX.^2;
V_SG_Square=SG_ZX.^2;
V_XDG_Square=XDG_ZX.^2; %JieDianShu(4)*length(W)
C_simpson=Unitlength/3;  %����ɭ��ʽϵ��
for ee=1:length(W) 
  V_sdg_Square=V_SDG_Square(:,ee);
  V_sg_Square=V_SG_Square(:,ee);
  V_xdg_Square=V_XDG_Square(:,ee);
  %%
  %��������ɭ��ʽ��ģ̬����
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
function F_generalize=Cal_Force_Transient(SDG_ZX,SG_ZX,XDG_ZX,W,F_EM,t)%���ϵ�����µ�������ĵ����(����̶����ڿڣ���պϻ�·)��Ӧ������
global  TotalUnit JieDianShu Unitlength
F_generalize=zeros(length(t),length(W));   %������
%��ڵ�������
F_JieDian=zeros(length(t),JieDianShu(4));
F_JieDian(:,1)=F_EM(:,1);
F_JieDian(:,JieDianShu(4))=F_EM(:,TotalUnit);
for ii=1:TotalUnit-1
  F_JieDian(:,ii+1)=(F_EM(:,ii)+F_EM(:,ii+1))*0.5;
end
% figure
% XX=((0:TotalUnit)*Unitlength)'
% plot(XX,F_JieDian(end,:))
% title('�ڵ���')
C_simpson=Unitlength/3;  %����ɭ��ʽϵ��

%�������
for ee=1:length(W) 
    V_sdg=SDG_ZX(:,ee);
%     V_sg=SG_ZX(:,ee);
    V_xdg=XDG_ZX(:,ee);
    %%
    %��������ɭ��ʽ��������Ӧ�Ĺ�������1/3����
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
function F_EM=Cal_Electromagnetic_Force(Current,t,h)%����ϵ�����µ�������ĵ����(����̶����ڿڣ���պϻ�·)
global L Unitlength TotalUnit   %TotalUnit=2800
miu0=4*pi*1e-7;  %��մŵ��� MLT^(-2)I^(-2)
XX=((0:(TotalUnit-1))*Unitlength)';  %TotalUnit�Σ����һ��(11.2-Unitlength),�ർ��ĩ��(��β)λ��
YY=L(4)-XX-Unitlength;           %TotalUnit�Σ����һ��0m,���ڿ�λ��
Current_Square=Current.^2*1e6;  %I^2 I��λ��KA!!!
F_EM_rr=zeros(length(t),length(XX));  %length(t)*TotalUnit
F_EM_ar=zeros(length(t),length(XX));
Constant1=miu0*(1/8)*(1/pi)*(1/h);

%ÿʱ�̵����ܵĵ����
%F_EM_rr_total=miu0/4/pi/h*Current_Square*(sqrt(4*h^2+L(4)^2)-2*h);
%F_EM_rr=F_EM_rr_total/L(4);  %��λ���ȵ������ܵ���� ????

for ii=1:length(XX)  %1:TotalUnit�Σ�ÿһ���µ���0-3ms�����ܵ�����ı仯(1-TotalUnit)
    F_EM_rr(:,ii)=Constant1*Current_Square*(sqrt(4*h^2+(Unitlength+YY(ii))^2)+sqrt(4*h^2+(Unitlength+XX(ii))^2) ...
      -sqrt(XX(ii)^2+4*h^2)-sqrt(YY(ii)^2+4*h^2));
end %F_EM_rr(:,i)��i�ε����
% figure
% plot(t,F_EM_rr(:,1))  %0-3s�µ�����ĩ��(�ڿ�)���ܵ����ʱ����������
% title('����β)���ܵ����ʱ����������');
% figure
% plot(XX,F_EM_rr(length(t),:)) %3ms�µ����������
% title('�µ����������')
for jj=1:length(XX)  %1-(TotalUnit-1)�Σ�  !ÿһ��!������µ���ĵ���������һ�β�����
    F_EM_ar(:,jj)=miu0*(1/4)*(1/pi)*Current_Square*(log(2*h^2+h*sqrt(4*h^2+(L(4)-XX(jj)-Unitlength)^2))+log(L(4)-XX(jj)) ...
      -log(L(4)-XX(jj)-Unitlength)-log(2*h^2+h*sqrt(4*h^2+(L(4)-XX(jj))^2)));
end  %F_EM_ar(:,i)��i�ε���� ���������
% figure
% plot(t,F_EM_ar(:,1))  %0-3s�µ�����ĩ��(��β)���ܵ������õ����ʱ����������
% title('2');
% figure
% plot(XX,F_EM_ar(length(t),:)) %3ms�µ����������
% title('������µ���������')
% F_EM=F_EM_rr+F_EM_ar;
F_EM=F_EM_rr;
% F_EM_mean=mean(F_EM,2);
% F_EM=repmat(F_EM_mean,1,2800);
F_EM=ones(3001,2800)*353.628821955585;
% figure
% plot(t,F_EM(:,1))
% title('����β)�����ܵ����ʱ����������');
% figure
% plot(XX,F_EM(length(t),:)) %3ms�µ����������
% title('�µ������������')
end
function Current=Cal_current(t) %�����������  ��λ��ǧ��
t11=(length(t)-1)/6+1;
t1=t(1:t11);
t22=(length(t)-1)/3;
t33=(length(t)-1)/2+2;
t3=t(t33:length(t));     
%%
% %����ʱ��12ms
% Current1=(300*sin(3140*t1/4000))'; %����������
% Current2=300*ones(t22,1);      %�㶨��    ʱ�䵥λת���ڴ��ѿ����ˣ�/1000)
% Current3=(300*exp(-(t3-6)/0.0005/4000))'; 
%%
% %����ʱ��3ms
Current1=300*sin(3140*t1/1000); %����������
Current2=300*ones(t22,1);      %�㶨��    ʱ�䵥λת���ڴ��ѿ����ˣ�/1000)     
Current3=300*exp(-(t3-1.5)/0.0005/1000);   %ָ���½���
%%
Current=[Current1;Current2;Current3];
% figure
% plot(t,Current)
% ylim([0,400])
% title('�������ʱ��仯����')
% xlabel('t/ms')
% ylabel('I/KA')
end
%%
function Test_M_Mode(M_Mode,F_generalize,W,SDG_ZX,SG_ZX,XDG_ZX,V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,t)%��֤ģ̬��������
global TotalUnit Unitlength JieDianShu
q_Mode=zeros(length(t),length(W));  %�����������length(t)*length(W)
 for ii=1:length(W) %���������
   q_Mode(:,ii)=F_generalize(ii)/M_Mode(ii)/W(ii)^2*(1-cos(W(ii)*t));  %q�Ľ���������ʼֵ(��Ϊ0)���
 end
 %ģ̬����
 V_SDG_StaticDisp_Test=SDG_ZX*q_Mode';
 V_SG_StaticDisp_Test=SG_ZX*q_Mode';
 V_XDG_StaticDisp_Test=XDG_ZX*q_Mode';
 
 %��ͼ
XX=(0:TotalUnit)*Unitlength;
figure
plot(XX,V_SDG_StaticDisp_Test(:,end),'--g',XX,V_SG_StaticDisp_Test(:,end),'-.b',XX,V_XDG_StaticDisp_Test(:,end),':r');
title(['t=' num2str(t(end)) 'msʱ�������λ��'])
legend('�ϵ���','���','�µ���')
% figure
% plot(XX,V_SDG_StaticDisp,'--g',XX,V_SG_StaticDisp,'-.b',XX,V_XDG_StaticDisp,':r')
% title('ģ̬���ӷ�������������µľ�λ��')
figure
plot(t,V_SDG_StaticDisp_Test(JieDianShu(4),:),'--g',t,V_SG_StaticDisp_Test(JieDianShu(4),:),'-.b',t,V_XDG_StaticDisp_Test(JieDianShu(4),:),':r');
legend('�ϵ���','���','�µ���')
% hold on
% plot([0,3],[V_SDG_StaticDisp(JieDianShu(4)),V_SDG_StaticDisp(JieDianShu(4))],'-g')
% hold on
% plot([0,3],[V_SG_StaticDisp(JieDianShu(4)),V_SG_StaticDisp(JieDianShu(4))],'-.b')
% hold on
% plot([0,3],[V_XDG_StaticDisp(JieDianShu(4)),V_XDG_StaticDisp(JieDianShu(4))],':r')
title('�ڿ�λ��ʱ������')
legend('�ϵ���','���','�µ���')
xlabel('t/ms')
ylabel('λ��/mm')
hold off
end
function [V_SDG_StaticDisp,V_SG_StaticDisp,V_XDG_StaticDisp,M_Mode,F_generalize]=Cal_StaticDisp(SDG_ZX,SG_ZX,XDG_ZX,W)%��ģ̬���ӷ�������Ӧ(��λ�ƣ�  length(W)�׵���
global m g  Unitlength JieDianShu TotalUnit
M_Mode=zeros(length(W),1);  %ģ̬����
F_generalize=zeros(length(W),1);   %������
q_t=zeros(length(W),1);     %��������
% for rr=1:length(W)   ��ģ��һ��
%     SDG_ZX(:,rr)=SDG_ZX(:,rr)/sqrt(sum(SDG_ZX(:,rr).^2));
%     SG_ZX(:,rr)=SG_ZX(:,rr)/sqrt(sum(SG_ZX(:,rr).^2));
%     XDG_ZX(:,rr)=XDG_ZX(:,rr)/sqrt(sum(XDG_ZX(:,rr).^2));
% end
V_SDG_Square=SDG_ZX.^2;
V_SG_Square=SG_ZX.^2;
V_XDG_Square=XDG_ZX.^2; %JieDianShu(4)*length(W)
C_simpson1=Unitlength/3;  %����ɭ��ʽϵ��
C_simpson2=Unitlength*2/45;
C_simpson3=7*Unitlength/17280;

%���������ģ̬�����͹�������
for ee=1:length(W) 
    V_sdg_Square=V_SDG_Square(:,ee);
    V_sg_Square=V_SG_Square(:,ee);
    V_xdg_Square=V_XDG_Square(:,ee);
    V_sdg=SDG_ZX(:,ee);
    V_sg=SG_ZX(:,ee);
    V_xdg=XDG_ZX(:,ee);
    %%
%     %ֱ����ģ̬�����͹�������������
%     M_Mode1=m(1)*sum(V_sdg_Square);
%     M_Mode2=m(2)*sum(V_sg_Square);
%     M_Mode3=m(3)*sum(V_xdg_Square);
%     F_generalize1=m(1)*g*sum(V_sdg);
%     F_generalize2=m(2)*g*sum(V_sg);
%     F_generalize3=m(3)*g*sum(V_xdg);
    %%
%     %ţ��-�ƴ�����(n=7)
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
%     %��������ɭ32/90����
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
    %��������ɭ��ʽ(1/3����ɭ����
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
%     %�������ι�ʽ
%     M_Mode1=0.004*m(1)*(V_sdg_Square(1)+V_sdg_Square(JieDianShu(4))+2*sum(V_sdg_Square(2:JieDianShu(4)-1)))/2;
%     M_Mode2=0.004*m(2)*(V_sg_Square(1)+V_sg_Square(JieDianShu(4))+2*sum(V_sg_Square(2:JieDianShu(4)-1)))/2;
%     M_Mode3=0.004*m(3)*(V_xdg_Square(1)+V_xdg_Square(JieDianShu(4))+2*sum(V_xdg_Square(2:JieDianShu(4)-1)))/2;
%     F_generalize1=0.004*m(1)*g*(SDG_ZX(1,ee)+SDG_ZX(JieDianShu(4),ee)+2*sum(SDG_ZX(2:JieDianShu(4)-1,ee)))/2;
%     F_generalize2=0.004*m(2)*g*(SG_ZX(1,ee)+SG_ZX(JieDianShu(4),ee)+2*sum(SG_ZX(2:JieDianShu(4)-1,ee)))/2;
%     F_generalize3=0.004*m(3)*(g)*(XDG_ZX(1,ee)+XDG_ZX(JieDianShu(4),ee)+2*sum(XDG_ZX(2:JieDianShu(4)-1,ee)))/2;
   %%
   M_Mode(ee)=M_Mode1+M_Mode2+M_Mode3;
%    F_generalize(ee)=-F_generalize1+F_generalize3; %�ϵ��������������µ�����������
%    F_generalize(ee)=F_generalize2; %ֻ���������
%    F_generalize(ee)=-F_generalize1; %ֻ���ϵ��츺����
%    F_generalize(ee)=F_generalize3; %ֻ���µ�������
   F_generalize(ee)=F_generalize1+F_generalize2+F_generalize3;
   q_t(ee)=F_generalize(ee)/M_Mode(ee)/W(ee)^2;
end

% % %��ģ̬������һ��  V'MV=1
% for rr=1:length(W)
%     SDG_ZX(:,rr)=SDG_ZX(:,rr)/sqrt(sum(SDG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
%     SG_ZX(:,rr)=SG_ZX(:,rr)/sqrt(sum(SG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
%     XDG_ZX(:,rr)=XDG_ZX(:,rr)/sqrt(sum(XDG_ZX(:,rr).^2))/M_Mode(rr)^0.5;
% end

% ģ̬����
V_SDG_StaticDisp=SDG_ZX*q_t;   %JieDianShu(4)*length(w) * length(w)*1
V_SG_StaticDisp=SG_ZX*q_t;
V_XDG_StaticDisp=XDG_ZX*q_t;

% ��ͼ
XX=(0:TotalUnit)*Unitlength;
figure
hold on
plot(XX,V_SDG_StaticDisp,'--g',XX,V_SG_StaticDisp,'-.b',XX,V_XDG_StaticDisp,':r')
%%%% legend('�ϵ��쾲λ��','��ܾ�λ��','�µ��쾲λ��');
%%%%  plot(,SDG_Coordinate,sdg_Height,'y:',SG_Coordinates,sg_Height,'-.',XDG_Coordinates,xdg_Height,'-')
% legend('�ϵ��쾲λ��','��ܾ�λ��','�µ��쾲λ��')%,'comsol�ϵ��쾲λ��','comsol��ܾ�λ��','comsol�µ��쾲λ��');
title('ģ̬���ӷ�������������µľ�λ��')
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('����βλ��/m')
ylabel('��λ��/mm')
hold off
end
%%
function HuaZhenXing(V1,V2,V3)  %������ ��Y)
global TotalUnit Unitlength
X=(0:TotalUnit)*Unitlength;
figure
hold on
plot(X,V1,'--g',X,V2,'-.b',X,V3,':r')
legend('�ϵ�������','�������','�µ�������');
set(gca,'LineWidth',3)
set(gca,'XLim',[0,12])
xlabel('X/m')
ylabel('�Ӷ�/mm')
hold off
end
function [SDG_ZX,SG_ZX,XDG_ZX,Z1_I_all]=Cal_ZX(W) %������ (Y)
global  JieDianShu Unitlength
Z1_I_all=zeros(12,length(W));
SDG_ZX=zeros(JieDianShu(4),length(W));   
SG_ZX=zeros(JieDianShu(4),length(W));
XDG_ZX=zeros(JieDianShu(4),length(W));  %ǰlength(W)�����ͣ�JieDianShu(4)=2801���ڵ� 1361+881+561-2=2801�� TotalUnit��
B=zeros((JieDianShu(1)-1)*12,12);  %�洢���ݾ��󣬽�Լ������
for aa=1:length(W)  %��aa�׹�������
    w=W(aa);
    [U_all,U1,U2,U3,U4,~,A_exp]=Cal_U(w);
    U_guiyi=[U_all([3,4,7,8,11,12],[1,2,5,6,9,10]);
                   1 1 1 1 1 1 ];   %7*6
    %[1 2 5 6 9 10]:y1,��1��y2,��2,y3,��3
    %[U_bar;[1...1]]��������ֵ�ֽ�����Z_bar����Ϊ����ֵ�ֽ�õ�������ֵ����S�����ȵ�
    Z_bar=(U_guiyi'*U_guiyi)\(U_guiyi')*[0 0 0 0 0 0 1]';% (6*7)*(7*6)*��6*7��*��7*1��
    %����ģ̬��������ʱ��һ������
    Z1_I=[Z_bar(1) Z_bar(2) 0 0 Z_bar(3) Z_bar(4) 0 0 Z_bar(5) Z_bar(6) 0 0]';
  
    Z1_I_all(:,aa)=Z1_I;
    SDG_ZX(1,aa)=Z1_I(1);
    SG_ZX(1,aa)=Z1_I(5);
    XDG_ZX(1,aa)=Z1_I(9);
    
%   ��������������ȫ����Yȫ���Ħ�ȫ����Mȫ����Q�����У��õ�ʱ���ÿһ��ժ������
%   X1=((1:(JieDianShu(1)-1))*0.08)';   %��Ԫ����0.08   67*1  ���һ��X1�Ĵ��ݾ������U1
%   U1_X1=Cal_U_elasticbeam(m,EI,w,k_Elasticbase,X1);       %(67*12)*12
%   ��һ�δ��ݾ���  A��12*12�ģ���ô��X1��67*1���������????
%   %Z1_X1=U1_X1(1:67,1:4)*Z1_I;      %  67*1  ��һ�������״̬ʸ��
%   SG_ZX(2:68,i)=U1_X1(1:67,1:4)*Z1_I;     %��һ��  �ڵ�ţ�2-68)

%   ѭ��
    for jj=2:JieDianShu(1)  %2-1361 JieDianShu(1)=1361
        X1=Unitlength*(jj-1);    %��Ԫ����0.004
        U1_X1=expm(A_exp*X1);    %  ���һ��X1�Ĵ��ݾ������U1
        B(12*(jj-2)+1:12*(jj-1),1:12)=U1_X1;  %�洢���ݾ���
        SDG_ZX(jj,aa)=U1_X1(1,:)*Z1_I;
        SG_ZX(jj,aa)=U1_X1(5,:)*Z1_I;   %2-1361
        XDG_ZX(jj,aa)=U1_X1(9,:)*Z1_I;
    end   %ѭ�������õ�������һ�δӵ�2���㵽��1361���������λ��ֵ
     
    Z3_I=U2*U1*Z1_I;  %��������; 12*1
    for dd=1:(JieDianShu(2)-1)    %JieDianShu(2)=881
        U3_X2=B(12*(dd-1)+1:12*dd,1:12);   %��ȡ���ݾ���  ���һ��U3_X2�Ĵ��ݾ������U3
        SDG_ZX(JieDianShu(1)+dd,aa)=U3_X2(1,:)*Z3_I;
        SG_ZX(JieDianShu(1)+dd,aa)=U3_X2(5,:)*Z3_I;
        XDG_ZX(JieDianShu(1)+dd,aa)=U3_X2(9,:)*Z3_I;  %1362-2241 
    end %ѭ�������õ������ڶ��δӵ�1362���㵽��2241���������λ��ֵ
    
    Z5_I=U4*U3*Z3_I;
    for bb=2:JieDianShu(3)      %JieDianShu(3)=561
        U5_X3=B(12*(bb-2)+1:12*(bb-1),1:12); %��ȡ���ݾ���  ���һ��U5_X3�Ĵ��ݾ������U5
        SDG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(1,:)*Z5_I;
        SG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(5,:)*Z5_I;
        XDG_ZX(JieDianShu(1)+JieDianShu(2)+bb-2,aa)=U5_X3(9,:)*Z5_I;  %2241-2801
    end  %ѭ�������õ����������δӵ�2241���㵽��2801���������λ��ֵ
end
end  
function w=Cal_delta(Pre,Step,n)   %�������Ƶ��(���ȣ�������������  (Y)
        w=zeros(n,1);  %��ǰ��n�׹���Ƶ�ʷ����ڴ��ַ 
        a=0 ;    %�㵽�ڼ��׹���Ƶ��
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
          %fun1*fun2<0   �ж�����û�б�Ҫ�ˣ����������һ���ж���������Ȼ����<0
                  x0=(x1+x2)/2;
                  fun0=Cal(x0);
                  if fun0*fun1>0
                     x1=x0;
                     fun1=fun0;
                  else
                     x2=x0;
                     fun2=fun0;
                  end
         end   %�˴�ѭ������ʱ��״̬Ϊabs(x2-x1)<Pre;fun1*fun2<0
         a=a+1;
         w(a)=(x1+x2)/2;  %�õ���a�׹���Ƶ��
         x1=x2;fun1=fun2;
         x2=x1+Step;fun2=Cal(x2);   %�ڴ˻����ϼ��������㣬����һ�׹���Ƶ��
        end
end  
function delta=Cal(w)  %�ֿ����(Y)
[U_all,~]=Cal_U(w);
U_bar=U_all([3,4,7,8,11,12],[1,2,5,6,9,10]); %6*6 ����߽���������������
delta=det(U_bar);
end
function [U_all,U1,U2,U3,U4,U5,A_exp]=Cal_U(w)  %ƴ�ܴ��ݾ���(Y)
global EI m k_Elasticbase k_Springsupport L
A_exp=Pin_exp(m,EI,w,k_Elasticbase(1),k_Elasticbase(2));
U1=expm(A_exp*L(1));
U2=Cal_Springsupport(k_Springsupport);
U3=expm(A_exp*L(2));
U4=Cal_Springsupport(k_Springsupport);
U5=expm(A_exp*L(3));      
U_all=U5*U4*U3*U2*U1; %�ܴ��ݾ���
end  
function TM_Springsupport=Cal_Springsupport(k)  %����Ԫ�����ӵ���֧�ŵĴ��ݾ���  (Y)
TM_Springsupport=eye(12,12);
TM_Springsupport(8,5)=-k;
end
function A_exp=Pin_exp(m,EI,w,k1,k2)         %΢�ַ������ϵ������A  (Y)
%k1 �ϵ�������ܼ�ֲ����Ըն�    k2  �µ�������ܼ�ֲ����Ըն�
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