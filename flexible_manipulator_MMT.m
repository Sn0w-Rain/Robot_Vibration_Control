% ****************************

% 作者: Zhihui Liu
% 程序简介: 这个程序用来进行计算弯曲的梁的自由振动特性
% 时间:2020-09-20 20:40:54

%****************************
tic
clc;
clear('all');
krxishu=30;%定义关节的关节刚度
curve=0;
% 0是直的，1是弯的
mode_num=2;
nlt=1;
% 0是不考虑非线性，1是考虑非线性的

element_num=50;%确定单元的数目
node_num=element_num+1;%这个是节点的数目
rho=2700;
% b=0.05;
b=0.0313219;
% h=0.005;
h=0.0048481;
A=b*h;
Ee=70e9;
Ii=b*h^3/12;


eta=0.5;
dt=1e-4;
tend=5;
tlist=0:dt:tend;
tnum=length(tlist);
tau0=0.5;
T0=4;




Rh=0.0;
N=100;

% kr=Ee*Ii*krxishu;
kr=624.6;

% Mp=5;
% Jp=0.02;
Mp=0;
Jp=0;
Jr=1.73e-5;
Jm=N^2*Jr;
Jh=0.000;


% f=@(x)sqrt(b^2-b^2/a^2*x.^2);
f=@(x)0.1*sin(pi*x);
xlist=linspace(0,1,node_num);
ylist=zeros(size(xlist));
% ylist=f(xlist);

x_start=xlist(1:end-1);
x_end=xlist(2:end);
x_mid=x_start+(x_end-x_start)/2;

y_start=ylist(1:end-1);
y_end=ylist(2:end);
y_mid=y_start+(y_end-y_start)/2;

element_lengths=((x_end - x_start).^2+(y_end - y_start).^2).^0.5;


middle_point_coord=[x_mid;y_mid];

x_cos=(x_end-x_start)./element_lengths;
x_sin=(y_end-y_start)./element_lengths;

y_cos=-x_sin;
y_sin=x_cos;

coor_R=@(i)[...
    x_cos(i) , x_sin(i);
    y_cos(i) , y_sin(i)];

Tmat=@(i)blkdiag(1,...
    [...
    x_cos(i) , x_sin(i) , 0;
    y_cos(i) , y_sin(i) , 0;
    0        , 0        , 1
    ],...
    [...
    x_cos(i) , x_sin(i) , 0;
    y_cos(i) , y_sin(i) , 0;
    0        , 0        , 1
    ]);
Kmat=zeros(3*node_num+2);
Mmat=zeros(3*node_num+2);
Cmat=zeros(3*node_num+2);
Gmat=zeros(3*node_num+2,1);


for i=1:element_num

    T=Tmat(i);
    rotation_R=coor_R(i);
    coor_in_local=rotation_R*middle_point_coord(:,i);

    xm=coor_in_local(1);
    ym=coor_in_local(2);
    elen=element_lengths(i);


    tmp_Kmat=...
    [0,0,0,0,0,0,0;0,A.*Ee.*elen.^(-1),0,0,(-1).*A.*Ee.*elen.^(-1),0, ...
      0;0,0,12.*Ee.*elen.^(-3).*Ii,6.*Ee.*elen.^(-2).*Ii,0,(-12).*Ee.* ...
      elen.^(-3).*Ii,6.*Ee.*elen.^(-2).*Ii;0,0,6.*Ee.*elen.^(-2).*Ii,4.* ...
      Ee.*elen.^(-1).*Ii,0,(-6).*Ee.*elen.^(-2).*Ii,2.*Ee.*elen.^(-1).* ...
      Ii;0,(-1).*A.*Ee.*elen.^(-1),0,0,A.*Ee.*elen.^(-1),0,0;0,0,(-12).* ...
      Ee.*elen.^(-3).*Ii,(-6).*Ee.*elen.^(-2).*Ii,0,12.*Ee.*elen.^(-3).* ...
      Ii,(-6).*Ee.*elen.^(-2).*Ii;0,0,6.*Ee.*elen.^(-2).*Ii,2.*Ee.* ...
      elen.^(-1).*Ii,0,(-6).*Ee.*elen.^(-2).*Ii,4.*Ee.*elen.^(-1).*Ii]; ...

    tmp_Mmat=...
          [(1/12).*A.*elen.^3.*rho+A.*elen.*rho.*xm.^2+A.*elen.*rho.*ym.^2,( ...
            -1/2).*A.*elen.*rho.*ym,(-1/10).*A.*elen.^2.*rho+(1/2).*A.*elen.* ...
            rho.*xm,(-1/120).*A.*elen.^3.*rho+(1/12).*A.*elen.^2.*rho.*xm,( ...
            -1/2).*A.*elen.*rho.*ym,(1/10).*A.*elen.^2.*rho+(1/2).*A.*elen.* ...
            rho.*xm,(-1/120).*A.*elen.^3.*rho+(-1/12).*A.*elen.^2.*rho.*xm;( ...
            -1/2).*A.*elen.*rho.*ym,(1/3).*A.*elen.*rho,0,0,(1/6).*A.*elen.* ...
            rho,0,0;(-1/10).*A.*elen.^2.*rho+(1/2).*A.*elen.*rho.*xm,0,(13/35) ...
            .*A.*elen.*rho,(11/210).*A.*elen.^2.*rho,0,(9/70).*A.*elen.*rho,( ...
            -13/420).*A.*elen.^2.*rho;(-1/120).*A.*elen.^3.*rho+(1/12).*A.* ...
            elen.^2.*rho.*xm,0,(11/210).*A.*elen.^2.*rho,(1/105).*A.*elen.^3.* ...
            rho,0,(13/420).*A.*elen.^2.*rho,(-1/140).*A.*elen.^3.*rho;(-1/2).* ...
            A.*elen.*rho.*ym,(1/6).*A.*elen.*rho,0,0,(1/3).*A.*elen.*rho,0,0;( ...
            1/10).*A.*elen.^2.*rho+(1/2).*A.*elen.*rho.*xm,0,(9/70).*A.*elen.* ...
            rho,(13/420).*A.*elen.^2.*rho,0,(13/35).*A.*elen.*rho,(-11/210).* ...
            A.*elen.^2.*rho;(-1/120).*A.*elen.^3.*rho+(-1/12).*A.*elen.^2.* ...
            rho.*xm,0,(-13/420).*A.*elen.^2.*rho,(-1/140).*A.*elen.^3.*rho,0,( ...
            -11/210).*A.*elen.^2.*rho,(1/105).*A.*elen.^3.*rho];


      tmp_Kmat=T.'*tmp_Kmat*T;
      tmp_Mmat=T.'*tmp_Mmat*T;



    for j=1:7
        for k=1:7
            if j==1
              jmap=2;
            else
              jmap=3*(i-1)+j+1;
            end

            if k==1
              kmap=2;
            else
              kmap=3*(i-1)+k+1;
            end

            Kmat(jmap,kmap)=Kmat(jmap,kmap)+tmp_Kmat(j,k);
            Mmat(jmap,kmap)=Mmat(jmap,kmap)+tmp_Mmat(j,k);
        end
    end

    if abs(i-element_num)<0.1
        T=Tmat(i);
        rotation_R=coor_R(i);
        coor_in_local=rotation_R*[x_mid(end);y_mid(end)];
        % coor_in_local=rotation_R*middle_point_coord(:,i);
        elen=element_lengths(i);
        xm=coor_in_local(1);
        ym=coor_in_local(2);



        extra_mass_Mmat=...
        [Jp+Mp.*(((1/2).*elen+xm).^2+ym.^2),0,0,0,(-1).*Mp.*ym,Mp.*((1/2) ...
          .*elen+xm),Jp;0,0,0,0,0,0,0;0,0,0,0,0,0,0;0,0,0,0,0,0,0;(-1).*Mp.* ...
          ym,0,0,0,Mp,0,0;Mp.*((1/2).*elen+xm),0,0,0,0,Mp,0;Jp,0,0,0,0,0, ...
          Jp];


        extra_mass_Mmat=T.'*extra_mass_Mmat*T;
        for j=1:7
            for k=1:7
                if j==1
                  jmap=2;
                else
                  jmap=3*(i-1)+j+1;
                end

                if k==1
                  kmap=2;
                else
                  kmap=3*(i-1)+k+1;
                end

                Mmat(jmap,kmap)=Mmat(jmap,kmap)+...
                extra_mass_Mmat(j,k);
            end
        end
    end


    if nlt==1
      tmp_Cmat=...
          -2*rho*A*[0,(1/12).*elen.^2+(-1/2).*elen.*xm,(-1/2).*elen.*ym,(-1/12).* ...
            elen.^2.*ym,(-1/12).*elen.^2+(-1/2).*elen.*xm,(-1/2).*elen.*ym,( ...
            1/12).*elen.^2.*ym;0,0,(7/20).*elen,(1/20).*elen.^2,0,(3/20).* ...
            elen,(-1/30).*elen.^2;0,(-7/20).*elen,0,0,(-3/20).*elen,0,0;0,( ...
            -1/20).*elen.^2,0,0,(-1/30).*elen.^2,0,0;0,0,(3/20).*elen,(1/30).* ...
            elen.^2,0,(7/20).*elen,(-1/20).*elen.^2;0,(-3/20).*elen,0,0,( ...
            -7/20).*elen,0,0;0,(1/30).*elen.^2,0,0,(1/20).*elen.^2,0,0];

      tmp_Gmat=...
        -rho*A*[0;(-1/12).*elen.*(elen+(-6).*xm);(1/2).*elen.*ym;(1/12).* ...
          elen.^2.*ym;(1/12).*elen.*(elen+6.*xm);(1/2).*elen.*ym;(-1/12).* ...
          elen.^2.*ym];

      tmp_Cmat=T.'*tmp_Cmat*T;
      tmp_Gmat=T.'*tmp_Gmat;

      for j=1:7
          for k=1:7
              if j==1
                jmap=2;
              else
                jmap=3*(i-1)+j+1;
              end

              if k==1
                kmap=2;
              else
                kmap=3*(i-1)+k+1;
              end
              Cmat(jmap,kmap)=Cmat(jmap,kmap)+tmp_Cmat(j,k);
          end
          Gmat(jmap,1)=Gmat(jmap,1)+tmp_Gmat(j,1);
      end

      if abs(i-element_num)<0.1
          T=Tmat(i);
          rotation_R=coor_R(i);
          coor_in_local=rotation_R*[x_mid(end);y_mid(end)];
          % coor_in_local=rotation_R*middle_point_coord(:,i);
          elen=element_lengths(i);
          xm=coor_in_local(1);
          ym=coor_in_local(2);



          extra_C_Mmat=...
            T.'*(-2*Mp*[0,0,0,0,(-1/2).*elen+(-1).*xm,(-1).*ym,0;0,0,0,0,0,0,0;0,0,0,0,0, ...
                      0,0;0,0,0,0,0,0,0;0,0,0,0,0,1,0;0,0,0,0,(-1),0,0;0,0,0,0,0,0,0])*T;

          extra_G_Mmat=...
            T.'*(-Mp*[0;0;0;0;(1/2).*elen+xm;ym;0]);

          for j=1:7
              for k=1:7
                  if j==1
                    jmap=2;
                  else
                    jmap=3*(i-1)+j+1;
                  end

                  if k==1
                    kmap=2;
                  else
                    kmap=3*(i-1)+k+1;
                  end

                  Cmat(jmap,kmap)=Cmat(jmap,kmap)+...
                  extra_C_Mmat(j,k);
              end
            Gmat(jmap,1)=Gmat(jmap,1)+extra_G_Mmat(j,1);
          end
      end
    end


end



% 对rotor进行处理
Mmat(1,1)=Mmat(1,1)+N^2*Jr;
% 对hub进行考虑
Mmat(2,2)=Mmat(2,2)+Jh;

%考虑关节柔性
Kmat(1:2,1:2)=Kmat(1:2,1:2)+kr*([1,-1].*[1,-1]');


Kmat(5,:)=[];
Kmat(:,5)=[];
Kmat(4,:)=[];
Kmat(:,4)=[];
Kmat(3,:)=[];
Kmat(:,3)=[];

Mmat(5,:)=[];
Mmat(:,5)=[];
Mmat(4,:)=[];
Mmat(:,4)=[];
Mmat(3,:)=[];
Mmat(:,3)=[];

Cmat(5,:)=[];
Cmat(:,5)=[];
Cmat(4,:)=[];
Cmat(:,4)=[];
Cmat(3,:)=[];
Cmat(:,3)=[];

Gmat(5,:)=[];
Gmat(4,:)=[];
Gmat(3,:)=[];

[hang,lie]=size(Kmat);

[V,D]=eig(Kmat,Mmat);
[eigvalue,ind]=sort(diag(D));
eigvec=V(:,ind);
freq=abs(sqrt(eigvalue)/2/pi);

modes=eigvec(:,1:mode_num);


Kmat=modes.'*Kmat*modes;
Mmat=modes.'*Mmat*modes;
Cmat=modes.'*Cmat*modes;
Gmat=modes.'*Gmat;


%%下面开始进行时域内的求解

 % 开始进行时域内的求解

%  @(t)modes'*[tau0.*((mod(t,1)-0.5)<=0)+-tau0.*((mod(t,1)-0.5)>0);zeros(hang-1,1)];
% tau=@(t)modes'*[tau0*sin(2*pi/T0*t).*(t<4);zeros(hang-1,1)];
tau=@(t)modes'*[sin(2*pi*t)*(t<4);zeros(hang-1,1)];
% tau=@(t)modes'*[tau0*t*(t<1)+tau0*(t>=1);zeros(hang-1,1)];
% tau=@(t)modes'*[((t-2)<0)*t+((t-2)>=0)*2;zeros(hang-1,1)];


alpha_m=(2*eta-1)/(eta+1);
alpha_f=eta/(1+eta);
gamma=1/2-alpha_m+alpha_f;
beta=1/4*(1-alpha_m+alpha_f)^2;
a1=(1-alpha_m)/(beta*dt^2);
a2=gamma*(1-alpha_f)/(beta*dt);
a3=(1-alpha_m)/(beta*dt);
a4=((1-alpha_f)*gamma-beta)/beta;
a5=(1-alpha_m-2*beta)/2/beta;
a6=(1-alpha_f)*(2*beta-gamma)*dt/2/beta;

d_mat=zeros(mode_num,tnum);
v_mat=zeros(mode_num,tnum);
a_mat=zeros(mode_num,tnum);
v_d_ratio=gamma/(beta*dt);
a_mat(:,1)=linsolve(Mmat,tau(0));
% 定义一个残差的函数
residual=@(t,k_a,k_v,k_d,k_dtheta,a0,v0,d0,dtheta0)tau(t)-(...
  Mmat*((1-alpha_m)*k_a+alpha_m*a0)+...
  Kmat*((1-alpha_f)*k_d+alpha_f*d0)+...
  (1-alpha_f)*(k_dtheta*Cmat*k_v+(k_dtheta)^2*Gmat)+...
  alpha_f*(dtheta0*Cmat*v0+(dtheta0)^2*Gmat)...
  );

iterative_num=zeros(1,tnum);
for i=2:tnum
  d0=d_mat(:,i-1);
  v0=v_mat(:,i-1);
  a0=a_mat(:,i-1);
  dtheta0=modes(2,:)*v0;
  t=(1-alpha_f)*tlist(i)+alpha_f*tlist(i-1);

  % 首先预估一个初步的解，在这个解的基础上继续进行迭代
  k_d=d0;
  k_a=((k_d-d0)/(dt^2)-v0/dt-(0.5-beta)*a0)/beta;
  k_v=v0+(1-gamma)*dt*a0+gamma*dt*k_a;

  k_dtheta=modes(2,:)*k_v;

  R=residual(t,k_a,k_v,k_d,k_dtheta,a0,v0,d0,dtheta0);
  counter=0;
  while norm(R)>1e-5
    counter=counter+1;
    if counter>100
      break;
      disp('This solution doesn''t convergent');
    end
    tangent_Mat=a1*Mmat+(1-alpha_f)*Kmat+...
      (1-alpha_f)*v_d_ratio*(k_dtheta*Cmat)+...
      (1-alpha_f)*v_d_ratio*(Cmat*k_v).*modes(2,:)+...
      (1-alpha_f)*2*k_dtheta*Gmat*v_d_ratio.*modes(2,:);
    increment_d=linsolve(tangent_Mat,R);
    k_d=k_d+increment_d;
    k_a=((k_d-d0)/(dt^2)-v0/dt-(0.5-beta)*a0)/beta;
    k_v=v0+(1-gamma)*dt*a0+gamma*dt*k_a;
    k_dtheta=modes(2,:)*k_v;
    R=residual(t,k_a,k_v,k_d,k_dtheta,a0,v0,d0,dtheta0);
  end

  iterative_num(i)=counter;

  d_mat(:,i)=k_d;
  v_mat(:,i)=k_v;
  a_mat(:,i)=k_a;

end
results=modes*d_mat;
% close all
hold all
tip_def=results(end-1,:);
plot(tlist,tip_def,'r--');
data=[tlist',tip_def'];

% theta=results(1,:);
% theta_h=results(2,:);
% x_pos=1*cos(theta_h)-tip_def.*sin(theta_h);
% y_pos=1*sin(theta_h)+tip_def.*cos(theta_h);
% figure
% plot(tlist,x_pos)
% hold on
% plot(tlist,y_pos)


WinOnTop
toc