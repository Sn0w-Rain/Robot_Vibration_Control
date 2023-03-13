% ****************************

% 作者: Zhihui Liu
% 程序简介: 这个程序用来进行计算弯曲的梁的自由振动特性
% 时间:2020-09-20 20:40:54

%****************************
% tic
% clc;
% clear('all');
% krxishu=10;

% element_num=100;%确定单元的数目
% node_num=element_num+1;%这个是节点的数目
% rho=2700;
% b=20e-3;
% h=4e-3;
% A=b*h;
% Ee=70e9;
% Ii=b*h^3/12;
% Mp=1;
% Rh=0.03;




% N=100;

% kr=Ee*Ii*krxishu;


% Jr=1e-5;
% Jm=N^2*Jr;
% Jh=2e-3;
% Jp=5e-3;


% % f=@(x)sqrt(b^2-b^2/a^2*x.^2);
% f=@(x)0.1*sin(pi*x);
% % xlist=linspace(0,1,node_num);
% xlist=linspace(Rh,1,node_num);
% % ylist=zeros(size(xlist));
% ylist=f(xlist);

tic
clc;
clear('all');
krxishu=10;
rho=2700;
b=0.02;
h=0.004;
A=b*h;
Ee=70e9;
Ii=b*h^3/12;
Mp=1;
Rh=0.00;
N=100;

kr=Ee*Ii*krxishu;

Jr=1.0e-5;
Jm=N^2*Jr;
Jh=0.002;
Jp=0.005;

freq_list=[];
omega_list=[];

element_num=50;
% element_num=20;%确定单元的数目
node_num=element_num+1;%这个是节点的数目

% f=@(x)sqrt(b^2-b^2/a^2*x.^2);
f=@(x)0.1*sin(pi*x);
xlist=linspace(Rh,1,node_num);
% ylist=zeros(size(xlist));
ylist=f(xlist);




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


for i=1:element_num

    T=Tmat(i);
    rotation_R=coor_R(i);
    coor_in_local=rotation_R*middle_point_coord(:,i);

    xm=coor_in_local(1);
    ym=coor_in_local(2);
    elen=element_lengths(i);
    % xm=x_mid(i);
    % ym=y_mid(i);

    % xc=x_cos(i);
    % yc=y_cos(i);

    % xs=x_sin(i);
    % ys=y_sin(i);

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

end



% 对rotor进行处理
Mmat(1,1)=Mmat(1,1)+N^2*Jr;
% 对hub进行考虑
Mmat(2,2)=Mmat(2,2)+Jh;

% 对负载质量进行处理
% Mmat(2,2)=Mmat(2,2)+Mp*tip_distance^2;
% Mmat(2,end-1)=Mmat(2,end-1)+Mp*tip_distance;
% Mmat(end-1,2)=Mmat(end-1,2)+Mp*tip_distance;
% Mmat(end-1,end-1)=Mmat(end-1,end-1)+Mp;


% % % 对负载惯性进行处理
% Mmat(2,2)=Mmat(2,2)+Jp;
% Mmat(2,end)=Mmat(2,end)+Jp;
% Mmat(end,2)=Mmat(end,2)+Jp;
% Mmat(end,end)=Mmat(end,end)+Jp;

%考虑关节柔性
Kmat(1:2,1:2)=Kmat(1:2,1:2)+kr*([1,-1].*[1,-1]');


Kmat(5,:)=[];
Kmat(:,5)=[];

Mmat(5,:)=[];
Mmat(:,5)=[];


Kmat(4,:)=[];
Kmat(:,4)=[];

Mmat(4,:)=[];
Mmat(:,4)=[];


Kmat(3,:)=[];
Kmat(:,3)=[];

Mmat(3,:)=[];
Mmat(:,3)=[];

Mmat(1,:)=[];
Mmat(:,1)=[];


Kmat(1,:)=[];
Kmat(:,1)=[];


[V,D]=eig(Kmat,Mmat);
[eigvalue,ind]=sort(diag(D));
eigvec=V(:,ind);
freq=abs(sqrt(eigvalue)/2/pi);

[hang,lie]=size(eigvec(1,:));
eigvec=[eigvec(1,:);zeros(3,lie);eigvec(2:end,:)];

omegas=freq*2*pi;

plot(xlist,ylist)
toc