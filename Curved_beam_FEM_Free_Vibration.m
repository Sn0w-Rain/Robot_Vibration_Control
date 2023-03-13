% ****************************

% 作者: Zhihui Liu
% 程序简介: 这个程序用来进行计算弯曲的梁的自由振动特性
% 时间:2020-09-20 20:40:54

%****************************
tic
clc;
clear('all');
node_num=100;
element_num=node_num-1;
rho=7850;
b=0.2;
h=0.1;
A=b*h;
Ee=200e9;
Ii=b*h^3/12;
Mp=1;



kr=Ee*Ii*50;






a=1;
b=0.3;

% f=@(x)sqrt(b^2-b^2/a^2*x.^2);
f=@(x)sin(pi*x);
xlist=linspace(0,5,node_num);
% ylist=zeros(size(xlist));
ylist=f(xlist);

% tip_distance=sqrt((xlist(end)+Rh)^2+ylist(end)^2);
tip_distance=1;

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

Tmat=@(i)blkdiag(...
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
Kmat=zeros(3*node_num);
Mmat=zeros(3*node_num);


for i=1:element_num

    T=Tmat(i);
    elen=element_lengths(i);
    xm=x_mid(i);
    ym=y_mid(i);

    xc=x_cos(i);
    yc=y_cos(i);

    xs=x_sin(i);
    ys=y_sin(i);


    tmp_Kmat=...
    [A.*Ee.*elen.^(-1),0,0,(-1).*A.*Ee.*elen.^(-1),0,0;0,12.*Ee.* ...
      elen.^(-3).*Ii,6.*Ee.*elen.^(-2).*Ii,0,(-12).*Ee.*elen.^(-3).*Ii, ...
      6.*Ee.*elen.^(-2).*Ii;0,6.*Ee.*elen.^(-2).*Ii,4.*Ee.*elen.^(-1).* ...
      Ii,0,(-6).*Ee.*elen.^(-2).*Ii,2.*Ee.*elen.^(-1).*Ii;(-1).*A.*Ee.* ...
      elen.^(-1),0,0,A.*Ee.*elen.^(-1),0,0;0,(-12).*Ee.*elen.^(-3).*Ii,( ...
      -6).*Ee.*elen.^(-2).*Ii,0,12.*Ee.*elen.^(-3).*Ii,(-6).*Ee.*elen.^( ...
      -2).*Ii;0,6.*Ee.*elen.^(-2).*Ii,2.*Ee.*elen.^(-1).*Ii,0,(-6).*Ee.* ...
      elen.^(-2).*Ii,4.*Ee.*elen.^(-1).*Ii];

  tmp_Mmat=...
          [(1/3).*A.*elen.*rho,0,0,(1/6).*A.*elen.*rho,0,0;0,(13/35).*A.* ...
            elen.*rho,(11/210).*A.*elen.^2.*rho,0,(9/70).*A.*elen.*rho,( ...
            -13/420).*A.*elen.^2.*rho;0,(11/210).*A.*elen.^2.*rho,(1/105).*A.* ...
            elen.^3.*rho,0,(13/420).*A.*elen.^2.*rho,(-1/140).*A.*elen.^3.* ...
            rho;(1/6).*A.*elen.*rho,0,0,(1/3).*A.*elen.*rho,0,0;0,(9/70).*A.* ...
            elen.*rho,(13/420).*A.*elen.^2.*rho,0,(13/35).*A.*elen.*rho,( ...
            -11/210).*A.*elen.^2.*rho;0,(-13/420).*A.*elen.^2.*rho,(-1/140).* ...
            A.*elen.^3.*rho,0,(-11/210).*A.*elen.^2.*rho,(1/105).*A.*elen.^3.* ...
            rho];

      tmp_Kmat=T.'*tmp_Kmat*T;
      tmp_Mmat=T.'*tmp_Mmat*T;

    for j=1:6
        for k=1:6
          jmap=j+(i-1)*3;
          kmap=k+(i-1)*3;
          Kmat(jmap,kmap)=Kmat(jmap,kmap)+tmp_Kmat(j,k);
          Mmat(jmap,kmap)=Mmat(jmap,kmap)+tmp_Mmat(j,k);
        end
    end

    if abs(i-element_num)<0.1
        xm=x_end(end);
        ym=y_end(end);
        extra_mass_Mmat=...
        [0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,Mp,0,0;0,0,0,0,Mp,0;0, ...
          0,0,0,0,0];
        extra_rotation_Mmat=...
        [0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0;0,0, ...
          0,0,0,Jp];

        extra_mass_Mmat=T.'*extra_mass_Mmat*T;
        extra_rotation_Mmat=T.'*extra_rotation_Mmat*T;
        for j=1:6
            for k=1:6

              jmap=j+(i-1)*3;
              kmap=k+(i-1)*3;
              Mmat(jmap,kmap)=Mmat(jmap,kmap)+...
              extra_mass_Mmat(j,k)+extra_rotation_Mmat(j,k);
            end
        end
    end

end



% 对rotor进行处理

% 对hub进行考虑


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







Kmat(3,:)=[];
Kmat(:,3)=[];

Mmat(3,:)=[];
Mmat(:,3)=[];


Kmat(2,:)=[];
Kmat(:,2)=[];

Mmat(2,:)=[];
Mmat(:,2)=[];

Kmat(1,:)=[];
Kmat(:,1)=[];

Mmat(1,:)=[];
Mmat(:,1)=[];
% Kmat=sparse(Kmat);
% Mmat=sparse(Mmat);


[V,D]=eig(Kmat,Mmat);
[eigvalue,ind]=sort(diag(D));
eigvec=V(:,ind);
freq=abs(sqrt(eigvalue)/2/pi);


% 提取有关弯曲振动的
aa=eigvec.'*Kmat*eigvec;
toc