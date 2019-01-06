clear
% declare constants
g=9.8;
R=6371000; % Earth radius (m)
omega=2*pi/(24*60*60); % Earth rotation angle velocity
% input data
OceanData=xlsread('woa13_2005-2012winter_density.csv');
size_OceanData=size(OceanData,1);
Dep=OceanData(1,3:49); %depth   1000 m at 49 
y=OceanData(2:size_OceanData,1); %latitude
x=OceanData(2:size_OceanData,2); %longitude

%find the region 173~179 E 28~35 S
j=0;
for i=1:size(y,1)
    if (x(i,1)<179.5 && x(i,1)>173.5)
        if (y(i,1)<-28 && y(i,1)>-35)  
        j=j+1;
        X(j,1)=x(i,1);
        Y(j,1)=y(i,1);
        % Calculation pressure at depth H
        p(j)=0;
        den(j,1)=OceanData(i+1,3);
        for k=2:size(Dep,2)-1
            %location i, at k-th depth
            den(j,k)=OceanData(i+1,k+2);
            % Numerical integration
            p(j)=0.5*den(j,1)*g*Dep(1,2);
            p(j)=p(j)+0.5*den(j,k)*g*(Dep(1,k+1)-Dep(1,k-1));
        end
        den(j,size(Dep,2))=OceanData(i+1,size(Dep,2)+2);
        p(j)=p(j)+0.5*den(j,size(Dep,2))*g*(Dep(1,size(Dep,2))-...
            Dep(1,size(Dep,2)-1))+1000*g*Dep(1,size(Dep,2));
        end
    end
end
clear OceanData
% Reconstruction x, y, p to 2D matrix
temp=X(1,1);
for i=2:size(X,1)
   if (X(i,1)==temp) 
       Nx=i-1;
       break; 
   end
end
Ny=size(X,1)/Nx;
for i=1:Ny
    for j=1:Nx
        geo_x(i,j)=X((i-1)*Nx+j,1);
        geo_y(i,j)=Y((i-1)*Nx+j,1);
        P(i,j)=p(1,(i-1)*Nx+j);
        Den_H(i,j)=den((i-1)*Nx+j,size(Dep,2))+1000;
    end
end

figure(10)
[C,h] = contourf(geo_x,-geo_y,P*10.0^(-5));
axis equal;
set(h,'LineColor','none')
set(gca,'Ydir','reverse');
title('Pressure (bar)','FontSize',15);
xlabel('E','FontSize',12);
ylabel('S','FontSize',12);
shading interp; 
colorbar;
colormap('jet');
Mat2VTK('Pressure_EAC_summer_dep1000.vtk',P,'binary');
% Estimate the horizontal velocity field at depth H
for i=2:size(geo_y,1)-1
    f=2*omega*sind(geo_y(i,1));
    for j=2:size(geo_x,2)-1
        % x, y: latitude in degree
        dx=(geo_x(i,j+1)-geo_x(i,j-1))*(R*cosd(geo_y(i,j))*pi/180);   
        dy=(geo_y(i+1,j)-geo_y(i-1,j))*(R*pi/180);  % unit: (m)
        % v=g/f (partial eta/partial x)
        v(i-1,j-1)=g/f*(P(i,j+1)-P(i,j-1))/dx/Den_H(i,j);
        % u=-g/f (partial eta/partial y)
        u(i-1,j-1)=-g/f*(P(i+1,j)-P(i-1,j))/dy/Den_H(i,j);
    end
end
figure(69)
quiver(geo_x(2:size(geo_y,1)-1,2:size(geo_x,2)-1),...
    geo_y(2:size(geo_y,1)-1,2:size(geo_x,2)-1),u,v,'LineWidth',1)
%title('Depth = 1000 (m), Summer','fontweight','bold','FontSize',15)
xlabel('Longitude','fontweight','bold','FontSize',12)
ylabel('Latitude','fontweight','bold','FontSize',12)
axis equal
axis([173.6 179.4 -35 -28])