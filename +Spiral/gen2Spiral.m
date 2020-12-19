function gen2Spiral()
%%this version can generate polygon having more than 4 sides.
% mex CdownSamplingUsingRealCoords.cpp
% mex CgetBoundariesFromLabel2.cpp
%% parameter setting
%   current unit: mm
f=50; % focal length
incidentAngle=0; %  illumination angle
lambda=1.75e-6; % wavelength
R=200e-3;   % radius
delta=R;    % separation of the two spirals
Nx=500;Ny=500;  % number of samplings for each sub-region
ringSampling=50;     % number of samplings for each zone
ringSamplingNum=Nx/ringSampling;
DownSamplingAccuCtrl=0.001; % downsampling accuracy
CoordsAccuCtrl=0.0001;  % boundary searching accuracy
filenamestr=['wave',num2str(lambda),'_d',num2str(delta)];   %   file name
filename=['C:\Users\T470P\Desktop\',filenamestr,'.gds'];    %   save path + file name
outputFile=OpenGDS(filename,filenamestr);
incidentAngle=incidentAngle/180*pi;
RingNum=(sqrt(f^2+(R+delta/2)^2)-f)/lambda-(sqrt(f^2+(delta/2)^2)-f)/lambda;    % number of zones
divnum=2*floor(RingNum/ringSamplingNum)+1;   %number of sub-region
dx=2*R/divnum;
dy=2*R/divnum;
Mx=divnum;
My=divnum;
fliped=divnum/2+1;
th=8;   %  threshold for generating hologram
%%  set unit to anstrom
db=10000000;
lambda=lambda*db;
delta=delta*db;
f=f*db;
R=R*db;
Rb=0*db;
dx=db*dx;
dy=db*dy;
offsetx=0*db;offsety=0*db;  %    offset of final coordinates
%%  generate initial boundary coordinates based coordinate transform
tic
ceout=cell(Mx,My);
parfor p=1:Mx
    for q=1:My
        ceout{p,q}=getCoords(lambda,delta,f,R,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th);
    end
end
toc

%% connect the shapes cut by dicing
Ixy=[];
uxy=[];
Ixy0=[];
sns=0;
for p=1:Mx
    for q=1:My
        Ixy0=[Ixy0,ceout{p,q}.Ixy0];
        Ixy=[Ixy,ceout{p,q}.Ixy];
        if ~isempty(ceout{p,q}.uxy)
            ceout{p,q}.uxy(:,4)=ceout{p,q}.uxy(:,4)+sns;
        end
        uxy=[uxy;ceout{p,q}.uxy];
        sns=length(Ixy);
    end
end

while ~isempty(uxy)
    t=uxy(1,4);
    temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)~=t);
    if isempty(temp)
        temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)==t);
        if isempty(temp)
            break;
        end
        uxy(temp,:)=[];
        temp3=find(uxy(:,4)==t);
        if isempty(temp3)
            Ixy0(end+1)=Ixy(t);
        end
        continue;
    end
    k=uxy(temp,4);
    Ixy(t).xr=[Ixy(t).xr;Ixy(k).xr];
    Ixy(t).yr=[Ixy(t).yr;Ixy(k).yr];
    uxy(temp,:)=[];
    uxy(1,:)=[];
    temp2=find(uxy(:,4)==k);
    if ~isempty(temp2)
        uxy(temp2,4)=t;
    end
    temp3=find(uxy(:,4)==t);
    if isempty(temp3)
        Ixy0(end+1)=Ixy(t);
    end
end

%% remove obscuration
for t=length(Ixy0):-1:1
    if ~all((Ixy0(t).xr).^2+(Ixy0(t).yr).^2>Rb.^2&(Ixy0(t).xr).^2+(Ixy0(t).yr).^2<=(R).^2)
        Ixy0(t)=[];
    end
end
%% generate GDS from each shape
tic
len=length(Ixy0);
for q=1:len
    cxy=[Ixy0(q).xr,Ixy0(q).yr];
    cxy=unique(cxy,'rows');
    cn=length(cxy);
    k=k+1;
    cx=cxy(:,1);
    cy=cxy(:,2);
    % find accurate coordinates and sort them
    B=CgetBoundariesFromLabel2(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
    % downsampling the coordinates
    Bs=CdownSamplingUsingRealCoords(B,lambda,delta,f,DownSamplingAccuCtrl);
    %% remove bad shapes
    for si=length(Bs):-1:1
        while 1
            num=length(Bs{si});
            angle=zeros(num,1);
            for js=1:num
                if js==1
                    angle(1)=calAngle(Bs{si}(end,1),Bs{si}(end,2),Bs{si}(1,1),Bs{si}(1,2),Bs{si}(2,1),Bs{si}(2,2));
                elseif js==num
                    angle(num)=calAngle(Bs{si}(num-1,1),Bs{si}(num-1,2),Bs{si}(num,1),Bs{si}(num,2),Bs{si}(1,1),Bs{si}(1,2));
                else
                    angle(js)=calAngle(Bs{si}(js-1,1),Bs{si}(js-1,2),Bs{si}(js,1),Bs{si}(js,2),Bs{si}(js+1,1),Bs{si}(js+1,2));
                end
            end
            da=diff([angle;angle(1)]);
            index=abs(da)>4;
            if sum(abs(angle)>3.1)>1
                index=index|abs(angle)>3.1;
            end
            Bs{si}(index,:)=[];
            [~,rss]=cart2pol(Bs{si}(:,1),Bs{si}(:,2));
            %control minimum unit
            %             if max(rss)-min(rss)<20000
            %                 Bs{si}=[];
            %             end
            if sum(index)==0||length(Bs{si})<3
                break;
            end
        end
    end
    %% generate boundary
    for ns=1:length(Bs)
        xy=Bs{ns};
        if isempty(xy)||size(xy,1)<3
            continue;
        end
        if (xy(1,1)^2+xy(1,2)^2)>R^2||(xy(1,1)^2+xy(1,2)^2)<Rb^2
            continue;
        end
        xy(end+1,:)=xy(1,:);
        xy(:,1)=xy(:,1)+offsetx;
        xy(:,2)=xy(:,2)+offsety;
        Np=size(xy,1)-1;%多边形顶点数
        CreateBoundary(outputFile,xy',Np);
    end
end
toc
%% finish GDS
CloseGDS(outputFile);
%% open generated GDS file
winopen(filename);


function stout=getCoords(lambda,delta,f,R,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th)
sx=(p-fliped)*dx;
sy=(q-fliped)*dy;
Rx=sx;
Rx(2)=Rx+dx;
Ry=sy;
Ry(2)=Ry+dy;
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y);
% rs=sqrt(r*R);
hd=delta/2/3;
rs=sqrt((r+hd)*(R+2*hd+2*sqrt(hd*(R+hd))))-sqrt(hd*(R+2*hd+2*sqrt(hd*(R+hd))));
% rs=sqrt((r+delta/2)*(R+delta+2*sqrt(delta/2*(R+delta/2))))-sqrt((delta/2)*(R+delta+2*sqrt(delta/2*(R+delta/2))));
xs=rs.*cos(TH);
ys=rs.*sin(TH);
stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th);

function stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th)
Ixy=[];
Ixy0=[];
uxy=[];
sns=1;
I=getIntensity(xs,ys,delta,f,lambda,incidentAngle);
[sr,sc]=size(xs);
I0=I;
deltas=0.1;
I(I<th-deltas)=0;
I(I>=th-deltas)=1;
[L,num]=bwlabel(I);
L0=L;
L0(I0>th+deltas)=0;
for i=1:num
    %up
    l=size(uxy,1);
    temp=find(L(1,:)==i);
    if ~isempty(temp)&&q~=1
        uxy(end+1,1:4)=[p,q-0.5,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p,q-0.5,temp(df),sns];
        end
    end
    %down
    temp=find(L(end,:)==i);
    if ~isempty(temp)&&q~=My
        uxy(end+1,1:4)=[p,q+0.5,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p,q+0.5,temp(df),sns];
        end
    end
    %left
    temp=find(L(:,1)==i);
    if ~isempty(temp)&&p~=1
        uxy(end+1,1:4)=[p-0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p-0.5,q,temp(df),sns];
        end
    end
    %right
    temp=find(L(:,end)==i);
    if ~isempty(temp)&&p~=Mx
        uxy(end+1,1:4)=[p+0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p+0.5,q,temp(df),sns];
        end
    end
    if l==size(uxy,1)
        [y,x]=find(L0==i);
        Ixy0(end+1).xr=xs((x-1)*sr+y);
        Ixy0(end).yr=ys((x-1)*sr+y);
    else
        sns=sns+1;
        [y,x]=find(L0==i);
        Ixy(end+1).xr=xs((x-1)*sr+y);
        Ixy(end).yr=ys((x-1)*sr+y);
    end
end
stout.Ixy=Ixy;
stout.Ixy0=Ixy0;
stout.uxy=uxy;
stout.sns=sns;

function angle=calAngle(x0,y0,x1,y1,x2,y2)
A=[x1-x0,y1-y0];
B=[x2-x1,y2-y1];
y=det([A;B]);
x=dot(A,B);
angle=atan2(y,x);

function I=getIntensity(xs,ys,delta,f,lambda,incidentAngle)
PI=3.14159265;
T=lambda/sin(atan(delta/f/2));
offset=f*tan(incidentAngle);
th1=atan2(ys,xs-delta/2);
th2=-atan2(ys,xs+delta/2);
I=(cos(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2)+th1)+...
    cos(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2)+th2)+...
    2*cos(2*PI*ys/lambda*sin(incidentAngle))).^2+...
    (-sin(2*PI/lambda*sqrt((xs-delta/2).^2+(ys-offset).^2+f.^2)+th1)+...
    -sin(2*PI/lambda*sqrt((xs+delta/2).^2+(ys-offset).^2+f.^2)+th2)+...
    2*sin(2*PI*ys/lambda*sin(incidentAngle))).^2;