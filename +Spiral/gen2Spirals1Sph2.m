function gen2Spirals1Sph2()
%%this version can generate polygon having more than 4 sides.
% mex CdownSamplingUsingRealCoords.cpp
% mex CgetBoundariesFromLabel2.cpp
%% parameter setting
%   current unit: mm
f=14.67; % focal length
incidentAngle=0; %  illumination angle
lambda=1.75e-4; % wavelength
R=165e-3;   % radius
delta=2*f/15*216e-3;    % separation of the two spirals
Nx=500;Ny=500;  % number of samplings for each sub-region
ringSampling=50;     % number of samplings for each zone
ringSamplingNum=Nx/ringSampling;
DownSamplingAccuCtrl=0.001; % downsampling accuracy
CoordsAccuCtrl=0.0001;  % boundary searching accuracy
filenamestr=['wave',num2str(lambda),'_d',num2str(delta)];   %   file name
filename=['D:\',filenamestr,'.gds'];    %   save path + file name
incidentAngle=incidentAngle/180*pi;
RingNum=(sqrt(f^2+(R+delta/2)^2)-f)/lambda-(sqrt(f^2+(delta/2)^2)-f)/lambda;    % number of zones
divnum=2*floor(RingNum/ringSamplingNum)+1;   %number of sub-region
dx=2*R/divnum;
dy=2*R/divnum;
Mx=divnum;
My=divnum;
fliped=divnum/2+1;
th=18;   %  threshold for generating hologram
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

%% extract coordinates from ceout(this part is slow)
tic
Ixy=[];% save diced boundary shapes
uxy=[];% used for connecting the diced boundary shapes
Ixy0=[];% save other shapes
sns=0; % number of each diced boundary shape
for p=1:Mx
    for q=1:My
        Ixy0=[Ixy0,ceout{p,q}.Ixy0];
        Ixy=[Ixy,ceout{p,q}.Ixy];
        if ~isempty(ceout{p,q}.uxy) % modify the number of the diced boundary shapes to global number
            ceout{p,q}.uxy(:,4)=ceout{p,q}.uxy(:,4)+sns;
        end
        uxy=[uxy;ceout{p,q}.uxy];
        sns=length(Ixy);
    end
end
toc
%% connect the boundary shapes by uxy (this part is slow)
tic
count=[];
while ~isempty(uxy)
    t=uxy(1,4);
    temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)~=t); % find boundary shapes belong to different sub-region
    if isempty(temp)
        temp=find(uxy(:,1)==uxy(1,1)&uxy(:,2)==uxy(1,2)&uxy(:,3)==uxy(1,3)&uxy(:,4)==t);% find error boundary shapes
        if isempty(temp)
            break;
        end
        uxy(temp,:)=[];% remove error boundary shapes
        temp3=find(uxy(:,4)==t);% make sure there is no shape marked number t
        if isempty(temp3)
            %             Ixy0(end+1)=Ixy(t);% transfer these shapes to Ixy0
            count(end+1)=t;
        end
        continue;
    end
    k=uxy(temp,4);
    Ixy(t).xr=[Ixy(t).xr;Ixy(k).xr];% merge same boundary shapes
    Ixy(t).yr=[Ixy(t).yr;Ixy(k).yr];
    uxy(temp,:)=[];
    uxy(1,:)=[];
    temp2=find(uxy(:,4)==k);
    if ~isempty(temp2)
        uxy(temp2,4)=t;% modify the number of the merged boundary shapes
    end
    temp3=find(uxy(:,4)==t);% make sure there is no shape marked number t
    if isempty(temp3)
        %         Ixy0(end+1)=Ixy(t);% transfer these shapes to Ixy0
        count(end+1)=t;
    end
end
Ixy0=[Ixy0,Ixy(count)];
toc
%% remove obscuration
len=length(Ixy0);
ts=zeros(1,len);
parfor t=1:len
    if ~all((Ixy0(t).xr).^2+(Ixy0(t).yr).^2>Rb.^2&(Ixy0(t).xr).^2+(Ixy0(t).yr).^2<=(R).^2)
        ts(t)=1;
    end
end
Ixy0(ts==1)=[];
%% generate GDS from each shape
tic
len=length(Ixy0); % number of total jobs
ds=floor(len/99); % number of each job
file=zeros(1,100);
for sp=0:19
    for p=sp*5+1:sp*5+5
        file(p)=fopen(['temp',num2str(p),'.dat'],'wb'); % create temporary file for save the processed coordinates (in order to perform CPU accelaration)
        if p~=100
            for q=(p-1)*ds+1:p*ds
                cxy=[Ixy0(q).xr,Ixy0(q).yr]; % coordinates
                cxy=unique(cxy,'rows');% remove redundant coordinates
                cn=length(cxy);
                cx=cxy(:,1);
                cy=cxy(:,2);
                % find accurate coordinates and sort them
                B=Spiral.CgetBoundariesFromLabel3(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
                % downsampling the coordinates
                Bs=Spiral.CdownSamplingUsingRealCoords(B,lambda,delta,f,DownSamplingAccuCtrl);
                %% remove bad shapes
                for si=length(Bs):-1:1
                    while 1
                        num=length(Bs{si});
                        angle=zeros(num,1);
                        for js=1:num % calculate the angle of three adjcent points
                            if js==1
                                angle(1)=Spiral.calAngle(Bs{si}(end,1),Bs{si}(end,2),Bs{si}(1,1),Bs{si}(1,2),Bs{si}(2,1),Bs{si}(2,2));
                            elseif js==num
                                angle(num)=Spiral.calAngle(Bs{si}(num-1,1),Bs{si}(num-1,2),Bs{si}(num,1),Bs{si}(num,2),Bs{si}(1,1),Bs{si}(1,2));
                            else
                                angle(js)=Spiral.calAngle(Bs{si}(js-1,1),Bs{si}(js-1,2),Bs{si}(js,1),Bs{si}(js,2),Bs{si}(js+1,1),Bs{si}(js+1,2));
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
                    if isempty(xy)||size(xy,1)<3 % remove line shapes and dot shapes
                        continue;
                    end
                    if (xy(1,1)^2+xy(1,2)^2)>R^2||(xy(1,1)^2+xy(1,2)^2)<Rb^2 % remove coordinates belong to the obscuration or over the radius
                        continue;
                    end
                    xy(end+1,:)=xy(1,:);
                    xy(:,1)=xy(:,1)+offsetx;
                    xy(:,2)=xy(:,2)+offsety;
                    Np=size(xy,1)-1;%多边形顶点数
                    Spiral.CreateBoundary(file(p),xy',Np); % generate gds shapes
                end
            end
        else % case p==100
            for q=99*ds+1:len
                cxy=[Ixy0(q).xr,Ixy0(q).yr];
                cxy=unique(cxy,'rows');
                cn=length(cxy);
                cx=cxy(:,1);
                cy=cxy(:,2);
                % find accurate coordinates and sort them
                B=Spiral.CgetBoundariesFromLabel3(cx,cy,cn,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl);
                % downsampling the coordinates
                Bs=Spiral.CdownSamplingUsingRealCoords(B,lambda,delta,f,DownSamplingAccuCtrl);
                %% remove bad shapes
                for si=length(Bs):-1:1
                    while 1
                        num=length(Bs{si});
                        angle=zeros(num,1);
                        for js=1:num
                            if js==1
                                angle(1)=Spiral.calAngle(Bs{si}(end,1),Bs{si}(end,2),Bs{si}(1,1),Bs{si}(1,2),Bs{si}(2,1),Bs{si}(2,2));
                            elseif js==num
                                angle(num)=Spiral.calAngle(Bs{si}(num-1,1),Bs{si}(num-1,2),Bs{si}(num,1),Bs{si}(num,2),Bs{si}(1,1),Bs{si}(1,2));
                            else
                                angle(js)=Spiral.calAngle(Bs{si}(js-1,1),Bs{si}(js-1,2),Bs{si}(js,1),Bs{si}(js,2),Bs{si}(js+1,1),Bs{si}(js+1,2));
                            end
                        end
                        da=diff([angle;angle(1)]);
                        index=abs(da)>4;
                        if sum(abs(angle)>3.1)>1
                            index=index|abs(angle)>3.1;
                        end
                        Bs{si}(index,:)=[];
                        [~,rss]=cart2pol(Bs{si}(:,1),Bs{si}(:,2));
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
                    Spiral.CreateBoundary(file(p),xy',Np);
                end
            end
        end
        fclose(file(p));
    end
end
toc
%% finish GDS
outputFile=Spiral.OpenGDS(filename,filenamestr);% initial gds
% transfer coordinates from the temporary files to gds file
for p=1:100
    file(p)=fopen(['temp',num2str(p),'.dat'],'rb');
    coords=fread(file(p),'int32','b');
    fclose(file(p));
    fwrite(outputFile,coords, 'int32','b');
end
Spiral.CloseGDS(outputFile);% finish gds
delete('*.dat');% delete temporay files
%% open generated GDS file
winopen(filename);


function stout=getCoords(lambda,delta,f,R,Nx,Ny,p,q,Mx,My,dx,dy,fliped,incidentAngle,th)
% generate sub-region coordinates
sx=(p-fliped)*dx;
sy=(q-fliped)*dy;
Rx=sx;
Rx(2)=Rx+dx;
Ry=sy;
Ry(2)=Ry+dy;
[x,y]=meshgrid(linspace(Rx(1),Rx(2),Nx*2+1),linspace(Ry(1),Ry(2),Ny*2+1));
[TH,r] = cart2pol(x,y);
% new coordinates transform
hd=delta/2/3;
rs=sqrt((r+hd)*(R+2*hd+2*sqrt(hd*(R+hd))))-sqrt(hd*(R+2*hd+2*sqrt(hd*(R+hd))));
xs=rs.*cos(TH);
ys=rs.*sin(TH);
stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th);

function stout=getBoundaries(xs,ys,delta,f,lambda,incidentAngle,p,q,Mx,My,th)
Ixy=[];% boundary shape
Ixy0=[];% normal shape
uxy=[];% used for connecting the diced boundary shapes
sns=1;% number of boundary shape
I=getIntensity(xs,ys,delta,f,lambda,incidentAngle); %generate hologram
[sr,sc]=size(xs);
I0=I;
deltas=0.1;
% generate rough zoneplate boundary coordinates
I(I<th-deltas)=0;
I(I>=th-deltas)=1;
[L,num]=bwlabel(I);
L0=L;
L0(I0>th+deltas)=0;
% distinguish and save boundary shapes
for i=1:num
    % top boundary
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
    % buttom boundary
    temp=find(L(end,:)==i);
    if ~isempty(temp)&&q~=My
        uxy(end+1,1:4)=[p,q+0.5,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p,q+0.5,temp(df),sns];
        end
    end
    %left boundary
    temp=find(L(:,1)==i);
    if ~isempty(temp)&&p~=1
        uxy(end+1,1:4)=[p-0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p-0.5,q,temp(df),sns];
        end
    end
    %right boundary
    temp=find(L(:,end)==i);
    if ~isempty(temp)&&p~=Mx
        uxy(end+1,1:4)=[p+0.5,q,temp(1),sns];
        if (temp(end)-temp(1))~=length(temp)-1
            df=find(diff(temp)~=1);
            df=df(1)+1;
            uxy(end+1,1:4)=[p+0.5,q,temp(df),sns];
        end
    end
    if l==size(uxy,1) % save other shapes
        [y,x]=find(L0==i);
        Ixy0(end+1).xr=xs((x-1)*sr+y);
        Ixy0(end).yr=ys((x-1)*sr+y);
    else %save boundary shapes
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
% calculate the angle between two vectors
A=[x1-x0,y1-y0];
B=[x2-x1,y2-y1];
y=det([A;B]);
x=dot(A,B);
angle=atan2(y,x);

function I=getIntensity(xs,ys,delta,f,lambda,incidentAngle)
% get the intensity of the hologram at (xs,ys)
PI=3.14159265;
T=lambda/sin(atan(delta/f/2));
offset=f*tan(incidentAngle);
th1=atan2(ys-delta/2*cos(pi/3),xs-delta/2*sin(pi/3));
th2=-atan2(ys-delta/2*cos(pi/3),xs+delta/2*sin(pi/3));
% function of the hologram
I=(cos(2*PI/lambda*sqrt((xs-delta/2*sin(pi/3)).^2+(ys-delta/2*cos(pi/3)).^2+f.^2)+th1)+...
    cos(2*PI/lambda*sqrt((xs+delta/2*sin(pi/3)).^2+(ys-delta/2*cos(pi/3)).^2+f.^2)+th2)+...
    cos(2*PI/lambda*sqrt((xs).^2+(ys+delta/2).^2+f.^2))+...
    3*cos(2*PI*ys/lambda*sin(incidentAngle))).^2+...
    (-sin(2*PI/lambda*sqrt((xs-delta/2*sin(pi/3)).^2+(ys-delta/2*cos(pi/3)).^2+f.^2)+th1)+...
    -sin(2*PI/lambda*sqrt((xs+delta/2*sin(pi/3)).^2+(ys-delta/2*cos(pi/3)).^2+f.^2)+th2)+...
    -sin(2*PI/lambda*sqrt((xs).^2+(ys+delta/2).^2+f.^2))+...
    3*sin(2*PI*ys/lambda*sin(incidentAngle))).^2;