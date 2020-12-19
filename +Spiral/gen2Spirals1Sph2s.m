tic
len=length(Ixy0);
ds=floor(len/99);
file=zeros(1,100);
for p=1:100
    file(p)=fopen(['temp',num2str(p),'dat'],'wb');
    if p~=100
        for q=(p-1)*ds+1:p*ds
            cxy=[Ixy0(q).xr,Ixy0(q).yr];
            cxy=unique(cxy,'rows');
            cn=length(cxy);
            k=k+1;
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
                Spiral.CreateBoundary(file(p),xy',Np);
            end
        end
    else
        for q=99*ds+1:len
            cxy=[Ixy0(q).xr,Ixy0(q).yr];
            cxy=unique(cxy,'rows');
            cn=length(cxy);
            k=k+1;
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
                Spiral.CreateBoundary(file(p),xy',Np);
            end
        end
    end
   fclose(file(p)); 
end
toc
%% finish GDS
outputFile=Spiral.OpenGDS(filename,filenamestr);
for p=1:100
    file(p)=fopen(['temp',num2str(p),'dat'],'rb');
    coords=fread(file(p),'int32','b');
    fclose(file(p));
    fwrite(outputFile,coords, 'int32','b');
end
Spiral.CloseGDS(outputFile);
delete('*.dat');