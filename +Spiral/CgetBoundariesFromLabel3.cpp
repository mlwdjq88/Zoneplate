#include <fstream>
#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>
#include <cstring>
#define PI 3.14159265

// get the maximum value
double dmax(double * A, long Ibegin, long Iend){
    double temp=A[Ibegin];
    for(long i=Ibegin+1; i<Iend; i++){
        if (temp<A[i]){
            temp=A[i];
        }
    }
    return temp;
}

// get the minimum value
double dmin(double * A, long Ibegin, long Iend){
    double temp=A[Ibegin];
    for(long i=Ibegin+1; i<Iend; i++){
        if (temp>A[i]){
            temp=A[i];
        }
    }
    return temp;
}

// find the specific number
long *dfind(double * A, long num,double B,long &Ns){
    long *Index;
    Index= (long*) malloc(num * sizeof(long));
    long k=0;
    for(long i=0; i<num; i++){
        if (A[i]==B){
            Index[k]=i;
            k++;
        }
    }
    Ns=k;
    return Index;
}
// from cart to polar coordinates
void cart2pol(double *xr,double *yr,double *ths,double *rs,long num,int angM){
    double maxN=8,maxP=8,minN=8,minP=8;
    for(long i=0;i<num;i++){
        ths[i]=atan2(yr[i],xr[i]);
        rs[i]=sqrt(pow(yr[i],2)+pow(xr[i],2));
        if (ths[i]>=0){
            if (maxP==8||ths[i]>maxP)
                maxP=ths[i];
            if (minP==8||ths[i]<minP)
                minP=ths[i];
        }
        else{
            if (maxN==8||ths[i]>maxN)
                maxN=ths[i];
            if (minN==8||ths[i]<minN)
                minN=ths[i];
        }
    }
    if (angM==1||(maxN!=8&&maxP!=8&&(2*PI+minN-maxP< minP-maxN))){
        for(long i=0;i<num;i++){
            if (ths[i]<0)
                ths[i]=ths[i]+2*PI;
        }
    }
}
// from cart to polar coordinates
void cart2polCir(double *xr,double *yr,double *ths,double *rs,long num,double *tr,int &angM){
    double maxN=8,maxP=8,minN=8,minP=8,r1=0,r2=0,r3=0,r4=0;
    for(long i=0;i<num;i++){
        ths[i]=atan2(yr[i],xr[i]);
        rs[i]=sqrt(pow(yr[i],2)+pow(xr[i],2));
        if (ths[i]>=0){
            if (maxP==8||ths[i]>maxP){
                maxP=ths[i];
                r3=rs[i];
            }
            if (minP==8||ths[i]<minP){
                minP=ths[i];
                r1=rs[i];
            }
        }
        else{
            if (maxN==8||ths[i]>maxN){
                maxN=ths[i];
                r4=rs[i];
            }
            if (minN==8||ths[i]<minN){
                minN=ths[i];
                r2=rs[i];
            }
        }
    }
    if (maxN!=8&&maxP!=8&&(2*PI+minN-maxP< minP-maxN)){
        angM=1;
        for(long i=0;i<num;i++){
            if (ths[i]<0)
                ths[i]=ths[i]+2*PI;
        }
        if ((minP-maxN)<20*PI/180){
            tr[0]=minP;
            tr[1]=r1;
            tr[2]=maxN;
            tr[3]=r4;
            
        }
    }
    else if(maxN!=8&&maxP!=8&&(2*PI+minN-maxP)<20*PI/180){
        tr[0]=minN;
        tr[1]=r2;
        tr[2]=maxP;
        tr[3]=r3;
    }
}
// from cart to polar coordinates (standard)
void cart2polStd(double *xr,double *yr,double *ths,double *rs,long num){
    for(long i=0;i<num;i++){
        ths[i]=atan2(yr[i],xr[i]);
        rs[i]=sqrt(pow(yr[i],2)+pow(xr[i],2));
    }
}
// from cart to polar coordinates (single point)
void cart2pol1(double xr,double yr,double &ths,double &rs){
    ths=atan2(yr,xr);
    rs=sqrt(pow(yr,2)+pow(xr,2));
}
// from  polar to cart coordinates
void pol2cart(double *ths,double *rs,double *xr,double *yr,long num){
    for(long i=0;i<num;i++){
        xr[i]=rs[i]*cos(ths[i]);
        yr[i]=rs[i]*sin(ths[i]);
    }
}
// from  polar to cart coordinates (singles point)
void pol2cart1(double ths,double rs,double &xr,double &yr){
    xr=rs*cos(ths);
    yr=rs*sin(ths);
}
// calculate the average value
double mean(double *A,long num){
    double s=0;
    for(long i=0;i<num;i++){
        s+=A[i];
    }
    s=s/(double)num;
    return s;
}

// remove NULL 
long removeNULL(double *x,long num){
    long k=0;
    if(num>10000){
        double *y;
        y= (double*) malloc(num * sizeof(double));
        for(int i=0;i<num;i++){
            if (x[i]!=NULL){
                y[k]=x[i];
                k++;
            }
        }
        for(int i=0;i<k;i++){
            x[i]=y[i];
        }
        free(y);
    }
    else {
        double y[10000];
        for(int i=0;i<num;i++){
            if (x[i]!=NULL){
                y[k]=x[i];
                k++;
            }
        }
        for(int i=0;i<k;i++){
            x[i]=y[i];
        }
    }
    return k;
}
// get the intensity of the hologram at (xs,ys)
double getIntensity(double xs,double ys,double delta,double f,double lambda,double incidentAngle){
    double I=0,delta1,delta2,T,offset;
    T=lambda/sin(atan(delta/f/2));
    offset=f*tan(incidentAngle);
    delta2=f*(tan(asin(sin(incidentAngle)+lambda/T)))-offset;
    delta1=offset-f*(tan(asin(sin(incidentAngle)-lambda/T)));
    //hologram function 
    I=pow(cos(2*PI/lambda*sqrt(pow(xs-delta/2*sin(PI/3),2)+pow(ys-delta/2*cos(PI/3),2)+pow(f,2))+atan2(ys-delta/2*cos(PI/3),xs-delta/2*sin(PI/3)))+cos(2*PI/lambda*sqrt(pow(xs+delta/2*sin(PI/3),2)+pow(ys-delta/2*cos(PI/3),2)+pow(f,2))-atan2(ys-delta/2*cos(PI/3),xs+delta/2*sin(PI/3)))+cos(2*PI/lambda*sqrt(pow(xs,2)+pow(ys+delta/2,2)+pow(f,2)))+3*cos(2*PI/lambda*ys*sin(incidentAngle)),2)\
           +pow(-sin(2*PI/lambda*sqrt(pow(xs-delta/2*sin(PI/3),2)+pow(ys-delta/2*cos(PI/3),2)+pow(f,2))+atan2(ys-delta/2*cos(PI/3),xs-delta/2*sin(PI/3)))-sin(2*PI/lambda*sqrt(pow(xs+delta/2*sin(PI/3),2)+pow(ys-delta/2*cos(PI/3),2)+pow(f,2))-atan2(ys-delta/2*cos(PI/3),xs+delta/2*sin(PI/3)))-sin(2*PI/lambda*sqrt(pow(xs,2)+pow(ys+delta/2,2)+pow(f,2)))+3*sin(2*PI/lambda*ys*sin(incidentAngle)),2);
    return I;
}

// get the accurate coordinates based on the rough coordinates
void getAccurateCoords(double *xr,double *yr,double *ths,double *rs,double x0,double y0,double delta,double f,double lambda,double th,double incidentAngle,double CoordsAccuCtrl,long &Ns){
    double *thi, *ri,rsmin,rsmax,rimin,rimax,*rs0,*r0,Ia,Ib,a,b,xs=0,ys=0,Im,s,*xr0,*yr0;
    long Ns0=Ns;
    int k;
    thi= (double*) malloc(Ns * sizeof(double));
    ri= (double*) malloc(Ns * sizeof(double));
    rs0= (double*) malloc(Ns * sizeof(double));
    r0= (double*) malloc(Ns * sizeof(double));
    xr0= (double*) malloc(Ns * sizeof(double));
    yr0= (double*) malloc(Ns * sizeof(double));
    for (long i=0;i<Ns;i++){
        xr0[i]=xr[i]-x0;
        yr0[i]=yr[i]-y0;
    }
    cart2polStd(xr0,yr0,thi,ri,Ns);
    for (long i=0;i<Ns;i++){
        xr0[i]=xr[i];
        yr0[i]=yr[i];
    }
    rsmin=dmin(rs,0,Ns);
    rsmax=dmax(rs,0,Ns);
    rimin=dmin(ri,0,Ns);
    rimax=dmax(ri,0,Ns);
    for (long i=0;i<Ns;i++){
        r0[i]=ri[i];// save original radius
    }
    // get the accurate coordinates based on the rough coordinates
    for (long i=0;i<Ns;i++){
        Ia=getIntensity(xr[i],yr[i],delta,f,lambda,incidentAngle)-th; // get inital instensity error
        if (fabs(Ia)<0.1){ // Judge whether the point belong to boundary
            a=ri[i]; // initial radius a
            b=a*(Ia/2+th)/th; // find a new radius b 
            pol2cart1(thi[i],b,xs,ys);
            Ib=getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th;//get new intensity error
            if (fabs(Ia)>fabs(Ib)) // whether the new error is smaller
                ri[i]=b;
            Im=fmin(fabs(Ia),fabs(Ib));
            k=0;
            // iterator for finding more accurate coordiantes minimizing the intensity error
            while(Im>CoordsAccuCtrl&&fabs(Ia-Ib)>0.0000001&&k<100){ 
                if (fabs(Ia)>fabs(Ib)){
                    s=(b-a)/(Ib-Ia)*Ib;
                    a=b-s;
                    pol2cart1(thi[i],a,xs,ys);
                    while(fabs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th)>fabs(Ia)&&a!=b){
                        s=s/2;
                        a=b-s;
                        pol2cart1(thi[i],a,xs,ys);
                    }
                    Ia=fabs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th);
                    ri[i]=a;
                }
                else{
                    s=(b-a)/(Ib-Ia)*Ia;
                    b=a-s;
                    pol2cart1(thi[i],b,xs,ys);
                    while(fabs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th)>fabs(Ib)&&a!=b){
                        s=s/2;
                        b=a-s;
                        pol2cart1(thi[i],b,xs,ys);
                    }
                    Ib=fabs(getIntensity(xs+x0,ys+y0,delta,f,lambda,incidentAngle)-th);
                    ri[i]=b;
                }
                Im=fmin(fabs(Ia),fabs(Ib));
                k++;
            }
        }
    }
    pol2cart(thi,ri,xr,yr,Ns);// back to cart 
    for (long i=0;i<Ns;i++){
        xr[i]=xr[i]+x0;
        yr[i]=yr[i]+y0;
        rs0[i]=sqrt(pow(xr[i],2)+pow(yr[i],2));
    }
    //remove worse accurate coordinates
    for (long i=0;i<Ns;i++){
        if (ri[i]>1.1*rimax||ri[i]<0.9*rimin||fabs(ri[i]-r0[i])/(rimax-rimin)>0.05||fabs(ri[i]-r0[i])/r0[i]>0.25||\
                (rs0[i]>(1.1*rsmax-0.1*rsmin))||(rs0[i]<(1.1*rsmin-0.1*rsmax))||\
                (fabs(getIntensity((xr[i]+xr0[i])/2.0,(yr[i]+yr0[i])/2.0,delta,f,lambda,incidentAngle)-th)>fabs(getIntensity(xr0[i],yr0[i],delta,f,lambda,incidentAngle)-th))){
            xr[i]=NULL;
            yr[i]=NULL;
        }
    }
    Ns=removeNULL(xr,Ns0);
    Ns=removeNULL(yr,Ns0);
    
    free(thi);
    free(ri);
    free(rs0);
    free(r0);
    free(xr0);
    free(yr0);
}

// sort coordinates
void sortCoods(double *xr,double *yr,double *ths,double *rs,long Ns){
    double thmin,thmax,rmin,rmax,ravg,*Sths,*Srs,sides,savgN=0,savgP=0;
    long *Imin,*Imax,k=0,dk=Ns-1,NI;
    thmin=dmin(ths,0,Ns);
    Imin=dfind(ths,Ns,thmin,NI);
    rmin=rs[Imin[0]];
    thmax=dmax(ths,0,Ns);
    Imax=dfind(ths,Ns,thmax,NI);
    rmax=rs[Imax[0]];
    ravg=mean(rs,Ns);
    Sths= (double*) malloc(Ns * sizeof(double));
    Srs= (double*) malloc(Ns * sizeof(double));
    if ((rmin-ravg)*(rmax-ravg)>0&&\
            ((fabs(rmin-ravg)/(ravg-dmin(rs,0,Ns))>1/3.0&&fabs(rmax-ravg)/(ravg-dmin(rs,0,Ns))>1/3.0)||\
            (fabs(rmin-ravg)/(dmax(rs,0,Ns)-ravg)>1/3.0&&fabs(rmax-ravg)/(dmax(rs,0,Ns)-ravg)>1/3.0)))
        rmin=ravg;
    // sort coordinates in the polar space 
    for (long i=0;i<Ns;i++){
        sides=rs[i]*(thmax-thmin)-ths[i]*(rmax-rmin)-rmin*thmax+rmax*thmin;
        if (sides<0){
            savgN=savgN+sides;
            Sths[k]=ths[i];
            Srs[k]=rs[i];
            if (k>0){
                for (long j=0;j<k;j++){
                    if (ths[i]<Sths[j]){
                        for (long p=k;p>j;p--){
                            Sths[p]=Sths[p-1];
                            Srs[p]=Srs[p-1];
                        }
                        Sths[j]=ths[i];
                        Srs[j]=rs[i];
                        break;
                    }
                }
                
            }
            k++;
        }
        else{
            savgP=savgP+sides;
            Sths[dk]=ths[i];
            Srs[dk]=rs[i];
            if (dk<Ns-1){
                for (long j=Ns-1;j>dk;j--){
                    if (ths[i]<Sths[j]){
                        for (long p=dk;p<j;p++){
                            Sths[p]=Sths[p+1];
                            Srs[p]=Srs[p+1];
                        }
                        Sths[j]=ths[i];
                        Srs[j]=rs[i];
                        break;
                    }
                }
                
            }
            dk--;
        }
    }
    if(k!=0)
        savgN=savgN/(k*1.0);
    if (k!=Ns-2)
        savgP=savgP/((Ns-2-k)*1.0);
    // if above sorting method is not good, then resorting the coordinates based on the gravity center 
    if(k<=1||k>=Ns-3||savgP==0||savgN/savgP<0.08||savgN/savgP>1/0.08){ 
        double x0=mean(xr,Ns);
        double y0=mean(yr,Ns);
        for(long i=0;i<Ns;i++){
            Sths[i]=atan2(yr[i]-y0,xr[i]-x0);
        }
        double temp;
        for(long i=0;i<Ns-1;i++){
            for(long j=0;j<Ns-1-i;j++){
                if(Sths[j]>Sths[j+1]){
                    temp=Sths[j];
                    Sths[j]=Sths[j+1];
                    Sths[j+1]=temp;
                    temp=xr[j];
                    xr[j]=xr[j+1];
                    xr[j+1]=temp;
                    temp=yr[j];
                    yr[j]=yr[j+1];
                    yr[j+1]=temp;
                }
            }
        }
    }
    else{
        pol2cart(Sths,Srs,xr,yr,Ns);
    }
    free(Sths);
    free(Srs);
    free(Imax);
    free(Imin);
}
// main function 
void getBoundariesFromLabel(double *xs,double *ys,double delta,double f,double lambda,double th,double incidentAngle,double CoordsAccuCtrl,double *ps,int ns,double* B,long *pB,int &nB){
    double T=lambda/sin(atan(delta/f/2)),cth,cr,cth2,cr2,tr[4];
    double *xr,*yr,*ths,*rs, th0, r0,x0,y0,thmin,thmax,l,div,*xm,*ym,*thm,*rm,flag;
    long sm=dmax(ps,0,ns);
   // std::cout << std::fixed <<sm<<  '\n';
    xr= (double*) malloc(sm * sizeof(double));
    yr= (double*) malloc(sm * sizeof(double));
    ths= (double*) malloc(sm * sizeof(double));
    rs= (double*) malloc(sm * sizeof(double));
    xm= (double*) malloc(sm * sizeof(double));
    ym= (double*) malloc(sm * sizeof(double));
    thm= (double*) malloc(sm * sizeof(double));
    rm= (double*) malloc(sm * sizeof(double));
    long Ns, Index=0,*Indexm,*Im,Nm,Ns0,countm;
    int SubN,angM,direc;
    nB=0;
    pB[0]=0;
    // ns number of shapes
    for (int i=0;i<ns;i++){
        tr[0]=0;
        tr[1]=0;
        tr[2]=0;
        tr[3]=0;
        direc=1;
        angM=0;
        Ns=(long)ps[i];
        if(nB>=19999){
            std::cout << std::fixed <<"allocated p array is too small"<<  '\n';
            break;
        }

        for(int j=0;j<Ns;j++){
                xr[j]=xs[Index+j];
                yr[j]=ys[Index+j];
        }
        Index=Index+Ns;
        cart2pol(xr,yr,ths,rs,Ns,angM);
        th0=mean(ths,Ns);
        r0=mean(rs,Ns);
        x0=(dmax(xr,0,Ns)+dmin(xr,0,Ns))/2;
        y0=(dmax(yr,0,Ns)+dmin(yr,0,Ns))/2;
        // get accurate coordinates
        getAccurateCoords(xr,yr,ths,rs,x0,y0,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl,Ns);
        if(Ns<3)
            continue;
        
        cart2polCir(xr,yr,ths,rs,Ns,tr,angM);
        thmin=dmin(ths,0,Ns);
        thmax=dmax(ths,0,Ns);
        //compare variance for finding the correct direction for sorting
        countm=0;
        for(int p=0;p<Ns;p++){
            if (ths[p]<=(thmin+(thmax-thmin)/4.0)){
                xm[countm]=rs[p];
                countm++;
            }
        }
        xm[countm]=mean(xm,countm);
        xm[countm+1]=0;
        for(int p=0;p<countm;p++)
            xm[countm+1]+=pow(xm[p]-xm[countm],2);
        xm[0]=xm[countm+1]/(countm*1.0);
        
        countm=0;
        for(int p=0;p<Ns;p++){
            if (ths[p]>=(thmax-(thmax-thmin)/4.0)){
                ym[countm]=rs[p];
                countm++;
            }
        }
        ym[countm]=mean(ym,countm);
        ym[countm+1]=0;
        for(int p=0;p<countm;p++)
            ym[countm+1]+=pow(ym[p]-ym[countm],2);
        
        ym[0]=ym[countm+1]/(countm*1.0);
        if(xm[0]<ym[0]){
            direc=0;//change angle direction
            cth=tr[2];
            cr=tr[3];
        }
        else{
            cth=tr[0];
            cr=tr[1];
        }
        
        // divide the long shape into multiple short shapes for sorting
        l=r0*(thmax-thmin);
        // dicing shape based on its length and zoneplate pitch
        if (dmax(rs,0,Ns)-dmin(rs,0,Ns)>2.5*lambda)
            div=T/5;
        else
            div=T/10;
        if (div>3*r0*PI/180)
            div=3*r0*PI/180;
        if (l>div){
            SubN=ceil(l/div);
            for(int j=0;j<SubN;j++){
                countm=0;
                if(direc==1){
                    flag=thmin+(thmax-thmin)*(j+1)/((double)SubN);
                   
                    for(int p=0;p<Ns;p++){
                        if (ths[p]<=flag){
                            xm[countm]=xr[p];
                            ym[countm]=yr[p];
                            countm++;
                        }
                    }
                    if(countm<5&&j!=(SubN-1)||countm<3)
                        continue;
                    for(int p=0;p<Ns;p++){
                        if (ths[p]<=flag){
                            xr[p]=NULL;
                            yr[p]=NULL;
                            ths[p]=NULL;
                        }
                    }
                }
                else{
                    flag=thmax-(thmax-thmin)*(j+1)/((double)SubN);
                    for(int p=0;p<Ns;p++){
                        if (ths[p]>=flag){
                            xm[countm]=xr[p];
                            ym[countm]=yr[p];
                            countm++;
                        }
                    }
                    if(countm<5&&j!=(SubN-1)||countm<3)
                        continue;
                    for(int p=0;p<Ns;p++){
                        if (ths[p]>=flag){
                            xr[p]=NULL;
                            yr[p]=NULL;
                            ths[p]=NULL;
                        }
                    }
                }
                
                Ns0=removeNULL(xr,Ns);
                Ns0=removeNULL(yr,Ns);
                Ns0=removeNULL(ths,Ns);
                Ns=Ns0;
                cart2pol(xm,ym,thm,rm,countm,angM);

                sortCoods(xm,ym,thm,rm,countm);
                cart2pol(xm,ym,thm,rm,countm,angM);
                if (direc==1)
                    Im=dfind(thm,countm,dmax(thm,0,countm),Nm);
                else
                    Im=dfind(thm,countm,dmin(thm,0,countm),Nm);
                xr[Ns]=xm[Im[0]];
                yr[Ns]=ym[Im[0]];
                ths[Ns]=thm[Im[0]];
                // finding two connection points  
                if (Im[0]==0){
                    if (getIntensity((xm[0]+xm[1])/2.0,(ym[0]+ym[1])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[0]+xm[countm-1])/2.0,(ym[0]+ym[countm-1])/2.0,delta,f,lambda,incidentAngle)){
                        xr[Ns+1]=xm[1];
                        yr[Ns+1]=ym[1];
                        ths[Ns+1]=thm[1];
                    }
                    else{
                        xr[Ns+1]=xm[countm-1];
                        yr[Ns+1]=ym[countm-1];
                        ths[Ns+1]=thm[countm-1];
                    }
                }
                else if(Im[0]==(countm-1)){
                    if (getIntensity((xm[countm-1]+xm[0])/2.0,(ym[countm-1]+ym[0])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[countm-1]+xm[countm-2])/2.0,(ym[countm-1]+ym[countm-2])/2.0,delta,f,lambda,incidentAngle)){
                        xr[Ns+1]=xm[0];
                        yr[Ns+1]=ym[0];
                        ths[Ns+1]=thm[0];
                    }
                    else{
                        xr[Ns+1]=xm[countm-2];
                        yr[Ns+1]=ym[countm-2];
                        ths[Ns+1]=thm[countm-2];
                    }
                }
                else{
                    if (getIntensity((xm[Im[0]]+xm[Im[0]+1])/2.0,(ym[Im[0]]+ym[Im[0]+1])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[Im[0]]+xm[Im[0]-1])/2.0,(ym[Im[0]]+ym[Im[0]-1])/2.0,delta,f,lambda,incidentAngle)){
                        xr[Ns+1]=xm[Im[0]+1];
                        yr[Ns+1]=ym[Im[0]+1];
                        ths[Ns+1]=thm[Im[0]+1];
                    }
                    else{
                        xr[Ns+1]=xm[Im[0]-1];
                        yr[Ns+1]=ym[Im[0]-1];
                        ths[Ns+1]=thm[Im[0]-1];
                    }
                }
                Ns=Ns+2;
                // if the shape is a circle
                if(j==0&&(cth!=0||cr!=0)){
                    Im=dfind(thm,countm,cth,Nm);
                    if(Nm==0)
                        Im=dfind(thm,countm,cth+2*PI,Nm);
                    if(Nm==0){
                        cth=0;
                        cr=0;
                    }
                    else{
                        cth=xm[Im[0]];
                        cr =ym[Im[0]];
                        if (Im[0]==0){
                            if (getIntensity((xm[0]+xm[1])/2.0,(ym[0]+ym[1])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[0]+xm[countm-1])/2.0,(ym[0]+ym[countm-1])/2.0,delta,f,lambda,incidentAngle)){
                                cth2=xm[1];
                                cr2 =ym[1];
                            }
                            else{
                                cth2=xm[countm-1];
                                cr2 =ym[countm-1];
                            }
                        }
                        else if(Im[0]==(countm-1)){
                            if (getIntensity((xm[countm-1]+xm[0])/2.0,(ym[countm-1]+ym[0])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[countm-1]+xm[countm-2])/2.0,(ym[countm-1]+ym[countm-2])/2.0,delta,f,lambda,incidentAngle)){
                                cth2=xm[0];
                                cr2 =ym[0];
                            }
                            else{
                                cth2=xm[countm-2];
                                cr2 =ym[countm-2];
                            }
                        }
                        else{
                            if (getIntensity((xm[Im[0]]+xm[Im[0]+1])/2.0,(ym[Im[0]]+ym[Im[0]+1])/2.0,delta,f,lambda,incidentAngle)>getIntensity((xm[Im[0]]+xm[Im[0]-1])/2.0,(ym[Im[0]]+ym[Im[0]-1])/2.0,delta,f,lambda,incidentAngle)){
                                cth2=xm[Im[0]+1];
                                cr2 =ym[Im[0]+1];
                            }
                            else{
                                cth2=xm[Im[0]-1];
                                cr2 =ym[Im[0]-1];
                            }
                        }
                    }     
                }
                    
                for(int k=0;k<countm;k++){
                    B[pB[nB]+k]=xm[k];
                    B[pB[nB]+k+countm]=ym[k];
                }
                pB[nB+1]=pB[nB]+2*countm;
                nB++;
                if(j==(SubN-1)&&(cth!=0||cr!=0)){
                    xr[Ns]=cth;
                    xr[Ns+1]=cth2;
                    yr[Ns]=cr;
                    yr[Ns+1]=cr2;
                    Ns=Ns+2;
                    cart2pol(xr,yr,ths,rs,Ns,angM);
                    sortCoods(xr,yr,ths,rs,Ns);// sort coordinates
                    for(int k=0;k<Ns;k++){
                        B[pB[nB]+k]=xr[k];
                        B[pB[nB]+k+Ns]=yr[k];
                    }
                    pB[nB+1]=pB[nB]+2*Ns;
                    nB++;
                }
            }
        }
        else{
            sortCoods(xr,yr,ths,rs,Ns);// sort coordinates
            for(int k=0;k<Ns;k++){
                B[pB[nB]+k]=xr[k];
                B[pB[nB]+k+Ns]=yr[k];
            }
            pB[nB+1]=pB[nB]+2*Ns;
            nB++;
        }
    }
    free(xr);
    free(yr);
    free(ths);
    free(rs);
    free(xm);
    free(ym);
    free(thm);
    free(rm);
    free(Im);
}


// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long N            = mxGetM(prhs[0]);
    double *xs        = mxGetPr(prhs[0]);
    double *ys        = mxGetPr(prhs[1]);
    double *ps           = mxGetPr(prhs[2]);
    int  ns           = mxGetN(prhs[2]);
    double delta      = mxGetScalar(prhs[3]);
    double f          = mxGetScalar(prhs[4]);
    double lambda     = mxGetScalar(prhs[5]);
    double th         = mxGetScalar(prhs[6]);
    double incidentAngle    = mxGetScalar(prhs[7]);
    double CoordsAccuCtrl   = mxGetScalar(prhs[8]);
    mxArray *Bs;
    int nB=0,rowLen;
    long pB[20000];
    double *pC;
    double *B;

    B=(double*) malloc(3*N * sizeof(double));
    getBoundariesFromLabel(xs,ys,delta,f,lambda,th,incidentAngle,CoordsAccuCtrl,ps,ns,B,pB,nB);// get accurate and sorted boundary coordinates
    plhs[0]=mxCreateCellMatrix(nB, 1);
    for(int i=0; i< nB ; i++){
        rowLen=(pB[i+1]-pB[i])/2;
        Bs=mxCreateDoubleMatrix(rowLen, 2,mxREAL);
        pC = mxGetPr(Bs);
        for(int j=0; j<rowLen; j++){
            pC[j]=B[pB[i]+j];
            pC[j+rowLen]=B[pB[i]+j+rowLen];
        }
        mxSetCell(plhs[0], i,Bs);
    }
    free(B);
}
