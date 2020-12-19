#include <stdio.h>
#include <mex.h>
#include <math.h>
#include <iostream>
#include <cstring>

//sum
int sum(int * A, int Ibegin, int Iend){
    int temp=A[Ibegin];
    for(int i=Ibegin+1; i<Iend; i++){
        temp+=A[i];
    }
    return temp;
}


// find the specific number
int *dfind(double * A, int num,double B,int &Ns){
    int *Index;
    Index= (int*) malloc(num * sizeof(int));
    int k=0;
    for(int i=0; i<num; i++){
        if (A[4*i+3]==B){
            Index[k]=i;
            k++;
        }
    }
    Ns=k;
    return Index;
}

int *dfind3(double * A, int num,int &Ns){
    int *Index;
    Index= (int*) malloc(num * sizeof(int));
    int k=0;
    for(int i=0; i<num; i++){
        if (A[4*i]==A[0]&&A[4*i+1]==A[1]&&A[4*i+2]==A[2]&&A[4*i+3]!=A[3]){
            Index[k]=i;
            k++;
        }
    }
    Ns=k;
    return Index;
}

int *dfind4(double * A, int num,int &Ns){
    int *Index;
    Index= (int*) malloc(num * sizeof(int));
    int k=0;
    for(int i=0; i<num; i++){
        if (A[4*i]==A[0]&&A[4*i+1]==A[1]&&A[4*i+2]==A[2]&&A[4*i+3]==A[3]){
            Index[k]=i;
            k++;
        }
    }
    Ns=k;
    return Index;
}

// remove NULL
int remove4NULL(double *x,int num){
    int k=0;
    if(num>10000){
        double *y;
        y= (double*) malloc(4*num * sizeof(double));
        for(int i=0;i<num;i++){
            if (x[4*i]!=NULL){
                y[4*k]=x[4*i];
                y[4*k+1]=x[4*i+1];
                y[4*k+2]=x[4*i+2];
                y[4*k+3]=x[4*i+3];
                k++;
            }
        }
        for(int i=0;i<k;i++){
            x[4*i]  = y[4*i];
            x[4*i+1]= y[4*i+1];
            x[4*i+2]= y[4*i+2];
            x[4*i+3]= y[4*i+3];
        }
        free(y);
    }
    else {
        double y[40000];
        for(int i=0;i<num;i++){
            if (x[4*i]!=NULL){
                y[4*k]=x[4*i];
                y[4*k+1]=x[4*i+1];
                y[4*k+2]=x[4*i+2];
                y[4*k+3]=x[4*i+3];
                k++;
            }
        }
        for(int i=0;i<k;i++){
            x[4*i]  = y[4*i];
            x[4*i+1]= y[4*i+1];
            x[4*i+2]= y[4*i+2];
            x[4*i+3]= y[4*i+3];
        }
    }
    return k;
}

// merge boundary shapes
void mergeBoundaryShapes(mxArray *Ixy,double *uxy,int num,int *count,int &Ns){
    int t,*Index,Nu=0,k;
    mxArray *It,*Ik,*I;
    double *s1,*s2,*s;
    while(num>0){
        t=uxy[3];
        Index=dfind3(uxy,num,Nu);
        if (Nu==0){
            Index=dfind4(uxy,num,Nu);
            if (Nu==0)
                break;
            for(int i=0;i<Nu;i++)
                uxy[4*Index[i]]=NULL;//remove error boundary shapes
            num=remove4NULL(uxy,num);
            //make sure there is no shape marked number t and save shape index
            Index=dfind(uxy,num,t,Nu);
            if(Nu==0){
                count[Ns]=t;
                Ns++;
            }
            continue;
        }
        k=uxy[4*Index[0]+3];
        //merge data to Ixy(t).xr
        It=mxGetField(Ixy, t-1, "xr");
        Ik=mxGetField(Ixy, k-1, "xr");
        I=mxCreateDoubleMatrix(mxGetM(It)+mxGetM(Ik),1, mxREAL);
        
        s1=mxGetPr(It);
        s2=mxGetPr(Ik);
        s=mxGetPr(I);
      
        
        for(int m=0;m<mxGetM(It);m++)
            s[m]=s1[m];
        for(int m=0;m<mxGetM(Ik);m++)
            s[m+mxGetM(It)]=s2[m];
        mxSetField(Ixy, t-1,"xr",I);
        //merge data to Ixy(t).yr
        It=mxGetField(Ixy, t-1, "yr");
        Ik=mxGetField(Ixy, k-1, "yr");
        I=mxCreateDoubleMatrix(mxGetM(It)+mxGetM(Ik),1, mxREAL);
        s1=mxGetPr(It);
        s2=mxGetPr(Ik);
        s=mxGetPr(I);
        for(int m=0;m<mxGetM(It);m++)
            s[m]=s1[m];
        for(int m=0;m<mxGetM(Ik);m++)
            s[m+mxGetM(It)]=s2[m];
        mxSetField(Ixy, t-1,"yr",I);
        uxy[0]=NULL;
        uxy[4*Index[0]]=NULL;
        num=remove4NULL(uxy,num);
        //modify the number of the merged boundary shapes
        if (num>0)  {
            Index=dfind(uxy,num,k,Nu);
            for(int i=0;i<Nu;i++)
                uxy[4*Index[i]+3]=t;
            //make sure there is no shape marked number t and save shape index
            Index=dfind(uxy,num,t,Nu);
            if (Nu==0){
                count[Ns]=t;
                Ns++;
            }
        }
        else{
            count[Ns]=t;
            Ns++;
        }
    }
}


// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int Mx       = mxGetScalar(prhs[1]);
    int My       = mxGetScalar(prhs[2]);
    int *Num0,*Num,*Nuxy,k0=0,k=0,ku=0,sns=0,Ns=0,*count0;
    Num0= (int*) malloc(Mx*My* sizeof(int));
    Num= (int*) malloc(Mx*My* sizeof(int));
    Nuxy= (int*) malloc(Mx*My* sizeof(int));
    mxArray *B,*C,*Ixy0,*Ixy,*Puxy;
    double *uxy,*uxy0,*s1,*s2,*count;
    
    for(int i=0; i< Mx*My ; i++){
        B=mxGetCell(prhs[0], i);
        C=mxGetField(B, 0, "Ixy0");
        Num0[i]=mxGetN(C);
        C=mxGetField(B, 0, "Ixy");
        Num[i]=mxGetN(C);
        C=mxGetField(B, 0, "uxy");
        Nuxy[i]=mxGetM(C);
    }
    const char *field_names[] = {"xr","yr"};
    //extract coordinates from ceout
    plhs[0]=mxCreateStructMatrix(1,sum(Num0,0,Mx*My),2,field_names);
    plhs[1]=mxCreateStructMatrix(1,sum(Num,0,Mx*My),2,field_names);
    count0= (int*) malloc(sum(Num,0,Mx*My)* sizeof(int));
    Puxy=mxCreateDoubleMatrix(4,sum(Nuxy,0,Mx*My), mxREAL);
    uxy=mxGetPr(Puxy);
    for(int i=0; i< Mx*My ; i++){
        B=mxGetCell(prhs[0], i);
        C=mxGetField(B, 0, "Ixy0");
        for(int j=0;j<Num0[i];j++){
            //transfer data to Ixy0.xr
            Ixy=mxGetField(C, j, "xr");
            Ixy0=mxCreateDoubleMatrix(mxGetM(Ixy),1, mxREAL);
            s1=mxGetPr(Ixy0);
            s2=mxGetPr(Ixy);
            for(int m=0;m<mxGetM(Ixy);m++){
                s1[m]=s2[m];
            }
            mxSetField(plhs[0], k0,"xr",Ixy0);
            //transfer data to Ixy0.yr
            Ixy=mxGetField(C, j, "yr");
            Ixy0=mxCreateDoubleMatrix(mxGetM(Ixy),1, mxREAL);
            s1=mxGetPr(Ixy0);
            s2=mxGetPr(Ixy);
            for(int m=0;m<mxGetM(Ixy);m++){
                s1[m]=s2[m];
            }
            mxSetField(plhs[0], k0,"yr",Ixy0);
            k0++;
        }
        C=mxGetField(B, 0, "Ixy");
        uxy0=mxGetPr(mxGetField(B, 0, "uxy"));
        for(int j=0;j<Num[i];j++){
            //transfer data to Ixy.xr
            Ixy=mxGetField(C, j, "xr");
            Ixy0=mxCreateDoubleMatrix(mxGetM(Ixy),1, mxREAL);
            s1=mxGetPr(Ixy0);
            s2=mxGetPr(Ixy);
            for(int m=0;m<mxGetM(Ixy);m++){
                s1[m]=s2[m];
            }
            mxSetField(plhs[1], k,"xr",Ixy0);
            //transfer data to Ixy.yr
            Ixy=mxGetField(C, j, "yr");
            Ixy0=mxCreateDoubleMatrix(mxGetM(Ixy),1, mxREAL);
            s1=mxGetPr(Ixy0);
            s2=mxGetPr(Ixy);
            for(int m=0;m<mxGetM(Ixy);m++){
                s1[m]=s2[m];
            }
            mxSetField(plhs[1], k,"yr",Ixy0);
            k++;
        }
        for(int j=0;j<Nuxy[i];j++){
            //transfer data to uxy
            uxy[4*ku]=uxy0[j];
            uxy[4*ku+1]=uxy0[j+Nuxy[i]];
            uxy[4*ku+2]=uxy0[j+2*Nuxy[i]];
            uxy[4*ku+3]=uxy0[j+3*Nuxy[i]]+sns;
            ku++;
        }
        sns+=Num[i];
    }
    
    //merge boundary shapes
    mergeBoundaryShapes(plhs[1],uxy,sum(Nuxy,0,Mx*My),count0,Ns);
    plhs[2]=mxCreateDoubleMatrix(1, Ns, mxREAL);
    count=mxGetPr(plhs[2]);
    for(int m=0;m<Ns;m++)
        count[m]=count0[m];
    free(Num0);
    free(Num);
    free(Nuxy);
    free(count0);
}
