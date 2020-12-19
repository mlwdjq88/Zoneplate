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

// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int Mx       = mxGetScalar(prhs[1]);
    int My       = mxGetScalar(prhs[2]);
    int *Num0,*Num,*Nuxy,k0=0,k=0,ku=0,sns=0;
    Num0= (int*) malloc(Mx*My* sizeof(int));
    Num= (int*) malloc(Mx*My* sizeof(int));
    Nuxy= (int*) malloc(Mx*My* sizeof(int));
    mxArray *B,*C,*Ixy0,*Ixy;
    double *uxy,*uxy0,*s1,*s2;
    
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
    plhs[0]=mxCreateStructMatrix(1,sum(Num0,0,Mx*My),2,field_names);
    plhs[1]=mxCreateStructMatrix(1,sum(Num,0,Mx*My),2,field_names);
    plhs[2]=mxCreateDoubleMatrix(4,sum(Nuxy,0,Mx*My), mxREAL);
    uxy=mxGetPr(plhs[2]);
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
    free(Num0);
    free(Num);
    free(Nuxy);
}
