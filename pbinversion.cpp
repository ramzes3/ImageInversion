//////////////////////////////////////////////////////////////////////////////////////
//  This part of the program is based on the old program written by Gustavo Garcia ///
//  The functions polar_b(), theta_f(), cubic() and ldist() are published with     ///
//  permission of Gustavo Garcia. For more information on the image inversion      ///
//  with pBasex please read the original paper:                                    ///
//  G. A. Garcia, L. Nahon, and I. Powis, “Two-dimensional charged particle image  ///
//  inversion using a polar basis function expansion,” Review of Scientific Instru-///
//  ments, vol. 75, no. 11, pp. 4989–4996, 2004.                                   ///
//////////////////////////////////////////////////////////////////////////////////////


#include "pbinversion.h"

#include "time.h"
#include <string.h>
#include <math.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// Qt support for read / write
#include <QFile>
#include <QDataStream>

void* p2Object;        // global variable which points to an arbitrary object


PBInversion::PBInversion()
{
    //////////////////// O L D   //////////////
    NR  = 256; // Radial binning
    NTH = 256; // Angular binning (180 degrees)
    NK  = 128; // Number of basis functions
    //////////////////// O L D   //////////////

    ///////////  N E W   //////////////////////////////
    //size_r = 256;
    //size_theta = 256;
    //size_basis = 128;

    ncolumns = 512; nrows = 512; // default cartesian array size
    S = NULL; V = NULL; U = NULL; ang = NULL;
    m_isloaded = 0;

    gsl_set_error_handler_off(); //  switch off default gsl error hanlder
}

PBInversion::~PBInversion()
{
    if(S) { gsl_matrix_free(V); gsl_matrix_free(U); gsl_vector_free(S); }
    if(ang) {
        for(int i=0;i<ang_size;i++) delete [] ang[i];
        delete [] ang;
    }
}

void PBInversion::run() // running a separate thread to creat basis set
{
    emit sent_respond( "thread is started" );
    write_forward(nl_v, odd_v);
    return;
    //load_matrices(nl_v, odd_v);
    //polar_b(data_v, Imrz_v, xc_v, yc_v, dr_v, nl_v, odd_v);
}

void PBInversion::polar_a(int nl,int odd)
{
    nl_v = nl; odd_v = odd;
}

void PBInversion::set_parameters(int nl,int odd)
{
    nl_v = nl; odd_v = odd;
}

/*
void PBInversion::Invert_Polar(double *input_array, double *output_array, int xc,int yc, int radius, int l_order, bool symmetric)
{
    gsl_vector *polar_array, *expansion_coeff;

    int n_leg; // number of Legandre polinomials
    if(symmetric) n_leg=l_order/2+1;
    else n_leg=l_order+1;

    int size_coef = size_basis*n_leg; // size of the array with basis coefficients

    polar_array=gsl_vector_calloc(size_r*size_theta);
    // copy input_array into polar gsl vector
    for (int i=0; i<size_theta*size_r; i++){gsl_vector_set(polar_array,i,input_array[i]);}

    /// Solving SVD to calculate basis coefficients ///
    emit sent_respond( "solving svd" );

    expansion_coeff=gsl_vector_calloc(size_coef);
    int error_code = gsl_linalg_SV_solve(U,V,S,polar_array,expansion_coeff);
    if (error_code) sent_error("SVD solve error");
    gsl_vector_free(polar_array);

    /// Extracting the PES and angular parameters ///


}
*/

void PBInversion::polar_b(double **data, double **Imrz,int xc,int yc,int dr,int nl,int odd)
{
    gsl_vector *polar,*coeff, *polar_large;
    int ix,iz,ir,it,ik,il,kmin,kmax,i,j;
    double pmax,psum,f,a0;
    double dx,dy,thta,rad;
    //double cubic(double);
    //double ldist(double , gsl_vector *, double [], int, int, int);
    double x,y,p; //P12,P34,
    int lowx,lowy;
    double fact, fact2;
    int N,NL;
    double *pl,leg,cpu0,cpu1,sum;

    cpu0=(double)clock()/(double)CLOCKS_PER_SEC;

    if(odd) NL=nl+1;
    else NL=nl/2+1;


    /* Re-initialise angular array */
    if(ang) {
        for(int k=0;k<ang_size;k++) delete [] ang[k];
        delete [] ang;
        ang = NULL;
    }

    ang = new double *[NL];
    for(int k=0;k<NL;k++) ang[k] = new double [NR];
    for(int ki=0;ki<NL;ki++){
        for(int kj=0;kj<NR;kj++){
            ang[ki][kj] = 0;
        }
    }
    ang_size = NL;

    //for(i=0;i<NL;i++)
      //memset(ang[i],NULL,NR*sizeof(double));

    N=NK*NL;

    emit sent_respond( "symmetrise image" );

    /* Symmetrise the image along the y axis since the Legendre        */
    /* polynomials are naturally symmetric along the polarisation axis */
    /* If the user wants to only use one half of the image, comment out*/
    /* this part and then choose the side in the expression below      */
    /* x=xc-rad*sin(thta) will use the left side and xc+... the right  */
    //The image has to be symmetrised before entering this function
    /*
    for(iz=0;iz<nrows;iz++){
      for(ix=xc;ix<ncolumns;ix++){
        if(ix<=(2*xc))
          data[iz][2*xc-ix]+=data[iz][ix];
      }
    }
*/
    emit sent_respond( "converting to polar coordinates" );
    /* Convert to polar co-ordinates */
    /* using bicubic interpolation   */

    polar=gsl_vector_calloc(NR*NTH);

    int size_l1, size_l2, scale1, scale2;
    scale1 = ceil( nrows / (double) NR ) ;
    size_l1 = NR * scale1;
    scale2 = ceil( ncolumns / (double) NTH ) ;
    size_l2 = NTH * scale2;
    polar_large = gsl_vector_calloc(size_l1*size_l2);

    if(ncolumns<=nrows) {fact=(double)ncolumns; fact2=(double)ncolumns;}
    else {fact=(double)nrows; fact2=(double)nrows;}


    fact=fact/2.0/(double)NR;
    //fact2=fact2/2.0/(double)nrows;
    fact2=fact2/2.0/(double)size_l1;

    // First converting cartesian matrix to similar size polar matrix (polar_large)
    for(ir=0; ir < size_l1; ir++){
      rad=(double)ir*fact2;
      for(it=0;it<size_l2;it++){
        thta=it/(double)(size_l2-1)*M_PI;
        x=xc-rad*sin(thta); // If the image hasn't been symmetrised use the sign to choose sides
        y=yc+rad*cos(thta);
        lowx=(int)x;
        lowy=(int)y;
        dx=x-lowx;
        dy=y-lowy;
        p=0.0;
        for (i=-1;i<=2;i++){
          for (j=-1;j<=2;j++){
            if((lowy+i)<nrows&&(lowx+j)<ncolumns&&(lowy+i)>=0&&(lowx+j)>=0){
                p+=data[lowy+i][lowx+j]*cubic(j-dx)*cubic(dy-i)/16.0;
            }
          }
        }
        gsl_vector_set(polar_large,ir*size_l2+it,p);
      }
    }

    //char work_file[255];
    //sprintf(work_file,"D:/Roman/Data/20130531/test1.dat");
    //save_gsl_vector(polar_large, work_file);

    // converting polar_large matrix into smaller matrix "polar" which matches the dimenssions of the pBasex basis functions

    for(ir=0; ir<NR; ir++){
      for(it=0;it<NTH;it++){
          p=0.0;
          for(i=0;i<scale1;i++){
            for(j=0;j<scale2;j++){
                p += gsl_vector_get(polar_large,(ir*scale1+i)*size_l2+(it*scale2+j))/(scale1*scale2);
            }
          }
          //p /= (scale1*scale2);
          gsl_vector_set(polar,ir*NTH+it,p);
      }
    }

    //char work_file[255];
    //sprintf(work_file,"D:/Roman/Data/20130531/test2.dat");
    //save_gsl_vector(polar, work_file);

    gsl_vector_free(polar_large);

    emit sent_respond( "solving svd" );

    /*Solve SVD to calculate coefficients*/

    coeff=gsl_vector_calloc(N);
    int error_code = gsl_linalg_SV_solve(U,V,S,polar,coeff);
    if (error_code) sent_error("SVD solve error");
    gsl_vector_free(polar);

    /* Construct the PES and calculate angular parameters */

    dr=(int)(dr/fact);

    pmax=0.0;
    for(ir=0;ir<dr;ir++){ // dr
      kmin=ir/2-5;
      if(kmin<0) kmin=0;
      kmax=ir/2+6;
      if(kmax>(NK-1)) kmax=NK-1;
      for(il=0;il<NL;il++){
        ang[il][ir]=0.0;
      }
      for(ik=kmin;ik<kmax;ik++){
        f=exp(-(ir-ik*2)*(ir-ik*2)/2.0);
        for(il=0;il<NL;il++)
          ang[il][ir]+=gsl_vector_get(coeff,ik*NL+il)*f;
      }
      a0=ang[0][ir];
      if(a0<=0){
        //a0=0;
        //for(il=0;il<NL;il++) ang[il][ir]=0.0;
      }
      else
          for(il=1;il<NL;il++) ang[il][ir]/=a0;

      //if(a0*ir*ir>pmax) pmax=a0*ir*ir;
    }

    psum = 0;
    for(ir=0;ir<dr;ir++) psum += ang[0][ir]*(ir*ir); // Normalise PES
    for(ir=0;ir<dr;ir++) ang[0][ir]*=(ir*ir)/psum; // Normalise PES  ang[0][ir]*=(ir*ir/pmax);

    /* Create cartesian image for display */

    pl=(double*)calloc(nl+1,sizeof(double)); // Allocate Legendre polynomials

    for(iz=0;iz<nrows;iz++){
      dy=(double)(iz-yc);
       for(ix=xc;ix<ncolumns;ix++){
         dx=(double)(ix-xc);
         rad=sqrt(dx*dx+dy*dy);
         rad/=fact;
         sum=0.0;
         if (rad<=dr){
           thta=theta_f(dx,dy);
           kmin=(int)(rad/2.0)-5;
           if(kmin<0) kmin=0;
           kmax=(int)(rad/2.0)+6;
           if(kmax>NK-1) kmax=NK-1;
           for(ik=kmin;ik<=kmax;ik++){
             f=exp(-(rad-ik*2)*(rad-ik*2)/2.0);
             leg=ldist(thta,coeff,pl,nl,odd,ik);
             leg*=f;
             sum+=leg;
           }
           if(sum>=0) Imrz[iz][ix]=sum;//*rad*rad;
           else Imrz[iz][ix]=0.0;
           Imrz[iz][2*xc-ix]=Imrz[iz][ix];
         }
         else Imrz[iz][ix]=0.0;
       }
    }
    cpu1=(double)clock()/(double)CLOCKS_PER_SEC;
    gsl_vector_free(coeff);
    free(pl);
    //printf("Done in %f seconds\n",cpu1-cpu0);
    char buffer[100];
    sprintf(buffer, "Done in %2.2f seconds",cpu1-cpu0 );
    emit sent_respond(buffer);

}


double PBInversion::theta_f(double x,double y)
{
    double thta;
    thta=atan(fabs(x)/fabs(y));
    if (x<0&&y>0) thta=2.0*M_PI-thta;
    if (x>0&&y<0) thta=M_PI-thta;
    if (x<0&&y<0) thta=M_PI+thta;
    if(x==0&&y>0) thta=0.0;
    if(x==0&&y<0) thta=M_PI;
    if(y==0&&x>0) thta=M_PI/2.0;
    if(y==0&&x<0) thta=3.0*M_PI/2.0;
    if(x==0&&y==0) thta=2.0*M_PI;
    return thta;
}

double PBInversion::cubic(double x)
{
    double p0,p1,p2,p3;
    if((x+2)>0) p0=(x+2)*(x+2)*(x+2);
    else p0=0.0;
    if((x+1)>0) p1=(x+1)*(x+1)*(x+1);
    else p1=0.0;
    if((x)>0) p2=(x)*(x)*(x);
    else p2=0.0;
    if((x-1)>0) p3=(x-1)*(x-1)*(x-1);
    else p3=0.0;
    return 1.0/6.0*(p0-4.0*p1+6.0*p2-4.0*p3);
}

double PBInversion::ldist(double X, gsl_vector *coeff, double *pl, int npl, int odd, int ik)
{
    double v=0.0;
    int i,index;
    int j;
    double twox,f2,f1,d,x;
    int NL;

    if(odd) NL=npl+1;
    else NL=npl/2+1;

    x = cos(X);
    pl[0]=1.0;
    pl[1]=x;
    if (npl >= 2) {
      twox=2.0*x;
      f2=x;
      d=1.0;
      for (j=2;j <= npl ;j++) {
        f1=d++;
        f2 += twox;
        pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
      }
    }
    for (i=0; i<=npl; i++){
      if((i%2)){
        if(odd)
          v += (double)(gsl_vector_get(coeff,ik*NL+i)*pl[i]);
      }
      else{
        if(odd)
          v += (gsl_vector_get(coeff,ik*NL+i)*pl[i]);
        else{
          index=i/2;
          v += (gsl_vector_get(coeff,ik*NL+index)*pl[i]);
        }
      }
    }
    return v;

}

void PBInversion::load_matrices(int l,int odd)
{
  FILE *fu; //*fs,*fu,*fv;
  char fileV[50],fileU[50],fileS[50],file[50],str[50]; //,alert[50];
  //static int L=0,O=0;
  int i, error_code;
  //static short LOADED=0;
  int N,M,NL;

  //if(l==L&&O==odd) /* Already loaded */
  //  return;
  emit sent_respond("set to zero");
  /* Set to zero */
  memset(fileV,NULL,50);
  memset(fileU,NULL,50);
  memset(fileS,NULL,50);
  memset(file,NULL,50);

  strcat(file,"");// Add root directory
  strcat(fileU,file);
  strcat(fileV,file);
  strcat(fileS,file);
  //memset(file,NULL,50);


  /* Describe matrices' names */
  for(i=0;i<=l;i++){
    if((i%2)){
      sprintf(str,"P%d",i);
      if(odd) strcat(file,str);
    }
    else{
      sprintf(str,"P%d",i);
      strcat(file,str);
    }
  }

  sprintf(str,"NR%dNTH%dNK%d",NR,NTH,NK);
  strcat(file,str);

  strcat(file,".dat"); // Add extension

  strcat(fileU,"U");
  strcat(fileV,"V");
  strcat(fileS,"S");

  strcat(fileU,file);
  strcat(fileV,file);
  strcat(fileS,file);

  if(odd) NL=l+1;
  else NL=l/2+1;

  N=NK*NL;
  M=NTH*NR;

  emit sent_respond("loading matrices");
  /* Load matrices */

  if(S){
    gsl_matrix_free(V); V = NULL;
    gsl_matrix_free(U); U = NULL;
    gsl_vector_free(S); S = NULL;
  }

  //fv=fopen(fileV,"rb");
  //if(fv==NULL){
 //   m_isloaded = 0;
 //   emit sent_respond("cannot open file"); return;
    //printf("Error from svd.c: Couldn't open %s \n",fileV);
    //sprintf(alert,"Could not find the basis set for l=%d, odd=%d, do you want to create it?\n(This may take a few hours)",l,odd);
    //if(fl_show_question(alert, 0)){
    //write_forward(l,odd);
      //fl_show_alert("PROGRAM NEEDS TO RESTART","","",0);
      //exit(0);
    //}
    //else
    //exit(0);
  //}

  emit sent_respond("loading v");
  V=gsl_matrix_alloc(N,N);
  error_code = this->read_gsl_matrix(V, fileV);
  if (error_code) sent_error("couldn't read the file fv");
  //fclose(fv);

  fu=fopen(fileU,"rb");
  if(fu==NULL){
    m_isloaded = 0;
    emit sent_respond("cannot open file u"); return;
    //printf("Error from svd.c: Couldn't open %s \n",fileU);
    exit(0);
  }

  U=gsl_matrix_alloc(M,N);
  //error_code = gsl_matrix_fread(fu,U);
  error_code = read_gsl_matrix(U, fileU);
  if (error_code) sent_error("couldn't read the file fu");
  //fclose(fu);
  emit sent_respond("loading u");
  //fs=fopen(fileS,"rb");
  //if(fs==NULL){
  //  m_isloaded = 0;
  //  emit sent_respond("cannot open file s"); return;
    //printf("Error from svd.c: Couldn't open %s \n",fileS);
 //   return; //exit(0);
  //}
  S=gsl_vector_alloc(N);
  error_code = read_gsl_vector(S, fileS);
  if (error_code) sent_error("couldn't read the file fv");
  //fclose(fs);
  emit sent_respond("loading s");
  m_isloaded = 1;
  /* Allocate angular space */
  //if(LOADED) delete [] ang; //fl_free_matrix(ang);

  //ang = new double*[NL]; // ang=fl_get_matrix(NL,NR,sizeof(double));
  // //for(int var=0;var<NL;var++) ang[var] = new double [NR];
  //if(ang==NULL){
  //  printf("Couldn't load angular matrix\n");
  //  exit(0);
  //}

  //LOADED=1;

  //L=l; // Set L and Odd values to avoid reloading matrices if not needed
  //O=odd;

  double result = gsl_matrix_get(U,3,3);
  char buffer[255]; // buffer for messages
  sprintf(buffer,"Result U %f ",result);
  emit sent_respond(buffer);
  sleep(3);
}

void PBInversion::write_forward(int l,int odd)
{
    gsl_set_error_handler_off(); //  switch off default gsl error hanlder
    int ir,it,ik,il,nl,basis_index,point_index,NL,N,M,i;
    double p,dx,dz,thta,rad;
    double cpu1,cpu0;
    gsl_matrix *Basis,*BasisT,*V,*X;
    BasisT = NULL; Basis = NULL;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(100000); // Allocate space for integration
    gsl_vector *W,*S;
    //FILE *fs,*fu,*fv;
    char fileV[50],fileU[50],fileS[50],file[50],str[50];

    if(odd) NL=l+1;
    else NL=l/2+1;
    N=NL*NK;
    M=NR*NTH;
    Basis=gsl_matrix_alloc(N,M); //Allocate Basis matrix
    if (Basis==NULL){
      sent_error("Couldn't allocate memory for basis functions\n");
      this->terminate();
    }
    emit sent_respond(" step 1"); // debugger;

    for(ik=0;ik<NK;ik++){ //  if not 0 change back
      for(il=0;il<NL;il++){
        basis_index=ik*NL+il;
        cpu0=(double)clock()/(double)CLOCKS_PER_SEC;
        if(odd) nl=il;
        else nl=il*2;
        for(ir=0;ir<NR;ir++){
          for(it=0;it<NTH;it++){
            point_index=ir*NTH+it;
            rad=(double)ir;
            thta=it/(double)(NTH-1)*M_PI;
            dx=rad*sin(thta);
            dz=rad*cos(thta);
            p=model_c(ik,nl,dx,dz,w);
            gsl_matrix_set(Basis,basis_index,point_index,p);
          }
        }
        cpu1=(double)clock()/(double)CLOCKS_PER_SEC;
        char buffer[255]; // buffer for messages
        sprintf(buffer,"Basis set %d,%d (index %d of %d) done in %f seconds",ik,il,basis_index,NK*NL,cpu1-cpu0);
        //printf("Basis set %d,%d (index %d of %d) done in %f seconds\n",ik,il,basis_index,NK*NL,cpu1-cpu0);
        emit sent_respond(buffer);
      }
    }
    gsl_integration_workspace_free(w);

    //printf("Decomposing matrix, please wait ...\n");
    emit sent_respond(" Decomposing matrix, please wait ..."); // debugger;

    /***********************/
    /* Decomposition stuff */
    /***********************/
    double result = gsl_matrix_get(Basis,3,3);
    char buffer[255]; // buffer for messages
    sprintf(buffer,"Result1 %f ",result);
    emit sent_error(buffer);
    sleep(3);

    /* Transpose it */

    BasisT=gsl_matrix_alloc(M,N);
    if (BasisT==NULL){
      sent_error("Couldn't allocate memory for basisT functions\n");
      gsl_matrix_free(Basis);
      return;
      //exit(0);
    }
    int error_code = gsl_matrix_transpose_memcpy(BasisT,Basis);
    if (error_code) sent_error("Transpose error");
    gsl_matrix_free(Basis);

    /*Allocate working space for decomposition */
    result = gsl_matrix_get(BasisT,3,3);
    sprintf(buffer,"Result2 %f ",result);
    emit sent_error(buffer);
    sleep(3);


    S=gsl_vector_alloc(N);
    V=gsl_matrix_alloc(N,N);
    W=gsl_vector_alloc(N);
    X=gsl_matrix_alloc(N,N);

    /* Decompose */
    error_code = gsl_linalg_SV_decomp_mod(BasisT,X,V,S,W);
    if (error_code) sent_error("decomposition error");
    gsl_matrix_free(X);
    gsl_vector_free(W);

    /* Describe matrices' names */

    memset(fileV,NULL,50);
    memset(fileU,NULL,50);
    memset(fileS,NULL,50);
    memset(file,NULL,50);

    for(i=0;i<=l;i++){
      if((i%2)){
        sprintf(str,"P%d",i);
        if(odd) strcat(file,str);
      }
      else{
        sprintf(str,"P%d",i);
        strcat(file,str);
      }
    }

    sprintf(str,"NR%dNTH%dNK%d",NR,NTH,NK);
    strcat(file,str);

    strcat(file,".dat"); // Add extension

    sprintf(fileU,"U");
    sprintf(fileV,"V");
    sprintf(fileS,"S");

    strcat(fileU,file);
    strcat(fileV,file);
    strcat(fileS,file);

    /* Write everything to file */

    //fu=fopen(fileU,"w");
    //error_code = gsl_matrix_fwrite(fu,BasisT);
    //if (error_code) sent_error("couldn't save the file fu");
    //fclose(fu);
    error_code = save_gsl_matrix(BasisT, fileU);
    if (error_code) sent_error("couldn't save the file fu");

    result = gsl_matrix_get(BasisT,3,3);
    sprintf(buffer,"Result2 %f ",result);
    emit sent_error(buffer);
    sleep(3);

    gsl_matrix_free(BasisT); BasisT = NULL;

    //fv=fopen(fileV,"w");
    error_code = save_gsl_matrix(V, fileV);
    sprintf(buffer,"couldn't save the file fv, Error %d ",error_code);
    if (error_code) sent_error(buffer);
    gsl_matrix_free(V); V = NULL;
    //fclose(fv);

    //fs=fopen(fileS,"w");
    error_code = save_gsl_vector(S, fileS);
    if (error_code) sent_error("couldn't save the file fs");
    gsl_vector_free(S); S = NULL;
    //fclose(fs);

    emit sent_respond(" basis functions are saved"); // debugger;
}

double PBInversion::model_f(double r,void *par)
{
  //double lpoly(double,double [],int);
  //double *pl = NULL;
  double *xz=(double*)par;
  //double theta_f(double,double);
  double thta=this->theta_f(r,xz[1]);
  double R=r/sin(thta), f;
  int k=(int)xz[2]*2; //change for different pixel widths
  int l=(int)xz[3];

  //pl=(double*)calloc(l+1,sizeof(double));
  //pl = new double[l+1];

  f=exp(-(R-k)*(R-k)/2.0);
  f*=lpoly(thta,l);
  //free(pl);
  //delete [] pl;

  return f*r/sqrt(r*r-xz[0]*xz[0]);
}

double PBInversion::model_c(int k,int l,double x,double z,gsl_integration_workspace *w)
{
  double par[4]={x,z,k,l};
  gsl_function f;
  //size_t neval;
  double error, result;
  p2Object = (double*) this;
  //f.function=&model_f(double, void*);
  f.function=PBInversion::Wrapper_To_Call_model_f;
  f.params=&par;
  int error_code; bool flag = 0;
  char buffer[255]; // buffer for messages
  for (int loop=0;loop<6;loop++){
    // apper integration limit was 500.0, it is changed now to NR+50.0
    error_code = gsl_integration_qag(&f,fabs(x),NR+50.0,0.00001,0.00005,100000,6-loop,w,&result,&error);
    if (error_code) {
        sprintf(buffer,"there is an error in %d,%d, attempt %d",k , l, loop);
        emit sent_error(buffer);
        flag = 1;
    }
    else {break;}
  }
  if (flag) {sprintf(buffer,"result: %f", result); emit sent_error(buffer);}
  return result;
}

/* Return Legendre polynomial at X */

double PBInversion::lpoly(double X, int npl)
{
  int j;
  double twox,f2,f1,d,x;

  double pl[npl+1];
  //pl = new double[npl+1];

  x=(double)cos(X);
  pl[0]=1.0;
  pl[1]=x;
  if (npl >= 2) {
    twox=2.0*x;
    f2=x;
    d=1.0;
    for (j=2;j <= npl ;j++) {
      f1=d++;
      f2 += twox;
      pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d;
    }
  }

  f1 = pl[npl];
  //delete [] pl; pl = NULL;
  return f1; //pl[npl];
}

void PBInversion::set_size(int size_x, int size_y)
{
    ncolumns = size_x;
    nrows    = size_y;
}

double PBInversion::Wrapper_To_Call_model_f(double r, void *par)
{
    PBInversion* mySelf = (PBInversion*) p2Object;
    // call member
    return mySelf->model_f(r, par);
}

int PBInversion::save_gsl_matrix(const gsl_matrix * matrix , char *filename)
{
    QFile file(filename);
    if(!file.open(QIODevice::WriteOnly)){ sent_error("couldn't open the file"); return -1;};
    QDataStream out(&file);
    for(int i=0; i < (matrix->size1 * matrix->size2);i++){
        out << (double) matrix->data[i];
    }
    file.close();
    return 0;
}

int PBInversion::read_gsl_matrix(const gsl_matrix * matrix, char *filename)
{
    QFile file(filename);
    if(!file.open(QIODevice::ReadOnly)){sent_error("couldn't open the file"); return -1;};
    QDataStream in(&file);    // read the data serialized from the file
    double a;
    for(int i=0; i < (matrix->size1 * matrix->size2);i++){
         in >> a;           // extract element
         matrix->data[i] = a;
    }
    file.close();
    return 0;
}

int PBInversion::save_gsl_vector(const gsl_vector * matrix , char *filename)
{
    QFile file(filename);
    if(!file.open(QIODevice::WriteOnly)){sent_error("couldn't open the file"); return -1;};;
    QDataStream out(&file);
    for(int i=0; i < (matrix->size);i++){
        out << (double) matrix->data[i];
    }
    file.close();
    return 0;
}

int PBInversion::read_gsl_vector(const gsl_vector * matrix, char *filename)
{
    QFile file(filename);
    if(!file.open(QIODevice::ReadOnly)){sent_error("couldn't open the file"); return -1;};;
    QDataStream in(&file);    // read the data serialized from the file
    double a;
    for(int i=0; i < (matrix->size);i++){
         in >> a;           // extract element
         matrix->data[i] = a;
    }
    file.close();
    return 0;
}

void PBInversion::generateTestVMI(const int dim, double fwhm, double r_centre, double b2){


    char fileout[255], buffer[255];
    sprintf(fileout,"VMItest.%s", "txt");
    ofstream outFile2(fileout, ios::out);
    if(!outFile2){ return; }

    //const int sizu = 4*dim*dim;
    //double* marray;
    //marray = new double [sizu];
    double var, sig;
    sig = sqrt(2)*fwhm/2.35;

    for(int x=-dim+1;x<dim;x++){
       for(int z=-dim+1;z<dim;z++){
           var = 0.0;
           for(int y=0;y<dim;y++){
               if ((x==0)&&(y==0)&&(z==0)) var += exp(-pow((sqrt(x*x + y*y + z*z) - r_centre)/sig,2));
               else var += exp(-pow((sqrt(x*x + y*y + z*z) - r_centre)/sig,2)) * (1 + b2/2*(3*pow(z/sqrt(x*x + y*y + z*z),2) - 1));
           }
           //marray[z+dim+(x+dim)*2*dim] = var;
           outFile2 << "\t" <<var;
       }
       //sprintf(buffer,"line: %d", x);
       //emit sent_respond(buffer);
       outFile2 << endl;
    }

    outFile2.close();
    //delete [] marray;

    /*
    LP2[x_, y_, z_] := b2/2*(3*(z/Sqrt[x^2 + y^2 + z^2])^2 - 1);
    LP0[x_, y_, z_] := 1;
    sig = 1; v0 = 100;
    gf3[x_, y_, z_] := Exp[-((Sqrt[x^2 + y^2 + z^2] - v0)/sig)^2];
    sf1[x_, y_, z_] := gf3[x, y, z]*(LP0[x, y, z] + LP2[x, y, z]);

    Proj[x_, z_] := 1000*N[2*Sum[sf1[x, y, z], {y, 0, 255}]];

    Proj[x_, z_] := 1000*N[2*Sum[sf1[x, y, z], {y, 0, 255}]];
    For[j = -255, j < 256, j++,
     List1 = Table[Proj[i, j], {i, -255, 255}];
     List2 = Append[List2, List1];
     Print[j];
     ]
*/




}
