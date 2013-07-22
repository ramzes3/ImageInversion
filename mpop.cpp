#include "mpop.h"

#include <math.h>
#include "time.h"

#include <fstream>

MPop::MPop()
{
    size_x = 512; size_y = 512; size_pes = 256;
    data_ini = NULL; data_inv = NULL; data_pes = NULL;
    PIXELS = NULL; AINC = NULL; LUT = NULL;
    Data = 0;

    load_arrays();
}

MPop::~MPop(){
    if(data_ini) {
        for(int i=0;i<size_y;i++) delete [] data_ini[i];
        delete [] data_ini;
    }
    if(data_inv) {
        for(int i=0;i<size_y;i++) delete [] data_inv[i];
        delete [] data_inv;
    }

    if(data_pes) delete data_pes;

    if(LUT) {
        for(int i=0;i<lut_rows;i++) { delete [] LUT[i]; }
        delete [] LUT;
    }
    if(PIXELS) {
        for(int i=0;i<pixels_rows;i++) { delete [] PIXELS[i]; }
        delete [] PIXELS;
    }
    if(AINC) delete [] AINC;

}

int MPop::invert(double *data_input, double *data_output, int s_x, int s_y, int rad, int p_x, int p_y)
{
    if(!PIXELS){emit sent_respond( "there is no PIXELS matrix loaded" ); return -1;}
    if(!AINC){emit sent_respond( "there is no AINC matrix loaded" ); return -1;}
    if(!data_input) {emit sent_respond( "there is no input data" ); return -1;}
    if(!data_output) {emit sent_respond( "output array hasnt been initialised" ); return -1;}
    size_x = s_x; size_y = s_y;

    // Initialising 2D arrays
    if(data_ini) {
        for(int i=0;i<size_y;i++) delete [] data_ini[i];
        delete [] data_ini;
    }
    if(data_inv) {
        for(int i=0;i<size_y;i++) delete [] data_inv[i];
        delete [] data_inv;
    }

    data_ini = new double *[size_x];
    for(i=0;i<size_x;i++) data_ini[i] = new double [size_y];
    data_inv = new double *[size_x];
    for(i=0;i<size_x;i++) data_inv[i] = new double [size_y];

    // fill initial array
    for (i=0;i<size_y;i++){
       for (j=0;j<size_x;j++){
           data_ini[j][i] = data_input[j + i*size_x];
           data_inv[j][i] = 0.0;
       }
   }

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    (1) CARTESIAN IMAGE CENTRING & 4 WAY FOLDING
    (zzx,zzy) is defined as the centre of the image
    The centre of a Cartesian pixel is defined as the centre of the pixel (NOT the bottom left vertex)

    Note: Only even dimension sizes can be used!
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    // find out size for a4 matrix
    int zz, zzx, zzy;
    if( (s_x/2) > p_x) zzx = p_x;
    else zzx = s_x - p_x;
    if( (s_y/2) > p_y) zzy = p_y;
    else zzy = s_y - p_y;
    if( zzx > zzy ) zz = zzy;
    else zz = zzx;
    // declare a4 matrix
    const int a4_size = zz;
    double a4[a4_size][a4_size], ira[a4_size][a4_size], iraraw[a4_size][a4_size], iradecon[a4_size][a4_size];    

    a4[0][0]=data_ini[p_x][p_y]*4;

    for(i=1;i<=zz-1;i++) {

       a4[i][0]=(data_ini[p_x+i][p_y]+data_ini[p_x-i][p_y])*2;

       for(j=1;j<=zz-1;j++) {

          a4[i][j]=data_ini[p_x+i][p_y+j]+data_ini[p_x-i][p_y+j]+data_ini[p_x+i][p_y-j]+data_ini[p_x-i][p_y-j];
          a4[0][j]=(data_ini[p_x][p_y+j]+data_ini[p_x][p_y-j])*2;

       }
    }

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    (2) CONVERSION FROM CARTESIAN TO POLAR COORDINATES: rp,alpha
    alpha is the angle and qp is the pixel designation in alpha. qp=0 -> alpha=0
    This conversion takes into account that each of the Cartesian pixels have a dimension 1x1
    This is also true for the polar pixels generated.
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    ira[0][0]=a4[0][0];

    for(rp=1;rp<=zz-1;rp++) {

    /*
    INFO:
    For a given radius (rp) the angular increment (ainc) between each polar pixel is calculated by first calculating the number of angular increments for that radius (nainc). Nainc will be given by the number of polar pixels (The number of polar pixels is calculated using (rp+1)*(pi/2) which is then rounded down using floor) at that radius minus one, as there are one less angular increments than there are polar pixels.
    */

       nainc=PIXELS[rp][1];
       ainc=AINC[rp];

       for(qp=0;qp<=nainc;qp++) {

    /*
    INFO:
    The angle at which a polar pixel lies (alpha) for a given radius (rp) is calculated using ainc*qp. The x and y coordinates of the polar pixel (xp and yp respectively) are then determined using xp=rp*sin(alpha) and yp=rp*cos(alpha) from basic trigonometry.

    The x and y coordinates of the Cartesian pixel (xc and yc respectively) for which the centre of the polar pixel (xp,yp) lies in is then determined by rounding xp and yp to the nearest integer (using int) to give xc and yc respectively. (Note: xc and yc are integers but xp and yp are floating point numbers).

    The absolute value of the difference between the the polar and Cartesian coordinates (i.e. abs(xc-xp) and abs(yc-yp)) is then subtracted from 1 (1 because the size of the Cartesian pixels is 1x1). This gives xd and yd respectively.
    */

          xp=rp*sin(ainc*qp);
          yp=rp*cos(ainc*qp);

          xc=int(xp);
          yc=int(yp);

          xd=1-fabs(xc-xp);
          yd=1-fabs(yc-yp);

    /*
    INFO:
    xd and yd multiplied together (i.e. xd*yd) provides the percentage overlap of the polar pixel (xp,yp) with the Cartesian pixel (xc,yc). The intensity of the cartesian pixel (xc,yc) in the a4 array (i.e. a4[xc][yc]) is then multiplied by this percentage overlap. This provides the first of four values to be added together to determine the polar pixel intensity (pint).

    If statements are used to determined which other three Cartesian pixels the polar pixel overlaps with. Firstly, the percentage of intensity that needs to be taken from each of these Cartesian pixels needs to be determined. Three moves in the Cartesian frame are possible and they each correspond to these percentages:

    1. x coord=stay the same, y coord=change, percentage=xd*(1-yd)
    2. x coord=change, y coord=stay the same, percentage=(1-xd)*yd
    3. x coord=change, y coord= change, precentage=(1-xd)*(1-yd)

    To determine the other three Cartesian pixels, a comparison between xp and xc and yp and yc is carried out. If polar>cartesian then +1 to the coordinate, if polar<cartesian then -1 from that coordinate.

    All four intensity values are summed to give the intensity of the polar pixel (pint).
    */

          if(xp>=xc && yp>=yc) {
             pint=((xd*yd)*a4[xc][yc])+((xd*(1-yd))*a4[xc][yc+1])+(((1-xd)*yd)*a4[xc+1][yc])+(((1-xd)*(1-yd))*a4[xc+1][yc+1]);
          }

          if(xp>=xc && yp<=yc) {
             pint=((xd*yd)*a4[xc][yc])+((xd*(1-yd))*a4[xc][yc-1])+(((1-xd)*yd)*a4[xc+1][yc])+(((1-xd)*(1-yd))*a4[xc+1][yc-1]);
          }

          if(xp<=xc && yp>=yc) {
             pint=((xd*yd)*a4[xc][yc])+((xd*(1-yd))*a4[xc][yc+1])+(((1-xd)*yd)*a4[xc-1][yc])+(((1-xd)*(1-yd))*a4[xc-1][yc+1]);
          }

          if(xp<=xc && yp<=yc) {
             pint=((xd*yd)*a4[xc][yc])+((xd*(1-yd))*a4[xc][yc-1])+(((1-xd)*yd)*a4[xc-1][yc])+(((1-xd)*(1-yd))*a4[xc-1][yc-1]);
          }

    /*
    INFO: The polar pixel intensity is placed in the appropriate polar pixel (i.e. ira[rp][qp]) in the polar matrix (ira).
    */

          ira[rp][qp]=ira[rp][qp]+pint;

       }
    }

    /*
    INFO: Saving of the polar image ira into iraraw
    */

    for(rp=0;rp<=zz-1;rp++) {

       nainc=PIXELS[rp][1];

       for(qp=0;qp<=nainc;qp++) {

          iraraw[rp][qp]=ira[rp][qp];

       }
    }

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    GENERATION OF THE CONVOLUTED PHOTOELECTRON SPECTRUM, PES I(R)
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    double PESR[a4_size], PESI[a4_size];
    for(i=0;i<=zz-1;i++) {

    nainc=PIXELS[i][1];
    PESR[i]=i;

       for(j=0;j<=nainc;j++) {

          PESI[i]=PESI[i]+ira[i][j];

       }
    }

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    PART II - POLAR ONION PEELING ALGORITHM

    PART II of the code carries out the deconvolution of the image via an onion peeling method. This is carried out via removal of the phi (f) contribution from the raw image (ira) at every radius (rp) to obtain a "slice" through the centre of the original 3D distribution of photofragments projected onto 2D detector surface i.e at f=0 only.

    This involves repeating steps (3) and (4) for every radius (rp). This is done by starting at the outermost radius (rp=zz-1) and incrementally decreasing rp until rp=0 is reached.
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    /*
    INFO:
    The deconvolution is carried out using an overall rp "for loop" that runs for every radius (rp) from zz-1 to 0.
    */
    double ratio;
    double B[3], A[4][4], Ain[4][4], Beta[3], Beta2[a4_size], Beta4[a4_size], PESId[a4_size], PESRd[a4_size];
    for(rp=zz-1;rp>=0;rp--) {

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    (3) OUTER RING FIT OF EXPERIMENTAL IMAGE (IRA) AT rp USING LINEAR REGRESSION

    Note: theta=alpha, rp=R
    Note: nainc+1 is the number of polar pixels at a given value of rp (i.e. for linear regression, the number of data points)

    For information on how the linear regression technique works in more details see:

    "Data Reduction and Error Analysis for Physical Sciences" 3rd Edition,
    by Philip R. Bevington & D. Keith Robinson

    Chapter 7 - Least Squares Fit to a Polynomial
    p122 - Matrix Solution(7.2)

    (more specifically for Legendre Polynomials, p132)
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    /*
    INFO:
    Zeroing of Matrices A and B before the outer ring fit begins for a given R
    */

    B[0]=0;
    B[1]=0;
    B[2]=0;

    for(i=0;i<=3;i++) {
       for(j=0;j<=3;j++) {
       A[i][j]=0;
       }
    }

    /*
    INFO:
    Calculation of Matrices A and B
    */

    nainc=PIXELS[rp][1];

    for(qp=0;qp<=nainc;qp++) {

       y=ira[rp][qp];
       xa=0.5*((3*pow(cos(AINC[rp]*qp),2))-1);
       xb=(1/8)*((35*pow(cos(AINC[rp]*qp),4))-(30*pow(cos(AINC[rp]*qp),2))+3);

       B[0]=B[0]+y;
       B[1]=B[1]+(y*xa);
       B[2]=B[2]+(y*xb);

       A[0][0]=A[0][0]+1;
       A[0][1]=A[0][1]+xa;
       A[0][2]=A[0][2]+xb;
       A[1][1]=A[1][1]+pow(xa,2);
       A[1][2]=A[1][2]+(xa*xb);
       A[2][2]=A[2][2]+pow(xb,2);

    }

    A[2][1]=A[1][2];
    A[2][0]=A[0][2];
    A[1][0]=A[0][1];

    /*
    INFO:
    Inversion of Matrix A
    */

    Ain[0][0]=(A[1][1]*A[2][2])-(A[1][2]*A[2][1]);
    Ain[0][1]=((A[1][0]*A[2][2])-(A[1][2]*A[2][0]))*-1;
    Ain[0][2]=(A[1][0]*A[2][1])-(A[1][1]*A[2][0]);
    Ain[1][1]=(A[0][0]*A[2][2])-(A[0][2]*A[2][0]);
    Ain[1][2]=((A[0][0]*A[2][1])-(A[0][1]*A[2][0]))*-1;
    Ain[2][2]=(A[0][0]*A[1][1])-(A[0][1]*A[1][0]);

    Ain[1][0]=Ain[0][1];
    Ain[2][0]=Ain[0][2];
    Ain[2][1]=Ain[1][2];

    /*
    INFO:
    Calculation of the Determinant of Matrix A
    */

    det=(A[0][0]*Ain[0][0])-(A[0][1]*(-Ain[0][1]))+(A[0][2]*Ain[0][2]);

    for(i=0;i<=3;i++) {
       for(j=0;j<=3;j++) {

       Ain[i][j]=Ain[i][j]/det;

       }
    }

    /*
    INFO:
    Calculation of the Beta co-efficient Matrix via mulitplication of the Matrix Ain and Matrix B:
    Beta[0] = N (Intensity Factor)
    Beta[1] = Beta2
    Beta[2] = Beta4
    */

    Beta[0]=(B[0]*Ain[0][0])+(B[1]*Ain[1][0])+(B[2]*Ain[2][0]);
    Beta[1]=((B[0]*Ain[0][1])+(B[1]*Ain[1][1])+(B[2]*Ain[2][1]))/Beta[0];
    Beta[2]=((B[0]*Ain[0][2])+(B[1]*Ain[1][2])+(B[2]*Ain[2][2]))/Beta[0];

    /*
    INFO:
    Generation of Beta Spectra for Beta2 and Beta4
    */

    Beta2[rp]=Beta[1];
    Beta4[rp]=Beta[2];

    Na=fabs(Beta[0]);

    if(Na==0) {

    PESId[rp]=0;
    PESRd[rp]=rp;

    continue;

    }


    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    (4) DECONVOLUTION

    (a) EXPAND NORMALISED LUT TO OBTAIN SIMULATED IMAGE
    (b) SELECT OUTER RING AT RP TO BE USED FOR IRADECON (EXPERIMENT OR SIMULATED)
    (c) SUBTRACT SIMULATED IMAGE FROM EXPERIMENTAL IMAGE FOR EACH PIXEL (qp) AT EACH RADIUS (R)
    (d) GENERATION OF DECONVOLUTED PES

    Note: For this simulated image and fit the values of x and y correspond to equivalent Cartesian coordinates
    Note: Cos(q)**2 = (R/rp)**2 * Cos(a)**2, where q=theta & a=alpha. This is due to alpha not being equal to theta at all values of phi (f).
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */

    PESRd[rp]=rp;

    for(R=0;R<=rp;R++) {

       pixels=PIXELS[R][0];
       it=(LUT[rp][R]*((nainc+1)/pixels));

       for(qp=0;qp<=pixels-1;qp++) {

          alpha=AINC[R]*qp;

          if(rp == 0) ratio = 0;
          else ratio = R/rp;
          fact=Beta[0]*it*(1+((Beta[1]*0.5)*(3*pow(ratio,2)*pow(cos(alpha),2)-1))+((Beta[2]*(1/8))*((35*pow(ratio,4)*pow(cos(alpha),4))-(30*pow(ratio,2)*pow(cos(alpha),2))+3)));

          if(Data==1 && R==rp) {  // simulated data

             PESId[R]=PESId[R]+fact*sqrt(R)*sin(alpha);
             iradecon[R][qp]=fact/sqrt(R);

          }

          else if(Data==0 && R==rp) { // raw data

             PESId[R]=PESId[R]+ira[R][qp]*sqrt(R)*sin(alpha);
             iradecon[R][qp]=ira[R][qp]/sqrt(R);

          }

          ira[R][qp]=ira[R][qp]-fact;

       }
    }

    }

    /*
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    END
    +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--
    */
    if(data_pes) delete data_pes;
    data_pes = new double[zz];
    // copy the inverted data
    size_pes = zz;
    char buffer[255];
    for(rp=zz-1;rp>=0;rp--) {
        data_pes[rp] = PESId[rp];
        sprintf(buffer, "value: %f", PESId[rp]);
        emit sent_error(buffer);

    }

    return 0;

}

int MPop::load_arrays()
{
    // prepare arrays
    if(LUT) {
        for(int i=0;i<lut_rows;i++) { delete [] LUT[i]; }
        delete [] LUT;
    }
    if(PIXELS) {
        for(int i=0;i<pixels_rows;i++) { delete [] PIXELS[i]; }
        delete [] PIXELS;
    }
    if(AINC) delete [] AINC;

    //  get file
    FILE* pFile; pFile = NULL;
    char filename[255], buffer[1024], symbol;

    int ncol, nrows; ncol = 0; nrows = 0;  // ncol for AINC file and nrows for Pixels file
    float var1; var1 = 0;
    int var2, error_code; var2 = 0;

    // Load AINC array
    sprintf(filename,"Ainc_LUT.txt");
    pFile = fopen (filename , "r");
    if (pFile == NULL) {return -1;}

    do {  // calculate number of columns
        error_code = fscanf(pFile, "%f", &var1);
        if(error_code == EOF) break;
        ncol++;
    }while(1);
    rewind (pFile); // bring posistion to the beginning of the file

    AINC = new float [ncol];  // allocate memory for AINC
    for(int i=0;i<ncol;i++) {  // read data to AINC array
        error_code = fscanf(pFile, "%f", &var1);
        AINC[i] = var1;
        //if(buffer == EOF) break;
    }
    fclose (pFile);

    // Load Pixels array
    sprintf(filename,"Pixels_LUT.txt");
    pFile = fopen (filename , "r");
    if (pFile == NULL) {return -1;}

    do {  // calculate number of rows
      error_code = fscanf(pFile, "%d", &var2);
      error_code = fscanf(pFile, "%d", &var2);
      nrows++;
    } while (error_code != EOF);
    rewind (pFile); // bring posistion to the beginning of the file

    PIXELS = new int *[nrows]; // Allocate memory for PIXELS
    for(int i=0;i<nrows;i++) {  PIXELS[i] = new int[2]; }

    for(int i=0;i<ncol;i++) {  // read data to PIXELS array
        error_code = fscanf(pFile, "%f", &var1);
        PIXELS[i][0] = var2;
        error_code = fscanf(pFile, "%f", &var1);
        PIXELS[i][1] = var1;
        //if(buffer == EOF) break;
    }
    fclose (pFile);

    pixels_rows = nrows;
    ainc_col = ncol;

// Next lines do work yet
    // Load LUT array
    sprintf(filename,"Basis_Set_Lin.txt");
    pFile = fopen (filename , "r");
    if (pFile == NULL) {return -1;}

    lut_rows = 0; lut_cols = 0;
    do {  // calculate number of rows
      symbol = getc(pFile);
      if ( (symbol == '\n') || (symbol == '\r') ) lut_rows++;
      //if ( ncol == 0 ) nrows++;
    } while (symbol != EOF);
    rewind (pFile); // bring posistion to the beginning of the file

    fpos_t pos;
    do {  // calculate number of columns
      fgetpos (pFile,&pos);
      symbol = getc(pFile);
      if( symbol == '\t') continue;
      if( symbol == ' ') continue;
      if( symbol == '\n') break;
      if( symbol == '\r') break;
      fsetpos (pFile,&pos);
      //fseek( pFile, pos, SEEK_SET);
      symbol = fscanf(pFile, "%f", &var1);
      lut_cols++;
      if( lut_cols > 10000 ) return 0;
    } while (1);

    LUT = new float *[lut_rows]; // Allocate memory for LUT
    for(int i=0;i<lut_rows;i++) {  LUT[i] = new float[lut_cols]; }

    for(int i=0;i<lut_rows;i++) {  // read data to LUT array
        for(int j=0;j<lut_cols;j++) {  // read data to LUT array
            error_code = fscanf(pFile, "%f", &var1);
            LUT[i][j] = var1;
        //if(buffer == EOF) break;
        }
    }
    fclose (pFile);
    return lut_rows;
}
