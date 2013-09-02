#include "mdata.h"

#include <stdio.h>
#include "math.h"
#include <inttypes.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <QFileDialog>

using namespace std;

std::string changeExtension(std::string name, const std::string& ext)
{
    std::size_t extPos = name.rfind('.') ;

    if ( extPos != std::string::npos )
        name.resize(extPos) ;

    return name + ext;
}

MData::MData()
{
    //array_2d = NULL;
    Nskip_lines = 0;
    t_flag = false; bkg_flag = false;
    size_x = 512; size_y = 512; // default array dimensions
    array_2D_initial = NULL; array_2D_processed = NULL; array_2D_result = NULL;
    array_1D_spectrum = NULL; array_coeffs = NULL;
    array_2D_tot_sig = NULL; array_2D_bkg = NULL;
    time_profile = NULL; delay_steps = 0;

    polar_2D_image = NULL;
    input_2d_array = NULL; output_2d_array = NULL;
}

MData::~MData()
{
    //if(array_2d) delete [] array_2d;
    if(array_2D_initial)    delete [] array_2D_initial;
    if(array_2D_processed)  delete [] array_2D_processed;
    if(array_2D_result)     delete [] array_2D_result;
    if(array_1D_spectrum)   delete [] array_1D_spectrum;
    if(array_coeffs)        delete [] array_coeffs;
    if(array_2D_tot_sig)    delete [] array_2D_tot_sig;
    if(array_2D_bkg)        delete [] array_2D_bkg;
    if(polar_2D_image)      delete [] polar_2D_image;
    if(time_profile)        delete [] time_profile;

    if(input_2d_array)  for(int i=0;i<size_y;i++) delete [] input_2d_array[i];
    if(output_2d_array) for(int i=0;i<size_y;i++) delete [] output_2d_array[i];
    delete [] input_2d_array;
    delete [] output_2d_array;
}

/* int MData::Load_2D_new(char* filename)
{
    //FILE* pFile; pFile = NULL;
    //int nrows, ncol;
    string line;

    ifstream pFile (filename);
    if (!pFile.is_open()) return 0;
    while(getline(pFile, line)){
        ;
    }
    return 1;
}
*/

int MData::Load_2D(char* filename, bool transpose_flag)
{
    FILE* pFile; pFile = NULL;
    int nrows, ncol;
    float var;
    char buffer;
    //char buffer [100];
    pFile = fopen (filename , "r");
    if (pFile == NULL) {return 0;}

    nrows = 0; ncol = 0;
    do {  // calculate number of rows
      buffer = getc(pFile);
      if ( buffer == '\n' ) nrows++;
      //if ( ncol == 0 ) nrows++;
    } while (buffer != EOF);
    rewind (pFile); // bring posistion to the beginning of the file

    var = 0;
    while (var != Nskip_lines){  // skip first Nlines
      buffer = getc(pFile);
      if ( buffer == '\n' ) var++;
    }

    fpos_t pos;    
    do {  // calculate number of columns
      fgetpos (pFile,&pos);
      buffer = getc(pFile);
      if( buffer == '\t') continue;
      if( buffer == ' ') continue;
      if( buffer == '\n') break;
      if( buffer == '\r') break;
      fsetpos (pFile,&pos);
      //fseek( pFile, pos, SEEK_SET);
      buffer = fscanf(pFile, "%f", &var);
      ncol++;
      if( ncol > 10000 ) return 0;
    } while (1);
    //nrows = nrows / ncol;
    if(transpose_flag) {size_y = ncol; size_x = nrows - Nskip_lines;}
    else {size_x = ncol; size_y = nrows - Nskip_lines;}

    rewind (pFile);
    nrows = 0;
    while (nrows != Nskip_lines) {  // skip first Nskip_lines
      buffer = getc(pFile);
      if ( buffer == '\n' ) nrows++;
    }

    if(array_2D_initial) delete [] array_2D_initial;
    array_2D_initial = new double [size_x*size_y];

    if(array_2D_processed) delete [] array_2D_processed;
    array_2D_processed = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_initial[i] = 0; array_2D_processed[i] = 0;} //array_2D_processed[i] = 0;

    const double bkglevel = 0;
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           fscanf(pFile, "%f", &var);
           buffer = getc(pFile);
           if(transpose_flag){
               array_2D_initial[i + j*size_y] = (double)var - bkglevel;
               array_2D_processed[i + j*size_y] = array_2D_initial[i + j*size_y];
           }
           else{

               array_2D_initial[j + i*size_x] = (double)var - bkglevel;
               array_2D_processed[j + i*size_x] = array_2D_initial[j + i*size_x];
           }
       } 
       //buffer = getc(pFile);
    }
    fclose (pFile);    

    //if(array_2D_result) {delete [] array_2D_result; array_2D_result = NULL;}
   // array_2D_result = new double [size_x*size_y];

    return 1;
}

int MData::Load_binary(char* filename, int depth, int width, int height, bool transpose_flag){
/*
    QFile file(filename);
    if(!file.open(QIODevice::ReadOnly)){ return -1;};;
    QDataStream in(&file);    // read the data serialized from the file

    size_x = width; size_y = height;
    if(array_2D_initial) delete [] array_2D_initial;
    array_2D_initial = new double [size_x*size_y];

    if(array_2D_processed) delete [] array_2D_processed;
    array_2D_processed = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_initial[i] = 0; array_2D_processed[i] = 0;} //array_2D_processed[i] = 0;

    double a;
    for(int i=0; i < (size_x*size_y);i++){
         in >> a;           // extract element
         array_2D_initial[i] = a;
         array_2D_processed[i] = array_2D_initial[i];
    }

*/

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly)) return 0;
    QByteArray var = file.readAll();

    if(transpose_flag){size_y = width; size_x = height;}
    else{ size_x = width; size_y = height; }
    if(array_2D_initial) delete [] array_2D_initial;
    array_2D_initial = new double [size_x*size_y];

    if(array_2D_processed) delete [] array_2D_processed;
    array_2D_processed = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_initial[i] = 0; array_2D_processed[i] = 0;} //array_2D_processed[i] = 0;

    unsigned int *emptyarray; emptyarray = new unsigned int [size_x*size_y];
    //int32_t *emptyarray; emptyarray = new int32_t [size_x*size_y];

    memcpy(emptyarray, var.data(), size_x*size_y*sizeof(int));

    if(transpose_flag) {
        for (int i=0;i<size_y;i++){
           for (int j=0;j<size_x;j++){
               array_2D_initial[i + j*size_y] = (double) emptyarray[j + i*size_x] / 100.0;
               array_2D_processed[i + j*size_y] = array_2D_initial[i + j*size_y];
           }
        }
    }
    else {
            for(int i=0; i<size_x*size_y; i++){ array_2D_initial[i] = (double) emptyarray[i] / 100.0; array_2D_processed[i] = (double) emptyarray[i] / 100.0;}
        }


    //if(array_2D_result) {delete [] array_2D_result; array_2D_result = NULL;}
    //array_2D_result = new double [size_x*size_y];

    file.close();
    delete [] emptyarray;

    return 1;
}

 int MData::Load_1D(char* filename)
 {
     FILE* pFile; pFile = NULL;
     char buffer;
     float var;
     //char buffer [100];
     pFile = fopen (filename , "r");
     if (pFile == NULL) {return 0;}

     delay_steps = 0;
     do {  // calculate number of rows
       buffer = getc(pFile);
       if ( buffer == '\n' ) delay_steps++;
     } while (buffer != EOF);
     rewind (pFile); // bring posistion to the beginning of the file

     if(time_profile) delete [] time_profile;
     time_profile = new double [delay_steps];
     for(int i=0; i<delay_steps; i++){
         fscanf(pFile, "%f", &var);
         fscanf(pFile, "%f", &var);
         buffer = getc(pFile);
         time_profile[i] = (double)var;
     }
     return 1;

 }

int MData::Copy_2D_ini_to_result()
{
    if( !array_2D_initial ) return -2;
    if(array_2D_result) {delete [] array_2D_result; array_2D_result = NULL;}
    array_2D_result = new double [size_x*size_y];

    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_result[j + i*size_x] = array_2D_initial[j + i*size_x];

       }
    }
    return 1;
}

int MData::Copy_2D_ini_to_processed()
{    
    if( !array_2D_initial ) return -2;
    if( !array_2D_processed ) array_2D_processed = new double [size_x*size_y];
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_processed[j + i*size_x] = array_2D_initial[j + i*size_x];

       }
    }
    return 1;
}

int MData::Add_result_ini_to_ini()
{
    if( !array_2D_result ) return -1;
    if( !array_2D_initial ) return -2;
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_result[j + i*size_x] += array_2D_initial[j + i*size_x];
           array_2D_initial[j + i*size_x] = array_2D_result[j + i*size_x];
           array_2D_processed[j + i*size_x] = array_2D_initial[j + i*size_x];
       }
    }
    return 1;
}

int MData::Load_histo(double* array_x, double*array_y, int size, double min_x, double max_x, double min_y, double max_y)
{
    if(array_2D_initial) delete [] array_2D_initial;
    array_2D_initial = new double [size_x*size_y];

    if(array_2D_processed) delete [] array_2D_processed;
    array_2D_processed = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_initial[i] = 0; array_2D_processed[i] = 0; }

    int posx, posy;
    for(int i=0; i<size; i++){
        posx = (int) size_x * ( array_x[i] - min_x )/( max_x - min_x );
        posy = (int) size_y * ( array_y[i] - min_y )/( max_y - min_y );
        if( posx > size_x || posx < 0 ) continue;
        if( posy > size_y || posy < 0 ) continue;
        array_2D_initial[posy + size_y * posx]++;
        array_2D_processed[posy + size_y * posx]++;
    }
    return 1;
}

int MData::GetsizeX(){return size_x;}
int MData::GetsizeY(){return size_y;}
double MData::Get_pix_value(int posx, int posy) {
    if(!array_2D_processed) return -1;
    if( (posx<size_x)&&(posy<size_y)&&(posx>-1)&&(posy>-1)) return array_2D_processed[posx + posy * size_x];
    else return 0;
}

double MData::find_max(double* input_array){
    double max = input_array[0];
    int size = size_x*size_y;
    for (int i = 0; i<size; i++){
        if(max < input_array[i]) max = input_array[i];
    }
    return max;
}

double MData::find_max(double* input_array, int size){
    double max = input_array[0];
    for (int i = 0; i<size; i++){
        if(max < input_array[i]) max = input_array[i];
    }
    return max;
}

void MData::normalise_2darray(double *input_array, double norm)
{
    int size = size_x*size_y;
    double max = find_max(input_array);
    for (int i = 0; i<size; i++){ input_array[i] = norm * input_array[i] / max; }
}

void MData::circle_mask(int xc, int yc, int radius )
{
    int rad2 = radius*radius;
    if( !array_2D_processed ) return;
    for (int i = 0; i<size_y; i++){
        for (int j = 0; j<size_x; j++){
            if ( ((i-yc)*(i-yc) + (j-xc)*(j-xc)) > rad2 ) array_2D_processed[ j + i * size_x] = 0;
        }
    }
}

void MData::transpose()
{
    if( !array_2D_processed ) return;
    if(size_x != size_y) return;
    if(array_2D_result) delete [] array_2D_result;
    array_2D_result = new double [size_x*size_y];

    for (int i = 0; i<size_y; i++){
        for (int j = 0; j<size_x; j++){
            array_2D_result[ i + j * size_y] = array_2D_processed[ j + i * size_x];
        }
    }
    for (int i = 0; i<size_y; i++){
        for (int j = 0; j<size_x; j++){
            array_2D_processed[ j + i * size_y] = array_2D_result[ j + i * size_x];
            array_2D_initial[ j + i * size_y] = array_2D_processed[ j + i * size_y];
        }
    }
}

double MData::TotalCounts(double *input_array, int posx, int posy, int diameter)
{
    double counts = 0;
    for (int i = 0; i<size_y; i++){
        for (int j = 0; j<size_x; j++){
            if( sqrt( pow(j-posx,2) + pow(i-posy,2)) > diameter ) continue;
            counts += input_array[ j + i * size_x];
        }
    }
    return counts;
}

double MData::Change_pix_value(int posx, int posy, double value)
{
    if( !array_2D_initial ) return -1;
    if ( (posx + posy * size_x) < size_x * size_y) {array_2D_initial[ posx + posy * size_x] = value; return (int) value;}
    else return -1;
}

double MData::GetSpectrCounts(int limit1, int limit2)
{
    if(!array_1D_spectrum) return -1;
    if( limit1 > size_spectrum || limit1 < 0 || limit2 > size_spectrum || limit2 < 0) return -1;
    double sum1;
    sum1 = 0;
    for(int i=limit1; i<limit2; i++){
        sum1 += array_1D_spectrum[i];
    }
    return sum1;
}

void MData::CartToPolar(int x_c, int y_c, int r_dim, int phi_dim)
{
    if( !array_2D_processed ) return;
    if(polar_2D_image) delete [] polar_2D_image;
    // polar image has only quarter of the bins it goes from 0 to
    size_r = r_dim; size_phi = phi_dim;
    polar_2D_image = new double [r_dim * phi_dim];

    if(size_x < size_y) scaling = size_x / (2*r_dim);
    else scaling = size_y / (2*r_dim);
    //scaling = sqrt(size_x * size_y) / (2*r_dim);

    int x1[2], y1[2]; //nodes
    double coeflin[4]; // linear interpolation coefficients
    double r_var, phi_var, x_var, y_var; // intermidiate variable
    // linear interplation
    for (int i=0;i<size_r;i++){
       r_var = i * scaling;
       for (int j=0;j<size_phi;j++){
           phi_var = M_PI*j/(size_phi - 1); x_var = r_var*cos(phi_var); y_var = r_var*sin(phi_var);
           x1[0] = floor(x_var) + x_c;
           x1[1] = ceil(x_var) + x_c;
           y1[0] = floor(y_var) + y_c;
           y1[1] = ceil(y_var) + y_c;
           coeflin[0] = 1-(ceil(x_var) - x_var) - (ceil(y_var) - y_var) + (ceil(x_var) - x_var) * (ceil(y_var) - y_var);
           coeflin[1] = (ceil(x_var) - x_var) - (ceil(x_var) - x_var) * (ceil(y_var) - y_var);
           coeflin[2] = (ceil(y_var) - y_var) - (ceil(x_var) - x_var) * (ceil(y_var) - y_var);
           coeflin[3] = (ceil(x_var) - x_var) * (ceil(y_var) - y_var);
           polar_2D_image[i + j * size_r] = 0;
           for(int k = 0; k<2;k++){
                for(int l = 0; l<2;l++){
                    if ( (x1[k] + y1[l] * size_x) < size_x*size_y ) polar_2D_image[j + i * size_phi] += array_2D_processed[ x1[k] + y1[l] * size_x] * coeflin[k + l*2];
                    else polar_2D_image[j + i * size_phi] += 0;
                }
           }
       }
    }
}

void MData::declare_matrices()
{
    if(input_2d_array) {
            for(int i=0;i<size_y;i++) delete [] input_2d_array[i];
            delete [] input_2d_array;
    }
    if(output_2d_array) {
            for(int i=0;i<size_y;i++) delete [] output_2d_array[i];
            delete [] output_2d_array;
    }

    input_2d_array = new double *[size_y];
    for(int i=0;i<size_y;i++) input_2d_array[i] = new double [size_x];

    output_2d_array = new double *[size_y];
    for(int i=0;i<size_y;i++) output_2d_array[i] = new double [size_x];

    if(array_2D_initial){
        for (int i=0;i<size_y;i++){
           for (int j=0;j<size_x;j++){
               input_2d_array[i][j] = (double) array_2D_processed[j + i*size_x];
               output_2d_array[i][j] = 0.0;
           }
       }
    }
}

void MData::copy_f_to_d(double** input_array)
{
    if(!input_array){return;}
    //if(!output_array){return;}

    if(array_2D_result) delete [] array_2D_result; array_2D_result = NULL;
    array_2D_result = new double [size_x*size_y];


    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_result[j + i*size_x] = (double) input_array[i][j];
       }
   }
}

void MData::copy_pes(double** input_array)
{
    if(!input_array) return;
    if(array_1D_spectrum) delete [] array_1D_spectrum; array_1D_spectrum = NULL;
    if(array_coeffs) delete [] array_coeffs; array_coeffs = NULL;

    array_1D_spectrum = new double [size_spectrum];
    array_coeffs = new double [size_spectrum * nL];


    for(int i=0; i<size_spectrum; i++){
        array_1D_spectrum[i] = input_array[0][i];
        for(int j=0; j<nL; j++){
             array_coeffs[j + i*nL] = input_array[j+1][i];
        }
    }

}

void MData::reload_data()
{
    if(!array_2D_initial) return;

    if(array_2D_processed) delete [] array_2D_processed;
    array_2D_processed = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_processed[i] = 0; }

    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_processed[j + i*size_x] = array_2D_initial[j + i*size_x];
       }
    }
}

void MData::SetBackground()
{
    if(!array_2D_initial) return;

    if(array_2D_bkg) delete [] array_2D_bkg;
    array_2D_bkg = new double [size_x*size_y];
    for(int i=0; i<size_x*size_y; i++){ array_2D_bkg[i] = 0; }

    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_bkg[j + i*size_x] = array_2D_initial[j + i*size_x];
       }
    }
    bkg_flag = true;
}

void MData::Subtract_bkg()
{
    if(!array_2D_initial) return;
    if(!array_2D_bkg) return;

    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_initial[j + i*size_x] = array_2D_initial[j + i*size_x] - array_2D_bkg[j + i*size_x];
           //array_2D_processed[j + i*size_x] = array_2D_initial[j + i*size_x];
       }
    }
}

void MData::Add_bkg(double division)
{
    if(!array_2D_initial) return;
    if(!array_2D_bkg) return;

    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_bkg[j + i*size_x] = (array_2D_initial[j + i*size_x] + array_2D_bkg[j + i*size_x])/division;
       }
    }
}

// copy the 1st (top left) quarter into all other
void MData::Get_quarter1(int x_c, int y_c)
{
    if(!array_2D_processed) return;
    if(x_c > size_x-1) return;
    if(y_c > size_y-1) return;

    // copy the 1st quarter into square number 2
    for (int i=0;i<y_c+1;i++){
       for (int j=x_c+1;j<size_x;j++){
           if( ((2*x_c-j) > size_x-1) || ((2*x_c-j) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + i*size_x];
       }
    }

    // copy 1st quarter into square number 3
    for (int i=y_c+1;i<size_y;i++){
       for (int j=0;j<x_c;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*y_c-i) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[j + (2*y_c-i)*size_x];
       }
    }

    // copy 1st quarter intosquare number 4
    for (int i=y_c+1;i<size_y;i++){
       for (int j=x_c;j<size_x;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*x_c-j) > size_x-1)) array_2D_processed[j + i*size_x] = 0;
           else {
               if( ((2*x_c-j) < 0)|| ((2*y_c-i) < 0))array_2D_processed[j + i*size_x] = 0;
               else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + (2*y_c-i)*size_x];
           }
       }
    }
}

// copy the 2nd (top right) quarter into all others
void MData::Get_quarter2(int x_c, int y_c)
{
    if(!array_2D_processed) return;
    if(x_c > size_x-1) return;
    if(y_c > size_y-1) return;

    // copy the 2nd quarter into square number 1
    for (int i=0;i<y_c+1;i++){  // copy
       for (int j=0;j<x_c+1;j++){
           if( ((2*x_c-j) > size_x-1) || ((2*x_c-j) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + i*size_x];
       }
    }


    // copy 2nd quarter into square number 3
    for (int i=y_c+1;i<size_y;i++){
       for (int j=0;j<x_c;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*x_c-j) > size_x-1)) array_2D_processed[j + i*size_x] = 0;
           else {
               if( ((2*x_c-j) < 0)|| ((2*y_c-i) < 0))array_2D_processed[j + i*size_x] = 0;
               else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + (2*y_c-i)*size_x];
           }
       }
    }

    // copy 2nd quarter intosquare number 4
    for (int i=y_c+1;i<size_y;i++){
       for (int j=x_c;j<size_x;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*y_c-i) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[j + (2*y_c-i)*size_x];
       }
    }
}

// copy the 3rd (bottom left) quarter into all other
void MData::Get_quarter3(int x_c, int y_c)
{
    if(!array_2D_processed) return;
    if(x_c > size_x-1) return;
    if(y_c > size_y-1) return;

    for (int i=0;i<y_c;i++){  // copy 3rd into 1st
       for (int j=0;j<x_c+1;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*y_c-i) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[j + (2*y_c-i)*size_x];
       }
    }

    for (int i=y_c;i<size_y;i++){ //copy 3rd into 4th
       for (int j=x_c+1;j<size_x;j++){
           if( ((2*x_c-j) > size_x-1) || ((2*x_c-j) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[2*x_c-j + i*size_x];
       }
    }

    for (int i=0;i<y_c;i++){ //copy 3rd into 2nd
       for (int j=x_c+1;j<size_x;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*x_c-j) > size_x-1)) array_2D_processed[j + i*size_x] = 0;
           else {
               if( ((2*x_c-j) < 0)|| ((2*y_c-i) < 0))array_2D_processed[j + i*size_x] = 0;
               else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + (2*y_c-i)*size_x];
           }
       }
    }
}

// copy the 4th (top left) quarter into all other
void MData::Get_quarter4(int x_c, int y_c)
{
    if(!array_2D_processed) return;
    if(x_c > size_x-1) return;
    if(y_c > size_y-1) return;


    // copy the 4th quarter into square number 1

    for (int i=0;i<y_c+1;i++){
       for (int j=0;j<x_c+1;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*x_c-j) > size_x-1)) array_2D_processed[j + i*size_x] = 0;
           else {
               if( ((2*x_c-j) < 0)|| ((2*y_c-i) < 0))array_2D_processed[j + i*size_x] = 0;
               else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + (2*y_c-i)*size_x];
           }
       }
    }

    // copy the 4th quarter into square number 2
    for (int i=0;i<y_c+1;i++){
       for (int j=x_c+1;j<size_x;j++){
           if( ((2*y_c-i) > size_y-1) || ((2*y_c-i) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[j + (2*y_c-i)*size_x];
       }
    }

    // copy 4th quarter intosquare number 3
    for (int i=y_c+1;i<size_y;i++){
       for (int j=0;j<x_c;j++){
           if( ((2*x_c-j) > size_x-1) || ((2*x_c-j) < 0)) array_2D_processed[j + i*size_x] = 0;
           else array_2D_processed[j + i*size_x] = array_2D_processed[(2*x_c-j) + i*size_x];
       }
    }
}

void MData::Symmetrise(int x_c, int y_c){
    if(!array_2D_processed) return;
    if(!array_2D_initial) return;
    if(array_2D_result) {delete [] array_2D_result; array_2D_result = NULL;}
    array_2D_result = new double [size_x*size_y];

    if(x_c > size_x-1) return;
    if(y_c > size_y-1) return;
    // declare result matrix
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_result[j + i*size_x] = array_2D_processed[j + i*size_x];
       }
    }

   // flip image over y coordinate and add to original
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           if ( ((-i+2*y_c)>-1) && ((-i+2*y_c)<size_y)) array_2D_result[j + i*size_x] += array_2D_processed[j + (-i+2*y_c)*size_x];
       }
    }
    // copy result matrix
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           array_2D_processed[j + i*size_x] = array_2D_result[j + i*size_x];
       }
    }
   // Flip image over x coordinate and add to already flipped over y coordinate
    for (int i=0;i<size_y;i++){
       for (int j=0;j<size_x;j++){
           if ( ((-j+2*x_c)>-1) && ((-j+2*x_c)<size_x) ) array_2D_processed[j + i*size_x] += array_2D_result[(-j+2*x_c) + (i)*size_x];
       }
    }

}

bool MData::SaveProcMatrix(char* filein)
{
    if(!array_2D_processed) {return 0;}
    //char fileout[255];
    //sprintf(fileout,"%s.%s", filein, "proc");
    string combined_str = changeExtension(filein, "_proc.dat" );
    ofstream outFile1(combined_str.c_str(), ios::out);
    if(!outFile1){ return 0; }

    for (int i=0;i<GetsizeY();i++){
       for (int j=0;j<GetsizeX();j++){
           outFile1 << "\t" << array_2D_processed[j + i*GetsizeX()];
       }
       outFile1 << endl;
    }
    outFile1.close();
    return 1;
}

bool MData::SaveInvMatrix(char* filein)
{
    if(!array_2D_result) {return 0;}
    //char fileout[255];
    //sprintf(fileout,"%s.%s", filein, "inv");
    string combined_str = changeExtension(filein, "_inv.dat" );
    ofstream outFile1(combined_str.c_str(), ios::out);
    if(!outFile1){ return 0; }

    for (int i=0;i<GetsizeY();i++){
       for (int j=0;j<GetsizeX();j++){
           outFile1 << "\t" << array_2D_result[j + i*GetsizeX()];
       }
       outFile1 << endl;
    }
    outFile1.close();
    return 1;
}

bool MData::SavePES(char* filein)
{
    if(!array_1D_spectrum) {return 0;}
    //char fileout[255];
    //sprintf(fileout,"%s.%s", filein, "pes");
    string combined_str = changeExtension(filein, "_pes.dat" );
    ofstream outFile1(combined_str.c_str(), ios::out);
    if(!outFile1){ return 0; }
    outFile1 << "px\t" << "I\t" << "px^2\t" << "I/px" << endl;
    outFile1 << 0 << "\t" << array_1D_spectrum[0] << "\t0\t" <<  "0" << endl;
    for(int i=1; i<size_spectrum; i++){
        outFile1 << i << "\t" << array_1D_spectrum[i] << "\t" << i*i << "\t" << array_1D_spectrum[i]/(float) i << endl;
    }
    outFile1.close();
    return 1;
}

bool MData::SaveAng(char* filein)
{
    if(!array_coeffs) {return 0;}
    //char fileout[255];
    //sprintf(fileout,"%s.%s", filein, "ang");
    string combined_str = changeExtension(filein, "_ang.dat" );
    ofstream outFile1(combined_str.c_str(), ios::out);
    if(!outFile1){ return 0; }

    for(int i=0; i<size_spectrum; i++){
        outFile1 << i;
        for(int j=0; j<nL; j++){
             outFile1 << "\t" << array_coeffs[j + i*nL];
        }
        outFile1 << endl;
    }
    outFile1.close();
    return 1;
}
