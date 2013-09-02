#ifndef MDATA_H
#define MDATA_H

#include "math.h"
//#include <gsl/gsl_matrix.h>

class MData
{
public:
    MData();
    ~MData();

    // Loading data and data access functions
    int Load_2D( char*, bool );  // Loading ASCII file 2D matrix
    //int Load_2D_new( char* );  // new version of the above but not finished also not required for now
    int Load_binary(char* , int, int, int, bool); // Load binary data
    int Load_histo( double*, double*, int, double, double, double, double );
    int Load_1D(char*); // loading normalisation time profile
    int GetsizeX();
    int GetsizeY();
    double Get_pix_value(int, int);
    void SetBackground();
    //double GetElement(int, int);
    double TotalCounts( double*, int, int, int );
    double GetSpectrCounts (int, int);

    // Data processing functions
    void transpose();
    void normalise_2darray( double*, double ); // normalisation on arbitrary number
    double find_max( double* ); // returns maximum from any double array of size = size_x*size_y
    double find_max( double*, int ); // returns maximum from any double array of any size
    void circle_mask( int, int, int ); // apply circle mask to the image
    void reload_data();  // copy initial data into processed data
    void Subtract_bkg(); // substrack background image
    void Add_bkg(double );   // Add two images and save in bkg file
    void Get_quarter3(int, int); //  Reflects the bottom left quarter of the image to all other quarters
    void Get_quarter1(int, int); // Reflects the top left quarter of the image to all other quarters
    void Get_quarter2(int, int); // Reflects the top right quarter of the image to all other quarters
    void Get_quarter4(int, int); // Reflects the bottom right quarter of the image to all other quarters

    void Symmetrise(int, int); // adds all quarters together
    void CartToPolar(int, int, int, int);   // convert cartesian 2D image to polar 2D image

    double Change_pix_value(int, int, double);



    // Simple data operations
    int Copy_2D_ini_to_processed();
    int Add_result_ini_to_ini();
    int Copy_2D_ini_to_result();

    // Saving result
    bool SavePES(char*);       // save photoelectron velocity and energy spectrum
    bool SaveProcMatrix(char* ); //save processed matrix
    bool SaveInvMatrix(char*); // save inverted matrix
    bool SaveAng(char*);       // save angular distribution

    // Main variables
    double      *array_2D_initial, *array_2D_processed, *array_2D_result;
    double      *array_2D_tot_sig, *array_2D_bkg;
    double      *polar_2D_image;
    double      *time_profile; // normailization of VMI images
    int         delay_steps; // number of time delays
    int         size_x, size_y;
    int         size_r, size_phi, scaling;

    double      *array_1D_spectrum, *array_coeffs;
    int         size_spectrum, nL;
    int         Nskip_lines; // amount of line to be skipped during loading the data frm file
    int         binary_depth;
    bool        t_flag, bkg_flag; //transpose flag

    // Additional variables
//
    //double      *array_2d;
    double       **input_2d_array, **output_2d_array; // double** used in pBasex
    double      max_element, min_element;
    //int         total_counts, useful_counts;


    //gsl_matrix *U,*V;

    void declare_matrices(); // allocate memory for **input_2d_array, **output_2d_array;
                             //and copy array_2d into input_2d_array
    void copy_f_to_d(double**); // copy any double** array to double* array_2D_result
    void copy_pes(double**); // copy photoelectron spectrum and beta parameters

    // constatns

    //static const double m_pi = 2*acos(0.0);

};

#endif // MDATA_H
