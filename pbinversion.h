#ifndef PBINVERSION_H
#define PBINVERSION_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
//#include <gsl/gsl_linalg.h>

#include <QThread>

class PBInversion : public QThread
{
    Q_OBJECT

public:
    PBInversion();
    ~PBInversion();

    void polar_a(int nl,int odd); // the function set variables nl_v = nl; odd_v = odd;
    void polar_b(double **data,double **Imrz,int xc,int yc,int dr,int nl,int odd); // main function of the image inversion
    double theta_f(double x,double y); // Obtains polar angle theta from the x, y coordinates
    double cubic(double); //
    double ldist(double X, gsl_vector *coeff, double *pl, int npl, int odd, int ik);
    void load_matrices(int l,int odd);
    void generateTestVMI(const int, double, double, double); // size x and y, gauss FWHM, centre position, beta2


    void write_forward(int l,int odd); // write basis function
    double model_f(double r,void *par); //
    double model_c(int k,int l,double x,double z,gsl_integration_workspace *w); //
    double lpoly(double X, int npl); //
    int save_gsl_matrix(const gsl_matrix * , char*);
    int read_gsl_matrix(const gsl_matrix * , char*);
    int save_gsl_vector(const gsl_vector * , char*);
    int read_gsl_vector(const gsl_vector * , char*);

    void set_size(int, int);
    void set_parameters(int, int);
    void Invert_Polar(double*, double*, int ,int ,int ,int ,bool ); // main function of the image inversion

    int  NR;  //  = 256; // Radial binning
    int  NTH; // = 256; // Angular binning (180 degrees)
    int  NK;  //  = 128; // Number of basis functions

    int size_r, size_theta, size_basis; // radial bining, angular binning (180 degrees) and number of basis functions

    int ncolumns, nrows, ang_size;
    bool m_isloaded;

    gsl_matrix *U,*V;
    gsl_vector *S;
    double **ang;
    double **data_v; double **Imrz_v;
    int xc_v, yc_v, dr_v, nl_v, odd_v;

    static double Wrapper_To_Call_model_f(double r,void *par);    

signals:

    void sent_respond(const char* );
    void sent_error(const char* );

protected:
    void run();

};

 //typedef  double (PBInversion::*model_f)(double r,void *par);

#endif // PBINVERSION_H
