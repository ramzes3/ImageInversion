#ifndef MPOP_H
#define MPOP_H

//#include <QThread>
#include <QObject>

class MPop : public QObject //public QThread
{
    Q_OBJECT

public:
    MPop();
    ~MPop();

    int invert(double*, double*, int, int, int, int, int);

    double **data_ini, **data_inv, *data_pes;
    int size_x, size_y, size_pes;

    bool Data; // choose simulated (1) or raw data (0)

    // Preload arrays
    int load_arrays();

signals:
    void sent_respond(char* );
    void sent_error(char* );

private:

    /*
    INPUT PARAMETERS AND DEFINITIONS
    Note: The numbers in brackets represent the sections that each of the parameters are used in.
    */

    int i; // For loop Increment (1,2,3,4)
    int j; // For loop Increment (1,2,3,4)
    int rp; // radius (2,3,4)
    int qp; // polar pixel number (2,3,4)
    int nainc; // number of angular increments for a radius rp (2,3,4)
    float ainc; // angular increment of a radius rp (2,3,4)
    float xp; // polar x coordinate for a given rp and alpha (2)
    float yp; // polar y coordinate for a given rp and alpha (2)
    int xc; // x coordinate of cartesian pixel for which the centre of the polar pixel (xp,yp) is in (2)
    int yc; // y coordinate of cartesian pixel for which the centre of the polar pixel (xp,yp) is in (2)
    float xd; // x percentage of cartesian pixel (xc,yc) that the polar pixel (xp,yp) occupies (2)
    float yd; // y percentage of cartesian pixel (xc,yc) that the polar pixel (xp,yp) occupies (2)
    float pint; // intensity of the polar pixel (rp,qp) (2)
    float xa; // Calculation of the second order legendre polynomial (3)
    float xb; // Calculation of the forth order legendre polynomial (3)
    float y; // Intensity of the ira matrix at ira[rp][qp] (3)
    float det; // Determinant of Matrix A (3)
    float Na; // Absolute value of N (Beta[0]) (3)
    int R; // radius of polar pixel extracted from LUT (3,4)
    int pixels; // number of polar pixels (qp) at a given radius (rp) (4)
    float it; // intensity of pixel (qp) in expanded normalised LUT (4)
    float alpha; // alpha (4)
    float fact; // scaled angular distribution at pixel (R,qp) from LUT (4)

    int **PIXELS; int pixels_rows;
    float *AINC, **LUT;  int ainc_col, lut_rows, lut_cols;

};

#endif // MPOP_H
