#ifndef COBOLDDATA_H
#define COBOLDDATA_H


#define NUM_CHANNELS 32
#define NUM_IONS 16

#include "Read_lmf/LMF_IO.h"

class LMF_IO;

class Cobolddata
{
public:
    Cobolddata(const char *LMF_Filename, int regime);
    ~Cobolddata();

    int   read_lmf_size(const char*);
    int   read_lmf(const char*, int);
    void  zeros(double*, int);

    //int                 iTDC[NUM_CHANNELS][NUM_IONS];
    //double      	dTDC[NUM_CHANNELS][NUM_IONS];
    unsigned int	number_of_hits[NUM_CHANNELS];

    double*       data_array; // structure: number of events -> number of channels -> number of hits
    int           array_size;

    double*       array_x;
    double*       array_y;
    int*          hits_histo[NUM_IONS];  // distribution of hits
    double        points_var[10];
    int           stat_counts;
    double*       array_tofx;  // not in use
    double*       array_tofy;  // not in use
    double*       array_corx;  // not in use
    double*       array_cory;  // not in use

    //double*       ch1; // channel  1
    //double*       ch2; // channel  2
    //double*       ch3; // channel  3
    //double*       ch4; // channel  4

    int           events_number;
    int           error_code;
    double        correction; // coefficient to correct ellipticity where y = y*correction

    //////////////////  Flags to choose data set and discard data points ////////////////
    int           datatype;

    bool          double_zero;
    bool          ToFselect;
    bool          ToFbandwidth; // if it is 1 you must set maxToF and minToF !!!!!!
    double        maxToF;
    double        minToF;
    bool          any_zero;
    bool          resorting;

};

#endif // COBOLDDATA_H
