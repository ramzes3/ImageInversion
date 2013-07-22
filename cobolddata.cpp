#include "cobolddata.h"
#include <math.h>

Cobolddata::Cobolddata(const char *LMF_Filename, int regime)
{
    // declaration default values
    data_array = NULL; array_size = 0; events_number = 0;
    array_x = NULL; array_y = NULL; array_tofx = NULL; array_tofy = NULL;
    error_code = 0;
    correction = 1;
    for(int i=0;i<NUM_IONS;i++)  hits_histo[i] = 0;


// check number of events in the lmf data file
    events_number = read_lmf_size(LMF_Filename);
    if(events_number == -1) {error_code = 1;}

    datatype = 1;

    double_zero = 1;
    ToFselect = 1;
    ToFbandwidth = 0;
    any_zero = 0;
    resorting = 0;

    //int error = read_lmf(LMF_Filename, regime);
}

Cobolddata::~Cobolddata()
{
    if(data_array){delete [] data_array;}
    if(array_x){delete [] array_x;}
    if(array_y){delete [] array_y;}
    if(array_tofx){delete [] array_tofx;}
    if(array_tofy){delete [] array_tofy;}
}

int Cobolddata::read_lmf_size(const char *LMF_Filename)
{
  LMF_IO * LMF = new LMF_IO(NUM_CHANNELS,NUM_IONS);

  if (!LMF->OpenInputLMF(LMF_Filename)) {
    delete LMF; return -1;
  }
  //cout << "Number of events = "<< LMF->uint64_Numberofevents <<endl;
  int n_events = LMF->uint64_Numberofevents;
  delete LMF;

  return n_events;
}

int Cobolddata::read_lmf(const char *LMF_Filename, int regime)
{
  unsigned int	event, counter, hits;
  unsigned __int32  number_of_hits[NUM_CHANNELS];
  memset(number_of_hits, 0, NUM_CHANNELS*4);
  __int32       iTDC[NUM_CHANNELS][NUM_IONS];
  double	dTDC[NUM_CHANNELS][NUM_IONS];


  LMF_IO * LMF = new LMF_IO(NUM_CHANNELS,NUM_IONS);

  if (!LMF->OpenInputLMF(LMF_Filename)) { delete LMF; return -1; }
  bool time_conversion_flag = false;

  if(array_x){delete [] array_x;}
  array_x = new double [events_number];
  if(array_y){delete [] array_y;}
  array_y = new double [events_number];
  //if(array_tofx){delete [] array_tofx;}
  //array_tofx = new double [events_number];
  //if(array_tofy){delete [] array_tofy;}
  //array_tofy = new double [events_number];

  //LMF->GetNumberOfChannels()
  //array_size = events_number * NUM_IONS * NUM_CHANNELS;
  //if(data_array){delete [] data_array;}
  //data_array = new double [array_size];
  //zeros(data_array, array_size);

  // Start reading event data:
  // ---------------------------------------
  double sumx, sumy, difx, dify, sumx_2, sumy_2;
  event = 0; counter = 0; stat_counts = 0;
  while(true) {
      if (LMF->ReadNextEvent())  {

          LMF->GetNumberOfHitsArray(number_of_hits);
          if (LMF->errorflag) {
            //LMF->GetErrorText(LMF->errorflag,error_text);
            error_code = LMF->errorflag;
            delete LMF; return false;
          }
          hits = number_of_hits[0];
          for(int j=1;j<4;j++) if(hits<number_of_hits[j]) hits = number_of_hits[j];
          hits_histo[hits]++;

          if (LMF->data_format_in_userheader==5) { LMF->GetTDCDataArray(&dTDC[0][0]);}
          else {LMF->GetTDCDataArray(&iTDC[0][0]);}

          ///////// START OF THE ROUTINE FOR POINTS SELSCTION
          int select = 0; // for(int k=0;k<4;k++){select[k]=0;}
          if(resorting && (hits>1)){  //f_multi_hit){
            double sumx_1, sumx_2, sumy_1, sumy_2, var;
            sumx_1 = iTDC[0][0] + iTDC[1][0];
            sumy_1 = iTDC[2][0] + iTDC[3][0];
            sumx_2 = iTDC[0][1] + iTDC[1][1];
            sumy_2 = iTDC[2][1] + iTDC[3][1];
            var = fabs(sumx_1 - sumy_1); select = 1;
            //if (fabs(sumx_1 - sumy_1) < 1000) {select = 1;}
            if (fabs(sumx_1 - sumy_2) < var) {var = fabs(sumx_1 - sumy_2); select = 2;}
            if (fabs(sumx_2 - sumy_1) < var) {var = fabs(sumx_2 - sumy_1); select = 3;}
            if (fabs(sumx_2 - sumy_2) < var) {var = fabs(sumx_2 - sumy_2); select = 4;}

            switch(select){
            case 1:
              sumx = iTDC[0][0] + iTDC[1][0];
              sumy = iTDC[2][0] + iTDC[3][0];
              difx = iTDC[0][0] - iTDC[1][0];
              dify = iTDC[2][0] - iTDC[3][0]; break;
            case 2:
              sumx = iTDC[0][0] + iTDC[1][0];
              sumy = iTDC[2][1] + iTDC[3][1];
              difx = iTDC[0][0] - iTDC[1][0];
              dify = iTDC[2][1] - iTDC[3][1]; break;
            case 3:
              sumx = iTDC[0][1] + iTDC[1][1];
              sumy = iTDC[2][0] + iTDC[3][0];
              difx = iTDC[0][1] - iTDC[1][1];
              dify = iTDC[2][0] - iTDC[3][0]; break;
            case 4:
              sumx = iTDC[0][1] + iTDC[1][1];
              sumy = iTDC[2][1] + iTDC[3][1];
              difx = iTDC[0][1] - iTDC[1][1];
              dify = iTDC[2][1] - iTDC[3][1]; break;

            }

          }
          else{
                sumx = iTDC[0][0] + iTDC[1][0];
                sumy = iTDC[2][0] + iTDC[3][0];
                difx = iTDC[0][0] - iTDC[1][0];
                dify = iTDC[2][0] - iTDC[3][0];
                sumx_2 = iTDC[0][1] + iTDC[1][1];
                sumy_2 = iTDC[2][1] + iTDC[3][1];
          }
/*
          if (LMF->data_format_in_userheader==5) {
              LMF->GetTDCDataArray(&dTDC[0][0]);
              sumx = dTDC[0][0] + dTDC[1][0];
              sumy = dTDC[2][0] + dTDC[3][0];
              difx = dTDC[0][0] - dTDC[1][0];
              dify = dTDC[2][0] - dTDC[3][0];
          }
          else {
              LMF->GetTDCDataArray(&iTDC[0][0]);
              sumx = iTDC[0][0] + iTDC[1][0];
              sumy = iTDC[2][0] + iTDC[3][0];
              difx = iTDC[0][0] - iTDC[1][0];
              dify = iTDC[2][0] - iTDC[3][0];
          }
*/
          ////////////////////////  END OF THE ROUTINE FOR POINTS SELSECTION /////////
          if (LMF->errorflag) {
            error_code = LMF->errorflag;
            delete LMF; return false;
          }

          bool flag = 1;

          if (difx == 0 && dify == 0) if(double_zero) flag = 0;

          if ((iTDC[0][0] == 0) || (iTDC[1][0] == 0) ||  (iTDC[2][0] == 0) ||  (iTDC[3][0] == 0)) if(any_zero) flag = 0;
          if (fabs(sumx - sumy) > 1000) if(ToFselect) flag = 0;

          if(ToFbandwidth) if(sumx > maxToF || sumy > maxToF || sumx < minToF || sumy < minToF) flag = 0;

          if(regime == 4){
              if(hits <2) flag =0;
              //if ((iTDC[0][1] == 0) || (iTDC[1][1] == 0) ||  (iTDC[2][1] == 0) ||  (iTDC[3][1] == 0)) flag = 0;
              //if ((iTDC[0][1] == 0) || (iTDC[1][1] == 0) ||  (iTDC[2][1] == 0) ||  (iTDC[3][1] == 0)) flag = 0;
          }

          double ToF2 = (sumx+sumy)/2;
          if( ( ToF2 < 3170 ) || (ToF2 > 3420) ) flag = 0;
          double rebin = 512/3000;
          if( sqrt(pow(difx*rebin + 4, 2) + pow(dify*rebin*correction + 2, 2)) > 230) flag = 0;

          if(flag)
          {
             switch (regime) {
                case 1:
                    array_x[counter] = difx;
                    array_y[counter] = dify*correction;
                    break;
                case 2:
                    array_x[counter] = sumx;
                    array_y[counter] = sumy;
                    break;
                case 3:
                    array_x[counter] = sqrt(difx*difx+dify*dify);
                    array_y[counter] = (sumx+sumy)/2;
                    break;
                case 4:
                    array_x[counter] = (sumx+sumy)/2;
                    array_y[counter] = (sumx_2+sumy_2)/2;
                    break;
                default :
                    array_x[counter] = difx;
                    array_y[counter] = dify;
                    break;
                }
             counter++; stat_counts++;
          }

          event++;
          if(event % (events_number/10) == 0) {
              int ix = (int) (event / (events_number/10)) - 1;
              points_var[ ix ] = stat_counts;
              stat_counts = 0;
          }
      } // end if(LMF->ReadNextEvent())
      if (LMF->errorflag) break;

  } // end while(true)
  array_size = counter;
  delete LMF;
  return 0;
}



void Cobolddata::zeros(double *array, int size)
{
    for(int i = 0; i<size;i++) array[i] = 0;
}
