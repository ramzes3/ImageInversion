#include "multifileinversion.h"
#include <QTime>
#include <QFileInfo>
#include <QDir>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

#include "mdata.h"
#include "pbinversion.h"

MultiFileInversion::MultiFileInversion(QObject *parent)
{
    files.clear();
    n_quadrant = 0;
}

void MultiFileInversion::run() { //Q_DECL_OVERRIDE {
    if (files.empty()) return;

    QString message;
    char filename[1024], fileout[1024], work_dir_c[1024];

    // Prepare output files
    QFileInfo info1(files.at(0));
    QByteArray var1 = info1.absolutePath().toLatin1();
    sprintf(work_dir_c,"%s", var1.data());
    sprintf(fileout,"%s/TotalCounts.%s", work_dir_c, "txt");
    ofstream outFile1(fileout, ios::out);
    if(!outFile1){ return; }

    sprintf(fileout,"%s/Summary_PES.%s", work_dir_c, "txt");
    ofstream outFile2(fileout, ios::out);
    if(!outFile2){ return; }

    sprintf(fileout,"%s/Summary_beta2.%s", work_dir_c, "txt");
    ofstream outFile3(fileout, ios::out);
    if(!outFile2){ return; }

    sprintf(fileout,"%s/Summary_beta4.%s", work_dir_c, "txt");
    ofstream outFile4(fileout, ios::out);
    if(!outFile2){ return; }
// //////////////////////////////////////////////////////

    /* expensive or blocking operation  */
    int error;
    int file_start, file_stop; file_start=0; file_stop = files.size(); //file_info.size();
    for(int k = file_start; k < file_stop; k++){ // main loop

        message = files.at(k);

        QByteArray var = files.at(k).toLatin1();// dir.toLatin1(); // convertion from QString to char*
        sprintf(filename,"%s",var.data()); // get path location
        switch (file_type) {
        case 1:
            error = mdata1->Load_2D(filename, mdata1->t_flag);
            break;
        case 2:
            error = mdata1->Load_binary(filename, mdata1->binary_depth, mdata1->size_x, mdata1->size_y, mdata1->t_flag);
            break;
        default:
            break;
        }
        emit resultReady(message, 100*(k-file_start + 1)/(file_stop-file_start + 1));

        if (!error){
            message = "Cannot load file";
            continue;
        }

        if(mdata1->bkg_flag && subtractbkg_flag) {
            mdata1->Subtract_bkg();
            mdata1->Copy_2D_ini_to_processed();
        }

        switch (n_quadrant)
        {
        case 1: mdata1->Get_quarter1(invers1->xc_v, invers1->yc_v);
            break;
        case 2: mdata1->Get_quarter2(invers1->xc_v, invers1->yc_v);
            break;
        case 3: mdata1->Get_quarter3(invers1->xc_v, invers1->yc_v);
            break;
        case 4: mdata1->Get_quarter4(invers1->xc_v, invers1->yc_v);
            break;
        default: mdata1->Symmetrise(invers1->xc_v, invers1->yc_v);
            break;
        }        

        outFile1 << filename << "\t Total counts: \t" << mdata1->TotalCounts(mdata1->array_2D_initial, invers1->xc_v, invers1->yc_v, 2*invers1->dr_v)<<endl;

        invertImages();

        outFile2 << filename << "\t";
        outFile3 << filename << "\t";
        outFile4 << filename << "\t";
        for (int j=0;j<mdata1->size_spectrum;j++) {
            outFile2 << "\t"<< mdata1->array_1D_spectrum[j];
            outFile3 << "\t" << mdata1->array_coeffs[0+j*mdata1->nL];
            outFile4 << "\t" << mdata1->array_coeffs[1+j*mdata1->nL];
        }
        outFile2 << endl; outFile3 << endl; outFile4 << endl;

        if(saving_flags[0]) mdata1->SaveProcMatrix(filename);
        if(saving_flags[1]) mdata1->SaveInvMatrix(filename);
        if(saving_flags[2]) mdata1->SavePES(filename);
        if(saving_flags[3]) mdata1->SaveAng(filename);

        emit resultReady(message, 100*(k-file_start + 1)/(file_stop-file_start + 1));
        //while (t1.elapsed() < 1350) {};

    }

    outFile1.close(); outFile2.close(); outFile3.close(); outFile4.close();

    message = "thread is finished";
    emit resultReady(message, 100);
}

void MultiFileInversion::initialize(QStringList filelist, MData *data, PBInversion *inv, int ftype, int nquad, bool *flags){

    for(int i=0;i<filelist.size();i++) files.append(filelist.at(i));
    mdata1 = data; invers1 = inv;
    n_quadrant = nquad;
    file_type = ftype;
    for(int i=0;i<4;i++) saving_flags[i] = flags[i];
    subtractbkg_flag = flags[4];
}

void MultiFileInversion::invertImages(){

    if(!mdata1) { return;}
    if(!invers1) { return;}

    mdata1->declare_matrices(); // prepare array in the pBasex format

    invers1->polar_b(mdata1->input_2d_array, mdata1->output_2d_array, invers1->xc_v, invers1->yc_v, invers1->dr_v, invers1->nl_v, invers1->odd_v);
    //delete invers;

    mdata1->copy_f_to_d(mdata1->output_2d_array);

    mdata1->size_spectrum = invers1->NR;
    mdata1->copy_pes(invers1->ang);


/*

    MData *mdata1; mdata1 = new MData();
    mdata1->Nskip_lines = Skip_lines;

    sprintf(filename,"%s/CS2images2_125000.dat.proc",work_dir_c);
    if (!mdata1->Load_2D(filename)) {message("cannot open bkg file for reading"); return;}
    SetBkg();

    char fileout[255];
    //sprintf(fileout,"%s/test.%s", work_dir_c, "txt");
    sprintf(fileout,"%s/test.%s", work_dir_c, "txt");
    ofstream outFile1(fileout, ios::out);
    if(!outFile1){ return; }

    sprintf(fileout,"%s/Summary_PES.%s", work_dir_c, "txt");
    ofstream outFile2(fileout, ios::out);
    if(!outFile2){ return; }

    sprintf(fileout,"%s/Summary_beta2.%s", work_dir_c, "txt");
    ofstream outFile3(fileout, ios::out);
    if(!outFile2){ return; }

    sprintf(fileout,"%s/Summary_beta4.%s", work_dir_c, "txt");
    ofstream outFile4(fileout, ios::out);
    if(!outFile2){ return; }

    QByteArray var = files.at(k).toLatin1(); //convertion from QString to char*
    sprintf(filename,"%s",var.data());

    if (!mdata1->Load_2D(filename)) {message("cannot open file for reading"); continue;}

    outFile1 << filename << "\t Total: \t"<< mdata1->TotalCounts(mdata1->array_2D_initial, xc, yc, dr) << "\t Bkg: \t" << mdata1->TotalCounts(mdata1->array_2D_bkg, xc, yc, dr)<<endl;

    mdata1->Subtract_bkg();
    mdata1->Copy_2D_ini_to_processed();

    switch (res)
    {
    case 1: mdata1->Get_quarter1(xc, yc);
        break;
    case 2: mdata1->Get_quarter2(xc, yc);
        break;
    case 3: mdata1->Get_quarter3(xc, yc);
        break;
    case 4: mdata1->Get_quarter4(xc, yc);
        break;
    default: mdata1->Symmetrise(xc, yc);
        break;
    }

    InvertImage1();

    outFile2 << filename << "\t";
    outFile3 << filename << "\t";
    outFile4 << filename << "\t";
    for (int j=0;j<mdata1->size_spectrum;j++) {
        outFile2 << "\t"<< mdata1->array_1D_spectrum[j];
        outFile3 << "\t" << mdata1->array_coeffs[0+j*mdata1->nL];
        outFile4 << "\t" << mdata1->array_coeffs[1+j*mdata1->nL];
    }
    outFile2 << endl; outFile3 << endl; outFile4 << endl;

    //if(!mdata1->SaveAng(filename)) {message("there is no angular distribution"); }

    */
}
