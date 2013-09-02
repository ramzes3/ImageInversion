#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include "math.h"

using namespace std;
//using namespace QtConcurrent;

#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
//#include <QThread>
//#include <qtconcurrentrun.h>

#include "mdata.h"
#include "mimageview.h"
#include "pbinversion.h"
#include "mpop.h"
//#include "cobolddata.h"
#include "mspectrum.h"
#include "multifileinversion.h"

#include <direct.h>

#include "options.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //Qt::WindowFlags flags = this->windowFlags();
    //this->setWindowFlags(flags | Qt::CustomizeWindowHint | Qt::WindowStaysOnTopHint);
    QString folder1 = QApplication::applicationDirPath(); // applicationFilePath();
    strncpy(current_folder, folder1.toLocal8Bit().constData(), 512);

    msignal_mapper =  new QSignalMapper(this); // required for exclusive choice of the inversion area ""circle" or "full"

    //this->setWindowFlags(Qt::Dialog);
    ui->mainToolBar->setVisible(0);
    //ui->mainToolBar->setShown(0);
    ui->progressBar->setHidden(1); // not in use right now
    ui->pushButton_do->setHidden(1);
    ui->pushButton_test->setHidden(1);
    ui->pushButton_convert->setHidden(1);
    ui->pushButton_centre->setHidden(1);
    this->setFixedSize(380, 317); // setting fixed size

// Setting up initial parameters ///
    file_status = false; // file is not loaded;
    sprintf(filename,"d:/"); // initialise string for storing folder pathway
    mdata1 = NULL; image_window = NULL; invers = NULL; pop = NULL; mspectr = NULL;
    Skip_lines = 0;   // initially don't skip anything
    area_flag = false;  // by default a whole area of the image is taken for the inversion
    transpose_flag = false; subtractbkg_flag = false;
    typeModel = NULL;
    lthickness = 1; color1 = QColor(Qt::green);
    /////////// !!!!!!!!!!!!!!!!!!!!!!! ////////////////////////
    // Initialise settings first //
    Skip_lines = 0; x_pos = 0; y_pos = 0; c_radius = 200;
    nL  = 2; n_r = 256; n_phi = 256; n_basisf = 128;

    bit_depth = 31; binary_width = 1024; binary_height = 1024; // seetings for loading binary file
    Lfile_type = 0; n_quadrant = 0;

    for(int i=0;i<5;i++) saving_flags[i]=0;
    saving_flags[5] = true; saving_flags[6] = true; saving_flags[7] = true; saving_flags[8] = true; // save all files by default
    LoadSettings(); // Load all the initial values from the file if file exists; !!!



// Connect signals to public slots
// Menu bar

    // Menu File
    connect( ui->actionQuit, SIGNAL( triggered() ), this, SLOT( close() ) );
    connect( ui->actionASCII, SIGNAL( triggered() ), this, SLOT( Load_Image_ASCII() ) );
    connect( ui->actionBinary, SIGNAL( triggered() ), this, SLOT( Load_Image_Binary() ) );

    //connect( ui->actionLMF, SIGNAL( triggered() ), this, SLOT( Load_Image_LMF() ) );
    connect( ui->actionAll, SIGNAL( triggered() ), this, SLOT( SaveAll() ) );
    connect( ui->actionSelected, SIGNAL( triggered() ), this, SLOT( SaveSelected() ) );
    connect( ui->actionImage, SIGNAL( triggered() ), this, SLOT( SaveImage() ) );


    // Menu Settings
    connect( ui->actionLoad_File, SIGNAL( triggered() ), this, SLOT( GetSkipLines() ) );
    // ...Selecting image area //
    connect(ui->actionFull, SIGNAL( triggered() ), msignal_mapper, SLOT(map()));
    msignal_mapper->setMapping(ui->actionFull, 0);
    connect(ui->actionCircle_only_2, SIGNAL( triggered() ), msignal_mapper, SLOT(map()));
    msignal_mapper->setMapping(ui->actionCircle_only_2, 1);
    connect(msignal_mapper, SIGNAL(mapped(int )), this, SLOT( CheckSelectedArea(int) ));
    connect( ui->actionSave_settings, SIGNAL( triggered() ), this, SLOT( SaveSettings() ) );
    connect( ui->actionLoad_Settings, SIGNAL( triggered() ), this, SLOT( LoadSettings() ) );
    connect( ui->actionSettings, SIGNAL( triggered() ), this, SLOT( OptionDialog() ) );

    // Menu Operations
    connect(ui->actionInvert_a_set, SIGNAL( triggered() ), this, SLOT( Analise() ));
    connect(ui->actionInvert_files, SIGNAL( triggered() ), this, SLOT( AnaliseSelected() ));
    connect( ui->actionGenerate_Basis_2, SIGNAL( triggered() ), this, SLOT( GenerateBasis() ) );
    connect( ui->actionSet_bkg, SIGNAL( triggered() ), this, SLOT( SetBkg() ) );
    connect( ui->actionShow_bkg, SIGNAL( triggered() ), this, SLOT( ShowBkg() ) );
    connect( ui->actionShow_image, SIGNAL( triggered() ), this, SLOT( ReloadImage() ) );
    connect( ui->actionSubtract_bkg, SIGNAL( triggered() ), this, SLOT( Subtract_bkg() ) );
    connect( ui->actionAdd_bkg, SIGNAL( triggered() ), this, SLOT( Add_bkg() ) );
    connect( ui->actionImgAdd, SIGNAL( triggered() ), this, SLOT( Add_Image_ASCII() ) );

    // Menu Help
    connect( ui->actionAbout, SIGNAL( triggered() ), this, SLOT( About() ) );


// Buttons
    connect( ui->pushButton_invert, SIGNAL( clicked() ), this, SLOT( InvertImage() ) );
    connect( ui->pushButton_Reload, SIGNAL( clicked() ), this, SLOT( ReloadImage() ) );
    connect( ui->pushButton_sym, SIGNAL( clicked() ), this, SLOT( Symmetrise()) );
    connect( ui->pushButton_polar, SIGNAL( clicked() ), this, SLOT( GetPolarImage()) );
    connect( ui->pushButton_getquarter, SIGNAL( clicked() ), this, SLOT( GetQuarter()) );

    connect( ui->pushButton_test, SIGNAL( clicked() ), this, SLOT( Test() ) );
    connect( ui->pushButton_do, SIGNAL( clicked() ), this, SLOT( DataManipulation()) );
    connect( ui->pushButton_convert, SIGNAL( clicked() ), this, SLOT( ConvertImages()) );
// Timers
    //timer_loading = new QTimer(this);
    //connect(timer_loading, SIGNAL(timeout()), this, SLOT(MUpdate()));

    DeclareMImageView();
}

MainWindow::~MainWindow()
{
    SaveSettings();

    delete ui;
    if(mdata1) delete mdata1;
    if(image_window) delete image_window;
    if(invers) {if (invers->isRunning()) invers->terminate(); delete invers;}
    if(pop) delete pop;
    if(mspectr) delete mspectr;
    if(typeModel) delete typeModel;
}

void MainWindow::Load_Image_ASCII()
{
    //if(mspectr) {delete mspectr; mspectr = NULL;}
    char buffer[255];
    getPath();
    if(!file_status){return;}
    //if(mdata1) delete mdata1;
    if(!mdata1) mdata1 = new MData();
    mdata1->Nskip_lines = Skip_lines;   

    if (!mdata1->Load_2D(filename, transpose_flag)) {message("cannot open file for reading"); return;}
    //if (!mdata1->Copy_2D_ini_to_processed()){message("cannot creat processed array"); return;}
    //mdata1->normalise_2darray(mdata1->array_2D_initial);
    double max_content = mdata1->find_max(mdata1->array_2D_initial);
    double counts = mdata1->TotalCounts(mdata1->array_2D_initial, mdata1->GetsizeX()/2, mdata1->GetsizeY()/2, mdata1->GetsizeY());
    sprintf (buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message( buffer );

    image_window->DrawArray( mdata1->array_2D_initial, mdata1->size_x, mdata1->size_y, max_content );
    image_window->update();   
    image_window->show();
}

void MainWindow::Load_Image_Binary(){

    char buffer[255];
    getPath();
    if(!file_status){return;}
    //if(mdata1) delete mdata1;
    if(!mdata1) mdata1 = new MData();
    mdata1->Nskip_lines = Skip_lines;

    if (!mdata1->Load_binary(filename, bit_depth, binary_width, binary_height, transpose_flag)) {message("cannot open file for reading"); return;}
    //if (!mdata1->Copy_2D_ini_to_processed()){message("cannot creat processed array"); return;}
    //mdata1->normalise_2darray(mdata1->array_2D_initial);
    double max_content = mdata1->find_max(mdata1->array_2D_initial);
    double counts = mdata1->TotalCounts(mdata1->array_2D_initial, mdata1->GetsizeX()/2, mdata1->GetsizeY()/2, mdata1->GetsizeY());
    sprintf (buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message( buffer );

    image_window->DrawArray( mdata1->array_2D_initial, mdata1->size_x, mdata1->size_y, max_content );
    image_window->update();
    image_window->show();
}

void MainWindow::Add_Image_ASCII()
{
    //if(mspectr) {delete mspectr; mspectr = NULL;}
    char buffer[255];
    getPath();
    if(!file_status){return;}
    //if(mdata1) delete mdata1;
    if(!mdata1) mdata1 = new MData();
    mdata1->Nskip_lines = Skip_lines;

    if (!mdata1->Copy_2D_ini_to_result()){message("cannot creat processed array"); return;}
    if (!mdata1->Load_2D(filename,transpose_flag)) {message("cannot open file for reading"); return;}
    if (!mdata1->Add_result_ini_to_ini()){message("cannot creat processed array"); return;}
    //mdata1->normalise_2darray(mdata1->array_2D_initial);
    double max_content = mdata1->find_max(mdata1->array_2D_initial);
    double counts = mdata1->TotalCounts(mdata1->array_2D_initial, mdata1->GetsizeX()/2, mdata1->GetsizeY()/2, mdata1->GetsizeY());
    sprintf (buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message( buffer );

    image_window->DrawArray( mdata1->array_2D_initial, mdata1->size_x, mdata1->size_y, max_content );
    image_window->update();
    image_window->show();
}

void MainWindow::InvertImage()
{
    if(!mdata1) {message("There is no data loaded"); return;}
    if(invers) delete invers;
    invers = new PBInversion();
    connect( invers, SIGNAL( sent_respond(const char*) ), this, SLOT( message(const char*)) );
    connect( invers, SIGNAL( sent_error(const char*) ), this, SLOT( merror(const char*)) );

    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();
    int dr = ui->spinBox_dr->value();
    // prepare image data
    if( area_flag == 1) mdata1->circle_mask( xc, yc, dr );
    mdata1->declare_matrices(); // prepare array in the pBasex format

    int nl = ui->spinBox_nL->value();
    int odd; if(ui->checkBox_odd->isChecked()) odd = true; else odd = false;

    //invers->polar_a(mdata1->input_2d_array, mdata1->output_2d_array, xc, yc, dr, nl, odd);
    SetBasisDim(); // change the dimesions of the basis set in PBasex in the invers class
    invers->set_size( mdata1->size_x, mdata1->size_y );
    //invers->start(QThread::NormalPriority);
    invers->load_matrices(nl,odd);    
    if(!invers->m_isloaded){message("can not load basex matricies"); return;}
    invers->polar_b(mdata1->input_2d_array, mdata1->output_2d_array, xc, yc, dr, nl, odd);
    //delete invers;
    mdata1->copy_f_to_d(mdata1->output_2d_array);
    if(odd) mdata1->nL = nl;
    else mdata1->nL = nl/2;
    mdata1->size_spectrum = invers->NR;
    mdata1->copy_pes(invers->ang);

    char buffer[100];
    double counts = (int) mdata1->TotalCounts(mdata1->array_2D_initial, xc, yc, dr);
    sprintf(buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message(buffer);

    DeclareMImageView();
    //mdata1->normalise_2darray(mdata1->array_2D_result);
    double max_content = mdata1->find_max(mdata1->array_2D_result);
    image_window->DrawArray(mdata1->array_2D_result, mdata1->size_x, mdata1->size_y, max_content);
    image_window->update();
    image_window->show();

    ShowSpectrum();
    //if(invers) {delete invers; invers = NULL;}
    disconnect( invers, SIGNAL( sent_respond(const char*) ), this, SLOT( message(const char*)) );
    disconnect( invers, SIGNAL( sent_error(const char*) ), this, SLOT( merror(const char*)) );
}

void MainWindow::Declear_pBasex(){

    if(!mdata1) {message("There is no data loaded"); return;}
    if(invers) delete invers;
    invers = new PBInversion();
    connect( invers, SIGNAL( sent_respond(const char*) ), this, SLOT( message(const char*)) );
    connect( invers, SIGNAL( sent_error(const char*) ), this, SLOT( merror(const char*)) );

    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();
    int dr = ui->spinBox_dr->value();
    // prepare image data
    if( area_flag == 1) mdata1->circle_mask( xc, yc, dr );

    int nl = ui->spinBox_nL->value();
    int odd; if(ui->checkBox_odd->isChecked()) odd = true; else odd = false;

    SetBasisDim(); // change the dimesions of the basis set in PBasex in the invers class
    invers->set_size( mdata1->size_x, mdata1->size_y );
    invers->xc_v = xc; invers->yc_v = yc; invers->dr_v = dr;
    invers->nl_v = nl; invers->odd_v = odd;
    //invers->start(QThread::NormalPriority);
    invers->load_matrices(nl,odd);
    if(!invers->m_isloaded){message("can not load basex matricies"); return;}
}

void MainWindow::InvertImage1(){
    if(!mdata1) {message("There is no data loaded"); return;}
    if(!invers) {message("pBasex class has not been decleared"); return;}

    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();
    int dr = ui->spinBox_dr->value();

    int nl = ui->spinBox_nL->value();
    int odd; if(ui->checkBox_odd->isChecked()) odd = true; else odd = false;

    mdata1->declare_matrices(); // prepare array in the pBasex format

    invers->polar_b(mdata1->input_2d_array, mdata1->output_2d_array, xc, yc, dr, nl, odd);
    //delete invers;

    mdata1->copy_f_to_d(mdata1->output_2d_array);
    if(odd) mdata1->nL = nl;
    else mdata1->nL = nl/2;
    mdata1->size_spectrum = invers->NR;
    mdata1->copy_pes(invers->ang);

    char buffer[100];
    double counts = (int) mdata1->TotalCounts(mdata1->array_2D_initial, xc, yc, dr);
    sprintf(buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message(buffer);

}

void MainWindow::ShowSpectrum()
{
    if(!mdata1) {message("There is no data loaded"); return;}
    if(!mdata1->array_1D_spectrum) {message("Image was not inverted"); return;}
    if(mspectr) delete mspectr;
    mspectr = new MSpectrum(this);
    mspectr->setWindowFlags(mspectr->windowFlags() | Qt::Window);
    mspectr->n_points = n_r;
    mspectr->SetData(mdata1->array_1D_spectrum);
    mspectr->show();
}

void MainWindow::GenerateBasis()
{
    if(invers) delete invers;
    invers = new PBInversion();
    connect( invers, SIGNAL( sent_respond(const char*) ), this, SLOT( message(const char*)) );
    connect( invers, SIGNAL( sent_error(const char*) ), this, SLOT( merror(const char*)) );

    int nl = ui->spinBox_nL->value();
    int odd; if(ui->checkBox_odd->isChecked()) odd = true; else odd = false;

    SetBasisDim(); // change the dimesions of the basis set
    //invers->write_forward(nl, odd);
    invers->polar_a(nl, odd);           // run basis generator in the thread
    invers->start(QThread::NormalPriority);   //  run basis generator in the thread
}

void MainWindow::Transpose()
{
    if(!mdata1) {message("There is no data loaded"); return;}
    if(!mdata1->array_2D_processed) {message("There is no data loaded"); return;}
    mdata1->transpose();

    DeclareMImageView();
    //mdata1->normalise_2darray(mdata1->array_2D_processed);
    double max_content = mdata1->find_max(mdata1->array_2D_processed);
    image_window->DrawArray(mdata1->array_2D_processed, mdata1->size_x, mdata1->size_y, max_content);
    image_window->update();
    image_window->show();
}

void MainWindow::AnaliseSelected()
{

    if(!mdata1) {message("There is no data loaded"); return;}
    //if(!invers) {message("pBasex class has not been decleared"); return;}
    QStringList files = QFileDialog::getOpenFileNames(
                            this,
                            "Select one or more files to join",
                            filename,
                            "Images (*.*)");
//  Preparations ///

    int nl = ui->spinBox_nL->value();
    int odd; if(ui->checkBox_odd->isChecked()) odd = true; else odd = false;
    if(odd) mdata1->nL = nl;
    else mdata1->nL = nl/2;

    Declear_pBasex();
    mdata1->binary_depth = bit_depth;
    mdata1->t_flag = transpose_flag;

// //////////////////////////////////

    ui->progressBar->setValue(0);
    ui->progressBar->setVisible(1);

    MultiFileInversion *task1 = new MultiFileInversion(this);
    //connect(task1, &MultiFileInversion::finished, this, &MainWindow::wakeup);
    //connect(task1, &MultiFileInversion::resultReady, this, &MainWindow::invProgress);
    //connect(task1, SIGNAL(finished(QPrivateSignal)), this, SLOT(wakeup()));
    connect(task1, SIGNAL(resultReady(QString, int)), this, SLOT(invProgress(QString, int)));
    saving_flags[4] = subtractbkg_flag;
    task1->initialize(files, mdata1, invers, Lfile_type, n_quadrant, saving_flags);
    task1->start(QThread::NormalPriority);
}

void MainWindow::wakeup(){
    ui->progressBar->setVisible(0); message("thread is fnished");
}

void MainWindow::invProgress(QString msend, int iter){
    QByteArray var = msend.toLatin1();
    message(var.data());

    if (iter == 100) ui->progressBar->setVisible(0);
    else ui->progressBar->setValue(iter);
}

void MainWindow::Analise()
{
    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();
    //int dr = ui->spinBox_dr->value();
//// Get working directory  ////////
    QString work_dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                    filename,
                                                    QFileDialog::ShowDirsOnly
                                                    | QFileDialog::DontResolveSymlinks);

    if(work_dir==NULL) {message("operation cancelled"); return;}
    char file_template[255], bufferc[255];
    QDir workDir(work_dir);
    QStringList dir_list = workDir.entryList(QDir::AllDirs);
    QByteArray buffer = work_dir.toLocal8Bit();
    const char *work_dir_c = buffer.data();

 //// Get file name template ////
    bool ok;
    QString file_temp;
    file_temp = QInputDialog::getText(this, tr("QInputDialog::getText()"),
                                         tr("File name template"), QLineEdit::Normal,
                                         tr("PumpProbe.txt"), &ok);
    if (ok && !file_temp.isEmpty()) {message("ok");}
    else {message("not ok"); return;}
    sprintf(file_template,"%s", file_temp.toLatin1().data());

    QApplication::setOverrideCursor(Qt::WaitCursor);

    if(mdata1) delete mdata1;
    mdata1 = new MData();
    mdata1->Nskip_lines = Skip_lines;

    char dirname[255];//, buffer[255];
    ui->progressBar->setValue(0);
    ui->progressBar->setVisible(1);
    int file_start, file_stop; file_start=0; file_stop = dir_list.size();
    for(int k = file_start; k < file_stop; k++){ // main loop

        ui->progressBar->setValue(100*(k-file_start + 1)/(file_stop-file_start + 1));

        QByteArray var = dir_list.at(k).toLatin1();// dir.toLatin1(); //  convertion from QString to char*
        sprintf(dirname,"%s",var.data()); // get path location

        sprintf(filename,"%s/%s/%s", work_dir_c, dirname, file_template);

        if (!mdata1->Load_2D(filename,transpose_flag)) {message("cannot open file for reading"); continue;}
        //mdata1->normalise_2darray(mdata1->array_2D_initial);
        // select quarter
        if(mdata1->bkg_flag && subtractbkg_flag) mdata1->Subtract_bkg();

        switch (n_quadrant)
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

        double max_content = mdata1->find_max(mdata1->array_2D_initial);

        image_window->DrawArray( mdata1->array_2D_initial, mdata1->size_x, mdata1->size_y, max_content);
        image_window->update();
        image_window->show();

        // Ciaran modifications
        sprintf(filename,"%s/%s/inversion", work_dir_c, dirname, file_template); // inversion folder path string
        mkdir(filename);
        sprintf(filename,"%s/%s/inversion/x%iy%i", work_dir_c, dirname, xc, yc); //output folder name
        mkdir(filename);

        switch (n_quadrant){
            case 0:
                sprintf(filename,"%s/%s/inversion/x%iy%i/whole", work_dir_c, dirname, xc, yc); //output folder
                break;
            default:
                sprintf(filename,"%s/%s/inversion/x%iy%i/q%i", work_dir_c, dirname, xc, yc, n_quadrant); //output folder
                break;
        }
        //

        sprintf(bufferc, "%s_orig", filename);
        SaveImage(bufferc);


        InvertImage();
        qApp->processEvents(QEventLoop::ExcludeUserInputEvents);

        if(saving_flags[2])
            if(!mdata1->SavePES(filename)) {message("there is no spectrum"); }
        if(saving_flags[3])
            if(!mdata1->SaveAng(filename)) {message("there is no angular distribution"); }
        if(saving_flags[1])
            if(!mdata1->SaveInvMatrix(filename)) {message("there is no inverted image"); }
        if(saving_flags[0])
            if(!mdata1->SaveProcMatrix(filename)) {message("there is no image image"); }

    } // main loop
    ui->progressBar->setVisible(0);
    if(mdata1) delete mdata1;
    mdata1 = NULL;

    QApplication::restoreOverrideCursor();
    return;
    //strcpy (current_folder,dirname);
}

void MainWindow::getPath()
{
    //char buffer[255]; // buffer for messages
    QStringList filters;
    filters << "Any files (*)"
            << "Text files (*.txt)"
            << "Data (*.dat)";

    QFileDialog newdialog(this);
    newdialog.setNameFilters(filters);
    newdialog.setReadOnly(true);

    //char* inidir;
    //inidir = "C:/Documents and Settings/Laser1/Desktop/Roman/Data/t_0_49_156mm6000counts";

    //QFileDialog* newdialog;
    //newdialog = new QFileDialog();
    QString path = newdialog.getOpenFileName(this, tr("Open File"), filename);
    //QString path = QFileDialog::getOpenFileName(this, tr("Open File"), inidir);
    //QString path = QString("C:/Documents and Settings/Laser1/Desktop/Roman/Data/t_0_49_156mm6000counts/PumpProbe.asc");
    if (!path.isEmpty()){ file_status = 1;} // file name has been choosen
    else{message("file name has not been choosen"); file_status = 0; return;}
    QByteArray var = path.toLatin1(); //  convertion from QString to char*
    sprintf(filename,"%s",var.data()); // get path location

    //delete newdialog;

}

void MainWindow::setCircle(int posx, int posy, int radius)
{
    ui->spinBox_xc->setValue( posx ); x_pos = posx;
    ui->spinBox_yc->setValue( posy ); y_pos = posy;
    ui->spinBox_dr->setValue(radius); c_radius = radius;
}

void MainWindow::ModifyPixel(int posx, int posy)
{
    if(!mdata1) return;
    if(posx < mdata1->size_x && posx>0 && posy < mdata1->size_y && posy>0) ;
    else return;
    bool ok; double res; char buffer[255];
    sprintf(buffer,"Do you want to change pixel (%d, %d): %f with the number below?",posx, posy, mdata1->Get_pix_value(posx, posy));
    res = QInputDialog::getDouble(this, tr("getDouble()"),
                                           tr(buffer), mdata1->Get_pix_value(posx+1, posy+1), 0, 100000000, 2, &ok);
    if (ok) {
        res = mdata1->Change_pix_value(posx, posy, res);
        ReloadImage();
        sprintf(buffer,"pixel (%d, %d) has been modified to: %f",posx, posy, res);
        message(buffer);
    }
    else message("change cancelled");
}

void MainWindow::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    /*
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
*/

    if(e->type()==QEvent::WindowStateChange&&isMinimized()){
        if(image_window) image_window->setWindowState(Qt::WindowMinimized);
        if(mspectr) mspectr->setWindowState(Qt::WindowMinimized);
    }
    else if(e->type()==QEvent::WindowStateChange&&isActiveWindow()){
        if(image_window) image_window->setWindowState(Qt::WindowNoState);
        if(mspectr) mspectr->setWindowState(Qt::WindowNoState);
    }

}

void MainWindow::SaveAll(){

    QFileDialog dialog(this);
    QString path = dialog.getSaveFileName(this, tr("Save File"), filename);

    if (!path.isEmpty()){} // file name has been choosen
    else{message("file name has not been choosen"); return;}
    QByteArray var = path.toLatin1(); //  convertion from QString to char*
    sprintf(filesave,"%s",var.data()); // get path location

    if(!mdata1) {message("there is no data"); return;}
    if(!mdata1->SavePES(filesave)) {message("there is no spectrum"); }
    if(!mdata1->SaveAng(filesave)) {message("there is no angular distribution"); }
    if(!mdata1->SaveInvMatrix(filesave)) {message("there is no inverted image"); }
    if(!mdata1->SaveProcMatrix(filesave)) {message("there is no image image"); }
    message("files successfully saved");
}


void MainWindow::SaveSelected(){

    if(!mdata1) {message("there is no data"); return;}

    QFileDialog dialog(this);
    QString path = dialog.getSaveFileName(this, tr("Save File"), filename);

    if (!path.isEmpty()){} // file name has been choosen
    else{message("file name has not been choosen"); return;}
    QByteArray var = path.toLatin1(); //  convertion from QString to char*
    sprintf(filesave,"%s",var.data()); // get path location

    /*
    //fp >> savePES >> saveANG >> savePROC >> saveINV;
    fp >> saving_flags[5] >> saving_flags[6] >> saving_flags[7] >> saving_flags[8];
    */

    if (saving_flags[2]){
        if(!mdata1->SavePES(filesave)) {message("there is no spectrum"); }
    }
    if (saving_flags[3]){
        if(!mdata1->SaveAng(filesave)) {message("there is no angular distribution"); }
    }
    if (saving_flags[0]){
        if(!mdata1->SaveProcMatrix(filesave)) {message("there is no image image"); }
    }
    if (saving_flags[1]){
        if(!mdata1->SaveInvMatrix(filesave)) {message("there is no inverted image"); }
    }
    message("files successfully saved");

}

void MainWindow::SaveImage()
{
    QStringList filters;
    filters << "Image (*.jpg)"
            << "Any files (*)";

    QFileDialog dialog(this);
    dialog.setNameFilters(filters);
    QString path = dialog.getSaveFileName(this, tr("Save File"), filename);

    if (!path.isEmpty()){} // file name has been choosen
    else{message("file name has not been choosen"); return;}
    QByteArray var = path.toLatin1(); //  convertion from QString to char*
    sprintf(filesave,"%s",var.data()); // get path location

    SaveImage(filesave);
}

void MainWindow::SaveImage(char *subname)
{
    if(!image_window) { message("nothing to save"); return; }
    if(!image_window->m_image) { message("nothing to save"); return; }

    char buffer[255];
    sprintf(buffer, "%s.jpg",subname);
    QString image_name(buffer);
    if(image_window->m_image->save(image_name, "JPEG", 100))
        message("image successfully saved");
    else message("can not save the image");
}

void MainWindow::ReloadImage()
{
    if(!mdata1) {message("There is no data to reload"); return;}
    mdata1->reload_data();
    char buffer[255];
    sprintf (buffer, "Array size: %d, %d", mdata1->GetsizeX(), mdata1->GetsizeY() );
    message( buffer );

    if(!mdata1->array_2D_processed) {message("There is no data to reload"); return;}
    double max_content = mdata1->find_max(mdata1->array_2D_processed);

    DeclareMImageView();

    image_window->DrawArray( mdata1->array_2D_processed, mdata1->size_x, mdata1->size_y, max_content );
    image_window->update();
    image_window->show();

    double counts = mdata1->TotalCounts(mdata1->array_2D_initial, mdata1->GetsizeX()/2, mdata1->GetsizeY()/2, mdata1->GetsizeY());
    sprintf (buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message( buffer );
}

void MainWindow::Symmetrise()
{
    // symmetrises image
    if(!mdata1) return;
    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();

    mdata1->Symmetrise(xc, yc);

    double max_content = mdata1->find_max(mdata1->array_2D_processed);

    image_window->DrawArray( mdata1->array_2D_processed, mdata1->size_x, mdata1->size_y, max_content);
    image_window->update();
    image_window->show();
    ///////// end of symmertrise the image
}

void MainWindow::GetQuarter(){

    if(!mdata1) {message("please open a file first"); return;}
    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();

    bool ok; int res; char buffer[255];
    sprintf(buffer,"Please choose your favourite quarter");
    res = QInputDialog::getInt(this, tr("GetQuarter"), tr(buffer), 1, 1, 4, 1, &ok);
    if (ok) {
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
        default: message("wrong choice, please try again"); return;
            break;
        }
    }
    else {message("action cancelled"); return;}

    double max_content = mdata1->find_max(mdata1->array_2D_processed);

    image_window->DrawArray( mdata1->array_2D_processed, mdata1->size_x, mdata1->size_y, max_content);
    image_window->update();
    image_window->show();

    message("How is it?");
}

void MainWindow::GetPolarImage()
{
    if(!mdata1) return;
    int xc = ui->spinBox_xc->value();
    int yc = ui->spinBox_yc->value();

    mdata1->CartToPolar(xc, yc, n_r, n_phi);

    double max_content = mdata1->find_max(mdata1->polar_2D_image, mdata1->size_phi * mdata1->size_r);

    image_window->DrawArray( mdata1->polar_2D_image, mdata1->size_r, mdata1->size_phi, max_content);
    image_window->update();
    image_window->show();

}

void MainWindow::Subtract_bkg()
{
    if(!mdata1) {message("There is no data to reload"); return;}
    mdata1->Subtract_bkg();
    mdata1->Copy_2D_ini_to_processed();
    ReloadImage();
}

void MainWindow::Add_bkg()
{
    if(!mdata1) {message("There is no data to reload"); return;}
    mdata1->Add_bkg(1);
    ShowBkg();
}

void MainWindow::GetSkipLines()
{
    bool ok;
    Skip_lines = QInputDialog::getInt(this, tr("QInputDialog::getInteger()"), tr("Number of lines"), Skip_lines, 0, 1000, 1, &ok);
}

void MainWindow::CheckSelectedArea(int moption)
{
    if(moption == 1) {
        area_flag = 1;
        if(ui->actionCircle_only_2->isChecked()) ui->actionFull->setChecked(0);
        else ui->actionCircle_only_2->setChecked(1);
        message("Circle area");
        return;
    }
    if(moption == 0) {
        area_flag = 0;
        if(ui->actionFull->isChecked()) ui->actionCircle_only_2->setChecked(0);
        else ui->actionFull->setChecked(1);
        message("Whole area");
        return;
    }
}

void MainWindow::SetBkg()
{
    if(!mdata1) {message("there is no data"); return;}
    mdata1->SetBackground();
    message("Background image is set");
}

void MainWindow::ShowBkg()
{
    if(!mdata1) {message("There is no data to reload"); return;}
    if(!mdata1->array_2D_bkg) {message("Please set bkg first"); return;}
    char buffer[255];
    sprintf (buffer, "Array size: %d, %d", mdata1->GetsizeX(), mdata1->GetsizeY() );
    message( buffer );

    double max_content = mdata1->find_max(mdata1->array_2D_bkg);

    image_window->DrawArray( mdata1->array_2D_bkg, mdata1->size_x, mdata1->size_y, max_content );
    image_window->update();
    image_window->show();

    double counts = mdata1->TotalCounts(mdata1->array_2D_bkg, mdata1->GetsizeX()/2, mdata1->GetsizeY()/2, mdata1->GetsizeY());
    sprintf (buffer, "Array size: %d, %d; Counts: %4.0f", mdata1->GetsizeX(), mdata1->GetsizeY(), counts );
    message( buffer );
}

void MainWindow::MUpdate()
{
    //char buffer[255], var[3];
    //var[0] = '.';
    //sprintf(buffer, "Loading...%s", var[0]);
    //message(buffer);
}

void MainWindow::SetBasisDim()
{
    if(invers){
        invers->NK = n_basisf;
        invers->NR = n_r;
        invers->NTH = n_phi;
    }
    else {message("There is no PBasex class localised"); return;}
}

void MainWindow::message(const char* message)
{
   ui->message_box->setText(message);
}

void MainWindow::merror(const char* message)
{
   error_list << tr(message);
   if(error_list.size() > 10) error_list.removeFirst();
   if (typeModel) delete typeModel;
   typeModel = new QStringListModel(error_list, this);
   ui->error_box->setModel(typeModel);
}

void MainWindow::About()
{
    QMessageBox mbox1;
    mbox1.setWindowTitle("About");
    mbox1.setStandardButtons(QMessageBox::Close);
    mbox1.setText("VMI Image inversion.");
    mbox1.setInformativeText("Program uses Qt GUI interface and compiled using Qt 5.0.2. \n \n"
                             "The program is provided for free with no warranty of any kind. "
                             "However, if you have any questions or suggestions please contact: Roman Spesyvtsev. \n \n"
                             "Current version of the program uses pBasex routine for image inversion, developed by G. A. Garcia, L. Nahon, and I. Powis");
    mbox1.exec();
}

void MainWindow::ConvertQstringToChar(QString qstring1, char* char_out){
    QByteArray var = qstring1.toLocal8Bit();
    strcpy(char_out, var.data());
    return;
}

void MainWindow::ConvertImages()
{    
    if(!mdata1) {message("There is no data loaded"); return;}
    double db_var, *ptr_var;

    getPath();
    if(!file_status){return;}
    if (!mdata1->Load_1D(filename)) {message("cannot open file for reading"); return;}

    QApplication::setOverrideCursor(Qt::WaitCursor);

    char work_dir_c[255],fileout[255], work_dir_c2[255], work_dir_c3[255], buffer[255];//, buffer[255];

    sprintf(work_dir_c,"D:/Roman/Data/20130829/VMIimages");
    //sprintf(work_dir_c2,"D:/Roman/Data/20130828/SecondRun");
    //sprintf(work_dir_c3,"D:/Roman/Data/20130828/ThirdRun");
/*
    QStringList files = QFileDialog::getOpenFileNames(
                            this,
                            "Select one or more files to join",
                            "D:/Musor/CS2Data",
                            "Images (*.*)");
*/
    QStringList files, files2, files3;

    sprintf(fileout,"%s/CS21000nmstep_-1500000.dat", work_dir_c);
    files << fileout;

    //sprintf(fileout,"%s/CS2_VMI_2_-1500000.dat", work_dir_c2);
    //files2 << fileout;

    //sprintf(fileout,"%s/CS2_VMI_3_-1500000.dat", work_dir_c3);
    //files3 << fileout;

    int delay = 295000;
    for(int i = 0; i<305;i++) {
        sprintf(fileout,"%s/CS21000nmstep_%d.dat", work_dir_c, delay);
        files << fileout;
        //sprintf(fileout,"%s/CS2_VMI_2_%d.dat", work_dir_c2, delay);
        //files2 << fileout;
        //sprintf(fileout,"%s/CS2_VMI_3_%d.dat", work_dir_c3, delay);
        //files3 << fileout;

        delay = delay + 1000;
    }



    if (!files.isEmpty()){} // file name has been choosen
    else{message("files have not been choosen"); return;}

    message("loading...");

    mdata1->Nskip_lines = 0;
    //sprintf(work_dir, "E:/Roman/Data/AnilineCopy/20120704");

    double max_content, total_counts;//, check_number;
    //var_array = new double [7*mdata1->size_y]; // load 5 images and calculate mean and stdev


    ui->progressBar->setValue(0);
    ui->progressBar->setVisible(1);

////////////////////////// A summary file  ////////////////////////////
/*
    sprintf(fileout,"%s/Summary_ALL_empty.%s", work_dir_c2, "txt");
    ofstream outFile1(fileout, ios::out);
    if(!outFile1){ message("cannot create summary file"); return; }
*/
///////////////////////////////

    //int xc = ui->spinBox_xc->value();
    //int yc = ui->spinBox_yc->value();
    //int dr = ui->spinBox_dr->value();

    int file_start, file_stop; file_start=0; file_stop = files.size();
    for(int k = file_start; k < file_stop; k++){ // main loop

        ui->progressBar->setValue(100*(k-file_start + 1)/(file_stop-file_start + 1));

        QByteArray var = files.at(k).toLatin1();// dir.toLatin1(); // convertion from QString to char*
        sprintf(filename,"%s",var.data()); // get path location
        //if (!mdata1->Load_2D(filename, transpose_flag)) {message("cannot open file for reading"); continue;}
        if (!mdata1->Load_binary(filename, 32, 1024, 1024, transpose_flag)) {message("cannot open file for reading 1"); continue;}
/*
        //  Adding another file
        if (!mdata1->Copy_2D_ini_to_result()){message("cannot creat result array"); continue;}
        var = files2.at(k).toLatin1(); // dir.toLatin1(); // convertion from QString to char*
        sprintf(filename,"%s",var.data()); // get path location
        //if (!mdata1->Load_2D(filename,transpose_flag)) {message("cannot open file for reading"); continue;}
        if (!mdata1->Load_binary(filename, 32, 1024, 1024, transpose_flag)) {message("cannot open file for reading 2"); continue;}
        if (!mdata1->Add_result_ini_to_ini()){message("cannot creat processed array"); continue;}
        // end of the adding file

        //  Adding another file        
        var = files3.at(k).toLatin1(); // dir.toLatin1(); // convertion from QString to char*
        sprintf(filename,"%s",var.data()); // get path location
        //if (!mdata1->Load_2D(filename,transpose_flag)) {message("cannot open file for reading"); continue;}
        if (!mdata1->Load_binary(filename, 32, 1024, 1024, transpose_flag)) {message("cannot open file for reading 2"); continue;}
        if (!mdata1->Add_result_ini_to_ini()){message("cannot creat processed array"); continue;}
        // end of the adding file
*/
        // Normalise array;
        mdata1->normalise_2darray(mdata1->array_2D_processed, mdata1->time_profile[k]);

        sprintf(fileout,"%s_sum", filename);
        if(!mdata1->SaveProcMatrix(filename)) {message("there is no image image"); }

        //total_counts = mdata1->TotalCounts(mdata1->array_2D_initial, xc, yc, dr); //mdata1->find_max(mdata1->array_2D_processed);
        //max_content = mdata2->array_2D_initial[k]*1000*(6561695/total_counts); //mdata1->find_max(mdata1->array_2D_processed);
        /*
        max_content = mdata1->find_max(mdata1->array_2D_processed)/20;

        image_window->DrawArray( mdata1->array_2D_processed, mdata1->size_x, mdata1->size_y, max_content);
        total_counts = (-35000+k*2000 + 7079.97)*200/30000;
        sprintf(buffer," %2.2ffs", total_counts);
        image_window->Write_text(buffer);
        image_window->update();
        image_window->show();

        sprintf(buffer, "%s/%d", work_dir_c, k);
        SaveImage(buffer);
        */
        qApp->processEvents(QEventLoop::ExcludeUserInputEvents);
    } // main loop


    ui->progressBar->setValue(100);
    message("done");

    QApplication::restoreOverrideCursor();
    return;
}

void MainWindow::DataManipulation()
{

    if(!mdata1) {message("There is no data loaded"); return;}


    //

    char work_dir_c[255],fileout[255],work_dir_c2[255], buffer[255];//, buffer[255];
    sprintf(work_dir_c,"D:/Roman/Data/20130531/CS2_electrons/SumImages/");
    //sprintf(work_dir_c2,"D:/Roman/Data/20130531/CS2_electrons/Inverted2");
/*
    QStringList files = QFileDialog::getOpenFileNames(
                            this,
                            "Select one or more files to join",
                            "D:/Musor/CS2Data",
                            "Images (*.*)");
    */
    QStringList files;

    sprintf(fileout,"%s/CS2images2_-1500000.dat.proc.ang", work_dir_c);
    files << fileout;

    int delay = -33000;
    for(int i = 0; i<83;i++) {
        sprintf(fileout,"%s/CS2images2_%d.dat.proc.ang", work_dir_c, delay);
        files << fileout;

        delay = delay + 2000;
    }



    if (!files.isEmpty()){} // file name has been choosen
    else{message("files have not been choosen"); return;}

    message("loading...");

    mdata1->Nskip_lines = 0;
    //sprintf(work_dir, "E:/Roman/Data/AnilineCopy/20120704");

    double *var_array;//, check_number;
    //var_array = new double [7*mdata1->size_y]; // load 5 images and calculate mean and stdev
    int array_size_m = 200*mdata1->size_y;

    var_array = new double [array_size_m];

    for(int i = 0; i<array_size_m;i++){var_array[i] =0;}


    ui->progressBar->setValue(0);
    ui->progressBar->setVisible(1);

////////////////////////// A summary file  ////////////////////////////
    sprintf(fileout,"%s/Summary_ALL_ANG.%s", work_dir_c, "txt");
    ofstream outFile1(fileout, ios::out);
    if(!outFile1){ message("cannot create summary file"); return; }

///////////////////////////////

    int file_start, file_stop; file_start=0; file_stop = files.size();
    for(int k = file_start; k < file_stop; k++){ // main loop

        ui->progressBar->setValue(100*(k-file_start + 1)/(file_stop-file_start + 1));

        QByteArray var = files.at(k).toLatin1();// dir.toLatin1(); // convertion from QString to char*
        sprintf(filename,"%s",var.data()); // get path location
        if (!mdata1->Load_2D(filename,transpose_flag)) {message("cannot open file for reading"); continue;}


        //for(int i = 0; i<mdata1->size_y;i++){var_array[i + k*mdata1->size_y] = mdata1->array_2D_initial[1 + i*mdata1->size_x];} // PES

        for(int i = 0; i<mdata1->size_y;i++){
            var_array[i + k*mdata1->size_y] = mdata1->array_2D_initial[1 + i*mdata1->size_x]; // beta 2
            //var_array[i + k*mdata1->size_y] = mdata1->array_2D_initial[2 + i*mdata1->size_x]; // beta 4
            //if(i%2) var_array[i/2 + k*mdata1->size_y] = mdata1->array_2D_initial[1 + i*mdata1->size_x] + mdata1->array_2D_initial[1 + (i+1)*mdata1->size_x]; // beta 2
        } // beta 2 and beta 4


////////////////////// O U T P U T     D A T A / ////////////////


        outFile1 << filename <<"\t";  // requierd for summary

        for(int i = 0; i<mdata1->size_y;i++){
        //outFile1 << i << "\t" << var_array[i + 5*mdata1->size_y]<< "\t"<< sqrt(var_array[i + 6*mdata1->size_y]/(5*4))<<endl;
            outFile1 << var_array[i + k*mdata1->size_y]<< "\t" ;
        }
        outFile1 << endl;
    } // main loop

    outFile1.close();
    delete [] var_array;

    ui->progressBar->setValue(100);
    message("done");
}


void MainWindow::Test()
{
    getPath();
    if(!file_status){return;}
    //if(mdata1) delete mdata1;
    if(!mdata1) mdata1 = new MData();

    if (!mdata1->Load_1D(filename)) {message("cannot open file for reading"); return;}


    if(mspectr) delete mspectr;
    mspectr = new MSpectrum(this);
    mspectr->setWindowFlags(mspectr->windowFlags() | Qt::Window);
    mspectr->n_points = mdata1->delay_steps;
    mspectr->SetData(mdata1->time_profile);
    mspectr->show();

}

void MainWindow::MDisable()
{
    //ui->menuBar;
}

void MainWindow::SaveSettings()
{
    char buffer[255];
    sprintf(buffer, "user.set");

    ofstream outFile1(buffer, ios::out);
    if(!outFile1){ message("can't open settings file"); return; }

    nL = ui->spinBox_nL->value();
    x_pos = ui->spinBox_xc->value(); y_pos = ui->spinBox_yc->value(); c_radius = ui->spinBox_dr->value();
    int r, g, b;
    color1.getRgb(&r, &g, &b);

    outFile1 << "\t" << Skip_lines << endl;
    outFile1 << "\t" << x_pos << "\t" << y_pos << "\t" << c_radius << endl;
    outFile1 << "\t" << nL << "\t" << n_r << "\t" << n_phi << "\t" << n_basisf << endl;
    outFile1 << "\t" << filename << endl;
    outFile1 << "\t" << current_folder << endl;
    outFile1 << "\t" << bit_depth << "\t" << binary_width << "\t" << binary_height << endl;
    outFile1 << "\t" << Lfile_type << "\t" << n_quadrant << endl;
    for (int i=0;i<4;i++) outFile1 << "\t" << saving_flags[i];    
    outFile1 << endl;
    outFile1 << "\t" << transpose_flag << "\t" << subtractbkg_flag<< endl;
    outFile1 << "\t" << lthickness << endl;
    outFile1 << "\t" << r << "\t" << g << "\t" << b << endl;
    //outFile1 << "\t" << PES << "\t" << ANG << "\t" << PROC << "\t" << INV << endl;
    outFile1.close();

    message("settings are saved");
}

void MainWindow::LoadSettings()
{
    char buffer[255];
    //char username = 'user';
    sprintf(buffer, "user.set");
    int r, g, b;

    ifstream fp;
    fp.open(buffer);
    if (!fp.is_open())
      {
        message("cannot open settings file");
        return;
      }
    fp >> Skip_lines;
    fp >> x_pos >> y_pos >> c_radius;
    fp >> nL >> n_r >> n_phi >> n_basisf;
    fp >> filename;
    fp >> current_folder;
    fp >> bit_depth >> binary_width >>binary_height;
    fp >> Lfile_type >> n_quadrant;
    for (int i=0;i<4;i++) fp >> saving_flags[i];
    fp >> transpose_flag >> subtractbkg_flag;
    fp >> lthickness;
    fp >> r >> g >> b;
    //fp >> savePES >> saveANG >> savePROC >> saveINV;
    fp.close();

    color1.setRed(r); color1.setGreen(g); color1.setBlue(b);

    ui->spinBox_nL->setValue(nL); ui->spinBox_dr->setValue(c_radius); ui->spinBox_xc->setValue(x_pos); ui->spinBox_yc->setValue(y_pos);
    message("settings are loaded");
    //x_pos, y_pos, c_radius;
    //int               nL, n_r, n_phi, n_basisf;
}

void MainWindow::OptionDialog()
{
    int mbasis[20];
    mbasis[0] = n_r; mbasis[1] = n_phi; mbasis[2] = n_basisf;
    mbasis[3] = bit_depth; mbasis[4] = binary_width; mbasis[5] = binary_height;
    mbasis[6] = Lfile_type; mbasis[7] = n_quadrant;
    mbasis[8] = Skip_lines;
    mbasis[9] = lthickness;
    bool other_flags[2]; other_flags[0] = transpose_flag; other_flags[1] = subtractbkg_flag;
    //Options                 (this, mbasis, saving_flags, other_flags, &color)
    Options *opt = new Options(this, mbasis, saving_flags, other_flags, &color1);
    bool ok = opt->exec() == QDialog::Accepted;
    delete opt;
    if (ok) {
        n_r = mbasis[0]; n_phi = mbasis[1]; n_basisf = mbasis[2]; message("new parameters are set");
        bit_depth = mbasis[3]; binary_width =  mbasis[4]; binary_height = mbasis[5];
        Lfile_type = mbasis[6]; n_quadrant = mbasis[7];
        Skip_lines = mbasis[8];
        transpose_flag = other_flags[0];
        subtractbkg_flag = other_flags[1];
        lthickness = mbasis[9];
        //for(int i=0;i<4;i++) saving_flags[i] = mbasis[7+i];
        if(image_window){image_window->set_circle_par(lthickness, color1);}
    }
    else message("settings cancelled");
}

void MainWindow::DeclareMImageView(){
    if(!image_window){
            // window for displaying the image
        image_window = new MImageView(this);
        image_window->setCircle(ui->spinBox_xc->value(), ui->spinBox_yc->value(),ui->spinBox_dr->value());
        image_window->set_circle_par(lthickness, color1);
        ///////////////////////////////////
        connect( image_window, SIGNAL( circleChanged(int, int, int) ), this, SLOT( setCircle(int, int, int) ) );
        connect( image_window, SIGNAL( ChangePixel(int, int) ), this, SLOT( ModifyPixel(int, int) ) );

        image_window->setWindowFlags(Qt::Window);
        // connecting spin boxes to image_window updates
        connect( ui->spinBox_xc, SIGNAL( valueChanged(int) ), image_window, SLOT( circlePXUpdate(int) ) );
        connect( ui->spinBox_yc, SIGNAL( valueChanged(int) ), image_window, SLOT( circlePYUpdate(int) ) );
        connect( ui->spinBox_dr, SIGNAL( valueChanged(int) ), image_window, SLOT( circleRUpdate(int) ) );


        //image_window->circlePXUpdate(ui->spinBox_xc->value());
        //image_window->circlePYUpdate(ui->spinBox_yc->value());
        //image_window->circleRUpdate(ui->spinBox_dr->value());
    }
}
