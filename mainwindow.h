#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui>
#include <QMainWindow>
#include <QTimer>
//#include <QThread>

class MData;
class MImageView;
class PBInversion;
class MPop;
class MSpectrum;

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
     QThread workerThread;
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
        void Load_Image_ASCII();  // open ASCII file, load image into MData object
        void Add_Image_ASCII(); // open ASCII file, load image into MData object and add it to existing one
        void Load_Image_Binary(); // opens binary images
        void getPath();          //  get a directory full file path and write it into filename variable
        void setCircle(int, int, int);  // slot which is used to update area of interest
        void ModifyPixel(int, int);  // modify pixel value of mdata1
        void InvertImage();  // Invert image currently using pBasex code
        void InvertImage1();  // Inverting image but the pBasex class has to be decleared in advance
        void Declear_pBasex(); // creating pBasex class
        void SetBasisDim();  // Set the dimensions of the PBasex basis
        void ShowSpectrum();  // Showing velocity or energy spectrum
        void GenerateBasis(); // Generate basis function for PBasex
        void message(const char*);     // write a message in the message line
        void merror(const char* ); // writing an error message in the error box
        //void message(string);     // write a message in the message line
        void Analise();         // automatic data analisys (of the spesific files in the folder)
        void AnaliseSelected(); // automatic analisys of selected files
        void Transpose();      // transpose image
        void SaveAll();         // save the result of image inversion from mdata1 class
        void SaveSelected(); // Saves selected files, chosen from settings menu
        void SaveImage();      // asking to choose file name then uses function void SaveImage(char*);
        void SaveImage(char*);      // saving image as jpg file from image_window
        void SaveSettings();    // saving imortant parameters in the settings file
        void LoadSettings();    // Load important parameters from the settings file
        void OptionDialog();    // start an option dialog for the program settings

        void DataManipulation(); // collecting data and averaging some random stuff
        void Symmetrise();  // symmetrise the image over all 4 quarters
        void GetPolarImage();
        void GetQuarter(); // taking only a single quarter of the image and copying it into the other quarters


        void GetSkipLines();   // get number of lines to be skipped from ASCII file
        void CheckSelectedArea(int ); // check if flag was cheched on full or circle area of the image
        void ReloadImage();  // copies original image into processed
        void SetBkg();
        void ShowBkg();
//        void ShowImage();
        void Subtract_bkg();
        void Add_bkg();

        void About();
        void MUpdate();
        void MDisable(); // disabling all the buttons of the main window

        void Test();
        void ConvertImages(); // converting ASCII data into jpg
        void ConvertQstringToChar(QString, char*);
        void wakeup(); //
        void invProgress(QString, int);

protected:
    void changeEvent(QEvent *e);

private:
    Ui::MainWindow *ui;

    void DeclareMImageView();  // Deaclaring MImageView object, repetitive action in the separate function

    char              filename[512], filesave[512], current_folder[512];
    //QString           current_folder;
    bool              file_status; // 1 - file is loaded, 0 - file is not loaded
    int               Skip_lines;  // number of line to skip in image ASCII file before reading it out
    int               x_pos, y_pos, c_radius;  // circular area of interest

    int               nL, n_r, n_phi, n_basisf;  // dimensions of the pBasex
    int               bit_depth, binary_width, binary_height; // seetings for loading binary file
    int               Lfile_type, n_quadrant;
    bool              saving_flags[9], transpose_flag, subtractbkg_flag;
    QColor            color1;
    int               lthickness;

    int               area_flag; // 0 - full area, 1 - circle area
    MData*            mdata1; //  class contains all the data and data loading functions
    MImageView*       image_window;  //displays 2D VMI image original or inverted
    PBInversion*      invers;        //  PBasex inversion procedures
    MPop*             pop;  // POP algorithm for inversion
    MSpectrum*        mspectr;       //  displays velocity or energy spectrum after inversion    
    QSignalMapper*    msignal_mapper;

    QTimer*           timer_loading;  // timer for loading the image
    QStringListModel *typeModel; // for displaying error messages
    QStringList       error_list; //list of the errors

};

#endif // MAINWINDOW_H
