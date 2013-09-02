#include "options.h"
#include "ui_options.h"
#include <QColorDialog>

Options::Options(QWidget *parent, int *par, bool *saving_flags, bool *flags2, QColor *color) :
    QDialog(parent),
    ui(new Ui::Options)
{    
    ui->setupUi(this);
    //ui->tab_5->setEnabled(0);
    ui->spinBox_r->setValue(par[0]);
    ui->spinBox_phi->setValue(par[1]);
    ui->spinBox_Rbasis->setValue(par[2]);

    ui->spinBox_bitdepth->setValue(par[3]);
    ui->spinBox_width->setValue(par[4]);
    ui->spinBox_height->setValue(par[5]);                

    switch (par[6])
    {
    case 1: ui->radioButton_ASCII->setChecked(1);
        break;
    case 2: ui->radioButton_binary->setChecked(1);
        break;
    default: ui->radioButton_ASCII->setChecked(1);
        break;
    }

    ui->spinBox_quadrant->setValue(par[7]);
    ui->spinBox_skiplines->setValue(par[8]);
    ui->doubleSpinBox_thickness->setValue(par[9]);
    //  checkboxes   
    //
    ui->checkBox_orig->setChecked(saving_flags[0]);
    ui->checkBox_inv->setChecked(saving_flags[1]);
    ui->checkBox_pes->setChecked(saving_flags[2]);
    ui->checkBox_betas->setChecked(saving_flags[3]);
    ui->checkBox_transpose->setChecked(flags2[0]);
    ui->checkBox_subtract_bkg->setChecked(flags2[1]);
    //

    ovar = &par[0];
    bool_var = &saving_flags[0];
    bool_var2 = &flags2[0];

    var_color = color;

    connect( ui->pushButton_getcolour, SIGNAL( clicked() ), this, SLOT( getColor()) );
}

Options::~Options()
{
    ovar[0] = ui->spinBox_r->value();
    ovar[1] = ui->spinBox_phi->value();
    ovar[2] = ui->spinBox_Rbasis->value();
    ovar[3] = ui->spinBox_bitdepth->value();
    ovar[4] = ui->spinBox_width->value();
    ovar[5] = ui->spinBox_height->value();

    if (ui->radioButton_ASCII->isChecked()){ ovar[6]=1;}
    else {ovar[6]=2;}
    ovar[7] = ui->spinBox_quadrant->value();
    ovar[8] = ui->spinBox_skiplines->value();
    ovar[9] = ui->doubleSpinBox_thickness->value();

    bool_var[0] = ui->checkBox_orig->isChecked();
    bool_var[1] = ui->checkBox_inv->isChecked();
    bool_var[2] = ui->checkBox_pes->isChecked();
    bool_var[3] = ui->checkBox_betas->isChecked();

    bool_var2[0] = ui->checkBox_transpose->isChecked();
    bool_var2[1] = ui->checkBox_subtract_bkg->isChecked();

    delete ui;
}

void Options::getColor(){

    QColor color = QColorDialog::getColor(var_color[0], this, "Choose a colour") ;
    if(color.isValid()) var_color->setRgb(color.rgb());

}

int Options::GetOptions()
{
    return 1;
}
