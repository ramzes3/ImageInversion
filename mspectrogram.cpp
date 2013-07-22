#include "mspectrogram.h"
#include "ui_mspectrogram.h"

MSpectrogram::MSpectrogram(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MSpectrogram)
{
    ui->setupUi(this);
}

MSpectrogram::~MSpectrogram()
{
    delete ui;
}
