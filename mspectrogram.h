#ifndef MSPECTROGRAM_H
#define MSPECTROGRAM_H

#include <QWidget>

namespace Ui {
    class MSpectrogram;
}

class MSpectrogram : public QWidget
{
    Q_OBJECT

public:
    explicit MSpectrogram(QWidget *parent = 0);
    ~MSpectrogram();

private:
    Ui::MSpectrogram *ui;
};

#endif // MSPECTROGRAM_H
