#ifndef MSPECTRUM_H
#define MSPECTRUM_H

#include <QWidget>

namespace Ui {
    class MSpectrum;
}

class MSpectrum : public QWidget {
    Q_OBJECT
public:
    MSpectrum(QWidget *parent = 0);
    ~MSpectrum();

    void SetData(double*);
    double* experimentalPoints;
    int n_points;

protected:
    void changeEvent(QEvent *e);
    void paintEvent(QPaintEvent *event);  // draw a simple spectrum graph

private:
    Ui::MSpectrum *ui;


};

#endif // MSPECTRUM_H
