#ifndef MULTIFILEINVERSION_H
#define MULTIFILEINVERSION_H

#include <QObject>
#include <QThread>
#include <QStringList>

class MData;
class PBInversion;

class MultiFileInversion : public QThread
{
    Q_OBJECT

public:
    explicit MultiFileInversion(QObject *parent = 0);


    void initialize(QStringList, MData* , PBInversion*, int, int, bool*);
    void invertImages();

signals:
    void resultReady(QString, int);
    void resultCurrent(int);

public slots:


private:

    QStringList files;
    MData *mdata1;
    PBInversion *invers1;
    //int x_c, y_c, dr;
    //int nl, odd;
    //int n_basisf, n_r, n_phi;
    int n_quadrant, file_type;
    bool saving_flags[4], subtractbkg_flag;

protected:
    void run();
};

#endif // MULTIFILEINVERSION_H
