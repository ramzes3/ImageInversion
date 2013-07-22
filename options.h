#ifndef OPTIONS_H
#define OPTIONS_H

#include <QDialog>

namespace Ui {
    class Options;
}

class Options : public QDialog
{
    Q_OBJECT

public:
    explicit Options(QWidget *parent = 0, int *par = NULL, bool *flags= NULL, bool *flags2= NULL, QColor *color = NULL );
    ~Options();

    static int GetOptions();    

private:
    Ui::Options *ui;
    int *ovar;
    bool *bool_var, *bool_var2;
    QColor *var_color;

private slots:
    void getColor();

};

#endif // OPTIONS_H
