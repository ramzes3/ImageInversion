#include "mspectrum.h"
#include "ui_mspectrum.h"

#include <QtGui>

MSpectrum::MSpectrum(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MSpectrum)
{
    ui->setupUi(this);
    experimentalPoints = NULL;
    n_points = 256;
}

MSpectrum::~MSpectrum()
{
    delete ui;
    if(experimentalPoints) delete [] experimentalPoints;
}

void MSpectrum::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    QSize widgetSize = this->size();
    //QRect frame( 0, geometry_yoffset, widgetSize.width(), widgetSize.height());
    //painter.drawImage(frame, m_image[0]);
    //char buffer[100];
    //sprintf(buffer,"%s","painted");
    //ui->lineEdit_bincontent->setText(buffer);
    int axis_offset = 10;
    painter.setPen(QPen(Qt::black, 1, Qt::SolidLine));
    painter.drawLine(QPoint(axis_offset,widgetSize.height()-axis_offset), QPoint(widgetSize.width()-axis_offset, widgetSize.height()-axis_offset));
    painter.drawLine(QPoint(axis_offset,widgetSize.height()-axis_offset), QPoint(axis_offset, axis_offset));
    //QPen pen(Qt::green, 1, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
    QPen pen(Qt::red, 3, Qt::SolidLine);
    painter.setPen(pen);
    //QPoint centre_p(5, 5), y_end(5, widgetSize.height()-5);
    //painter.drawLine(centre_p, y_end);
    QPoint st_point, nd_point;
    st_point = QPoint(axis_offset, widgetSize.height()-axis_offset - (widgetSize.height()-2*axis_offset)*experimentalPoints[0]);
    for(int i=1;i<n_points;i++){
        //experimentPoint = QPoint(sin(i),i);
        nd_point = QPoint(axis_offset + i*(widgetSize.width()-2*axis_offset)/n_points, axis_offset + (widgetSize.height()-2*axis_offset) - (widgetSize.height()-2*axis_offset)*experimentalPoints[i]);
        painter.drawPoint(st_point);
        painter.drawLine(st_point, nd_point);
        st_point = nd_point;
    }

}

void MSpectrum::SetData(double* data)
{
    if(experimentalPoints) delete [] experimentalPoints;
    experimentalPoints = new double [n_points];
    double max_value = 0;
    for(int i=0;i<n_points;i++){
        experimentalPoints[i] = data[i];
        if(max_value < experimentalPoints[i]) max_value = experimentalPoints[i];
    }
    for(int i=0;i<n_points;i++){
        experimentalPoints[i] /= max_value;
    }
}

void MSpectrum::changeEvent(QEvent *e)
{
    QWidget::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}
