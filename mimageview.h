#ifndef MMIMAGEVIEW_H
#define MIMAGEVIEW_H

#include <QWidget>

class QLabel;

namespace Ui {
    class MImageView;
}

class MImageView : public QWidget {
    Q_OBJECT
public:
    MImageView(QWidget *parent = 0);
    ~MImageView();

    void DrawArray( double*, int, int, double);
    void Write_text(const char*);
    void set_circle_par(double, QColor);

    QImage          *m_image;

public slots:
    void setCircle(int, int, int); // initialize circle dimensions
    void circlePXUpdate(int); // update of x position of the circle centre
    void circlePYUpdate(int);   // update of y position of the circle centre
    void circleRUpdate(int); // update of the circle radius
    void ReDraw();  // Redraw and image when something changed
    void resizeWidget(int);  // resize the window but keep the aspect ratio
    void Zoom_out();

signals:
    void circleChanged(int, int, int);
    void ChangePixel(int, int);



protected:
    void changeEvent(QEvent *e);
    void mouseMoveEvent(QMouseEvent *event);  // tracking mouse movements on the widget
    void mousePressEvent(QMouseEvent *event);  // function is called when mouse is pressedin the window
    void mouseDoubleClickEvent( QMouseEvent * event );// used to change pixel value
    void paintEvent(QPaintEvent *event);  // reimplementation of the widget paint even to draw a circle / ellipse
    void resizeEvent(QResizeEvent * event);  // update variables when widget size have been changed

private:
    Ui::MImageView *ui;
    QRgb colourMap(double intensity, double max_value);
    void update_scaling(); // update parameters when widget size is changes
    //double           pixel_scale; // ratio of the x/y pixel size of the monitor.

    QLabel          *imageLabel;
    double          *array_view;

    int              geometry_yoffset;
    int              cursor_x, cursor_y;
    int              el_radius_x, el_radius_y, c_radius;
    double           el_thickness;
    QColor           el_color;
    double           x_scale, y_scale, aspect_ratio;
    double           m_value;
    //double           tot_counts;
    int              image_width;
    int              image_height;
    QSize            m_widgetSize;

    bool             circle_st; // 0 no circle, 1 - there is circle
};

#endif // MIMAGEVIEW_H
