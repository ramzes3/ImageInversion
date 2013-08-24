#include "mimageview.h"
#include "ui_mimageview.h"

#include "math.h"
#include <algorithm>

#include <QtGui>
#include <QLabel>

MImageView::MImageView(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MImageView)
{
    ui->setupUi(this);
    geometry_yoffset = 26; // offset from the top of the window to avoid overlap between image and LineEdit bars //
    //this->setFixedSize(512, 538); // setting widget size
    Qt::WindowFlags flags = this->windowFlags();
    this->setWindowFlags(flags | Qt::CustomizeWindowHint | Qt::WindowStaysOnBottomHint);
    //Qt::WindowFlags flag = Qt::WindowStaysOnBottomHint;
    //this->setWindowFlags(flag);
    // declaring pointers as zero objects
    m_image = NULL; array_view = NULL;
    circle_st = false;
    aspect_ratio = 1;
    x_scale = 1; y_scale = 1;
    el_color = QColor(Qt::green);
    el_thickness = 1;

    //QRect geom1 = QRect(0, 0, 500, (500 + geometry_yoffset));
    //this->setGeometry(geom1);


    //char buffer[100];
    //sprintf(buffer,"%2.2f",pixel_scale);
    //ui->lineEdit_bincontent->setText(buffer);

    // imageLabel is a container for displaying 2D distribution //
    imageLabel = new QLabel(this);
    imageLabel->setMouseTracking(true);
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setScaledContents(true);

    imageLabel->setPalette(QPalette(QColor(250, 250, 250)));
    imageLabel->setAutoFillBackground(true);

    QRect LabelG; LabelG.moveTop(geometry_yoffset);
    imageLabel->setGeometry(LabelG);

    // loadading a colour map for 2D image
    //char colormap_file[] = "/home/roman/Prog/QtBased/MImageViewer/colourmap.py"; // Not a good solution. It has to be changed
    //char colormap_file[] = "colourmap.py"; // Not a good solution. It has to be changed
    //load_color_map(colormap_file, &newmap);

    //connect( ui->horizontalSlider_bottom, SIGNAL( sliderReleased() ), this, SLOT( ReDraw() ) );
    connect( ui->horizontalSlider_top, SIGNAL( sliderReleased() ), this, SLOT( ReDraw() ) );
    connect( ui->horizontalSlider_bottom, SIGNAL( valueChanged(int) ), this, SLOT( resizeWidget(int) ) );
    // set default position of the slide bar
    //ui->horizontalSlider_bottom->setValue(150);
    // unzoom button
    ui->pushButton_unzoom->setIcon(QIcon("zoom_out.bmp"));
    connect( ui->pushButton_unzoom, SIGNAL( clicked() ), this, SLOT( Zoom_out() ) );
}

MImageView::~MImageView()
{
    delete ui;
    if( m_image ) delete m_image;
    if(imageLabel) delete imageLabel;
    if(array_view) delete [] array_view;
}

// Fill 2D histogram with appropriate colourmap and paint it in imageLabel canvas
void MImageView::DrawArray(double* array, int size_x, int size_y, double max_value)
{
    //QPainter painter(this);
    if(!array) return;
    if (x_scale == 1 ) this->setFixedSize(size_x, size_y + geometry_yoffset);
    //QRect geom1 = QRect(this->geometry().x(), this->geometry().y(), size_x, (size_y + geometry_yoffset));
    //this->setGeometry(geom1);
    image_width  = size_x;
    image_height = size_y;

    m_value = max_value;
    aspect_ratio = (double) size_y / (double) size_y;

    //this->size().scale(image_width, image_height + geometry_yoffset, Qt::AspectRatioMode);

    if (m_image) delete m_image;
    m_image = new QImage(size_x, size_y, QImage::Format_RGB32);
    if (array_view) delete [] array_view;
    array_view = new double [size_x * size_y];

    for (int i=0;i<size_y;i++){  //  y coordinate / vertical
        for (int j=0;j<size_x;j++){ // x coordinate / horizontal
            array_view[j + i*size_x] = array[j + i*size_x];
        }
    }

    imageLabel->setPixmap(QPixmap::fromImage(m_image[0]));
    m_widgetSize = this->size() - QSize(0, geometry_yoffset);
    //QSize labelSize = imageLabel->size();
    x_scale = (double) (m_widgetSize.width() ) / (double) (image_width);
    y_scale = (double) (m_widgetSize.height() ) / (double) (image_height);
    //Q_ASSERT(imageLabel->pixmap());
    //imageLabel->resize(widgetSize);
    //x_scale = (double) (size_x - 1) / (double) (imageLabel->size().width() - 1);
    //y_scale = (double) (size_y - 1) / (double) (imageLabel->size().height() - 1);
    ReDraw();
}

QRgb MImageView::colourMap(int intensity, int max_value){
    /* Implementation of JET colour map from http://www.metastine.com/?p=7
     *
     */
    double norm  = intensity/(double) max_value;
    double fourValue = 4 * norm;
    int red   = std::min(fourValue - 1.5, -fourValue + 4.5)*255;
    int green = std::min(fourValue - 0.5, -fourValue + 3.5)*255;
    int blue  = std::min(fourValue + 0.5, -fourValue + 2.5)*255;

    if (red < 0)   red = 0;
    if (green < 0) green = 0;
    if (blue < 0)  blue = 0;

    return qRgb(red, green, blue); // qt rgb type
}

void MImageView::ReDraw()
{
    if (!array_view) return;
    //double cvalue, m_value_new;
    //int pos_bott = ui->horizontalSlider_bottom->value();
    int pos_top = ui->horizontalSlider_top->value(); //ui->horizontalSlider_top->value();
    int intensity ;
    //m_value = 2;
    QRgb value;
    for (int i=0;i<image_height;i++){  //  y coordinate / vertical
        for (int j=0;j<image_width;j++){ // x coordinate / horizontal
            //cvalue = log10(10 + (array_view[j + i*image_width] / m_value)*100);
            //m_value_new = log10(110);
            intensity = (int)( 255.0 * array_view[j + i*image_width] / m_value);
            //int intensity = (int)( 255.0 * cvalue / m_value_new);
            intensity = intensity * ( 255.0 / pos_top );
            if ( intensity > 255.0 || intensity < 0.0) value = qRgb(189, 149, 39);
            else {value = colourMap(intensity, m_value);}

            m_image->setPixel( QPoint(j,i), value );
        }
    }
    //Write_text("fdfddf");
    imageLabel->setPixmap(QPixmap::fromImage(m_image[0]));
    //QSize widgetSize = this->size() - QSize(0, geometry_yoffset);
    //QSize labelSize = imageLabel->size();
    //x_scale = (double) (size_x - 1) / (double) (widgetSize.width() - 1);
    //y_scale = (double) (size_y - 1) / (double) (widgetSize.height() - 1);
    Q_ASSERT(imageLabel->pixmap());
    update();
}

void MImageView::Write_text(const char* text){

    QPainter p(m_image);
    p.setPen(QPen(Qt::white));
    p.setFont(QFont("Times", 60, QFont::Bold));
    p.drawText(m_image->rect(),Qt::AlignBottom, text);
}

void MImageView::mouseMoveEvent(QMouseEvent *event)
{
    // checking if mouse is within the image frame
    int posx = ceil( (event->x()) / x_scale );
    if( (posx < 0) || (posx > image_width-1) ) return;

    int posy = ceil( (event->y() - geometry_yoffset) / y_scale );
    if( (posy < 0) || (posy > image_height-1) ) return;


    char buffer[100];
    sprintf(buffer,"%d", posx);
    ui->lineEdit_xpos->setText(buffer);

    sprintf(buffer,"%d", posy);
    ui->lineEdit_ypos->setText(buffer);

    if(m_image){
        //QPoint m_point(posx, posy);
        //QColor m_color(m_image->pixel(m_point));
        //sprintf(buffer,"%d, %d, %d", m_color.red(), m_color.green(), m_color.blue());
        sprintf(buffer,"%2.2f",array_view[ posx + image_width*posy ]);
        ui->lineEdit_bincontent->setText(buffer);
    }

    if(event->buttons() == Qt::RightButton) {
       // int im_height = this->width() - geometry_yoffset;
        el_radius_x = abs( event->x() - cursor_x); el_radius_y = abs( event->y() - cursor_y);
        int radius = (int) sqrt( pow(el_radius_x,2) + pow(el_radius_y, 2) );
        circleRUpdate(radius / sqrt(y_scale*x_scale));
       // int radius = (int) sqrt( pow(el_radius_x,2) + pow(el_radius_y, 2) );
        //if(cursor_x > this->width/2)
         //   if( (cursor_x + radius) > this->width() )
        circle_st = true;
        emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );
        update();
    }

}

void MImageView::mousePressEvent(QMouseEvent *event)
{
    // checking if mouse is within the image frame
    int posx = (double) (event->x()) / x_scale;
    if( (posx < 0) || (posx > image_width-1) ) return;

    int posy = (double) (event->y() - geometry_yoffset) / y_scale;
    if( (posy < 0) || (posy > image_height-1) ) return;

    //char buffer[100];
    //sprintf(buffer,"%s","piff puff");
    //ui->lineEdit_bincontent->setText(buffer);

    // only with a right mouse button pressed you can draw circle
    //if (event->button() == Qt::LeftButton) { cursor_x = event->x(); cursor_y = event->y();circleRUpdate(c_radius); update();}
    if (event->button() == Qt::RightButton) {
        cursor_x = event->x(); cursor_y = event->y(); circleRUpdate(ceil(c_radius/sqrt(y_scale*x_scale)));
        emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );
    }
}

void MImageView::mouseDoubleClickEvent ( QMouseEvent * event )
{
    if(event->buttons() == Qt::LeftButton)  emit ChangePixel(ceil(event->x()/x_scale), ceil((event->y() - geometry_yoffset)/y_scale));
}

void MImageView::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    if(!m_image) return;
    QSize widgetSize = this->size() - QSize(0, geometry_yoffset);
    QRect frame( 0, geometry_yoffset, widgetSize.width(), widgetSize.height());
    painter.drawImage(frame, m_image[0]);

    if(circle_st){
        QPen pen(el_color, el_thickness, Qt::DashDotLine, Qt::RoundCap, Qt::RoundJoin);
        painter.setPen(pen);
        QPoint centre_p(cursor_x, cursor_y);
        painter.drawEllipse( centre_p, c_radius, c_radius );
    }
}

void MImageView::resizeWidget(int zoom_level)
{
    // recover original cursor position
    circle_st = false;
    //c_radius = c_radius/sqrt(y_scale*x_scale);
    //cursor_x = cursor_x * x_scale; cursor_y = (cursor_y-geometry_yoffset ) * y_scale + geometry_yoffset;


    // definining minimum image width by looking at the most right widget element
    double min_width = ui->horizontalSlider_top->geometry().x() + ui->horizontalSlider_top->width(); //move it as more general definition;
    double min_scale = min_width/ (double) image_width;
    // calculating scaling factor
    double scale = min_scale + (double) zoom_level / (150.0 / (1-min_scale));
    //x_scale = scale; y_scale = scale;
    //resizing widget
    this->setFixedSize(image_width*scale, (image_height*scale + geometry_yoffset));

    // read out size of the widget and find out scaling factor
    m_widgetSize = this->size() - QSize(0, geometry_yoffset);
    x_scale = (double) (m_widgetSize.width() ) / (double) (image_width);
    y_scale = (double) (m_widgetSize.height() ) / (double) (image_height);

    // new cursor scaling
    //cursor_x = cursor_x / x_scale; cursor_y = (cursor_y-geometry_yoffset ) / y_scale + geometry_yoffset;
    //c_radius = c_radius * sqrt(y_scale*x_scale);
    // update circled area
    //emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );
    //circleRUpdate(ceil(c_radius/sqrt(y_scale*x_scale)));
}

void MImageView::resizeEvent(QResizeEvent *event)
{
    //double aspect_ratio = event->oldSize().width() / event->size().width();

    //QSize newSize(event->size().width(), (event->oldSize().height() / aspect_ratio));
    //QSize widgetSize = event->size() - QSize(0, geometry_yoffset);
    //QSize labelSize = imageLabel->size();
    //if (aspect_ratio) this->resize(newSize);

    //x_scale = (double) (image_width ) / (double) (widgetSize.width() );
    //y_scale = (double) (image_height ) / (double) (widgetSize.height() );

    //char buffer[100];
    //sprintf(buffer,"%2.2f", x_scale);
    //ui->lineEdit_bincontent->setText(buffer);
}

/*
void MImageView::load_color_map(char* filename, colormap* var)
{
    FILE* pFile; pFile = NULL;
    char buffer, str[255];
    pFile = fopen (filename , "r");
    if (pFile == 0) {return;}
    buffer = fscanf(pFile, "%s", str);
    strcpy(var->title, str);
    int col1, col2, col3;
    for (int i=0;i<256;i++){
        do { buffer = getc(pFile); } while (buffer != '\(');
        buffer = fscanf(pFile, "%d", &col1);
        buffer = getc(pFile); buffer = getc(pFile);
        buffer = fscanf(pFile, "%d", &col2);
        buffer = getc(pFile); buffer = getc(pFile);
        buffer = fscanf(pFile, "%d", &col3);
        var->red[i] = col1; var->green[i] = col2; var->blue[i] = col3;
        var->value[i] = qRgb(col1, col2, col3);
        buffer = getc(pFile);
        buffer = getc(pFile);
    }
    fclose (pFile);
}
*/

void MImageView::setCircle(int posx, int posy, int radius){
    cursor_x = posx*x_scale;
    cursor_y = posy*y_scale + geometry_yoffset;
    c_radius = radius * sqrt(y_scale*x_scale);
}

void MImageView::circlePXUpdate(int posx)
{
    cursor_x = posx*x_scale;
    circle_st = true;

    if( (c_radius + cursor_x) > (image_width*x_scale) ) {c_radius = image_width*x_scale - cursor_x;  emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}
    if( (cursor_x - c_radius) < 0 ) {c_radius = cursor_x;  emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}

    update();
}
void MImageView::circlePYUpdate(int posy)
{
    cursor_y = posy*y_scale + geometry_yoffset;    
    circle_st = true;
    if( (c_radius + cursor_y) > (image_height*y_scale + geometry_yoffset) ) {c_radius = image_height*y_scale + geometry_yoffset - cursor_y; emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}
    if( (cursor_y - c_radius) < geometry_yoffset) {c_radius = cursor_y - geometry_yoffset; emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}

    update();
}
void MImageView::circleRUpdate(int radius)
{
    c_radius = radius * sqrt(y_scale*x_scale);
    circle_st = true;
    if( (c_radius + cursor_x) > (image_width*x_scale) ) {c_radius = image_width*x_scale - cursor_x;  emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}
    if( (cursor_x - c_radius) < 0 ) {c_radius = cursor_x;  emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}
    if( (c_radius + cursor_y) > (image_height*y_scale + geometry_yoffset) ) {c_radius = image_height*y_scale + geometry_yoffset - cursor_y; emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}
    if( (cursor_y - c_radius) < geometry_yoffset) {c_radius = cursor_y - geometry_yoffset; emit circleChanged( ceil(cursor_x/x_scale), ceil((cursor_y - geometry_yoffset)/y_scale) , ceil(c_radius/sqrt(y_scale*x_scale)) );}

    update();
}

void MImageView::update_scaling()
{
    QSize widgetSize = this->size() - QSize(0, geometry_yoffset);
    //QSize labelSize = imageLabel->size();
    x_scale = (double) (widgetSize.width()) / (double) (image_width);
    y_scale = (double) (widgetSize.height()) / (double) (image_height);

    char buffer[100];
    sprintf(buffer,"%f", x_scale);
    ui->lineEdit_bincontent->setText(buffer);
}

void MImageView::Zoom_out(){ ui->horizontalSlider_bottom->setValue(150); }

void MImageView::set_circle_par(double thickness, QColor mcolor){
   el_thickness = thickness;
   el_color = mcolor;
}

void MImageView::changeEvent(QEvent *e)
{
    QWidget::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    case QEvent::Resize:
        update_scaling();
        break;
    default:
        break;
    }
}
