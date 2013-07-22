# -------------------------------------------------
# Project created by QtCreator 2010-05-15T18:59:18
# -------------------------------------------------
TARGET = Iminv

# debug:TARGET = /Iminv/debug/Iminv
TEMPLATE = app\

QT += widgets

 #LIBS += -L/usr/local/lib \
 #-lgsl \
 #-lgslcblas \
 #-lm
LIBS += C:\GSL\lib\libgslcblas.dll.a
LIBS += C:\GSL\lib\libgsl.dll.a
INCLUDEPATH += C:\GSL\include

#LIBS += C:\cygwin\lib\libgsl.dll.a
#LIBS += C:\cygwin\lib\libgslcblas.dll.a
#INCLUDEPATH += C:\cygwin\usr

# debug:TARGET = C:\Documents and Settings\Laser1\Desktop\Roman\GUI\ImageInvertion_2\Iminv\debug
SOURCES += main.cpp \
    mainwindow.cpp \
    mimageview.cpp \
    mdata.cpp \
    pbinversion.cpp \
    mspectrum.cpp \
    mspectrogram.cpp \
    mpop.cpp \
    options.cpp \
    multifileinversion.cpp
HEADERS += mainwindow.h \
    mimageview.h \
    mdata.h \
    pbinversion.h \
    mspectrum.h \
    mspectrogram.h \
    mpop.h \
    options.h \
    multifileinversion.h
FORMS += mainwindow.ui \
    mimageview.ui \
    mspectrum.ui \
    mspectrum.ui \
    mspectrogram.ui \
    dialog1.ui \
    options.ui
