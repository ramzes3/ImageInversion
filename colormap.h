#ifndef COLORMAP_H
#define COLORMAP_H

#include <QRgb>

#include <iostream>
#include <string>
#include <sstream>

//const int rgb_binning = 255;

struct colormap {
  char title[255];
  int red[256], green[256], blue[256];
  QRgb value[256];
};
/*
void load_color_map(char* filename, colormap* var)
{
    FILE* pFile; pFile = NULL;
    char buffer, str[255];
    pFile = fopen (filename , "r");
    if (pFile == 0) {return;}
    buffer = fscanf(pFile, "%s", str);
    strcpy(var->title, str);
    int col1, col2, col3;
    for (int i=0;i<255;i++){
        do { buffer = getc(pFile); } while (buffer != '\(');
        buffer = fscanf(pFile, "%d", &col1);
        buffer = getc(pFile); buffer = getc(pFile);
        buffer = fscanf(pFile, "%d", &col2);
        buffer = getc(pFile); buffer = getc(pFile);
        buffer = fscanf(pFile, "%d", &col3);
        var->red[i] = col1; var->green[i] = col2; var->blue[i] = col3;
        var->value[i] = qRgb(col3, col2, col3);
        buffer = getc(pFile);
        buffer = getc(pFile);
    }
    fclose (pFile);
}

*/
#endif // COLORMAP_H
