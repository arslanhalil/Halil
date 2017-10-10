/*
 * Image.h - header file of the Image library which defines 
 *           a new class "Image" and the associated member functions
 *
 * Note: 
 *   This is a simple C++ library for image processing.
 *   This library can only read in PGM/PPM format images. 
 */

#ifndef IMAGE_H
#define IMAGE_H

#include <iostream>
using namespace std;

#define PGMRAW   1                     // magic number is 'P5'
#define PPMRAW   2                     // magic number is 'P6'
#define PGMASCII 3                     // magic number is 'P2'
#define PPMASCII 4                     // magic number is 'P3'
#define GRAY     10                    // gray-level image
#define BINARY   11                    // binary image

#define PI       3.1415926

// coefficients for D4 wavelet transform
#define C0 0.4829629131445341
#define C1 0.8365163037378079
#define C2 0.2241438680420134
#define C3 -0.1294095225512604

class Image {
  friend ostream & operator<<(ostream &, Image &);
 public:
  // constructors and destructor
  Image();                             // default constructor 
  Image(int,                           // constructor with row
    int);                              // column (for grayscale image)
  Image(int,                           // constructor with row
    int,                               // column
    int);                              // type (use PGMRAW, PPMRAW, 
                                       //           PGMASCII, PPMASCII)
  Image(const Image &);                // copy constructor 
  ~Image();                            // destructor 

  // create an image
  void createImage();                  // create an image, parameters all set
  void createImage(int,                // create an image with row
           int,                        // column
           int t=PGMRAW);              // and type, default is PGMRAW
  void initImage(float init=0.0);      // initiate the pixel value of an image
                                       // the default is 0.0

  // get and set functions
  int getRow() const;                  // get row # / the height of the image 
  int getCol() const;                  // get col # / the width of the image 
  int getChannel() const;              // get channel number of the image
  int getType() const;                 // get the image type 
  int getMaximum() const;              // get the maximum pixel value
  int getMinimum() const;              // get the mininum pixel value
  int getSetmax() const;               // get the indicator
  Image getRed() const;                // get red channel
  Image getGreen() const;              // get green channel
  Image getBlue() const;               // get blue channel
  int* getHist() const;                // get the histogram of an image

  void setRow(int);                    // set row number 
  void setCol(int);                    // set column number 
  void setChannel(int);                // set the number of channel
  void setType(int);                   // set the image type
  void setMaximum(int);                // set the maximum pixel value
  void setRed(Image &);                // set the red channel
  void setGreen(Image &);              // set the green channel
  void setBlue(Image &);               // set the blue channel

  void readImage(char *);              // read image from a file
  void writeImage(char *);             // write image to a file

  // matrix manipulations
  Image transpose();                   // transpose of an image
  Image inverse();                     // inverse of an image
  Image subImage(int, int, int, int);  // clip an image 
                                       // (startrow, startcol, endrow, endcol) 
  // operator overloading functions
  float & operator()(int, 
             int, 
             int k = 0) const;             // operator overloading, default k=0
  const Image operator=(const Image &);    // = operator overloading
  Image operator+(const Image &) const;    // overloading + operator
  Image operator-(const Image &) const;    // overloading - operator
  Image operator*(const Image &) const;    // overloading pixelwise *
  Image operator/(const Image &) const;    // overloading pixelwise division
  Image operator->*(const Image &) const;  // overloading ->* operator 
                                           // (matrix multiplication) 
  
  // logic operators
  Image AND(Image &);                  // logic AND
  Image OR(Image &);                   // logic OR
  Image NOT();                         // logic NOT
  Image XOR(Image &);                  // logic XOR
  int XOR(int, int);                   // logic XOR

  // Fourier transforms
  void fft(Image &,                    // forward FT with Fourier spectrum
       Image &);                       // phase
  void ifft(Image &,                   // inverse FT with Fourier spectrum
        Image &);                      // phase
  void fftifft(Image &,                // called by fft() or ifft()
           Image &,                    // phase
           int);                       // 1: perform fft, -1: perform ifft

  // 1-D wavelet transform
  Image wt(int, int);
  void daub4(int, int);
  Image wt2(int, int);
  void daub42(int, int);
  // point-based image enhancement
  Image logtran();                     // log transformation

  // color model convertion
  Image RGB2HSI();                     // convert color image from RGB to HSI 
                                       // (H-1st chan, S-2nd, I-3rd) 
  Image HSI2RGB();                     // convert color image from HSI to RGB

  // morphological operators
  // the origin is specified as (int origRow, int origCol)
  // foreground pixels are white (greater than 0)
  Image erode(Image &, int, int);      // erosion given the origin of se
  Image dilate(Image &, int, int);     // dilation given the origin of se
  Image open(Image &, int, int);       // opening given the origin of se
  Image close(Image &, int, int);       // closing given the origin of se

 private:
  int row;                  // number of rows / height 
  int col;                  // number of columns / width 
  int channel;              // nr of channels (1 for gray, 3 for color)
  int type;                 // image type (PGM, PPM, etc.)
  int maximum;              // the maximum pixel value
  int setmax;               // indicates if users want to set their own maximum
  float *image;             // image buffer
  void outofMemory();       // output out of memory message

  // the following four functions are used by inverse()
  int findPivot(Image &, int);        // find the row with max abs value in that col
  void switchRow(Image &, int, int);  // switch two rows
  void dividePivot(Image &, int);     // divide that row with the element in that col
  void eliminate(Image &, int);       // eliminate the following columns
};

#endif
