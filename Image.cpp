/*
 * Image.cpp - the image library which implements
 *             the member functions defined in Image.h
 */

#include "Image.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
using namespace std;

#define PI 3.1415926

// default constructor
Image::Image()
{
  createImage(0, 0);
}

// constructor for grayscale images
Image::Image(int r, int c)
{
  createImage(r, c);
}

// constructor for grayscale/color images
Image::Image(int r, int c, int t)
{
  createImage(r, c, t);
}

// copy constructor
Image::Image(const Image &img)
{
  int i, j, k;

  row = img.getRow();
  col = img.getCol();
  type = img.getType();
  if (type == PGMRAW || type == PGMASCII)
    channel = 1;
  else if (type == PPMRAW || type == PPMASCII)
    channel = 3;
  else
    cout << "Undefined image type!\n";

  createImage(row, col, type);             // allocate memory
  //  maximum = img.getMaximum();
  //  setmax = img.getSetmax();
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
    image[(i*col+j)*channel+k] = img(i, j, k);
}

// destructor
Image::~Image()
{
  if (image)
    delete [] image;       // free the image buffer
}

// allocate memory for the image
void Image::createImage()
{
  if (type == PGMRAW || type == PGMASCII)
    channel = 1;
  else if (type == PPMRAW || type == PPMASCII)
    channel = 3;
  else
    cout << "Undefined image type!\n";
  setmax = 0;
  maximum = 255;

  image = (float *) new float [row * col * channel];
  if (!image) 
    outofMemory();

  initImage();
}

// allocate memory for the image
void Image::createImage(int nr, int nc, int ntype)
{
  row = nr;
  col = nc;
  type = ntype;
  if (type == PGMRAW || type == PGMASCII)
    channel = 1;
  else if (type == PPMRAW || type == PPMASCII)
    channel = 3;
  else
    cout << "Undefined image type!\n";
  setmax = 0;
  maximum = 255;

  image = (float *) new float [row * col * channel];
  if (!image) 
    outofMemory();

  initImage();
}

// initialize the image
void Image::initImage(float init)
{
  int i;

  for (i=0; i<row*col*channel; i++)
    image[i] = init;
}

// get the image row number
int Image::getRow() const
{
  return row;
}

// get the image column number
int Image::getCol() const
{
  return col;
}

// get the image channel number
int Image::getChannel() const
{
  return channel;
}

// get the image type
int Image::getType() const
{
  return type;
}

// get the maximum pixel value of the image
int Image::getMaximum() const
{
  int tmax = 0;
  int i, j, k;

  if (channel == 3) {
    cout << "getMaximum: "
     << "Can only return maximum pixel value of gray-scale image.\n";
    exit(3);
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      if (tmax < image[i*col+j])
        tmax = (int)image[i*col+j];
  
  return tmax;
}

// get the minimum pixel value of the image
int Image::getMinimum() const
{
  float tmin = 255.0;
  int i, j; 

  if (channel == 3) {
    cout << "getMinimum: "
     << "Can only return minimum pixel value of gray-scale image.\n";
    exit(3);
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      if (tmin > image[i*col+j])
        tmin = image[i*col+j];

  return (int)tmin;
}

// get the indicator
int Image::getSetmax() const
{
  return setmax;
}
      
// get the red channel
Image Image::getRed() const 
{
  Image temp;
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  temp.createImage(row, col);

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i, j) = image[(i*col+j)*channel];

  return temp;
}

// get the green channel
Image Image::getGreen() const
{
  Image temp;
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  temp.createImage(row, col);
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i, j) = image[(i*col+j)*channel+1];

  return temp;
}

// get the blue channel
Image Image::getBlue() const
{
  Image temp;
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  temp.createImage(row, col);
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      temp(i, j) = image[(i*col+j)*channel+2];
  
  return temp;
}

// get the histogram
int * Image::getHist() const
{
  int i, j, k;
  int *hist;
  int maxi, mini;
  
  if (channel == 3) {
    cout << "GETHIST: can't handle color images\n";
    exit(3);
  }
  
  // after histogram equalization, the intensity level might not be within 255
  maxi = getMaximum();
  mini = getMinimum();
  if (mini<0 || maxi>255) {
    cout << "getHist: "
         << "The intensity value is outside [0, 255], rescale.\n";  

    if (maxi == 0.0 && mini == 0.0) {
      for (i=0; i<row; i++)
        for (j=0; j<col; j++)
          image[i*col+j] = 0.0;
    }
    else if (maxi == mini) {
      for (i=0; i<row; i++)
        for (j=0; j<col; j++)
          image[i*col+j] = 255.0;
    }
    else {
      for (i=0; i<row; i++)
        for (j=0; j<col; j++)
          image[i*col+j] = (image[i*col+j]-mini) * 255.0 / (maxi-mini); 
    }
    
    maxi = getMaximum();
  }
  
  // allocate memory according to the new maximum
  hist = (int *) new int [maxi];
  
  for (k=0; k<maxi; k++)
    hist[k] = 0;
    
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      hist[(int)image[(i*col+j)]]++;
      
  return hist;
}

// set the row number
void Image::setRow(int r)
{
  row = r;
}

// set the column number
void Image::setCol(int c)
{
  col = c;
}

// set the image type
void Image::setType(int t)
{
  type = t;
  if (t == PGMRAW || t == PGMASCII)
    channel = 1;
  else if (t == PPMRAW || t == PPMASCII)
    channel = 3;
  else
    cout << "Undefined image type!\n";
}

void Image::setChannel(int c)
{
  channel = c;
}

// set the maximum pixel value of the image
void Image::setMaximum(int m)
{
  maximum = m;
  setmax = 1;
}

// set the red channel
void Image::setRed(Image &red) 
{
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      image[(i*col+j)*channel] = red(i, j);
}

// set the green channel
void Image::setGreen(Image &green) 
{
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      image[(i*col+j)*channel+1] = green(i, j);
}

// set the red channel
void Image::setBlue(Image &blue) 
{
  int i, j;

  if (channel != 3) {
    cout << "This is not a color image\n";
    exit(3);
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      image[(i*col+j)*channel+2] = blue(i, j);
}

// read image from a file                     
void Image::readImage(char *fname)
{
  ifstream ifp;
  char dummy[80];
  unsigned char *img;
  int i, j, k;
  int tmp;

  ifp.open(fname, ios::in | ios::binary);

  if (!ifp) {
    cout << "Can't read image: " << fname << endl;
    exit(1);
  }

  // identify image format
  ifp.getline(dummy, 80, '\n');

  if (dummy[0] == 'P' && dummy[1] == '5') {
    type = PGMRAW;
    channel = 1;
  }
  else if (dummy[0] == 'P' && dummy[1] == '6') {
    type = PPMRAW;
    channel = 3;
  }
  else if (dummy[0] == 'P' && dummy[1] == '2') {
    type = PGMASCII;
    channel = 1;
  }
  else if (dummy[0] == 'P' && dummy[1] == '3') {
    type = PPMASCII;
    channel = 3;
  }
  else {
    cout << "Can't identify image format." << endl;
  }

  // skip the comments
  ifp.getline(dummy, 80, '\n');

  while (dummy[0] == '#') {
    ifp.getline(dummy, 80, '\n');
  }

  // read the row number and column number
  sscanf(dummy, "%d %d", &col, &row);

  // read the maximum pixel value
  ifp.getline(dummy, 80, '\n');
  sscanf(dummy, "%d", &maximum); 

  // read the image data
  img = (unsigned char *) new unsigned char [row * col * channel];
  if (!img) 
    outofMemory();

  // convert the data type from unsigned char to float
  image = (float *) new float [row * col * channel];
  if (!image) 
    outofMemory();

  // added capability to process ASCII format as well
  if (type == PGMRAW || type == PPMRAW) {
    ifp.read((char *)img, (row * col * channel * sizeof(unsigned char)));
    for (i=0; i<row; i++)
      for (j=0; j<col; j++)
        for (k=0; k<channel; k++)
          image[(i*col+j)*channel+k] = (float)img[(i*col+j)*channel+k];
  }
  else {   // ASCII formats
    for (i=0; i<row; i++)
      for (j=0; j<col; j++)
        for (k=0; k<channel; k++) {
          ifp >> tmp;
          image[(i*col+j)*channel+k] = (float)tmp;
        }
  }

  ifp.close();
  delete [] img;
}

// write image buffer to a file
void Image::writeImage(char *fname)
{
  ofstream ofp;
  char dummy[80];
  int i, j, k;
  float *maxi, *mini;
  unsigned char *img;

  ofp.open(fname, ios::out | ios::binary);

  if (!ofp) {
    cout << "Can't write image: " << fname << endl;
    exit(1);
  }

  // Write the format ID
  switch (type) {
  case PGMRAW:
    ofp << "P5" << endl;
    break;
  case PPMRAW:
    ofp << "P6" << endl;
    break;
  case PGMASCII:
    ofp << "P2" << endl;
    break;
  case PPMASCII:
    ofp << "P3" << endl;
    break;
  otherwise:
    cout << "Can't identify image type\n";
  }

  ofp << col << " " << row << endl;

  maxi = (float *) new float [channel];
  mini = (float *) new float [channel];
  for (k=0; k<channel; k++) {
    maxi[k] = image[0];
    mini[k] = image[0];
  }

  // if user has already set the maximum, then skip the checking
  if (!setmax) {
    // look for the maximum and minimum pixel value
    for (i=0; i<row; i++) {
      for (j=0; j<col; j++) {
        for (k=0; k<channel; k++) {
          if (image[(i*col+j)*channel+k] > maxi[k])
            maxi[k] = image[(i*col+j)*channel+k];
          if (image[(i*col+j)*channel+k] < mini[k])
            mini[k] = image[(i*col+j)*channel+k];
        }
      }
    }

    // if out-of-range, rescale
    // Modified: if not out-of-range, still rescale
    // Fixed equal max and min problem
    for (k=0; k<channel; k++) {
      maximum = 255;
      if (maxi[k] == 0.0 && mini[k] == 0.0) {
        for (i=0; i<row; i++)
          for (j=0; j<col; j++)
            image[(i*col+j)*channel+k] = 0.0;
      }
      else if (maxi[k] == mini[k]) {
        for (i=0; i<row; i++)
          for (j=0; j<col; j++)
            image[(i*col+j)*channel+k] = 255.0;
      }
      else {
        for (i=0; i<row; i++)
          for (j=0; j<col; j++)
            image[(i*col+j)*channel+k] = 
              (image[(i*col+j)*channel+k] - mini[k]) * 255.0 / 
              (maxi[k]-mini[k]); 
      }
    }
  }
  else {  // if the maximum has been set, any intensity that is
          // larger than the maximum would be set as maximum
    for (i=0; i<row; i++)
      for (j=0; j<col; j++)
        for (k=0; k<channel; k++)
          if (image[(i*col+j)*channel+k] > maximum)
            image[(i*col+j)*channel+k] = maximum;
  }
  
  ofp << maximum << endl;

  // convert the image data type back to unsigned char
  img = (unsigned char *) new unsigned char [row * col * channel];
  if (!img) 
    outofMemory();

  if (type == PGMRAW || type == PPMRAW) {
    for (i=0; i<row; i++)
      for (j=0; j<col; j++)
        for (k=0; k<channel; k++) 
          img[(i*col+j)*channel+k] = (unsigned char)image[(i*col+j)*channel+k];
    ofp.write((char *)img, (row * col * channel * sizeof(unsigned char)));
  }
  else {   // ASCII format
    for (i=0; i<row; i++)
      for (j=0; j<col; j++)
        for (k=0; k<channel; k++) 
          ofp << (unsigned int)image[(i*col+j)*channel+k] << endl; 
  }

  ofp.close();
  delete [] img;
}

// overloading () operator
float & Image::operator()(int i, int j, int k) const
{
  return image[(i*col+j)*channel+k];
}

// overloading = operator
const Image Image::operator=(const Image& img)
{
  int i, j, k;

  row = img.getRow();
  col = img.getCol();
  channel = img.getChannel();
  type = img.getType();
  createImage(row, col, type);             // allocate memory
  //  maximum = img.getMaximum();
  //  setmax = img.getSetmax();

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        image[(i*col+j)*channel+k] = img(i, j, k);

  return *this;
}

// overloading + operator
Image Image::operator+(const Image& img) const
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do addition\n";
    exit(3);
  }

  temp.createImage(row, col, type);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = image[(i*col+j)*channel+k] + img(i, j, k);

  return temp;
}

// overloading - operator
Image Image::operator-(const Image &img) const
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do subtraction\n";
    exit(3);
  }

  temp.createImage(row, col, type);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = image[(i*col+j)*channel+k] - img(i, j, k);

  return temp;
}

// overloading * operator
Image Image::operator*(const Image &img) const
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, "
         << "can't do multiplication\n";
    exit(3);
  }

  temp.createImage(row, col, type);
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = image[(i*col+j)*channel+k] * img(i, j, k);

  return temp;
}

// overloading / operator
Image Image::operator/(const Image &img) const
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do division\n";
    exit(3);
  }

  temp.createImage(row, col, type);
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = image[(i*col+j)*channel+k] / img(i, j, k);

  return temp;
}

// overloading ->* operator
Image Image::operator->*(const Image &img) const
{
  int i, j, k, m, nr, nc, nchan, nt;
  Image temp;
  float tmp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (col != nr || nt != type) {
    cout << "Image size is not consistent, can't do multiplication\n";
    exit(3);
  }

  temp.createImage(row, nc, type);
  
  for (k=0; k<channel; k++)
    for (i=0; i<row; i++)
      for (j=0; j<nc; j++) {
        tmp = 0.0;
        for (m=0; m<col; m++) {
          tmp += image[(i*col+m)*channel+k] * img(m, j, k);
          temp(i, j, k) = tmp;
        }
      }

  return temp;
}

// overloading the << operator
ostream & operator<<(ostream &out, Image &img)
{
  int i, j;
  
  for (i=0; i<img.getRow(); i++) {
    for (j=0; j<img.getCol(); j++)
      out << setw(4) << img(i,j) << ' ';
    out << endl;
  }

  return out; 
}

// logic operator AND
Image Image::AND(Image &img)
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do logic AND\n";
    exit(3);
  }

  temp.createImage(row, col, type);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = (int)image[(i*col+j)*channel+k] && (int)img(i, j, k);

  return temp;
}

// logic operator OR
Image Image::OR(Image &img)
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do logic OR\n";
    exit(3);
  }

  temp.createImage(row, col, type);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = (int)image[(i*col+j)*channel+k] || (int)img(i, j, k);

  return temp;
}

// logic operator NOT
Image Image::NOT()
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  temp.createImage(row, col, type);             

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = !(int)image[(i*col+j)*channel+k];

  return temp;
}

// logic operator XOR
Image Image::XOR(Image &img)
{
  int i, j, k, nr, nc, nchan, nt;
  Image temp;

  nr = img.getRow();
  nc = img.getCol();
  nchan = img.getChannel();
  nt = img.getType();

  if (nr != row || nc != col || nchan != channel || nt != type) {
    cout << "Images are not of the same size or type, can't do logic XOR\n";
    exit(3);
  }

  temp.createImage(row, col, type);             
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i, j, k) = XOR((int)image[(i*col+j)*channel+k], 
                  (int)img(i, j, k));

  return temp;
}    

// logic operator XOR
int Image::XOR(int a, int b)
{
  return (!a && b) || (a && !b);
}

// output out of memory error
void Image::outofMemory()
{
  cout << "Out of memory!\n";
  exit(1);
}

// return a clip of an image
Image Image::subImage(int startrow, int startcol, int endrow, int endcol)
{
  int i, j, k;
  Image temp;

  if (startrow > endrow || startcol > endcol || 
      endrow > row -1 || endcol > col - 1) {
    cout << "Check size\n";
    exit(3);
  }

  temp.createImage(endrow-startrow+1, endcol-startcol+1, type);

  for (i=startrow; i<=endrow; i++)
    for (j=startcol; j<=endcol; j++)
      for (k=0; k<channel; k++)
        temp(i-startrow, j-startcol, k) = image[(i*col+j)*channel+k];

  return temp;
}

// forward FFT
void Image::fft(Image &mag, Image &phase)
{
  fftifft(mag, phase, 1);
}

// inverse FFT
void Image::ifft(Image &mag, Image &phase)
{
  fftifft(mag, phase, -1);
}

// fast Fourier transform
void Image::fftifft(Image &mag, Image &phase, int scale)
{
  int N, M, L;
  int i, j, k, m, n, p, q, sign, index;
  float *realv, *imagv;
  float *real, *imag;            // the real and imaginary part of the fft
  float temp1, temp2, theta;
  int *loc;                      // the reordered index

  // allocate memory for real and imag
  real = (float *) new float [row * col];
  if (!real) {
    outofMemory();
  }
  imag = (float *) new float [row * col];
  if (!imag) { 
    outofMemory();
  }
  if (row != col) {
    cout << "The image has to be square and the dimension power of 2. " 
         << "The current image is not square.\n";
    exit(3);
  }
  
  if (pow(2.0, log((double)row)/log(2.0)) != (double)row) {
    cout << "The image dimension is not power of 2.\n";
    exit(3);
  }

  N = row; 
  M = N/2;
  L = (int)ceil(log((double)N)/log(2.0));
  realv = (float *) new float [N];
  if (!realv)
    outofMemory();
  imagv = (float *) new float [N];
  if (!imagv)
    outofMemory();
  loc = (int *) new int [N];
  if (!loc)
    outofMemory();

  // ordering the 1D vector for FFT
  loc[0] = 0;
  k = 1;
  m = 1;
  n = M;
  for (i=0; i<L; i++) {
    for (j=0; j<pow(2.0,(double)i); j++) {
      loc[k] = loc[k-m]+n;
      k++;
    }
    m <<= 1;
    n >>= 1;
  }

  // 1D FFT - 2D row transform
  for (i=0; i<row; i++) {
    for (j=0; j<col; j++) {
      sign = (int)pow(-1.0, (double)loc[j]);
      if (scale == 1) {
        realv[j] = sign * image[i*col+loc[j]];
        imagv[j] = sign * 0;
      }
      else {
        temp1 = mag(i, loc[j]);
        temp2 = tan(phase(i, loc[j]));
        realv[j] = sqrt(temp1 * temp1 / (1+temp2*temp2));
        if (phase(i, loc[j]) > PI/2 || phase(i, loc[j]) < -PI/2)
          realv[j] = -realv[j];
        imagv[j] = realv[j] * temp2;
      }
    }
    for (p=0; p<L; p++) {      // each level
      n = (int)pow(2.0,(double)p);
      for (q=0; q<n; q++) {
        theta = 2 * PI * q / (float)(2*n);
        for (k=q; k<N; k+=2*n) {
          temp1 = realv[k+n] * cos(theta) + imagv[k+n] * sin(theta) * scale;
          temp2 = -realv[k+n] * sin(theta) * scale + imagv[k+n] * cos(theta);
          realv[k+n] = realv[k] - temp1;
          realv[k] = realv[k] + temp1;
          imagv[k+n] = imagv[k] - temp2;
          imagv[k] = imagv[k] + temp2;
        }
      }
    }
    for (j=0; j<col; j++) {
      index = i*col + j;
      real[index] = realv[j];
      imag[index] = imagv[j];
    }    
  }

  // 1D FFT - 2D column transform
  for (j=0; j<col; j++) {
    for (i=0; i<row; i++) {
      sign = (int)pow(-1.0, (double)loc[i]);  // translate image to the center
      if (scale == 1) {
        realv[i] = sign * real[loc[i]*col+j];   // reorder
        imagv[i] = sign * imag[loc[i]*col+j];
      }
      else {
        realv[i] = real[loc[i]*col+j];
        imagv[i] = imag[loc[i]*col+j];
      }
    }
    for (p=0; p<L; p++) {    
      n = (int)pow(2.0,(double)p);
      for (q=0; q<n; q++) {
        theta = 2 * PI * q / (2*n);
        for (k=q; k<N; k+=2*n) {
          temp1 = realv[k+n] * cos(theta) + imagv[k+n] * sin(theta) * scale;
          temp2 = -realv[k+n] * sin(theta) * scale + imagv[k+n] * cos(theta);
          realv[k+n] = realv[k] - temp1;        // F(u+M)
          realv[k] = realv[k] + temp1;          // F(u)
          imagv[k+n] = imagv[k] - temp2;
          imagv[k] = imagv[k] + temp2;
        }
      }
    }
    for (i=0; i<row; i++) {
      index = i*col + j;
      real[index] = realv[i];
      imag[index] = imagv[i];
    }
  }

  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      index = i*col + j;
      if (scale == 1) {
        mag(i, j) = sqrt(real[index]*real[index] + imag[index]*imag[index]);
        phase(i, j) = atan2(imag[index], real[index]);
      }
      else {
        image[index] = sqrt(real[index]*real[index] + 
                imag[index]*imag[index])/(N*N);
      }
    }

  delete [] realv;
  delete [] imagv;
  delete [] real;
  delete [] imag;
  delete [] loc;
}

// log transformation
Image Image::logtran()
{
  Image temp;
  int i, j, k;

  temp.createImage(row, col, type);

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(i,j,k) = log(1.0+fabs(image[((i*col+j)*channel+k)]));

  return temp;
}

// convert the current color image from RGB to HSI model
Image Image::RGB2HSI()
{
  Image temp;        // a 3-channel image holds H, S, and I components
  float r, g, b, theta, mini;
  int i, j, index;

  if (channel != 3) {
    cout << "RGB2HSI: This is not a color image\n";
    exit(3);
  }

  temp.createImage(row, col, type);

  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      index = (i*col+j)*channel;
      r = image[index];
      g = image[index+1];
      b = image[index+2];
      mini = r;
      if (mini > g) 
        mini = g;
      if (mini > b)
        mini = b;
      theta = acos((0.5*((r-g)+(r-b))) / sqrt((r-g)*(r-g)+(r-b)*(g-b)+0.001));
      temp(i,j,0) = (g >= b) ? theta : (2*PI-theta);        // H component
      temp(i,j,2) = (r + g + b) / sqrt(3.0);                // I component
      if (temp(i,j,2)>0)
        temp(i,j,1) = 1 - sqrt(3.0) * mini / temp(i, j, 2);     // S component
      else
        temp(i,j,1) = 0;
    }

  return temp;
}

// convert the current color image from HSI to RGB model
Image Image::HSI2RGB()
{
  Image temp;        // a 3-channel image holds R, G, B components
  float hue, saturation, intensity, theta;
  int i, j, index;

  if (channel != 3) {
    cout << "HSI2RGB: This is not a color image\n";
    exit(3);
  }

  temp.createImage(row, col, type);

  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      index = (i*col+j)*channel;
      hue = image[index];
      saturation = image[index+1];
      intensity = image[index+2];
      if (hue < 2*PI/3) {
        temp(i,j,0) = intensity * 
          (1 + saturation * cos(hue)/cos(PI/3-hue)) / sqrt(3.0);
        temp(i,j,2) = intensity * (1-saturation) / sqrt(3.0);
        temp(i,j,1) = sqrt(3.0) * intensity - temp(i,j,0) - temp(i,j,2);
      }
      else if (hue < 4*PI/3) {
        temp(i,j,1) = intensity *
          (1 + saturation * cos(hue-2*PI/3)/cos(PI-hue)) / sqrt(3.0);
        temp(i,j,0) = intensity * (1 - saturation) / sqrt(3.0);
        temp(i,j,2) = sqrt(3.0) * intensity - temp(i,j,0) - temp(i,j,1);
      }
      else {
        temp(i,j,2) = intensity * 
          (1 + saturation * cos(hue-4*PI/3) / cos(5*PI/3 - hue)) / sqrt(3.0);
        temp(i,j,1) = intensity * (1 - saturation) / sqrt(3.0);
        temp(i,j,0) = sqrt(3.0) * intensity - temp(i,j,1) - temp(i,j,2);
      }
    }

  return temp;
}

// compute transpose of an image
Image Image::transpose()
{
  int i, j, k;
  Image temp;

  temp.createImage(col, row, type);

  for (i=0; i<row; i++)
    for (j=0; j<col; j++)
      for (k=0; k<channel; k++)
        temp(j, i, k) = image[(i*col+j)*channel+k];

  return temp;
}

// compute the inverse of an image using Gauss-Jordan elimination
Image Image::inverse()
{
  int i, j, m, k, pivot;
  Image temp, img;

  if (row != col) {
    cout << "Can't compute the inverse of a non-square matrix\n";
    exit(3);
  }

  temp.setRow(row);
  temp.setCol(2*col);
  temp.setType(type);
  temp.createImage();

  // construct the A | I matrix
  k = 0;
  for (i=0; i<row; i++) {
    for (j=0; j<col; j++)
      temp(i, j) = image[i*col+j];
    for (j=col; j<2*col; j++)
      temp(i, j) = 0.0;
    temp(i, col+k++) = 1;
  }

  pivot = 0;
  
  while (pivot < row) {
    m = findPivot(temp, pivot);    // partial pivoting

    if (m != pivot) {
      switchRow(temp, pivot, m);
    }
    else if (temp(pivot, pivot) == 0) {
      cout << "This is a singular matrix, noninvertible\n";
      exit(3);
    }

    if (temp(pivot, pivot) != 1.0) {
      dividePivot(temp, pivot);
    }
      
    eliminate(temp, pivot);

    pivot++;
  }

  img = temp.subImage(0, col, row-1, 2*col-1);

  return img;
}

// find the maximum absolute value in the pivotal column
int Image::findPivot(Image &img, int pivot)
{
  float maxi;
  int i, index;

  maxi = abs(img(pivot, pivot));
  index = pivot;

  for (i=pivot+1; i<row; i++)
    if (abs(img(i, pivot)) > maxi) {
      index = i;
      maxi = abs(img(i, pivot));
    }

  return index;
}

// switch two rows pivot <--> m
void Image::switchRow(Image &img, int pivot, int m)
{
  int i, j;
  float tmp;

  for (j=pivot; j<2*col; j++) {
    tmp = img(pivot, j);
    img(pivot, j) = img(m, j);
    img(m, j) = tmp;
  }
}

// divide the pivotal row by pivot
void Image::dividePivot(Image &img, int pivot)
{
  float scale;
  int j;

  scale = img(pivot, pivot);

  for (j=pivot; j<2*col; j++) 
    img(pivot, j) /= scale;
}

// eliminate elements on the pivotal column
void Image::eliminate(Image &img, int pivot)
{
  int i, j;
  float scale;

  for (i=pivot+1; i<row; i++) {
    scale = img(i, pivot);
    for (j=pivot; j<2*col; j++) 
      img(i, j) -= img(pivot, j) * scale;
  }

  for (i=pivot-1; i>=0; i--) {
    scale = img(i, pivot);
    for (j=pivot; j<2*col; j++) 
      img(i, j) -= img(pivot, j) * scale;
  }
}

// 1-D wavelet transform by Xiaoling Wang
// reference: Numerical Recipe in C
//   n - the length of the 1-D signal
//   isign - wavelet transform when isign=1
//           inverse wavelet transform when isign=-1
Image Image::wt(int n, int isign)
{
  int nn, i, j, ncor, count;
  float max;
  Image temp1, temp2,temp3;// temp2 rows temp3 for colunmns

  // check if the length of n is an integer power of 2
  // if not, then cut the rest of the image
  count = 0;
  while (n>=2) {
    if (n%2==0) {
      n /= 2;
      count++;
    }
    else  
      n--;
  }
 
  ncor = (int)pow(2.0,(double)count);
  //cout<<"ncor="<<ncor<<endl;
  if (ncor<4) 
    exit(1);

  temp1.createImage(1,ncor);
  temp2.createImage(row,ncor);
  temp3.createImage(ncor,row); 
  //cout<<"row="<<row<<endl;

  for (i=0; i<row; i++) {
    temp1 = subImage(i,0,i,ncor-1);

    if (isign>=0) {
      for (nn=ncor; nn>=ncor; nn>>=1){ 
         //cout<<"nn="<<nn<<endl;
         temp1.daub4(nn,isign);
        }  
    }

    else {
      //cout<<"ncor="<<ncor<<endl;
      for (nn=ncor; nn<=ncor; nn<<=1){
        //cout<<"nn2="<<nn<<endl;
        temp1.daub4(nn,isign);
      }
    }

    for (j=0; j<ncor; j++){
      temp2(i,j) = temp1(0,j);
      //cout<<"row="<<i<<"column="<<j<<temp2(i,j)<<endl;
    }
  }

  
  return temp2;
}

// Daubechies 4-coefficient wavelet filter
//   called by wt( )
//   n - length of the data vector, (power of 2)
//   isign - apply wavelet filter when isign=1
//           apply the transpose of wavelet filter when isign=-1
void Image::daub4(int n, int isign)
{
  Image wksp;
  int nh, nh1, i, j;

  if (n<4) 
    exit(1);

  wksp.createImage(1,n);
  nh1 = (nh=n>>1) + 1;
  //cout<<"nh="<<nh<<"  "<<"nh1="<< nh1 <<endl;
  if (isign>=0) {
    for (i=0,j=0; j<n-3; j+=2,i++) {
      wksp(0,i) = C0*image[j]+C1*image[j+1]+C2*image[j+2]+C3*image[j+3];
      //cout<<i<<"wksp(0,i)"<<wksp(0,i)<<endl;
      wksp(0,i+nh) = C3*image[j]-C2*image[j+1]+C1*image[j+2]-C0*image[j+3];
      //cout<<i+nh<<"wksp(0,i+nh)"<<wksp(0,i+nh)<<endl;
    }
    wksp(0,i) = C0*image[n-2]+C1*image[n-1]+C2*image[0]+C3*image[1];
    //cout<<i<<"wksp(0,i)"<<wksp(0,i)<<endl;
    wksp(0,i+nh) = C3*image[n-2]-C2*image[n-1]+C1*image[0]-C0*image[1];
    //cout<<i+nh<<"wksp(0,i)"<<wksp(0,i+nh)<<endl;
  }
  else {
    wksp(0,0) = C2*image[nh-1]+C1*image[n-1]+C0*image[0]+C3*image[nh1-1];
    wksp(0,1) = C3*image[nh-1]-C0*image[n-1]+C1*image[0]-C2*image[nh1-1];
    for (i=0,j=2; i<nh-1; i++) {
      //cout<<"i="<<i<<endl;  
      wksp(0,j++) = C2*image[i]+C1*image[i+nh]+C0*image[i+1]+C3*image[i+nh1];
      wksp(0,j++) = C3*image[i]-C0*image[i+nh]+C1*image[i+1]-C2*image[i+nh1];
    }
  }

  for (i=0; i<n; i++){ 
    image[i] = wksp(0,i);
  } 
}

// dilate image with structuring element
Image Image::dilate(Image &se, int origRow, int origCol)
{
  Image temp, seref;
  int nr, nc, nchan, nt, i, j, m, n, sum;

  nr = se.getRow();
  nc = se.getCol();
  nchan = se.getChannel();
  nt = se.getType();

  if (channel > 1 || nchan > 1) {
    cout << "This function can only dilate binary images\n";
    exit(3);
  }
  if (nt != type) {
    cout << "Can't dilate two images of different types\n";
    exit(3);
  }
  if (origRow >= nr || origCol >= nc || origRow < 0 || origCol < 0) {
    cout << "The origin of se needs to be inside se\n";
    exit(3);
  }
  if (!(int)se(origRow, origCol)) {
    cout << "The origin of se should be a foreground pixel\n";
    exit(3);
  }

  temp.createImage(row, col);
  
  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      if ((int)image[i*col+j]) {       
        // if the foreground pixel is not zero, then fill in the pixel
        // covered by the s.e.
        for (m=0; m<nr; m++)
          for (n=0; n<nc; n++) {
            if ((i-origRow+m) >= 0 && (j-origCol+n) >=0 && 
                (i-origRow+m) < row && (j-origCol+n) < col)
              if (!(int)temp(i-origRow+m, j-origCol+n))
                temp(i-origRow+m, j-origCol+n) = ( (int)se(m,n) ? 1 : 0 );
          }
      }
    }

  return temp;
}

// erode image with structuring element
Image Image::erode(Image &se, int origRow, int origCol)
{
  Image temp;
  int nr, nc, nchan, nt, i, j, m, n, sum, count;

  nr = se.getRow();
  nc = se.getCol();
  nchan = se.getChannel();
  nt = se.getType();

  if (channel > 1 || nchan > 1) {
    cout << "This function can only dilate binary images\n";
    exit(3);
  }
  if (nt != type) {
    cout << "Can't dilate two images of different types\n";
    exit(3);
  }
  if (origRow >= nr || origCol >= nc || origRow < 0 || origCol < 0) {
    cout << "The origin of se needs to be inside se\n";
    exit(3);
  }
  if (!(int)se(origRow, origCol)) {
    cout << "The origin of se should be a foreground pixel\n";
    exit(3);
  }

  temp.createImage(row, col);
  
  sum = 0;
  for (i=0; i<nr; i++)
    for (j=0; j<nc; j++)
      if ((int)se(i, j))
        sum++;

  for (i=0; i<row; i++)
    for (j=0; j<col; j++) {
      if ((int)image[i*col+j]) {     // if the foreground pixel is 1
        count = 0;
        for (m=0; m<nr; m++)
          for (n=0; n<nc; n++) {
            if ((i-origRow+m) >= 0 && (j-origCol+n) >=0 && 
                (i-origRow+m) < row && (j-origCol+n) < col)
              count += (int)image[(i-origRow+m)*col+(j-origCol+n)] && (int)se(m, n);
          }
        if (sum == count)
          temp(i, j) = 1;
      }
    }

  return temp;
}

// opening
Image Image::open(Image &se, int origRow, int origCol)
{
  Image temp1, temp2;
  int i, j;

  temp1 = erode(se, origRow, origCol);
  temp2 = temp1.dilate(se, origRow, origCol);

  return temp2;
}

// closing
Image Image::close(Image &se, int origRow, int origCol)
{
  Image temp1, temp2;
  int i, j;

  temp1 = dilate(se, origRow, origCol);
  temp2 = temp1.erode(se, origRow, origCol);

  return temp2;
}













