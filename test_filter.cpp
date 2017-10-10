
#include "Image.h"  
#include "IpFunctions.h"       
#include <iostream>
#include <cstdlib>
using namespace std;


#define Usage "./test  input-img output-img \n"

int main(int argc, char **argv)
{
  Image img1, img2, filter,img3,filter2,img4,img5,img6,img7,img8,img9,img10,img11,img12;
	int rds=1;
	int nc,nr,ntype,nchan;
	img1.readImage(argv[1]);
	nc = img1.getCol();
	nr = img1.getRow();
	ntype = img1.getType();
	nchan = img1.getChannel();
  // objects are initialized with 0 rows, 0 cols, and type 1. there is no data.

  if (argc < 3) {
    cout << Usage;
    exit(3);
  }

  // read in the input image
	
	
	
	filter.createImage(3, 3, ntype);
	filter (-1,-1) = -1; //laplacian filter values
	filter (-1, 0) = 0;
	filter (-1, 1) = 1;
	filter ( 0,-1) = -2;
	filter ( 0, 0) =0;
	filter ( 0, 1) = 2;
	filter ( 1,-1) = -1;
	filter ( 1, 0) = 0;
	filter ( 1, 1) = 1;
	
	int filter_size= filter.getCol();

	filter2.createImage(3, 3, ntype);
	filter2 (-1,-1) = 1; //laplacian filter values
	filter2 (-1, 0) = 2;
	filter2 (-1, 1) = 1;
	filter2 ( 0,-1) = 0;
	filter2 ( 0, 0) =0;
	filter2 ( 0, 1) = 0;
	filter2 ( 1,-1) = -1;
	filter2 ( 1, 0) = -2;
	filter2 ( 1, 1) = -1;
	
	int filter2_size= filter2.getCol();
	
	
	
	img2=filterImage(img1, filter)+filterImage(img1, filter2);
	img6=cannyImage(img2);//gradient magnitude part
	img3=filterImage(img1,filter2);
	img4=filterImage(img1,filter);
	img8=canny2Image(img3, img4);
	img7=canny3Image(img6,img8);
	img5=canny4Image(img7,150);
	img9=canny4Image(img7,80);
	img10=img9+img5;
	img11=canny5Image(img10);
	//img12=img11+img10;



	//img3=filterImage(img1, filter);
	//img3=cannyImage(img3);

	//img4=filterImage(img1, filter2);
	//img4=cannyImage(img4);
	//img5=img4/img3;
	//img7=canny2Image(img5);//theta part gradient
	//img8=canny4Image(img6);
	
		

	
  // write image2 to the output image
  img11.writeImage(argv[2]);

  return 0;
}

