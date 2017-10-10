
#include "Image.h"
#include "IpFunctions.h" 
#include <iostream>
#include <cmath>

using namespace std;
using std::cout;
		
Image filterImage(Image &inimg, Image &filter )
{
	Image outimg;	
	int nr, nc, nchan, ntype;
	int i, j, k,a;
	int rds=1;
	int filter_size= filter.getCol();
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 

	
	for (i=1; i<nr; i++)
		for (j=1; j<nc; j++)
		{
			for (k=0; k<nchan; k++)
			a=inimg(i-1,j-1,k)*filter(-1,-1)+inimg(i,j-1,k)*filter(0,-1)+inimg(i+1,j-1,k)*filter(1,-1)+inimg(i-1,j,k)*filter(-1,0)+inimg(i,j,k)*filter(0,0)+inimg(i+1,j,k)*filter(1,0)+inimg(i-1,j+1,k)*filter(-1,1)+inimg(i,j+1,k)*filter(0,1)+inimg(i+1,j+1,k)*filter(1,1);
			outimg(i,j,k)=a*a;
      	}

	return outimg;
}

Image cannyImage(Image &inimg)
{
	Image outimg;	
	int nr, nc, nchan, ntype;
	int i, j, k;
	
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 

	
	for (i=1; i<nr; i++)
		for (j=1; j<nc; j++)
		{
			for (k=0; k<nchan; k++)
			outimg(i,j,k)=sqrt(inimg(i,j,k));
      	}

	return outimg;
}
Image canny2Image(Image &inimg , Image &inimg2)
{
	Image outimg;	
	int nr, nc, nchan, ntype;
	int i, j, k;
	
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 

	
	for (i=1; i<nr; i++)
		for (j=1; j<nc; j++)
		{
			for (k=0; k<nchan; k++)
			outimg(i,j,k)=atan2(inimg2(i,j,k),inimg(i,j,k))*180/3.14;
			
      	}

	return outimg;
}
Image canny3Image(Image &inimg,Image &inimg2)
{
	Image outimg;
	Image img,img2;
			
	int nr, nc, nchan, ntype;
	int i, j, k;
	int f,g;
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  	

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 
	
	
	for (i=1; i<nr; i++){
		for (j=1; j<nc; j++)
		{
			for (k=0; k<nchan; k++){
				
				if((inimg2(i,j,k)>=0) && (inimg2(i,j,k)<45)){ 
					if((inimg(i,j,k)>=inimg(i,j-1,k)) && (inimg(i,j,k)>=inimg(i,j+1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
						outimg(i,j,k)=0;
					
									
				}
				if((inimg2(i,j,k)>=45) && (inimg2(i,j,k)<90)){
					if((inimg(i,j,k)>=inimg(i-1,j-1,k)) && (inimg(i,j,k)>=inimg(i+1,j+1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
									
				}
				if((inimg2(i,j,k)>=90) && (inimg2(i,j,k)<135)){

					if((inimg(i,j,k)>=inimg(i-1,j,k)) && (inimg(i,j,k)>=inimg(i+1,j,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}
				if((inimg2(i,j,k)>=135) && (inimg2(i,j,k)<180)){

					if((inimg(i,j,k)>=inimg(i-1,j-1,k)) && (inimg(i,j,k)>=inimg(i+1,j+1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}
				if((inimg2(i,j,k)>= 180) && (inimg2(i,j,k)<225)){

					if((inimg(i,j,k)>=inimg(i,j-1,k)) && (inimg(i,j,k)>=inimg(i,j+1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}
				if((inimg2(i,j,k)>= 225) && (inimg2(i,j,k)<270)){

					if((inimg(i,j,k)>=inimg(i+1,j-1,k)) && (inimg(i,j,k)>=inimg(i-1,j+1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}
				if((inimg2(i,j,k)>= 270) && (inimg2(i,j,k)<315)){

					if((inimg(i,j,k)>=inimg(i+1,j,k)) && (inimg(i,j,k)>=inimg(i-1,j,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}
				if((inimg2(i,j,k)>= 315) && (inimg2(i,j,k)<=360)){

					if((inimg(i,j,k)>=inimg(i+1,j+1,k)) && (inimg(i,j,k)>=inimg(i-1,j-1,k))){
						outimg(i,j,k)=inimg(i,j,k);
					}
					else
					outimg(i,j,k)=0;
				}





			}
		}				
      	}

	return outimg;
}
Image canny4Image(Image &inimg,int T)
{
	Image outimg;	
	int nr, nc, nchan, ntype;
	int i, j, k;
	Image max;
	
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 

	
	for (i=1; i<nr; i++){
		for (j=1; j<nc; j++){
			for (k=0; k<nchan; k++){

				if(inimg(i,j,k)<T){
					outimg(i,j,k)=0;}
				else
					outimg(i,j,k)=255;
				
			}	
		}
						
      	}

	return outimg;
}

Image canny5Image(Image &inimg)
{
	Image outimg;	
	int nr, nc, nchan, ntype;
	int i, j, k;
	Image max;
	
	nc = inimg.getCol();
	nr = inimg.getRow();
	ntype = inimg.getType();
	nchan = inimg.getChannel();
  

  	
  	outimg.createImage(nr, nc, ntype);
	
	outimg.initImage(0); 

	
	for (i=1; i<nr; i++){
		for (j=1; j<nc; j++){
			for (k=0; k<nchan; k++){

				if((inimg(i,j,k)>=inimg(i-1,j-1,k))&&(inimg(i,j,k)>=inimg(i+1,j+1,k))&&(inimg(i,j,k)>=inimg(i-1,j,k))&&(inimg(i,j,k)>=inimg(i-1,j+1,k))&&(inimg(i,j,k)>=inimg(i,j-1,k))&&(inimg(i,j,k)>=inimg(i,j+1,k))&&(inimg(i,j,k)>=inimg(i+1,j-1,k))&&(inimg(i,j,k)>=inimg(i+1,j,k))){
					outimg(i,j,k)=0;}
				else
					outimg(i,j,k)=255;
				
			}	
		}
						
      	}

	return outimg;
}

