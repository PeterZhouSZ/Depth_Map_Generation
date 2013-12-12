#include "math.h"
#include "cv.h"
#include "highgui.h"
#include "FlyCapture2.h"
#include <iostream>
#include <fstream>
#include <string>
#include <direct.h>

using namespace std;
using namespace cv;
using namespace  FlyCapture2;
#define NUM_IMAGES 13
#define imageW 1920
#define imageH 1440

CvMat* mx1 = cvCreateMat( imageH, imageW, CV_32F );
CvMat* my1 = cvCreateMat( imageH, imageW, CV_32F );
CvMat* mx2 = cvCreateMat( imageH, imageW, CV_32F );
CvMat* my2 = cvCreateMat( imageH, imageW, CV_32F );

IplImage* capImg[NUM_IMAGES + 2] = {0};
//IplImage* capImgR[NUM_IMAGES] = {0};
IplImage* slPat[NUM_IMAGES/2 + 1] = {0};

IplImage* indexMatL = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);	
IplImage* indexMatR = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);	
IplImage* codeL = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_16U, 1);
IplImage* codeR = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_16U, 1);
IplImage* depthL = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_32F, 1);
IplImage* depthR = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_32F, 1);



void convertParams(){
	CvMat* leftCamera = cvCreateMat(3, 3, CV_64F);
	CvMat* rightCamera = cvCreateMat(3, 3, CV_64F);
	
	CvMat* D1 = cvCreateMat(1, 5, CV_64F);
	CvMat* D2 = cvCreateMat(1, 5, CV_64F);
	CvMat* R = cvCreateMat(3, 3, CV_64F);
	CvMat* T = cvCreateMat(3, 1, CV_64F);
	
	CvMat* R1 = cvCreateMat(3, 3, CV_64F);
	CvMat* R2 = cvCreateMat(3, 3, CV_64F);
	CvMat* P1 = cvCreateMat(3, 4, CV_64F);
	CvMat* P2 = cvCreateMat(3, 4, CV_64F);
	CvMat* Q = cvCreateMat(4, 4, CV_64F);

	CvSize imageSize = cvSize(imageW,imageH);
	

	leftCamera = (CvMat*)cvLoad("leftCamera.xml");
	rightCamera = (CvMat*)cvLoad("rightCamera.xml");

	D1 = (CvMat*)cvLoad("D1.xml");
	D2 = (CvMat*)cvLoad("D2.xml");
	
	cvSetReal2D(R,0, 0, 1.0000); cvSetReal2D(R,0, 1, 0.0009); cvSetReal2D(R,0, 2, 0.0044);  
	cvSetReal2D(R,1, 0, -0.0007); cvSetReal2D(R,1, 1, 0.9991); cvSetReal2D(R,1, 2,  -0.0419);  
	cvSetReal2D(R,2, 0,-0.0045); cvSetReal2D(R,2, 1, 0.0419); cvSetReal2D(R,2, 2, 0.9991);  
	cvSetReal2D(T,0, 0, -243.01258); cvSetReal2D(T, 1, 0, -1.40936); cvSetReal2D(T, 2, 0, -42.02111);  

	/*
	R	= (CvMat*)cvLoad("Rotation.xml");
	printf("%f %f %f\n", cvGetReal2D(R,0, 0), cvGetReal2D(R,0, 1),cvGetReal2D(R,0, 2));
	printf("%f %f %f\n", cvGetReal2D(R,1, 0), cvGetReal2D(R,1, 1),cvGetReal2D(R,1, 2));
	printf("%f %f %f\n", cvGetReal2D(R,2, 0), cvGetReal2D(R,2, 1),cvGetReal2D(R,2, 2));
	*/
	//T =	(CvMat*)cvLoad("Translation.xml");
	
	cvStereoRectify( leftCamera, rightCamera, D1, D2, cvSize(imageW,imageH),
					 R, T, R1, R2, P1, P2, Q, 0, -1, cvSize(0,0) );  /*CV_CALIB_ZERO_DISPARITY*/
	cvInitUndistortRectifyMap(leftCamera, D1, R1, P1, mx1, my1);
	cvInitUndistortRectifyMap(rightCamera, D2, R2, P2, mx2, my2);
	
	cvSave("map1x.xml", mx1);
	cvSave("map1y.xml", my1);

	cvSave("map2x.xml", mx2);
	cvSave("map2y.xml", my2);
	cvSave("reProj.xml", Q);
}

// 
int Dis(int a, int b, int grey){
	int counter = 0; 
	if (grey) {
		
		int sumAB = a^b;
		for (int i = 0; i < (NUM_IMAGES/2); i++){
			int base = pow(2.00,i);
		
			if (sumAB & base) {
				counter++;		
			}
		}

		if (counter > NUM_IMAGES/2){
			printf("Error! %d, %d \n", a, b);
		}
	}
	else
		counter = abs(a - b);
	return counter;
}	

void loadImage(CvSize size,int grey,int depthInd){
	char s[128];
	//IplImage * colorImg;
	for(int i = 0; i < NUM_IMAGES + 2; i++){
		capImg[i] = cvCreateImage(size,IPL_DEPTH_8U,1);
		if(i % 2) {
			if ( i == NUM_IMAGES + 1)
				sprintf(s,"../images/clean_imageR_%d.png",depthInd);  //clean_imageR_rect_%d
			else	
				sprintf(s,"../images/%d_resultR_%d.png",depthInd, i/2 + 1);   //right_rect_%d.jpg
		}
		else{
			if ( i == NUM_IMAGES)
				sprintf(s,"../images/clean_imageL_%d.png",depthInd);
			else
				sprintf(s,"../images/%d_resultL_%d.png",depthInd, i/2 + 1);
		}
		capImg[i] = cvLoadImage(s,0);
		//capImg[i] = cvCloneImage(colorImg);
		//cvCvtColor(colorImg, capImg[i], CV_RGB2GRAY);

		IplImage *tmp = cvCloneImage(capImg[i]);
		if(i % 2)
			cvRemap(capImg[i], tmp, mx2, my2, CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS);
		else
			cvRemap(capImg[i], tmp, mx1, my1, CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS);

		cvCopyImage(tmp,capImg[i]);
		//cvNamedWindow("debug");
		//cvShowImage("debug",capImg[i]);
		//cvWaitKey(-1);
		cvReleaseImage(&tmp);
		//cvDestroyWindow("debug");

	}
	//cvReleaseImage(&colorImg);
}

void initilizePattern(bool grey){
	char s[128];
	for(int i=0; i < NUM_IMAGES; i++)
		{
			if (grey) 
				//sprintf(s,"../images/%d.png", i+1);
				sprintf(s,"%d.png", i+1);
				
			else 
				sprintf(s,"../result/Color_Debruijin_1024_sub4.png");
				
			slPat[i]=cvLoadImage(s);
			if(!slPat[i]){
				printf("Could not load image file: %s\n",s);
				exit(0);
			}
	}

}

int RunSingleCamera( PGRGuid guid, char* filename)
{
    const int k_numImages = 1;

    Error error;
    Camera cam;

    // Connect to a camera
    error = cam.Connect(&guid);

    // Get the camera information
    CameraInfo camInfo;
    error = cam.GetCameraInfo(&camInfo);
      

    // Start capturing images
    error = cam.StartCapture();

    Image rawImage;    
    for ( int imageCnt=0; imageCnt < k_numImages; imageCnt++ )
    {                
        // Retrieve an image
        error = cam.RetrieveBuffer( &rawImage );
        printf( "Grabbed image %d\n", imageCnt );
        // Create a converted image
        Image convertedImage;

        // Convert the raw image
        error = rawImage.Convert( PIXEL_FORMAT_MONO8, &convertedImage );

        error = convertedImage.Save( filename );
    }            

    // Stop capturing images
    error = cam.StopCapture(); 

    // Disconnect the camera
    error = cam.Disconnect();
    return 0;
}

void captureImage(bool grey, int depthInd){
	char s[128]; 
	int index = 0;
	Error error;
	BusManager busMgr;

	//setting the projector
	printf("Starting capturing images:....\n");
	cvStartWindowThread();
	cvNamedWindow("mainWin", CV_WINDOW_NORMAL);
	cvMoveWindow("mainWin",3000,0);
	cvSetWindowProperty("mainWin", CV_WND_PROP_FULLSCREEN, CV_WINDOW_FULLSCREEN);
	cvWaitKey(500);  
	sprintf(s, "../images/depth/%d", depthInd);
	mkdir(s);
	//CreateDirectory((LPCTSTR) s, NULL);
	for (int i = 0; i < NUM_IMAGES; i = i + 1) {
	//for the stereo version
		cvShowImage("mainWin", slPat[i] ); //index
		cvWaitKey(100);
		for (int cap = 0; cap < 2; cap++) {	
			PGRGuid guid;
			error = busMgr. GetCameraFromIndex(cap, &guid);
			if (cap) {
				sprintf(s,"../images/depth/%d/R_%02d.png", depthInd, i + 1);
			}
			else{
				sprintf(s,"../images/depth/%d/L_%02d.png", depthInd, i + 1);
			}
			RunSingleCamera(guid, s);
		}
		cvWaitKey(3);
	}
	//cvWaitKey(-1);
	cvDestroyWindow("mainWin");
	//loadImage(cvSize(imageW,imageH), grey, depthInd);
}
/*
void captureImage(bool grey,int depthInd){ 
	 char s[128];

	 CvCapture* captureL = cvCaptureFromCAM(0);
	 CvCapture* captureR = cvCaptureFromCAM(1);

	 cvSetCaptureProperty(captureL, CV_CAP_PROP_FRAME_WIDTH, imageW);
	 cvSetCaptureProperty(captureL, CV_CAP_PROP_FRAME_HEIGHT, imageH);

	 cvSetCaptureProperty(captureR, CV_CAP_PROP_FRAME_WIDTH, imageW);
	 cvSetCaptureProperty(captureR, CV_CAP_PROP_FRAME_HEIGHT, imageH);

     printf("Initializing images.... \n");

	 for (int i = 0; i < NUM_IMAGES + 1; i = i + 2){
		 if (grey) {
			capImg[i] = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);
			capImg[i + 1] = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);
		 }
		 else{
			capImg[i] = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 3);
			capImg[i + 1] = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 3);
		 }
	 }


	 int index = 0;
	 for(int i = 0; i < NUM_IMAGES + 1; i = i + 2, index ++)
	 {
		printf("Capturing image number: %d\n", i);
		IplImage * tL, *tR;
		IplImage *tmpImgL, *tmpImgR;
		IplImage *tmpL, *tmpR;
		// show the image
		cvShowImage("mainWin", slPat[index] );
		cvWaitKey(100);

		for (int q = 0; q < 10; ++q) {
		// capture the image
		cvGrabFrame(captureL);
		cvGrabFrame(captureL);
		cvGrabFrame(captureL);
		cvGrabFrame(captureL);
		cvGrabFrame(captureL);

		cvGrabFrame(captureR); 
		cvGrabFrame(captureR);
		cvGrabFrame(captureR);
		cvGrabFrame(captureR);
		cvGrabFrame(captureR);
		cvWaitKey(100);

		tL = cvRetrieveFrame(captureL);   
		tR = cvRetrieveFrame(captureR);


		if (q == 0){
			tmpImgL = cvCloneImage(tL);
			tmpImgR = cvCloneImage(tR);
		} else {
			tmpL = cvCloneImage(tL);
			tmpR = cvCloneImage(tR);
			cvAddWeighted(tmpL, 1.00/(q+1.00), tmpImgL, q/(q+1.00), 0.00, tmpImgL);
			cvAddWeighted(tmpR, 1.00/(q+1.00), tmpImgR, q/(q+1.00), 0.00, tmpImgR);
			cvShowImage("show", tmpL);
			cvWaitKey(15);
			cvReleaseImage(&tmpL);
			cvReleaseImage(&tmpR);
			}
		}
		if (grey) {
			IplImage *tmp = cvCloneImage(capImg[0]);
			cvCvtColor(tmpImgL, capImg[i], CV_RGB2GRAY);
			cvCvtColor(tmpImgR, capImg[i + 1], CV_RGB2GRAY);
			
			cvRemap(capImg[i], tmp, mx1, my1, CV_INTER_CUBIC+CV_WARP_FILL_OUTLIERS);
			cvCopyImage(tmp,capImg[i]);
			
			cvRemap(capImg[i + 1], tmp, mx2, my2, CV_INTER_CUBIC+CV_WARP_FILL_OUTLIERS);
			cvCopyImage(tmp,capImg[i + 1]);
			cvReleaseImage(&tmp);
		}
		else{
			cvCopyImage(tmpImgL,capImg[i]);
			cvCopyImage(tmpImgR,capImg[i + 1]);
		}
		// write the image
		if (i == NUM_IMAGES) 
			sprintf(s,"../images/clean_imageL_%d.png",depthInd);
		else
			sprintf(s,"../images/resultL_%d.png",i/2 + 1);
		cvSaveImage(s,tmpImgL);

		if (i == NUM_IMAGES + 1) 
			sprintf(s,"../images/clean_imageR_%d.png",depthInd);
		else
			sprintf(s,"../images/resultR_%d.png",i/2 + 1);
		cvSaveImage(s,tmpImgR);	

		cvReleaseImage(&tmpImgL);
		cvReleaseImage(&tmpImgR);
		cvWaitKey(30);

		}
		cvReleaseCapture( &captureL );
		cvReleaseCapture( &captureR );
		printf("Done capture...\n");
		cvDestroyWindow("mainWin");
}
*/
void codeImageGreyNew(int strPos, int endPos, bool svCode, bool visCode, int codeThre, int visThre){
	char s[128];
	IplImage *tmpW, *tmpB, *unCode, *unCodePre, *temp_double, *code;
	IplImage* invCode;
	int index = 1;
	printf("Starting code calculation...\n");
	
	//Initialize the code
	tmpW = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_8U,1);
	tmpB = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_8U,1);
	unCode = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_8U,1);
	unCodePre = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_8U,1);
	code = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_16U,1);
	temp_double = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_16U, 1);
	invCode = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);
	cvSetZero(unCode); cvSetZero(unCodePre); cvSetZero(code);
	cvAddWeighted(unCodePre,1, unCode, 0, 1, unCodePre);
	
	for (int i = strPos; i< endPos; i = i + 2){
		cvSub(capImg[endPos + 2],capImg[i],invCode);
		cvAddS(invCode,cvRealScalar(codeThre),tmpW);
		cvSubS(invCode,cvRealScalar(codeThre),tmpB);
		cvCmp(capImg[i], tmpW, tmpW, CV_CMP_GT);
		cvCmp(capImg[i], tmpB, tmpB, CV_CMP_LT);
			
		cvOr(tmpW,tmpB,unCode);
		cvAnd(unCode, unCodePre, unCodePre);
		cvAnd(tmpW,unCodePre,tmpW);
		cvAnd(tmpB,unCodePre,tmpB);

		cvConvert(tmpW, temp_double);
		cvAddWeighted(code, 1.00, temp_double, (pow(2.00, (NUM_IMAGES/2 - index)/255.00 )), 0, code); 
		index = index + 1;
	}
	cvMul(code, unCodePre, code);
	printf("Done...\n");
	cvReleaseImage(&tmpW);
	cvReleaseImage(&tmpB);
	cvReleaseImage(&unCode);
	cvReleaseImage(&unCodePre);
	cvReleaseImage(&temp_double);
	//Save the code
	if (svCode){
		FILE *f;
		if (strPos%2 )
			f = fopen("../result/codeR.xml","wt");
		else
			f = fopen("../result/codeL.xml","wt");
	
		for (int i = 0;i<code->height;i++){
			for (int j = 0;j<code->width;j++){
				fprintf(f,"%f\n",cvGetReal2D(code,i,j));
				}
		}
		fclose(f);
	}

	// Visualize the code map
	if (visCode){
		IplImage *temp = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);
		cvNamedWindow("code");
		cvNormalize( code, temp, 0, 256, CV_MINMAX );
		cvShowImage("code", temp);
		cvWaitKey(-1);
		cvDestroyWindow("code");
		cvReleaseImage(&temp);
	}

	// Copy the result to corresponding matrix
	if (strPos%2 == 0){
		sprintf(s,"../result/codeL.png");
		cvCopyImage(code,codeL);
	}
	else{
		sprintf(s,"../result/codeR.png");
		cvCopyImage(code,codeR);
	}
	cvReleaseImage(&code);
}

void imagesc(IplImage *img, char * str, float low, float high){
	
	cvAbs(img, img);
	cvAddWeighted(img, 1.00/(high-low), img, 0.00, -low/(high-low), img);
	cvShowImage(str, img);
	cvSaveImage("../result/depth.png",img);

}

void calDepth(IplImage *codeMapL, IplImage* codeMapR, int lr, int grey, int depthInd){
	
	CvSize size = cvSize(imageW, imageH);
	IplImage* depthMap = cvCreateImage(size, IPL_DEPTH_32F, 1);
	IplImage* indexMat;
	printf("Calculating Depth Map...");
	cvNamedWindow("depth");
	FILE *depth;
	int low, high;
	char s[128];
	if (lr == 0){
		sprintf(s,"../result/resultL_%d.xml",depthInd);
		low = 100; high = 500;  //150, 50
	}
	else {
		sprintf(s,"../result/resultR_%d.xml",depthInd);
		low = 125; high = 200; //50, 150
	}
	depth = fopen(s,"wt");

	for (int y = 0; y < codeMapL->height; y++){
		for (int xL = 0; xL < codeMapL->width; xL++){
			if (cvGetReal2D(codeL,y,xL) != 0){				// &&  cvGetReal2D(indexMat,y,xL)/255 == 1
				int intensityL = cvGetReal2D(codeMapL,y,xL);
				int intensityR;
				int diffPre, diffCur,posDiffPre = 0; 
				diffPre = 256;
				for(int xR = max(xL - high, 0); xR< max(xL - low, 0);xR++){ //
					if( cvGetReal2D(codeR,y,xR) != 0) {
						intensityR = cvGetReal2D(codeMapR,y,xR);
						diffCur = Dis(intensityR,  intensityL, grey); //grey
						if(diffPre > diffCur){
							diffPre = diffCur;
							posDiffPre = (xR - xL);
						}
					}
				}
				cvSetReal2D(depthMap,y,xL,posDiffPre);
				fprintf(depth,"%d\n",posDiffPre);
			}
			else{
				cvSetReal2D(depthMap,y,xL,0);
				fprintf(depth,"%d\n",0);
			}

		}

		//printf("Line %d\n",y);
	}
	fclose(depth);
	

	if (lr == 0)
		cvCopyImage(depthMap,depthL);
	else
		cvCopyImage(depthMap,depthR);

	cvAbs(depthMap, depthMap);
	
	cvErode(depthMap, depthMap, NULL, 2);
	cvDilate(depthMap, depthMap, NULL, 2);
	
	//	cvSmooth(depthMap, depthMap, CV_MEDIAN, 5,5);
	imagesc(depthMap,"depth",40, 80);
	cvWaitKey(-1);

	CvMat* Q = cvCreateMat(4,4,CV_32F);
	IplImage *xC = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_32F,1);
	IplImage *yC = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_32F,1);
	IplImage *zC = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_32F,1);

	Q = (CvMat*)cvLoad("reProj.xml");
	CvMat* pCloud = cvCreateMat(imageH,imageW,CV_32FC3);

	cvReprojectImageTo3D(depthMap, pCloud,Q, 1);	
	cvSplit(pCloud,xC, yC, zC, NULL);
	if (lr == 0) 
		sprintf(s,"../result/pCloudL_%d.txt",depthInd);
	
	else
		sprintf(s,"../result/pCloudR_%d.txt",depthInd);
	
	//cvSave(s, pCloud);
	FILE *f = fopen(s,"wt");
	for(int i = 0; i<imageH; i++) {
		for(int j = 0; j<imageW; j++){
			fprintf(f,"%f\n",cvGetReal2D(xC,i,j));
			fprintf(f,"%f\n",cvGetReal2D(yC,i,j));
			fprintf(f,"%f\n",cvGetReal2D(zC,i,j));
		}
	}
	fclose(f);

	cvDestroyWindow("depth");
	cvReleaseImage(&depthMap);
	cvReleaseImage(&indexMat);
	cvReleaseImage(&xC);
	cvReleaseImage(&yC);
	cvReleaseImage(&zC);
	cvReleaseMat(&Q);
	cvReleaseMat(&pCloud);
}


void sl(bool option,int index){
	// Structured Light
	if(option){
		loadImage(cvSize(imageW,imageH),1,index);
	}
	else{
		initilizePattern(1);
		captureImage(1,index);
		}
 
	//int thresholdL = 40;
	//int thresholdR = 50;
	//codeImageGreyNew(0, NUM_IMAGES - 2, 1, 1, 3.00, thresholdL);
	//codeImageGreyNew(1, NUM_IMAGES - 1, 1, 1, 3.00, thresholdR);
	
	
	cvReleaseImage(capImg);
	cvReleaseImage(slPat);
	//calDepth(codeL, codeR, 0, 1, index);
	//calDepth(codeR, codeL, 1, 1, index);
	printf("Done..\n");
}	


void rmShadow(float threshold){
	int  disL, disR, diff;
	FILE *f = fopen("../result/testDepth.xml","wt");
	FILE *f1 = fopen("../result/testDiff.xml","wt");
	IplImage* depthMap = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_32F, 1);
	IplImage* diffMap = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_32F, 1);
	for(int i = 0; i < imageH;i++){
		for(int j = 0; j < imageW; j++){
			disL = cvGetReal2D(depthL,i,j);
			disR = cvGetReal2D(depthR,i,j+disL);
			diff = abs(disL) - (disR);
			fprintf(f1,"%f\n",diff);
			cvSetReal2D(diffMap,i,j,diff);
			if (abs(diff)  <= threshold){
				cvSetReal2D(depthMap,i,j,abs(disL));
				fprintf(f,"%f\n",abs(disL));
			}
			else{
				cvSetReal2D(depthMap,i,j,0);
				fprintf(f,"%f\n",0);
			}
		}

	}
	fclose(f);
	fclose(f1);
	cvNamedWindow("depth");
	//cvNormalize(diffMap,diffMap,0, 256,CV_MINMAX);
	//cvShowImage("depth",diffMap);
	//cvWaitKey(-1);
	imagesc(depthMap,"depth",10,40);
	cvWaitKey(-1);

}

int main(int numImages){
	printf("Start structured light scanning....\n");
	//convertParams();
	
	//sl(0,13);
	int index = 1;
	char userInput;
	while (index <= numImages) {
		cin >> userInput;
		
		if (userInput == 'c'){
			sl(0,index);
			index++;
		}
		
	}
	
	
	// for removing test
	/*
	float temp;
	while (1) {
	printf("Input the threshold you want to use to remove the shadow (0  for exit): \n");
	cin >> temp;
	if (temp == 0)
		break;
	rmShadow(temp);
	}
	*/
	cvReleaseImage(&codeL);
	cvReleaseImage(&codeR);
	cvReleaseImage(&depthL);
	cvReleaseImage(&depthR);
	cvReleaseImage(&indexMatL);
	cvReleaseImage(&indexMatR);
	cvReleaseMat(&mx1);
	cvReleaseMat(&mx2);
	cvReleaseMat(&my1);
	cvReleaseMat(&my2);
	
}


/*
void codeImageGrey(bool lr,bool svCode, bool visCode, int threshold){
	char s[128];
	CvSize size = cvSize(imageW,imageH);
	IplImage * code;
	if (lr == 0)
		code = cvCloneImage(codeL);
	else
		code = cvCloneImage(codeR);

	cvNamedWindow("debug");
	IplImage * tmp,*minImage,*maxImage,*half;
	minImage = cvCreateImage(size,capImgL[0]->depth,capImgL[0]->nChannels);
	maxImage = cvCreateImage(size,capImgL[0]->depth,capImgL[0]->nChannels);
	cvSetZero(code);
	printf("Computing min/max for each pixel....\n");
	if (lr == 0){
		tmp = cvCloneImage(capImgL[0]);
		minImage = cvCloneImage(capImgL[0]);
		maxImage = cvCloneImage(capImgL[0]);
		for (int i =0; i< NUM_IMAGES; i++){
			//cvRemap(capImgL[i], tmp, mx1, my1, CV_INTER_CUBIC+CV_WARP_FILL_OUTLIERS);
			//cvCopyImage(tmp, capImgL[i]);
			cvMin(capImgL[i], minImage, minImage);
			cvMax(capImgL[i], maxImage, maxImage);
		}

	}
	else {
		tmp = cvCloneImage(capImgR[0]);
		minImage = cvCloneImage(capImgR[0]);
		maxImage = cvCloneImage(capImgR[0]);
		for (int i =0; i< NUM_IMAGES; i++){
			//cvRemap(capImgR[i], tmp, mx2, my2, CV_INTER_CUBIC+CV_WARP_FILL_OUTLIERS);
			//cvCopyImage(tmp, capImgR[i]);
			cvMin(capImgR[i], minImage, minImage);
			cvMax(capImgR[i], maxImage, maxImage);
		}
	}
	printf("Done with code calculation....\n");

	half = cvCreateImage(size, minImage->depth, minImage->nChannels);
	cvAddWeighted(maxImage, 0.5, minImage, 0.5, 0, half);  
	
	IplImage* temp = cvCreateImage(size, IPL_DEPTH_8U, 1);
	IplImage* temp_double = cvCreateImage(size, code->depth, 1);
	

	printf("Starting code calculation...\n");
	for (int i = 0; i < NUM_IMAGES; i++){	 
		if (lr == 0)
			cvCmp(capImgL[i], half, temp, CV_CMP_GT);
		else
			cvCmp(capImgR[i], half, temp, CV_CMP_GT);
		cvConvert(temp, temp_double);
		cvAddWeighted(code, 1.00, temp_double, (pow(2.00, NUM_IMAGES - i - 1)/255.00), 0, code); //pow(2.00, -i)/255.00
		cvShowImage("debug", temp_double);
		cvWaitKey(30);
	}
	cvDestroyWindow("debug");

	printf("Done with code calculation....\n");
	// For saving the code and debuging
	if (svCode){
		FILE *f;
		if (lr == 0)
			f = fopen("../result/codeL.xml","wt");
		else
			f = fopen("../result/codeR.xml","wt");
	
		for (int i = 0;i<code->height;i++){
			for (int j = 0;j<code->width;j++){
				fprintf(f,"%i\n",cvGetReal2D(code,i,j));
				}
		}
		fclose(f);
	}

	// Visualize the code map
	if (visCode){
		cvNamedWindow("code");
		cvNormalize( code, temp, 0, 256, CV_MINMAX );
		cvShowImage("code", temp);
		cvWaitKey(-1);
		cvDestroyWindow("code");
	}

	//calculate index matrix
	IplImage* ind = cvCreateImage(cvSize(imageW,imageH), IPL_DEPTH_8U, 1);
	cvSetZero(ind);
	IplImage* subImage = cvCreateImage(size, minImage->depth, minImage->nChannels);
	cvAddWeighted(maxImage,1, minImage,-1, 0, subImage);

		for (int i = 0; i<size.height; i++)
			for(int j = 0; j<size.width; j++)
				if (cvGetReal2D(subImage,i,j) > threshold)
					cvSetReal2D(ind,i,j, 255);
	if (lr == 0) 	
		cvCopyImage(ind,indexMatL);
	else
		cvCopyImage(ind,indexMatR);

	cvReleaseImage(&subImage);
	cvReleaseImage(&ind);
	
	//Saving the results
	if (lr == 0){
		sprintf(s,"../result/codeL.png");
		cvCopyImage(code,codeL);
	}
	else{
		sprintf(s,"../result/codeR.png");
		cvCopyImage(code,codeR);
	}
	cvSaveImage(s, code);
	

	cvReleaseImage(&tmp);
	cvReleaseImage(&temp);
	cvReleaseImage(&temp_double);
	cvReleaseImage(&minImage);
	cvReleaseImage(&maxImage);
	cvReleaseImage(&half);
	cvReleaseImage(&code);
	
}

*/