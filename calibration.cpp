#include "cv.h"
#include "cxmisc.h"
#include "highgui.h"
#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "FlyCapture2.h"

#define imageW  1920
#define imageH 1440

using namespace std;
using namespace cv;
using namespace  FlyCapture2;
//CvMat *mapL = cvCreateMat( 480, 640, CV_8U);
//CvMat *mapR = cvCreateMat( 480, 640, CV_8U);

//
// Given a list of chessboard images, the number of corners (nx, ny)
// on the chessboards, and a flag: useCalibrated for calibrated (0) or
// uncalibrated (1: use cvStereoCalibrate(), 2: compute fundamental
// matrix separately) stereo. Calibrate the cameras and display the
// rectified results along with the computed disparity images.
//
void takePhoto(int);

//void sl(void);
//int gray2bin(int);

static void
	StereoCalib(int nx, int ny, int useUncalibrated, bool photoTake)
{	
	//#define imageW 640
	//#define imageH 480
	int displayCorners = 1;
	int showUndistorted = 1;
	bool isVerticalStereo = false;//OpenCV can handle left-right
	//or up-down camera arrangements
	const int maxScale = 1;
	const float squareSize = 0.035f; //Set this to your actual square size
	const int numImages = 40;

	printf("Square Size: %.3f\n", squareSize);

	if (photoTake) 
		takePhoto(numImages);


	//FILE* f = fopen(imageList, "rt");

	int j, lr, nframes, n = nx*ny, N = 0;
	vector<string> imageNames[2];
	vector<CvPoint3D32f> objectPoints;
	vector<CvPoint2D32f> pointsOld[2];
	vector<int> npoints;
	vector<int> active[2];
	vector<CvPoint2D32f> temp(n);
	CvSize imageSize = {0,0};
	// ARRAY AND VECTOR STORAGE:
	double M1[3][3], M2[3][3], D1[5], D2[5];
	double R[3][3], T[3], E[3][3], F[3][3];
	CvMat _M1 = cvMat(3, 3, CV_64F, M1 );
	CvMat _M2 = cvMat(3, 3, CV_64F, M2 );
	CvMat _D1 = cvMat(1, 5, CV_64F, D1 );
	CvMat _D2 = cvMat(1, 5, CV_64F, D2 );

	CvMat _R = cvMat(3, 3, CV_64F, R );
	CvMat _T = cvMat(3, 1, CV_64F, T );
	CvMat _E = cvMat(3, 3, CV_64F, E );
	CvMat _F = cvMat(3, 3, CV_64F, F );
	if( displayCorners )
		cvNamedWindow( "corners", 1 );

	// READ IN THE LIST OF CHESSBOARDS:
	/* if( !f )
	{
	fprintf(stderr, "can not open file %s\n", imageList );
	return;
	}
	*/
	printf("Starting main loop\n");
	for(int i=0; i< 2*numImages ;i++)
	{
		char str[128];
		printf("Iteration %d\n", i);
		//char buf[1024];
		int count = 0, result=0;
		lr = i % 2;
		vector<CvPoint2D32f>& pts = pointsOld[lr];
		//if( !fgets( buf, sizeof(buf)-3, f ))
		//   break;
		//size_t len = strlen(buf);
		//while( len > 0 && isspace(buf[len-1]))
		//   buf[--len] = '\0';
		//if( buf[0] == '#')
		//  continue;
		if (lr)
			sprintf(str,"../images/stereo/right_%d.png",i/2 + 1);
		else
			sprintf(str,"../images/stereo/left_%d.png",i/2 + 1);

		IplImage* img = cvLoadImage( str, 0 );
		if( !img )
			break;
		imageSize = cvGetSize(img);
		imageNames[lr].push_back(str);
		//FIND CHESSBOARDS AND CORNERS THEREIN:
		for( int s = 1; s <= maxScale; s++ )
		{
			IplImage* timg = img;
			if( s > 1 )
			{
				timg = cvCreateImage(cvSize(img->width*s,img->height*s),
					img->depth, img->nChannels );
				cvResize( img, timg, CV_INTER_CUBIC );
			}
			result = cvFindChessboardCorners( timg, cvSize(nx, ny),
				&temp[0], &count,
				CV_CALIB_CB_ADAPTIVE_THRESH |
				CV_CALIB_CB_NORMALIZE_IMAGE);
			if( timg != img )
				cvReleaseImage( &timg );
			if( result || s == maxScale )
				for( j = 0; j < count; j++ )
				{
					temp[j].x /= s;
					temp[j].y /= s;
				}
				if( result )
					break;
		}
		N = pts.size();
		pts.resize(N + n, cvPoint2D32f(0,0));
		active[lr].push_back(result);
		//assert( result != 0 );
		if( result )
		{
			//Calibration will suffer without subpixel interpolation
			cvFindCornerSubPix( img, &temp[0], count,
				cvSize(11, 11), cvSize(-1,-1),
				cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS,
				30, 0.01) );
			copy( temp.begin(), temp.end(), pts.begin() + N );
		}
		if( 0 ) //displayCorners
		{
			printf("%s\n", str);
			IplImage* cimg = cvCreateImage( imageSize, 8, 3 );
			cvCvtColor( img, cimg, CV_GRAY2BGR );
			cvDrawChessboardCorners( cimg, cvSize(nx, ny), &temp[0],
				count, result );
			cvShowImage( "corners", cimg );
			cvReleaseImage( &cimg );
			int c = cvWaitKey(1000);
			if( c == 27 || c == 'q' || c == 'Q' ) //Allow ESC to quit
				exit(-1);
		}
		else
			putchar('.');


		cvReleaseImage( &img );
	}

	vector<CvPoint2D32f> points[2];
	vector<string> imageNamesNew[2];
	nframes = 0;
	char str[128];
	for (int fr = 0; fr < numImages; ++fr){
			if ((active[0][fr] == 0) | (active[1][fr] == 0))
				printf("frame %d was bad\n", fr);
			else
			{
				nframes = nframes + 1;
				printf("frame  %d was good\n", fr);
				if (1) {
					vector<CvPoint2D32f>& pts = points[0];
					vector<CvPoint2D32f>& ptsOld = pointsOld[0];
					N = pts.size();
					pts.resize(N + n, cvPoint2D32f(0,0));
					copy(ptsOld.begin()+fr*n, ptsOld.begin()+fr*n+n, pts.begin()+N);
					
					IplImage* cimg = cvCreateImage( imageSize, 8, 3 );
					sprintf(str, "../images/stereo/left_%d.png",fr + 1);
					IplImage *left = cvLoadImage(str, 0);
					cvCvtColor( left, cimg, CV_GRAY2BGR );
					cvDrawChessboardCorners( cimg, cvSize(nx, ny), &pts[nframes*n-n], 54, 1);
					cvShowImage( "corners", cimg );
					cvReleaseImage( &cimg );
					cvWaitKey(200);

					sprintf(str, "../images/stereo/left_%d.png",fr + 1);
					left = cvLoadImage(str);
					sprintf(str,"../images/stereo/Idleft_%d.png",nframes);
					cvSaveImage(str,left);
				}
				if (1) {
					vector<CvPoint2D32f>& pts = points[1];
					vector<CvPoint2D32f>& ptsOld = pointsOld[1];
					N = pts.size();
					pts.resize(N + n, cvPoint2D32f(0,0));
					copy(ptsOld.begin()+fr*n, ptsOld.begin()+fr*n+n, pts.begin()+N);

					IplImage* cimg = cvCreateImage( imageSize, 8, 3 );
					sprintf(str, "../images/stereo/right_%d.png",fr + 1);
					IplImage *left = cvLoadImage(str, 0);
					cvCvtColor( left, cimg, CV_GRAY2BGR );
					cvDrawChessboardCorners( cimg, cvSize(nx, ny), &pts[nframes*n-n], 54, 1);
					cvShowImage( "corners", cimg );
					cvReleaseImage( &cimg );
					cvWaitKey(200);



					sprintf(str, "../images/stereo/right_%d.png",fr + 1);
					left = cvLoadImage(str);
					sprintf(str,"../images/stereo/Idright_%d.png",nframes);
					cvSaveImage(str,left);
					
				}



			}

	}
	//fclose(f);
	printf("DONE!!! found %d frames!\n", nframes);
	// HARVEST CHESSBOARD 3D OBJECT POINT LIST:
	
	FILE * fid1 = fopen("leftPoints.txt", "w");
	FILE * fid2 = fopen("rightPoints.txt", "w");

	for (int i = 0; i < nframes; ++i){
		for (int j = 0; j < n; ++j){
			
			CvPoint2D32f ptsx = points[0][i*n+j];
			fprintf(fid1, "%f %f ", ptsx.x, ptsx.y);
			ptsx = points[1][i*n+j];
			fprintf(fid2, "%f %f ", ptsx.x, ptsx.y);
		}
		fprintf(fid1, "\n");
		fprintf(fid2, "\n");
	}
	fclose(fid1);
	fclose(fid2);
	
	objectPoints.resize(nframes*n);
	for( int i = 0; i < ny; i++ )
		for( j = 0; j < nx; j++ )
			objectPoints[i*nx + j] =
			cvPoint3D32f(i*squareSize, j*squareSize, 0);
	for( int i = 1; i < nframes; i++ )
		copy( objectPoints.begin(), objectPoints.begin() + n,
		objectPoints.begin() + i*n );
	npoints.resize(nframes,n);
	N = nframes*n;
	CvMat _objectPoints = cvMat(1, N, CV_32FC3, &objectPoints[0] );
	CvMat _imagePoints1 = cvMat(1, N, CV_32FC2, &points[0][0] );
	CvMat _imagePoints2 = cvMat(1, N, CV_32FC2, &points[1][0] );
	CvMat _npoints = cvMat(1,npoints.size(), CV_32S, &npoints[0] );
	cvSetIdentity(&_M1);
	cvSetIdentity(&_M2);
	cvZero(&_D1);
	cvZero(&_D2);

	// CALIBRATE THE STEREO CAMERAS
	printf("Running stereo calibration ...");
	fflush(stdout);
	cvStereoCalibrate( &_objectPoints, &_imagePoints1,
		&_imagePoints2, &_npoints,
		&_M1, &_D1, &_M2, &_D2,
		imageSize, &_R, &_T, &_E, &_F,
		cvTermCriteria(CV_TERMCRIT_ITER+
		CV_TERMCRIT_EPS, 500, 1e-6),
		CV_CALIB_FIX_FOCAL_LENGTH);

	printf(" done\n");

	// SAVE TO FILE CAMERA MATRIXES AND DISTORTION COEFFICIENTS --> FILES
	/*
	FileStorage fs("leftCamera.xml",FileStorage::WRITE);
	fs<<"leftCamera"<<_M1;
	fs.release(); 

	*/


	FILE *f_out = fopen("leftCamera.txt","wt");
	fprintf( f_out, "%e %e %e %e %e %e %e %e %e\n",
		(double)cvGet2D(&_M1,0,0).val[0],(double)cvGet2D(&_M1,0,1).val[0],(double)cvGet2D(&_M1,0,2).val[0],
		(double)cvGet2D(&_M1,1,0).val[0],(double)cvGet2D(&_M1,1,1).val[0],(double)cvGet2D(&_M1,1,2).val[0],
		(double)cvGet2D(&_M1,2,0).val[0],(double)cvGet2D(&_M1,2,1).val[0],(double)cvGet2D(&_M1,2,2).val[0] );

	fprintf( f_out, "%e %e %e %e %e \n",
		(double)cvGet2D(&_D1,0,0).val[0],(double)cvGet2D(&_D1,0,1).val[0],(double)cvGet2D(&_D1,0,2).val[0],
		(double)cvGet2D(&_D1,0,3).val[0],(double)cvGet2D(&_D1,0,4).val[0] );
	fclose( f_out );

	f_out = fopen("rightCamera.txt","wt");
	fprintf( f_out, "%e %e %e %e %e %e %e %e %e\n",
		(double)cvGet2D(&_M2,0,0).val[0],(double)cvGet2D(&_M2,0,1).val[0],(double)cvGet2D(&_M2,0,2).val[0],
		(double)cvGet2D(&_M2,1,0).val[0],(double)cvGet2D(&_M2,1,1).val[0],(double)cvGet2D(&_M2,1,2).val[0],
		(double)cvGet2D(&_M2,2,0).val[0],(double)cvGet2D(&_M2,2,1).val[0],(double)cvGet2D(&_M2,2,2).val[0] );

	fprintf( f_out, "%e %e %e %e %e \n",
		(double)cvGet2D(&_D2,0,0).val[0],(double)cvGet2D(&_D2,0,1).val[0],(double)cvGet2D(&_D2,0,2).val[0],
		(double)cvGet2D(&_D2,0,3).val[0],(double)cvGet2D(&_D2,0,4).val[0] );
	fclose( f_out );

	f_out = fopen("fundamentalMatrix.txt","wt");
	fprintf( f_out, "%e %e %e\n%e %e %e\n%e %e %e\n",
		(double)cvGet2D(&_F,0,0).val[0],(double)cvGet2D(&_F,0,1).val[0],(double)cvGet2D(&_F,0,2).val[0],
		(double)cvGet2D(&_F,1,0).val[0],(double)cvGet2D(&_F,1,1).val[0],(double)cvGet2D(&_F,1,2).val[0],
		(double)cvGet2D(&_F,2,0).val[0],(double)cvGet2D(&_F,2,1).val[0],(double)cvGet2D(&_F,2,2).val[0] );
	fclose( f_out );

	f_out = fopen("essentialMatrix.txt","wt");
	fprintf( f_out, "%e %e %e\n%e %e %e\n%e %e %e\n",
		(double)cvGet2D(&_E,0,0).val[0],(double)cvGet2D(&_E,0,1).val[0],(double)cvGet2D(&_E,0,2).val[0],
		(double)cvGet2D(&_E,1,0).val[0],(double)cvGet2D(&_E,1,1).val[0],(double)cvGet2D(&_E,1,2).val[0],
		(double)cvGet2D(&_E,2,0).val[0],(double)cvGet2D(&_E,2,1).val[0],(double)cvGet2D(&_E,2,2).val[0] );
	fclose( f_out );

	f_out = fopen("RT.txt","wt");
	fprintf( f_out, "%e %e %e %e %e %e %e %e %e\n",
		(double)cvGet2D(&_R,0,0).val[0],(double)cvGet2D(&_R,0,1).val[0],(double)cvGet2D(&_R,0,2).val[0],
		(double)cvGet2D(&_R,1,0).val[0],(double)cvGet2D(&_R,1,1).val[0],(double)cvGet2D(&_R,1,2).val[0],
		(double)cvGet2D(&_R,2,0).val[0],(double)cvGet2D(&_R,2,1).val[0],(double)cvGet2D(&_R,2,2).val[0] );

	fprintf( f_out, "%e %e %e\n",
		(double)cvGet2D(&_T,0,0).val[0],(double)cvGet2D(&_T,1,0).val[0],(double)cvGet2D(&_T,2,0).val[0] );
	fclose( f_out );


	// CALIBRATION QUALITY CHECK
	// because the output fundamental matrix implicitly
	// includes all the output information,
	// we can check the quality of calibration using the
	// epipolar geometry constraint: m2^t*F*m1=0
	vector<CvPoint3D32f> lines[2];
	points[0].resize(N);
	points[1].resize(N);
	_imagePoints1 = cvMat(1, N, CV_32FC2, &points[0][0] );
	_imagePoints2 = cvMat(1, N, CV_32FC2, &points[1][0] );
	lines[0].resize(N);
	lines[1].resize(N);
	CvMat _L1 = cvMat(1, N, CV_32FC3, &lines[0][0]);
	CvMat _L2 = cvMat(1, N, CV_32FC3, &lines[1][0]);
	//Always work in undistorted space
	cvUndistortPoints( &_imagePoints1, &_imagePoints1,
		&_M1, &_D1, 0, &_M1 );
	cvUndistortPoints( &_imagePoints2, &_imagePoints2,
		&_M2, &_D2, 0, &_M2 );
	cvComputeCorrespondEpilines( &_imagePoints1, 1, &_F, &_L1 );
	cvComputeCorrespondEpilines( &_imagePoints2, 2, &_F, &_L2 );
	double avgErr = 0;
	for( int i = 0; i < N; i++ )
	{
		double err = fabs(points[0][i].x*lines[1][i].x +
			points[0][i].y*lines[1][i].y + lines[1][i].z)
			+ fabs(points[1][i].x*lines[0][i].x +
			points[1][i].y*lines[0][i].y + lines[0][i].z);
		avgErr += err;
	}
	printf( "avg err = %g\n", avgErr/(nframes*n) );
	//COMPUTE AND DISPLAY RECTIFICATION
	if( showUndistorted )
	{
		CvMat* mx1 = cvCreateMat( imageSize.height,
			imageSize.width, CV_32F );
		CvMat* my1 = cvCreateMat( imageSize.height,
			imageSize.width, CV_32F );
		CvMat* mx2 = cvCreateMat( imageSize.height,
			imageSize.width, CV_32F );
		CvMat* my2 = cvCreateMat( imageSize.height,
			imageSize.width, CV_32F );
		CvMat* img1r = cvCreateMat( imageSize.height,
			imageSize.width, CV_8U );
		CvMat* img2r = cvCreateMat( imageSize.height,
			imageSize.width, CV_8U );
		CvMat* disp = cvCreateMat( imageSize.height,
			imageSize.width, CV_16S );
		CvMat* vdisp = cvCreateMat( imageSize.height,
			imageSize.width, CV_8U );
		CvMat* pair;
		double R1[3][3], R2[3][3], P1[3][4], P2[3][4];
		CvMat _R1 = cvMat(3, 3, CV_64F, R1);
		CvMat _R2 = cvMat(3, 3, CV_64F, R2);
		// IF BY CALIBRATED (BOUGUET'S METHOD)
		if( useUncalibrated == 0 )
		{
			CvMat *Q = cvCreateMat(4,4,CV_64F);

			CvMat _P1 = cvMat(3, 4, CV_64F, P1);
			CvMat _P2 = cvMat(3, 4, CV_64F, P2);
			cvStereoRectify( &_M1, &_M2, &_D1, &_D2, imageSize,
				&_R, &_T,
				&_R1, &_R2, &_P1, &_P2, Q,
				0,-1,cvSize(0,0)/*CV_CALIB_ZERO_DISPARITY*/ );
			cvSave("reProj.xml", Q);
			isVerticalStereo = fabs(P2[1][3]) > fabs(P2[0][3]);
			//Precompute maps for cvRemap()
			cvInitUndistortRectifyMap(&_M1,&_D1,&_R1,&_P1,mx1,my1);
			cvInitUndistortRectifyMap(&_M2,&_D2,&_R2,&_P2,mx2,my2);
		}
		//OR ELSE HARTLEY'S METHOD
		else if( useUncalibrated == 1 || useUncalibrated == 2 )
			// use intrinsic parameters of each camera, but
			// compute the rectification transformation directly
			// from the fundamental matrix
		{
			double H1[3][3], H2[3][3], iM[3][3];
			CvMat _H1 = cvMat(3, 3, CV_64F, H1);
			CvMat _H2 = cvMat(3, 3, CV_64F, H2);
			CvMat _iM = cvMat(3, 3, CV_64F, iM);
			//Just to show you could have independently used F
			if( useUncalibrated == 2 )
				cvFindFundamentalMat( &_imagePoints1,
				&_imagePoints2, &_F);
			cvStereoRectifyUncalibrated( &_imagePoints1,
				&_imagePoints2, &_F,
				imageSize,
				&_H1, &_H2, 3);
			cvInvert(&_M1, &_iM);
			cvMatMul(&_H1, &_M1, &_R1);
			cvMatMul(&_iM, &_R1, &_R1);
			cvInvert(&_M2, &_iM);
			cvMatMul(&_H2, &_M2, &_R2);
			cvMatMul(&_iM, &_R2, &_R2);
			//Precompute map for cvRemap()
			cvInitUndistortRectifyMap(&_M1,&_D1,&_R1,&_M1,mx1,my1);

			cvInitUndistortRectifyMap(&_M2,&_D1,&_R2,&_M2,mx2,my2);
		}
		else
			assert(0);
		//cvNamedWindow( "rectified", 1 );
		// RECTIFY THE IMAGES AND FIND DISPARITY MAPS
		if( !isVerticalStereo )
			pair = cvCreateMat( imageSize.height, imageSize.width*2,
			CV_8UC3 );
		else
			pair = cvCreateMat( imageSize.height*2, imageSize.width,
			CV_8UC3 );
		//Setup for finding stereo corrrespondences
		printf("Running the Structured Light scanning....\n");

		cvSave("map1x.xml", mx1);
		cvSave("map1y.xml", my1);

		cvSave("map2x.xml", mx2);
		cvSave("map2y.xml", my2);

		IplImage* leftSrc = cvLoadImage("../images/stereo/left_1.png",0);
		IplImage* rightSrc = cvLoadImage("../images/stereo/right_1.png",0);

		IplImage *tmp = cvCreateImage(cvSize(imageW,imageH),IPL_DEPTH_8U,1);

		cvRemap(leftSrc, tmp, mx2, my2, CV_INTER_CUBIC + CV_WARP_FILL_OUTLIERS); 
		//cvNormalize(tmp, tmp, 0, 256, CV_MINMAX );
		cvNamedWindow("left");
		cvShowImage("left",tmp);
		cvWaitKey(-1);
		cvSaveImage("../images/stereo/rec_lef.png",tmp);
		cvRemap(rightSrc, tmp, mx1, my1, CV_INTER_CUBIC + CV_WARP_FILL_OUTLIERS); //
		cvSaveImage("../images/stereo/rec_rig.png",tmp);
		//cvNormalize(tmp, tmp, 0, 256, CV_MINMAX );
		cvShowImage("left",tmp);
		cvWaitKey(-1);
		cvDestroyWindow("left");
		cvReleaseMat( &mx1 );
		cvReleaseMat( &my1 );
		cvReleaseMat( &mx2 );
		cvReleaseMat( &my2 );
		cvReleaseMat( &img1r );
		cvReleaseMat( &img2r );
		cvReleaseMat( &disp );
		printf("Done!\n");
	}

}

int RunSingleCamera( PGRGuid guid, char* filename, int lr, int numImage)
{
	printf("Begina Capturing the image:....\n");
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
		if (lr)
			printf( "Grabbed right image %d\n", numImage);
		else
			printf( "Grabbed left image %d\n", numImage);
		// Create a converted image
		Image convertedImage;

		// Convert the raw image
		error = rawImage.Convert( PIXEL_FORMAT_MONO8, &convertedImage );

		// Create a unique filename
		//        char filename[512];
		//       sprintf( filename, "FlyCapture2Test-%u-%d.png", camInfo.serialNumber, imageCnt );

		// Save the image. If a file format is not passed in, then the file
		// extension is parsed to attempt to determine the file format.
		error = convertedImage.Save( filename );
	}            

	// Stop capturing images
	error = cam.StopCapture(); 

	// Disconnect the camera
	error = cam.Disconnect();
	return 0;
}
void takePhoto(int num){
	char s[128]; 
	int trigger;

	Error error;
	BusManager busMgr;
	int im = 0;
	cvNamedWindow("Show");
	while (im < num) {
		char userInput; 
		cin>>userInput;
		//userInput = cvWaitKey(1);
		if (userInput == 'n') {
			for (int cap = 0; cap < 2; cap++) {
				if (cap)
					sprintf(s,"../images/stereo/right_%d.png",im + 1);
				else
					sprintf(s,"../images/stereo/left_%d.png",im + 1);
				PGRGuid guid;
				error = busMgr. GetCameraFromIndex(cap, &guid);
				RunSingleCamera(guid, s, cap, im);

			}
			
			IplImage* img = cvLoadImage(s);
			cvShowImage("Show",img);
			//cvWaitKey(2000);
			printf("Done..\n");
			im++;
		}
			
		
	}
	/*
	int im = 0; char s[128];
	cvNamedWindow( "Camera l", 1 );
	cvNamedWindow( "Camera r", 1 );
	CvCapture* captureL = cvCaptureFromCAM(0);
	CvCapture* captureR = cvCaptureFromCAM(1);

	cvSetCaptureProperty(captureL, CV_CAP_PROP_FRAME_WIDTH, imageW);
	cvSetCaptureProperty(captureL, CV_CAP_PROP_FRAME_HEIGHT, imageH);

	cvSetCaptureProperty(captureR, CV_CAP_PROP_FRAME_WIDTH, imageW);
	cvSetCaptureProperty(captureR, CV_CAP_PROP_FRAME_HEIGHT, imageH);

	IplImage* imL, *imR;
	while(1){


	IplImage* srcLeft = cvQueryFrame( captureL );
	IplImage* srcRight = cvQueryFrame( captureR );

	cvShowImage( "Camera l", srcLeft );
	cvShowImage( "Camera r", srcRight);


	// capture the image

	cvGrabFrame(captureL); //Get to last frame in buffer 
	cvGrabFrame(captureR); //Get to last frame in buffer 



	if ( cvWaitKey(10) == 27 ) {
	imL =  cvRetrieveFrame(captureL);   
	imR =  cvRetrieveFrame(captureR);   

	sprintf(s,"../images/stereo/left_%d.png",im + 1);
	cvSaveImage(s,imL);

	sprintf(s,"../images/stereo/right_%d.png",im + 1);
	cvSaveImage(s,imR);

	im++;
	}

	if(im >= num) 
	break;
	}


	cvDestroyWindow("Camera l");
	cvDestroyWindow("Camera r");
	*/
}



int main(void)
{
	//CvMat leftCamMat, rightCamMat;
	//CvMat leftDisCoe, rightDisCoe;
	//StereoCalib("stereo_calib.txt", 6, 7, 1, &leftCamMat, &leftDisCoe, &rightCamMat, &rightDisCoe );
	StereoCalib( 8, 6, 0, 0);

	return 0;
}
