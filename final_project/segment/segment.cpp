/*
 * segment.cpp
 *
 *  Created on: Nov 10, 2012
 *      Author: thewbp
 */

#include "segment.h"

#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "segment-image.h"

using namespace std;

void groupSegments(ImageSegment *segment, IplImage *image)
{
	segment->segMat = cv::Mat::zeros(image->height, image->width, CV_32SC1);
	segment->nPixels =  (int *)calloc(segment->nSeg, sizeof(int));

	int *hash_table = (int *)calloc(segment->imSize.width * segment->imSize.height, sizeof(int));

	int gid = 1;
	for (int y = 0; y < image->height; y++)
	{
		uchar *ptr = (uchar *)(image->imageData + y * image->widthStep);
		for (int x = 0; x < image->width; x++)
		{
			int comp = 10000 * ptr[3 * x] + 100 * ptr[3 * x + 1] + ptr[3 * x + 2];
			if (hash_table[comp] == 0)
			{
				hash_table[comp] = gid++;
			}
			segment->segMat.at<int>(y, x) = hash_table[comp] - 1;
			segment->nPixels[hash_table[comp] - 1]++;
		}
	}

	free(hash_table);
}

void getAdjMat(ImageSegment *segment)
{
	int s1, s2, s3;
	segment->adjMat = cv::Mat::zeros(segment->nSeg,segment->nSeg, CV_8UC1);
	for (int y = 0; y < segment->imSize.height - 1; y++)
	{
		for (int x = 0; x < segment->imSize.width - 1; x++)
		{
			s1 = segment->segMat.at<int>(y, x);
			s2 = segment->segMat.at<int>(y, x + 1);
			s3 = segment->segMat.at<int>(y + 1, x);
			segment->adjMat.at<unsigned char>(s1, s2) = 1;
			segment->adjMat.at<unsigned char>(s2, s1) = 1;
			segment->adjMat.at<unsigned char>(s1, s3) = 1;
			segment->adjMat.at<unsigned char>(s3, s1) = 1;
		}
	}
}

void segmentImage(char *inputName, char *outputName, float sigma, int k, int minSize, ImageSegment *segment)
{
	cout << "Loading input image." << endl;
	image<rgb> *input = loadPPM(inputName);

	cout << "Processing..." << endl;
	int num_ccs;
	image<rgb> *seg = segment_image(input, sigma, k, minSize, &num_ccs);
	savePPM(seg, outputName);

	printf("Got %d components!\n", num_ccs);

	// image name
	strcpy(segment->imName, outputName);
	// size
	segment->imSize = cvSize(seg->width(), seg->height());
	// segments
	segment->nSeg = num_ccs;

	IplImage *img = cvLoadImage(outputName);

	groupSegments(segment, img);
	getAdjMat(segment);

//	int n = 0;
//	for (int i = 0; i < 122; i++)
//	{
//		for (int j = 0; j < 122; j++)
//		{
//			if (segment->adjMat.at<unsigned char>(i, j) == 1)
//				n++;
//		}
//	}
//
//	cout << "total adj = " << n << endl;

	// just to show the segmentation(use hashed rgb color)
//	cvNamedWindow("Segmented Image", CV_WINDOW_AUTOSIZE);
//	cvShowImage("Segmented Image", img);
//	cvWaitKey(0);
//	cvReleaseImage(&img);
//	cvDestroyWindow("Segmented Image");

	cvReleaseImage(&img);
}

void releaseImageSegment(ImageSegment *segment)
{

}
