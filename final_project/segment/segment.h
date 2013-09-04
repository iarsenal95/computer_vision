/*
 * segment.h
 *
 *  Created on: Nov 10, 2012
 *      Author: thewbp
 */

#ifndef SEGMENT_H_
#define SEGMENT_H_

#include <string.h>

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

typedef struct _ImageSegment
{
	char 	imName[32];	// image name
	CvSize 	imSize;		// image size
	cv::Mat	segMat;	// a matrix with the segment distribution
	int		nSeg;		// number of segments
	int		*nPixels;	// number of pixels in each segment
	cv::Mat	adjMat;	// the adjacent matrix of segments
}ImageSegment;

void segmentImage(char *inputName, char *outputName, float sigma, int k, int minSize, ImageSegment *segment);

#endif /* SEGMENT_H_ */
