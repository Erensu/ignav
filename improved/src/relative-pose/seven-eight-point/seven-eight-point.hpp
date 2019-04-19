#ifndef ONE_POINT_HPP
#define ONE_POINT_HPP

#include <opencv2/opencv.hpp>

int cvFindFundamentalMat( const CvMat* points1, const CvMat* points2,
                          CvMat* fmatrix, int method,
                          double param1, double param2, CvMat* mask );
int cvFindHomography( const CvMat* objectPoints, const CvMat* imagePoints,
                      CvMat* __H, int method, double ransacReprojThreshold,
                      CvMat* mask );

cv::Mat cv::findHomography( InputArray _points1, InputArray _points2,
                            OutputArray _mask, int method,
                            double ransacReprojThreshold );

cv::Mat cv::findFundamentalMat( InputArray _points1, InputArray _points2,
                                int method, double param1, double param2,
                                OutputArray _mask );

cv::Mat cv::findFundamentalMat( InputArray _points1, InputArray _points2,
                                OutputArray _mask, int method, double param1,
                                double param2 );
#endif