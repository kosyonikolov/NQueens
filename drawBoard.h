#ifndef DRAW_BOARD_H
#define DRAW_BOARD_H

#ifdef DRAW_SOLUTION_BOARD
#include <opencv2/opencv.hpp>

cv::Mat drawBoard(const int * queenColIdx, const int n);
#endif

#endif