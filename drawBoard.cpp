#include "drawBoard.h"

#ifdef DRAW_SOLUTION_BOARD

cv::Mat drawBoard(const int * queenColIdx, const int n)
{
    constexpr int SQUARE_SIZE = 21;
    constexpr float LINE_THICKNESS = 1.0f;

    const int imageSize = SQUARE_SIZE * n;

    cv::Mat result(imageSize, imageSize, CV_8UC3);

    // make checkboard
    for (int y = 0; y < SQUARE_SIZE; y++)
    {
        for (int x = 0; x < SQUARE_SIZE; x++)
        {
            cv::Rect roi(x * SQUARE_SIZE, y * SQUARE_SIZE, SQUARE_SIZE, SQUARE_SIZE);
            const bool black = ((x ^ y) & 1) == 0;
            const cv::Scalar color = black ? cv::Scalar(0, 0, 0) : cv::Scalar(255, 255, 255);
            cv::rectangle(result, roi, color, cv::FILLED);
        }
    }

    // draw attack lines
    if (true)
    {
        for (int yq = 0; yq < SQUARE_SIZE; yq++)
        {
            const int xq = queenColIdx[yq];

            const int y = yq * SQUARE_SIZE + SQUARE_SIZE / 2;
            const int x = xq * SQUARE_SIZE + SQUARE_SIZE / 2;

            cv::line(result, cv::Point(x - SQUARE_SIZE * n, y - SQUARE_SIZE * n), 
                            cv::Point(x + SQUARE_SIZE * n, y + SQUARE_SIZE * n),
                            cv::Scalar(0, 0, 255), LINE_THICKNESS);

            cv::line(result, cv::Point(x - SQUARE_SIZE * n, y + SQUARE_SIZE * n), 
                            cv::Point(x + SQUARE_SIZE * n, y - SQUARE_SIZE * n),
                            cv::Scalar(0, 0, 255), LINE_THICKNESS);
        }
    }    
    
    // draw queens
    for (int yq = 0; yq < SQUARE_SIZE; yq++)
    {
        const int xq = queenColIdx[yq];

        const int y = yq * SQUARE_SIZE + SQUARE_SIZE / 2;
        const int x = xq * SQUARE_SIZE + SQUARE_SIZE / 2;

        cv::circle(result, cv::Point2f(x, y), SQUARE_SIZE / 2 - 2, cv::Scalar(128, 128, 128), cv::FILLED);
    }

    return result;
}

#endif