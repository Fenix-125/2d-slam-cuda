//
// Created by fenix on 27.12.20.
//

#ifndef INC_2D_SLAM_CUDA_SCANMATCHER_H
#define INC_2D_SLAM_CUDA_SCANMATCHER_H

#include <vector>
#include <cmath>
#include <algorithm>

#include <slam_integrate_header.h>
#include <utils/OccupancyGrid.h>

class ScanMatcher {
    OccupancyGrid og;
    double searchRadius, searchHalfRad, scanSigmaInNumGrid, moveRSigma, maxMoveDeviation,
            turnSigma, missMatchProbAtCoarse, coarseFactor;

    ScanMatcher(OccupancyGrid og, double searchRadius, double searchHalfRad,
                double scanSigmaInNumGrid, double moveRSigma, double maxMoveDeviation,
                double turnSigma, double missMatchProbAtCoarse, double coarseFactor) :
            og{og},
            searchRadius{searchRadius},
            searchHalfRad{searchHalfRad},
            scanSigmaInNumGrid{scanSigmaInNumGrid},
            moveRSigma{moveRSigma},
            maxMoveDeviation{maxMoveDeviation},
            turnSigma{turnSigma},
            missMatchProbAtCoarse{missMatchProbAtCoarse},
            coarseFactor{coarseFactor} {}

    std::vector<double> frameSearchSpace(double estimatedX, double estimatedY, double unitLength, double sigma,
                                         double missMatchProbAtCoarse) {
        double maxScanRadius = 1.1 * og.lidarMaxRange + searchRadius;
        double_size_t2 tmp_dd2;
        size_t2 xRangeList{estimatedX - maxScanRadius, estimatedX + maxScanRadius};
        size_t2 yRangeList{estimatedY - maxScanRadius, estimatedY + maxScanRadius};

        auto idxEndX = static_cast<size_t>((xRangeList.second - xRangeList.first) / unitLength);
        auto idxEndY = static_cast<size_t>((yRangeList.second - yRangeList.first) / unitLength);
        matrix2d_double searchSpace = np::ones(idxEndY + 1u, idxEndX + 1u);
        for_each(searchSpace.begin(), searchSpace.end(),
                 [&](std::vector<double> &v) { for (auto &el : v) { el *= log(missMatchProbAtCoarse); }});

        og.checkAndExapndOG(xRangeList, yRangeList);
        tmp_dd2 = og.convertRealXYToMapIdx(xRangeList, yRangeList);
        size_t2 xRangeListIdx = tmp_dd2.first, yRangeListIdx = tmp_dd2.second;
        matrix2d_double ogMap = og.occupancyGridVisited[yRangeListIdx.first : yRangeListIdx[1],
                                        xRangeListIdx.first: xRangeListIdx[1]] /
                                og.occupancyGridTotal[yRangeListIdx.first : yRangeListIdx[1],
                                        xRangeListIdx.first: xRangeListIdx[1]];
    }
};


#endif //INC_2D_SLAM_CUDA_SCANMATCHER_H
