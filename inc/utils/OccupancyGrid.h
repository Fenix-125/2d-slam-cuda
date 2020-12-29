//
// Created by fenix on 27.12.20.
//

#ifndef INC_2D_SLAM_CUDA_OCCUPANCYGRID_H
#define INC_2D_SLAM_CUDA_OCCUPANCYGRID_H

#include <slam_integrate_header.h>

#include <utility>
#include <cmath>
#include <vector>
#include <cinttypes>

class OccupancyGrid {
    struct spokesGrid_t {
        size_t xGrid, yGrid, bearingIdxGrid, rangeIdxGrid;
    };

    struct angel_t {
        double radByX, radByY, radByR;
    };
public:
    size_t mapXLength, mapYLength, unitGridSize, numSamplesPerRev, lidarMaxRange, wallThickness, numSpokes, spokesStartIdx;
    size_t2 initXY  /* first = x */, mapXLim, mapYLim;
    double angularStep, lidarFOV, radByX, radByY, radByR;
    std::vector<std::vector<double>> OccupancyGridX, OccupancyGridY;
    std::vector<std::vector<size_t>> occupancyGridVisited, occupancyGridTotal;

    OccupancyGrid(size_t mapXLength, size_t mapYLength, size_t2 initXY, double unitGridSize, double lidarFOV,
                  size_t numSamplesPerRev, size_t lidarMaxRange, size_t wallThickness) {
        size_t xNum = mapXLength / unitGridSize;
        size_t yNum = mapYLength / unitGridSize;
        auto x = np::linspace(-xNum * unitGridSize / 2., xNum * unitGridSize / 2., xNum + 1u);
        np::add_1d(x, initXY.first);
        auto y = np::linspace(-xNum * unitGridSize / 2., xNum * unitGridSize / 2., yNum + 1u);
        np::add_1d(y, initXY.second);

        auto tmp_mesh = np::meshgrid(x, y);
        OccupancyGridX = std::move(tmp_mesh.first);
        OccupancyGridY = std::move(tmp_mesh.second);
        this->occupancyGridVisited = np::ones<size_t>(xNum + 1u, yNum + 1u);
        occupancyGridTotal = np::ones<size_t>(xNum + 1u, yNum + 1u);
        np::mul_2d(occupancyGridTotal, 2);
        this->unitGridSize = unitGridSize;
        this->lidarFOV = lidarFOV; // TODO: check type
        this->lidarMaxRange = lidarMaxRange;
        this->wallThickness = wallThickness;
        mapXLim = size_t2{OccupancyGridX[0][0], OccupancyGridX[0][-1]};
        mapYLim = size_t2{OccupancyGridY[0][0], OccupancyGridY[-1][0]};
        this->numSamplesPerRev = numSamplesPerRev;
        angularStep = lidarFOV / numSamplesPerRev;
        numSpokes = static_cast<size_t>(round(2 * M_PI / angularStep));
        const auto tmp_sp_grid = spokesGrid();
        const auto tmp_angel = itemizeSpokesGrid(tmp_sp_grid);
        radByX = tmp_angel.radByX;
        radByY = tmp_angel.radByY;
        radByR = tmp_angel.radByR;
        // theta= 0 is x direction. spokes=0 is y direction, spokesStartIdx is the first ray of lidar scan direction.
        // spokes increase counter-clockwise
        spokesStartIdx = static_cast<size_t>(((numSpokes / 2u - this->numSamplesPerRev) / 2u) % numSpokes);
    }

    spokesGrid_t spokesGrid() {
        spokesGrid_t res{};
        // 0th ray is at south, then counter-clock wise increases. Theta 0 is at east.
        size_t numHalfElem = lidarMaxRange / unitGridSize;
        auto bearingIdxGrid = np::zeros<size_t>(2u * numHalfElem + 1u, 2u * numHalfElem + 1u);
        auto x = np::linspace(-lidarMaxRange, lidarMaxRange, 2 * numHalfElem + 1);
        auto y = np::linspace(-lidarMaxRange, lidarMaxRange, 2 * numHalfElem + 1);
        auto tmp_mesh = np::meshgrid(x, y);
        auto xGrid = std::move(tmp_mesh.first), yGrid = std::move(tmp_mesh.second);
        for (size_t j = 0; j < bearingIdxGrid.size(); ++j) {
            for (size_t i = numHalfElem + 1u; i < 2u * numHalfElem + 1u; ++i) {
                bearingIdxGrid[j][i] = static_cast<size_t>(round((M_PI / 2u + atan(static_cast<double>(yGrid[j][i])
                                                                                   /
                                                                                   static_cast<double>(xGrid[j][i]))) /
                                                                 M_PI / 2u * numSpokes - 0.5));
            }
        }
//        bearingIdxGrid[:, 0: numHalfElem] = np.fliplr(np.flipud(bearingIdxGrid))[:, 0: numHalfElem] + int(self.numSpokes / 2)
//        bearingIdxGrid[numHalfElem + 1: 2 * numHalfElem + 1, numHalfElem] = int(self.numSpokes / 2)
//        rangeIdxGrid = np.sqrt(xGrid**2 + yGrid**2)
//        return xGrid, yGrid, bearingIdxGrid, rangeIdxGrid
        return spokesGrid_t{};
    }

    angel_t itemizeSpokesGrid(spokesGrid_t sp_grid) {
        return angel_t{};
    }

    void checkAndExapndOG(size_t2 xRangeList, size_t2 yRangeList) {

    }

    double_size_t2 convertRealXYToMapIdx(size_t2 xRangeList, size_t2 yRangeList) {

    }
};


#endif //INC_2D_SLAM_CUDA_OCCUPANCYGRID_H
