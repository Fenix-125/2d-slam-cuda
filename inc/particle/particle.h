//
// Created by fenix on 10.12.20.
//

#ifndef INC_2D_SLAM_CUDA_PARTICLE_H
#define INC_2D_SLAM_CUDA_PARTICLE_H

#include <vector>
#include <cinttypes>

namespace ucu {
    struct map_params_t {
        size_t initMapXLength, initMapYLength, initXY, unitGridSize;
        double lidarFOV, lidarMaxRange;
        uint32_t numSamplesPerRev;
        double wallThickness;
    };

    struct stat_params_t {
        double scanMatchSearchRadius,
            scanMatchSearchHalfRad,
            scanSigmaInNumGrid,
            moveRSigma,
            maxMoveDeviation,
            turnSigma,
            missMatchProbAtCoarse,
            coarseFactor;
    };

    class particle {
//        map_m
//        sm_m
//        xTrajectory_m
//        yTrajectory_m
//        weight_m
        particle(ogParameters, smParameters) {
            map_m = OccupancyGrid(ogParameters);
            sm_m = ScanMatcher(map_m, smParameters);
            xTrajectory_m = [];
            yTrajectory_m = [];
            weight_m = 1.;
        }


        updateEstimatedPose(currentRawReading) {
                estimatedTheta = self.prevMatchedReading['theta'] + currentRawReading['theta'] - self.prevRawReading['theta']
                estimatedReading = {'x': self.prevMatchedReading['x'], 'y': self.prevMatchedReading['y'],
                            'theta': estimatedTheta,
                            'range': currentRawReading['range']}

                dx, dy = currentRawReading['x'] - self.prevRawReading['x'], currentRawReading['y'] - self.prevRawReading['y']
                estMovingDist = math.sqrt(dx ** 2 + dy ** 2)
                rawX, rawY, prevRawX, prevRawY = currentRawReading['x'], currentRawReading['y'], self.prevRawReading['x'], \
                                         self.prevRawReading['y']
                rawXMove, rawYMove = rawX - prevRawX, rawY - prevRawY
                rawMove = math.sqrt((rawX - prevRawX) ** 2 + (rawY - prevRawY) ** 2)

                if rawMove > 0.3:
                if self.prevRawMovingTheta is not None:
                if rawYMove > 0:
                rawMovingTheta = math.acos(rawXMove / rawMove)  # between -pi and +pi
                else:
                rawMovingTheta = -math.acos(rawXMove / rawMove)
                rawTurnTheta = rawMovingTheta - self.prevRawMovingTheta
                estMovingTheta = self.prevMatchedMovingTheta + rawTurnTheta
                else:
                if rawYMove > 0:
                rawMovingTheta = math.acos(rawXMove / rawMove)  # between -pi and +pi
                else:
                rawMovingTheta = -math.acos(rawXMove / rawMove)
                estMovingTheta = None
                else:
                rawMovingTheta = None
                estMovingTheta = None

                return estimatedReading, estMovingDist, estMovingTheta, rawMovingTheta
        }

                def getMovingTheta(self, matchedReading):
        x, y, theta, range = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
        prevX, prevY = self.xTrajectory[-1], self.yTrajectory[-1]
        xMove, yMove = x - prevX, y - prevY
                move = math.sqrt(xMove ** 2 + yMove ** 2)
        if move != 0:
        if yMove > 0:
        movingTheta = math.acos(xMove / move)
        else:
        movingTheta = -math.acos(xMove / move)
        else:
        movingTheta = None
        return movingTheta

                def update(self, reading, count):
        if count == 1:
        self.prevRawMovingTheta, self.prevMatchedMovingTheta = None, None
        matchedReading, confidence = reading, 1
        else:
        currentRawReading = reading
        estimatedReading, estMovingDist, estMovingTheta, rawMovingTheta = self.updateEstimatedPose(
                currentRawReading)
        matchedReading, confidence = self.sm.matchScan(estimatedReading, estMovingDist, estMovingTheta, count,
                                                       matchMax=False)
        self.prevRawMovingTheta = rawMovingTheta
        self.prevMatchedMovingTheta = self.getMovingTheta(matchedReading)
        self.updateTrajectory(matchedReading)
        self.map.updateOccupancyGrid(matchedReading)
        self.prevMatchedReading, self.prevRawReading = matchedReading, reading
        self.weight *= confidence

                def updateTrajectory(self, matchedReading):
        x, y = matchedReading['x'], matchedReading['y']
        self.xTrajectory.append(x)
        self.yTrajectory.append(y)

        def plotParticle(self):
        plt.figure(figsize=(19.20, 19.20))
        plt.scatter(self.xTrajectory[0], self.yTrajectory[0], color='r', s=500)
        colors = iter(cm.rainbow(np.linspace(1, 0, len(self.xTrajectory) + 1)))
        for i in range(len(self.xTrajectory)):
        plt.scatter(self.xTrajectory[i], self.yTrajectory[i], color=next(colors), s=35)
        plt.scatter(self.xTrajectory[-1], self.yTrajectory[-1], color=next(colors), s=500)
        plt.plot(self.xTrajectory, self.yTrajectory)
        self.map.plotOccupancyGrid([-13, 20], [-25, 7], plotThreshold=False)

    };
}


#endif //INC_2D_SLAM_CUDA_PARTICLE_H
