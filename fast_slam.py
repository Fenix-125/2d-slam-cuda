import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from slam.occupy_grid import occupy_grid
from slam.scan_matcher import scan_matcher
import math
import copy


class particle_filter:
    def __init__(self, numParticles, ogParameters, smParameters):
        self.numParticles = numParticles
        self.particles = []
        self.p_init(ogParameters, smParameters)
        self.step = 0
        self.prevMatchedReading = None
        self.prevRawReading = None
        self.particlesTrajectory = []

    def p_init(self, ogParameters, smParameters):
        for _ in range(self.numParticles):
            p = particle(ogParameters, smParameters)
            self.particles.append(p)

    def p_update(self, reading, count):
        for i in range(self.numParticles):
            self.particles[i].update(reading, count)

    # if sum of variances is too big the weight are decided by Sum of all particle weights
    def p_balance(self):
        self.p_normalize()
        variance = 0
        for i in range(self.numParticles):
            variance += (self.particles[i].weight - 1 / self.numParticles) ** 2

        print(f"Variance {variance}")
        if variance > (self.numParticles ** 2 - self.numParticles - .000000000000001) / (self.numParticles ** 2):
            return True
        else:
            return False

    def p_normalize(self):
        weightSum = 0
        for i in range(self.numParticles):
            weightSum += self.particles[i].weight
        for i in range(self.numParticles):
            self.particles[i].weight /= weightSum

    def resample(self):
        weights = np.zeros(self.numParticles)
        tempParticles = []
        for i in range(self.numParticles):
            weights[i] = self.particles[i].weight
            tempParticles.append(copy.deepcopy(self.particles[i]))
        resampledParticlesIdx = np.random.choice(np.arange(self.numParticles), self.numParticles, p=weights)
        for i in range(self.numParticles):
            self.particles[i] = copy.deepcopy(tempParticles[resampledParticlesIdx[i]])
            self.particles[i].weight = 1 / self.numParticles


class particle:
    def __init__(self, ogParameters, smParameters):
        initMapXLength, initMapYLength, initXY, unitGridSize, lidarFOV, lidarMaxRange, numSamplesPerRev, \
        wallThickness = ogParameters
        scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, moveRSigma, maxMoveDeviation, \
        turnSigma, missMatchProbAtCoarse, coarseFactor = smParameters
        p_map = occupy_grid(initMapXLength, initMapYLength, initXY, unitGridSize, lidarFOV, numSamplesPerRev,
                              lidarMaxRange, wallThickness)
        sm = scan_matcher(p_map, scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, moveRSigma,
                         maxMoveDeviation, turnSigma, missMatchProbAtCoarse, coarseFactor)
        self.map = p_map
        self.sm = sm
        self.xTrajectory = []
        self.yTrajectory = []
        self.weight = 1

    def update_est_pos(self, currentRawReading):
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

    def get_move_direction(self, matchedReading):
        x, y, theta, _ = matchedReading['x'], matchedReading['y'], matchedReading['theta'], matchedReading['range']
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
            estimatedReading, estMovingDist, estMovingTheta, rawMovingTheta = self.update_est_pos(
                currentRawReading)
            matchedReading, confidence = self.sm.matchScan(estimatedReading, estMovingDist, estMovingTheta, count,
                                                           matchMax=False)
            self.prevRawMovingTheta = rawMovingTheta
            self.prevMatchedMovingTheta = self.get_move_direction(matchedReading)
        self.update_trajectory(matchedReading)
        self.map.updateOccupancyGrid(matchedReading)
        self.prevMatchedReading, self.prevRawReading = matchedReading, reading
        self.weight *= confidence

    def update_trajectory(self, matchedReading):
        x, y = matchedReading['x'], matchedReading['y']
        self.xTrajectory.append(x)
        self.yTrajectory.append(y)

    def plot(self):
        plt.figure(figsize=(19.20, 19.20))
        plt.scatter(self.xTrajectory[0], self.yTrajectory[0], color='r', s=500)
        colors = iter(cm.rainbow(np.linspace(1, 0, len(self.xTrajectory) + 1)))
        for i in range(len(self.xTrajectory)):
            plt.scatter(self.xTrajectory[i], self.yTrajectory[i], color=next(colors), s=35)
        plt.scatter(self.xTrajectory[-1], self.yTrajectory[-1], color=next(colors), s=500)
        plt.plot(self.xTrajectory, self.yTrajectory)
        self.map.plotOccupancyGrid([-13, 20], [-25, 7], plotThreshold=False)


def process_sensor_data(pf, sensorData, plotTrajectory=True):
    count = 0
    plt.figure(figsize=(19.20, 19.20))
    for time in sorted(sensorData.keys()):
        count += 1
        print(count)
        pf.p_update(sensorData[time], count)
        if pf.p_balance():
            pf.resample()
            print("resample")

        plt.figure(figsize=(19.20, 19.20))
        maxWeight = -1
        for particle in pf.particles:
            if maxWeight < particle.weight:
                maxWeight = particle.weight
                bestParticle = particle
                plt.plot(particle.xTrajectory, particle.yTrajectory)

        xRange, yRange = [-13, 20], [-25, 7]
        ogMap = bestParticle.map.occupancyGridVisited / bestParticle.map.occupancyGridTotal
        xIdx, yIdx = bestParticle.map.convertRealXYToMapIdx(xRange, yRange)
        ogMap = ogMap[yIdx[0]: yIdx[1], xIdx[0]: xIdx[1]]
        ogMap = np.flipud(1 - ogMap)
        plt.imshow(ogMap, cmap='gray', extent=[xRange[0], xRange[1], yRange[0], yRange[1]])
        plt.savefig('Output/' + str(count).zfill(3) + '.png')
        plt.close()

    maxWeight = 0
    for particle in pf.particles:
        particle.plotParticle()
        if maxWeight < particle.weight:
            maxWeight = particle.weight
            bestParticle = particle
    bestParticle.plotParticle()


def readJson(jsonFile):
    with open(jsonFile, 'r') as f:
        input_j = json.load(f)
        return input_j['map']


def main():
    initMapXLength = 50  # in Meters
    initMapYLength = 50  # in Meters
    unitGridSize = 0.02  # in Meters
    lidarFOV = np.pi  # in Meters
    lidarMaxRange = 10  # in Meters

    scanMatchSearchRadius = 1.4
    scanMatchSearchHalfRad = 0.25
    scanSigmaInNumGrid = 2
    wallThickness = 5 * unitGridSize
    moveRSigma = 0.1
    maxMoveDeviation = 0.25
    turnSigma = 0.3
    missMatchProbAtCoarse = 0.15
    coarseFactor = 5  # TODO: What is coarseFactor?

    with open("DataSet/PreprocessedData/intel_gfs", 'r') as f:
        sensorData = json.load(f)['map']
    # {
    #     "976052890.244111": {
    #         "range": [
    #             1.09,
    #             ...
    #         ], # 180 components
    #         "theta": -0.463373,
    #         "x": 0.698,
    #         "y": -0.015
    #     },
    #     ...
    # }

    # Get how many points per revolution; angel view discrete parts count
    numSamplesPerView = len(sensorData[list(sensorData)[0]]['range'])
    initXY = sensorData[sorted(sensorData.keys())[0]]
    numParticles = 10

    #              [ map width    , map height   ,start_pos, map cell size, view angel, max dist view, angel parts count
    ogParameters = [initMapXLength, initMapYLength, initXY, unitGridSize, lidarFOV, lidarMaxRange, numSamplesPerView,
                    wallThickness]
    smParameters = [scanMatchSearchRadius, scanMatchSearchHalfRad, scanSigmaInNumGrid, moveRSigma, maxMoveDeviation,
                    turnSigma,
                    missMatchProbAtCoarse, coarseFactor]
    pf = particle_filter(numParticles, ogParameters, smParameters)
    process_sensor_data(pf, sensorData, plotTrajectory=True)


if __name__ == '__main__':
    main()
