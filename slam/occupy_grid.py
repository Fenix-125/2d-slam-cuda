import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class occupy_grid:
    def __init__(self, mapXLength, mapYLength, initXY, unitGridSize, lidarFOV, numSamplesPerRev, lidarMaxRange,
                 wallThickness):
        xNum = int(mapXLength / unitGridSize)
        yNum = int(mapYLength / unitGridSize)
        x = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=xNum + 1) + initXY['x']
        y = np.linspace(-xNum * unitGridSize / 2, xNum * unitGridSize / 2, num=yNum + 1) + initXY['y']
        self.occupy_by_x, self.occupy_by_y = np.meshgrid(x, y)
        self.visited = np.ones((xNum + 1, yNum + 1))
        self.total = 2 * np.ones((xNum + 1, yNum + 1))
        self.unitGridSize = unitGridSize
        self.lidarFOV = lidarFOV
        self.lidarMaxRange = lidarMaxRange
        self.wallThickness = wallThickness
        self.mapXLim = [self.occupy_by_x[0, 0], self.occupy_by_x[0, -1]]
        self.mapYLim = [self.occupy_by_y[0, 0], self.occupy_by_y[-1, 0]]
        self.numSamplesPerRev = numSamplesPerRev
        self.angularStep = lidarFOV / numSamplesPerRev
        self.numSpokes = int(np.rint(2 * np.pi / self.angularStep))
        xGrid, yGrid, bearingIdxGrid, rangeIdxGrid = self.init_grid()
        radByX, radByY, radByR = self.itemize(xGrid, yGrid, bearingIdxGrid, rangeIdxGrid)
        self.radByX = radByX
        self.radByY = radByY
        self.radByR = radByR
        self.spokesStartIdx = int(((self.numSpokes / 2 - self.numSamplesPerRev) / 2) % self.numSpokes)

    # INIT GRID
    def init_grid(self):
        # 0th ray is at south, then counter-clock wise increases. Theta 0 is at east.
        numHalfElem = int(self.lidarMaxRange / self.unitGridSize)
        bearingIdxGrid = np.zeros((2 * numHalfElem + 1, 2 * numHalfElem + 1))
        x = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        y = np.linspace(-self.lidarMaxRange, self.lidarMaxRange, 2 * numHalfElem + 1)
        xGrid, yGrid = np.meshgrid(x, y)
        bearingIdxGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] = np.rint((np.pi / 2 + np.arctan(
            yGrid[:, numHalfElem + 1: 2 * numHalfElem + 1] / xGrid[:, numHalfElem + 1: 2 * numHalfElem + 1]))
                                                                          / np.pi / 2 * self.numSpokes - 0.5).astype(
            int)
        bearingIdxGrid[:, 0: numHalfElem] = np.fliplr(np.flipud(bearingIdxGrid))[:, 0: numHalfElem] + int(
            self.numSpokes / 2)
        bearingIdxGrid[numHalfElem + 1: 2 * numHalfElem + 1, numHalfElem] = int(self.numSpokes / 2)
        rangeIdxGrid = np.sqrt(xGrid ** 2 + yGrid ** 2)
        return xGrid, yGrid, bearingIdxGrid, rangeIdxGrid

    # todo: EXPLAIN
    def itemize(self, xGrid, yGrid, bearingIdxGrid, rangeIdxGrid):
        # Due to discretization, later theta added could lead to up to 1 deg discretization error
        radByX = []
        radByY = []
        radByR = []
        for i in range(self.numSpokes):
            idx = np.argwhere(bearingIdxGrid == i)
            radByX.append(xGrid[idx[:, 0], idx[:, 1]])
            radByY.append(yGrid[idx[:, 0], idx[:, 1]])
            radByR.append(rangeIdxGrid[idx[:, 0], idx[:, 1]])
        return radByX, radByY, radByR

    # EXPAND THE MAP IN TO NEEDED DIRECTION
    def fit_size(self, position, axis):
        gridShape = self.visited.shape
        if axis == 0:
            insertion = np.ones((int(gridShape[0] / 5), gridShape[1]))
            if position == 0:
                x = self.occupy_by_x[0]
                y = np.linspace(self.mapYLim[0] - int(gridShape[0] / 5) * self.unitGridSize, self.mapYLim[0],
                                num=int(gridShape[0] / 5), endpoint=False)
            else:
                x = self.occupy_by_x[0]
                y = np.linspace(self.mapYLim[1] + self.unitGridSize,
                                self.mapYLim[1] + (int(gridShape[0] / 5)) * self.unitGridSize,
                                num=int(gridShape[0] / 5), endpoint=False)
        else:
            insertion = np.ones((gridShape[0], int(gridShape[1] / 5)))
            if position == 0:
                y = self.occupy_by_y[:, 0]
                x = np.linspace(self.mapXLim[0] - int(gridShape[1] / 5) * self.unitGridSize, self.mapXLim[0],
                                num=int(gridShape[1] / 5), endpoint=False)
            else:
                y = self.occupy_by_y[:, 0]
                x = np.linspace(self.mapXLim[1] + self.unitGridSize,
                                self.mapXLim[1] + (int(gridShape[1] / 5)) * self.unitGridSize,
                                num=int(gridShape[1] / 5), endpoint=False)
        self.visited = np.insert(self.visited, [position], insertion, axis=axis)
        self.total = np.insert(self.total, [position], 2 * insertion, axis=axis)
        xv, yv = np.meshgrid(x, y)
        self.occupy_by_x = np.insert(self.occupy_by_x, [position], xv, axis=axis)
        self.occupy_by_y = np.insert(self.occupy_by_y, [position], yv, axis=axis)
        self.mapXLim[0] = self.occupy_by_x[0, 0]
        self.mapXLim[1] = self.occupy_by_x[0, -1]
        self.mapYLim[0] = self.occupy_by_y[0, 0]
        self.mapYLim[1] = self.occupy_by_y[-1, 0]

    # DECODE THE NEEDED DIRECTION AND EXPAND GREED
    def expandOccupancyGrid(self, expandDirection):
        gridShape = self.visited.shape
        if expandDirection == 1:
            self.fit_size(0, 1)
        elif expandDirection == 2:
            self.fit_size(gridShape[1], 1)
        elif expandDirection == 3:
            self.fit_size(0, 0)
        else:
            self.fit_size(gridShape[0], 0)

    # GET FRO X Y COORDINATES X Y INDEXES
    def convertRealXYToMapIdx(self, x, y):
        # mapXLim is (2,) array for left and right limit, same for mapYLim
        xIdx = (np.rint((x - self.mapXLim[0]) / self.unitGridSize)).astype(int)
        yIdx = (np.rint((y - self.mapYLim[0]) / self.unitGridSize)).astype(int)
        return xIdx, yIdx

    # CHACK IF EXPANSION OF THE MAP IS NEEDED
    def checkMapToExpand(self, x, y):
        if any(x < self.mapXLim[0]):
            return 1
        elif any(x > self.mapXLim[1]):
            return 2
        elif any(y < self.mapYLim[0]):
            return 3
        elif any(y > self.mapYLim[1]):
            return 4
        else:
            return -1

    # EXPAND THE MAP IF NEEDED
    def checkAndExapndOG(self, x, y):
        """check x, y (vector points) are inside OG. If not, expand OG."""
        expandDirection = self.checkMapToExpand(x, y)
        while expandDirection != -1:
            self.expandOccupancyGrid(expandDirection)
            expandDirection = self.checkMapToExpand(x, y)

    # UPDATE THE MAP DUE TO ODOMETERY
    def updateOccupancyGrid(self, reading, dTheta=0, update=True):
        x, y, theta, rMeasure = reading['x'], reading['y'], reading['theta'], reading['range']
        theta += dTheta
        rMeasure = np.asarray(rMeasure)
        spokesOffsetIdxByTheta = int(np.rint(theta / (2 * np.pi) * self.numSpokes))
        emptyXList, emptyYList, occupiedXList, occupiedYList = [], [], [], []
        for i in range(self.numSamplesPerRev):
            spokeIdx = int(np.rint((self.spokesStartIdx + spokesOffsetIdxByTheta + i) % self.numSpokes))
            xAtSpokeDir = self.radByX[spokeIdx]
            yAtSpokeDir = self.radByY[spokeIdx]
            rAtSpokeDir = self.radByR[spokeIdx]
            if rMeasure[i] < self.lidarMaxRange:
                emptyIdx = np.argwhere(rAtSpokeDir < rMeasure[i] - self.wallThickness / 2)
            else:
                emptyIdx = []
            occupiedIdx = np.argwhere(
                (rAtSpokeDir > rMeasure[i] - self.wallThickness / 2) & (
                        rAtSpokeDir < rMeasure[i] + self.wallThickness / 2))
            xEmptyIdx, yEmptyIdx = self.convertRealXYToMapIdx(x + xAtSpokeDir[emptyIdx], y + yAtSpokeDir[emptyIdx])
            xOccupiedIdx, yOccupiedIdx = self.convertRealXYToMapIdx(x + xAtSpokeDir[occupiedIdx],
                                                                    y + yAtSpokeDir[occupiedIdx])
            if update:
                self.checkAndExapndOG(x + xAtSpokeDir[occupiedIdx], y + yAtSpokeDir[occupiedIdx])
                if len(emptyIdx) != 0:
                    self.total[yEmptyIdx, xEmptyIdx] += 1
                if len(occupiedIdx) != 0:
                    self.visited[yOccupiedIdx, xOccupiedIdx] += 2
                    self.total[yOccupiedIdx, xOccupiedIdx] += 2
            else:
                emptyXList.extend(x + xAtSpokeDir[emptyIdx])
                emptyYList.extend(y + yAtSpokeDir[emptyIdx])
                occupiedXList.extend(x + xAtSpokeDir[occupiedIdx])
                occupiedYList.extend(y + yAtSpokeDir[occupiedIdx])
        if not update:
            return np.asarray(emptyXList), np.asarray(emptyYList), np.asarray(occupiedXList), np.asarray(occupiedYList)

    def plot(self, xRange=None, yRange=None, plotThreshold=True):
        if xRange is None or xRange[0] < self.mapXLim[0] or xRange[1] > self.mapXLim[1]:
            xRange = self.mapXLim
        if yRange is None or yRange[0] < self.mapYLim[0] or yRange[1] > self.mapYLim[1]:
            yRange = self.mapYLim
        ogMap = self.visited / self.total
        xIdx, yIdx = self.convertRealXYToMapIdx(xRange, yRange)
        ogMap = ogMap[yIdx[0]: yIdx[1], xIdx[0]: xIdx[1]]
        ogMap = np.flipud(1 - ogMap)
        plt.imshow(ogMap, cmap='gray', extent=[xRange[0], xRange[1], yRange[0], yRange[1]])
        plt.show()
        if plotThreshold:
            ogMap = ogMap >= 0.5
            plt.matshow(ogMap, cmap='gray', extent=[xRange[0], xRange[1], yRange[0], yRange[1]])
            plt.show()
