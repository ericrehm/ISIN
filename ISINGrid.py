#!/usr/bin/env python
"""GCmatchup6.py: Matchup BIO chl data with CCI pixels in NetCDF4 format."""

# imports
import numpy as np

class ISINGrid:

    # Class variables, e.g. ISINGrid.SEAWIFS 
    RE = 6378.137    # Radius of the earth
    SEAWIFS = 2160   # Row count for SeaWiFS
    MODIS   = 4320   # Row count for MODIS, GlobColour, CCI, VIIRS

    # ISINGrid constructor
    def __init__(self, rc):    

        # Instance variables
        self.rowCount = 0
        self.deltaLat = 0
        self.totalBinCount = 0

        if rc > 0:
            self.rowCount = rc
        else:
            self.rowCount = self.MODIS  # 4 km MODIS grid

        self.binSize = np.pi * self.RE / self.rowCount;    
        self.deltaLat = 180.0 / self.rowCount

        self.lats       = np.zeros(self.rowCount)
        self.rowLength  = np.zeros(self.rowCount)
        self.deltaLons  = np.zeros(self.rowCount)
        self.binOffsets = np.zeros(self.rowCount)
        
        #         /* effective size (km and deg) of a cell in latitudinal direction */
        #         /* computes the number of cells in longitudinal direction for each rowIndex */
        #         /* and the bin effective longitudinal size */
        #         lats = new double [rowCount];
        #         deltaLons = new double [rowCount];
        #         rowLength = new int [rowCount];
        #         binOffsets = new int [rowCount];

        # Constructor code
        binCount = 0
        for i in range(0,self.rowCount):
            self.lats[i] = -90.0 + (i + 0.5) * self.deltaLat
            self.rowLength[i] =  np.floor(0.5 + 2.0 * self.rowCount * np.cos(self.lats[i]*np.pi/180.))
            self.deltaLons[i] = 360.0 / self.rowLength[i]
            if i == 0:
                self.binOffsets[i] = 0
            else:
                self.binOffsets[i] = self.binOffsets[i - 1] + self.rowLength[i - 1]
            binCount = binCount + self.rowLength[i]
            
        self.totalBinCount = binCount

    # Class methods
        
    def getLat(self, rowIndex):
        return self.lats[rowIndex]

    def getDeltaLon(self, rowIndex):
        return self.deltaLons[rowIndex]

    def getRowLength(self, rowIndex):
        return self.rowLength[rowIndex]
    
    def getBinOffset(self, rowIndex):
        return self.binOffsets[rowIndex]
    
    def getRowIndex(self, lat):
        rowIndex = np.floor((lat + 90.0)*(self.rowCount)/180)
        rowIndex = min(rowIndex, self.rowCount)
        return rowIndex
    
    def getRowIndex2(self, binIndex):
        g = np.where(binIndex >= self.binOffsets)[0]
        rowIndex = g[-1]
        return rowIndex
    
    #     /**
    #      * Gets the zero-based column index in the global ISIN grid for the given zero-based row index and longitude.
    #      *
    #      * @param rowIndex the zero-based row index in the range 0...{@link #getRowCount()}-1
    #      * @param lon      the longitude in the range 0...360 degree or -180...+180
    #      *
    #      * @return the zero-based column index only if the longitude is in the range 0...360 degree, otherwise undefined
    #      *
    #      */
    def getColIndex(self, rowIndex, lon):
        if (lon >= 180):
            lon = lon - 360
        colIndex =  np.floor((lon + 180.0) * self.getRowLength(rowIndex)/360.0)
        return colIndex
        
    def getBinIndexLon(self, rowIndex, lon):
        colIndex = self.getColIndex(rowIndex, lon)
        binIndex = self.getBinOffset(rowIndex) + colIndex
        return binIndex
    
    def getBinIndex(self, lat, lon):
        rowIndex = self.getRowIndex(lat)
        binIndex = self.getBinIndexLon(rowIndex, lon)
        return binIndex

    def getGridPoint(self, binIndex):
        rowIndex = self.getRowIndex2(binIndex)
        if rowIndex == -1:
            colIndex = -1
        else:
            colIndex = binIndex - self.binOffsets[rowIndex]
        p = np.zeros(2)
        p[0] = self.getLat(rowIndex)
        p[1] = ((360.0*(colIndex + 0.5)/self.rowLength[rowIndex])-180.0);
        return p
        
if __name__ == '__main__':

    lat = 60
    lon = -20
    idx = 5543444
    print "lat, lon = ", lat, lon
    
    igrid = ISINGrid(ISINGrid.SEAWIFS)
    # igrid = ISINGrid(ISINGrid.MODIS)
    print "rowCount = %d (%d)" % (igrid.rowCount, 2160)
    print "totbins  = %d (%d)" % (igrid.totalBinCount, 5940422)
    print "deltaLat = %f (%f)" % (igrid.deltaLat, 0.0833333)
    print "binSize  = %f (%f) km" % ( igrid.binSize, 9.27662)

    jrow  = igrid.getRowIndex(lat)
    jrow2 = igrid.getRowIndex2(idx)
    icol  = igrid.getColIndex(jrow2, lon)
    print "(jrow, icol) = (%d, %d) ((%d, %d)), idx = %d (%d)" % (jrow, icol, 1800, 958, idx, 5543444)
    print "jrow2 = %d" % jrow2

    idx2 = igrid.getBinIndex(lat, lon)
    print "idx2 = %d (%d) from (lat, lon) = (%f, %f)" % (idx2, idx, lat, lon)

    p = igrid.getGridPoint(idx2)
    print "%f N, %f E  (60.0417 N, -20.0278 E) from idx2 = %d" % (p[0], p[1], idx2)
