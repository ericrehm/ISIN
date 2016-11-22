//
//  ISINGrid.cpp
//  Teledetection2
//
//  Created by Eric Rehm on 2014-04-09.
//  Copyright (c) 2014 Takuvik. All rights reserved.
//
//  Based on:
//      Java package org.esa.beam.dataio.merisl3  ISINGrid.java,
//      Copyright (C) 2010 Brockmann Consult GmbH (info@brockmann-consult.de)
//
//  Note that the original Brockmann Java code had errors which are corrected here.
//

#include "Arctic.h"
#include "ISINGridAPI.h"
#include "ISINGrid.h"
#include "math.h"
#include "string.h"
using namespace std;


//
// ISINGrid Class definition
//

class ISINGrid {
    
private:
    int rowCount;
    double binSize;
    double deltaLat;
    double *lats;
    double *deltaLons;
    int *rowLength;
    int *binOffsets;
    int totalBinCount;
    int binOffset45N;
    
public:
    
    /**
     * Earth equatorial radius in km.
     */
    static const double RE;
    /**
     * Default number of latitude rows.
     */
    static const int DEFAULT_ROW_COUNT; // number of latitude rows
    
    
    ISINGrid(int rowCount) {
        this->rowCount = rowCount;
        binSize = PI * RE / rowCount;
        deltaLat = 180.0 / rowCount;
        
        /* effective size (km and deg) of a cell in latitudinal direction */
        /* computes the number of cells in longitudinal direction for each rowIndex */
        /* and the bin effective longitudinal size */
        lats = new double [rowCount];
        deltaLons = new double [rowCount];
        rowLength = new int [rowCount];
        binOffsets = new int [rowCount];
        binOffset45N = -1;
        
        int binCount = 0;
        for (int i = 0; i < rowCount; i++) {
            lats[i] = -90.0 + (i + 0.5) * deltaLat;
            rowLength[i] = (int) (0.5 + 2.0 * rowCount * cos(lats[i]*PI/180.));
            deltaLons[i] = 360.0 / rowLength[i];
            binOffsets[i] = (i == 0) ? 0 : binOffsets[i - 1] + rowLength[i - 1];
            binCount += rowLength[i];
            if ((lats[i] >= 45.0) && (binOffset45N == -1)) {
                binOffset45N = binOffsets[i];
            }
        }
        totalBinCount = binCount;
    }
    
    int getRowCount() {
        return rowCount;
    }
    
    double getBinSize() {
        return binSize;
    }
    
    double getDeltaLat() {
        return deltaLat;
    }
    
    int getTotalBinCount() {
        return totalBinCount;
    }
    
    double getLat(int rowIndex) {
        return lats[rowIndex];
    }
    
    double getDeltaLon(int rowIndex) {
        return deltaLons[rowIndex];
    }
    
    int getRowLength(int rowIndex) {
        return rowLength[rowIndex];
    }
    
    int getBinOffset(int rowIndex) {
        return binOffsets[rowIndex];
    }
    
    int getRowIndex(int binIndex) {
        int *binOffsets = this->binOffsets;
        int totalBinCount = this->totalBinCount;
        
        int iL = 0;
        int iU = rowCount - 1;
        int iM;
        
        if (binIndex < 0 || binIndex >= totalBinCount) {
            return -1;
        }
        
        // Binary search for rowIndex
        do {
            // todo - check if this is an optimization or not
            //            if (binIndex == binOffsets[iL]) {
            //                return iL;
            //            }
            //            if (binIndex == binOffsets[iU]) {
            //                return iU;
            //            }
            iM = (iL + iU - 1) / 2;
            if (binIndex < binOffsets[iM]) {
                iU = iM;
            } else if (binIndex >= binOffsets[iM + 1]) {
                iL = iM + 1;
            } else {
                return iM;
            }
        } while (iL != iU);
        
        return iL;
    }
    
    int getRowIndex(double lat) {
        int rowIndex = ((lat + 90.0)*rowCount/180);
        rowIndex = min(rowIndex, rowCount-1);
        return rowIndex;
    }
    
    
    GridPoint getGridPoint(int binIndex, GridPoint p) {
        int rowIndex = getRowIndex(binIndex);
        int colIndex = (rowIndex == -1) ? -1 : (binIndex - binOffsets[rowIndex]);
        p.lat = getLat(rowIndex);
        p.lon = ((360.0*(colIndex + 0.5)/getRowLength(rowIndex))-180.0);
        return p;
    }

    GridPoint getGridPoint45N(int binIndex45N, GridPoint p) {
        int binIndex = binIndex45N += binOffset45N;
        return getGridPoint(binIndex, p);
    }
    
    /**
     * Gets the zero-based column index in the global ISIN grid for the given zero-based row index and longitude.
     *
     * @param rowIndex the zero-based row index in the range 0...{@link #getRowCount()}-1
     * @param lon      the longitude in the range 0...360 degree or -180...+180
     *
     * @return the zero-based column index only if the longitude is in the range 0...360 degree, otherwise undefined
     *
     */
    int getColIndex(int rowIndex, double lon) {
        if (lon >= 180) lon -= 360;
        return (int) ((lon + 180.0) * getRowLength(rowIndex)/360.0);
    }
    
    /**
     * Gets the zero-based bin index in the global ISIN grid for the given zero-based row index and longitude.
     *
     * @param rowIndex the zero-based row index in the range 0...{@link #getRowCount()}-1
     * @param lon      the longitude in the range 0...360 degree or -180...+180
     *
     * @return the zero-based bin index only if the longitude is in the range 0...360 degree, otherwise undefined
     *
     */
    int getBinIndex(int rowIndex, double lon) {
        int colIndex = getColIndex(rowIndex, lon);
        return getBinOffset(rowIndex) + colIndex;
    }
    
    /**
     * Gets the zero-based bin index in the global ISIN grid for the given zero-based row index and longitude.
     *
     * @param lat      the latitude in the range of -90 ... 90
     * @param lon      the longitude in the range 0...360 degree or -180...+180
     *
     * @return the zero-based bin index 
     *
     */
    int getBinIndex(double lat, double lon) {
        int rowIndex = getRowIndex(lat);
        return getBinIndex(rowIndex, lon);
    }

    int getBinIndex(GridPoint p) {
        return getBinIndex(p.lat, p.lon);
    }
    
    /**
     * Gets the zero-based bin index in the Arctic ISIN grid for the given zero-based row index and longitude.
     *
     * @param rowIndex the zero-based row index in the range 0...{@link #getRowCount()}-1
     * @param lon      the longitude in the range 0...360 degree or -180...+180
     *
     * @return the zero-based bin index only if the longitude is in the range 0...360 degree, otherwise undefined
     *
     */
    int getBinIndex45N(int rowIndex, double lon) {
        return getBinIndex(rowIndex, lon) - binOffset45N;
    }
    
    /**
     * Gets the zero-based bin index in the Arctic ISIN grid for the given zero-based row index and longitude.
     *
     * @param lat      the latitude in the range of -90 ... 90
     * @param lon      the longitude in the range 0...360 degree or -180...+180
     *
     * @return the zero-based bin index
     *
     */
    int getBinIndex45N(double lat, double lon) {
        return getBinIndex(lat, lon) - binOffset45N;
    }
    
    int getBinIndex45N(GridPoint p) {
        return getBinIndex(p.lat, p.lon) - binOffset45N;
    }


#ifdef XXXX
    /**
     * Detects the row count from the product name.
     *
     * @param productName the name of the L3 product
     *
     * @return the row count
     *
     * From original BEAM MERIS Java code ... not yet ported.
     */
    static int detectRowCount(string productName) {
        Pattern p = Pattern.compile(".*_(\\d{4})x(\\d{4})_.*");
        Matcher m = p.matcher(productName);
        if (m.matches() && m.groupCount() == 2) {
            String binSize1 = m.group(1);
            String binSize2 = m.group(2);
            if (binSize1.equals(binSize2)) {
                int binSize = Integer.parseInt(binSize1);
                return (int) round(PI * RE * 1000/ binSize);
            }
        }
        return DEFAULT_ROW_COUNT;
    }
#endif
};

const double ISINGrid::RE = 6378.137;  // Radius of the earth
const int ISINGrid::DEFAULT_ROW_COUNT = 2160; // number of latitude rows

//
// External C API to ISINGrid class
//

API void *ISINGridCreate(int rowCount) {
  return (void *) new ISINGrid(rowCount);
};

API int ISINGridGetBinIndex(void *grid, GridPoint p) {
  return ((ISINGrid *)grid)->getBinIndex(p);
};

API int  ISINGridGetBinIndex45N(void *grid, GridPoint p) {
  return ((ISINGrid *)grid)->getBinIndex45N(p);
};

API GridPoint ISINGridGetGridPoint(void *grid, int binIndex, GridPoint p) {
  return ((ISINGrid *)grid)->getGridPoint(binIndex, p);
};

API GridPoint ISINGridGetGridPoint45N(void *grid, int binIndex, GridPoint p) {
  return ((ISINGrid *)grid)->getGridPoint45N(binIndex, p);
};

API void ISINGridDelete(void *grid) {
    delete (ISINGrid *)grid;
};



//
// Test routine
//
// Rename ISINGrid_Test to main
//
int test_main(int argc, char **argv) {
    
    int irow, rowCount, numbin[SZLAT_SEAWIFS], basebin[SZLAT_SEAWIFS], totbins;
    float latbin[SZLAT_SEAWIFS];
    
    float lat, lon, lon2, rlat, rlon;
    int jrow, icol, idx;
    
    cout << "***** ISINGrid ******  \n";

    lat = 60;
    lon = -20;
    
    if (argc > 2) {
        lat = atof(argv[1]);
        lon = atof(argv[2]);
    }
    printf("%8.4f N, %9.4f E\n", lat, lon);
    
    rowCount = SZLAT_SEAWIFS;
    
    //
    // Reference code based on FORTRAN Code published in NASA Tech Memo 104566, p. 64
    // Vol 32, Level-3 SeaWiFS Data Products: Spatial and Temporal Binning Algorithms
    // Appendix A: Equal-area gridding scheme for SeaWiFS Binned Data
    //
    printf("***** NASA Code ******\n");

    cout << "rowCount = " << rowCount << endl;

    basebin[0] = 0;
    for (irow = 0; irow < rowCount; irow++) {
        latbin[irow] = ((irow + 0.5) *180.0/rowCount) - 90.0;
        numbin[irow] = (int) (2*rowCount*cos(latbin[irow]*PI/180) + 0.5);
        if (irow > 0) basebin[irow] = basebin[irow-1] + numbin[irow-1];
        //        printf("(%d) : %8.4f  %7d %4d\n", irow, latbin[irow], basebin[irow], numbin[irow]);
    }
    totbins = basebin[rowCount-1] + numbin[rowCount-1];
    printf("totbins  = %d\n", totbins);

    // lat, lon to bin number
    jrow = (int) ((lat + 90.0)*rowCount/180);
    jrow = min(jrow, rowCount-1);
    lon2 = (lon > 180 ? lon -= 360 : lon);
    icol = (int) ((lon2 + 180.0) * numbin[jrow]/360);
    icol = min(icol, numbin[jrow]-1);
    idx  = basebin[jrow] + icol;
    
    printf("(jrow, icol) = (%d, %d), idx = %d\n", jrow, icol, idx);
    printf("(%d) : %8.4f  %7d %4d\n", jrow, latbin[jrow], basebin[jrow], numbin[jrow]);
    printf("%8.4F N, %9.4f E ==> %d\n", lat, lon, idx - basebin[SZLAT_45_SEAWIFS]);

    rlat = latbin[jrow];
    rlon = (float)((360.0*(icol + 0.5)/numbin[jrow])-180.0);
    printf("%8.4f N, %9.4f E\n", rlat, rlon);
    
    // Now compare with this ISINGrid implementation
    
    printf("***** ISINGrid Code ******\n");
    
    ISINGrid *igrid = new ISINGrid((int) SZLAT_SEAWIFS);
    
    cout << "rowCount = " << igrid->getRowCount() << endl;
    cout << "totbins  = " << igrid->getTotalBinCount() << endl;
    cout << "deltaLat =     " << igrid->getDeltaLat() << endl;
    cout << "binSize  =     " << igrid->getBinSize() << " km" << endl;

    jrow = igrid->getRowIndex(lat);
    int jrow2 = igrid->getRowIndex(idx);
    icol = igrid->getColIndex(jrow2, (double)lon);

    // row, lon or lat,lon to bin number
    idx = igrid->getBinIndex(jrow, lon);
    int idx2 = igrid->getBinIndex(lat, lon);
    printf("(jrow, icol) = (%d, %d), idx = %d:%d, %d:%d\n", jrow, icol, idx, idx2, idx - basebin[SZLAT_45_SEAWIFS], igrid->getBinIndex45N(lat, lon));
    printf("(%d:%d) : %8.4f  %7d %4d\n",
           jrow, jrow2, igrid->getLat(jrow), igrid->getBinOffset(jrow), igrid->getRowLength(jrow));
    printf("%8.4f N, %9.4f deltaLon\n", igrid->getLat(jrow), igrid->getDeltaLon(jrow));
    
    // Test ISINGrip GridPoint routines
    GridPoint p;
    p.lat = 0;
    p.lon = 0;
    p = igrid->getGridPoint(idx, p);
    printf("%8.4f N, %9.4f E\n", p.lat, p.lon);
    
    idx2 = igrid->getBinIndex45N(lat, lon);      // idx2 is the index in our Arctic Grid
    printf("%8.4f N, %9.4f E ==>    idx2 = %d\n", lat, lon, idx2);
    p.lat = 0;
    p.lon = 0;
    p = igrid->getGridPoint45N(idx2, p);
    printf("%8.4f N, %9.4f E\n", p.lat, p.lon);

    printf("bye\n");
    exit(0);
    
}
