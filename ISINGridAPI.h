//
//  ISINGridAPI.h 
//  Teledetection2
//
//  Created by Eric Rehm on 2014-04-11.
//  Copyright (c) 2013 Takuvik. All rights reserved.
//

#ifndef Teledetection2_ISINGridAPI_h
#define Teledetection2_ISINGridAPI_h

#ifdef __cplusplus
#   define API extern "C"
#else
#   define API
#endif

// ISINGrid C-language API

typedef struct GridPoint {
    double lon;
    double lat;
} GridPoint;


API void *ISINGridCreate(int rowCount);
API int ISINGridGetBinIndex(void *grid, GridPoint p);
API int  ISINGridGetBinIndex45N(void *grid, GridPoint p);
API GridPoint ISINGridGetGridPoint(void *grid, int binIndex, GridPoint p);
API GridPoint ISINGridGetGridPoint45N(void *grid, int binIndex, GridPoint p);
API void ISINGridDelete(void *grid);

#endif
