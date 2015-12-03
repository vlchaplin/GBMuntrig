/*
 *  Gbm_PositionIO.h
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */


#ifndef GBMPOSITION_H
#define GBMPOSITION_H

#include <math.h>
#include <string.h>
#include <iostream>
#include <iomanip>

#include <vector>
#include <iterator>
#include <algorithm>
extern "C" {
#include "fitsio.h"
};
#include "PHA_IO.hh"
#include "gbm_geometry.hh"
#include "algorithms_t.h"

using namespace std;

class GbmPosReader : public PHAReader {

	public:
	
	GbmPosReader(){};
	
	virtual int getTimes( double *& tstarts, double *& tstops, long * ntimes, bool alloc=true ) {return -1;};
	virtual int getTimes( vector<double>& tstart, vector<double>& tstop, long startRow=0, long endRow=0 ) {return -1;};
	
	virtual int getScAttPos( double * times, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows = 1 ) {return -1;};
	virtual int getScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows = 1 ) {return -1;};
	virtual int getScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow=0, long endRow=0 ) {return -1;};
	long findRow( double time, vector<double> * starts = NULL, vector<double> * ends = NULL );

};

class GbmPosWriter : public PHAWriter {

	public:
	
	GbmPosWriter(){};
	
	virtual int writeScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long rownum ) {return -1;};
	virtual int writeScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow=0, long endRow=0 ) {return -1;};

};


class TrigdatReader : public GbmPosReader {

	private:
	
	template<typename T>
	void bindKeysToMem( OGIP_PHA_misc<T>& keys ) {
		
		FitsKeyLocation links[] = { {"TRIGTIME",TDOUBLE,1,(void*)&keys.tzero},
									{"TSTART",TDOUBLE, 1,(void*)&keys.tmin}, 
									{"TSTOP",TDOUBLE,  1,(void*)&keys.tmax}
									};
		
		this->fileToMemLink.assign( links, links+5 );
	};


	


	public:
	TrigdatReader() {};
	
	
	virtual int getTimes( double *& tstarts, double *& tstops, long * ntimes, bool alloc=true  );
	virtual int getTimes( vector<double>& tstart, vector<double>& tstop, long startRow=0, long endRow=0 );
	virtual int getScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows = 1 );
	virtual int getScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow=0, long endRow=0 );
	//virtual long findRow( double time );
	//int getScAttPos(long firsRow, long lastRow);


};

class TrigdatWriter : public GbmPosWriter {

	private:

/*	
	template<typename T>
	void bindKeysToMem( OGIP_PHA_misc<T>& keys ) {
		
		FitsKeyLocation links[] = { {"TRIGTIME",TDOUBLE,1,(void*)&keys.tzero},
									{"TSTART",TDOUBLE, 1,(void*)&keys.tmin}, 
									{"TSTOP",TDOUBLE,  1,(void*)&keys.tmax}
									};
		
		this->fileToMemLink.assign( links, links+5 );
	};*/

	public:
	TrigdatWriter() {};
	
	int writeRateData( long detNum, double * rates, int nChan, long firstRow=1, long numRows=1 );
	int writeTimes( double * tstart, double * tstop, long firstRow=1, long numRows=1 );
	virtual int writeScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long startRow=1, long numRows=0 );
	virtual int writeScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow=1, long numRows=0 ) {return -1;};

};

class PoshistReader : public GbmPosReader {

	private:
	
	template<typename T>
	void bindKeysToMem( OGIP_PHA_misc<T>& keys ) {
		
		FitsKeyLocation links[] = { {"TRIGTIME",TDOUBLE,1,(void*)&keys.tzero},
									{"TSTART",TDOUBLE, 1,(void*)&keys.tmin}, 
									{"TSTOP",TDOUBLE,  1,(void*)&keys.tmax}
									};
		
		this->fileToMemLink.assign( links, links+5 );
	};


	string attExtName;


	public:
	PoshistReader() { attExtName="GLAST POS HIST"; };
	
	virtual int getTimes( vector<double>& tstart, vector<double>& tstop, long startRow=0, long endRow=0 );
	virtual int getScAttPos( double * tstart, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows=1 );
	virtual int getScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows=1 ) {
		
		int status = this->getScAttPos( tstart, quats, pos, rownum, numrows );
		if ( status ) return status;
		
		for(long i=0; i<numrows; i++) tstop[i] = tstart[i];
		return status;
	};
	virtual int getScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow=0, long endRow=0 );
	//virtual long findRow( double time );

};


class NewPosReader {

	private:
	
	NewPosReader() {};
	
	public:
	
	static GbmPosReader * type(string type) {
		if (type == "trigdat") {
			return new TrigdatReader;
		} else if (type == "poshist") {
			return new PoshistReader;
		} else {
			return NULL;
		}
	};


};




#endif