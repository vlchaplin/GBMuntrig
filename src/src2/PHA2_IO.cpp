/*
 *  PHA2_IO.cpp
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#ifndef PHA_2_IODEFINED
#define PHA_2_IODEFINED

#include "PHA2_IO.hh"
#include <typeinfo>

//#define DEBUG_PHA_SEARCH

int PHA2Reader::binSortedIntersection(long& rowMIN, long& rowMAX, double& searchTime, FxbReadStat& stat ) {

	if ( stat.status || rowMIN > rowMAX || rowMIN < 1 ) return -1;

	long median = ( rowMIN+rowMAX ) / 2;
	double mtime;
	
	if ( rowMAX - rowMIN == 0 ) {
		#ifdef DEBUG_PHA_SEARCH
		cout << "equal " << rowMIN << endl;
		#endif
		return 0;
	}
	if ( rowMAX - rowMIN <= 4 ) {
		double minRow, residue, minResidue;
		rowMAX++;
		
		ffgcv( stat.fptr, TDOUBLE, search_Col, rowMAX, 1, 1, 0, &mtime, NULL, &(stat.status) );
		
		minResidue = abs(searchTime - mtime);
		minRow = rowMAX;
		
		while ( rowMIN <= rowMAX && stat.status == 0 )
		{
			ffgcv( stat.fptr, TDOUBLE, search_Col, --rowMAX, 1, 1, 0, &mtime, NULL, &(stat.status) );
			
			residue = abs(searchTime - mtime);
			
			if ( minResidue > residue ) {
				minResidue = residue;
				minRow = rowMAX;
			}
		}
		
		rowMIN = minRow;
		rowMAX = minRow;

		return 0;
	}
	
	ffgcv( stat.fptr, TDOUBLE, search_Col, median, 1, 1, 0, &mtime, NULL, &(stat.status) );
	bool branchLower = mtime > searchTime;
	bool branchUpper = mtime < searchTime;
	bool equal = abs(searchTime - mtime) < 0.064;
	
	if (equal) {
		while ( mtime > searchTime && stat.status == 0 )
		{
			ffgcv( stat.fptr, TDOUBLE, search_Col, --median, 1, 1, 0, &mtime, NULL, &(stat.status) );
		}
		#ifdef DEBUG_PHA_SEARCH
		char temp[100];
		sprintf( temp, "{ %10ld - %10ld } -> %10ld ,  %16f, %16f ", rowMIN, rowMAX, median, mtime, searchTime - mtime); 
		cout << "equal " << temp << endl;
		#endif
		rowMIN=median;
		return 0;
	} 
	else if (branchUpper) {
		#ifdef DEBUG_PHA_SEARCH
		char temp[100];
		sprintf( temp, "{ %10ld - %10ld } -> %10ld ,  %16f, %16f ", rowMIN, rowMAX, median, mtime, searchTime - mtime); 
		cout << "upper " << temp << endl;
		#endif
		rowMIN = median+1;
		return binSortedIntersection( rowMIN, rowMAX, searchTime, stat );
	}
	else if (branchLower) {	
		#ifdef DEBUG_PHA_SEARCH
		char temp[100];
		sprintf( temp, "{ %10ld - %10ld } -> %10ld ,  %16f, %16f ", rowMIN, rowMAX, median, mtime, searchTime - mtime ); 
		cout << "lower " << temp << endl;
		#endif
		rowMAX = median;
		return binSortedIntersection( rowMIN, rowMAX, searchTime, stat );
	} 
	
	return -1;
}

int PHA2Reader::findRow( double time, long& row ) {

	FxbReadStat * fxbCntl = this->getFxb();
	if ( ! fxbCntl->okForOps() ) return fxbCntl->status;


	double tstart,tstop;
	long R0,R1;
	int rc;
	R0 = 1;
	R1 = fxbCntl->n_rows;
	
	ffgcv(fxbCntl->fptr, TDOUBLE, search_Col, R0, 1, 1, 0, &tstart, NULL, &(fxbCntl->status) );
	ffgcv(fxbCntl->fptr, TDOUBLE, search_Col, R1, 1, 1, 0, &tstop, NULL,  &(fxbCntl->status) );

	
	if ( time < tstart || time > tstop ) 
		rc = 1;
	else
		rc = this->binSortedIntersection( R0,R1,time, *fxbCntl );
	
	if ( ! rc ) row = R0;
	
	if ( fxbCntl->status ) fxbCntl->errmsg();
	
	return rc;
};

int PHA2Reader::getSpectra( double * dataTable, double * exposure, long startrow, long numrows )
{

	FxbReadStat * fxbCntl = this->getFxb();
	
	if ( ! fxbCntl->okForOps() ) return fxbCntl->status;
	
	int nChan;
	if ( fxbCntl->move2hdu( (char*)"EBOUNDS" ) ) return fxbCntl->status;
	nChan = fxbCntl->n_rows;
	
	if ( fxbCntl->move2hdu( (char*)"SPECTRUM" ) ) return fxbCntl->status;

	//fftscl( fxbCntl->fptr, 1, 1.0, 2*32768.0, &(fxbCntl->status) );
	//ffflus( fxbCntl->fptr, &(fxbCntl->status) );
	int z=0;
//	cout << "here" << endl;
//	ffuky( fxbCntl->fptr, TINT, (char*)"TZERO1", &z, NULL, &(fxbCntl->status) );
//	ffflus( fxbCntl->fptr, &(fxbCntl->status) );
//	size_t index=0;
	long max_rows = fxbCntl->n_rows;
	long n_eff_rows = fxbCntl->n_eff_rows;
	long q=0;
	long thisrow = startrow;
	int dataCol,exp_c;
	if ( startrow > max_rows ) return -1;
	
    ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"EXPOSURE",    &exp_c, &(fxbCntl->status) );
	if ( ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"COUNTS",  &dataCol, &(fxbCntl->status) ) ) 
	{
		fxbCntl->status=0;
		if ( ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"RATE",    &dataCol, &(fxbCntl->status) ) )
			return fxbCntl->status;
	}
	
	if ( n_eff_rows > numrows ) n_eff_rows = numrows;
	
	long lastrow = startrow+numrows-1;
	
	if ( lastrow  > max_rows ) lastrow = max_rows;

	//cout << n_eff_rows << endl;
	//n_eff_rows=1;
	//dataCol=1;

	while ( fxbCntl->status == 0 && thisrow <= lastrow ) {
	
		if ( (thisrow + n_eff_rows - 1) > lastrow ) n_eff_rows = lastrow - thisrow + 1;
		
		ffgcv(fxbCntl->fptr, TDOUBLE, dataCol, thisrow, 1, n_eff_rows*nChan, 0, (dataTable + nChan*q), NULL, &(fxbCntl->status) );
        if ( exposure != NULL )
            ffgcv(fxbCntl->fptr, TDOUBLE, exp_c, thisrow, 1, n_eff_rows, 0, (exposure + q), NULL, &(fxbCntl->status) );
		
		thisrow+=n_eff_rows;
		q+=n_eff_rows;

	}
	
	if ( fxbCntl->status ) fxbCntl->errmsg();
	
	return (fxbCntl->status);
}

int PHA2Reader::getTimes( double * tstart, double * tstop, long startrow, long numrows )
{

	FxbReadStat * fxbCntl = this->getFxb();
	
	if ( ! fxbCntl->okForOps() ) return fxbCntl->status;

//	size_t index=0;
	long n_eff_rows = fxbCntl->n_eff_rows;
	long thisrow = startrow;
	long max_rows = fxbCntl->n_rows;
	long lastrow = startrow+numrows-1;
	if ( startrow > max_rows ) return -1;
	if ( lastrow  > max_rows ) lastrow = max_rows;
	
	long q=0;
	int exp_c, tstart_c, tstop_c;
	ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"EXPOSURE",    &exp_c, &(fxbCntl->status) );
	ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"TIME",  &tstart_c, &(fxbCntl->status) );
	ffgcno(fxbCntl->fptr, CASEINSEN, (char *)"ENDTIME",    &tstop_c, &(fxbCntl->status) );
	
	if ( n_eff_rows > numrows ) n_eff_rows = numrows;

	while ( fxbCntl->status == 0 && thisrow <= lastrow ) {
	
		if ( (thisrow + n_eff_rows - 1) > lastrow ) n_eff_rows = lastrow - thisrow + 1;
		
		//cout << thisrow << ", " << lastrow << ", " << n_eff_rows << ", " << endl;
		
		if ( tstart!= NULL) ffgcv(fxbCntl->fptr, TDOUBLE, tstart_c, thisrow, 1, n_eff_rows, 0, tstart + q, NULL, &(fxbCntl->status));
		if ( tstop != NULL) ffgcv(fxbCntl->fptr, TDOUBLE, tstop_c , thisrow, 1, n_eff_rows, 0, tstop  + q, NULL, &(fxbCntl->status));
		
		thisrow+=n_eff_rows;
		q+=n_eff_rows;
	}
	
	return (fxbCntl->status);
}



#endif