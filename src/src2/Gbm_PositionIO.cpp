/*
 *  Gbm_PositionIO.cpp
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#include "Gbm_PositionIO.h"

long GbmPosReader::findRow( double time, vector<double> * sTimes, vector<double> * eTimes ) {
	
	long row = -1;
	
	if ( sTimes == NULL && eTimes == NULL ) {
		sTimes = new vector<double>;
		eTimes = new vector<double>;
		if ( this->getTimes( *sTimes, *eTimes ) == 0 ) {
			row = findNearestElement( *sTimes, time );		
		}
		delete sTimes;
		delete eTimes;
	} else if ( sTimes->size()*eTimes->size() == 0 ) {
	
		if ( this->getTimes( *sTimes, *eTimes ) == 0 ) {
			row = findNearestElement( *sTimes, time );		
		}
		
	} else if ( sTimes->size() > 0 ) {
		row = findNearestElement( *sTimes, time );	
	}
	
	
	
	return row;
};


int TrigdatReader::getTimes( double *& tstarts, double *& tstops, long * ntimes, bool alloc )
{
	this->setOpRows(0);
	
	if ( this->stat.move2hdu((char*)"EVNTRATE" ) ) return this->stat.status;
	
    
    
	*ntimes = stat.n_rows;
	
    if ( alloc ) {
        tstarts = new double[stat.n_rows];
        tstops = new double[stat.n_rows];
        cout << "new() alloc in " << __FILE__ << ":" << __LINE__ << endl;
    }
	
	//double * ptr1 = &(tstarts[0]);
	//double * ptr2 = &(tstops[0]);
	long i=1;
    
    size_t q=0;
    
    long neffrows = stat.n_eff_rows;
    long nrows = *ntimes;
    if ( neffrows > nrows ) neffrows = nrows;
    neffrows=1;
	
	while( stat.status == 0 && i <= nrows ) {
	
		if ( (i + neffrows - 1) > nrows) neffrows = nrows - i + 1;
	
		if ( tstarts != NULL) ffgcv(stat.fptr, TDOUBLE, 1, i, 1, neffrows,0, (tstarts+q),NULL,&(this->stat.status) );
		if ( tstops  != NULL) ffgcv(stat.fptr, TDOUBLE, 2, i, 1, neffrows,0, (tstops+q), NULL,&(this->stat.status) );
	

	//	cout << i << " " << tstarts[q] << "," << tstops[q] << endl;
		
		if (stat.status != 0) {
            this->errmsg();
            break;
        }
        
        this->setOpRows(i);
        q+=neffrows;
		i+=neffrows;

	};
	
	//tstarts = ptr1;
	//tstops = ptr2;
	
	return stat.status;

}

int TrigdatReader::getTimes( vector<double>& tstart, vector<double>& tstop, long startRow, long endRow )
{
	this->setOpRows(0);
	
	if ( this->stat.move2hdu((char*)"EVNTRATE" ) ) return this->stat.status;
	
	long stride = stat.n_eff_rows;
	long n_rows = stat.n_rows;
	
	if ( n_rows <=0 ) return -1;
	
	if ( endRow > 0 && startRow > 0 ) {
		n_rows = endRow - startRow + 1;
	} else {
		endRow = n_rows;
		startRow = 1;
	}
	
	tstart.reserve(n_rows);
	tstop.reserve(n_rows);
	
	double * tstarts = new double[stride];
	double * tstops = new double[stride];
	
	long i=startRow;
	
	while( stat.status == 0 && i <= endRow ) {
	
		if ( (i + stride - 1) > endRow) stride = endRow - i + 1;
	
		ffgcv(stat.fptr, TDOUBLE, 1, i, 1, stride,0, tstarts,NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 2, i, 1, stride,0, tstops, NULL,&(this->stat.status) );
	
		for (int k=0;k<stride;k++) {
			tstart.push_back( tstarts[k] );
			tstop.push_back( tstops[k] );
		}
		
		if (stat.status != 0) break;
		
		this->setOpRows(i);
        i+=stride;
	}
	
	delete [] tstarts;
	delete [] tstops;
	
	return stat.status;

}

int TrigdatReader::getScAttPos( double * tstart, double * tstop, geom::gquat * quat, geom::x3 * pos, long rownum, long numrows )
{
	this->setOpRows(0);
	
	if ( this->stat.move2hdu((char*)"EVNTRATE" ) ) return this->stat.status;
	if (stat.status != 0) return this->stat.status;
	
	double scattit[4];
	double ecipos[3];
	
	long i=rownum;
	
	if ( i > stat.n_rows  || i < 1 ) return -1;
	
	ffgcv(stat.fptr, TDOUBLE, 1, i, 1, 1,0, tstart,NULL,&(this->stat.status) );
	ffgcv(stat.fptr, TDOUBLE, 2, i, 1, 1,0, tstop, NULL,&(this->stat.status) );
	ffgcv(stat.fptr, TDOUBLE, 3, i, 1, 4,0, scattit, NULL,&(this->stat.status) );
	ffgcv(stat.fptr, TDOUBLE, 4, i, 1, 3,0, ecipos, NULL,&(this->stat.status) );
	
	quat->set(scattit);
	pos->set(ecipos);
	
	this->setOpRows(1);
	
	return stat.status;
}

int TrigdatReader::getScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow, long endRow )
{

	this->setOpRows(0);

	if ( this->stat.move2hdu((char*)"EVNTRATE" ) ) return this->stat.status;
	
	long stride = stat.n_eff_rows;
	long n_rows = stat.n_rows;
	
	if ( endRow > 0 && startRow > 0 ) {
		n_rows = endRow - startRow + 1;
	} else {
		endRow = n_rows;
		startRow = 1;
	}
	
	tstart.reserve(n_rows);
	tstop.reserve(n_rows);
	quats.reserve(n_rows);
	pos.reserve(n_rows);
	
	double * tstarts = new double[stride];
	double * tstops = new double[stride];
	double * scattit = new double[4*stride];
	double * ecipos3 = new double[3*stride];
	
	long i=startRow;
	
	while( stat.status == 0 && i <= endRow ) {
	
		if ( (i + stride - 1) > endRow) stride = endRow - i + 1;
	
		ffgcv(stat.fptr, TDOUBLE, 1, i, 1, stride,0, tstarts,NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 2, i, 1, stride,0, tstops, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 3, i, 1, 4*stride,0, scattit, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 4, i, 1, 3*stride,0, ecipos3, NULL,&(this->stat.status) );
		
	
		for (int k=0;k<stride;k++) {
			tstart.push_back( tstarts[k] );
			tstop.push_back( tstops[k] );
			
			quats.push_back( geom::gquat( scattit + k ) );
			pos.push_back( geom::x3( ecipos3 + k ) );
		}
				
		
		
		if (stat.status != 0) break;
		
		this->setOpRows(i);
        i+=stride;

	};
	delete [] tstarts;
	delete [] tstops;
	delete [] scattit;
	delete [] ecipos3;
	
	return stat.status;
}


int PoshistReader::getTimes( vector<double>& tstart, vector<double>& tstop, long startRow, long endRow )
{
	
	this->setOpRows(0);
	
	if ( this->stat.move2hdu((char*)attExtName.c_str()) ) {
		cout << " at " << __func__ << endl;
		return this->stat.status;
	}
	
	long stride = stat.n_eff_rows;
	long n_rows = stat.n_rows;
	
	if ( endRow > 0 && startRow > 0 ) {
		n_rows = endRow - startRow + 1;
	} else {
		endRow = n_rows;
		startRow = 1;
	}
	
	tstart.reserve(n_rows);
	tstop.reserve(n_rows);
	
	double * tstarts = new double[stride];
	
	long i=startRow;
	
	while( stat.status == 0 && i <= endRow ) {
	
		if ( (i + stride - 1) > endRow) stride = endRow - i + 1;
	
		ffgcv(stat.fptr, TDOUBLE, 1, i, 1, stride,0, tstarts,NULL,&(this->stat.status) );
	
		for (int k=0;k<stride;k++) tstart.push_back( tstarts[k] );

		
		if (stat.status != 0) break;
		
		this->setOpRows(i);
        i+=stride;

	};
	
	std::copy( tstart.begin(), tstart.end()-1, tstop.begin() );
	//The last entry has tstart=tstop
	tstop.push_back( tstarts[stride-1] );
	
	
	delete [] tstarts;
	
	return stat.status;

}

int PoshistReader::getScAttPos( double * times, geom::gquat * quat, geom::x3 * pos, long rownum, long numrows )
{
	FxbReadStat * fstruct = this->getFxb();
	
	this->setOpRows(0);
	
	if (fstruct->status != 0) return fstruct->status;
	if (! fstruct->okForOps() ) return -1;
	
	if (numrows <= 0) { numrows = 1; }
	
	long row = rownum;
	long lastrow = row+numrows-1;
	long buffsz = fstruct->n_eff_rows;
	long buff_j;
	size_t k=0;


	if ( lastrow > fstruct->n_rows ) lastrow = fstruct->n_rows;
	if ( lastrow-row+1 > buffsz ) buffsz = lastrow-row+1;
	//if ( row+lastrow-1 < buffsz ) buffsz = row+lastrow-1;
	
	double * scattit = new double[4*buffsz];
	double * ecipos  = new double[3*buffsz];
		
	while ( row <= lastrow && fstruct->status == 0 ) {
		
		if ( row+buffsz-1 > lastrow ) buffsz = lastrow-row+1;
		
		//cout << "reading " << row << " + "<< row+buffsz << endl;
	
		ffgcv(fstruct->fptr, TDOUBLE, 1, row, 1, buffsz,0, times + k,NULL,&(fstruct->status) );		
		ffgcv(fstruct->fptr, TDOUBLE, 2, row, 1, buffsz,0, scattit, NULL, &(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 3, row, 1, buffsz,0, scattit+1*buffsz, NULL,&(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 4, row, 1, buffsz,0, scattit+2*buffsz, NULL,&(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 5, row, 1, buffsz,0, scattit+3*buffsz, NULL,&(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 9, row, 1, buffsz,0, ecipos, NULL, &(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 10,row, 1, buffsz,0, ecipos+1*buffsz, NULL, &(fstruct->status) );
		ffgcv(fstruct->fptr, TDOUBLE, 11,row, 1, buffsz,0, ecipos+2*buffsz, NULL, &(fstruct->status) );
	
		if ( fstruct->status ) break;
		
		for (buff_j=0; buff_j<buffsz; buff_j++) {
			//cout << setprecision(16) << times[buff_j+k] << " : " << k + buff_j << " : " << scattit[buff_j] << "," << scattit[buff_j+1*buffsz] << "," << scattit[buff_j +2*buffsz] << "," << scattit[buff_j +3*buffsz ] << endl;
			
			(quat + k + buff_j)->set( scattit + buff_j, buffsz );
			(pos  + k + buff_j)->set( ecipos  + buff_j, buffsz );
		}
		
		row+=buffsz;
		k+=buffsz;
		this->setOpRows(k);
	}
	
	delete [] scattit;
	delete [] ecipos;
	
	return fstruct->status;
}

/*
int PoshistReader::getScAttPos( double * times, double * quaternions, double * xyz_vecs, long rownum, long numrows )
{

	FxbReadStat * fstruct = this->getFxb();
	
	if (fstruct->status != 0) return fstruct->status;
	if (! fstruct->okForOps() ) return -1;
	
	if (numrows <= 0) { numrows = 1; }
	
	long row = rownum;
	long lastrow = row+numrows-1;
	long buffsz = fstruct->n_eff_rows;
	size_t k=0;
	
	if ( row+lastrow-1 > buffsz ) buffsz = row+lastrow-1;
	
	double * scattit = new double[4*buffsz];
	double * ecipos  = new double[3*buffsz];
	
	while ( row <= lastrow && fstruct->status == 0 ) {
		
		cout << row << endl;
		
		if ( row+buffsz-1 > lastrow ) buffsz = lastrow-row+1;
	
		ffgcv(stat.fptr, TDOUBLE, 1, row, 1, 1,0, times, NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 2, row, 1, 1,0, scattit + k*buffsz    , NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 3, row, 1, 1,0, scattit + (k+1)*buffsz, NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 4, row, 1, 1,0, scattit + (k+2)*buffsz, NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 5, row, 1, 1,0, scattit + (k+3)*buffsz, NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 9, row, 1, 1,0, ecipos  + k*buffsz    , NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 10,row, 1, 1,0, ecipos  + (k+1)*buffsz, NULL, &(fstruct->status) );
		ffgcv(stat.fptr, TDOUBLE, 11,row, 1, 1,0, ecipos  + (k+2)*buffsz, NULL, &(fstruct->status) );
		
		k += buffsz;
		row += buffsz;
	}
	
	
	return fstruct->status;
}
*/

int PoshistReader::getScAttPos( vector<double>& tstart, vector<double>& tstop, vector<geom::gquat>& quats, vector<geom::x3>& pos, long startRow, long endRow )
{
	this->setOpRows(0);
	
	if ( this->stat.move2hdu( (char*)attExtName.c_str()) ) return this->stat.status;
	
	long stride = stat.n_eff_rows;
	long n_rows = stat.n_rows;
	
	if ( endRow > 0 && startRow > 0 ) {
		n_rows = endRow - startRow + 1;
	} else {
		endRow = n_rows;
		startRow = 1;
	}
	
	tstart.reserve(n_rows);
	tstop.reserve(n_rows);
	quats.reserve(n_rows);
	pos.reserve(n_rows);
	
	double * tstarts = new double[stride];
	double * scattit = new double[4*stride];
	double * ecipos = new double[3*stride];
	
	long i=startRow;
	
	while( stat.status == 0 && i <= endRow ) {
	
		if ( (i + stride - 1) > endRow) stride = endRow - i + 1;
	
		ffgcv(stat.fptr, TDOUBLE, 1, i, 1, stride,0, tstarts,NULL,&(this->stat.status) );
		
		ffgcv(stat.fptr, TDOUBLE, 2, i, 1, stride,0, scattit, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 3, i, 1, stride,0, scattit+1*stride, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 4, i, 1, stride,0, scattit+2*stride, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 5, i, 1, stride,0, scattit+3*stride, NULL,&(this->stat.status) );
		
		ffgcv(stat.fptr, TDOUBLE,  9, i, 1, stride,0, ecipos, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 10, i, 1, stride,0, ecipos+1*stride, NULL,&(this->stat.status) );
		ffgcv(stat.fptr, TDOUBLE, 11, i, 1, stride,0, ecipos+2*stride, NULL,&(this->stat.status) );
		
	
		for (int k=0;k<stride;k++) {
			tstart.push_back( tstarts[k] );
			
			quats.push_back( geom::gquat( scattit[k], scattit[k+1*stride], scattit[k+2*stride], scattit[k+3*stride] ) );
			pos.push_back( geom::x3( ecipos[k], ecipos[k+1*stride], ecipos[k+2*stride]) );
		}
		
		if (stat.status != 0) break;
		
		this->setOpRows(i);
        i+=stride;

	};
	
	
	std::copy( tstart.begin(), tstart.end()-1, tstop.begin() );
	//The last entry has tstart=tstop
	tstop.push_back( tstarts[stride-1] );
	
	delete [] tstarts;
	delete [] scattit;
	delete [] ecipos;
	
	return stat.status;
}

int TrigdatWriter::writeTimes( double * tstart, double * tstop, long firstRow, long numRows )
{
	FxbReadStat * fstruct = this->getFxb();
	
	this->setOpRows(0);
	if (fstruct->status != 0) return fstruct->status;
	if (! fstruct->okForOps() ) return -1;
	if (numRows <= 0) { numRows = 1; }
	
	ffgnrw(fstruct->fptr,&(fstruct->n_rows),&(fstruct->status) );
	ffgrsz(fstruct->fptr,&(fstruct->n_eff_rows),&(fstruct->status) );
	    
	long row = firstRow;
	long lastrow = row+numRows-1;
	long buffsz = fstruct->n_eff_rows;
	size_t k=0;

	if ( row+lastrow-1 < buffsz ) buffsz = row+lastrow-1;
    
    buffsz=1;
	
	while ( row <= lastrow && fstruct->status == 0 ) {
		
		//cout << row << "," << lastrow << "," << buffsz << "," << fstruct->status  << endl;
		
		if ( row+buffsz-1 > lastrow ) buffsz = lastrow-row+1;
		
		if (tstart != NULL) ffpcl( fstruct->fptr, TDOUBLE, 1, row, 1, buffsz, tstart+k, &(fstruct->status) );
		if (tstop != NULL)  ffpcl( fstruct->fptr, TDOUBLE, 2, row, 1, buffsz, tstop +k, &(fstruct->status) );
		
		if ( row == 1 ) {
		//	ffflus( fstruct->fptr, &(fstruct->status) );
		}
		
		row += buffsz;
		k += buffsz;
		
		this->setOpRows(k);
	}
	//ffflus( fstruct->fptr, &(fstruct->status) );
	ffpdat(fstruct->fptr, &(fstruct->status) );
	ffpcks(fstruct->fptr, &(fstruct->status) );

	return fstruct->status;
}


int TrigdatWriter::writeRateData( long detNum, double * rates, int nChan, long firstRow, long numRows )
{
	FxbReadStat * fstruct = this->getFxb();
	
	this->setOpRows(0);
	
	if (fstruct->status != 0) return fstruct->status;
	if (! fstruct->okForOps() ) return -1;

	if (numRows <= 0) return -1;
	
	long row = firstRow;
	long lastRow = firstRow + numRows - 1;
	long q=0;
	
    long currentDims[2];
	long axDims[2] = { nChan, 14 };
	
	while ( row <= lastRow && fstruct->status == 0 ) {
	
		//if ( row+buffsz-1 > lastrow ) buffsz = lastrow-row+1;
		
		//cout << nChan*(detNum) + 1 << " , " << row << " " << (rates + nChan*q)[0] << endl;
		
		ffpcl( fstruct->fptr, TDOUBLE, 5, row, nChan*(detNum) + 1, nChan, rates + nChan*q, &(fstruct->status) );
		
		
		if ( row == 1 ) {
			ffflus( fstruct->fptr, &(fstruct->status) );
		}

		row++;
		q++;
		this->setOpRows(q);
		
		
	}

//	cout << "in here" << endl;
    
    int ndims;

	ffuky(fstruct->fptr, TINT, (char *)"DETCHANS", &nChan, NULL, &(fstruct->status) );
    
    ffgtdm(fstruct->fptr, 5, 2, &ndims, (long int*)&currentDims, &(fstruct->status) );
    if (ndims == 2 && (currentDims[0] != axDims[0] || currentDims[1] != axDims[1])) {
        char tdim[20];
        sprintf(tdim, "(%ld,%ld)", axDims[0], axDims[1]);
        ffuky(fstruct->fptr, TSTRING, "TDIM5", tdim, NULL, &(fstruct->status));
        //ffptdm(fstruct->fptr, 5, 2, axDims, &(fstruct->status) );
        cout << "Updated TDIM5" << endl;
    }
	ffpdat(fstruct->fptr, &(fstruct->status) );
	ffpcks(fstruct->fptr, &(fstruct->status) );

//	cout << "in here2" << endl;

	return fstruct->status;

}


int TrigdatWriter::writeScAttPos( double * tstart, double * tstop, geom::gquat * quats, geom::x3 * pos, long rownum, long numrows ) 
{
	FxbReadStat * fstruct = this->getFxb();
	
	this->setOpRows(0);
	
	if (fstruct->status != 0) return fstruct->status;
	if (! fstruct->okForOps() ) return -1;
	
	if (numrows <= 0) { numrows = 1; }
	
	long row = rownum;
	long lastrow = row+numrows-1;
	//long buffsz = fstruct->n_eff_rows;
	//long buffsz = 1;
	size_t k=0;
	
//	if ( row+lastrow-1 < buffsz ) buffsz = row+lastrow-1;
	
	while ( row <= lastrow && fstruct->status == 0 ) {
		
		//cout << row << "," << lastrow << "," << buffsz << endl;
		
//		if ( row+buffsz-1 > lastrow ) buffsz = lastrow-row+1;
		
		if (tstart != NULL) ffpcl( fstruct->fptr, TDOUBLE, 1, row, 1, 1, tstart+k, &(fstruct->status) );
		if (tstop != NULL)  ffpcl( fstruct->fptr, TDOUBLE, 2, row, 1, 1, tstop +k, &(fstruct->status) );
		ffpcl( fstruct->fptr, TDOUBLE, 3, row, 1, 4, quats[k].q, &(fstruct->status) );
		ffpcl( fstruct->fptr, TDOUBLE, 4, row, 1, 3,   pos[k].x, &(fstruct->status) );
		
		if ( row == 1 ) {
			ffflus( fstruct->fptr, &(fstruct->status) );
		}
		
		k++;
		row++;
		
		this->setOpRows(k);
	}
	
	//ffflus(fstruct->fptr, &(fstruct->status) );
	
	ffpdat(fstruct->fptr, &(fstruct->status) );
	ffpcks(fstruct->fptr, &(fstruct->status) );
	
	if (fstruct->status != 0) fstruct->errmsg();
	
	return fstruct->status;
};
