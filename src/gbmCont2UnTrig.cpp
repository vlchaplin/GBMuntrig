/*
 *  gbmCont2UnTrig.cpp
 *  UntrigData
 *
 *
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#include "Gbm_PositionIO.h"
#include "PHA2_IO.hh"
#include "GeoTransform.h"
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include "spoccExeUtilities.h"
#include "DBStringUtilities.hh"

//#define UNTRIG_TEMPLATE (char*)"/gbm/Codes/Untrig_proj/fits_template/untrigdat.txt"

#ifdef TRIGDATATEMPLATE
#define UNTRIG_TEMPLATE TRIGDATATEMPLATE
#endif

#ifndef UNTRIG_TEMPLATE
#define UNTRIG_TEMPLATE getenv((char*)"UNTRIG_TEMPLATE")
#endif

#ifndef DEFAULT_PRETRIG_SECONDS
#define DEFAULT_PRETRIG_SECONDS 200.0
#endif

#ifndef DEFAULT_POSTTRIG_SECONDS
#define DEFAULT_POSTTRIG_SECONDS 600.0
#endif

#ifndef DEFAULT_FULLRES_START
#define DEFAULT_FULLRES_START -10.0
#endif

#ifndef DEFAULT_FULLRES_END
#define DEFAULT_FULLRES_END 10.0
#endif

#ifndef DEFAULT_COARSE_RESOLUTION
#define DEFAULT_COARSE_RESOLUTION 4.096  //seconds
#endif

#ifndef DEFAULT_OBJECT
#define DEFAULT_OBJECT "UN"
#endif

#ifndef UNTRIG_EXE
#define UNTRIG_EXE "gbmCont2Untrig"
#endif

#ifndef PACKAGE_NAME
#define PACKAGE_NAME UNTRIG_EXE
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "xcode" 
#endif

inline double nearest( double val, double interval ) {

	double val0 = interval*floor( val / interval );
	double val1 = val0 + interval;
	
	if ( abs(val0 - val) < abs(val1 - val) ) return val0;
	else return val1;
	
};

using namespace std;
using namespace geom;

template<typename T>
inline string val2string(T value)
{
    stringstream temp;
    temp << value;
    return temp.str();
};

void _print_usage() {
    int nameW=10;
    int descW=80;
    cout << endl;
    cout << NewArgument("Usage: "+string(UNTRIG_EXE), "[-met -mjd -utc] trigtime [-fmt utc_format] -ctime [*.pha] -ph poshist.fit -alg 'search_algorithm'  [--help --version]  ", 25, 125) << endl;
    cout << NewArgument("", "[-pre sec -post sec] [-fstart time -fend time]", 25, 75) << endl;
    cout << NewArgument("", "[-coarse sec -obj 'object_string' -v]", 25, 75) << endl << endl;
    cout << NewArgument("Argument", "Value description", nameW, descW) << endl;
    cout << setfill('-') << setw(nameW+descW+2) << "" << endl;
    cout << setfill(' ');
    cout << NewArgument("Time arg:", "Trigger time (T0) in MET, MJD, or UTC time (see below).", nameW, descW) << endl;
    cout << NewArgument("-met", "", nameW, descW) << endl;
    cout << NewArgument("-mjd", "", nameW, descW) << endl ;
    cout << NewArgument("-utc", "[-fmt string] to optionally give posix UTC format code", nameW, descW) << endl << endl;
    cout << NewArgument("-ph", "GBM poshist file containing data for the given time range" , nameW, descW) << endl;
    cout << NewArgument("-ctime", "Complete list of GBM CTIME data files that contain data for the given time", nameW, descW) << endl;
    cout << NewArgument("-alg", "USER_ALG keyword value. Use an informative string such as the name search program", nameW, descW) << endl;

    cout << endl;
    cout << "Optional Arguments:" << endl;
    cout << NewArgument("-pre", "Number of seconds to include before the event (positive number)", nameW, descW) << endl;
    cout << NewArgument("", "Default = "+val2string(DEFAULT_PRETRIG_SECONDS), nameW, descW) << endl;
    cout << NewArgument("-post", "Number of seconds to include after the event (positive number)", nameW, descW) << endl;
    cout << NewArgument("", "Default = "+val2string(DEFAULT_POSTTRIG_SECONDS), nameW, descW) << endl;
    cout << NewArgument("-fstart", "Start time (relative to T0) of full time resolution", nameW, descW) << endl;
    cout << NewArgument("", "Default = "+val2string(DEFAULT_FULLRES_START), nameW, descW) << endl;
    cout << NewArgument("-fend", "End time (relative to T0) of full time resolution", nameW, descW) << endl;
    cout << NewArgument("", "Default = "+val2string(DEFAULT_FULLRES_END), nameW, descW) << endl;
    cout << NewArgument("-coarse", "Coarse time resolution in seconds", nameW, descW) << endl;
    cout << NewArgument("-obj", "OBJECT kewyord value. Leave blank for the default:"+val2string(DEFAULT_OBJECT)+"yymmddfff", nameW, descW) << endl << endl;
    cout << NewArgument("-v", "Be verbose with output", nameW, descW) << endl;
    cout << NewArgument("--help", "Print this help and exit", nameW, descW) << endl;
    cout << NewArgument("--version", "Print version number and exit", nameW, descW) << endl;
    printSingleTimeHelp();
    
    cout << setfill('-') << setw(nameW+descW+2) << "" << endl;
    cout << setfill(' ');
    cout << "Example:" << endl;
    cout << string(UNTRIG_EXE) << " -utc '2009-10-05 13:30:05.256' -alg 'SearchProg v0.0' -ctime glg_ctime*091005*.pha -ph glg_poshist_all_091005_v00.fit" << endl << endl;
    
}


int main ( int argc, char ** argv )
{

	char ** argitr;
	double UnTrigTime, trigLowerMargin, trigUpperMargin, fullResStart, fullResEnd, coarseResolution;
	string timeSys;
	vector<string> pha_files;
	string testCtimeFile;
    string outputUntrigFile;
    string objectString;
    string algString;
    string posFile;
    
    bool verboseMode=0;
	
	int detNum;
    char detName[20];
	
	/* Parse all arguments */
	
	if ( getCmdOption( (char**)argv, (char**)argv + argc,"--help",0) ) {
		_print_usage();
		return 0;
	}
    if ( getCmdOption( (char**)argv, (char**)argv + argc,"--version",0) ) {
		cout << PACKAGE_VERSION << endl;
		return 0;
	}
	if ( getTimeArgMET((char**)argv, (char**)(argv + argc), UnTrigTime, timeSys, NULL, argc) )
	{
        cout << __FUNCTION__ << endl;
        cout << " * Event MET = " << setprecision(16) << UnTrigTime << endl;
	
	} else {
		cout << " ***** Error parsing time argument" << endl;
		return 1;
	}
	
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-v",0);
	if (argitr != NULL) {
		verboseMode=true;
	}
    
    argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-ph",0);
	if (argitr == NULL) {
		cout << "Specify poshist file(s) containing data for the given time range" << endl;
		return 1;
	} else {
        posFile = *(argitr+1);
    }
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-ctime",0);
	if (argitr == NULL) {
		cout << "Specify data files with the '-ctime file1.pha file2.pha ...' argument" << endl;
		return 1;
	}
	int haveDataMask=0;
	argitr++;
	while( (argitr != (char**)argv + argc) && strncmp(*argitr, (char*)"-", 1) ) {
	
		testCtimeFile = string(*argitr);
		argitr++;
		
		PHA2Reader * tempPhaReader = new PHA2Reader;  
		tempPhaReader->setFile( testCtimeFile );
		tempPhaReader->open(0,READWRITE);
		if ( tempPhaReader->read_key( (char*)"DETNAM", TSTRING, detName ) ) 
		{
			cout << "Unable to read DETNAM keyword" << endl;
            tempPhaReader->errclose();
			delete tempPhaReader;
			continue;
		}
		detNum = gbmDetname2Num( detName );
		
		haveDataMask |= (int)pow( 2.0, (int)detNum );
		delete tempPhaReader;
		
		pha_files.push_back(testCtimeFile);
	}
	
	if ( haveDataMask < 2047 )
	{
		for (detNum=0; detNum < 12; detNum++ ) {
			if ( haveDataMask & (int)pow( 2.0, (int)detNum ) ) continue;
			
			cout << "No data for detector: " << gbmDetLongname(detNum) << endl;
		}
		//return 1;
	}
    
    argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-alg",1);
	if (argitr != NULL)
	{
        algString = string( *(argitr+1) );
		cout << " * USER_ALG = " << algString << endl;
	} else {
        cout << "-alg option required" << endl << endl;
        return 1;
    }
    argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-obj",1);
	if (argitr != NULL)
	{
        objectString = string( *(argitr+1) );
		cout << " * OBJECT = " << objectString << endl;
	} 
    
    if (objectString.length() == 0) {
        cout << " * Using default object name" << endl;
    }
	
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-coarse",1);
	if (argitr == NULL)
	{
		cout << " * Coarse resolution set to the default " << setprecision(8) << DEFAULT_COARSE_RESOLUTION << " seconds" << endl;
		coarseResolution = DEFAULT_COARSE_RESOLUTION;
	} else {
		coarseResolution = atof( *(argitr+1) );
		if ( coarseResolution < 1.024 ) {
			cout << endl << " ***** Not a valid value for coarse time resolution (need real number >= 1.024) " << endl;
			return 1;
		} else if ( coarseResolution > 20.0 ) {
			cout << endl << " ***** Coarse resolution must be less than 20 seconds " << endl;
			return 1;
		}
		
		coarseResolution = nearest( coarseResolution, 0.256 ) ;
		
		cout << " * Coarse resolution set to " << coarseResolution << " seconds (multiple of 256ms)" << endl;		
	}
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-pre",1);
	if (argitr == NULL)
	{
		cout << " * Pre-trigger interval to the default: " << setprecision(8) << DEFAULT_PRETRIG_SECONDS << " seconds" << endl;
		trigLowerMargin = DEFAULT_PRETRIG_SECONDS;
	} else {
		trigLowerMargin = atof( *(argitr+1) );
		if ( trigLowerMargin < DEFAULT_PRETRIG_SECONDS ) {
			//cout << endl << "Error: Pre-trigger interval must be at least " << setprecision(8) << DEFAULT_PRETRIG_SECONDS << " seconds"<< endl;
			//return 1;
		} else {
			cout << " * Pre-trigger interval set to " << setprecision(8) << trigLowerMargin << endl;
		}
	}
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-post",1);
	if (argitr == NULL)
	{
		cout << " * Post-trigger interval to the default: " << setprecision(8) << DEFAULT_POSTTRIG_SECONDS << " seconds" << endl;
		trigUpperMargin = DEFAULT_POSTTRIG_SECONDS;
	} else {
		trigUpperMargin = atof( *(argitr+1) );
		if ( trigUpperMargin < DEFAULT_POSTTRIG_SECONDS ) {
			//cout << endl << "Error: Post-trigger interval must be at least " << setprecision(8) << DEFAULT_POSTTRIG_SECONDS << " seconds" << endl;
			//return 1;
		} else {
			cout << " * Post-trigger interval set to " << setprecision(8) << trigUpperMargin << endl;
		}
	}
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-fstart",1);
	if (argitr == NULL)
	{
		cout << " * '-fstart' set to the default " << setprecision(8) << DEFAULT_FULLRES_START << " seconds" << endl;
		fullResStart = DEFAULT_FULLRES_START;
	} else {
		fullResStart = atof( *(argitr+1) );
		if ( fullResStart < -trigLowerMargin ) {
		
			fullResStart = -trigLowerMargin;
			//cout << endl << "Error: Post-trigger interval must be at least " << setprecision(8) << DEFAULT_POSTTRIG_SECONDS << " seconds" << endl;
			//return 1;
		} else {
			cout << " * Full time resolution begins at " << setprecision(8) << fullResStart << endl;
		}
	}
	argitr = getCmdOption( (char**)argv, (char**)argv + argc,"-fend",1);
	if (argitr == NULL)
	{
		cout << " * '-fend' set to the default " << setprecision(8) << DEFAULT_FULLRES_END << " seconds" << endl;
		fullResEnd = DEFAULT_FULLRES_END;
	} else {
		fullResEnd = atof( *(argitr+1) );
		if ( fullResEnd > trigUpperMargin ) {
		
			fullResEnd = trigUpperMargin;
			//cout << endl << "Error: Post-trigger interval must be at least " << setprecision(8) << DEFAULT_POSTTRIG_SECONDS << " seconds" << endl;
			//return 1;
		} else {
			cout << " * Full time resolution ends at " << setprecision(8) << fullResEnd << endl;
		}
	}

    
    char dateAbbrev[7];
    char unNum[20];
    char trigDetBits[15] = "11111111111100";
    int frOfDay = glast_met_integer_dayfraction(UnTrigTime, 3);
    
    glast_met_strftime(UnTrigTime, dateAbbrev, (char*)"%y%m%d");
    sprintf(unNum, (char*)"%s%03d", dateAbbrev, frOfDay);
    outputUntrigFile = "glg_trigdat_all_un"+string(unNum)+"_v00.fit";
	
    if ( objectString.length() == 0 ) {
        objectString = string(DEFAULT_OBJECT)+string(unNum);
    }
    
    cout << "Interpolated file name: " << outputUntrigFile << endl << endl;
    


//	string posFile = "/gbm/fastcopy/archive/level1/poshist/glg_poshist_all_091005_v00.fit";
//	testCtimeFile = "/gbm/fastcopy/archive/level1/continuous/2009-10-05/glg_ctime_na_091005_v00.pha";
	//string testCtimeFile = "glg_ctime_na_091005_v00.pha";

	/* Variable declarations */
	GbmPosReader * posReader;
	FxbReadStat * fxbControl;
	FxbReadStat * trigFxbOut;
	fitsfile * fptr;
	int status;
	double tstart;
	double tstop;
    double newFileStart;
    double newFileEnd;
	
	bool needPrevDay=false;
	bool needNextDay=false;
	
	   
	TrigdatWriter trigInit,trigWrite;
	TrigdatReader trigRead;

	/* variables for poshist reading, and interpolated Quats, XYZ vectors */
	double * inputPosTimes;
	geom::gquat * inputQUATs;
	geom::x3 * inputXYZs;
	geom::gquat * outQUATs;
	geom::x3 * outXYZs;
	//-------------------------
	
	/* variables to hold PHA data */ 	
	double * inputStarts;
	double * inputStops;
	double * inputCounts;
    double * inputExposure;
	//
	
	/* For PHA rebinning algorithm */
	double liveTime, targetBinSize, segmentBegin, segmentRange;
	double * outputExposure;
	double * outputCounts;
	double * outputTimes;
	double * outputStops;
	
	/* miscellaneous; iterators */
	long ri, bi, ck;
	long fineTimeStartIdx1, fineTimeEndIdx1, maxBins, triggerBin;
	long ctimeWindowStart, ctimeWindowEnd, NumCtimeRows, ctimeRow, fullResStartRow, fullResEndRow, numFullResRows, nChan;
	long poshistWindowStart, poshistWindowEnd, poshistEvRow, NumPosRows;
	char text[120];
	
	
	trigInit.setFile(outputUntrigFile);
	
	trigFxbOut = trigInit.getFxb();
	trigFxbOut->template_init( (char*)(trigInit.getFile()).c_str(), UNTRIG_TEMPLATE );
	
	if ( trigInit.status() ) {
		cerr << "Unable to create UNTRIG file: " << trigInit.getFile() << endl;
		trigInit.errmsg();
		return 1;
	} else {
		//trigWrite.flush();
		trigInit.close();
		//trigWrite.move2hdu((char*)"EVNTRATE");
	}
    
	trigWrite.setFile(outputUntrigFile);

    /****************************************************
    
     Loop over PHA files.  Read and bin the spectra, and write the data segment to the output file.
    
    ****************************************************/
	for( size_t fn=0; fn < pha_files.size(); fn++ ) 
	{
		cout << "Reading:"<<endl<<"------------------------" << endl << pha_files[fn] << endl << endl;
		
		testCtimeFile = pha_files[fn];
		PHA2Reader phaReader;  
		phaReader.setFile( testCtimeFile );
		phaReader.open(1,READWRITE);
		phaReader.errmsg();
		if ( phaReader.move2hdu((char*)"EBOUNDS") ) {
			cout << "Error accessing EBOUNDS" << endl;
			//phaReader.errclose();
			//delete phaReader;
			return 1;
		} else
		nChan = phaReader.n_rows();
		
		if ( phaReader.move2hdu((char*)"SPECTRUM") ) {
			cout << "Error accessing SPECTRUM" << endl;
			//phaReader.errclose();
			//delete phaReader;
			return 1;
		} 
		
		if ( phaReader.findRow( UnTrigTime, ctimeRow ) ) 
		{
			cout << endl;
			cout << "Trigger time not found in file: " << phaReader.getFile() << endl;
			//phaReader.close();
			//delete phaReader;
			return 1;
		}
		if ( phaReader.findRow( UnTrigTime - trigLowerMargin, ctimeWindowStart ) ) 
		{
			cout << endl;
			cout << "Start time not found in file. " << endl;
			cout << "Using the first row" << endl;
			ctimeWindowStart = 1;
		}
		if ( phaReader.findRow( UnTrigTime + trigUpperMargin, ctimeWindowEnd ) ) 
		{
			cout << endl;
			cout << "End time not found in file. " << endl;
			cout << "Using the last row" << endl;
			ctimeWindowEnd = phaReader.n_rows();
			//return 1;
		} else {
			ctimeWindowEnd++;
		}
		
		phaReader.read_key( (char*)"DETNAM", TSTRING, detName );
		detNum = gbmDetname2Num( detName );
		
	//	phaReader.move2hdu((char*)"GTI");
	//	phaReader.search_Col = 1;
		phaReader.findRow( UnTrigTime, ctimeRow );
		cout << "Event Row: " << ctimeRow << endl;

		phaReader.findRow( UnTrigTime + fullResStart, fullResStartRow );
		phaReader.findRow( UnTrigTime + fullResEnd, fullResEndRow );
		
		phaReader.getTimes( &fullResStart, NULL, fullResStartRow, 1 );
		phaReader.getTimes( NULL, &fullResEnd,  fullResEndRow, 1 );
		phaReader.getTimes( &tstart, NULL, 1, 1 );
		phaReader.getTimes( NULL, &tstop, phaReader.n_rows(), 1 );
		
		if ( ctimeWindowStart <= 1 || (tstart > UnTrigTime - trigLowerMargin) ) {
			if ( tstart > UnTrigTime - trigLowerMargin ) needPrevDay = true;
			ctimeWindowStart = 1;
		}
		if ( ctimeWindowEnd > phaReader.n_rows() || (tstop < UnTrigTime + trigUpperMargin) ) {
			if ( tstop < UnTrigTime + trigUpperMargin ) needNextDay = true;
			ctimeWindowEnd = phaReader.n_rows();
		}
		
		NumCtimeRows = ctimeWindowEnd - ctimeWindowStart + 1 ;
		numFullResRows = fullResEndRow - fullResStartRow + 1;
		
		inputStarts = new double[ NumCtimeRows+1 ];
		inputStops =  new double[ NumCtimeRows+1 ];
        inputExposure = new double[ NumCtimeRows+1 ];
		inputCounts = new double[ nChan*(NumCtimeRows+1) ] ;

		phaReader.getTimes( inputStarts, inputStops, ctimeWindowStart, NumCtimeRows );
		phaReader.getSpectra( inputCounts, inputExposure, ctimeWindowStart, NumCtimeRows );
		if ( phaReader.status() ) {
			phaReader.errclose();
			//delete phaReader;
			return 1;
		} else {
			phaReader.close();
			//delete phaReader;
		}	
		/*Start of the PHA data selection*/
		segmentBegin = inputStarts[0];
		/*Length of the PHA data selection*/
		segmentRange = inputStops[NumCtimeRows-1] - inputStarts[0];
		
		
		
		fineTimeStartIdx1 = fullResStartRow - ctimeWindowStart;
		fineTimeEndIdx1 = fullResEndRow - ctimeWindowStart;
		
		// Maximum number bins possible
		maxBins = floor( segmentRange / coarseResolution ) - floor( ( inputStops[fineTimeEndIdx1] - inputStarts[fineTimeStartIdx1] ) / coarseResolution ) + numFullResRows;
		
		outputExposure = new double[ maxBins+1 ];
		outputCounts = new double[ nChan*(maxBins+1) ];
		outputTimes = new double[ maxBins+1 ];
		outputStops = new double[ maxBins+1 ];
		
		for ( bi=0; bi < maxBins; bi++ ) {
			
			for ( ck=0; ck < nChan; ck++ ) outputCounts[bi*nChan + ck] = 0;
			
			outputTimes[bi] = -1;
			outputExposure[bi] = 0;
		}
        cout << detName << " Ctime rows:" << endl;
		cout << setw(10) << "Start: " << ctimeWindowStart << endl;
		cout << setw(10) << "Event: " << ctimeRow << endl;
		cout << setw(10) << "End: " << ctimeWindowEnd << endl;
		cout << setw(10) << "Num: " << NumCtimeRows << endl;
		
		
		/* Loop over spectra.  If a bin is in the coarse-time region (user specified) then add bins
		until the desired livetime is approximated.  Bins in the full-time region are merely copied to
		the output spectrum. 
		*/
		ri=0;
		bi=0;
        triggerBin=-1;
		while ( ri < NumCtimeRows && bi < maxBins )
		{
		
			//sprintf( text, (char*)"%9ld : [%6ld]", ri, bi );
			//cout << text << endl;
		
		
			/* Full resolution region */
			if (ri >= fineTimeStartIdx1 && 
				ri <= fineTimeEndIdx1
				) 
			{
				outputTimes[bi] = inputStarts[ri];
				outputStops[bi] = inputStops[ri];
                outputExposure[bi] += inputExposure[ri];
				for ( ck=0; ck < nChan; ck++ )
					outputCounts[bi*nChan + ck] += inputCounts[ri*nChan + ck];
                
                
                if ( UnTrigTime >= outputTimes[bi] && UnTrigTime <= outputStops[bi] ) triggerBin = bi;
                
				ri++;
				bi++;
				
                if ( verboseMode ) {
                    sprintf( text, (char*)"%9ld : [%6ld] %10s %16f - %16f", ri, bi-1, "bi++, ri++", outputTimes[bi-1] - UnTrigTime, outputStops[bi-1] - UnTrigTime);
                    cout << text << endl;
                }
                
				continue;
			}		
			
			
			
			/* Coarse resolution region */
			liveTime = 0.0;
			targetBinSize = coarseResolution;

			outputTimes[bi] = inputStarts[ri];
			liveTime += inputStops[ri] - inputStarts[ri];
			while ( liveTime < (targetBinSize + 0.001) && 
					ri < NumCtimeRows && 
				   (ri < fineTimeStartIdx1 || ri > fineTimeEndIdx1) 
				) 
			{
				
				for ( ck=0; ck < nChan; ck++ )
					outputCounts[bi*nChan + ck] += inputCounts[ri*nChan + ck];
					
				outputExposure[bi] += inputExposure[ri];
				ri++;
				liveTime += inputStops[ri] - inputStarts[ri];
			}
			
			outputStops[bi] = inputStops[ri-1];
			bi++;
			
            if ( verboseMode ) {
                sprintf( text, (char*)"%9ld : [%6ld] %10s [%16f : %16f]", ri, bi-1, "ri++", outputTimes[bi-1] - UnTrigTime, outputStops[bi-1] - UnTrigTime);
                cout << text << endl;
            }
		}
        
        for (bi=0;bi<maxBins;bi++)
            for ( ck=0; ck < nChan; ck++ ) 
                outputCounts[bi*nChan + ck] /= outputExposure[bi];
		
		/* maxBins now becomes the total number of bins */
		maxBins = bi;
        cout << maxBins << " output bins" << endl;
		
		long orow=1;
		
        trigRead.setFile(trigWrite.getFile());
        trigRead.open(1);
		//trigRead.reopen( &trigWrite );
		trigRead.move2hdu((char*)"EVNTRATE");
	
		if ( trigRead.status() ) {
			cerr << "Error re-opening" << endl;
			trigRead.errclose();
		} else {
			orow = trigRead.findRow( outputTimes[0] ) + 1;
			trigRead.close();
			
			if ( orow <= 0 ) orow = 1;
		}
        trigWrite.open(1, READWRITE);
        trigWrite.move2hdu((char*)"EVNTRATE");
		trigWrite.writeTimes( outputTimes, outputStops, orow, maxBins );
        trigWrite.writeRateData( detNum, outputCounts, nChan, orow, maxBins );
        trigWrite.flush();
        
        if ( triggerBin != -1 ) {
            trigWrite.move2hdu((char*)"TRIGRATE");
            trigWrite.writeTimes( &outputTimes[triggerBin], &outputStops[triggerBin], 1 );
            trigWrite.writeRateData( detNum, &outputCounts[triggerBin*nChan], nChan, 1 );
            
            trigWrite.move2hdu((char*)"MAXRATES");
            trigWrite.writeTimes( &outputTimes[triggerBin], &outputStops[triggerBin], 1 );
            trigWrite.writeRateData( detNum, &outputCounts[triggerBin*nChan], nChan, 1 );
        }

        if ( trigWrite.status() ) trigWrite.errmsg(); 
        else {
            cout << "Wrote data for " << gbmDetLongname( detNum ) << endl << endl;
        }
		
		trigWrite.close();
		
		delete [] outputExposure;
		delete [] outputCounts;
		delete [] outputTimes;
		delete [] outputStops;
				
		delete [] inputStarts;
		delete [] inputStops;
		delete [] inputCounts;
        delete [] inputExposure;
		
	//	phaReader.close();
	//	delete phaReader;
		
	}
	
	cout << "end" << endl;
	//trigWrite.close();
	//trigRead.close();
	
/*
	Open the POSHIST file
*/	
	
	posReader = NewPosReader::type((char*)"poshist");	
	posReader->setFile(posFile);
	status = posReader->open(0);
	fxbControl = posReader->getFxb();
	fptr = fxbControl->fptr;
	
	fxbControl->CloseOnError=0;
	
	status = fxbControl->status;
	if ( status != 0 ) { 
		fxbControl->errclose();
		cout << "Error opening " << posReader->getFile() << endl;
		delete posReader; 
		return 0;
	}

	vector<double> tempStr;
	vector<double> tempEnd;
	poshistWindowStart = posReader->findRow( UnTrigTime - trigLowerMargin - 2, &tempStr, &tempEnd );
	poshistWindowEnd = posReader->findRow( UnTrigTime + trigUpperMargin + 2, &tempStr, &tempEnd );
	poshistEvRow = posReader->findRow( UnTrigTime, &tempStr, &tempEnd );
	NumPosRows = poshistWindowEnd - poshistWindowStart + 1;
	
	tempStr.clear();
	tempEnd.clear();
	
	inputPosTimes = new double[NumPosRows];
	inputQUATs = new geom::gquat[NumPosRows];
	inputXYZs = new geom::x3[NumPosRows];
	
	posReader->getScAttPos( inputPosTimes, inputQUATs, inputXYZs, poshistWindowStart, NumPosRows );
	if ( posReader->status() ) {
		posReader->errclose();
		return 1;
	}
	
/***************
*
* Position /attitude interpolation
*
****************/	

    trigRead.setFile( trigWrite.getFile() );
    trigRead.open();
    trigRead.move2hdu((char*)"EVNTRATE");
	maxBins = trigRead.n_rows(); 
    
	
	outputTimes = new double[maxBins+1];
	outputStops = new double[maxBins+1];
	outQUATs = new geom::gquat[maxBins+1];
	outXYZs = new geom::x3[maxBins+1];
    
    if ( trigRead.getTimes(outputTimes, outputStops, &maxBins, false) ) 
    {
        cout << endl;
        cout << "Error re-opening untrigdat" << endl;
        trigRead.errclose();
        return 1;
    } else {
        trigRead.close();
    }
	
	cout << endl;
	cout << "Interpolating " << posReader->operationRows() << " poshist rows to " << maxBins << endl << endl;
    
	bi = 0;
	double frame_sep;
	double binMid;
	double unit_time=0.0;
	
	long k=0;
	
	while (bi < maxBins) {
		
		binMid = 0.5*(outputStops[bi] + outputTimes[bi]);
		
		//if ( binMid < inputPosTimes[k] ) 
		
		while( !( binMid >= inputPosTimes[k] && binMid <= inputPosTimes[k+1]) && k < posReader->operationRows()-1 ) 
		{
            //sprintf( text, (char*)"%10ld , %16f : {%16f , %16f}", k, binMid , inputPosTimes[k] - UnTrigTime, inputPosTimes[k+1] - UnTrigTime  );
            //cout  << text << " " << endl;
            k++;
		}
		
		if ( k == posReader->operationRows() ) {
			cout << "Timing error here" << endl;
			break;
		}

		if ( !( binMid >= inputPosTimes[k] && binMid <= inputPosTimes[k+1]) ) {
			
			//sprintf( text, (char*)"Err- %10ld {%6ld : %6ld}, %16f ne {%16f , %16f }", bi, k, k+1, binMid - UnTrigTime, inputPosTimes[k] - UnTrigTime, inputPosTimes[k+1] - UnTrigTime  );
			//cout  << text << " " << endl;
			
			k++;
			continue;
		}
		
		frame_sep = inputPosTimes[k+1] - inputPosTimes[k];

		unit_time = ( binMid - inputPosTimes[k] ) / frame_sep;

		
		outQUATs[bi].slerp( unit_time, inputQUATs[k], inputQUATs[k+1] );
		outXYZs[bi].linterp( unit_time, inputXYZs[k], inputXYZs[k+1] );
		
		//sprintf( text, (char*)"+++ %10ld {%6ld : %6ld}, %16f  e {%16f , %16f }", bi, k, k+1, binMid - UnTrigTime, inputPosTimes[k] - UnTrigTime, inputPosTimes[k+1] - UnTrigTime  );
		//cout  << text << " " << endl;
		
		
		bi++;
	}

    if ( bi == maxBins ) cout << "No issues" << endl;
    else {
        cout << "There is likely part of an SAA passage in your time interval selection" << endl;
        cout << "Attempting to proceed anyway" << endl << endl;
    }
	
    //If there was an SAA during the time interval, the number of predicted rows will be greater than 
    //what was actually obtained. These would show up as elements where the bin time == -1.0, and are
    //at the end of the data block.  We will just remove them.
    
    triggerBin=-1;
	for ( bi=0; bi < maxBins; bi++ ) {
		outXYZs[bi] /= 1000.0;
        if ( UnTrigTime >= outputTimes[bi] && UnTrigTime <= outputStops[bi] ) triggerBin = bi;
        
        //cout << setw(8) << bi << setw(20) << outputTimes[bi]-UnTrigTime << endl;
        
        if ( outputTimes[bi] < 0 ) 
        {
            
            cout << "Removing excess rows due to SAA passage, " << bi+1 << " -> " << maxBins << endl;
            
            trigWrite.open(1, READWRITE);
            trigWrite.move2hdu("EVNTRATE");
            FxbReadStat * fstruct = trigWrite.getFxb();
            ffdrow(fstruct->fptr, bi+1, maxBins - (bi), &(fstruct->status) );
            trigWrite.up_checksum_date();
            trigWrite.errclose();
            cout << endl;
            
            maxBins = bi; //this loop will terminate
        }
	}
    
    if ( triggerBin == -1 ) {
        cout << "Warning - unable to find the 'trigger' (T0) time bin" << endl;
        cout << "The given T0 is probably in an SAA passage. The nearest time bin will be used." << endl;
        cout << "Keywords like RA_SCZ and GEO_LONG/GEO_LAT are for this time." << endl;
        
        cout << endl;
        cout << "If you see a FITS error 402 and this message: 'Error in ffd2e: double value is a NaN or INDEF'," << endl;
        cout << "this is because the derived TSTART or TSTOP is during an SAA outage.  Re-run your time selection with a larger " << endl;
        cout << "pre- or post- burst interval, so that some data before or after the SAA is captured, using the '-pre' or '-post' argument." << endl;
        
        double minDiff=UnTrigTime;
        
        for ( bi=0; bi < maxBins; bi++ ) 
        {
            if ( abs(UnTrigTime - outputTimes[bi]) > minDiff ) {
                minDiff=abs(UnTrigTime - outputTimes[bi]);
                triggerBin=bi;
            }
        }
        
        cout << endl << " --> Row number " << triggerBin+1 << endl << endl;
        
    }
    
    char date_fmt[] = "%Y-%m-%dT%H:%M:%S";
    char date_obs[30];
    char date_end[30];
    string creator= string(PACKAGE_NAME) + " v" + string(PACKAGE_VERSION);
    newFileStart = outputTimes[0];
    newFileEnd = outputStops[maxBins-1];
    
    glast_met_strftime(newFileStart, date_obs, date_fmt);
	glast_met_strftime(newFileEnd  , date_end, date_fmt);
    
    x3 xaxis(1,0,0);
    x3 zaxis(0,0,1);
    double x_ra,x_dec,z_ra,z_dec;
    
    cart2sphere( outQUATs[triggerBin].trot(xaxis) , x_ra, x_dec, 1);
    cart2sphere( outQUATs[triggerBin].trot(zaxis) , z_ra, z_dec, 1);
    
    vGeographic geo_pos = gei2geo(UnTrigTime, outXYZs[triggerBin]);
    
    geo_pos.lat /= DTOR;
    geo_pos.lon /= DTOR;
    geo_pos.lat = nearest( geo_pos.lat , 0.001);
    geo_pos.lon = nearest( geo_pos.lon , 0.001);
    
    z_ra = nearest(z_ra, 0.001);
    z_dec = nearest(z_dec, 0.001);
    x_ra = nearest(x_ra, 0.001);
    x_dec = nearest(x_dec, 0.001);
    
    trigWrite.open(1, READWRITE);
    for (int hdu_num=1; hdu_num <=6; hdu_num++) 
    {
        trigWrite.move2hdu(hdu_num);

        trigWrite.write_key((char*)"TSTART", TDOUBLE, (void*)&newFileStart);
        trigWrite.write_key((char*)"TSTOP", TDOUBLE, (void*)&newFileEnd);
        trigWrite.write_key((char*)"TRIGTIME", TDOUBLE, (void*)&UnTrigTime);
        trigWrite.write_key((char*)"DATE-OBS", TSTRING, date_obs);
        trigWrite.write_key((char*)"DATE-END", TSTRING, date_end);
        trigWrite.write_key((char*)"OBJECT", TSTRING, (void*)objectString.c_str() );                        
        trigWrite.write_key((char*)"CREATOR", TSTRING, (void*)creator.c_str() );   
        trigWrite.write_key((char*)"DET_MASK", TSTRING, trigDetBits ); 
        if (hdu_num==1) {
            trigWrite.write_key((char*)"FILENAME", TSTRING, (void*)outputUntrigFile.c_str() );
            trigWrite.write_key((char*)"USER_ALG", TSTRING, (void*)algString.c_str() ); 
            trigWrite.write_key((char*)"GEO_LONG", TDOUBLE, &geo_pos.lon);
            trigWrite.write_key((char*)"GEO_LAT", TDOUBLE, &geo_pos.lat);
            trigWrite.write_key((char*)"RA_SCX", TDOUBLE, &x_ra);
            trigWrite.write_key((char*)"DEC_SCX", TDOUBLE, &x_dec);
            trigWrite.write_key((char*)"RA_SCZ", TDOUBLE, &z_ra);
            trigWrite.write_key((char*)"DEC_SCZ", TDOUBLE, &z_dec);
        }
        trigWrite.up_checksum_date();
    }
    
	trigWrite.move2hdu((char*)"EVNTRATE");
	if ( trigWrite.writeScAttPos( NULL, NULL, outQUATs, outXYZs, 1, maxBins ) ) trigWrite.errmsg();	
    
    if ( triggerBin != -1 ) 
    {
        trigWrite.move2hdu((char*)"TRIGRATE");
        if ( trigWrite.writeScAttPos( NULL, NULL, &outQUATs[triggerBin], &outXYZs[triggerBin]) ) trigWrite.errmsg();
        
        
        trigWrite.move2hdu((char*)"MAXRATES");
        if ( trigWrite.writeScAttPos( NULL, NULL, &outQUATs[triggerBin], &outXYZs[triggerBin]) ) trigWrite.errmsg();	
    }
    
    trigWrite.move2hdu((char*)"OB_CALC");
    ffirow( trigWrite.getFxb()->fptr, 0, 1, &(trigWrite.getFxb()->status) );
    
    trigWrite.up_checksum_date();
    
    cout << endl << "----------------" << endl;
    if ( trigWrite.getFxb()->status == 0) cout << "Wrote: " << trigWrite.getFile() << endl << endl;
    
	trigWrite.close();
	
	delete posReader;
    
	delete [] outputTimes;
	delete [] outputStops;
	delete [] outQUATs;
	delete [] outXYZs;
	
	delete [] inputPosTimes;
	delete [] inputQUATs;
	delete [] inputXYZs;
	
	
	return 0;
	
	
};