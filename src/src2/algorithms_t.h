/*
 *  algorithms_t.h
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */

#include <vector>

#ifndef AGLORITHM_T_H
#define AGLORITHM_T_H


using namespace std;

template<class T, class V>
inline size_t findNearestElement( T points, V toValue, size_t maxPoint) {
	
	size_t m=0;
	size_t n=maxPoint;
	
	if ( toValue < points[m] ) return -1;
	if ( toValue > points[n] ) return -1;
	
	return nearestBinStepGeneric( points, toValue, m, n );

};

template<class T, class V>
inline size_t nearestBinStepGeneric( T& points, V& toValue, size_t& lowBx, size_t& upBx ) {

	size_t dx = upBx - lowBx;
	
	if (dx == 0 || upBx == 0) return upBx;
	
	long mdx = ( lowBx+upBx ) / 2;

	V mtime = points[mdx];
	
	if (toValue < mtime) {
		//Branch to the lower half
		upBx = mdx;
		return nearestBinStepGeneric( points, toValue, lowBx, upBx );
	}
	else if (toValue > mtime) {
		lowBx = mdx+1;
		return nearestBinStepGeneric( points, toValue, lowBx, upBx );
	} else {
		
		//if ( (toValue - mtime) > (points[mdx+1] - toValue) ) return mdx+1;
	
		return mdx;
	}
	
};


template<class T>
inline long nearestBinStep( vector<T>& points, T& toValue, long& lowBx, long& upBx ) {
    
	long dx = upBx - lowBx;
	
	if (dx == 0 || upBx == 0) return upBx;
	
	long mdx = ( lowBx+upBx ) / 2;
    
	T mtime = points[mdx];
	
	if (toValue < mtime) {
		//Branch to the lower half
		upBx = mdx;
		return nearestBinStep( points, toValue, lowBx, upBx );
	}
	else if (toValue > mtime) {
		lowBx = mdx+1;
		return nearestBinStep( points, toValue, lowBx, upBx );
	} else {
		
		//if ( (toValue - mtime) > (points[mdx+1] - toValue) ) return mdx+1;
        
		return mdx;
	}
	
};

template<class T>
inline long findNearestElement( vector<T>& points, T toValue) {
	
	long m=0;
	long n=points.size()-1;
	
	if ( toValue < points[m] ) return -1;
	if ( toValue > points[n] ) return -1;
	
	return nearestBinStep( points, toValue, m, n );

};

#endif