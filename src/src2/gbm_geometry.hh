/*
 *  gbm_geometry.hh
 *
 *  Created by Vandiver Chaplin
 *  The University of Alabama in Huntsville
 *
 */
 
#ifndef GBM_GEOMETRY
#define GBM_GEOMETRY 1
 
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <iostream>

#include "DynMatrix.h"
 
using namespace std;

int gbmDetname2Num(char *);
string gbmDetShortname(int);
string gbmDetShortname(char *);
string gbmDetLongname(int det);

namespace geom {
	extern double PI;
	extern double DTOR;
	
	inline double deg2rad(double& d) { return d*DTOR; };
	inline double rad2deg(double& r) { return r/DTOR; };
	
	class x3 {
		public:
		
		double x[3];
		
		x3(){
			x[0] = x[1] = x[2] = 0.0;
		};
		x3(double q0, double q1, double q2) {
			x[0] = q0;
			x[1] = q1;
			x[2] = q2;
		};
		x3(double * vect) {
			set(vect);
		};
		
		x3& set(double * vect, size_t step=1) {
			x[0] = *vect;
			x[1] = *(vect+1*step);
			x[2] = *(vect+2*step);
			return *this;
		};
		
		double& x0() { return x[0]; };
		double& x1() { return x[1]; };
		double& x2() { return x[2]; };
		
		double mag() {
			return sqrt(x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
		};
		friend double dot(x3& ths, x3& that) {
			return ths.x[0] * that.x[0] + ths.x[1] * that.x[1] + ths.x[2] * that.x[2];
		};
		/*
		x3 operator+(x3& that) {
			return x3(x[0] + that.x[0], x[1] + that.x[1], x[2] + that.x[2]);
		};
		x3 operator-(x3& that) {
			return x3(x[0] - that.x[0], x[1] - that.x[1], x[2] - that.x[2]);
		};
		x3 operator*(x3& that) {
			return x3(x[0] * that.x[0], x[1] * that.x[1], x[2] * that.x[2]);
		};
		x3 operator+(double scalar) {
			return x3(x[0] + scalar, x[1] + scalar, x[2] + scalar);
		};
		x3 operator-(double scalar) {
			return x3(x[0] - scalar, x[1] - scalar, x[2] - scalar);
		};
		x3 operator*(double scalar) {
			return x3(x[0] * scalar, x[1] * scalar, x[2] * scalar);
		};
		x3 operator/(double scalar) {
			return x3(x[0] / scalar, x[1] / scalar, x[2] / scalar);
		};
		*/
		
		x3& operator+=(x3& that) {
			x[0] += that.x[0];
			x[1] += that.x[1];
			x[2] += that.x[2];
			return *this;
		};
		x3& operator-=(x3& that) {
			x[0] -= that.x[0];
			x[1] -= that.x[1];
			x[2] -= that.x[2];
			return *this;
		};
		x3& operator=(x3& that) {
			this->set( that.x );
			return *this;
		};
		x3& operator=(double scalar) {
			x[0] = x[1] = x[2] = scalar;
			return *this;
		};
		x3& operator+=(double scalar) {
			x[0] += scalar;
			x[1] += scalar;
			x[2] += scalar;
			return *this;
		};
		x3& operator-=(double scalar) {
			x[0] -= scalar;
			x[1] -= scalar;
			x[2] -= scalar;
			return *this;
		};
		x3& operator*=(double scalar) {
			x[0] *= scalar;
			x[1] *= scalar;
			x[2] *= scalar;
			return *this;
		};
		x3& operator/=(double scalar) {
			x[0] /= scalar;
			x[1] /= scalar;
			x[2] /= scalar;
			return *this;
		};
		template<typename T>
		x3& operator=(vector<T>& vec3) {
			x[0] = vec3[0];
			x[1] = vec3[1];
			x[2] = vec3[2];
			return *this;
		};
		
		x3& linterp( double& unit_time, x3& v1, x3& v2 ) {
		
			x3 temp1 = v1;
			x3 temp2 = v2;
		
			temp1*=(1.0 - unit_time);
			temp2*=unit_time;
			temp2+=temp1;
			this->set( temp2.x );
			return *this;
		};
		
		double& operator[](int i) {
			switch (i % 3) {
				case 0: return x[0];
				case 1: return x[1];
				case 2: return x[2];
				default: return x[0];
			};
		};
		
		friend ostream& operator<<(ostream& os, x3& data)
		{
			os << "(" << data.x[0] << "," << data.x[1] << "," << data.x[2] << ")";
			return os;
		}
		
	};

	template<typename T>
	inline T degsep( T * pnt1, T * pnt2 ) {
		
		vector<T> p1,p2;
		double dot = vdot( sphere2cart(p1, pnt1[0], pnt1[1]), sphere2cart(p2, pnt2[0], pnt2[1]), 0 );
		return acos(dot)/DTOR ;
	};

	template<typename T> 
	inline T vdot(const vector<T>& a,const vector<T>& b, unsigned int cos=0) {
		size_t i,N;
		N = a.size();
		if (b.size() < N) N = b.size();
		
		i=0;
		T dot = 0;
		if (cos == 0) {
			while (i < N) { dot += a[i]*b[i]; i+=1; }
		}
		else {
			T norm1 = 0;
			T norm2 = 0;
			while (i < N) {
				dot += a[i]*b[i];
				norm1 += pow(a[i],2);
				norm2 += pow(b[i],2);
				i++;
			}
			
			//cout << "dot = " << dot << "," << norm1 << "," << norm2 << endl;
			
			dot /= (sqrt(norm1)*sqrt(norm2));
			
			
			if (dot > 1.0) dot = 1.0;
			else if (dot < -1.0) dot = -1.0;
			
			if ( cos > 1 ) dot = acos(dot);
		}
		
		
		return dot;
	};

	/* q = xyquadrant(x,y):
	q : 
			 +y
			  |
			2 | 1 
	   -x ----|---- +x
			3 | 4 
			 -y
	*/
	template<typename T>
	inline int xyquadrant( T& x, T& y) {
				
		int i,j,q,jdim;
		if (x > 0) j=1; else j=0;
		if (y < 0) i=1; else i=0;

		jdim=2;
		if (i % 2 == 0) q = i*jdim + (jdim - j-1); 
		else q = i*jdim + j;

		return q+1;
	};

	template<typename T>
	inline vector<T>& vec3(vector<T>& v,T v0,T v1,T v2) {
		v.resize(3);
		v[0] = v0;
		v[1] = v1;
		v[2] = v2;
		return v;
	};

	template<typename T>
	inline vector<T>& sphere2cart (vector<T>& v, T& az, T& el) {
		
		double x,y,z;
		x = cos(DTOR*az)*cos(DTOR*el);
		y = sin(DTOR*az)*cos(DTOR*el);
		z = sin(DTOR*el);
		
		return vec3( v, x,y,z);
	};
	inline x3& sphere2cart (x3& v, double& az, double& el, int u=0) {

		if (u == 1) {
			v[0] = cos(DTOR*az)*cos(DTOR*el);
			v[1] = sin(DTOR*az)*cos(DTOR*el);
			v[2] = sin(DTOR*el);
		} else {
			v[0] = cos(az)*cos(el);
			v[1] = sin(az)*cos(el);
			v[2] = sin(el);
		}
		
		return v;
	};
	inline void cart2sphere (x3& v, double& az, double& el, int u=0) {
		double r = v.mag();
		az = atan( v[1] / v[0] );
		el = asin( v[2] / r );
		if (v[0] < 0) az += PI;
		
		if (az < 0) az += (PI+PI);
		
        
		if (u == 1) {
			az /= DTOR;
			el /= DTOR;
		};
		
	};
	inline void cart2sphere (x3& v, double& r, double& az, double& el, int u=0) {
		r = v.mag();
		az = atan( v[1] / v[0] );
		el = asin( v[2] / r );
		if (v[0] < 0) az += PI;
		
		if (az < 0) az += (PI+PI);
		
		if (u == 1) {
			az /= DTOR;
			el /= DTOR;
		};
		
	};
	
	class Rmatrix {
		private:
		double _A[3][3];
		
		public:
		
		friend class quaternion;
		
		Rmatrix(){};
		Rmatrix(double g00, double g11, double g22) {
			this->zero();
			this->setdiag(g00,g11,g22);
		};
		
		void setdiag(double g00, double g11, double g22) {
			(*this)(0,0) = g00;
			(*this)(1,1) = g11;
			(*this)(2,2) = g22;
		};
		void zero() {
			size_t i,m,n;
			i=0;
			while(i<9) {
				m = i / 3;
				n = i % 3;
				(*this)(m,n) = 0;
				i++;
			}
		};
		void identity() {
			this->zero();
			this->setdiag(1.0, 1.0, 1.0);
		};
		
		Rmatrix& tr() {
			Rmatrix tmp  = *this;
			size_t i,m,n;
			i=0;
			while(i<9) {
				m = i / 3;
				n = i % 3;
				(*this)(m,n) = tmp(n,m);
				i++;
			}
			return *this;
		};
		
		double& operator() (size_t& m, size_t& n) {
			return _A[m][n];
		};
		double& operator() (int m, int n) {
			return _A[m][n];
		};
		Rmatrix& operator= (double scalar) {
			size_t i,j;
			i=0;
			while(i < 3) {
				j=0;
				while(j < 3) {
					(*this)(i,j) = scalar;
					j++;
				}
				i++;
			}
			
			return *this;
		};
		virtual Rmatrix& operator= (Rmatrix& that) {
			size_t i,j;
			i=0;
			while(i < 3) {
				j=0;
				while(j < 3) {
					(*this)(i,j) = that(i,j);
					j++;
				}
				i++;
			}
			
			return *this;
		};
		Rmatrix operator* (Rmatrix& that) {			
			Rmatrix tmp;
			
			size_t i,j,q;
			i=0;
			while(i < 3) {
				j=0;
				while(j < 3)  {
					q=0;
					tmp(i,j)=0;
					while(q < 3) {
						tmp(i,j) += (*this)(i,q)*that(q,j);
						q++;
					}
					j++;
				}
				
				i++;
			}
			
			return tmp;
		};
		
		//Note about self-assignment:
		/*
		A *= B is equivalent to A = B*A
		
		This operation is the reverse of typical self-assignment operators,
		where A *= B usually is equivalent as A = A * B.
		
		This is to take advantage of
		familiar concepts from matrix algebra where B 'operates' on A,
		and results in a mutated A.
				
		Using self-assignment allows
		mutation with less memory overhead, as the construction
		A = B*A requires holding the result in a temporary allocation
		and using the copy constructor =.
		
		If Pn is the nth product of n+1 matrices, R[j], then the INefficient way to
		calculate Pn is:
		
		P1 = R1*R0;
		P2 = R2*P1 = R2*(R1*R0);
		P3 = R3*P2 = R3*(R2*(R1*R0));
		.
		.
		Pn = Rn * (Rn-1 * (Rn-2 * ... (R3 * (R2 * (R1 * R0))) ... );
		
		This would result in allocating a temporary 3x3 matrix n times, and
		copying it n times.
		
		An equivalent formula which does not allocate any uneccesary memory is:
				
		Pn = ( ... (((R0*=R1) *= R2) *= R3) ... ) *= Rn-2) *= Rn-1) *= Rn );
		
		Now, R0 has been mutated n times.  It's original value is lost because it's memory
		now holds the elements of Pn, the nth product.
		
		In C++, this composition could be written in a loop, for n+1 matrices.
		If we don't want to effect R0, then simply initializing Pn = R0 before
		calculating the product let's us use *= operations on the variable Pn.
		
		
		Rmatrix R[n+1]; //the Rn matrices
		Rmatrix Pn;   //to hold product of the Rn matrices
		...code to define R[n]...
		
		Pn = R0; //initialize
		
		for(i=1; i<n+1; i++) {
			Pn *= R[i];
		}
		
		now Pn equals the compound operation of all R matrices, and it only required an overhead of one 3x3 matrix
		instead of an overhead of n+1 3x3 matrices.
		
		It's much better to use '*=' than '*' for matrix multiplication.
		*/
		
		Rmatrix& operator*= (Rmatrix& that) {
			double ell[3];
			size_t i,j,q;
			i=0;
			while(i < 3) {
				j=0;
				while(j < 3)  {
					q=0;
					ell[j]=0;
					while(q < 3) {
						ell[j] += (*this)(i,q)*that(q,j);
						q+=1;
					}
					j+=1;
				}
				
				j=0;
				while(j < 3) {
					(*this)(i,j) = ell[j];
					j+=1;
				}
				i+=1;
			}
			
			return *this;
		};
		
		template<typename T>
		vector<T> rot(vector<T>& v) {
			size_t i,q;
			vector<T> vprime(3);
			i=0;
			while(i<3) {
				vprime[i]=0;
				q=0;
				while(q < 3) vprime[i] += (*this)(i,q)*v[q++];
				i=i+1;
			}
			return vprime;
		};
		
		x3& rot(x3& point) {
			x3 tmp = point;
			size_t i,q;
			i=0;
			while(i<3) {
				point[i]=0;
				q=0;
				while(q < 3) {
					point[i] += (*this)(i,q)*tmp[q];
					q+=1;
				}
				i=i+1;
			}
			return point;
		};
	};

	extern Rmatrix I3;
	
	
	
	class gquat {
		public:
		
		double q[4];
		
		gquat(){};
		gquat(double qx, double qy, double qz, double qs) {
			q[0] = qx;
			q[1] = qy;
			q[2] = qz;
			q[3] = qs;
		};
		gquat(double * vect, size_t step=1) {
			set(vect, step);
		};
		
		double& qx() { return q[0]; };
		double& qy() { return q[1]; };
		double& qz() { return q[2]; };
		double& qs() { return q[3]; };
		
		friend double dot( gquat& q1, gquat& q2 )
		{
			double dots=0;
			for( size_t i=0; i < 4; i++ ) dots += q1.q[i]*q2.q[i];
			return dots;
		};
		
		static gquat& slerp( gquat& newQuat, double& unit_time, gquat& q1, gquat& q2 )
		{
			double angle = acos( dot( q1, q2 ) / (q1.norm() * q2.norm()) );
			gquat temp1 = q1;
			gquat temp2 = q2;
			
			temp1 *= sin( angle*(1 - unit_time) );
			temp2 *= sin( angle*unit_time );
			
			newQuat = temp1;
			newQuat += temp2;
			newQuat /= angle;
		
			return newQuat;
		};
		gquat& slerp( double& unit_time, gquat& q1, gquat& q2 )
		{
			
			double angle = acos( dot( q1, q2 ) / (q1.norm() * q2.norm()) );
			gquat temp1 = q1;
			gquat temp2 = q2;
			
			temp1 *= sin( angle*(1 - unit_time) );
			temp2 *= sin( angle*unit_time );
			
			(*this) = temp1;
			(*this) += temp2;
			(*this) /= angle;
		
			return (*this);
		};
		gquat& linterp( double& unit_time, gquat& q1, gquat& q2 )
		{
			gquat temp1 = q1;
			gquat temp2 = q2;
			
			temp1 *= ( 1 - unit_time );
			temp2 *= unit_time;
			
			(*this) = temp1;
			(*this) += temp2;
			(*this) /= (*this).norm();
		
			return (*this);
		};
		
		double norm() {
			return sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
		};
		
		gquat& operator=(gquat& that) {
			for (size_t i=0;i<4;i++) q[i] = that.q[i];
			return *this;
		};
		gquat& operator+=(gquat& that) {
			for (size_t i=0;i<4;i++) q[i] += that.q[i];
			return *this;
		};
		gquat& operator-=(gquat& that) {
			for (size_t i=0;i<4;i++) q[i] -= that.q[i];
			return *this;
		};
	
		gquat& operator/=(double scalar) {
			for (size_t i=0;i<4;i++) q[i] /= scalar;
			return *this;
		};
		
		gquat& operator*=(double scalar) {
			for (size_t i=0;i<4;i++) q[i] *= scalar;
			return *this;
		};
		
		gquat& set( double * vect, size_t step=1 ) {
			for(int i=0;i<=3; i++) q[i] = *(vect+i*step);
			return *this;
		};
		
        //left-handed rotation
		x3& rot(x3& v) {
			Rmatrix m;
			this->rm( m );
			m.rot( v );
			return v;
		};
        //right-handed rotation (inefficient)
        x3& trot(x3& v) {
			Rmatrix m;
			this->rm( m );
            m.tr();
			m.rot( v );
			return v;
		};
		
		Rmatrix& rm(Rmatrix& m) {
			//left-handed matrix form:
			double qx = q[0];
			double qy = q[1];
			double qz = q[2];
			double qs = q[3];
			
			m(0,0) = qx*qx - qy*qy - qz*qz + qs*qs;
			m(1,1) = -qx*qx + qy*qy - qz*qz + qs*qs;
			m(2,2) = -qx*qx - qy*qy + qz*qz + qs*qs;
			
			m(0,1) = 2 * ( qx*qy + qz*qs );
			m(0,2) = 2 * ( qx*qz - qy*qs );
			
			m(1,0) = 2 * ( qx*qy - qz*qs );
			m(1,2) = 2 * ( qx*qs + qy*qz );
			
			m(2,0) = 2 * ( qx*qz + qy*qs );
			m(2,1) = 2 * ( qy*qz - qx*qs );
			
			return m;
		};
		
		
	};
	
	
	

};

using namespace geom;

class GBM_Geometry {

	private:
	DynMatrix<double> normals;
	
	
	public:
	typedef DynMatrix<double>::size_type size_type;
	
	GBM_Geometry();
	
	
	const vector<double>& detNormal( size_type det );
	vector<double> copyDetNormal( size_type det );
	
	double ang2det(size_type det, double s_x, double s_y, double s_z, int cos=2);
	double ang2det(size_type det, vector<double>& s, int cos=2);
	double ang2det(size_type det, const vector<double>& s, int cos=2);
	double ang2det(size_type det, double s[2], int cos=2 );
	
		
};
class Rz : public Rmatrix {
	
	private:
	double angle;
	
	public:
	Rz(double angle) {
		this->set(angle);
	};
	
	Rz() { Rz(0.0); };
	
	Rz& set(double& angle) {
		this->angle = angle;
		(*this)(0,0) = cos(angle);
		(*this)(0,1) = sin(angle);
		(*this)(0,2) = 0.0;
		
		(*this)(1,0) = -(*this)(0,1);
		(*this)(1,1) = (*this)(0,0);
		(*this)(1,2) = 0.0;
		
		(*this)(2,0) = 0.0;
		(*this)(2,1) = 0.0;
		(*this)(2,2) = 1.0;
		return *this;
	};
	
	virtual Rz& operator=(Rmatrix& base) {
		this->Rmatrix::operator=(base);
		return *this;
	};
	
};

class Ry : public Rmatrix {
	
	private:
	double angle;
	
	public:
	Ry(double angle) {
		this->set(angle);
	};
	Ry() { Ry(0.0); };
	
	Ry& set(double& angle) {
		this->angle = angle;
		(*this)(0,0) = cos(angle);
		(*this)(0,1) = 0.0;
		(*this)(0,2) = -sin(angle);
		
		(*this)(1,0) = 0.0;
		(*this)(1,1) = 1.0;
		(*this)(1,2) = 0.0;
		
		(*this)(2,0) = sin(angle);
		(*this)(2,1) = 0.0;
		(*this)(2,2) = cos(angle);
		return *this;
	};
	
	virtual Ry& operator=(Rmatrix& base) {
		this->Rmatrix::operator=(base);
		return *this;
	};
	
};

class Rzyz : public Rmatrix {
	
	private:
	double a1;
	double a2;
	double a3;
	
	public:
	Rzyz() {
		
	
	};
	Rzyz(double z_rot, double yp_rot, double zpp_rot) {
		
		this->set(z_rot, yp_rot, zpp_rot);
	};
	
	Rzyz& set(double& z_rot, double& yp_rot, double& zpp_rot) {
		a1 = z_rot;
		a2 = yp_rot;
		a3 = zpp_rot;
		double c1 = cos(z_rot);
		double s1 = sin(z_rot);
		double c2 = cos(yp_rot);
		double s2 = sin(yp_rot);
		double c3 = cos(zpp_rot);
		double s3 = sin(zpp_rot);
		
		(*this)(0,0) = c1*c2*c3 - s1*s3;
		(*this)(0,1) = s1*c2*c3 + c1*s3;
		(*this)(0,2) = -s2*c3;
		
		(*this)(1,0) = -c1*c2*s3 - s1*c3;
		(*this)(1,1) = -s1*c2*s3 + c1*c3;
		(*this)(1,2) = s2*s3;
		
		(*this)(2,0) = c1*s2;
		(*this)(2,1) = s1*s2;
		(*this)(2,2) = c2;
		return *this;
	};
	
	virtual Rzyz& operator=(Rmatrix& base) {
		this->Rmatrix::operator=(base);
		return *this;
	};

};










#endif