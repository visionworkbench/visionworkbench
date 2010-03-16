#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vw/Image.h>
#include <vw/Image/PixelMath.h>
#include <vw/Image/PixelMask.h>
#include <vw/Image/MaskViews.h>
#include <vw/FileIO.h>
#include <math.h>
#include <time.h>

//using namespace std;
//using namespace vw;
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs=boost::filesystem;

#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <boost/math/distributions/chi_squared.hpp>

#include <vw/Core.h>
#include <vw/Image.h>
#include <vw/FileIO.h>
#include <vw/Cartography.h>
#include <vw/Math.h>
using namespace std;
using namespace vw;
using namespace vw::math;
using namespace vw::cartography;
#include <vw/Photometry/Reconstruct.h>
#include <vw/Photometry/Outlier.h>
#include <vw/Photometry/Exposure.h>
#include <vw/Photometry/Misc.h>

void weight_images(std::vector<std::string> output_files, std::vector<std::string> input_files) {
        for (unsigned i = 0; i < input_files.size(); ++i) {
                if (  fs::exists(output_files[i]) ) {
                        std::cout << " " << i << "th weight image: " << output_files[i] << " is preserved." << std::endl;
                } else {
                        std::cout << " " << i << "th uniform weight image: " << output_files[i] << " is created." << std::endl;
                        GeoReference geo;
                        read_georeference(geo, input_files[i]);
                        DiskImageView<PixelMask<PixelGray<uint8> > > image(input_files[i]);
                        ImageView<PixelMask<PixelGray<float> > > weight(image.cols(), image.rows());
                        weight = copy_mask(apply_mask(weight,0), image);
                        fill(weight,1);

                        write_georeferenced_image(output_files[i], weight, geo, TerminalProgressCallback("{Core}","Processing:"));
                }
        }
}

// Preserve output_files if output_files is modified later than input_files and exposure_times.
float save_weight_images(std::vector<std::string> output_files, std::vector<std::string> realexp_files,
                                                 std::vector<std::string> realavr_files, std::vector<std::string> realvar_files, char * exp_time_file) {
        GeoReference geo;
        using boost::math::chi_squared;

        struct stat exp_time_stat;
        Vector<float> exposure_times = load_exposure_times(exp_time_file, realavr_files.size());
        stat(exp_time_file, &exp_time_stat);
        struct stat output_stat, realexp_stat, realavr_stat, realvar_stat;
        for (int i = 0; i < output_files.size(); ++i) {
                read_georeference(geo, realexp_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > realexp(realexp_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > realavr(realavr_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > realvar(realvar_files[i]);
                ImageView<PixelMask<PixelGray<float> > > weight(realexp.cols(), realexp.rows());
                weight = copy_mask(apply_mask(weight,0), realexp);

                if ( !stat(output_files[i].c_str(), &output_stat) ) {
                        stat(realexp_files[i].c_str(), &realexp_stat);
                        stat(realavr_files[i].c_str(), &realavr_stat);
                        stat(realvar_files[i].c_str(), &realvar_stat);
                        if ( difftime(output_stat.st_mtime, exp_time_stat.st_mtime) > 0 && difftime(output_stat.st_mtime, realavr_stat.st_mtime) > 0  && difftime(output_stat.st_mtime, realvar_stat.st_mtime) > 0 ) {
                                std::cout << " " << i << "th weight image: " << output_files[i] << " is preserved." << std::endl;
                                continue;
                        }
                }

                std::cout << " Constructing " << i << "th weight image: " << output_files[i] << "." << std::endl;
                ImageView<PixelMask<PixelGray<float> > > ratio = 2*realavr/realvar;
                ImageView<PixelMask<PixelGray<float> > > radia = realexp/exposure_times(i);
                for (unsigned x=0; x<weight.cols(); ++x)
                        for (unsigned y=0; y<weight.rows(); ++y)
                                if ( is_valid(weight(x,y)) ) {
                                        if ( realvar(x,y) == 0 ) {
                                                std::cout << x << "," << y << ": p " << " w " << realvar(x,y) << " ";
                                                weight(x,y) = 1;
                                        } else {
                                                chi_squared dist(ratio(x,y)*realavr(x,y));
                                                float p = 2*cdf(dist, ratio(x,y)*radia(x,y));
                                                weight(x,y) = (p>1) ? 2-p: p;
                                                //                                                if ( weight(x,y) < 1.0/512 ) {
                                                //                                                        std::cout << x << "," << y << ": p " << p << ": z " << ratio(x,y)*radia(x,y) << " dof " << ratio(x,y)*realavr(x,y) << std::endl;
                                                //                                                        weight(x,y) = 1.0/512;
                                                //                                                }
                                        }
                                }

                write_georeferenced_image(output_files[i], weight, geo, TerminalProgressCallback("{Core}","Processing:"));
         }
        return 1;
}

#include <boost/math/distributions/normal.hpp> // for normal_distribution

// Preserve output_files if output_files is modified later than input_files and exposure_times.
float save_normal_images(std::vector<std::string> output_files, std::vector<std::string> realexp_files,
                                                 std::vector<std::string> realrad_files, char * exp_time_file) {

        GeoReference geo;
        using boost::math::normal; // typedef provides default type is float.
        using boost::math::gamma_distribution;

        struct stat output_stat, realexp_stat, realrad_stat, exp_time_stat;
        Vector<float> t = load_exposure_times(exp_time_file, realrad_files.size());
        stat(exp_time_file, &exp_time_stat);
        for (int i = 0; i < output_files.size(); ++i) {
                if ( !stat(output_files[i].c_str(), &output_stat) ) {
                        stat(realexp_files[i].c_str(), &realexp_stat);
                        stat(realrad_files[i].c_str(), &realrad_stat);
                        if ( difftime(output_stat.st_mtime, exp_time_stat.st_mtime) > 0 && difftime(output_stat.st_mtime, realrad_stat.st_mtime) > 0 ) {
                                std::cout << " " << i << "th normal image: " << output_files[i] << " is preserved." << std::endl;
                                continue;
                        }
                }

                std::cout << " Constructing " << i << "th normal image: " << output_files[i] << "." << std::endl;

                DiskImageView<PixelMask<PixelGray<float> > > R(realrad_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > X(realexp_files[i]);
                ImageView<PixelMask<PixelGray<float> > > Rt = R*t[i];
                ImageView<PixelMask<PixelGray<float> > > Z(X.cols(), X.rows());
                Z = copy_mask(apply_mask(Z,0), X);

                for (unsigned x=0; x < R.cols(); ++x)
                        for (unsigned y=0; y < R.rows(); ++y)
                                if ( is_valid(X(x,y)) ) {
                                        normal s;
                                        gamma_distribution<float> p(X(x,y)+1,1/Rt(x,y));
                                        float q = cdf(p,1);
                                        if ( 1-q < SIGNIFICANCE_LEVEL ) Z(x,y) = quantile(s, 1-SIGNIFICANCE_LEVEL);
                                        else if ( q < SIGNIFICANCE_LEVEL ) Z(x,y) = quantile(s, SIGNIFICANCE_LEVEL);
                                        else Z(x,y) = quantile(s, cdf(p, 1));
                                }
                read_georeference(geo, realexp_files[i]);

                write_georeferenced_image(output_files[i], Z, geo, TerminalProgressCallback("{Core}","Processing:"));
         }
        return 1;
}

std::vector<float> real_quadratic_roots(float a, float b, float c, float lb = 0, float ub = 1) {
        std::vector<float> roots;
        if ( a == 0 ) {
                if ( b != 0 ) {
                        float x = -c/b;
                        if (x > lb && x < ub) roots.push_back(x);
                }
        } else {
                float d = b*b-4*a*c;
                if ( d > 0 ) {
                        d = sqrt(d);
                        float x1 = ( -b + ( (b>0) ? -d : d ) ) / (2*a);
                        float x2 = c/a/x1;
                        if (x1 > lb && x1 < ub) roots.push_back(x1);
                        if (x2 > lb && x2 < ub) roots.push_back(x2);
                }
                if ( d == 0 ) {
                        float x = -b/2/a;
                        if (x > lb && x < ub) roots.push_back(x);
                }
        }

        if ( roots.size() > 0 ) std::cout << ".";
        return roots;
}

float weight_likelihood(Matrix<float> T, float w0, float e) {
        float n=T(0,0);
        float Tx=T(0,1);
        float Tx2=T(0,2);
        float Tw=T(1,0);
        float Twx=T(1,1);
        float T2=T(1,2);

        float b0 = Tw*T2+Tw*w0*e*e+w0*T2-Twx*Twx-2*Twx*w0*e;

        if ( b0 == 0 || Tw == 0 ) {
                std::cout << "Singular !!!" << std::endl;
                std::cout << "w0" << w0 << ", e" << e <<" ";
                return numeric_limits<float>::max();
        } else {
                float A0 = (Tw+w0)*(Tw+w0)*Tx2-2*Tx*(Tw+w0)*(Twx+w0*e)+n*(Twx+w0*e)*(Twx+w0*e);
                return(A0/b0*(1-1/Tw));
        }
}

// Preserve output_files if output_files is modified later than input_files and exposure_times.
Matrix<float> total_weighted_statistics(std::vector<std::string> output_files, std::vector<std::string> input_files) {
        Matrix<float> T(2,3);
        struct stat output_stat, input_stat;
        for (int i = 0; i < output_files.size(); ++i) {
                DiskImageView<PixelMask<PixelGray<float> > > X(input_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > W(output_files[i]);
                ImageView<PixelMask<PixelGray<float> > > WX = W*X;
                ImageView<PixelMask<PixelGray<float> > > X2 = X*X;
                ImageView<PixelMask<PixelGray<float> > > WX2 = W*X2;

                for (unsigned x=0; x < W.cols(); ++x)
                        for (unsigned y=0; y < W.rows(); ++y) {
                                ++T(0,0);
                                T(0,1) += X(x,y);
                                T(0,2) += X2(x,y);
                                T(1,0) += W(x,y);
                                T(1,1) += WX(x,y);
                                T(1,2) += WX2(x,y);
                        }
        }

        std::cout << "Exposure & Weight Statistics: " << T << std::endl;
        return T;
}


#include  <complex>

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^3 + p[1] x^2...+p[3]=0
// x=r[1][k] + i r[2][k]  k=1,...,3
// Assumes 0<arctan(x)<pi/2 for x>0

int CubicRoots( double p[5], double r[3][5] )
{
        double s,t,b,c,d;
        int k;

        if( p[0] != 1. )
        {
                for(k = 1; k < 4; k++ ) { p[k] = p[k]/p[0]; }
                p[0] = 1.;
        }
        s = p[1]/3.0;
        t = s*p[1];
        b = 0.5*( s*( t/1.5 - p[2] ) + p[3] );
        t = ( t - p[2] )/3.0;
        c = t*t*t;
        d = b*b - c;

        if( d >= 0. )
        {
                d = std::pow( (std::sqrt(d) + std::fabs(b) ), 1.0/3.0 );

                if( d != 0. )
                {
                        if( b > 0. ) { b = -d; }
                        else         { b =  d; }
                        c =  t/b;
                }
                d       =  std::sqrt(0.75)*(b - c);
                r[2][2] =  d;
                b       =  b + c;
                c       = -0.5*b-s;
                r[1][2] =  c;

                if( ( b > 0. &&  s <= 0. ) || ( b < 0. && s > 0. ) )
                {
                        r[1][1] =  c;
                        r[2][1] = -d;
                        r[1][3] =  b - s;
                        r[2][3] =  0;
                }
                else
                {
                        r[1][1] =  b - s;
                        r[2][1] =  0.;
                        r[1][3] =  c;
                        r[2][3] = -d;
                }
        }              // end of 2 equal or complex roots
        else
        {
                if( b == 0. ) { d =  std::atan(1.0)/1.5; }
                else          { d =  std::atan( std::sqrt(-d)/std::fabs(b) )/3.0; }

                if( b < 0. )  { b =  std::sqrt(t)*2.0; }
                else          { b = -2.0*std::sqrt(t); }

                c =  std::cos(d)*b;
                t = -std::sqrt(0.75)*std::sin(d)*b - 0.5*c;
                d = -t - c - s;
                c =  c - s;
                t =  t - s;

                if( std::fabs(c) > std::fabs(t) ) { r[1][3] = c; }
                else
                {
                        r[1][3] = t;
                        t       = c;
                }
                if( std::fabs(d) > std::fabs(t) ) { r[1][2] = d; }
                else
                {
                        r[1][2] = t;
                        t       = d;
                }
                r[1][1] = t;

                for(k = 1; k < 4; k++ ) { r[2][k] = 0.; }
        }
        return 0;
}

std::vector<double> RealCubicRoots(double p[5], double lb = 0, double ub = 1)
{
        int k;
        double r[3][5];
        double s,t,b,c,d;
        std::vector<double> roots;

        if( p[0] != 1. )
        {
                for(k = 1; k < 4; k++ ) { p[k] = p[k]/p[0]; }
                p[0] = 1.;
        }
        s = p[1]/3.0;
        t = s*p[1];
        b = 0.5*( s*( t/1.5 - p[2] ) + p[3] );
        t = ( t - p[2] )/3.0;
        c = t*t*t;
        d = b*b - c;

        if( d >= 0. )
        {
                d = std::pow( (std::sqrt(d) + std::fabs(b) ), 1.0/3.0 );

                if( d != 0. )
                {
                        if( b > 0. ) { b = -d; }
                        else         { b =  d; }
                        c =  t/b;
                }
                d       =  std::sqrt(0.75)*(b - c);
                r[2][2] =  d;
                b       =  b + c;
                c       = -0.5*b-s;
                r[1][2] =  c;

                if( ( b > 0. &&  s <= 0. ) || ( b < 0. && s > 0. ) )
                {
                        r[1][1] =  c;
                        r[2][1] = -d;
                        r[1][3] =  b - s;
                        r[2][3] =  0;
                }
                else
                {
                        r[1][1] =  b - s;
                        r[2][1] =  0.;
                        r[1][3] =  c;
                        r[2][3] = -d;
                }
        }              // end of 2 equal or complex roots
        else
        {
                if( b == 0. ) { d =  std::atan(1.0)/1.5; }
                else          { d =  std::atan( std::sqrt(-d)/std::fabs(b) )/3.0; }

                if( b < 0. )  { b =  std::sqrt(t)*2.0; }
                else          { b = -2.0*std::sqrt(t); }

                c =  std::cos(d)*b;
                t = -std::sqrt(0.75)*std::sin(d)*b - 0.5*c;
                d = -t - c - s;
                c =  c - s;
                t =  t - s;

                if( std::fabs(c) > std::fabs(t) ) { r[1][3] = c; }
                else
                {
                        r[1][3] = t;
                        t       = c;
                }
                if( std::fabs(d) > std::fabs(t) ) { r[1][2] = d; }
                else
                {
                        r[1][2] = t;
                        t       = d;
                }
                r[1][1] = t;

                for(k = 1; k < 4; k++ ) { r[2][k] = 0.; }
        }
        for (k = 1; k < 4; ++k)
                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) roots.push_back(r[1][k]);
        //        for (k = 1; k < 5; ++k)
        //                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) std::cout << k << " " << r[1][k] << " ";
        //        for (k = 0; k < roots.size(); ++k)
        //                std::cout << k << "r " << roots[k] << " ";
        return roots;
}

//////////////////////////////////////////////////////////////////////////////
//
// Array r[3][5]  p[5]
// Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
// x=r[1][k] + i r[2][k]  k=1,...,4

std::vector<double> QuarticRoots(double p[5], double lb = 0, double ub = 1)
{
        std::vector<double> roots;

        double r[3][5];
        double a0, a1, a2, a3, y1;
        double R2, D2, E2, D, E, R = 0.;
        double a, b, c, d, ds;

        double reRoot[4];
        int k, noRoots, noReRoots = 0;

        for( k = 0; k < 4; k++ ) { reRoot[k] = numeric_limits<double>::max(); }

        if( p[0] != 1.0 )
        {
                for( k = 1; k < 5; k++) { p[k] = p[k]/p[0]; }
                p[0] = 1.;
        }
        a3 = p[1];
        a2 = p[2];
        a1 = p[3];
        a0 = p[4];

        // resolvent cubic equation cofs:

        p[1] = -a2;
        p[2] = a1*a3 - 4*a0;
        p[3] = 4*a2*a0 - a1*a1 - a3*a3*a0;

        noRoots = CubicRoots(p,r);

        for( k = 1; k < 4; k++ )
        {
                if( r[2][k] == 0. ) // find a real root
                {
                        noReRoots++;
                        reRoot[k] = r[1][k];
                }
                else reRoot[k] = numeric_limits<double>::max(); // kInfinity;
        }
        y1 = numeric_limits<double>::max(); // kInfinity;
        for( k = 1; k < 4; k++ )
        {
                if ( reRoot[k] < y1 ) { y1 = reRoot[k]; }
        }

        R2 = 0.25*a3*a3 - a2 + y1;
        b  = 0.25*(4*a3*a2 - 8*a1 - a3*a3*a3);
        c  = 0.75*a3*a3 - 2*a2;
        a  = c - R2;
        d  = 4*y1*y1 - 16*a0;

        if( R2 > 0.)
        {
                R = std::sqrt(R2);
                D2 = a + b/R;
                E2 = a - b/R;

                if( D2 >= 0. )
                {
                        D       = std::sqrt(D2);
                        r[1][1] = -0.25*a3 + 0.5*R + 0.5*D;
                        r[1][2] = -0.25*a3 + 0.5*R - 0.5*D;
                        if ( r[1][1] > lb && r[1][1] < ub ) roots.push_back(r[1][1]);
                        if ( r[1][2] > lb && r[1][2] < ub ) roots.push_back(r[1][2]);
                }
                if( E2 >= 0. )
                {
                        E       = std::sqrt(E2);
                        r[1][3] = -0.25*a3 - 0.5*R + 0.5*E;
                        r[1][4] = -0.25*a3 - 0.5*R - 0.5*E;
                        if ( r[1][3] > lb && r[1][3] < ub ) roots.push_back(r[1][3]);
                        if ( r[1][4] > lb && r[1][4] < ub ) roots.push_back(r[1][4]);
                }
        }
        else if( R2 < 0.)
        {
                R = std::sqrt(-R2);
                complex<double> CD2(a,-b/R);
                complex<double> CD = std::sqrt(CD2);

                r[1][1] = -0.25*a3 + 0.5*real(CD);
                r[1][2] = -0.25*a3 - 0.5*real(CD);
                r[2][1] =  0.5*R + 0.5*imag(CD);
                r[2][2] =  0.5*R - 0.5*imag(CD);
                complex<double> CE2(a,b/R);
                complex<double> CE = std::sqrt(CE2);

                r[1][3] = -0.25*a3 + 0.5*real(CE);
                r[1][4] = -0.25*a3 - 0.5*real(CE);
                r[2][3] =  -0.5*R + 0.5*imag(CE);
                r[2][4] =  -0.5*R - 0.5*imag(CE);

                for (k = 1; k < 5; ++k)
                        if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) roots.push_back(r[1][k]);
        }
        else // R2=0 case
        {
                if(d >= 0.)
                {
                        D2 = c + std::sqrt(d);
                        E2 = c - std::sqrt(d);

                        if( D2 >= 0. )
                        {
                                D       = std::sqrt(D2);
                                r[1][1] = -0.25*a3 + 0.5*R + 0.5*D;
                                r[1][2] = -0.25*a3 + 0.5*R - 0.5*D;
                                if ( r[1][1] > lb && r[1][1] < ub ) roots.push_back(r[1][1]);
                                if ( r[1][2] > lb && r[1][2] < ub ) roots.push_back(r[1][2]);
                        }
                        if( E2 >= 0. )
                        {
                                E       = std::sqrt(E2);
                                r[1][3] = -0.25*a3 - 0.5*R + 0.5*E;
                                r[1][4] = -0.25*a3 - 0.5*R - 0.5*E;
                                if ( r[1][3] > lb && r[1][3] < ub ) roots.push_back(r[1][3]);
                                if ( r[1][4] > lb && r[1][4] < ub ) roots.push_back(r[1][4]);
                        }
                }
                else
                {
                        ds = std::sqrt(-d);
                        complex<double> CD2(c,ds);
                        complex<double> CD = std::sqrt(CD2);

                        r[1][1] = -0.25*a3 + 0.5*real(CD);
                        r[1][2] = -0.25*a3 - 0.5*real(CD);
                        r[2][1] =  0.5*R + 0.5*imag(CD);
                        r[2][2] =  0.5*R - 0.5*imag(CD);

                        complex<double> CE2(c,-ds);
                        complex<double> CE = std::sqrt(CE2);

                        r[1][3] = -0.25*a3 + 0.5*real(CE);
                        r[1][4] = -0.25*a3 - 0.5*real(CE);
                        r[2][3] =  -0.5*R + 0.5*imag(CE);
                        r[2][4] =  -0.5*R - 0.5*imag(CE);
                        for (k = 1; k < 5; ++k)
                                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) roots.push_back(r[1][k]);
                }
        }
        //        for (k = 1; k < 5; ++k)
        //                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) roots.push_back(r[1][k]);
        //        for (k = 1; k < 5; ++k)
        //                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) std::cout << k << " " << r[1][k] << " ";
        //        for (k = 0; k < roots.size(); ++k)
        //                std::cout << k << "r " << roots[k] << " ";
        return roots;
}

std::vector<double> RealQuarticRoots( double p[5], double lb = 0, double ub = 1)
{
        std::vector<double> roots;
        double r[3][5];
        double a0, a1, a2, a3, y1;
        double R2, D2, E2, D, E, R = 0.;
        double a, b, c, d, ds;

        double reRoot[4];
        int k, noRoots, noReRoots = 0;

        for( k = 0; k < 4; k++ ) { reRoot[k] = numeric_limits<double>::max(); }

        if( p[0] != 1.0 )
        {
                for( k = 1; k < 5; k++) { p[k] = p[k]/p[0]; }
                p[0] = 1.;
        }
        a3 = p[1];
        a2 = p[2];
        a1 = p[3];
        a0 = p[4];

        // resolvent cubic equation cofs:

        p[1] = -a2;
        p[2] = a1*a3 - 4*a0;
        p[3] = 4*a2*a0 - a1*a1 - a3*a3*a0;

        noRoots = CubicRoots(p,r);

        for( k = 1; k < 4; k++ )
        {
                if( r[2][k] == 0. ) // find a real root
                {
                        noReRoots++;
                        reRoot[k] = r[1][k];
                }
                else reRoot[k] = numeric_limits<double>::max(); // kInfinity;
        }
        y1 = numeric_limits<double>::max(); // kInfinity;
        for( k = 1; k < 4; k++ )
        {
                if ( reRoot[k] < y1 ) { y1 = reRoot[k]; }
        }

        R2 = 0.25*a3*a3 - a2 + y1;
        b  = 0.25*(4*a3*a2 - 8*a1 - a3*a3*a3);
        c  = 0.75*a3*a3 - 2*a2;
        a  = c - R2;
        d  = 4*y1*y1 - 16*a0;

        if( R2 > 0.)
        {
                R = std::sqrt(R2);
                D2 = a + b/R;
                E2 = a - b/R;

                if( D2 >= 0. )
                {
                        D       = std::sqrt(D2);
                        r[1][1] = -0.25*a3 + 0.5*R + 0.5*D;
                        r[1][2] = -0.25*a3 + 0.5*R - 0.5*D;
                        r[2][1] = 0.;
                        r[2][2] = 0.;
                }
                else
                {
                        D       = std::sqrt(-D2);
                        r[1][1] = -0.25*a3 + 0.5*R;
                        r[1][2] = -0.25*a3 + 0.5*R;
                        r[2][1] =  0.5*D;
                        r[2][2] = -0.5*D;
                }
                if( E2 >= 0. )
                {
                        E       = std::sqrt(E2);
                        r[1][3] = -0.25*a3 - 0.5*R + 0.5*E;
                        r[1][4] = -0.25*a3 - 0.5*R - 0.5*E;
                        r[2][3] = 0.;
                        r[2][4] = 0.;
                }
                else
                {
                        E       = std::sqrt(-E2);
                        r[1][3] = -0.25*a3 - 0.5*R;
                        r[1][4] = -0.25*a3 - 0.5*R;
                        r[2][3] =  0.5*E;
                        r[2][4] = -0.5*E;
                }
        }
        else if( R2 < 0.)
        {
                R = std::sqrt(-R2);
                std::complex<double> CD2(a,-b/R);
                std::complex<double> CD = std::sqrt(CD2);

                r[1][1] = -0.25*a3 + 0.5*real(CD);
                r[1][2] = -0.25*a3 - 0.5*real(CD);
                r[2][1] =  0.5*R + 0.5*imag(CD);
                r[2][2] =  0.5*R - 0.5*imag(CD);
                std::complex<double> CE2(a,b/R);
                std::complex<double> CE = std::sqrt(CE2);

                r[1][3] = -0.25*a3 + 0.5*real(CE);
                r[1][4] = -0.25*a3 - 0.5*real(CE);
                r[2][3] =  -0.5*R + 0.5*imag(CE);
                r[2][4] =  -0.5*R - 0.5*imag(CE);
        }
        else // R2=0 case
        {
                if(d >= 0.)
                {
                        D2 = c + std::sqrt(d);
                        E2 = c - std::sqrt(d);

                        if( D2 >= 0. )
                        {
                                D       = std::sqrt(D2);
                                r[1][1] = -0.25*a3 + 0.5*R + 0.5*D;
                                r[1][2] = -0.25*a3 + 0.5*R - 0.5*D;
                                r[2][1] = 0.;
                                r[2][2] = 0.;
                        }
                        else
                        {
                                D       = std::sqrt(-D2);
                                r[1][1] = -0.25*a3 + 0.5*R;
                                r[1][2] = -0.25*a3 + 0.5*R;
                                r[2][1] =  0.5*D;
                                r[2][2] = -0.5*D;
                        }
                        if( E2 >= 0. )
                        {
                                E       = std::sqrt(E2);
                                r[1][3] = -0.25*a3 - 0.5*R + 0.5*E;
                                r[1][4] = -0.25*a3 - 0.5*R - 0.5*E;
                                r[2][3] = 0.;
                                r[2][4] = 0.;
                        }
                        else
                        {
                                E       = std::sqrt(-E2);
                                r[1][3] = -0.25*a3 - 0.5*R;
                                r[1][4] = -0.25*a3 - 0.5*R;
                                r[2][3] =  0.5*E;
                                r[2][4] = -0.5*E;
                        }
                }
                else
                {
                        ds = std::sqrt(-d);
                        std::complex<double> CD2(c,ds);
                        std::complex<double> CD = std::sqrt(CD2);

                        r[1][1] = -0.25*a3 + 0.5*real(CD);
                        r[1][2] = -0.25*a3 - 0.5*real(CD);
                        r[2][1] =  0.5*R + 0.5*imag(CD);
                        r[2][2] =  0.5*R - 0.5*imag(CD);

                        std::complex<double> CE2(c,-ds);
                        std::complex<double> CE = std::sqrt(CE2);

                        r[1][3] = -0.25*a3 + 0.5*real(CE);
                        r[1][4] = -0.25*a3 - 0.5*real(CE);
                        r[2][3] =  -0.5*R + 0.5*imag(CE);
                        r[2][4] =  -0.5*R - 0.5*imag(CE);
                }
        }
        for (k = 1; k < 5; ++k)
                if ( r[2][k] == 0 && r[1][k] > lb && r[1][k] < ub ) roots.push_back(r[1][k]);
        return roots;
}

//
//
//////////////////////////////////////////////////////////////////////////////

/*
 float optimal_weight(Matrix<float> T, float w, float x) {
 float n=T(0,0);
 float tx1=T(0,1);
 float tx2=T(0,2);
 float sw1=T(1,0);
 float swx=T(1,1);
 float wx2=T(1,2);

 float x2 = x*x;
 float x3 = x*x2;
 float x4 = x*x3;

 float sw12 = sw1*sw1;
 float sw13 = sw1*sw12;
 float sw14 = sw1*sw13;

 float swx2 = swx*swx;
 float swx3 = swx*swx2;

 float d0 = (tx2*sw12-2*tx1*sw1*swx+n*swx2+(sw1-1)*(2*tx2*sw1-2*tx1*swx-2*tx1*x*sw1+2*n*x*swx))*sw1*(wx2*sw1-swx2)-(sw1-1)*(tx2*sw12-2*tx1*sw1*swx+n*swx2)*(wx2*sw1-swx2+sw1*(x2*sw1+wx2-2*x*swx));
 float d1 = 4*tx1*sw1*swx3-6*tx2*sw12*swx2+4*tx2*sw13*wx2+4*tx1*sw13*swx*x2-4*n*swx2*x2*sw12-2*n*swx2*wx2*sw1-2*tx2*sw14*x2+4*tx2*sw13*x*swx+4*sw12*n*x*swx*wx2-4*tx1*swx*wx2*sw1+4*tx1*x*sw12*wx2+4*tx1*x*sw1*swx2+2*tx2*sw13*x2+2*n*swx2*wx2-4*tx2*sw12*x*swx-4*tx1*sw12*swx*x2+2*tx2*sw1*swx2-4*n*x*swx3-2*sw12*n*x2*wx2+4*sw1*n*x2*swx2+2*sw13*n*x2*wx2-8*tx1*x*sw13*wx2;
 float d2 = 6*tx2*wx2*sw12-3*n*x2*swx2+2*tx1*x*swx2-6*sw1*tx2*swx2+tx2*sw12*x2-n*swx2*wx2+sw13*n*x4-2*tx1*x3*sw13-2*tx1*swx*wx2-sw12*n*x4-5*sw1*n*x2*swx2+6*sw1*tx1*x*swx2+5*n*x2*wx2*sw12-12*tx1*x*wx2*sw12-2*n*x2*wx2*sw1+2*tx1*x*wx2*sw1+4*n*x3*swx*sw1-2*tx1*swx*x2*sw1-2*tx2*sw1*x*swx+6*tx1*swx*x2*sw12+2*n*x*swx*wx2-2*sw12*n*x3*swx+2*tx1*swx3+tx2*swx2+2*n*x*swx*wx2*sw1;
 float d3 = 2*(tx2-2*tx1*x+n*x2)*(2*wx2*sw1-swx2+x2*sw12-2*sw1*x*swx);
 float d4 = (tx2-2*tx1*x+n*x2)*(x2*sw1+wx2-2*x*swx);

 double p[5] = {d4, d3, d2, d1, d0};
 float lb = (sw1 < 1) ? 1-sw1: MIN_WEIGHT;
 std::vector<double> roots = RealQuarticRoots(p,lb);
 roots.push_back(1);
 roots.push_back(lb);

 float max_w, max_r = -numeric_limits<float>::max();
 for (unsigned k=0; k < roots.size(); k++) {
 w = roots[k];

 float tw1 = w+sw1;
 float twx = w*x+swx;
 float b = tw1*(tw1*(w*x2+wx2)-twx*twx);
 if ( b != 0 ) {
 float a = (tw1-1)*(tw1*tw1*tx2-2*tw1*tx1*twx+n*twx*twx);
 if ( max_r < (a/b) ) {
 max_w = w;
 }
 }
 }
 if ( max_w > lb  && max_w < 1 ) std::cout << " w: " << max_w << " ";

 return max_w;
 }
 */
#include <boost/math/distributions/fisher_f.hpp>

float optimal_weight(Matrix<float> T, float w, float x) {
        using boost::math::fisher_f;

        float n=T(0,0);
        float tx1=T(0,1);
        float tx2=T(0,2);
        float sw1=T(1,0);
        float swx=T(1,1);
        float wx2=T(1,2);

        float x2 = x*x;
        float x3 = x*x2;
        float x4 = x*x3;

        float sw12 = sw1*sw1;
        float sw13 = sw1*sw12;
        float sw14 = sw1*sw13;

        float swx2 = swx*swx;
        float swx3 = swx*swx2;

        float d0 = (tx2*sw12-2*tx1*sw1*swx+n*swx2+(sw1-1)*(2*tx2*sw1-2*tx1*swx-2*tx1*x*sw1+2*n*x*swx))*sw1*(wx2*sw1-swx2)-(sw1-1)*(tx2*sw12-2*tx1*sw1*swx+n*swx2)*(wx2*sw1-swx2+sw1*(x2*sw1+wx2-2*x*swx));
        float d1 = 4*tx1*sw1*swx3-6*tx2*sw12*swx2+4*tx2*sw13*wx2+4*tx1*sw13*swx*x2-4*n*swx2*x2*sw12-2*n*swx2*wx2*sw1-2*tx2*sw14*x2+4*tx2*sw13*x*swx+4*sw12*n*x*swx*wx2-4*tx1*swx*wx2*sw1+4*tx1*x*sw12*wx2+4*tx1*x*sw1*swx2+2*tx2*sw13*x2+2*n*swx2*wx2-4*tx2*sw12*x*swx-4*tx1*sw12*swx*x2+2*tx2*sw1*swx2-4*n*x*swx3-2*sw12*n*x2*wx2+4*sw1*n*x2*swx2+2*sw13*n*x2*wx2-8*tx1*x*sw13*wx2;
        float d2 = 6*tx2*wx2*sw12-3*n*x2*swx2+2*tx1*x*swx2-6*sw1*tx2*swx2+tx2*sw12*x2-n*swx2*wx2+sw13*n*x4-2*tx1*x3*sw13-2*tx1*swx*wx2-sw12*n*x4-5*sw1*n*x2*swx2+6*sw1*tx1*x*swx2+5*n*x2*wx2*sw12-12*tx1*x*wx2*sw12-2*n*x2*wx2*sw1+2*tx1*x*wx2*sw1+4*n*x3*swx*sw1-2*tx1*swx*x2*sw1-2*tx2*sw1*x*swx+6*tx1*swx*x2*sw12+2*n*x*swx*wx2-2*sw12*n*x3*swx+2*tx1*swx3+tx2*swx2+2*n*x*swx*wx2*sw1;
        float d3 = 2*(tx2-2*tx1*x+n*x2)*(2*wx2*sw1-swx2+x2*sw12-2*sw1*x*swx);
        float d4 = (tx2-2*tx1*x+n*x2)*(x2*sw1+wx2-2*x*swx);

        double pq[5] = {d4, d3, d2, d1, d0};
        float lb = (sw1 < 1) ? 1-sw1: MIN_WEIGHT;
        if ( w+sw1-1 < 1 ) {
                //                std::cout << " w " << w << " < lb" << lb << " return 1 ";
                w = 1;
                return w;
        }
        std::vector<double> quartic_roots = RealQuarticRoots(pq,lb);
        quartic_roots.push_back(1);
        quartic_roots.push_back(lb);

        float alpha = 0.95;
        float f = quantile(fisher_f(n, w+sw1-1), alpha);
        //        std::cout << " f(" << n << "," << w+sw1-1 << ")=" << f << " ";

        d0 = (sw1-1)*(tx2*sw12-2*tx1*sw1*swx+n*swx2)-n*f*sw1*(-swx2+wx2*sw1);
        d1 = -n*f*sw12*x2+((sw1-1)*(-2*sw1*tx1+2*n*swx)+2*n*f*sw1*swx)*x+tx2*sw12-2*tx1*sw1*swx+n*swx2+(sw1-1)*(2*sw1*tx2-2*tx1*swx)-n*f*(-swx2+wx2*sw1)-n*f*wx2*sw1;
        d2 = (-n+sw1*n-n*f*sw1)*x2+(-4*sw1*tx1+2*n*swx+2*tx1+2*n*f*swx)*x+3*sw1*tx2-2*tx1*swx-n*f*wx2-tx2;
         d3 = tx2-2*tx1*x+n*x2;

        double pc[5] = {d3, d2, d1, d0, 0};
        std::vector<double> cubic_roots = RealCubicRoots(pc,lb);
        for (unsigned k=0; k < quartic_roots.size(); k++) {
                float w = quartic_roots[k];

                float tw1 = w+sw1;
                float twx = w*x+swx;
                float b = tw1*(tw1*(w*x2+wx2)-twx*twx);
                if ( b != 0 ) {
                        float a = (tw1-1)*(tw1*tw1*tx2-2*tw1*tx1*twx+n*twx*twx);
                        if ( (a/b) > n*f ) cubic_roots.push_back(w);
                }
        }

        float min_w, min_r = numeric_limits<float>::max();
        if ( cubic_roots.size() > 0 ) {
                for (unsigned k=0; k < cubic_roots.size(); k++) {
                        float w = cubic_roots[k];

                        float tw1 = w+sw1;
                        float twx = w*x+swx;
                        float b = tw1*(tw1-1);
                        if ( b != 0 ) {
                                float a = (w*x2+wx2)-twx*twx;
                                float r = a/b;
                                if ( min_r > r ) {
                                        min_r = r;
                                        min_w = w;
                                }
                        }
                }
        } else {
                min_w = w;
                //                std::cout << " no root: " << w << " ";
        }
        return min_w;
}

// Preserve output_files if output_files is modified later than input_files and exposure_times.
void save_inlier_images(std::vector<std::string> output_files, std::vector<std::string> input_files, std::vector<std::string> index_files, char * weight_file) {
        struct stat output_stat, input_stat;
        Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);

        for (int i = 0; i < output_files.size(); ++i) {
                std::cout << " Constructing " << i << "th weight image: " << output_files[i] << "." << std::endl;

                GeoReference geo;
                read_georeference(geo, input_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > X(input_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > W(output_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);

                ImageView<PixelMask<PixelGray<float> > > TX1=X;
                ImageView<PixelMask<PixelGray<float> > > TX2=X*X;
                ImageView<PixelMask<PixelGray<float> > > SW1(index.cols(),index.rows());
                ImageView<PixelMask<PixelGray<float> > > SWX(index.cols(),index.rows());
                ImageView<PixelMask<PixelGray<float> > > WX2(index.cols(),index.rows());
                SW1 = copy_mask(apply_mask(SW1,0), X);
                SWX = copy_mask(apply_mask(SWX,0), X);
                WX2 = copy_mask(apply_mask(WX2,0), X);
                fill(SW1,0);
                fill(SWX,0);
                fill(WX2,0);
                ImageView<PixelMask<PixelGray<float> > > V=W;

                int jmin = (i>4) ? i-4: 0;
                for (int j = jmin; j<i; ++j) {
                        uint8 mask = 0x01 << j-i+4;
                        ImageView<PixelMask<PixelGray<float> > > X = interpolate_image(input_files[j], index_files[i], mask);
                        ImageView<PixelMask<PixelGray<float> > > W = interpolate_image(output_files[j], index_files[i], mask);
                        TX1 += X;
                        TX2 += X*X;
                        SW1 += W;
                        SWX += W*X;
                        WX2 += SWX*X;
                }

                int jmax = ((i+5)<input_files.size()) ? i+5: input_files.size();
                for (int j = i+1; j<jmax; ++j) {
                        uint8 mask = 0x10 << j-i-1;
                        ImageView<PixelMask<PixelGray<float> > > X = interpolate_image(input_files[j], index_files[i], mask);
                        ImageView<PixelMask<PixelGray<float> > > W = interpolate_image(output_files[j], index_files[i], mask);
                        TX1 += X;
                        TX2 += X*X;
                        SW1 += W;
                        SWX += W*X;
                        WX2 += SWX*X;
                }

                for (unsigned x=0; x<index.cols(); ++x)
                        for (unsigned y=0; y<index.rows(); ++y) {
                                if ( is_valid(index(x,y)) ) {
                                        if ( inverse_weight(index(x,y)) < 3 ) {
                                                V(x,y) = 1;
                                        } else {
                                                float some_data[6] = {inverse_weight(index(x,y)), TX1(x,y), TX2(x,y), SW1(x,y), SWX(x,y), WX2(x,y)};
                                                Matrix<float, 2, 3> T(some_data);
                                                //                                        std::cout << T << std::endl;
                                                V(x,y) = optimal_weight(T,W(x,y),X(x,y));
                                        }
                                }
                        }

                write_georeferenced_image(output_files[i], V, geo, TerminalProgressCallback("{Core}","Processing:"));
         }
}


// Preserve output_files if output_files is modified later than input_files and exposure_times.
void save_outlier_images(std::vector<std::string> output_files, std::vector<std::string> input_files, std::vector<std::string> index_files, char * weight_file) {
        struct stat output_stat, input_stat;
        Vector<uint8> inverse_weight = load_inverse_weight(weight_file, 256);
        Matrix<float> T = total_weighted_statistics(output_files, input_files);

        float n=T(0,0);
        float Tx=T(0,1);
        float Tx2=T(0,2);
        float Tw=T(1,0);
        float Twx=T(1,1);
        float T2=T(1,2);


        for (int i = 0; i < output_files.size(); ++i) {
                std::cout << " Constructing " << i << "th weight image: " << output_files[i] << "." << std::endl;

                GeoReference geo;
                read_georeference(geo, input_files[i]);
                DiskImageView<PixelMask<PixelGray<uint8> > > index(index_files[i]);
                DiskImageView<PixelMask<PixelGray<float> > > X(input_files[i]);
                ImageView<PixelMask<PixelGray<float> > > W(X.cols(),X.rows());
                W = copy_mask(apply_mask(W,0), X);
                float lb;
                for (unsigned x=0; x<X.cols(); ++x)
                        for (unsigned y=0; y<X.rows(); ++y) {
                                if ( is_valid(X(x,y)) ) {
                                        if ( inverse_weight(index(x,y)) < 3 ) {
                                                W(x,y) = 1;
                                        } else {
                                                float e = X(x,y);
                                                float e2 = e*e;
                                                float Tw2=Tw*Tw;
                                                float Tw3=Tw*Tw2;
                                                float Twx2=Twx*Twx;
                                                float Twx3=Twx*Twx2;
                                                float D0 = (-Tx2*Tw3+2*Twx*Tx*Tw2-n*Twx2*Tw)*e2+(-2*T2*Tx*Tw2-2*Tx*Twx2*Tw+2*Twx*Tx2*Tw2+2*n*T2*Twx*Tw)*x+Tw2*Tx2*T2-2*Tw*Tx2*Twx2+2*Twx3*Tx-n*Twx2*T2;
                                                float D1 = 2*(Tw*T2-Twx2)*(Tx2-2*x*Tx+n*e2);
                                                float D2 = (Tx2-2*e*Tx+n*e2)*(T2-2*e*Twx+e2*Tw);
                                                //                                                if ( D0 == 0 ) std::cout << "(" << x << "," << y << ") " << D0 << " " << D1 << " " << D2 << " ";
                                                std::vector<float> roots = real_quadratic_roots(D2,D1,D0,MIN_WEIGHT);
                                                //                                                std::cout << roots.size() << " ";
                                                roots.push_back(MIN_WEIGHT);
                                                W(x,y) = 1;
                                                float max_likelihood = weight_likelihood(T,1,e);
                                                for (unsigned k=0; k < roots.size(); k++)
                                                        if ( max_likelihood < weight_likelihood(T,roots[k],e) ) W(x,y) = roots[k];
                                        }
                                }
                        }
                write_georeferenced_image(output_files[i], W, geo, TerminalProgressCallback("{Core}","Processing:"));
         }
}


