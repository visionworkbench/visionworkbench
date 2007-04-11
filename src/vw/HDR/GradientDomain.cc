
// STL
#include <vector>
#include <list>
#include <sstream>
#include <iostream>

// Vision Workbench
#include <vw/ImageView.h>
#include <vw/ImageViewBase.h>
#include <vw/FileIO.h>
#include <vw/Filter.h>
#include <vw/ImageAlgorithms.h>
#include <vw/ImageMath.h>
#include <vw/ImageStatistics.h>
#include <vw/Transform.h>

using namespace std;

namespace vw { namespace HDR {

		//  Solve a PDE using successive overrelaxation
		void SuccessiveOverrelaxation( ImageView<double> &U, const ImageView<double> &RHS, unsigned iterations) {
      
			cout << "Solving using successive overrelaxation: " << endl;
      
			float anormf = 0.0, alpha = 1.0;
			const float EPSILON = 1.0e-12;
			
			// The asymptotic spectral radius of the gauss-seidel method is
			// approximately 1-pi^2/N^2 for a square NxN grid.  (This was
			// pulled from "The Nature of Mathematical Modeling" by Neil
			// Gershenfeld.
			int32 x_max = RHS.rows();
			int32 y_max = RHS.cols(); 
			float rho_gs = 1.0 - 6.28/((x_max<y_max) ? x_max : y_max);
			
			// Compute initial norm of residual and terminate iteration when
			// norm has been reduced by a factor EPSILON.
			fill(U,0.0);
			for (int32 i = 0; i < RHS.cols(); i++) {
				for (int32 j = 0; j < RHS.rows(); j++) {
					anormf += fabs( RHS(i,j) );
				}
			}
			
			//Assumes initial u is zero.
			for (unsigned n = 0; n < iterations; n++) {
				float anorm = 0.0;
				int jsw = 0;
				for (unsigned ipass = 0; ipass < 2; ipass++, jsw=1-jsw) {
					int isw = jsw;
					for (int32 j = 0; j < RHS.rows(); j++, isw = 1-isw) {
						for (int32 i = isw; i < RHS.cols(); i += 2) {
							int w, n, e, s;
							w = (i == 0 ? 0 : i-1);
							n = (j == 0 ? 0 : j-1);
							s = (j+1 == RHS.rows() ? j : j+1);
							e = (i+1 == RHS.cols() ? i : i+1);
							
							double resid = U(e,j) + U(w,j) + U(i,n) + U(i,s) - 4.0 * U(i,j) + RHS(i,j);
							anorm += fabs(resid);
							U(i,j) -= alpha * resid / -4.0;
						}
					}
					// The overrelaxation factor is computed using Chebychev
					// acceleration according to the technique in Numercial
					// Methods.
					alpha = ( n==0 && ipass==0 ? 
										1.0 / (1.0 - 0.5 * rho_gs * rho_gs) :
										1.0 / (1.0 - 0.25 * rho_gs * rho_gs * alpha));
				}
				
				if( n%50 == 0 || n==1)
					cout << "\tSOR -- Iteration: " << n << "\tAnorm: " << anorm << "  Alpha: " << alpha << "\n";
				if (anorm < EPSILON * anormf ) {
					cout << "\tSOR -- solved.\n";
					return;
				}
			}
			cout << "\tSOR -- ITERATIONS exceeded\n";
		}
		
		void relax_jacobi( ImageView<double> &U, const ImageView<double> &RHS ) {
			
			double h = 1.0 / (RHS.rows() - 1);
			double h2 = h*h;
			
			for( int32 i=0 ; i<RHS.cols() ; i++ ) {
				for( int32 j=0 ; j<RHS.rows() ; j++) {
					int w, n, e, s;
					w = (i == 0 ? 0 : i-1);
					n = (j == 0 ? 0 : j-1);
					s = (j+1 == RHS.rows() ? j : j+1);
					e = (i+1 == RHS.cols() ? i : i+1);
          
					// Note that the sign of h2*RHS(i,j) depends on the sign
					// of the laplacian that was used to compute the RHS
					// image.  Change the final term in this series from
					// positive to negative to accomodate a change in the sign
					// of the laplacian.
					U(i,j) = 0.25*( U(e,j) + U(w,j) + U(i,n) + U(i,s) + h2*RHS(i,j) );
				}
			} 
		}
		
		// Solve the boundary value problem using Gauss-seidel red/black
		// relaxation.
		void relax_gauss_seidel( ImageView<double> &U, const ImageView<double> &RHS ) {
			
			int isw, jsw = 0;
			
			//    double h = 1.0 / (RHS.rows() - 1);
			double h = 1.0;
			double h2 = h*h;

			for ( unsigned ipass = 0; ipass < 2; ipass++, jsw=1-jsw ) {
				isw = jsw;
				for( int32 j=0 ; j<RHS.rows() ; j++, isw = 1-isw) {
					for( int32 i=isw ; i<RHS.cols() ; i+=2 ) {
						int w, n, e, s;
						w = (i == 0 ? 0 : i-1);
						n = (j == 0 ? 0 : j-1);
						s = (j+1 == RHS.rows() ? j : j+1);
						e = (i+1 == RHS.cols() ? i : i+1);
          
						U(i,j) = 0.25*( U(e,j) + U(w,j) + U(i,n) + U(i,s) + h2*RHS(i,j) );
					}
				}
			} 
		}  
  
		void compute_residual( ImageView<double> &res, ImageView<double> &U, const ImageView<double> &RHS ) {

			res.set_size(RHS.cols(), RHS.rows());

			int n = U.rows();
			double h = 1.0 / (n - 1);
			double h2i = 1 / (h*h);

			for (int32 i = 0;i < RHS.cols();i++) {
				for (int32 j = 0;j < RHS.rows();j++) {
					int w, n, e, s;
					w = (i == 0 ? 0 : i-1);
					n = (j == 0 ? 0 : j-1);
					s = (j+1 == RHS.rows() ? j : j+1);
					e = (i+1 == RHS.cols() ? i : i+1);
	
					res(i,j) = h2i*(U(e,j)+U(w,j)+U(i,n)+U(i,s)-4.0*U(i,j))-RHS(i,j);
				}
			}
		}

		void MGV( ImageView<double> U, const ImageView<double> RHS, int j) {
			ostringstream stream;
			stream << "level-" << j;
    
			cout << "MGV: Level " << j << "\n";
			const int NPRE=5, NPOST=5;
			int jpost, jpre;

			if (j == 0) {
				// Compute the "exact" solution if we have reached the coarsest
				// level.
				cout << "--> computing base level... ";
				for (int n = 0; n < 100; n++)
					relax_gauss_seidel(U, RHS);
				cout << "done.\n";
			} else {
				ImageView<double> residual;

				for (jpre=0; jpre<NPRE;jpre++) {
					relax_gauss_seidel(U,RHS);
					compute_residual(residual,U,RHS);
					//	cout << "PRE " << jpre << ": Mean residual: " << mean_channel_value(residual) << "\n";
				}
				residual = -1 * resample(residual, 0.5, int(U.cols()), int(U.rows()));   // Restrict and change sign

				// V is initially zero at each iteration
				ImageView<double> v(residual.cols(), residual.rows());
				fill(v, 0.0);
      
				// Recurse: compute the error: v.
				MGV(v,residual,j-1);

				// Add in the error term
				U += v;

				for (jpost=0; jpost<NPOST;jpost++) {
					relax_gauss_seidel(U,RHS);
				}

			}
		}    

#define MIN_MULTIGRID_DIMENSION 128

		void MultiGrid( ImageView<double> &u, const ImageView<double> &rhs) {
    
			cout << "MultiGrid: Starting\n";
    
			// Compute the depth of the multigrid and compute the subsampled
			// versions of the right hand side of the poisson equation. 
			int nn = rhs.rows() < rhs.cols() ? rhs.rows() : rhs.cols();
			int nlevels = 1;
			while( nn > MIN_MULTIGRID_DIMENSION ) { 
				nn >>= 1;
				nlevels++;
			}
			cout << "MultiGrid: Using " << nlevels << " levels.\n";

			// The pyramid is built such that the first element is the most
			// subsampled (coarse) and the last element is the least
			// subsampled (fine).
			std::vector<ImageView<double> > rho(nlevels);
			rho[nlevels - 1] = copy(rhs);
			for (int n = nlevels-2; n >= 0; n--) {
				rho[n] = resample(rho[n+1], 0.5);  // restrict
			}
        
			// Solve the bottom level of the pyramid "exactly"
			ImageView<double> uj(rho[0].cols(), rho[0].rows());
			for (int n = 0; n < 100; n++) {
				if (n%20 == 0) { std::cout << "\tRelax: " << n << "\n"; }
				relax_gauss_seidel(uj, rho[0]);
			}
			ImageView<double> temp = normalize(uj);
			write_image("ubottom2.jpg", temp);
      
			for (int j = 1; j<nlevels; j++) {
				cout << "MGV: Processing level " << j << "\n";
				uj = resample(uj, 2.0);            // interpolate
				for (int jcycle = 0; jcycle < 1; jcycle++) {
					MGV(uj, rho[j], j);
				}
			}
			u = uj;
			cout << "MultiGrid: Complete.\n";
		}  


		class LuminanceGradientPyramid {
			std::vector<std::pair<ImageView<double>, ImageView<double> > > _levels;
			ImageView<double> _x_gradient;
			ImageView<double> _y_gradient;
			ImageView<double> _original;
			

		public:
			
			LuminanceGradientPyramid(ImageView<double> base_image, double sigma = 2.0, unsigned smallest_dimension = 32) {
				
				std::vector<double> identity_kernel(1);                                identity_kernel[0] = 1;
				std::vector<double> central_difference_kernel(3);
				central_difference_kernel[0] = -1;  central_difference_kernel[1] = 0;  central_difference_kernel[2] = 1;
				
				std::cout << "Exp -- Min: " << min_channel_value(base_image) << "   "
									<< "Max: " << max_channel_value(base_image) << "\n";
				ImageView<double> level_image = base_image;
				_original = level_image;
				
				std::cout << "Log -- Min: " << min_channel_value(_original) << "   "
									<< "Max: " << max_channel_value(_original) << "\n";
				
				// Build the gaussian pyramid and compute the gradient at each
				// level.  The gradient is computed using a central difference
				// kernel that is weighted according to the subsampling level.
				std::cout << "Building Gaussian Gradient Pyramid.\n";
				int level = 0;
				while ((level_image.cols() >= smallest_dimension) && (level_image.rows() >= smallest_dimension)) {
					std::cout << "\tLevel " << level << ": [" << level_image.cols() << " x " << level_image.rows() << "]\n";
					ImageView<double> x_gradient = separable_convolution_filter( level_image, 
																																			 central_difference_kernel, 
																																			 identity_kernel ) / pow(2,level+1);
					ImageView<double> y_gradient = separable_convolution_filter( level_image, 
																																			 identity_kernel,
																																			 central_difference_kernel) / pow(2,level+1);
					level++;
					
					// Add the gradient to the vector, blur, and then subsample by
					// a factor of 2.
					_levels.push_back(std::pair<ImageView<double>, ImageView<double> >(x_gradient, y_gradient));
					level_image = subsample(gaussian_filter(level_image, sigma), 2);
				}

				_x_gradient = _levels[0].first;
				_y_gradient = _levels[0].second;
			}
			
			
			
			void attenuate_gradient() {
				std::cout << "Computing attenuation factor.\n";
				
				ImageView<double> x_gradient = _levels[_levels.size()-1].first;
				ImageView<double> y_gradient = _levels[_levels.size()-1].second;	
				ImageView<double> magnitude = hypot(x_gradient, y_gradient);
				
				// These two factors affect the strength of the weighting (see
				// Section 4 of the paper referred to in the comments at the
				// top of this file for more information on these factors).
				// Although these are hardwired for now, they should
				// eventually be broken out as function arguments.
				double beta = 0.98;
				double alpha = 0.01 * mean_channel_value(magnitude);
				std::cout << "Alpha: " << alpha << ".\n";
				ImageView<double> attenuation = alpha / magnitude * pow(magnitude / alpha, beta); 

				for (int i = _levels.size() - 2; i >= 0; i--) {
					x_gradient = _levels[i].first;
					y_gradient = _levels[i].second;
					magnitude = hypot(x_gradient, y_gradient);
					
					ImageView<double> phi = alpha / magnitude * pow(magnitude / alpha, beta); 
					attenuation = phi * resample(attenuation, 2.0f, int(phi.cols()), int(phi.rows()));
					
					//	std::ostringstream filename_prefix;
					//	filename_prefix << "attenuation-level-" << i;
					// 	write_image(filename_prefix.str() + ".jpg", normalize(-1*attenuation));
					// 	write_image(filename_prefix.str() + "-mag.jpg", normalize(-1*magnitude));
				}
				std::cout << "Appyling attenuation.\n";
				_x_gradient *= attenuation;
				_y_gradient *= attenuation;
			}

			
			ImageView<double> reconstruct_image() {
				
				// Compute the divergence of the gradient, div G = dG/dx + dG/dy
				//  using a backward difference approximation.
				ImageView<double> divG(_x_gradient.cols(), _x_gradient.rows());
				for( int32 i = 0 ; i < _x_gradient.cols() ; i++ ) {
					for( int32 j = 0 ; j < _y_gradient.rows() ; j++ ) {
						divG(i,j) = _x_gradient(i,j) + _y_gradient(i,j);
						if( i > 0 ) divG(i,j) -= _x_gradient(i-1,j  );
						if( j > 0 ) divG(i,j) -= _y_gradient(i,  j-1);
					}
				}
				
				// Now, solve the poisson equation to find an original image
				// that matches our newly scaled gradient function as closely as
				// possible (in a least squares sense)
				ImageView<double> U(divG.cols(), divG.rows());
				fill(U,0.0);
				
				//      SuccessiveOverrelaxation( U, divG, 1000);
				for (int n = 0; n < 1000; n++) {
					if (n%20 == 0) { std::cout << "\tRelax: " << n << "\n"; }
					relax_gauss_seidel(U,divG);
				}
				//      MultiGrid( U, divG);
				return U;
			}

			
			void save_gradient(std::string prefix) {
				std::ostringstream filename_prefix;
				filename_prefix << "gradient-post";
				
				ImageView<double> output = normalize(abs(_x_gradient));
				write_image(prefix + "-x.jpg", output);

				output = normalize(abs(_y_gradient));
				write_image(prefix + "-y.jpg", output);
			}
			
			void save_pyramid() {
				std::cout << "Saving Gaussian Gradient Pyramid.\n";
				
				typedef std::vector<std::pair<ImageView<double>, ImageView<double> > >::iterator iterator_type;
				int level = 0;
				for (iterator_type iter = _levels.begin(); iter != _levels.end(); ++iter) {
					ImageView<double> x_gradient = normalize(abs((*iter).first));
					ImageView<double> y_gradient = normalize(abs((*iter).second));
					
					std::ostringstream filename_prefix;
					filename_prefix << "pyramid-level-" << level++ << "-";
					write_image(filename_prefix.str() + "x.jpg", x_gradient);
					write_image(filename_prefix.str() + "y.jpg", y_gradient);
				}
			}
		};
		
		namespace tonemap {
			class GadientDomain {
				
				template <class ChannelT>
				void operator() (ImageViewBase<PixelRGB<ChannelT> > &image) {
					
					ImageView<PixelHSV<double> > hsv_image = image;
					ImageView<double> Lin = select_channel(hsv_image, 2);
					
					LuminanceGradientPyramid pyramid(Lin);
					//      pyramid.save_pyramid();
					//      pyramid.save_gradient("gradient-pre");
					pyramid.attenuate_gradient();
					//      pyramid.save_gradient("gradient-post");
					ImageView<double> Lout = normalize(pyramid.reconstruct_image());
				image = pow(image/Lin, 0.9)*Lout;
				}
			}; 
		}// namespace tonemap
		

		//! Tone map an image using tone mapping method MethodT.
		template <class ImageT, class MethodT>
		ImageT tone_map(ImageViewBase<ImageT> &image, MethodT &method) {
			return method(image.impl());
		}


}} // namespace vw::HDR
