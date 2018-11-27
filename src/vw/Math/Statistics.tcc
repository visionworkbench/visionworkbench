// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file Statistics.tcc


/// Compute the mean of an std container.
/// - It us up to the user to verify values.size() > 0
template <typename T>
double mean(T const& values) {

  double result = 0;
  typename T::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter)
    result += static_cast<double>(*iter);

  return result / static_cast<double>(values.size());
}

/// Compute the standard deviation of an std container given the mean.
template <typename T>
double standard_deviation(T const& values, double mean_value) {

  double result = 0, diff = 0;
  typename T::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter) {
    diff = static_cast<double>(*iter) - mean_value;
    result += diff * diff;
  }

  // No error checking done here as 
  return sqrt(result / static_cast<double>(values.size()-1));
}

/// Find the median of an std container
template <typename T>
double median(T const& values) {
  
  // Copy the input list to a vector
  const size_t num_elements = values.size();
  std::vector<double> v(num_elements);
  typename T::const_iterator iter = values.begin();
  for (size_t i=0; i<num_elements; ++i) {
    v[i] = *iter;
    ++iter;
  }

  // A copy is needed so we can preserve the order of the input list.  
  std::sort(v.begin(), v.end());
  size_t center = num_elements/2;
  if (num_elements % 2 == 1) // Odd length case
    return v[center];
  else  // Even length case
    return (v[center] + v[center-1]) / 2.0;
}



//---------------------------------------------------------------------------------

template <class ContainerT>
ContainerT MeanFunctor::operator() (std::vector<ContainerT> const& points) const {
  ContainerT result = points[0]; // to resize container if necessary
  size_t num_points = points.size();
  size_t dimensions = points[0].size();
  size_t last       = points[0].size() - 1;
  if (m_homogeneous)
    dimensions--;

  for (size_t i = 0; i < dimensions; ++i)
    result[i] = 0;
  if (m_homogeneous)
    result[last] = 1;

  if (m_homogeneous) {
    for (size_t i = 0; i < num_points; ++i)
      for (size_t j = 0; j < dimensions; ++j)
        result[j] += points[i][j] / points[i][last];
  }
  else {
    for (size_t i = 0; i < num_points; ++i)
      for (size_t j = 0; j < dimensions; ++j)
        result[j] += points[i][j];
  }

  for (size_t i = 0; i < dimensions; ++i)
    result[i] /= num_points;

  return result;
}


template <class ContainerT>
ContainerT StandardDeviationFunctor::operator() (std::vector<ContainerT> const& points) const {
  ContainerT result = points[0]; // to resize container if necessary
  ContainerT temp = points[0]; // to resize container if necessary
  MeanFunctor mean_func(m_homogeneous);
  ContainerT mean = mean_func(points);
  unsigned num_points = points.size();
  unsigned dimensions = points[0].size();
  unsigned last = points[0].size() - 1;
  if (m_homogeneous)
    dimensions--;

  for (unsigned int i = 0; i < dimensions; ++i)
    result[i] = 0;
  if (m_homogeneous)
    result[last] = 1;

  if (m_homogeneous) {
    for (unsigned i = 0; i < num_points; ++i)
      for (unsigned int j = 0; j < dimensions; ++j) {
        temp[j] = points[i][j] / points[i][last] - mean[j];
        result[j] += temp[j] * temp[j];
      }
  }
  else {
    for (unsigned i = 0; i < num_points; ++i)
      for (unsigned int j = 0; j < dimensions; ++j) {
        temp[j] = points[i][j] - mean[j];
        result[j] += temp[j] * temp[j];
      }
  }

  for (unsigned int i = 0; i < dimensions; ++i) {
    result[i] /= num_points;
    result[i] = sqrt(result[i]);
  }

  return result;
}

template<class T>
void find_outlier_brackets(std::vector<T> const& p,
                           double pct, double outlier_factor,
                           double & b, double & e){
  if (pct > 0.5)
    pct = 1.0 - pct;
  if (pct < 0)
    pct = 0;
    
  b = 0.0; e = 0.0; // initialize
  std::vector<T> q = p;
  std::sort(q.begin(), q.end());

  int len = q.size();
  if (len <= 0) return;

  b = q[0]; e = q[len-1];
  if (len <= 3) return; // too few points for analysis
    
  int bn = int(round(pct*len));
  int en = int(round((1.0-pct)*len))-1;

  // Adjust the bounds just in case
  bn = std::min(len-1, std::max(bn, 0));
  en = std::min(len-1, std::max(en, 0));
  
  double Q1 = q[bn];
  double Q3 = q[en];
  b = Q1 - outlier_factor*(Q3-Q1);
  e = Q3 + outlier_factor*(Q3-Q1);

  return;
}


// ------------------------------------------------------------------
//   CDFAccumulator

template <class ValT>
template <class InputIterator1, class InputIterator2>
double CDFAccumulator<ValT>::trapezoidal_rule( InputIterator1 first1, InputIterator1 last1,
                                               InputIterator2 first2 ) {
  // We will never acess last1 .. remember that
  double sum = 0.0;
  while ( first1 + 1 != last1 ) {
    sum += 0.5 * ( *(first1+1) - *first1 ) * ( *first2 + *(first2+1) );
    first1++;
    first2++;
  }
  return sum;
}

template <class ValT>
void CDFAccumulator<ValT>::pdf_differentiation( std::vector<double> const& quantile,
                                                std::vector<double> const& cdf,
                                                std::vector<double> & output_pdf ) {
  output_pdf[0] = 0;
  output_pdf[output_pdf.size()-1] = 0;

  for ( size_t i = 1; i < quantile.size() - 1; i++ ) {
    // Trying to get a centered derivative
    output_pdf[i] = 0;
    double quantile_diff = quantile[i] - quantile[i-1];
    if ( fabs(quantile_diff) > 1e-3 )
      output_pdf[i] = (cdf[i] - cdf[i-1])/quantile_diff;
    quantile_diff = quantile[i+1]-quantile[i];
    if ( fabs(quantile_diff) > 1e-3 )
      output_pdf[i] += (cdf[i+1]-cdf[i])/quantile_diff;
    output_pdf[i] /= 2;
  }
}

template <class ValT>
void CDFAccumulator<ValT>::correct_pdf_integral_to_1( std::vector<double> const& quantile,
                                                      std::vector<double> & pdf,
                                                      double additional_scalar ) {

  double pdf_error = trapezoidal_rule( quantile.begin(), quantile.end(), pdf.begin() );

  const double MAX_PDF_ERROR = 0.01;
  if ( pdf_error >= 1.0 - MAX_PDF_ERROR ) {
    double correction = (1.0 / pdf_error) * additional_scalar;
    for ( size_t i = 0; i < pdf.size(); i++ )
      pdf[i] *= correction;
  } else {
    vw_throw( MathErr() << "CDFMathAccumulator: pdf_error < 0.99: " << pdf_error );
  }
}


template <class ValT>
void CDFAccumulator<ValT>::resize( size_t buffersize, size_t quantiles ) {
  VW_ASSERT(quantiles > 0, LogicErr() << "Cannot have 0 quantiles");
  m_buffer_idx = m_num_samples = 0;
  m_sample_buf.resize( buffersize );

  m_q0 =  std::numeric_limits<double>::max();
  m_qm = -std::numeric_limits<double>::max();

  m_num_quantiles = quantiles;
  if ( !(quantiles%2) )
    m_num_quantiles++;
  m_quantile.resize(m_num_quantiles);
  m_cdf.resize(m_num_quantiles);

  // Setting a generic cdf to start things off, where 80% of the
  // distribution is in the middle third.
  size_t third = m_num_quantiles/3;
  size_t third2 = third*2;
  double slope = 10.0 / double(third);
  double first_tertile_gain = 1.0 - slope;

  // Filling middle
  for ( size_t j = third; j <= third2; j++ )
    m_cdf[j] = 0.8*(double(j-third)/double(third2-third))+0.1;
  // Filling first tertile
  for ( ssize_t j = third-1; j >= 0; j-- )
    m_cdf[j] = first_tertile_gain*m_cdf[j+1];
  // Filling third tertile
  for ( size_t j = third2+1; j < m_num_quantiles; j++ )
    m_cdf[j] = 1.0 - first_tertile_gain*(1.0-m_cdf[j-1]);
}

template <class ValT>
void CDFAccumulator<ValT>::update() {
  // Early exit if an update already happened
  if (!m_buffer_idx)
    return;

  size_t jd=0, jq=1;
  double target, told=0, tnew=0, qold, qnew;
  std::vector<double> m_new_quantile(m_num_quantiles);
  std::sort( m_sample_buf.begin(),
             m_sample_buf.begin()+m_buffer_idx ); // For partial updates
  // Setting to global min and max;
  qold = qnew = m_quantile[0] = m_new_quantile[0] = m_q0;
  m_quantile.back() = m_new_quantile.back() = m_qm;
  // .. then setting comparable probabilities
  m_cdf[0] = std::min(0.5/(m_buffer_idx+m_num_samples),
                      0.5*m_cdf[1]);
  m_cdf.back() = std::max(1-0.5/(m_buffer_idx+m_num_samples),
                          0.5*(1+m_cdf[m_num_quantiles-2]));
  // Looping over target probability values for interpolation
  for ( size_t iq = 1; iq < m_num_quantiles-1; iq++ ) {
    target = (m_num_samples+m_buffer_idx)*m_cdf[iq];
    if ( tnew < target ) {
      while (1) {
        // Locating a succession of abscissa-ordinate pairs
        // (qnew,tnew) that are the discontinuities of value or
        // slope, breaking to perform an interpolation as we cross
        // each target.
        if ( jq < m_num_quantiles &&
             ( jd >= m_buffer_idx ||
               m_quantile[jq] < m_sample_buf[jd] ) ) {
          // Found slope discontinuity from old CDF.
          qnew = m_quantile[jq];
          tnew = jd + m_num_samples*m_cdf[jq++];
          if ( tnew >= target ) break;
        } else {
          // Found value discontinuity from batch data CDF.
          qnew = m_sample_buf[jd];
          tnew = told;
          if ( m_quantile[jq] > m_quantile[jq-1] )
            tnew += m_num_samples*(m_cdf[jq]-m_cdf[jq-1])*
              (qnew-qold)/(m_quantile[jq]-m_quantile[jq-1]);
          jd++;
          if ( tnew >= target ) break;
          told = tnew++;
          qold = qnew;
          if ( tnew >= target ) break;
        }
        told = tnew;
        qold = qnew;
      }
    }
    // Performing new interpolation
    if ( tnew == told )
      m_new_quantile[iq] = 0.5*(qold+qnew);
    else
      m_new_quantile[iq] = qold + (qnew-qold)*(target-told)/(tnew-told);
    told = tnew;
    qold = qnew;
  }
  // Reset'n
  m_quantile = m_new_quantile;
  m_num_samples += m_buffer_idx;
  m_buffer_idx = 0;
}

template <class ValT>
void CDFAccumulator<ValT>::operator()( ValT const& arg ) {
  // Assimilate, We are the Borg, your data is my data!
  m_sample_buf[m_buffer_idx++] = arg;
  if ( arg < m_q0 )
    m_q0 = arg; // stretch cdf?
  if ( arg > m_qm )
    m_qm = arg;
  if ( m_buffer_idx == m_sample_buf.size() )
    update(); // merge cdf?
}

template <class ValT>
void CDFAccumulator<ValT>::operator()( CDFAccumulator<ValT>& other ) {
  
  // If this is an empty object we need duplicate the other CDF.
  // -> This not be the ideal behavior but the code below will throw an exception.
  if (m_num_samples == 0) {
    duplicate(other);
    return;
  }
  
  update();
  other.update();

  // Work out the new range of m_q0 and m_qm
  if ( m_qm < other.m_qm )
    m_qm = other.m_qm;
  if ( m_q0 > other.m_q0 )
    m_q0 = other.m_q0;

  // Generate PDF of both distribution functions by differentiating
  std::vector<double> pdf1(m_num_quantiles), pdf2(other.m_num_quantiles);
  pdf_differentiation( m_quantile, m_cdf, pdf1 );
  pdf_differentiation( other.m_quantile, other.m_cdf, pdf2 );

  // Resample both PDFs to higher number of quantiles that
  // intersects both ranges.

  // Defining quantile range
  std::vector<double> new_quantile( m_quantile.size() );
  new_quantile[0] = m_q0;
  new_quantile.back() = m_qm;
  double distance = m_qm - m_q0;
  std::transform( m_cdf.begin() + 1, m_cdf.end() - 1, new_quantile.begin() + 1,
                  std::bind1st( std::multiplies<double>(), distance ) );
  std::transform( new_quantile.begin() + 1, new_quantile.end() - 1,
                  new_quantile.begin() + 1,
                  std::bind1st( std::plus<double>(), m_q0 ) );

  // Resampling PDFs to our new quantile range so that we may
  // later add these 2 PDFs together.
  std::vector<double> pdf1_resample( m_quantile.size(), 0.0 ), pdf2_resample( m_quantile.size(), 0.0 );
  {
    size_t q1_i0 = 0, q2_i0 = 0;
    for ( size_t i = 0; i < m_quantile.size(); i++ ) {

      // Finding indexes that are just past our target quantile
      while ( q1_i0 < m_quantile.size() &&
              m_quantile[q1_i0] < new_quantile[i] ) {
        q1_i0++;
      }
      while ( q2_i0 < other.m_quantile.size() &&
              other.m_quantile[q2_i0] < new_quantile[i] ) {
        q2_i0++;
      }

      if ( q1_i0 != 0 && q1_i0 < m_quantile.size() ) {
        // The index represents represents
        // the next greatest matching quantile.
        pdf1_resample[i] = pdf1[q1_i0-1] + 
                           (pdf1[q1_i0] - pdf1[q1_i0-1]) * ( new_quantile[i] - m_quantile[q1_i0-1] ) / (m_quantile[q1_i0] - m_quantile[q1_i0-1]);
      }
      if ( q2_i0 != 0 && q2_i0 < other.m_quantile.size() ) {
        // The index represents represents
        // the next greatest matching quantile.
        pdf2_resample[i] = pdf2[q2_i0-1] + 
                           (pdf2[q2_i0] - pdf2[q2_i0-1]) * ( new_quantile[i] - other.m_quantile[q2_i0-1] ) / (other.m_quantile[q2_i0] - other.m_quantile[q2_i0-1]);
      }
    }
  }

  // Correct the problem where the resampled PDF doesn't actually sum to 1
  correct_pdf_integral_to_1( new_quantile, pdf1_resample,
                             double(m_num_samples) / double(m_num_samples+other.m_num_samples) );
  correct_pdf_integral_to_1( new_quantile, pdf2_resample,
                             double(other.m_num_samples) / double(m_num_samples+other.m_num_samples) );

  std::transform( pdf1_resample.begin(), pdf1_resample.end(),
                  pdf2_resample.begin(), pdf1_resample.begin(),
                  std::plus<double>() );


  // Integrate out a CDF
  std::vector<double> new_cdf( new_quantile.size() );
  new_cdf[0] = 0.0;
  for ( size_t i = 0; i < new_quantile.size() - 1; i++ ) {
    new_cdf[i+1] = new_cdf[i] + trapezoidal_rule( new_quantile.begin() + i,
                                                  new_quantile.begin() + i + 2,
                                                  pdf1_resample.begin() + i );
  }

  // Resample CDF to match our current rig
  m_quantile[0] = m_q0;
  m_quantile.back() = m_qm;
  {
    size_t cdf_index = 0;
    for ( size_t i = 1; i < m_quantile.size() - 1; i++ ) {
      while ( cdf_index < new_cdf.size() && new_cdf[cdf_index] < m_cdf[i] ) {
        cdf_index++;
      }
      if ( cdf_index == 0 ) {
        m_quantile[i] = new_quantile[cdf_index];
      } else if ( cdf_index == new_cdf.size() ) {
        m_quantile[i] = new_quantile.back();
      } else {
        m_quantile[i] = new_quantile[cdf_index-1] + 
                        (m_cdf[i] - new_cdf[cdf_index-1]) * (new_quantile[cdf_index]-new_quantile[cdf_index-1]) / ( new_cdf[cdf_index] - new_cdf[cdf_index-1]);
      }
    }
  }
  m_num_samples += other.m_num_samples;
}

template <class ValT>
void CDFAccumulator<ValT>::duplicate(CDFAccumulator<ValT> const& other) {
  
  m_num_quantiles = other.m_num_quantiles;
  m_buffer_idx    = other.m_buffer_idx;
  m_num_samples   = other.m_num_samples;

  m_cdf.resize(other.m_cdf.size());
  for (size_t i=0; i<m_cdf.size(); ++i)
    m_cdf[i] = other.m_cdf[i];
  m_sample_buf.resize(other.m_sample_buf.size());
  for (size_t i=0; i<m_sample_buf.size(); ++i)
    m_sample_buf[i] = other.m_sample_buf[i];
  m_quantile.resize(other.m_quantile.size());
  for (size_t i=0; i<m_quantile.size(); ++i)
    m_quantile[i] = other.m_quantile[i];

  m_q0 = other.m_q0;
  m_qm = other.m_qm;
}

template <class ValT>
ValT CDFAccumulator<ValT>::quantile( double const& arg ) const {
  double q;

  // if ( m_buffer_idx > 0 ) update();
  size_t jl=0,jh=m_num_quantiles-1,j;
  while ( jh - jl > 1 ) {
    j = (jh+jl)>>1;
    if ( arg > m_cdf[j] ) jl=j;
    else jh=j;
  }
  j = jl;
  q = m_quantile[j]+(m_quantile[j+1]-m_quantile[j])*(arg-m_cdf[j])/(m_cdf[j+1]-m_cdf[j]);

  // Keeping estimate in CDF
  return std::max(m_quantile[0],std::min(m_quantile.back(),q));
}

template <class ValT>
ValT CDFAccumulator<ValT>::approximate_mean( float const& stepping ) const {
  ValT   mean = 0;
  size_t count = 0;
  for ( float i = stepping; i < 1+stepping; i+=stepping ) {
    count++;
    mean += ( quantile(i) + quantile(i-stepping) ) / 2.0;
  }
  return mean / ValT(count);
}

template <class ValT>
ValT CDFAccumulator<ValT>::approximate_stddev( float const& stepping ) const {
  ValT mean = approximate_mean(stepping);
  ValT stddev = 0;
  size_t count = 0;
  for ( float i = stepping; i < 1+stepping; i+=stepping ) {
    count++;
    stddev += pow((quantile(i) + quantile(i-stepping))/2-mean,2);
  }
  return sqrt(stddev/ValT(count));
}


//--------------------------------------------------------------------------
// Class Histogram
inline
void Histogram::initialize(size_t num_bins, double min_value, double max_value) {

  if (num_bins == 0)
    vw_throw(ArgumentErr() << "Can't create a Histogram with zero bins!");

  m_num_bins    = num_bins;
  m_max_bin     = num_bins-1;
  m_num_values  = 0; 
  m_min_value   = min_value;
  m_max_value   = max_value;
  m_bin_centers.resize(num_bins);
  m_bin_values.resize (num_bins);

  
  // Compute the bin centers
  m_range     = max_value - min_value;
  m_bin_width = m_range / num_bins;
  
  for (size_t i=0; i<num_bins; ++i) {
    m_bin_centers[i] = min_value + i*m_bin_width + m_bin_width/2.0;
    m_bin_values [i] = 0;
  }

}
inline
void Histogram::add_value(double value, bool saturate) {
  // Compute the bin
  int bin = static_cast<int>(round( m_max_bin * ( (value - m_min_value)/m_range ) ));
  // Handle out of range input values
  if (bin < 0) {
    if (saturate)
      bin = 0;
    else
      return;
  }
  if (bin >= m_num_bins) {
    if (saturate)
      bin = m_max_bin;
    else
      return;
  }
  m_bin_values[bin] += 1.0;
  
  ++m_num_values;
}
inline
void Histogram::add_value_no_check(double value) {
  // Compute the bin
  int bin = static_cast<int>(round( m_max_bin * ( (value - m_min_value)/m_range ) ));
  m_bin_values[bin] += 1.0;    
  ++m_num_values;
}
inline
size_t Histogram::get_percentile(double percentile) const {

  // Verify the input percentile is in the legal range
  if ((percentile < 0) || (percentile > 1.0)) {
    vw_throw(ArgumentErr() << "get_histogram_percentile: illegal percentile request: " << percentile << "\n");
  }

  double sum = static_cast<double>(m_num_values);

  double running_percentile = 0;
  for (int i=0; i<m_num_bins; ++i) {
    double this_percent = get_bin_value(i) / sum;
    running_percentile += this_percent;
    if (running_percentile >= percentile)
      return i;
  }
  vw_throw(LogicErr() << "get_histogram_percentile: Illegal histogram encountered!");
}
inline
void Histogram::write_to_disk(std::string const& path) const {
  std::ofstream f(path.c_str());
  f << "# Histogram_Value, Histogram_Count";
  for (int i=0; i<m_num_bins; ++i) {
    f << get_bin_center(i) << ", " << get_bin_value(i) << std::endl;
  }
  f.close();
}

