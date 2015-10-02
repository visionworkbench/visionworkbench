package Statistics::Descriptive;

use strict;
use warnings;

##This module draws heavily from perltoot v0.4 from Tom Christiansen.

require 5.00404;  ##Yes, this is underhanded, but makes support for me easier
		  ##Not only that, but it's the latest "safe" version of
		  ##Perl5.  01-03 weren't bug free.
use vars (qw($VERSION $Tolerance));

$VERSION = '3.0100';

$Tolerance = 0.0;

package Statistics::Descriptive::Sparse;

use vars qw(%fields);
use Carp;

sub _make_accessors
{
    my ($pkg, $methods) = @_;

    no strict 'refs';
    foreach my $method (@$methods)
    {
        *{$pkg."::".$method} =
            do {
                my $m = $method;
                sub {
                    my $self = shift;

                    if (@_)
                    {
                        $self->{$m} = shift;
                    }
                    return $self->{$m};
                };
            };
    }

    return;
}

sub _make_private_accessors
{
    my ($pkg, $methods) = @_;

    no strict 'refs';
    foreach my $method (@$methods)
    {
        *{$pkg."::_".$method} =
            do {
                my $m = $method;
                sub {
                    my $self = shift;

                    if (@_)
                    {
                        $self->{$m} = shift;
                    }
                    return $self->{$m};
                };
            };
    }

    return;
}

##Define the fields to be used as methods
%fields = (
  count			=> 0,
  mean			=> 0,
  sum			=> 0,
  sumsq			=> 0,
  min			=> undef,
  max			=> undef,
  mindex		=> undef,
  maxdex		=> undef,
  sample_range		=> undef,
  variance => undef,
  );

__PACKAGE__->_make_accessors( [ grep { $_ ne "variance" } keys(%fields) ] );
__PACKAGE__->_make_accessors( ["_permitted"] );
__PACKAGE__->_make_private_accessors(["variance"]);

sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  my $self = {
    %fields,
  };
  bless ($self, $class);
  $self->_permitted(\%fields);
  return $self;
}

sub _is_permitted
{
    my $self = shift;
    my $key = shift;

    return exists($self->_permitted()->{$key});
}

sub add_data {
  my $self = shift;  ##Myself
  my $oldmean;
  my ($min,$mindex,$max,$maxdex,$sum,$sumsq,$count);
  my $aref;

  if (ref $_[0] eq 'ARRAY') {
    $aref = $_[0];
  }
  else {
    $aref = \@_;
  }

  ##If we were given no data, we do nothing.
  return 1 if (!@{ $aref });

  ##Take care of appending to an existing data set
  
  if (!defined($min = $self->min()))
  {
      $min = $aref->[$mindex = 0];
  }
  else
  {
      $mindex = $self->mindex();
  }

  if (!defined($max = $self->max()))
  {
      $max = $aref->[$maxdex = 0];
  }
  else
  {
      $maxdex = $self->maxdex();
  }

  $sum = $self->sum();
  $sumsq = $self->sumsq();
  $count = $self->count();

  ##Calculate new mean, sumsq, min and max;
  foreach ( @{ $aref } ) {
    $sum += $_;
    $sumsq += $_**2;
    $count++;
    if ($_ >= $max) {
      $max = $_;
      $maxdex = $count-1;
    }
    if ($_ <= $min) {
      $min = $_;
      $mindex = $count-1;
    }
  }

  $self->min($min);
  $self->mindex($mindex);
  $self->max($max);
  $self->maxdex($maxdex);
  $self->sample_range($max - $min);
  $self->sum($sum);
  $self->sumsq($sumsq);
  $self->mean($sum / $count);
  $self->count($count);
  ##indicator the value is not cached.  Variance isn't commonly enough
  ##used to recompute every single data add.
  $self->_variance(undef());
  return 1;
}

sub standard_deviation {
  my $self = shift;  ##Myself
  return undef if (!$self->count());
  return sqrt($self->variance());
}

##Return variance; if needed, compute and cache it.
sub variance {
  my $self = shift;  ##Myself
  my $div = @_ ? 0 : 1;
  my $count = $self->count();
  if ($count < 1 + $div) {
      return 0;
  }

  if (!defined($self->_variance())) {
    my $variance = ($self->sumsq()- $count * $self->mean()**2);

    # Sometimes due to rounding errors we get a number below 0.
    # This makes sure this is handled as gracefully as possible.
    #
    # See:
    #
    # https://rt.cpan.org/Public/Bug/Display.html?id=46026
    if ($variance < 0)
    {
        $variance = 0;
    }

    $variance /= $count - $div;

    $self->_variance($variance);
  }
  return $self->_variance();
}

##Clear a stat.  More efficient than destroying an object and calling
##new.
sub clear {
  my $self = shift;  ##Myself
  my $key;

  return if (!$self->count());
  while (my($field, $value) = each %fields) {
    $self->{$field} = $value;
  }
}

1;

package Statistics::Descriptive::Full;

use Carp;

use POSIX ();

use vars qw(@ISA $a $b %fields);

@ISA = qw(Statistics::Descriptive::Sparse);

##Create a list of fields not to remove when data is updated
%fields = (
  _permitted => undef,  ##Place holder for the inherited key hash
  data       => undef,  ##Our data
  presorted  => undef,  ##Flag to indicate the data is already sorted
  _reserved  => undef,  ##Place holder for this lookup hash
);

__PACKAGE__->_make_private_accessors(
    [qw(data frequency geometric_mean harmonic_mean 
        least_squares_fit median mode
       )
    ]
);
__PACKAGE__->_make_accessors([qw(presorted _reserved _trimmed_mean_cache)]);

sub _clear_fields
{
    my $self = shift;

    # Empty array ref for holding data later!
    $self->_data([]);
    $self->_reserved(\%fields);
    $self->presorted(0);
    $self->_trimmed_mean_cache(+{});

    return;
}

##Have to override the base method to add the data to the object
##The proxy method from above is still valid
sub new {
  my $proto = shift;
  my $class = ref($proto) || $proto;
  # Create my self re SUPER
  my $self = $class->SUPER::new();  
  bless ($self, $class);  #Re-anneal the object
  $self->_clear_fields();
  return $self;
}

sub _is_reserved
{
    my $self = shift;
    my $field = shift;

    return exists($self->_reserved->{$field});
}

sub _delete_all_cached_keys
{
    my $self = shift;

    KEYS_LOOP:
    foreach my $key (keys %{ $self }) { # Check each key in the object
        # If it's a reserved key for this class, keep it
        if ($self->_is_reserved($key) || $self->_is_permitted($key))
        {
            next KEYS_LOOP;
        }
        delete $self->{$key};          # Delete the out of date cached key
    }
    return;
}

##Clear a stat.  More efficient than destroying an object and calling
##new.
sub clear {
    my $self = shift;  ##Myself
    my $key;

    if (!$self->count())
    {
        return;
    }

    $self->_delete_all_cached_keys();
    $self->SUPER::clear();
    $self->_clear_fields();
}

sub add_data {
  my $self = shift;
  my $aref;

  if (ref $_[0] eq 'ARRAY') {
    $aref = $_[0];
  }
  else {
    $aref = \@_;
  }
  $self->SUPER::add_data($aref);  ##Perform base statistics on the data
  push @{ $self->_data() }, @{ $aref };
  ##Clear the presorted flag
  $self->presorted(0);

  $self->_delete_all_cached_keys();

  return 1;
}

sub get_data {
  my $self = shift;
  return @{ $self->_data() };
}

sub sort_data {
  my $self = shift;

  if (! $self->presorted())
  {
      ##Sort the data in descending order
      $self->_data([ sort {$a <=> $b} @{$self->_data()} ]);
      $self->presorted(1);
      ##Fix the maxima and minima indices
      $self->mindex(0);
      $self->maxdex($#{$self->_data()});
  }

  return 1;
}

sub percentile {
  my $self = shift;
  my $percentile = shift || 0;
  ##Since we're returning a single value there's no real need
  ##to cache this.

  ##If the requested percentile is less than the "percentile bin
  ##size" then return undef.  Check description of RFC 2330 in the
  ##POD below.
  my $count = $self->count();
  return undef if $percentile < 100 / $count;

  $self->sort_data();
  my $num = $count*$percentile/100;
  my $index = &POSIX::ceil($num) - 1;
  my $val = $self->_data->[$index];
  return wantarray
    ? ($val, $index)
    : $val
    ;
}

sub _calc_new_median
{
    my $self = shift;
    my $count = $self->count();

    ##Even or odd
    if ($count % 2)
    {   
        return $self->_data->[($count-1)/2];
    }
    else
    {
        return
        (
            ($self->_data->[($count)/2] + $self->_data->[($count-2)/2] ) / 2
        );
    }
}

sub median {
    my $self = shift;

    ##Cached?
    if (! defined($self->_median()))
    {
        $self->sort_data();
        $self->_median($self->_calc_new_median());
    }
    return $self->_median();
}

sub quantile {
    my ( $self, $QuantileNumber ) = @_;

    unless ( defined $QuantileNumber and $QuantileNumber =~ m/^0|1|2|3|4$/ ) {
       carp("Bad quartile type, must be 0, 1, 2, 3 or 4\n");
       return;
    }
    
    $self->sort_data();

    return $self->_data->[0] if ( $QuantileNumber == 0 );

    my $count = $self->count();

    return $self->_data->[ $count - 1 ] if ( $QuantileNumber == 4 );

    my $K_quantile = ( ( $QuantileNumber / 4 ) * ( $count - 1 ) + 1 );
    my $F_quantile = $K_quantile - POSIX::floor($K_quantile);
    $K_quantile = POSIX::floor($K_quantile);

    # interpolation
    my $aK_quantile     = $self->_data->[ $K_quantile - 1 ];
    return $aK_quantile if ( $F_quantile == 0 );
    my $aKPlus_quantile = $self->_data->[$K_quantile];
    
    # Calcul quantile
    my $quantile = $aK_quantile
      + ( $F_quantile * ( $aKPlus_quantile - $aK_quantile ) );

    return $quantile;
}

sub _real_calc_trimmed_mean
{
    my $self = shift;
    my $lower = shift;
    my $upper = shift;

    my $lower_trim = int ($self->count()*$lower); 
    my $upper_trim = int ($self->count()*$upper); 
    my ($val,$oldmean) = (0,0);
    my ($tm_count,$tm_mean,$index) = (0,0,$lower_trim);

    $self->sort_data();
    while ($index <= $self->count() - $upper_trim -1)
    {
        $val = $self->_data()->[$index];
        $oldmean = $tm_mean;
        $index++;
        $tm_count++;
        $tm_mean += ($val - $oldmean) / $tm_count;
    }

    return $tm_mean;
}

sub trimmed_mean
{
    my $self = shift;
    my ($lower,$upper);
    #upper bound is in arg list or is same as lower
    if (@_ == 1)
    {
        ($lower,$upper) = ($_[0],$_[0]);
    }
    else
    {
        ($lower,$upper) = ($_[0],$_[1]);
    }

    ##Cache
    my $thistm = join ':',$lower,$upper;
    my $cache = $self->_trimmed_mean_cache();
    if (!exists($cache->{$thistm}))
    {
        $cache->{$thistm} = $self->_real_calc_trimmed_mean($lower, $upper);
    }

    return $cache->{$thistm};
}

sub _test_for_too_small_val
{
    my $self = shift;
    my $val = shift;

    return (abs($val) <= $Statistics::Descriptive::Tolerance);
}

sub _calc_harmonic_mean
{
    my $self = shift;

    my $hs = 0;

    foreach my $item ( @{$self->_data()} )
    {
        ##Guarantee that there are no divide by zeros
        if ($self->_test_for_too_small_val($item))
        {
            return;
        }

        $hs += 1/$item;
    }

    if ($self->_test_for_too_small_val($hs))
    {
        return;
    }

    return $self->count()/$hs;
}

sub harmonic_mean
{
    my $self = shift;

    if (!defined($self->_harmonic_mean()))
    {
        $self->_harmonic_mean(scalar($self->_calc_harmonic_mean()));
    }

    return $self->_harmonic_mean();
}

sub mode
{
    my $self = shift;

    if (!defined ($self->_mode()))
    {
        my $mode = 0;
        my $occurances= 0;
        my $flag = 1;

        my %count;

        foreach my $item (@{ $self->_data() })
        {
            $count{$item}++;
            $flag = 0 if ($count{$item} > 1);
        }

        #Distribution is flat - no mode exists
        if ($flag)
        {
            return undef;
        }

        foreach my $val (keys %count)
        {
            if ($count{$val} > $occurances)
            {
                $occurances = $count{$val};
                $mode = $val;
            }
        }

        $self->_mode($mode);
    }

    return $self->_mode();
}

sub geometric_mean {
    my $self = shift;

    if (!defined($self->_geometric_mean()))
    {
        my $gm = 1;
        my $exponent = 1/$self->count();

        for my $val (@{ $self->_data() })
        {
            if ($val < 0)
            {
                return undef;
            }
            $gm *= $val**$exponent;
        }

        $self->_geometric_mean($gm);
    }

    return $self->_geometric_mean();
}

sub frequency_distribution_ref
{
    my $self = shift;
    my @k = ();
    # Must have at least two elements
    if ($self->count() < 2)
    {
        return undef;
    }

    if ((!@_) && (defined $self->_frequency()))
    {
        return $self->_frequency()
    }

    my %bins;
    my $partitions = shift;

    if (ref($partitions) eq 'ARRAY')
    {
        @k = @{ $partitions };
        return undef unless @k;  ##Empty array
        if (@k > 1) {
            ##Check for monotonicity
            my $element = $k[0];
            for my $next_elem (@k[1..$#k]) {
                if ($element > $next_elem) {
                    carp "Non monotonic array cannot be used as frequency bins!\n";
                    return undef;
                }
                $element = $next_elem;
            }
        }
        %bins = map { $_ => 0 } @k;
    }
    else
    {
        return undef unless $partitions >= 1;
        my $interval = $self->sample_range() / $partitions;
        foreach my $idx (1 .. ($partitions-1))
        {
            push @k, ($self->min() + $idx * $interval);
        }

        $bins{$self->max()} = 0;

        push @k, $self->max();
    }

    ELEMENT:
    foreach my $element (@{$self->_data()})
    {
        foreach my $limit (@k)
        {
            if ($element <= $limit)
            {
                $bins{$limit}++;
                next ELEMENT;
            }
        }
    }

    return $self->_frequency(\%bins);
}

sub frequency_distribution {
    my $self = shift;

    my $ret = $self->frequency_distribution_ref(@_);

    if (!defined($ret))
    {
        return undef;
    }
    else
    {
        return %$ret;
    }
}

sub least_squares_fit {
  my $self = shift;
  return () if $self->count() < 2;

  ##Sigma sums
  my ($sigmaxy, $sigmax, $sigmaxx, $sigmayy, $sigmay) = (0,0,0,0,$self->sum);
  my ($xvar, $yvar, $err);

  ##Work variables
  my ($iter,$y,$x,$denom) = (0,0,0,0);
  my $count = $self->count();
  my @x;

  ##Outputs
  my ($m, $q, $r, $rms);

  if (!defined $_[1]) {
    @x = 1..$self->count();
  }
  else {
    @x = @_;
    if ( $self->count() != scalar @x) {
      carp "Range and domain are of unequal length.";
      return ();
    }
  }
  foreach $x (@x) {
    $y = $self->_data->[$iter];
    $sigmayy += $y * $y;
    $sigmaxx += $x * $x;
    $sigmaxy += $x * $y;
    $sigmax  += $x;
    $iter++;
  }
  $denom = $count * $sigmaxx - $sigmax*$sigmax;
  return ()
    unless abs( $denom ) > $Statistics::Descriptive::Tolerance;

  $m = ($count*$sigmaxy - $sigmax*$sigmay) / $denom;
  $q = ($sigmaxx*$sigmay - $sigmax*$sigmaxy ) / $denom;

  $xvar = $sigmaxx - $sigmax*$sigmax / $count;
  $yvar = $sigmayy - $sigmay*$sigmay / $count;

  $denom = sqrt( $xvar * $yvar );
  return () unless (abs( $denom ) > $Statistics::Descriptive::Tolerance);
  $r = ($sigmaxy - $sigmax*$sigmay / $count )/ $denom;

  $iter = 0;
  $rms = 0.0;
  foreach (@x) {
    ##Error = Real y - calculated y
    $err = $self->_data->[$iter] - ( $m * $_ + $q );
    $rms += $err*$err;
    $iter++;
  }

  $rms = sqrt($rms / $count);
  
  $self->_least_squares_fit([$q, $m, $r, $rms]);

  return @{ $self->_least_squares_fit() };
}

1;

package Statistics::Descriptive;

##All modules return true.
1;

__END__

=head1 NAME

Statistics::Descriptive - Module of basic descriptive statistical functions.

=head1 SYNOPSIS

  use Statistics::Descriptive;
  $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(1,2,3,4); $mean = $stat->mean();
  $var  = $stat->variance();
  $tm   = $stat->trimmed_mean(.25);
  $Statistics::Descriptive::Tolerance = 1e-10;

=head1 DESCRIPTION

This module provides basic functions used in descriptive statistics.
It has an object oriented design and supports two different types of
data storage and calculation objects: sparse and full. With the sparse
method, none of the data is stored and only a few statistical measures
are available. Using the full method, the entire data set is retained
and additional functions are available.

Whenever a division by zero may occur, the denominator is checked to be
greater than the value C<$Statistics::Descriptive::Tolerance>, which
defaults to 0.0. You may want to change this value to some small
positive value such as 1e-24 in order to obtain error messages in case
of very small denominators.

Many of the methods (both Sparse and Full) cache values so that subsequent
calls with the same arguments are faster.

=head1 METHODS

=head2 Sparse Methods

=over 5

=item $stat = Statistics::Descriptive::Sparse->new();

Create a new sparse statistics object.

=item $stat->clear();

Effectively the same as

  my $class = ref($stat);
  undef $stat;
  $stat = new $class;

except more efficient.

=item $stat->add_data(1,2,3);

Adds data to the statistics variable. The cached statistical values are 
updated automatically.

=item $stat->count();

Returns the number of data items.

=item $stat->mean();

Returns the mean of the data.

=item $stat->sum();

Returns the sum of the data.

=item $stat->variance();

Returns the variance of the data.  Division by n-1 is used.

=item $stat->standard_deviation();

Returns the standard deviation of the data. Division by n-1 is used.

=item $stat->min();

Returns the minimum value of the data set.

=item $stat->mindex();

Returns the index of the minimum value of the data set.

=item $stat->max();

Returns the maximum value of the data set.

=item $stat->maxdex();

Returns the index of the maximum value of the data set.

=item $stat->sample_range();

Returns the sample range (max - min) of the data set.

=back

=head2 Full Methods

Similar to the Sparse Methods above, any Full Method that is called caches
the current result so that it doesn't have to be recalculated.  In some
cases, several values can be cached at the same time.

=over 5

=item $stat = Statistics::Descriptive::Full->new();

Create a new statistics object that inherits from
Statistics::Descriptive::Sparse so that it contains all the methods
described above.

=item $stat->add_data(1,2,4,5);

Adds data to the statistics variable.  All of the sparse statistical
values are updated and cached.  Cached values from Full methods are
deleted since they are no longer valid.  

I<Note:  Calling add_data with an empty array will delete all of your
Full method cached values!  Cached values for the sparse methods are
not changed>

=item $stat->get_data();

Returns a copy of the data array.

=item $stat->sort_data();

Sort the stored data and update the mindex and maxdex methods.  This
method uses perl's internal sort.

=item $stat->presorted(1);

=item $stat->presorted();

If called with a non-zero argument, this method sets a flag that says
the data is already sorted and need not be sorted again.  Since some of
the methods in this class require sorted data, this saves some time.
If you supply sorted data to the object, call this method to prevent
the data from being sorted again. The flag is cleared whenever add_data
is called.  Calling the method without an argument returns the value of
the flag.

=item $x = $stat->percentile(25);

=item ($x, $index) = $stat->percentile(25);

Sorts the data and returns the value that corresponds to the
percentile as defined in RFC2330:

=over 4

=item

For example, given the 6 measurements:

-2, 7, 7, 4, 18, -5

Then F(-8) = 0, F(-5) = 1/6, F(-5.0001) = 0, F(-4.999) = 1/6, F(7) =
5/6, F(18) = 1, F(239) = 1.

Note that we can recover the different measured values and how many
times each occurred from F(x) -- no information regarding the range
in values is lost.  Summarizing measurements using histograms, on the
other hand, in general loses information about the different values
observed, so the EDF is preferred.

Using either the EDF or a histogram, however, we do lose information
regarding the order in which the values were observed.  Whether this
loss is potentially significant will depend on the metric being
measured.

We will use the term "percentile" to refer to the smallest value of x
for which F(x) >= a given percentage.  So the 50th percentile of the
example above is 4, since F(4) = 3/6 = 50%; the 25th percentile is
-2, since F(-5) = 1/6 < 25%, and F(-2) = 2/6 >= 25%; the 100th
percentile is 18; and the 0th percentile is -infinity, as is the 15th
percentile.

Care must be taken when using percentiles to summarize a sample,
because they can lend an unwarranted appearance of more precision
than is really available.  Any such summary must include the sample
size N, because any percentile difference finer than 1/N is below the
resolution of the sample.

=back

(Taken from:
I<RFC2330 - Framework for IP Performance Metrics>,
Section 11.3.  Defining Statistical Distributions.
RFC2330 is available from:
L<http://www.ietf.org/rfc/rfc2330.txt> .)

If the percentile method is called in a list context then it will
also return the index of the percentile.

=item $x = $stat->quantile($Type);

Sorts the data and returns estimates of underlying distribution quantiles based on one 
or two order statistics from the supplied elements.

This method use the same algorithm as Excel and R language (quantile B<type 7>).

The generic function quantile produces sample quantiles corresponding to the given probabilities.

B<$Type> is an integer value between 0 to 4 :

  0 => zero quartile (Q0) : minimal value
  1 => first quartile (Q1) : lower quartile = lowest cut off (25%) of data = 25th percentile
  2 => second quartile (Q2) : median = it cuts data set in half = 50th percentile
  3 => third quartile (Q3) : upper quartile = highest cut off (25%) of data, or lowest 75% = 75th percentile
  4 => fourth quartile (Q4) : maximal value

Exemple : 

  my @data = (1..10);
  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@data);
  print $stat->quantile(0); # => 1
  print $stat->quantile(1); # => 3.25
  print $stat->quantile(2); # => 5.5
  print $stat->quantile(3); # => 7.75
  print $stat->quantile(4); # => 10


=item $stat->median();

Sorts the data and returns the median value of the data.

=item $stat->harmonic_mean();

Returns the harmonic mean of the data.  Since the mean is undefined
if any of the data are zero or if the sum of the reciprocals is zero,
it will return undef for both of those cases.

=item $stat->geometric_mean();

Returns the geometric mean of the data.

=item $stat->mode();

Returns the mode of the data. 

=item $stat->trimmed_mean(ltrim[,utrim]);

C<trimmed_mean(ltrim)> returns the mean with a fraction C<ltrim> 
of entries at each end dropped. C<trimmed_mean(ltrim,utrim)> 
returns the mean after a fraction C<ltrim> has been removed from the
lower end of the data and a fraction C<utrim> has been removed from the
upper end of the data.  This method sorts the data before beginning
to analyze it.

All calls to trimmed_mean() are cached so that they don't have to be
calculated a second time.

=item $stat->frequency_distribution_ref($partitions);

=item $stat->frequency_distribution_ref(\@bins);

=item $stat->frequency_distribution_ref();

C<frequency_distribution_ref($partitions)> slices the data into
C<$partition> sets (where $partition is greater than 1) and counts the
number of items that fall into each partition. It returns a reference to
a hash where the keys are the numerical values of the
partitions used. The minimum value of the data set is not a key and the
maximum value of the data set is always a key. The number of entries
for a particular partition key are the number of items which are
greater than the previous partition key and less then or equal to the
current partition key. As an example,

   $stat->add_data(1,1.5,2,2.5,3,3.5,4);
   $f = $stat->frequency_distribution_ref(2);
   for (sort {$a <=> $b} keys %$f) {
      print "key = $_, count = $f->{$_}\n";
   }

prints

   key = 2.5, count = 4
   key = 4, count = 3

since there are four items less than or equal to 2.5, and 3 items
greater than 2.5 and less than 4.

C<frequency_distribution_refs(\@bins)> provides the bins that are to be used
for the distribution.  This allows for non-uniform distributions as
well as trimmed or sample distributions to be found.  C<@bins> must
be monotonic and contain at least one element.  Note that unless the
set of bins contains the range that the total counts returned will
be less than the sample size.

Calling C<frequency_distribution_ref()> with no arguments returns the last
distribution calculated, if such exists.

=item my %hash = $stat->frequency_distribution($partitions);

=item my %hash = $stat->frequency_distribution(\@bins);

=item my %hash = $stat->frequency_distribution();

Same as C<frequency_distribution_ref()> except that returns the hash clobbered
into the return list. Kept for compatibility reasons with previous
versions of Statistics::Descriptive and using it is discouraged.

=item $stat->least_squares_fit();

=item $stat->least_squares_fit(@x);

C<least_squares_fit()> performs a least squares fit on the data,
assuming a domain of C<@x> or a default of 1..$stat->count().  It
returns an array of four elements C<($q, $m, $r, $rms)> where

=over 4

=item C<$q and $m>

satisfy the equation C($y = $m*$x + $q).

=item C<$r>

is the Pearson linear correlation cofficient.

=item C<$rms>

is the root-mean-square error.

=back

If case of error or division by zero, the empty list is returned.

The array that is returned can be "coerced" into a hash structure
by doing the following:

  my %hash = ();
  @hash{'q', 'm', 'r', 'err'} = $stat->least_squares_fit();

Because calling C<least_squares_fit()> with no arguments defaults
to using the current range, there is no caching of the results.

=back

=head1 REPORTING ERRORS

I read my email frequently, but since adopting this module I've added 2
children and 1 dog to my family, so please be patient about my response
times.  When reporting errors, please include the following to help
me out:

=over 4

=item *

Your version of perl.  This can be obtained by typing perl C<-v> at
the command line.

=item *

Which version of Statistics::Descriptive you're using.  As you can
see below, I do make mistakes.  Unfortunately for me, right now
there are thousands of CD's with the version of this module with
the bugs in it.  Fortunately for you, I'm a very patient module
maintainer.

=item *

Details about what the error is.  Try to narrow down the scope
of the problem and send me code that I can run to verify and
track it down.

=back

=head1 AUTHOR

Current maintainer:

Shlomi Fish, L<http://www.shlomifish.org/> , C<shlomif@cpan.org>

Previously:

Colin Kuskie

My email address can be found at http://www.perl.com under Who's Who
or at: http://search.cpan.org/author/COLINK/.

=head1 REFERENCES

RFC2330, Framework for IP Performance Metrics

The Art of Computer Programming, Volume 2, Donald Knuth.

Handbook of Mathematica Functions, Milton Abramowitz and Irene Stegun.

Probability and Statistics for Engineering and the Sciences, Jay Devore.

=head1 COPYRIGHT

Copyright (c) 1997,1998 Colin Kuskie. All rights reserved.  This
program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

Copyright (c) 1998 Andrea Spinelli. All rights reserved.  This program
is free software; you can redistribute it and/or modify it under the
same terms as Perl itself.

Copyright (c) 1994,1995 Jason Kastner. All rights
reserved.  This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 LICENSE

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=head1 REVISION HISTORY

=over 4

=item v2.3

Rolled into November 1998

Code provided by Andrea Spinelli to prevent division by zero and to
make consistent return values for undefined behavior.  Andrea also
provided a test bench for the module.

A bug fix for the calculation of frequency distributions.  Thanks to Nick
Tolli for alerting this to me.

Added 4 lines of code to Makefile.PL to make it easier for the ActiveState
installation tool to use.  Changes work fine in perl5.004_04, haven't
tested them under perl5.005xx yet.

=item v2.2

Rolled into March 1998.

Fixed problem with sending 0's and -1's as data.  The old 0 : true ? false
thing.  Use defined to fix.

Provided a fix for AUTOLOAD/DESTROY/Carp bug.  Very strange.

=item v2.1

August 1997

Fixed errors in statistics algorithms caused by changing the
interface.

=item v2.0

August 1997

Fixed errors in removing cached values (they weren't being removed!)
and added sort_data and presorted methods.

June 1997

Transferred ownership of the module from Jason to Colin.

Rewrote OO interface, modified function distribution, added mindex,
maxdex.

=item v1.1

April 1995

Added LeastSquaresFit and FrequencyDistribution.

=item v1.0 

March 1995

Released to comp.lang.perl and placed on archive sites.

=item v.20

December 1994

Complete rewrite after extensive and invaluable e-mail 
correspondence with Anno Siegel.

=item v.10

December 1994

Initital concept, released to perl5-porters list.

=back

=cut
