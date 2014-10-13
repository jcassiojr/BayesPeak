//A number of functions from the GNU Scientific Library
//=====================================================


#include <math.h>
#include <limits.h>
#define GSL_DBL_EPSILON        2.2204460492503131e-16

//RNG
//---

#define N_SHUFFLE 32
#define N_DIV (1 + 2147483562/N_SHUFFLE)

typedef struct
  {
    unsigned long int x;
    unsigned long int y;
    unsigned long int n;
    unsigned long int shuffle[N_SHUFFLE];
  }
ran2_state_t;

void ran2_set (const ran2_state_t* vstate, unsigned long int s)
{
	static const long int m1 = 2147483563, a1 = 40014, q1 = 53668, r1 = 12211;
	//static const long int m2 = 2147483399, a2 = 40692, q2 = 52774, r2 = 3791;
	//JMC: not sure what this second line is for? alternative (m, a, q, r) values?

  ran2_state_t *state = (ran2_state_t *) vstate;
  int i;

  if (s == 0)
    s = 1;      /* default seed is 1 */

  state->y = s;

  for (i = 0; i < 8; i++)
    {
      long int h = s / q1;
      long int t = a1 * (s - h * q1) - h * r1;
      if (t < 0)
        t += m1;
      s = t;
    }

  for (i = N_SHUFFLE - 1; i >= 0; i--)
    {
      long int h = s / q1;
      long int t = a1 * (s - h * q1) - h * r1;
      if (t < 0)
        t += m1;
      s = t;
      state->shuffle[i] = s;
    }

  state->x = s;
  state->n = s;

  return;
}

unsigned long int ran2_get (const ran2_state_t* vstate)
{
	static const long int m1 = 2147483563, a1 = 40014, q1 = 53668, r1 = 12211;
	static const long int m2 = 2147483399, a2 = 40692, q2 = 52774, r2 = 3791;


  ran2_state_t *state = (ran2_state_t *) vstate;

  const unsigned long int x = state->x;
  const unsigned long int y = state->y;

  long int h1 = x / q1;
  long int t1 = a1 * (x - h1 * q1) - h1 * r1;

  long int h2 = y / q2;
  long int t2 = a2 * (y - h2 * q2) - h2 * r2;

  if (t1 < 0)
    t1 += m1;

  if (t2 < 0)
    t2 += m2;

  state->x = t1;
  state->y = t2;

  {
    unsigned long int j = state->n / N_DIV;
    long int delta = state->shuffle[j] - t2;
    if (delta < 1)
      delta += m1 - 1;
    state->n = delta;
    state->shuffle[j] = t1;
  }

  return state->n;
}


double ran2_get_double (const ran2_state_t* vstate)
{
  float x_max = 1 - 1.2e-7f ; /* Numerical Recipes version of 1-FLT_EPS */

  float x = ran2_get (vstate) / 2147483563.0f ;
 
  if (x > x_max) 
    return x_max ;
  
  return x ;
}

unsigned long int
gsl_rng_uniform_int (const ran2_state_t* vstate, unsigned long int n)
{
  unsigned long int offset = 0;
  unsigned long int range = 2147483563;
  unsigned long int scale;
  unsigned long int k;

  if (n > range || n == 0) 
    {
      //fail miserably
    }

  scale = range / n;

  do
    {
      k = ((ran2_get(vstate)) - offset) / scale;
    }
  while (k >= n);

  return k;
}


//LOG GAMMA FIXME
//---------

/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
 -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
 -176.61502916214059906584551354,
  12.507343278686904814458936853,
 -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};


/* Lanczos method for real x > 0;
 * gamma=7, truncated at 1/(z+8) 
 * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
 */

#define LogRootTwoPi_  0.9189385332046727418

static double lngammafn(double x)
{
  int k;
  double Ag;
  double term1, term2;

  x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag = lanczos_7_c[0];
  for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

  /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
  term1 = (x+0.5)*log((x+7.5)/M_E);
  term2 = LogRootTwoPi_ + log(Ag);
  return (term1 + (term2 - 7.0));
  //result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
  //result->err += GSL_DBL_EPSILON * fabs(result->val);
}


//LOG FACTORIAL
//-------------

double logfact(int x)
{
	static double calc[170];
	if(x <= 1) return 0;
	if(x < 170)
	{
		if(calc[x] > 0)
		{
			return calc[x];
		}
		else
		{
			return calc[x] = lngammafn(x + 1);
		}
	}
	return lngammafn(x+1);
}


//NORMAL DEVIATES 
//---------------


double
gsl_ran_gaussian_ratio_method (const ran2_state_t * r, const double sigma)
{
  double u, v, x, y, Q;
  const double s = 0.449871;    /* Constants from Leva */
  const double t = -0.386595;
  const double a = 0.19600;
  const double b = 0.25472;
  const double r1 = 0.27597;
  const double r2 = 0.27846;

  do                            /* This loop is executed 1.369 times on average  */
    {
      /* Generate a point P = (u, v) uniform in a rectangle enclosing
         the K+M region v^2 <= - 4 u^2 log(u). */

      /* u in (0, 1] to avoid singularity at u = 0 */
      u = 1 - ran2_get_double (r);

      /* v is in the asymmetric interval [-0.5, 0.5).  However v = -0.5
         is rejected in the last part of the while clause.  The
         resulting normal deviate is strictly symmetric about 0
         (provided that v is symmetric once v = -0.5 is excluded). */
      v = ran2_get_double (r) - 0.5;

      /* Constant 1.7156 > sqrt(8/e) (for accuracy); but not by too
         much (for efficiency). */
      v *= 1.7156;

      /* Compute Leva's quadratic form Q */
      x = u - s;
      y = fabs (v) - t;
      Q = x * x + y * (a * y - b * x);

      /* Accept P if Q < r1 (Leva) */
      /* Reject P if Q > r2 (Leva) */
      /* Accept if v^2 <= -4 u^2 log(u) (K+M) */
      /* This final test is executed 0.012 times on average. */
    }
  while (Q >= r1 && (Q > r2 || v * v > -4 * u * u * log (u)));

  return sigma * (v / u);       /* Return slope */
}



//NORMAL CDF FIXME
//---------------

#ifndef M_1_SQRT2PI
#define M_1_SQRT2PI (M_2_SQRTPI * M_SQRT1_2 / 2.0)
#endif

#define SQRT32 (4.0 * M_SQRT2)

#define GAUSS_EPSILON  (GSL_DBL_EPSILON / 2)
#define GAUSS_XUPPER (8.572)
#define GAUSS_XLOWER (-37.519)
#define GAUSS_SCALE (16.0)

static double
get_del (double x, double rational)
{
  double xsq = 0.0;
  double del = 0.0;
  double result = 0.0;

  xsq = floor (x * GAUSS_SCALE) / GAUSS_SCALE;
  del = (x - xsq) * (x + xsq);
  del *= 0.5;

  result = exp (-0.5 * xsq * xsq) * exp (-1.0 * del) * rational;

  return result;
}

/*
 * Normal cdf for fabs(x) < 0.66291
 */
static double
gauss_small (const double x)
{
  unsigned int i;
  double result = 0.0;
  double xsq;
  double xnum;
  double xden;

  const double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
  };
  const double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
  };

  xsq = x * x;
  xnum = a[4] * xsq;
  xden = xsq;

  for (i = 0; i < 3; i++)
    {
      xnum = (xnum + a[i]) * xsq;
      xden = (xden + b[i]) * xsq;
    }

  result = x * (xnum + a[3]) / (xden + b[3]);

  return result;
}

/*
 * Normal cdf for 0.66291 < fabs(x) < sqrt(32).
 */
static double
gauss_medium (const double x)
{
  unsigned int i;
  double temp = 0.0;
  double result = 0.0;
  double xnum;
  double xden;
  double absx;

  const double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
  };
  const double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
  };

  absx = fabs (x);

  xnum = c[8] * absx;
  xden = absx;

  for (i = 0; i < 7; i++)
    {
      xnum = (xnum + c[i]) * absx;
      xden = (xden + d[i]) * absx;
    }

  temp = (xnum + c[7]) / (xden + d[7]);

  result = get_del (x, temp);

  return result;
}

/*
 * Normal cdf for 
 * {sqrt(32) < x < GAUSS_XUPPER} union { GAUSS_XLOWER < x < -sqrt(32) }.
 */
static double
gauss_large (const double x)
{
  int i;
  double result;
  double xsq;
  double temp;
  double xnum;
  double xden;
  double absx;

  const double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
  };
  const double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
  };

  absx = fabs (x);
  xsq = 1.0 / (x * x);
  xnum = p[5] * xsq;
  xden = xsq;

  for (i = 0; i < 4; i++)
    {
      xnum = (xnum + p[i]) * xsq;
      xden = (xden + q[i]) * xsq;
    }

  temp = xsq * (xnum + p[4]) / (xden + q[4]);
  temp = (M_1_SQRT2PI - temp) / absx;

  result = get_del (x, temp);

  return result;
}



double
gsl_cdf_ugaussian_Q (const double x)
{
  double result;
  double absx = fabs (x);

  if (absx < GAUSS_EPSILON)
    {
      result = 0.5;
      return result;
    }
  else if (absx < 0.66291)
    {
      result = gauss_small (x);

      if (x < 0.0)
        {
          result = fabs (result) + 0.5;
        }
      else
        {
          result = 0.5 - result;
        }

      return result;
    }
  else if (absx < SQRT32)
    {
      result = gauss_medium (x);

      if (x < 0.0)
        {
          result = 1.0 - result;
        }

      return result;
    }
  else if (x > -(GAUSS_XLOWER))
    {
      result = 0.0;
      return result;
    }
  else if (x < -(GAUSS_XUPPER))
    {
      result = 1.0;
      return result;
    }
  else
    {
      result = gauss_large (x);

      if (x < 0.0)
        {
          result = 1.0 - result;
        }

    }

  return result;
}


//GAMMA DEVIATES FIXME
//---------------

double
gsl_ran_gamma (const ran2_state_t * r, const double a, const double b)
{
	double gsl_ran_gaussian_ziggurat (const ran2_state_t* r, const double sigma);
  /* assume a > 0 */

  if (a < 1)
    {
      //double u = gsl_rng_uniform_pos (r);
      double u = 1.0 - ran2_get_double (r);
      return gsl_ran_gamma (r, 1.0 + a, b) * pow (u, 1.0 / a);
    }

  {
    double x, v, u;
    double d = a - 1.0 / 3.0;
    double c = (1.0 / 3.0) / sqrt (d);

    while (1)
      {
        do
          {
            x = gsl_ran_gaussian_ziggurat (r, 1.0);
            v = 1.0 + c * x;
          }
        while (v <= 0);

        v = v * v * v;
        u = 1.0 - ran2_get_double (r);

        if (u < 1 - 0.0331 * x * x * x * x) 
          break;

        if (log (u) < 0.5 * x * x + d * (1 - v + log (v)))
          break;
      }
    
    return b * d * v;
  }
}


//FIXME

//UINT_MAX defn should be ok

//double
//gsl_ran_gamma_knuth (const ran2_state_t* r, const double a, const double b)
//{
//  /* assume a > 0 */
//  unsigned int na = floor (a);

//  if(a >= UINT_MAX) 
//    {
//      return b * (gamma_large (r, floor (a)) + gamma_frac (r, a - floor (a)));
//    }
//  else if (a == na)
//    {
//      return b * gsl_ran_gamma_int (r, na);
//    }
//  else if (na == 0)
//    {
//      return b * gamma_frac (r, a);
//    }
//  else
//    {
//      return b * (gsl_ran_gamma_int (r, na) + gamma_frac (r, a - na)) ;
//    }
//}


//BETA DEVIATES FIXME
//---------------

double
gsl_ran_beta (const ran2_state_t* r, const double a, const double b)
{
  double x1 = gsl_ran_gamma (r, a, 1.0);
  double x2 = gsl_ran_gamma (r, b, 1.0);

  return x1 / (x1 + x2);
}


//FOR SOME REASON WE NEED GAUSSIAN ZIGGURAT FOR THE GAMMA FUNCTION
//---------------------------------------------------------------

/* position of right-most step */
#define PARAM_R 3.44428647676

/* tabulated values for the heigt of the Ziggurat levels */
static const double ytab[128] = {
  1, 0.963598623011, 0.936280813353, 0.913041104253,
  0.892278506696, 0.873239356919, 0.855496407634, 0.838778928349,
  0.822902083699, 0.807732738234, 0.793171045519, 0.779139726505,
  0.765577436082, 0.752434456248, 0.739669787677, 0.727249120285,
  0.715143377413, 0.703327646455, 0.691780377035, 0.68048276891,
  0.669418297233, 0.65857233912, 0.647931876189, 0.637485254896,
  0.62722199145, 0.617132611532, 0.607208517467, 0.597441877296,
  0.587825531465, 0.578352913803, 0.569017984198, 0.559815170911,
  0.550739320877, 0.541785656682, 0.532949739145, 0.524227434628,
  0.515614886373, 0.507108489253, 0.498704867478, 0.490400854812,
  0.482193476986, 0.47407993601, 0.466057596125, 0.458123971214,
  0.450276713467, 0.442513603171, 0.434832539473, 0.427231532022,
  0.419708693379, 0.41226223212, 0.404890446548, 0.397591718955,
  0.390364510382, 0.383207355816, 0.376118859788, 0.369097692334,
  0.362142585282, 0.355252328834, 0.348425768415, 0.341661801776,
  0.334959376311, 0.328317486588, 0.321735172063, 0.31521151497,
  0.308745638367, 0.302336704338, 0.29598391232, 0.289686497571,
  0.283443729739, 0.27725491156, 0.271119377649, 0.265036493387,
  0.259005653912, 0.253026283183, 0.247097833139, 0.241219782932,
  0.235391638239, 0.229612930649, 0.223883217122, 0.218202079518,
  0.212569124201, 0.206983981709, 0.201446306496, 0.195955776745,
  0.190512094256, 0.185114984406, 0.179764196185, 0.174459502324,
  0.169200699492, 0.1639876086, 0.158820075195, 0.153697969964,
  0.148621189348, 0.143589656295, 0.138603321143, 0.133662162669,
  0.128766189309, 0.123915440582, 0.119109988745, 0.114349940703,
  0.10963544023, 0.104966670533, 0.100343857232, 0.0957672718266,
  0.0912372357329, 0.0867541250127, 0.082318375932, 0.0779304915295,
  0.0735910494266, 0.0693007111742, 0.065060233529, 0.0608704821745,
  0.056732448584, 0.05264727098, 0.0486162607163, 0.0446409359769,
  0.0407230655415, 0.0368647267386, 0.0330683839378, 0.0293369977411,
  0.0256741818288, 0.0220844372634, 0.0185735200577, 0.0151490552854,
  0.0118216532614, 0.00860719483079, 0.00553245272614, 0.00265435214565
};

/* tabulated values for 2^24 times x[i]/x[i+1],
 * used to accept for U*x[i+1]<=x[i] without any floating point operations */
static const unsigned long ktab[128] = {
  0, 12590644, 14272653, 14988939,
  15384584, 15635009, 15807561, 15933577,
  16029594, 16105155, 16166147, 16216399,
  16258508, 16294295, 16325078, 16351831,
  16375291, 16396026, 16414479, 16431002,
  16445880, 16459343, 16471578, 16482744,
  16492970, 16502368, 16511031, 16519039,
  16526459, 16533352, 16539769, 16545755,
  16551348, 16556584, 16561493, 16566101,
  16570433, 16574511, 16578353, 16581977,
  16585398, 16588629, 16591685, 16594575,
  16597311, 16599901, 16602354, 16604679,
  16606881, 16608968, 16610945, 16612818,
  16614592, 16616272, 16617861, 16619363,
  16620782, 16622121, 16623383, 16624570,
  16625685, 16626730, 16627708, 16628619,
  16629465, 16630248, 16630969, 16631628,
  16632228, 16632768, 16633248, 16633671,
  16634034, 16634340, 16634586, 16634774,
  16634903, 16634972, 16634980, 16634926,
  16634810, 16634628, 16634381, 16634066,
  16633680, 16633222, 16632688, 16632075,
  16631380, 16630598, 16629726, 16628757,
  16627686, 16626507, 16625212, 16623794,
  16622243, 16620548, 16618698, 16616679,
  16614476, 16612071, 16609444, 16606571,
  16603425, 16599973, 16596178, 16591995,
  16587369, 16582237, 16576520, 16570120,
  16562917, 16554758, 16545450, 16534739,
  16522287, 16507638, 16490152, 16468907,
  16442518, 16408804, 16364095, 16301683,
  16207738, 16047994, 15704248, 15472926
};

/* tabulated values of 2^{-24}*x[i] */
static const double wtab[128] = {
  1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08,
  3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
  3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08,
  4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
  5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08,
  5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
  5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08,
  6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
  6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08,
  6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
  7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08,
  7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
  7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08,
  8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
  8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08,
  8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
  9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08,
  9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
  9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07,
  1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
  1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07,
  1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
  1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07,
  1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
  1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07,
  1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
  1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07,
  1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
  1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07,
  1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
  1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07,
  1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};


double
gsl_ran_gaussian_ziggurat (const ran2_state_t* r, const double sigma)
{
  unsigned long int i, j;
  int sign;
  double x, y;

  const unsigned long int range = 2147483563;
  const unsigned long int offset = 0;

  while (1)
    {
      if (range >= 0xFFFFFFFF)
        {
          unsigned long int k = ran2_get(r) - offset;
          i = (k & 0xFF);
          j = (k >> 8) & 0xFFFFFF;
        }
      else if (range >= 0x00FFFFFF)
        {
          unsigned long int k1 = ran2_get(r) - offset;
          unsigned long int k2 = ran2_get(r) - offset;
          i = (k1 & 0xFF);
          j = (k2 & 0x00FFFFFF);
        }
      else
        {
          i = gsl_rng_uniform_int (r, 256); /*  choose the step */
          j = gsl_rng_uniform_int (r, 16777216);  /* sample from 2^24 */
        }

      sign = (i & 0x80) ? +1 : -1;
      i &= 0x7f;

      x = j * wtab[i];

      if (j < ktab[i])
        break;

      if (i < 127)
        {
          double y0, y1, U1;
          y0 = ytab[i];
          y1 = ytab[i + 1];
          U1 = ran2_get_double (r);
          y = y1 + (y0 - y1) * U1;
        }
      else
        {
          double U1, U2;
          U1 = 1.0 - ran2_get_double (r);
          U2 = ran2_get_double (r);
          x = PARAM_R - log (U1) / PARAM_R;
          y = exp (-PARAM_R * (x - 0.5 * PARAM_R)) * U2;
        }

      if (y < exp (-0.5 * x * x))
        break;
    }

  return sign * sigma * x;
}

