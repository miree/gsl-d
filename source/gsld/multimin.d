import std.stdio;
import std.math;
import gsl.gsl_multimin;
import gsl.gsl_vector;
import gsl.gsl_deriv;

/+ F function to minimize
 +/
struct Fs(F) 
{
	F f;
	double h;	// step size during derivative calculation
}

struct LenPtr(T) { typeof(T[0].length) len; T* ptr; }
union Conv(T) {	LenPtr!T lenptr; T array[]; }
Conv!double conv_global;
// convert gsl_vector to D array
double[] conv(gsl_vector *x)
{
	conv_global.lenptr.len = x.size;
	conv_global.lenptr.ptr = gsl_vector_const_ptr(x,0);
	return conv_global.array;
}

// function for the gsl
extern(C) double min_f_pp(F) (gsl_vector *x, void *params)
{
	auto fs = cast(Fs!F *) params;
	double xval[] = conv(x);
	return fs.f(xval);
}

// helper struct to build a gsl-derivative-calculation compatible function
// F function type
struct Dms(F) // drivative structure
{
	F f;						// function
	double xval[];            	// parameter array
	uint n; 					// the n-th parameter will be modified
};

// gsl derivative of a function F with respect to a certain parameter
extern(C) double min_deriv_f(F)(double x, void *p)
{
	auto ds = cast(Dms!F *) p;
	double tmp = ds.xval[ds.n];				// save n-th parameter
	ds.xval[ds.n] = x;						// overwrite n-th parameter with x
	double f = ds.f(ds.xval);
	ds.xval[ds.n] = tmp;	                // put the n-th parameter back on its place
	return f;
}

extern(C) void min_df_pp(F)(gsl_vector *x, void *params, gsl_vector *g)
{
	auto fs = cast(Fs!F *) params;
	double xval[] = conv(x);
	Dms!(F) ds = {fs.f, xval, 0};
	
	foreach(i, xi; xval)
	{
		ds.f = fs.f;
		ds.n = cast(uint)(i);
		gsl_function gslf = {&min_deriv_f!F, &ds};
		double result, abserr;
		gsl_deriv_central(&gslf, ds.xval[i], fs.h, &result, &abserr);
		gsl_vector_set(g, i, result);
	}
}

extern(C) void min_fdf_pp(F)(gsl_vector *x, void *params, double *f, gsl_vector *g)
{
	*f = min_f_pp!F(x,params);
	min_df_pp!F(x,params,g);
}

struct Multimin(F)
{
	Fs!F fs;
	gsl_multimin_function_fdf gslf;
	gsl_multimin_fdfminimizer *s;
	gsl_vector *x;
	ulong n;
	this(F f, double x_start[], double stepsize = 0.01, double tolerance = 1e-4, double hstep = 1e-10)
	{
		this.n = x_start.length;
		this.fs = Fs!F(f,hstep);
		this.gslf.n = x_start.length;
		this.gslf.f = &min_f_pp!F;
		this.gslf.df = &min_df_pp!F;
		this.gslf.fdf = &min_fdf_pp!F;
		this.gslf.params = cast(void*)&fs;
		this.x = gsl_vector_alloc(x_start.length);
		foreach(i,xi; x_start) gsl_vector_set(x,i,xi);
		this.s = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, n);
		gsl_multimin_fdfminimizer_set(s,&gslf,x,stepsize,tolerance);
	}
	~this()
	{
		gsl_vector_free(x);
		gsl_multimin_fdfminimizer_free(s);
	}
	
	int run(int steps = 50, double epsilon = 1e-5)
	{
		int iter = 0;
		int status;
		do
		{
			status = gsl_multimin_fdfminimizer_iterate(s);
			if (status)
				break;
			status = gsl_multimin_test_gradient(s.gradient, epsilon);
		}
		while(status == GSL_CONTINUE && iter < steps);
		return iter;
	}
	@property double[] x_min()
	{
		auto result = new double[n];
		foreach(i,ref xi;result)
			xi = gsl_vector_get(s.x,i);
		return result;
	}
}
