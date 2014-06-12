import std.stdio;
import multifit_nlin;

double eval(double x, double f, double s)
{
	return (x-f)/s;
}
double funa(double x, double[] ps)
{
	//writeln(ps);
	double result = 0;
	double pow = 1;
	foreach(p; ps)
	{
		result += p * pow;
		pow *= x;
	}
	return result;
}
void main()
{
	double[6] params = 0.1;
	
	Dp!double[] data;
	data ~= Dp!double(1,0);
	data ~= Dp!double(2,1);
	data ~= Dp!double(3,2);
	data ~= Dp!double(4,1);
	data ~= Dp!double(5,2);
	data ~= Dp!double(6,0);
	data ~= Dp!double(7,1);
	data ~= Dp!double(8,2);
	data ~= Dp!double(9,1);
	data ~= Dp!double(10,2);
	
	auto fit = MultifitNlin!(double,typeof(&funa))
				(&funa, data, params, true);
	fit.run();

	foreach(line; fit.result_covar)
	{
		writeln(line);
	}
}


