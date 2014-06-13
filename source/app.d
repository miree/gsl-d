// Simple example of how to use the wrapper structs.
// Copyright (C) 2014  Michael Reese
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// 
// questions & comments: michaelate at gmail.com

import std.stdio;
import std.math;
import multifit_nlin;
import multimin;

double eval(double x, double f, double s)
{
	return (x-f)/s;
}
double funa(double x, double[] ps)
{
	double result = 0;
	double pow = 1;
	foreach(p; ps)
	{
		result += p * pow;
		pow *= x;
	}
	return result;
}
void multifit_nlin_test()
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
//----------------------------------------------------------------------

//double min_fun(double x[])
//{
//	double result = 0;
//	foreach(xi;x) result += sin(xi+1)^^2;
//	return result+1;
//}
//void multimin_test()
//{
//	double x_start[] = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
//	auto min = Multimin!(typeof(&min_fun))
//			(&min_fun, x_start);
//	min.run();
//	writeln("min: ", min.x_min);
//}

//----------------------------------------------------------------------
void main()
{
	multifit_nlin_test();
	//multimin_test();
}


