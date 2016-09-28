#include "util.hpp"
#include "geom.hpp"
#include "data.hpp"

Data	Data::model[M];
float	Data::shrink;
float	Data::bondCA, Data::bondPP, Data::distH = 1.0;
float	Data::scalein, Data::scaleout;
float	Data::delay = 0.1; // sec.s
int	Data::depth, Data::nmodels;
int	Data::norun = 0, Data::noview = 0, Data::putpdb = 0;
int	Data::frame, Data::frozen = -1, Data::hidden = 0, Data::tinker = 0, Data::unweight = 0;
Vec	Data::focus = Vec(0,0,0);
float	Data::Eratio[E] = { 0.1,1.0/10.0,2.0/10.0,3.0/10.0,4.0/10.0,5.0/10.0,6.0/10.0,7.0/10.0,8.0/10.0,9.0/10.0,
                            1.0,10.0/9.0,10.0/8.0,10.0/7.0,10.0/6.0,10.0/5.0,10.0/4.0,10.0/3.0,10.0/2.0,10.0};
