#include "csm/csm_all.h"


gsl_vector * vector_from_array(unsigned int n, double *x) {
	gsl_vector * v = gsl_vector_alloc(n);
	unsigned int i;
	for(i=0;i<n;i++)
		gvs(v,i,x[i]);

	return v;
}

void vector_to_array(const gsl_vector*v, double*x){
	int i;
	for(i=0;i<v->size();i++)
		x[i] = gvg(v,i);
}

void ominus(const gsl_vector*x, gsl_vector*res) {
	double c = cos(gvg(x,2));
	double s = sin(gvg(x,2));
    gvs(res,0,  -c*gvg(x,0) - s*gvg(x,1));
    gvs(res,1,   s*gvg(x,0) - c*gvg(x,1));
	gvs(res,2,  -gvg(x,2));
}


const char* gsl_friendly_pose(gsl_vector*v) {
	return friendly_pose(v->data());
}

