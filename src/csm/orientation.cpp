#include "csm/csm_all.h"
#include "egsl/egsl_macros.h"

#include <math.h>
#include <vector>


void find_neighbours(LDP ld, int i, int max_num, std::vector<int>& indexes, size_t*num_found);
void filter_orientation(double theta0, double rho0, size_t n,
 	const std::vector<double>& thetas, const std::vector<double>& rhos, double *alpha, double*cov0_alpha );

/** Requires the "cluster" field to be set */
void ld_compute_orientation(LDP ld, int size_neighbourhood, double sigma) {
	int i;
    for(i=0; i<ld->nrays; i++){
		if(!ld_valid_ray(ld,i) || (ld->cluster[i] == -1)) {
			ld->alpha[i] = std::numeric_limits<double>::quiet_NaN(); //GSL_NAN;
			ld->cov_alpha[i] = std::numeric_limits<double>::quiet_NaN(); //GSL_NAN;
			ld->alpha_valid[i] = 0;
			continue;
		}
		
        // int neighbours[size_neighbourhood*2];
		std::vector<int> neighbours(size_neighbourhood*2, 0);
		size_t num_neighbours;
		find_neighbours(ld, i, size_neighbourhood, neighbours, &num_neighbours);

        if(0 == num_neighbours) {
			ld->alpha[i] = std::numeric_limits<double>::quiet_NaN(); //GSL_NAN;
			ld->cov_alpha[i] = std::numeric_limits<double>::quiet_NaN(); //GSL_NAN;
			ld->alpha_valid[i] = 0;
			continue;
		}

		//double thetas[num_neighbours];
		std::vector<double> thetas(num_neighbours, 0.0);
		//double readings[num_neighbours];
		std::vector<double> readings(num_neighbours, 0.0);
        size_t a = 0;
        for(a=0; a<num_neighbours; a++) {
			thetas[a] = ld->theta[neighbours[a]];
			readings[a] = ld->readings[neighbours[a]];
		}
		
        double alpha = 42, cov0_alpha = 32;
        filter_orientation(ld->theta[i], ld->readings[i], num_neighbours,
            thetas, readings, &alpha, &cov0_alpha);
#ifndef WINDOWS
		if(std::isnan(alpha)) {
#else
		if(_isnan(alpha)) {
#endif
			ld->alpha[i] = std::numeric_limits<double>::quiet_NaN();
			ld->cov_alpha[i] = std::numeric_limits<double>::quiet_NaN();
			ld->alpha_valid[i] = 0;
		} else {  
			ld->alpha[i] = alpha;
			ld->cov_alpha[i] = cov0_alpha * square(sigma);
			ld->alpha_valid[i] = 1;
		}
	}
}

/** A very cool algorithm for finding the orientation */
void filter_orientation(double theta0, double rho0, size_t n,
    const std::vector<double>& thetas, const std::vector<double>& rhos, double *alpha, double*cov0_alpha ) {

    egsl_push();
    /* Y = L x + R epsilon */
    val Y = zeros(n, 1);
    val L = ones(n, 1);
    val R = zeros(n ,n+1);

    for(size_t i=0; i<n; ++i) {
        *egsl_atmp(Y, i, 0)   = (rhos[i] - rho0) / (thetas[i] - theta0);
        *egsl_atmp(R, i, 0)   =               -1 / (thetas[i] - theta0);
        *egsl_atmp(R, i, i+1) =               +1 / (thetas[i] - theta0);
    }

    val eRinv = inv(m(R, tr(R)));
    val vcov_f1 = inv(m3(tr(L), eRinv, L));
    val vf1 =   m4(vcov_f1, tr(L), eRinv, Y);

    double cov_f1 = *egsl_atmp(vcov_f1, 0, 0);
    double f1 = *egsl_atmp(vf1, 0, 0);

    *alpha = theta0 - atan(f1/rho0);

    if(cos(*alpha)*cos(theta0) + sin(*alpha)*sin(theta0) > 0)
        *alpha = *alpha + M_PI;

    double dalpha_df1  = rho0 / (square(rho0) + square(f1));
    double dalpha_drho = -f1 /  (square(rho0) + square(f1));

    *cov0_alpha	= square(dalpha_df1) * cov_f1 + square(dalpha_drho);

#ifndef WINDOWS
    if(std::isnan(*alpha)) {
#else
    if(_isnan(*alpha)) {
#endif
        egsl_print("Y", Y);
        egsl_print("L", L);
        egsl_print("R", R);
        egsl_print("eRinv", eRinv);
        egsl_print("vcov_f1", vcov_f1);

        printf("   f1 = %f cov =%f \n", f1, cov_f1);
        printf("   f1/rho = %f \n", f1/rho0);
        printf("   atan = %f \n", atan(f1/rho0));
        printf("   theta0= %f \n", theta0);
    }

    egsl_pop();
}

/* indexes: an array of size "max_num*2" */
void find_neighbours(LDP ld, int i, int max_num, std::vector<int>& indexes, size_t*num_found) {
	*num_found = 0;
	
	int up = i; 
    while ((up+1 <= i+max_num) && (up+1 < ld->nrays) && ld_valid_ray(ld,up+1)
			&& (ld->cluster[up+1] == ld->cluster[i])) {
        up++;
		indexes[(*num_found)++] = up;
	}
	int down = i; 
    while ((down >= i-max_num) && (down-1 >= 0) && ld_valid_ray(ld,down-1) &&
			(ld->cluster[down-1] == ld->cluster[i])) {
        down--;
		indexes[(*num_found)++] = down;
	}
}


