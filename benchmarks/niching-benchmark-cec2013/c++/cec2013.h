/******************************************************************************
 * Version: 1.0
 * Last modified on: 21 January, 2013 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/
#ifndef __CEC2013_NICHING_COMPETITION_H__
#define __CEC2013_NICHING_COMPETITION_H__
#include <vector>
#include "cfunction.h"

/* Interface for all benchmarks */
class CEC2013
{
	//non-copyable
	CEC2013(const CEC2013 &);
	CEC2013& operator=(const CEC2013 &);
public:
	CEC2013(const int &nofunc);
	~CEC2013();

	tFitness evaluate(const double *x);
	tFitness evaluate(const std::vector<double> &x);
	tFitness get_fitness_goptima() const { return fopt_[nfunc_-1]; };

	double get_lbound(const int &n) const;
	double get_ubound(const int &n) const;
	
	int get_dimension() const { return dimensions_[nfunc_-1]; }
	int get_no_goptima() const { return nopt_[nfunc_-1]; }
	int how_many_goptima(const std::vector< std::vector<double> > &pop, 
		std::vector< std::vector<double> > &seeds,
		const double &accuracy);
	int count_goptima(const std::vector<std::vector<double> > &pop,
		   		      const std::vector<double> &fit,
					  double accuracy);

	double get_rho() const { return rho_[nfunc_-1]; }
	int get_maxfes() const { return maxfes_[nfunc_-1]; }
private:
	/*Specific function dimensions for the competition */
	int nfunc_;
	CFunction *cfunc_;
	int dimensions_[20];
	tFitness fopt_[20];
	double rho_[20];
	double nopt_[20];
	long maxfes_[20];

	void init_vars_();
};

void get_seeds(const std::vector< std::vector<double> > &pop, std::vector< std::vector<double> > &seeds,
		const double &radius);
double get_eudist(const std::vector<double> &v1, const std::vector<double> &v2);
void print_pop(std::vector< std::vector<double> > pop);
#endif
