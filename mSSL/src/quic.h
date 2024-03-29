// quic.h
// [[Rcpp::depends(RcppArmadillo)]]

#ifndef QUIC_H
# define QUIC_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <RcppArmadillo.h>

#define EPS (double(2.22E-16))

using namespace Rcpp;
using namespace arma;
using namespace std;

/*
  Functions related to QUIC
  More or less copied exactly from the original QUIC source code
*/

//The original one
/*
extern "C" {
    void dpotrf_(char* uplo, ptrdiff_t* n, double* A, ptrdiff_t* lda,
		 ptrdiff_t* info);
    void dpotri_(char* uplo, ptrdiff_t* n, double* A, ptrdiff_t* lda,
		 ptrdiff_t* info);
}
*/

typedef struct
{
	int i;
	int j;
} ushort_pair_t;

int IsDiag(int q, arma::mat &A)
{
	int flag = 0;
	for (int k = 0, i = 0; i < q; i++, k += q)
	{
		for (int j = 0; j < i; j++)
		{
			if (A(k + j) != 0.0)
			{
				flag = 0;
				;
			}
			else
			{
				flag = 1;
				;
			}
		}
	}
	return flag;
}


void CoordinateDescentUpdate(
    const int &q,                       // dimension
    const arma::mat &S,           // Y^TY/n
    const arma::mat &Rho,         // penalty
    const arma::mat &Omega,       // precision matrix, X
    const arma::mat &Sigma,       // covariance matrix, W=X^-1
    arma::mat &U,                 // DW
    arma::mat &D,                 // Newton directions, to be updated
    int i, int j,                 // looks like coordinates
    double &normD, double &diffD) // something related to the direction
{
    // calculating a
    double a = Sigma(i, j) * Sigma(i, j);                       // this is the W_ij^2
    if (i != j)
    {
        a += Sigma(i, i) * Sigma(j, j); // W_ij^2+W_iiW_jj when it is not diagonal
    }                                                          // not diagonal

    double ainv = 1.0 / a; // multiplication is cheaper than division

    // calculating b
    double b = S(i, j) - Sigma(i, j) +
               as_scalar(Sigma.row(i) * U.col(j));     // this is from GLASSO
               

    // calculating c
    double c = Omega(i, j) + D(i, j);

    // the soft threshold function related
    double l = Rho(i, j) * ainv;
    double f = b * ainv;
    double mu;
    normD -= fabs(D(i, j));
    // calculate the direction mu
    if (c > f)
    {
        mu = -f - l;
        if (c + mu < 0.0)
        {
            mu = -c;
            D(i, j) = -Omega(i, j);
        }
        else
        {
            D(i, j) += mu;
        }
    }
    else
    {
        mu = -f + l;
        if (c + mu > 0.0)
        {
            mu = -c;
            D(i, j) = -Omega(i, j);
        }
        else
        {
            D(i, j) += mu;
        }
    }
    diffD += fabs(mu);
    normD += fabs(D(i, j));
    D(j, i) = D(i, j); // make symetric
    // updating U, Q is not gonna change here
    if (mu != 0.0)
    {
        U.row(i) += mu * Sigma.row(j);
        if (i != j)
        {
            U.row(j) += mu * Sigma.row(i);
        }
    }
}


double DiagNewton(int q, const arma::mat &S,
				  const arma::mat &Rho, const arma::mat &Omega,
				  const arma::mat &Sigma, arma::mat &D)
{
	for (int iq = 0, i = 0; i < q; i++, iq += q)
	{
		for (int jq = 0, j = 0; j < i; j++, jq += q)
		{
			int ij = iq + j;
			double a = Sigma(i,i) * Sigma(j,j);
			double ainv = 1.0 / a; // multiplication is cheaper than division
			double b = S(i,j);
			double l = Rho(i,j) * ainv;
			double f = b * ainv;
			double mu;
			double x = -b * ainv;
			if (0 > f)
			{
				mu = -f - l;
				x -= l;
				if (mu < 0.0)
				{
					mu = 0.0;
					D(i,j) = -Omega(i,j);
				}
				else
				{
					D(i,j) += mu;
				}
			}
			else
			{
				mu = -f + l;
				if (mu > 0.0)
				{
					mu = 0.0;
					D(i,j) = -Omega(i,j);
				}
				else
				{
					D(i,j) += mu;
				}
			}
		}
	}
	double logdet = 0.0;
	double l1normX = 0.0;
	double trSX = 0.0;
	for (int i = 0, k = 0; i < q; i++, k += (q + 1))
	{
		logdet += log(Omega(k));
		l1normX += fabs(Omega(k)) * Rho(k);
		trSX += Omega(k) * S(k);
		double a = Sigma(k) * Sigma(k);
		double ainv = 1.0 / a; // multiplication is cheaper than division
		double b = S(k) - Sigma(k);
		double l = Rho(k) * ainv;
		double c = Omega(k);
		double f = b * ainv;
		double mu;
		if (c > f)
		{
			mu = -f - l;
			if (c + mu < 0.0)
			{
				D(k) = -Omega(k);
				continue;
			}
		}
		else
		{
			mu = -f + l;
			if (c + mu > 0.0)
			{
				D(k) = -Omega(k);
				continue;
			}
		}
		D(k) += mu;
	}
	D += D.t();
	D.diag()/=2.;
	double fX = -logdet + trSX + l1normX;
	return fX;
}

double projLogDet(int q, const arma::mat &S,
				  arma::mat &Sigma, arma::mat prW, const arma::mat &Rho)
{
	// The computed Sigma does not satisfy |Sigma - S| .< Rho.  Project it.
	for (int i = 0; i < q; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double tmp = Sigma(i,j);
			if (S(i,j) - Rho(i,j) > tmp)
				tmp = S(i,j) - Rho(i,j);
			if (S(i,j) + Rho(i,j) < tmp)
				tmp = S(i,j) + Rho(i,j);
			prW(i,j) = tmp;
			prW(j,i) = tmp;
		}
	}
	//  ptrdiff_t info = 0;
	//ptrdiff_t p0 = q;
	//dpotrf_((char*) "U", &p0, prW, &p0, &info);
	bool info;
	//int p0 = q;
	//dpotrf_((char *)"U", &p0, prW, &p0, &info);

	info = chol(prW, prW);

	if (prW.n_rows == 0)
		return 1e+15;
	double logdet = 0.0;
	for (int i = 0, k = 0; i < q; i++, k += (q + 1))
		logdet += log(prW(k));
	logdet += logdet;
	return logdet;
}


cube my_quic(const int &q, 
const arma::mat &S, const arma::mat &Rho, 
const double &tol, 
int &max_iter_quic)
{
	srand(1);
    

	int maxNewtonIter = max_iter_quic;
	double cdSweepTol = 0.05;
	int max_lineiter = 20;
	double fX = 1e+15;
	double fX1 = 1e+15;
	double fXprev = 1e+15;
	double sigma = 0.001;
	bool info;

	arma::mat Omega(q, q, fill::eye);
	arma::mat Sigma(q, q, fill::eye);

	arma::mat D(q, q, fill::zeros);
	arma::mat U(q, q, fill::zeros);


	ushort_pair_t *activeSet = (ushort_pair_t *)malloc(q * (q + 1) / 2 * sizeof(ushort_pair_t));

	double l1normX = 0.0;
	double trSX = 0.0;
	double logdetX = 0.0;

	for (int i = 0, k = 0; i < q; i++, k += q)
	{
		for (int j = 0; j < i; j++)
		{
			l1normX += Rho(k + j) * fabs(Omega(k + j));
			trSX += Omega(k + j) * S(k + j);
		}
	}
	l1normX *= 2.0;
	trSX *= 2.0;
	for (int i = 0, k = 0; i < q; i++, k += (q + 1))
	{
		l1normX += Rho(k) * fabs(Omega(k));
		trSX += Omega(k) * S(k);
	}

	int NewtonIter = 1;
	for (; NewtonIter <= maxNewtonIter; NewtonIter++)
	{
		double normD = 0.0;
		double diffD = 0.0;
		double subgrad = 1e+15;
		

		if (NewtonIter == 1 && IsDiag(q, Omega))
		{
			D *= 0.;
			fX = DiagNewton(q, S, Rho, Omega, Sigma, D);
		}
		else
		{
			
			int numActive = 0;
			U *= 0.;
			D *= 0.;
			subgrad = 0.0;

			for (int k = 0, i = 0; i < q; i++, k += q)
			{
				for (int j = 0; j <= i; j++)
				{
					double g = S(k + j) - Sigma(k + j);
					if (Omega(k + j) != 0.0 || fabs(g) > Rho(k + j))
					{
						activeSet[numActive].i = i;
						activeSet[numActive].j = j;
						numActive++;
						if (Omega(k + j) > 0)
						{
							g += Rho(k + j);
						}
						else if (Omega(k + j) < 0)
						{
							g -= Rho(k + j);
						}
						else
						{
							g = fabs(g) - Rho(k + j);
						}
						subgrad += fabs(g);
					}
				}
			}
			
			//Rcout << "Newton iteration" <<  NewtonIter << endl;
			//Rcout << "Active set size" << numActive << endl;
			//Rcout << "subgradient = " << subgrad << "l1-norm of Omega = " << l1normX << endl;

			for (int cdSweep = 1; cdSweep <= 1 + NewtonIter / 3; cdSweep++)
			{
				diffD = 0.0;
				for (int i = 0; i < numActive; i++)
				{
					int j = i + rand() % (numActive - i);
					int k1 = activeSet[i].i;
					int k2 = activeSet[i].j;
					activeSet[i].i = activeSet[j].i;
					activeSet[i].j = activeSet[j].j;
					activeSet[j].i = k1;
					activeSet[j].j = k2;
				}
				for (int l = 0; l < numActive; l++)
				{
					int i = activeSet[l].i;
					int j = activeSet[l].j;
					CoordinateDescentUpdate(q, S, Rho, Omega, Sigma, U, D, i, j, normD, diffD);
				}
				if (diffD <= normD * cdSweepTol)
					break;
			}
		}
		if (fX == 1e+15)
		{
			//ptrdiff_t info = 0;
			//ptrdiff_t p0 = q;
			
			int p0 = q;
			U = Omega;
			
			//dpotrf_((char *)"U", &p0, U, &p0, &info);
			chol(U,U);
			
			if (U.n_rows == 0)
			{
				// lack of positive definiteness
				//iter = -1;
				free(activeSet);
				U = Omega; 
				//free(U);
				//free(D);
				//return;
			}
			for (int i = 0, k = 0; i < q; i++, k += (q + 1))
			{
				logdetX += log(U(k));
			}
			logdetX *= 2.0;
			fX = (trSX + l1normX) - logdetX;
		}
		
		double trgradgD = 0.0;
		for (int i = 0, k = 0; i < q; i++, k += q)
		{
			for (int j = 0; j < i; j++)
			{
				trgradgD += (S(k + j) - Sigma(k + j)) * D(k + j);
			}
		}
		trgradgD *= 2.0;
		for (int i = 0, k = 0; i < q; i++, k += q)
		{
			trgradgD += (S(k) - Sigma(k)) * D(k);
		}
		
		/*
		originally, I had omitted lines 373 - 375 from QUIC.cpp
		so we only updated trgradgD with the diagonal
		for(int i = 0, k = 0; i < q; i++, k += (q+1)){
			trgradgD += (S[k] - Sigma[k])*D[k];
		}
		*/
		double alpha = 1.0;
		double l1normXD = 0.0;
		double fX1prev = 1e+15;
		for (int lineiter = 0; lineiter < max_lineiter; lineiter++)
		{
			double l1normX1 = 0.0;
			double trSX1 = 0.0;
			//Sigma = Omega;
			Sigma = Omega + alpha * D; 
			arma::mat temp = Sigma;
			
			
			for (int i = 0, k = 0; i < q; i++, k += q)
			{
				for (int j = 0; j < i; j++)
				{
					int ij = k + j;
					//Sigma[ij] = Omega[ij] + D[ij] * alpha;
					l1normX1 += fabs(Sigma(ij)) * Rho(ij);
					trSX1 += Sigma(ij) * S(ij);
				}
			}
			
			l1normX1 *= 2.0;
			trSX1 *= 2.0;
			for (int i = 0, k = 0; i < q; i++, k += (q + 1))
			{
				//Sigma[k] = D[k] * alpha + Omega[k];
				l1normX1 += fabs(Sigma(k)) * Rho(k);
				trSX1 += Sigma(k) * S(k);
			}
			
			//ptrdiff_t info = 0;
			//ptrdiff_t p0 = q;
			int p0 = q;
			chol(Sigma,Sigma);
			if (Sigma.n_rows==0)
			{
				alpha *= 0.5;
				Sigma = Omega;
				continue;
			}

			double logdetX1 = 0.0;
			for (int i = 0, k = 0; i < q; i++, k += (q + 1))
			{
				logdetX1 += log(Sigma(k));
			}
			logdetX1 *= 2.0;
			fX1 = (trSX1 + l1normX1) - logdetX1;
			if (alpha == 1.0)
			{
				l1normXD = l1normX1;
			}
			if (fX1 <= fX + alpha * sigma * (trgradgD + l1normXD - l1normX) || normD == 0)
			{
				fXprev = fX;
				fX = fX1;
				l1normX = l1normX1;
				logdetX = logdetX1;
				trSX = trSX1;
				Omega = temp;
				break;
			}
			if (fX1prev < fX1)
			{
				fXprev = fX;
				l1normX = l1normX1;
				logdetX = logdetX1;
				trSX = trSX1;
				Omega = temp;
				break;
			}
			fX1prev = fX1;
			alpha *= 0.5;
		}
		
		
		// compute Sigma = Omega^-1
		
		//ptrdiff_t info;
		//ptrdiff_t p0 = q;
		//int info = 0;
		int p0 = q;
		inv_sympd(Sigma, Omega);
		

		
		// check convergence
		//Rcout << "alpha "<<alpha<< " l1X " <<l1normX << " fX " << fX <<" relative change " << fabs((fX - fXprev) / fX )<< endl;
		
		if (subgrad * alpha >= l1normX * tol && (fabs((fX - fXprev) / fX) >= EPS))
			continue;
		if (NewtonIter == maxNewtonIter)
		{
			Rcout << "hit max newton Iter within QUIC" << endl;
		}
		break;
	}
	//opt = fX;
	//	double logdetW = projLogDet(q, S, Sigma, U, Rho);
	//double gap = -logdetW - q - logdetX + trSX + l1normX;
	//dGap = gap;
	//iter = NewtonIter;

	free(activeSet);
	//free(U);
	//free(D);
	
	cube results(q, q, 2);
	results.slice(0) = Omega;
	results.slice(1) = Sigma;
	return results;

	//double elapsedTime = (clock() - timeBegin)/CLOCKS_PER_SEC;
	//if(cputime != NULL) cputime = elapsedTime;
}

#endif
