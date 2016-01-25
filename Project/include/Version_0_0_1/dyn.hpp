
#include <thrust/device_vector.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/odeint/external/thrust/thrust.hpp>
#include <boost/numeric/odeint.hpp>

#include "TestEnv.hpp"
#include "config.hpp"
#include "sys_base_functor.hpp"


#pragma once
namespace blas = boost::numeric::ublas;
using namespace std;
using namespace boost::numeric::odeint;
using namespace blas;


const int nqFirst = 3;
const int nqLast = 6;
const int nMsize = nstate * nstate;
const int nNsize = nstate * ncontr;
const int NEQ = (nvar + 1) * nstate;
const value_type eps = 1e-2;

typedef unbounded_array<value_type> UARRAY;

typedef blas::vector< value_type > vector_type;
typedef blas::matrix< value_type > matrix_type;
typedef thrust::device_vector< value_type > state_type;

//typedef thrust::detail::normal_iterator< thrust::device_ptr < const value_type > > device_ptr_iter;
//typedef array<value_type, nsubhesse> SUBHESSEARRAY;
typedef matrix<value_type, column_major, vector_type> BlasMatrix;
typedef matrix<value_type, row_major, JACOBIARRAY> MatrixJ;
//typedef matrix<value_type, row_major, SUBHESSEARRAY> MatrixSubH;




struct state
{
	state()
	{
		/* r */
		r[0] = 0;
		r[1] = 0;
		r[2] = 0;
		/* q */
		q[0] = 1;
		q[1] = 0;
		q[2] = 0;
		q[3] = 0;
		/* v */
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
		/* omega */
		omega[0] = 0;
		omega[1] = 0;
		omega[2] = 0;
		/* u */
		u[0] = 0;
		u[1] = 0;
		u[2] = 0;
		u[3] = 0;

	}

	state(const state& d)
	{
		copy(d);	
	}
	void operator=(const state& d)
	{
		copy(d);
	}
	void copy(const state& d)
	{
		/* r */
		r[0] = d.r[0];
		r[1] = d.r[1];
		r[2] = d.r[2];
		/* q */
		q[0] = d.q[0];
		q[1] = d.q[1];
		q[2] = d.q[2];
		q[3] = d.q[3];
		/* v */
		v[0] = d.v[0];
		v[1] = d.v[1];
		v[2] = d.v[2];
		/* omega */
		omega[0] = d.omega[0];
		omega[1] = d.omega[1];
		omega[2] = d.omega[2];
		/* u */
		u[0] = d.u[0];
		u[1] = d.u[1];
		u[2] = d.u[2];
		u[3] = d.u[3];
	}
	
	value_type r[3];
	value_type q[4];
	value_type v[3];
	value_type omega[3];
	value_type u[4];
};

class dyn : public TestEnv
{
	
	value_type m_g;
	value_type m_m;
	value_type m_kT;
	value_type m_kQ;
	value_type m_d;
	value_type m_IM;
	value_type m_Iges[3];
public:

	
	struct sys_functor : sys_base_functor
	{
		
		sys_functor(value_type g, value_type m, value_type kT,
			value_type kQ, value_type d, value_type IM, 
			value_type* Iges) : sys_base_functor(g, m, kT, kQ, d, IM, Iges)
		{
	
		}
		

	template <class Tuple>
		__host__ __device__ 
		void operator()(Tuple t)
		{
			value_type r[3];
			value_type q[4];
			value_type v[3];
			value_type omega[3];
			value_type u[4];

			FUNCARRAY rhs;
			JACOBIARRAY pd;

			state inp = thrust::get<0>(t);
			matrix<value_type> M = thrust::get<1>(t);
			matrix<value_type> N = thrust::get<2>(t);

			memcpy(r, inp.r, 3 * sizeof(value_type));
			memcpy(q, inp.q, 4 * sizeof(value_type));
			memcpy(v, inp.v, 3 * sizeof(value_type));
			memcpy(omega, inp.omega, 3 * sizeof(value_type));
			memcpy(u, inp.u, 4 * sizeof(value_type));


			func(r, q, v, 
			  omega, u, rhs
			);

			jacobi(r, q, v, 
				omega, u, pd
			);

			MatrixJ J(nstate, nvar, pd);
			matrix<value_type> Jx = subrange(J, 0, nstate, 0, nstate);
			matrix<value_type> Ju = subrange(J, 0, nstate, nstate, nvar);

			matrix<value_type> Nout0 = prod(Jx, N);
			matrix<value_type, column_major> Mout = prod(Jx, M);
			matrix<value_type, column_major> Nout = Nout0 + Ju;

			thrust::get<3>(t) = rhs;
			thrust::get<4>(t) = Mout.data();
			thrust::get<5>(t) = Nout.data();
		}
	};

	struct jac_functor : sys_base_functor
	{
		jac_functor(value_type g, value_type m, value_type kT,
			value_type kQ, value_type d, value_type IM, 
			value_type* Iges) : sys_base_functor(g, m, kT, kQ, d, IM, Iges)
		{

	
		}
		__device__
		void calc(const HESSEARRAY& arrHesse, const UARRAY& Msp, matrix_type &Mi)
		{
			
			for(int k = 0; k < nstate; k++)
			{
				for(int m = 0; m < nstate; m++)
				{
					value_type tmp = 0;
					for(int l = 0; l < nstate; l++)
					{
						int index = k * nvar * nvar + l * nvar + m;
						tmp += arrHesse[index] * Msp[l];
					}
					Mi(k, m) = tmp;
				}
			}
		}

		__device__
		void getNx1(const HESSEARRAY& arrHesse, matrix_type& Nx1, int ind)
		{
			const int l = ind;
			for(int k = 0; k < nstate; k++)
			{
				for(int m = 0; m < nstate; m++)
				{
					int index = k * nvar * nvar + l * nvar + m;
					Nx1(k, m) = arrHesse[index];
				}		

			}
		}

		template <class Tuple>
			__host__ __device__
		void operator()(Tuple t)
		{
			value_type r[3];
			value_type q[4];
			value_type v[3];
			value_type omega[3];
			value_type u[4];

			JACOBIARRAY pd;
			HESSEARRAY pdd;

			const int nmx = nstate * nstate;
			const int irowMxstart = nstate;
			const int irowNxstart = nstate + nmx;

			state inp = thrust::get<0>(t);
			matrix<value_type> M = thrust::get<1>(t);
			matrix<value_type> N = thrust::get<2>(t);

			memcpy(r, inp.r, 3 * sizeof(value_type));
			memcpy(q, inp.q, 4 * sizeof(value_type));
			memcpy(v, inp.v, 3 * sizeof(value_type));
			memcpy(omega, inp.omega, 3 * sizeof(value_type));
			memcpy(u, inp.u, 4 * sizeof(value_type));

			jacobi(r, q, v, 
				omega, u, pd
			);

			hesse(r, q, v, 
				omega, u, pdd
			);

			MatrixJ J(nstate, nvar, pd);
			matrix_type Jx = subrange(J, 0, nstate, 0, nstate);


			matrix_type dfdy(NEQ, NEQ);

			for(int i = 0; i < nstate; i++)
			{
				const int startindexJx =  i*nstate;
				const int endindexJx = (i+1) * nstate;

				subrange(dfdy, startindexJx, endindexJx, startindexJx, endindexJx) = Jx;


				const int startindexM = irowMxstart + i*nstate;
				const int endindexM = irowMxstart + (i+1) * nstate;
				matrix_type Mx(nstate, nstate);
				matrix_type Msp = subrange(M, 0, nstate, i, i+1);

				calc(pdd, Msp.data(), Mx);

				subrange(dfdy, startindexM, endindexM, 0, nstate ) = Mx;

				if(i < ncontr)
				{
					const int startindexN = irowNxstart + i*nstate;
					const int endindexN = irowNxstart + (i+1) * nstate;

					matrix_type Nx(nstate, nstate);
					matrix_type Nx0(nstate, nstate);
					matrix_type Nx1(nstate, nstate);
					matrix_type Nsp = subrange(N, 0, nstate, i, i+1);
					
					calc(pdd, Nsp.data(), Nx0);
					getNx1(pdd, Nx1, i + nstate);
					Nx = Nx0 + Nx1;
					subrange(dfdy, startindexN, endindexN, 0, nstate) = Nx;
				}

			}
			 
			for(int i = nstate ; i < (nvar +1); i++)
			{
				const int startindexJx =  i*nstate;
				const int endindexJx = (i+1) * nstate;

				subrange(dfdy, startindexJx, endindexJx, startindexJx, endindexJx) = Jx;
			}
			thrust::get<3>(t) = dfdy;
		}
	};

	dyn(size_t N, const vector_type &u, value_type g, value_type m, value_type kT, 
				 value_type kQ, value_type d, value_type IM,
				 value_type* Iges) : mN(N), m_g(g), m_m(m), m_kT(kT),
				 m_kQ(kQ), m_d(d), m_IM(IM)
	{
		this->m_Iges[0] = Iges[0];
		this->m_Iges[1] = Iges[1];
		this->m_Iges[2] = Iges[2];
		
		mu = state_type(ncontr * N);
		thrust::copy(u.begin(), u.end(), mu.begin());
	}

	void pre(const vector_type& x, thrust::device_vector<state>& host_vec, 
		thrust::device_vector<BlasMatrix>& hostM,
		thrust::device_vector<BlasMatrix>& hostN
	)
	{
		int startIndex = 0;
		for(int i = 0; i < mN; i++)
		{
			/*
			vector_type::const_iterator iter 
				= x.begin() + i * (nstate + nMsize + nNsize) ;
			*/

			const int ctrIndex = ncontr * i;
			vector_type iter = subrange(x, startIndex, startIndex + nstate);

//			 r 	
			state d1;
			d1.r[0] = iter[0];
			d1.r[1] = iter[1];
			d1.r[2] = iter[2];
//			 q 
			d1.q[0] = iter[3];
			d1.q[1] = iter[4];
			d1.q[2] = iter[5];
			d1.q[3] = iter[6];
//			v 
			d1.v[0] = iter[7];
			d1.v[1] = iter[8];
			d1.v[2] = iter[9];
//			omega 
			d1.omega[0] = iter[10];
			d1.omega[1] = iter[11];
			d1.omega[2] = iter[12];
//			u 
			d1.u[0] = mu[ctrIndex];
			d1.u[1] = mu[ctrIndex + 1];
			d1.u[2] = mu[ctrIndex + 2];
			d1.u[3] = mu[ctrIndex + 3];

			startIndex += nstate;
			host_vec[i] = d1;

			vector_type tmpM = subrange(x, startIndex, startIndex + nMsize);
			startIndex = startIndex + nMsize;

			vector_type tmpN = subrange(x, startIndex, startIndex + nNsize);
			startIndex = startIndex + nNsize;
			
			BlasMatrix m(nstate, nstate, tmpM);
			BlasMatrix n(nstate, ncontr, tmpN);

			hostM[i] = BlasMatrix(nstate, nstate, tmpM);
			hostN[i] = BlasMatrix(nstate, ncontr, tmpN);
		}
	}


	void fex(const vector_type& x, vector_type& dxdt)
	{

		thrust::device_vector<state> host_vec(mN);
		thrust::device_vector<BlasMatrix> hostM(mN);
		thrust::device_vector<BlasMatrix> hostN(mN);
		thrust::device_vector<FUNCARRAY> hostFuncArray(mN);
		thrust::device_vector<UARRAY> hostMArray(mN);
		thrust::device_vector<UARRAY> hostNArray(mN);

		pre(x, host_vec, hostM, hostN);
		
		thrust::for_each(
                thrust::make_zip_iterator( 
                	thrust::make_tuple( 
                		boost::begin(host_vec),
                		boost::begin(hostM),
                		boost::begin(hostN),
                		boost::begin(hostFuncArray),
                		boost::begin(hostMArray),
                		boost::begin(hostNArray)
                	)
                ),
                thrust::make_zip_iterator( 
                	thrust::make_tuple( 
                		boost::end(host_vec),
                		boost::end(hostM),
                		boost::end(hostN),
                		boost::end(hostFuncArray),
                		boost::end(hostMArray),
                		boost::end(hostNArray)
                	)
                ),
                sys_functor(m_g, m_m, m_kT,
			m_kQ, m_d, m_IM, 
			m_Iges)
        );

 		int Inddxdt = 0;
        for(int i = 0; i < mN; i++)
        {
        	FUNCARRAY arrFunc = hostFuncArray[i];
        	UARRAY MArray = hostMArray[i];
        	UARRAY NArray = hostNArray[i];

        	for(int j = 0; j < nstate; j++) // 13
        		dxdt[Inddxdt++] = arrFunc[j];
        	for(int j = 0; j < nMsize; j++)
        		dxdt[Inddxdt++] = MArray[j];//13x13
        	for(int j = 0; j < nNsize; j++)
        		dxdt[Inddxdt++] = NArray[j];//13x4
        }
	}

	void jac(const vector_type& x, matrix_type& J)
	{
		

		thrust::device_vector<state> host_vec(mN);
		thrust::device_vector<BlasMatrix> hostM(mN);
		thrust::device_vector<BlasMatrix> hostN(mN);
		thrust::device_vector<matrix_type> hostPD(mN);


		pre(x, host_vec, hostM, hostN);

		thrust::for_each(
                thrust::make_zip_iterator( 
                	thrust::make_tuple( 
                		boost::begin(host_vec),
                		boost::begin(hostM),
                		boost::begin(hostN), 
                		boost::begin(hostPD)
                	)
                ),
                thrust::make_zip_iterator( 
                	thrust::make_tuple( 
                		boost::end(host_vec),
                		boost::end(hostM),
                		boost::end(hostN),
                		boost::end(hostPD)
                	)
                ),
                jac_functor(m_g, m_m, m_kT,
			m_kQ, m_d, m_IM, 
			m_Iges)
        );
        
		
        for(int i = 0; i < mN; i++)
        {
        	const int startindex =  i*NEQ;
			const int endindex = (i+1) * NEQ;
			matrix_type pd = hostPD[i];
        	subrange(J, startindex, endindex, startindex, endindex) = pd;
        }
	}
protected:
	size_t mN;
private:
	state_type mu;

 	static void generateData(int timepoints, const vector_type& vec, vector_type& state0, vector_type& u0)
	{
		int indx = 0;
		int indu = 0;

		matrix<value_type, row_major, vector_type> identMatrix =
				identity_matrix<value_type>(nstate);
		matrix<value_type, row_major, vector_type> zeroMatrix = 
				zero_matrix<value_type> (ncontr, nstate);

		vector_type arrIdentMatrix = identMatrix.data();
		vector_type arrZeroMatrix = zeroMatrix.data();
		
		for(int j = 0; j < timepoints; j++)
		{
			for(int i = 0; i < nstate; i++)
			{
				int tmpInd = (nvar * j) + i;
				state0[indx++] = vec[tmpInd];
			}

			for(int i = 0; i < nMsize; i++)
				state0[indx++] = arrIdentMatrix[i];
			for(int i = 0; i < nNsize; i++)
				state0[indx++] = arrZeroMatrix[i];
			
			for(int i = nstate; i < nvar; i++)
			{
				int tmpInd = (nvar * j) + i;
				u0[indu++] = vec[tmpInd];
			}
			
		}
	}

public:
/*
	static void TestJacobiFunctionality(int nintervall, value_type g,
	value_type m,
	value_type kT,
	value_type kQ,
	value_type d,
	value_type IM,
	value_type* Iges
	)
	{
		cout << "Test: Jacobi Functionality" << endl;
		int ntimepoints = (nintervall + 1);
		int nvec = ntimepoints * nvar;
		int nallcontr = ncontr * ntimepoints;
		vector_type vec(nvec);
		vector_type state0(NEQ * ntimepoints);
		vector_type u0(nallcontr);


		generator(nvec, vec);

		generateData(ntimepoints, vec, state0, u0);
		dyn Dynamic(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);
		
		matrix_type J(NEQ * ntimepoints, NEQ * ntimepoints);
		Dynamic.jac(state0, J);

		cout << "J: " << J << endl;
	}
*/
	static void TestFuncFunctionality(int nintervall, value_type g,
	value_type m,
	value_type kT,
	value_type kQ,
	value_type d,
	value_type IM,
	value_type* Iges
	)
	{
		cout << "Test: Func Functionality" << endl;
		int ntimepoints = (nintervall + 1);
		int nvec = ntimepoints * nvar;
		int nallcontr = ncontr * ntimepoints;

		INITARRAY matlabData;

		vector_type u0(nallcontr);
		u0[0] = 10000;
		u0[1] = 10000;
		u0[2] = 10000;
		u0[3] = 10000;

		dyn Dynamic(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);
		sys_functor sys(g, m,  kT,
			kQ,  d,  IM, Iges);

	
		sys.Data(matlabData);
		vector_type state0(NEQ);
		vector_type dxdt(NEQ);
		vector_type result(NEQ);

		for(int i = 0; i < NEQ; i++)
			state0[i] = matlabData[i];
		int j = 0;
		for(int i = NEQ; i < 2* NEQ; i++)
		{
			result[j] = matlabData[j];
			j++;
		}
		
		matrix_type J(NEQ * ntimepoints, NEQ * ntimepoints);
		matrix_type numJ(NEQ * ntimepoints, NEQ * ntimepoints);
		

		Dynamic.fex(state0, dxdt);
		Dynamic.jac(state0, J);
		
		//numDiff_nD(int n, dyn& Dynamic, const vector_type& vec_old, matrix_type& numJ)
		numDiff_nD(NEQ * ntimepoints, Dynamic, state0, numJ);
		matrix_type diffJ = numJ - J;
		value_type normJ = norm_inf(diffJ);
		cout << "Norm - Diff: " << normJ << endl;
	}

	static void TestJacobi(int nintervall, value_type g,
	value_type m,
	value_type kT,
	value_type kQ,
	value_type d,
	value_type IM,
	value_type* Iges
	)
	{
		cout << "Test: Jacobi" << endl;
		int ntimepoints = (nintervall + 1);
		int nallcontr = ntimepoints * ncontr;
		vector_type state0(NEQ * ntimepoints);
		vector_type u0(nallcontr);

		generator(NEQ * ntimepoints, state0);
		generator(nallcontr, u0);

		
		dyn Dynamic(ntimepoints, u0, g, m, kT, kQ, d, IM, Iges);
		matrix_type J(NEQ * ntimepoints, NEQ * ntimepoints);
		matrix_type numJ(NEQ * ntimepoints, NEQ * ntimepoints);
		Dynamic.jac(state0, J);

		numDiff_nD(NEQ * ntimepoints, Dynamic, state0, numJ);
		matrix_type diffJ = numJ - J;
		value_type normJ = norm_inf(diffJ);
		cout << "Norm - Diff: " << normJ << endl;
	}

	static void numDiff_nD(int n, dyn& Dynamic, const vector_type& vec_old, matrix_type& numJ)
	{

		for(int j = 0; j < n; j++)
		{
			
			vector_type vec_p  = vec_old;
			vector_type vec_n  = vec_old;
			
			vector_type dxdtp(n);
			vector_type dxdtn(n);


			vec_p[j] = vec_old[j] + eps;
			vec_n[j] = vec_old[j] - eps;
			
			Dynamic.fex(vec_p, dxdtp);
			Dynamic.fex(vec_n, dxdtn);

			
			for (int i = 0; i < n; i++)
			{
				value_type d = (dxdtp[i] - dxdtn[i]) /(2.0 * eps);
				numJ(i, j) = d;
			}
			
		}
	}

};