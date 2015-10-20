

class Solver
{
	private: 
		int timepoints;
	protected:
		double* y0;
		double tend, t0, h;
		void(*func)(int neq, double t, double* y);


	protected:
		int neq, nelm;
	public:
		Solver(int neq, double t0, double tend, double* y0, double h, void(*func)(int neq, double t, double* y))
		{
			this->init(neq, t0, tend, y0, h, func);
		}
		~Solver(){}
	public:
		void init(int neq, double t0, double tend, double* y0, double h, void(*func)(int neq, double t, double* y));
		
		int getMemorySize();
		int getMemorySizeT();
		int solve(double *tout, double* yout);
	protected:
		virtual void step(double tlast, double* ylast, double h) = 0;

};

class ForwEuler : public Solver
{
	public:
		ForwEuler(int neq, double t0, double tend, double* y0, double h, void(*func)(int neq, double t, double* y)) 
		: Solver(neq, t0, tend, y0, h, func) {
		}
		~ForwEuler(){}
	protected:
		void step(double tlast, double* ylast, double h);
};

class RungeKutta : public Solver
{
	public:
		RungeKutta(int neq, double t0, double tend, double* y0, double h, void(*func)(int neq, double t, double* y)) 
		: Solver(neq, t0, tend, y0, h, func) {
		}
		~RungeKutta(){}
	protected:
		void step(double tlast, double* ylast, double h);
};