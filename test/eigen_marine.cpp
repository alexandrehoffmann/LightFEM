#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>

#include <LightFEM/Mesh.hpp>
#include <LightFEM/Analysis.hpp>
#include <LightFEM/Expression.hpp>

#include <chrono>
#include <iostream>

int main()
{	
	constexpr size_t degree = 4;
	
	Eigen::initParallel();
	
	//////////////////////////////////////////////////////////////////////
	////  We want to solve the following PDE                          ////
	////        $\int_\Omega \ddot{u}v + \nabla u.\nabla v dx$        ////
	////               $+ \int\partial\Omega  \dot{u}v$               ////
	////                 $= \int_\Omega v\delta_0 dx$                 ////
	//// With a P7 Finite element method                              ////
	//////////////////////////////////////////////////////////////////////

	// Initialize mesh from a Gmsh file
	// In this mesh we defined three sub domains:
	// * bathy (y > -70m)
	// * shallow (-70m > y > -2km)
	// * deep (-2km > y > -5km)
	// and two boundaries:
	// * dirichlet
	// * neumann
	GmshMesh mesh("marine.msh");
	// Initialize function space P7 defined on the mesh 
	PnFunctionSpace Vh(&mesh, degree);
	// Initialize source term, a delta function
	SamplingFunction delta(&mesh, NodeWorld(1.0, -0.02));
	// Initialize wave velocity
	ConstFactory cst = ConstFactory(&mesh);
	ScalarField c(&mesh, [](const double, const double y) -> double
	{
		if (y > -2.0)
		{
			return 1.5;
		}
		return 3.5;
	});
	// Print wave velocity using Gnuplot file format.
	// the filed can be displayed with set view 0,0; splot "err_1.dat" u 1:2:3 w pm3d
	printFunction("c.dat", c);
	
	//// We define our time discretization and our Newmark scheme 
	const double t0 = 0.0;
	const double t1 = 3.0;
	const double dt = 0.001;
	const size_t nt = (t1-t0) / dt;

	const double gamma = 0.5;
	const double beta = 0.5*gamma;
	
	// our source term can be written as $f(t,x) = s(t)\delta_0(x)$ where $f(t)$ is a ricker wavelet.
	const double f_M = 3.5;
	auto ricker = [f_M](const double t) -> double { return (1.0 - 2.0*M_PI*M_PI*f_M*f_M*t*t)*std::exp(-(M_PI*f_M*t)*(M_PI*f_M*t)); };
	
	//////////////////////////////////////////////////////////////////////
	////               Our equation can be written as:                ////
	//// A(a_{k+1}, w) = (f,w) - a_C(v_k + dt(1 - gamma)a_k,w)        ////
	////                 - a_K((uk + dt vk + dt dt(0.5 - beta)ak, w)  ////
	////       v_{k+1} = v_{k} + (1 - gamma)dt ak + gamma dt a_{k+1}  ////
	////       u_{k+1} = u_k + dt v_k + (0.5 - beta) dt dt a_k        ////
	////                 + beta dt dt a_{k+1}                         ////
	//////////////////////////////////////////////////////////////////////
	
	BilinearForm formA(&Vh, &Vh, [&mesh, &c, dt, gamma, beta, &cst](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, u*v / (c*c) + cst(beta*dt*dt)*inner(grad(u), grad(v)))
			+ boundaryIntegral(mesh, "neumann", cst(gamma*dt)*u*v / c); // the first argument is the domain, the second argument is an expression 
	});
	formA.setIdentityOnBoundary({"dirichlet"});
	BilinearForm formC(&Vh, &Vh, [&mesh, &c](const TrialFunction& u, const TestFunction& v) -> double
	{
		return boundaryIntegral(mesh, "neumann", u*v / c);  
	});
	formC.setZeroOnBoundary({"dirichlet"}, false, true);
	formC.pruneNullEntries();
	BilinearForm formK(&Vh, &Vh, [&mesh](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, inner(grad(u), grad(v))); // the first argument is the domain, the second argument is an expression 
	});
	formK.setZeroOnBoundary({"dirichlet"}, false, true);
	formK.pruneNullEntries();
	BilinearForm formM(&Vh, &Vh, [&mesh, c](const TrialFunction& u, const TestFunction& v) -> double // used to compute a_0
	{
		return integral(mesh, u*v / (c*c));
	});
	formM.setIdentityOnBoundary({"dirichlet"});
	// create a linear form $l(v) := \int_\Omega v\delta_0 dx$
	LinearForm l(&Vh, [&mesh, &delta](const TestFunction& v) -> double
	{
		return integral(mesh, v, delta); // the first argument is the domain, the second argument is an expression, the third argument is a measure (it can be a weighted sampling function or a quadrature) 
	});
	l.setZeroOnBoundary({"dirichlet"});
	//////////////////////////////////////////////////////////////////////
	////              We discretize our PDE using eigen               ////
	////         Note that any Linear algebra could be used.          ////
	//////////////////////////////////////////////////////////////////////
	
	Eigen::SparseMatrix< double > A(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	A.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formA.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formA.getCoefData() + formA.getNCoefs()));
	
	Eigen::SparseMatrix< double > M(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	M.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formM.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formM.getCoefData() + formM.getNCoefs()));
	
	Eigen::SparseMatrix< double, Eigen::RowMajor > C(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	C.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formC.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formC.getCoefData() + formC.getNCoefs()));
	
	Eigen::SparseMatrix< double, Eigen::RowMajor > K(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	K.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formK.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formK.getCoefData() + formK.getNCoefs()));
	
	Eigen::Map<const Eigen::VectorXd> f(l.getCoefData(), l.getNCoefs()); // source term f_i = $l(v_i)$
	
	Eigen::SparseLU< Eigen::SparseMatrix< double > > LLt_A(A);
	Eigen::SparseLU< Eigen::SparseMatrix< double > > LLt_M(M);
	
	//// We prepare arrays to store our solutions on the boundaries 
	
	Eigen::DiagonalMatrix<double, -1> R(Vh.getNBasisFunction());
	for (size_t globId=0;globId<Vh.getNBasisFunction();++globId)
	{
		R.diagonal()[globId] = (Vh.isIdOnBoundary(globId)) ? 0.0 : 1.0;
	}
	std::vector<size_t> locToGlobId;
	for (size_t globId=0;globId<Vh.getNBasisFunction();++globId)
	{
		if (Vh.isIdOnBoundary(globId))
		{
			locToGlobId.push_back(globId);
		}
	}
	std::vector< std::vector< double > > ub(nt+1);
	std::vector< std::vector< double > > vb(nt+1);
	std::vector< std::vector< double > > ab(nt+1);
	
	//// We initialize our scheme
	
	Eigen::VectorXd uk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); // wavefield $u(t_k)$
	Eigen::VectorXd vk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); // wavefield velocity $\dot{u}(t_k)$
	Eigen::VectorXd ak = LLt_M.solve(ricker(0)*f); // wavefield acceleration $\ddot{u}(t_k)$
	
	std::vector< double > duration(nt+1);
	
	size_t idx = 0;
	for (size_t k=0;k<nt+1;++k)
	{
		auto start = std::chrono::high_resolution_clock::now();
		
		//const double tk = dt*k + t0;
		const double tkp1 = dt*(k+1) + t0;
		
		ub[k].resize(locToGlobId.size());
		vb[k].resize(locToGlobId.size());
		ab[k].resize(locToGlobId.size());
		for (size_t locId=0;locId<locToGlobId.size();++locId)
		{
			ub[k][locId] = uk[locToGlobId[locId]];
			vb[k][locId] = vk[locToGlobId[locId]];
			ab[k][locId] = ak[locToGlobId[locId]];
		} 
		
		if (k % 10 == 0) // print the function every 100 iterations
		{
			// create a ScalarField in Vh from a list of coefficients 
			// at this point the ScalarField hasn't been discretized yet. 
			// it must be discretized before being printed, integrated or used in an expression
			FiniteElementScalarField u(&Vh, reinterpret_cast< Scalar* >(uk.data()));
			
			printFunction("u_" + std::to_string(idx) + ".dat", u.discretize());
			++idx;
		}
		if (k < nt)
		{
			const Eigen::VectorXd akp1 = LLt_A.solve( ricker(tkp1)*f - C*(vk + dt*(1.0 - gamma)*ak) - K*(uk + dt*vk + dt*dt*(0.5 - beta)*ak) );
			const Eigen::VectorXd vkp1 = vk + (1.0 - gamma)*dt*ak + gamma*dt*akp1;
			const Eigen::VectorXd ukp1 = uk + dt*vk + (0.5 - beta)*dt*dt*ak + beta*dt*dt*akp1;
			
			ak = akp1;
			vk = vkp1;
			uk = ukp1;
		}
		
		auto stop = std::chrono::high_resolution_clock::now();
		duration[k] = double(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())*1.0e-6;
		
		std::cout << "Iteration " << k << "/ " << nt << " computed in " << duration[k] << " seconds" << std::endl;
	}
	std::cout << "Average time for a single iteration: " << std::reduce(std::begin(duration), std::end(duration)) / duration.size() << std::endl;
	
	return EXIT_SUCCESS;
}
