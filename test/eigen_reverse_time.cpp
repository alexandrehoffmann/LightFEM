#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>

#include <LightFEM/Mesh.hpp>
#include <LightFEM/Analysis.hpp>
#include <LightFEM/Expression.hpp>

#include <chrono>
#include <iostream>

int main()
{	
	constexpr size_t degree = 7;
	
	Eigen::initParallel();
	
	//////////////////////////////////////////////////////////////////////
	////  We want to solve the following PDE                          ////
	////        $\int_\Omega \ddot{u}v + \nabla u.\nabla v dx$        ////
	////               $+ \int\partial\Omega  \dot{u}v$               ////
	////                 $= \int_\Omega v\delta_0 dx$                 ////
	//// With a P7 Finite element method                              ////
	//////////////////////////////////////////////////////////////////////

	// Initialize mesh from a Gmsh file
	GmshMesh mesh("circle.msh");
	// Initialize function space P7 defined on the mesh 
	PnFunctionSpace Vh(&mesh, degree);
	// Initialize source term, a delta function
	SamplingFunction delta(&mesh, NodeWorld(0.0, 0.0));
	// Initialize wave velocity
	ConstFactory cst = ConstFactory(&mesh);
	//Const c = cst(0.5);	
	ScalarField c(&mesh, [](const double x, const double y) -> double
	{
		if (std::sqrt((x - 0.25)*(x - 0.25) + (y - 0.25)*(y - 0.25)) < 0.1)
		{
			return 1.0;
		}
		if (std::sqrt((x + 0.25)*(x + 0.25) + (y + 0.25)*(y + 0.25)) < 0.1)
		{
			return 0.25;
		}
		return 0.5;
	});
	// Print wave velocity using Gnuplot file format.
	// the filed can be displayed with set view 0,0; splot "err_1.dat" u 1:2:3 w pm3d
	printFunction("c.dat", c);
	
	std::cout << "max(c) = " << max(c) << std::endl;
	const double dx = Vh.getHMin();
	std::cout << "h = " << dx << std::endl;
	std::cout << "theoritical dt = " << 0.5*min(cst(dx) / c ) << std::endl;
	
	//// We define our time discretization and our Newmark scheme 
	const double t0 = 0.0;
	const double t1 = 3.0;
	const double dt = 0.001;
	const size_t nt = (t1-t0) / dt;

	//const double gamma = 0.5;
	//const double beta = 0.5*gamma;
	const double gamma = 0.5;
	const double beta = 0.0*gamma;
	
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
			+ boundaryIntegral(mesh, cst(gamma*dt)*u*v / c); // the first argument is the domain, the second argument is an expression 
	});
	BilinearForm formC(&Vh, &Vh, [&mesh, &c](const TrialFunction& u, const TestFunction& v) -> double
	{
		return boundaryIntegral(mesh, u*v / c);  
	});
	BilinearForm formK(&Vh, &Vh, [&mesh](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, inner(grad(u), grad(v))); // the first argument is the domain, the second argument is an expression 
	});
	
	BilinearForm formM(&Vh, &Vh, [&mesh, c](const TrialFunction& u, const TestFunction& v) -> double // used to compute a_0
	{
		return integral(mesh, u*v / (c*c));
	});
	// create a linear form $l(v) := \int_\Omega v\delta_0 dx$
	LinearForm l(&Vh, [&mesh, &delta](const TestFunction& v) -> double
	{
		return integral(mesh, v, delta); // the first argument is the domain, the second argument is an expression, the third argument is a measure (it can be a weighted sampling function or a quadrature) 
	});
	
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
	
	//////////////////////////////////////////////////////////////////////
	////  We solve the reversed time equation:                        ////
	////              (\ddot{u_r}, w) + a(u_r,w) = (f,w)              ////
	////  With:                                                       ////
	////                       $u_r(t) = u(T-t)$ on $\partial\Omega$  ////
	////                       $u_r(0) = u(T)$                        ////
	////                 $\dot{u_r}(0) = -\dot{u}(T)$                 ////
	////  The equation is re-written as :                             ////
	////            (\ddot{u^h_r},w) + a(u^h_r,w) = (f,w)             ////
	////                                           - (\ddot{u^b_r},w) ////
	////                                           - a(u^b_r,w)       ////
	////  With:                                                       ////
	////              $u^h_r(t) = 0$ on $\partial\Omega$              ////
	////              $u^h_r(0) = u(T) - u^b_r(0)$                    ////
	////       $\dot{u^h_r}(0) = -\dot{u}(T) + \dot{u^b_r}(0)$        ////
	////                                                              ////
	////  Our equation can be written as:                             ////
	//// A(a_{k+1}, w) = (f,w)                                        ////
	////                 - a_K((uk + dt vk + dt dt(0.5 - beta)ak, w)  ////
	////       v_{k+1} = v_{k} + (1 - gamma)dt ak + gamma dt a_{k+1}  ////
	////       u_{k+1} = u_k + dt v_k + (0.5 - beta) dt dt a_k        ////
	////                 + beta dt dt a_{k+1}                         ////
	//////////////////////////////////////////////////////////////////////
	
	BilinearForm formAr(&Vh, &Vh, [&mesh, &c, dt, beta, &cst](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, u*v / (c*c) + cst(beta*dt*dt)*inner(grad(u), grad(v)));
	}); 
	formAr.setIdentityOnBoundary();
	BilinearForm formKr(&Vh, &Vh, [&mesh](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, inner(grad(u), grad(v))); 
	});
	formKr.setZeroOnBoundary(false, true); // nullyfy all rows from Kr i.e. we force v to live in H_0^1(\Omega) 
	BilinearForm formMr(&Vh, &Vh, [&mesh, c](const TrialFunction& u, const TestFunction& v) -> double // used to compute a_0
	{
		return integral(mesh, u*v / (c*c));
	});
	formMr.setIdentityOnBoundary(); // we know a_0 on the boundaries
	
	//////////////////////////////////////////////////////////////////////
	////              We discretize our PDE using eigen               ////
	////         Note that any Linear algebra could be used.          ////
	//////////////////////////////////////////////////////////////////////
	
	Eigen::SparseMatrix< double > Ar(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	Ar.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formAr.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formAr.getCoefData() + formAr.getNCoefs()));
	
	Eigen::SparseMatrix< double > Mr(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	Mr.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formMr.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formMr.getCoefData() + formMr.getNCoefs()));
	
	Eigen::SparseMatrix< double, Eigen::RowMajor > Kr(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	Kr.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(formKr.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(formKr.getCoefData() + formKr.getNCoefs()));
	
	Eigen::SparseLU< Eigen::SparseMatrix< double > > LLt_Ar(Ar);
	Eigen::SparseLU< Eigen::SparseMatrix< double > > LLt_Mr(Mr);
	
	//// We initialize our scheme
	
	Eigen::VectorXd ub_final = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); 
	Eigen::VectorXd ab_final = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); 
	
	for (size_t locId=0;locId<locToGlobId.size();++locId)
	{
		ub_final[locToGlobId[locId]] = uk[locToGlobId[locId]];
		ab_final[locToGlobId[locId]] = ak[locToGlobId[locId]];
	}
	
	Eigen::VectorXd uhk =   R*uk; // wavefield $u(t_k)$
	Eigen::VectorXd vhk = -(R*vk); // wavefield velocity $\dot{u}(t_k)$
	
	Eigen::VectorXd ahk = LLt_Mr.solve(R*(ricker(t1)*f - M*ab_final - Kr*(ub_final + uhk)));
	
	idx=0;
	for (size_t k=0;k<nt+1;++k)
	{
		auto start = std::chrono::high_resolution_clock::now();
		
		const double tkp1 = t1 - dt*(k+1);
		
		Eigen::VectorXd ubk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); 
		Eigen::VectorXd vbk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); 
		Eigen::VectorXd abk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); 
		for (size_t locId=0;locId<locToGlobId.size();++locId)
		{
			ubk[locToGlobId[locId]] = ub[nt-k][locId];
			vbk[locToGlobId[locId]] = -vb[nt-k][locId];
			abk[locToGlobId[locId]] = ab[nt-k][locId];
		}
		if (k % 10 == 0) // print the function every 100 iterations
		{
			// create a ScalarField in Vh from a list of coefficients 
			// at this point the ScalarField hasn't been discretized yet. 
			// it must be discretized before being printed, integrated or used in an expression
			
			const Eigen::VectorXd urk = uhk + ubk;
			const Eigen::VectorXd vrk = vhk + vbk;
			const Eigen::VectorXd ark = ahk + abk;
			
			FiniteElementScalarField u(&Vh, reinterpret_cast< const Scalar* >(urk.data()));
			//FiniteElementScalarField dot_u(&Vh, reinterpret_cast< const Scalar* >(vrk.data()));
			//FiniteElementScalarField ddot_u(&Vh, reinterpret_cast< const Scalar* >(ark.data()));
			
			printFunction("ur_" + std::to_string(idx) + ".dat", u.discretize());
			//printFunction("vr_" + std::to_string(idx) + ".dat", dot_u.discretize());
			//printFunction("ar_" + std::to_string(idx) + ".dat", ddot_u.discretize());
			++idx;
		}
		
		const Eigen::VectorXd ahkp1 = LLt_Ar.solve( R*(ricker(tkp1)*f - M*abk - Kr*(ubk + uhk + dt*vhk + dt*dt*(0.5 - beta)*ahk)) );
		const Eigen::VectorXd vhkp1 = vhk + (1.0 - gamma)*dt*ahk + gamma*dt*ahkp1;
		const Eigen::VectorXd uhkp1 = uhk + dt*vhk + (0.5 - beta)*dt*dt*ahk + beta*dt*dt*ahkp1;
		
		ahk = ahkp1;
		vhk = vhkp1;
		uhk = uhkp1;
		
		auto stop = std::chrono::high_resolution_clock::now();
		duration[k] = double(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())*1.0e-6;
		
		std::cout << "Iteration " << k << "/ " << nt << " computed in " << duration[k] << " seconds" << std::endl;
	}
	
	std::cout << "Average time for a single iteration: " << std::reduce(std::begin(duration), std::end(duration)) / duration.size() << std::endl;
	
	return EXIT_SUCCESS;
}
