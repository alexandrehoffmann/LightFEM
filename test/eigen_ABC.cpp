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
	//ConstFactory cst = ConstFactory(&mesh);
	//Const c = cst(1.0);	
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
	// create a bilinear form $a_K(u,v) := \int_\Omega \nabla u.\nabla v dx$
	BilinearForm aK(&Vh, &Vh, [&mesh](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, inner(grad(u), grad(v))); // the first argument is the domain, the second argument is an expression 
	});
	// create a bilinear form $a_M(u,v) := \int_\Omega uv dx$
	BilinearForm aM(&Vh, &Vh, [&mesh, &c](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, u*v / (c*c)); 
	});
	// create a bilinear form $a_C(u,v) := \int_\partial\Omega uv dx$
	BilinearForm aC(&Vh, &Vh, [&mesh, &c](const TrialFunction& u, const TestFunction& v) -> double
	{
		return boundaryIntegral(mesh, u*v / c);  
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
	
	// mass matrix M_ij = $m(u_j,v_i)$
	Eigen::SparseMatrix< double, Eigen::RowMajor > M(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	M.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(aM.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(aM.getCoefData() + aM.getNCoefs()));
	
	// stiffness matrix K_ij = $a(u_j,v_i)$
	Eigen::SparseMatrix< double, Eigen::RowMajor > K(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	K.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(aK.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(aK.getCoefData() + aK.getNCoefs()));
	
	// damping matrix C_ij = $c(u_j,v_i)$
	Eigen::SparseMatrix< double, Eigen::RowMajor > C(Vh.getNBasisFunction(), Vh.getNBasisFunction()); 
	C.setFromTriplets(reinterpret_cast< const Eigen::Triplet<double>* >(aC.getCoefData()), reinterpret_cast< const Eigen::Triplet<double>* >(aC.getCoefData() + aC.getNCoefs()));
	
	Eigen::Map<const Eigen::VectorXd> f(l.getCoefData(), l.getNCoefs()); // source term f_i = $l(v_i)$
	
	Eigen::VectorXd uk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); // wavefield $u(t_k)$
	Eigen::VectorXd vk = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); // wavefield velocity $\dot{u}(t_k)$
	Eigen::VectorXd ak = Eigen::VectorXd::Zero(Vh.getNBasisFunction()); // wavefield acceleration $\ddot{u}(t_k)$
	
	// our source term can be written as $f(t,x) = s(t)\delta_0(x)$ where $f(t)$ is a ricker wavelet.
	const double f_M = 3.5;
	auto ricker = [f_M](const double t) -> double { return (1.0 - 2.0*M_PI*M_PI*f_M*f_M*t*t)*std::exp(-(M_PI*f_M*t)*(M_PI*f_M*t)); };
	
	//////////////////////////////////////////////////////////////////////
	////  We solve the following semi-discretized ODE:                ////
	////                $M\ddot{u} + Cdot{u} + Ku = f$                ////
	////  with a Newmark scheme                                       ////
	//////////////////////////////////////////////////////////////////////
	
	// define our time discretization and our Newmark scheme
	const double t0 = 0.0;
	const double t1 = 3.0;
	const double dt = 0.001;
	const size_t nt = (t1-t0) / dt;

	const double gamma = 0.5;
	
	Eigen::SparseLU<Eigen::SparseMatrix< double, Eigen::RowMajor > > LLt(M + gamma*dt*C);
	
	std::vector< double > duration(nt+1);
	
	std::ofstream out("norm_u.dat");
	
	size_t idx=0;
	for (size_t k=0;k<nt+1;++k)
	{
		auto start = std::chrono::high_resolution_clock::now();
		
		const double t = dt*k + t0;
		
		out << t << " " << uk.dot(M*uk) << std::endl;
		
		// create a ScalarField in Vh from a list of coefficients 
		// at this point the ScalarField hasn't been discretized yet. 
		// it must be discretized before being printed, integrated or used in an expression
		FiniteElementScalarField u(&Vh, reinterpret_cast< Scalar* >(uk.data()));
		if (k % 10 == 0) // print the function every 100 iterations
		{
			printFunction("u_" + std::to_string(idx) + ".dat", u.discretize());
			++idx;
		}
		
		const Eigen::VectorXd akp1 = LLt.solve( ricker(t)*f - C*(vk + dt*(1.0 - gamma)*ak) - K*(uk + dt*vk + 0.5*dt*dt*ak) );
		const Eigen::VectorXd vkp1 = vk + (1.0 - gamma)*dt*ak + gamma*dt*akp1;
		const Eigen::VectorXd ukp1 = uk + dt*vk + 0.5*dt*dt*ak;
		
		ak = akp1;
		vk = vkp1;
		uk = ukp1;
		
		auto stop = std::chrono::high_resolution_clock::now();
		duration[k] = double(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())*1.0e-6;
		
		std::cout << "Iteration " << k << "/ " << nt << " computed in " << duration[k] << " seconds" << std::endl;
	}
	
	std::cout << "Average time for a single iteration: " << std::reduce(std::begin(duration), std::end(duration)) / duration.size() << std::endl;
	
	return EXIT_SUCCESS;
}
