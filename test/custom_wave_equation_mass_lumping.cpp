#include <LightFEM/Mesh.hpp>
#include <LightFEM/Analysis.hpp>
#include <LightFEM/Expression.hpp>

#include <chrono>

#include <MySparseMatrix.hpp>

int main()
{
	constexpr size_t degree = 7;
	
	//////////////////////////////////////////////////////////////////////
	////  We want to solve the following PDE                          ////
	////       $\int_\Omega \ddot{u}v/c + \nabla u.\nabla v dx$       ////
	////                 $= \int_\Omega v\delta_0 dx$                 ////
	//// With a P7 Finite element method                              ////
	//////////////////////////////////////////////////////////////////////

	// Initialize mesh from a Gmsh file
	GmshMesh mesh("circle.msh");
	// Initialize function space P7 defined on the mesh 
	PnFunctionSpace Vh(&mesh, degree);
	// Initialize a 8 node Gauss-Lobatto-Legendre quadrature.
	// With this quadrature, P7 polynomials are orthogonal to each other.
	GaussLobattoQuadrature quadrature(&mesh, degree+1);
	// Initialize source term, a delta function
	SamplingFunction delta(&mesh, NodeWorld(0.0, 0.0));
	// Initialize wave velocity
	ScalarField c(&mesh, [](const double x, const double y) -> double
	{
		if (std::sqrt((x - 0.25)*(x - 0.25) + (y - 0.25)*(y - 0.25)) < 0.1)
		{
			return 2.0;
		}
		if (std::sqrt((x + 0.25)*(x + 0.25) + (y + 0.25)*(y + 0.25)) < 0.1)
		{
			return 0.5;
		}
		return 1.0;
	});
	// Print wave velocity using Gnuplot file format.
	// the filed can be displayed with set view 0,0; splot "err_1.dat" u 1:2:3 w pm3d
	printFunction("c.dat", c);
	// create a bilinear form $a(u,v) := \int_\Omega \nabla u.\nabla v dx$
	BilinearForm a(&Vh, &Vh, [&mesh](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, inner(grad(u), grad(v))); // the first argument is the domain, the second argument is an expression 
	});
	// create a bilinear form $m(u,v) := \int_\Omega uv/c dx$
	BilinearForm m(&Vh, &Vh, [&mesh, &c, &quadrature](const TrialFunction& u, const TestFunction& v) -> double
	{
		return integral(mesh, u*v / (c*c), quadrature); // Here we use our Gauss-Lobatto-Legendre quadrature to compute our expression
		// u and v are interpolating functions defined on a 8 node Gauss-Lobatto-Legendre quadrature, $\int_\Omega u_i v_j$ = $w_i\delta_ij$ where
		// $w_i$ is a weight and $\delta_ij = 1$ if i=j 0 if i!=j
	});
	// create a linear form $l(v) := \int_\Omega v\delta_0 dx$
	LinearForm l(&Vh, [&mesh, &delta](const TestFunction& v) -> double
	{
		return integral(mesh, v, delta); // the first argument is the domain, the second argument is an expression, the third argument is a measure (it can be a weighted sampling function or a quadrature) 
	});
	// enforces Homogeneous Dirichlet boundary condition on the domain
	a.setIdentityOnBoundary();
	m.setZeroOnBoundary();
	l.setZeroOnBoundary();
	
	//////////////////////////////////////////////////////////////////////
	////We discretize our PDE using the STL and a custom SparseMatrix ////
	////         Note that any Linear algebra could be used.          ////
	//////////////////////////////////////////////////////////////////////
	
	// mass matrix M_ij = $m(u_j,v_i)$
	std::vector< double > M(Vh.getNBasisFunction());
	#pragma omp parallel for
	for (size_t idx=0;idx<m.getNCoefs();++idx)
	{
		const auto& [i,j,mij] = m.getCoef(idx);
		if (i == j) { M[i] = mij; }
	}
	
	// stiffness matrix K_ij = $a(u_j,v_i)$
	MySparseMatrix<double> K(Vh.getNBasisFunction(), Vh.getNBasisFunction());
	std::vector<MySparseMatrix<double>::Entry> Kcoefs(a.getNCoefs());
	for (size_t idx=0;idx<a.getNCoefs();++idx)
	{
		const auto& [i,j,Kij] = a.getCoef(idx);
		Kcoefs.push_back(std::make_tuple(i, j, Kij));
	}
	K.setFromTriplets(std::cbegin(Kcoefs), std::cend(Kcoefs));
	
	std::vector< double > f(l.getCoefData(), l.getCoefData() + l.getNCoefs()); // source term f_i = $l(v_i)$
	
	// Compute the inverse of our mass matrix
	std::vector< double > invSemM(Vh.getNBasisFunction());
	#pragma omp parallel for
	for (size_t i=0;i<Vh.getNBasisFunction();++i)
	{
		invSemM[i] = 1.0 / M[i];
	}
	
	std::vector< double > uk(Vh.getNBasisFunction()); // wavefield $u(t_k)$
	std::vector< double > vk(Vh.getNBasisFunction()); // wavefield velocity $\dot{u}(t_k)$
	std::vector< double > ak(Vh.getNBasisFunction()); // wavefield acceleration $\ddot{u}(t_k)$
	
	std::vector< double > ukp12(Vh.getNBasisFunction()); // $uk + dt*vk + 0.5*dt*dt*ak$
	
	// our source term can be written as $f(t,x) = s(t)\delta_0(x)$ where $f(t)$ is a ricker wavelet.
	const double f_M = 3.5;
	auto ricker = [f_M](const double t) -> double { return (1.0 - 2.0*M_PI*M_PI*f_M*f_M*t*t)*std::exp(-(M_PI*f_M*t)*(M_PI*f_M*t)); };
	
	//////////////////////////////////////////////////////////////////////
	////  We solve the following semi-discretized ODE:                ////
	////                     $M\ddot{u} + Ku = f$                     ////
	////  with a Newmark scheme                                       ////
	//////////////////////////////////////////////////////////////////////
	
	// define our time discretization and our Newmark scheme
	const double t0 = 0.0;
	const double t1 = 1.25;
	const double dt = 0.0001;
	const size_t nt = (t1-t0) / dt;

	const double gamma = 0.5;
	
	std::vector< double > duration(nt+1);
	
	size_t idx=0;
	for (size_t k=0;k<nt+1;++k)
	{
		auto start = std::chrono::high_resolution_clock::now();
		
		const double t = dt*k + t0;
		
		// create a ScalarField in Vh from a list of coefficients 
		// at this point the ScalarField hasn't been discretized yet. 
		// it must be discretized before being printed, integrated or used in an expression
		FiniteElementScalarField u(&Vh, reinterpret_cast< Scalar* >(uk.data()));
		if (k % 100 == 0) // print the function every 100 iterations
		{
			printFunction("std_sem_u_" + std::to_string(idx) + ".dat", u.discretize());
			++idx;
		}
		#pragma omp parallel for
		for (size_t i=0;i<Vh.getNBasisFunction();++i)
		{
			ukp12[i] = uk[i] + dt*vk[i] +  + 0.5*dt*dt*ak[i];
		}
		#pragma omp parallel for
		for (size_t i=0;i<Vh.getNBasisFunction();++i)
		{
			const double Kukp12i = K.apply(ukp12, i);
			
			const double akp1i = invSemM[i]*(ricker(t)*f[i] - Kukp12i);
			const double vkp1i = vk[i] + (1.0 - gamma)*dt*ak[i] + gamma*dt*akp1i;
			const double ukp1i = uk[i] + dt*vk[i] + 0.5*dt*dt*ak[i];
			
			uk[i] = ukp1i;
			vk[i] = vkp1i;
			ak[i] = akp1i;
		}
		
		auto stop = std::chrono::high_resolution_clock::now();
		duration[k] = double(std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count())*1.0e-6;
		
		std::cout << "Iteration " << k << "/ " << nt << " computed in " << duration[k] << " seconds" << std::endl;
	}
	
	std::cout << "Average time for a single iteration: " << std::reduce(std::begin(duration), std::end(duration)) / duration.size() << std::endl;

	return EXIT_SUCCESS;
}
