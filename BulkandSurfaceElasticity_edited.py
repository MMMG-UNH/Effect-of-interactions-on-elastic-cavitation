"""
Code for bulk and surface elasticity models in FEniCS.

This module defines the Bulk_SurfaceElasticity class, which sets up the finite element
formulation for an incompressible neo-Hookean solid with optional surface energy at the
cavity boundary. It handles mesh loading, function space construction, kinematics,
constitutive relations, and the weak form assembly for use in axisymmetric simulations.

The writeResults() function at the bottom handles exporting displacement, pressure,
and stress fields to XDMF files for post-processing in ParaView.
"""

# Fenics-related packages
from dolfin import *

class Bulk_SurfaceElasticity:
    def __init__(self, elastic_incompressible, bulk_elastic_model,
                 surface_elastic_model, bulkparams, surfaceparams):
        """
        Store model type flags and material parameters.

        Parameters
        ----------
        elastic_incompressible : int
            Set to 1 for incompressible material (uses mixed u-p formulation)
        bulk_elastic_model : int
            Set to 1 for neo-Hookean bulk model (currently the only supported option)
        surface_elastic_model : int
            Set to 1 for isotropic surface tension model (currently the only supported option)
        bulkparams : list
            [Gshear, Kbulk] for incompressible neo-Hookean:
                Gshear — shear modulus (sets the stress scale; dimensionless results use Gshear=1)
                Kbulk  — penalty bulk modulus (set large, e.g. 1e12*Gshear, to enforce incompressibility)
        surfaceparams : list
            [gambar] — dimensionless surface energy parameter (set to 0 to neglect surface tension)
        """
        self.elastic_incompressible = elastic_incompressible
        self.bulk_elastic_model = bulk_elastic_model
        self.surface_elastic_model = surface_elastic_model

        if(bulk_elastic_model==1):
            if(elastic_incompressible==1):
                assert len(bulkparams) == 2, "Incorrect number of parameters for bulk elastic model"
                self.Gshear  = bulkparams[0]
                # Kbulk acts as a penalty parameter to enforce near-incompressibility
                self.Kbulk   = bulkparams[1]
            else:
                assert len(bulkparams) == 2, "Incorrect number of parameters for bulk elastic model"
                self.Gshear  = bulkparams[0]
                # lambdaparam is the first Lame parameter for compressible formulation
                self.lambdaparam = bulkparams[1]
        else:
            assert False, "Invalid choice for bulk elastic model"

        if(surface_elastic_model==1):
            assert len(surfaceparams) == 1, "Incorrect number of parameters for surface elastic model"
            # gambar: dimensionless surface energy coefficient at the cavity surface
            self.gambar = surfaceparams[0]
        else:
            assert False, "Invalid choice for surface elastic model"


    def setup(self, meshname_string, simulation_type, surface_energy_boundaries):
        """
        Load the mesh, set up integration measures, and build the weak form.

        Parameters
        ----------
        meshname_string : str
            Base name of the mesh files (without extension), e.g. '2Cavity_yc1.10_finest'
        simulation_type : int
            1 = full 3D, 2 = axisymmetric, 3 = plane strain
        surface_energy_boundaries : int
            Physical group tag of the cavity surface boundary (where surface energy is applied)
        """
        self.mesh = Mesh()
        # mesh files are expected in the 'meshes/' subdirectory
        meshname = "meshes/" + meshname_string + ".xdmf"
        facet_meshname = "meshes/facet_" + meshname_string + ".xdmf"
        self.simulation_type = simulation_type
        self.surface_energy_boundaries = surface_energy_boundaries

        # read the triangular domain mesh from the XDMF file
        with XDMFFile(meshname) as infile:
            infile.read(self.mesh)

        # store reference (undeformed) spatial coordinates for use in axisymmetric integrals
        self.x = SpatialCoordinate(self.mesh)

        # read mesh boundary facets (line elements with physical group tags)
        if (simulation_type == 1):
            mvc = MeshValueCollection("size_t", self.mesh, 2)
        else:
            mvc = MeshValueCollection("size_t", self.mesh, 1)
        # read the facet physical group tags from the facet XDMF file
        with XDMFFile(facet_meshname) as infile:
            infile.read(mvc, "name_to_read")

        # store facet tags as a MeshFunction for use in boundary condition and measure definitions
        self.facets = cpp.mesh.MeshFunctionSizet(self.mesh, mvc)

        # ds: boundary integration measure, tagged by physical group IDs from the .geo file
        self.ds = Measure('ds', domain=self.mesh, subdomain_data=self.facets)
        # dsm and dxm: integration measures with explicit quadrature settings for accuracy
        self.dsm= self.ds(metadata={"quadrature_degree": 2, "quadrature_scheme": "default"})
        self.dxm = dx(metadata={"quadrature_degree": 2, "quadrature_scheme": "default"})

        # facet normal vector in the reference configuration
        if(simulation_type==1):
            self.Facet_normal = FacetNormal(self.ds)
        else:
            self.Facet_normal=FacetNormal(self.mesh)

        self.setUpFunctionSpaces()
        self.setUpWeakform()

    def setUpFunctionSpaces(self):
        """
        Define the mixed function space (displacement + pressure) and trial/test functions.

        Uses the Taylor-Hood element pair: quadratic (P2) for displacement,
        linear (P1) for pressure — a stable inf-sup satisfying combination for
        incompressible elasticity.
        """
        # P2 vector element for displacement
        self.U2 = VectorElement("Lagrange", self.mesh.ufl_cell(), 2)
        # P1 scalar element for pressure (Lagrange multiplier for incompressibility)
        self.P1 = FiniteElement("Lagrange", self.mesh.ufl_cell(), 1)

        if(self.elastic_incompressible==1):
            # mixed Taylor-Hood element: (u, p) coupled in a single function space
            TH = MixedElement([self.U2, self.P1])
            self.ME = FunctionSpace(self.mesh, TH)
        else:
            # displacement-only space for compressible formulation
            self.ME = FunctionSpace(self.mesh, self.U2)

        # test function for the weak form residual
        w_test = TestFunction(self.ME)
        if(self.elastic_incompressible==1):
            # split into displacement and pressure test functions
            self.u_test, self.p_test = split(w_test)
        else:
            self.u_test = w_test

        # trial function used by FEniCS automatic differentiation to compute the Jacobian
        self.dw = TrialFunction(self.ME)

        # the actual solution function (holds current displacement and pressure DOFs)
        self.w = Function(self.ME)
        if(self.elastic_incompressible==1):
            # split into displacement field u and pressure field p
            self.u, self.p = split(self.w)
        else:
            self.u = self.w

        # copy of the solution from the previous load step (used for incremental updates)
        self.w_old = Function(self.ME)
        if(self.elastic_incompressible==1):
            u_old, p_old = split(self.w_old)
        else:
            u_old = self.w_old

    '''''''''''''''''''''
    Some subroutines
    '''''''''''''''''''''
    # VECTOR FIELD GRADIENT

    def axi_grad_vector(self, u):
        """
        Compute the 3x3 displacement gradient tensor for axisymmetric problems (r, z, theta).

        In axisymmetry, the out-of-plane (theta) component of the gradient is u_r / r.
        A conditional is used to avoid division by zero at the symmetry axis (r=0).
        """
        grad_u = grad(u)
        grad_u33_exp =  u[0]/self.x[0]
        # at r=0 the hoop strain u_r/r is replaced by its limit (0) to avoid singularity
        grad_u33 = conditional(lt(self.x[0], DOLFIN_EPS), 0.0, grad_u33_exp)
        return as_tensor([[grad_u[0,0], grad_u[0,1], 0],
                    [grad_u[1,0], grad_u[1,1], 0],
                    [0, 0, grad_u33]])

    def pe_grad_vector(self, u):
        """
        Compute the 3x3 displacement gradient tensor for plane strain problems.

        The out-of-plane strain component is zero by assumption.
        """
        grad_u = grad(u)
        return as_tensor([[grad_u[0,0], grad_u[0,1], 0],
                    [grad_u[1,0], grad_u[1,1], 0],
                    [0, 0, 0]])

    def grad_vector(self, u):
        """
        Dispatch to the correct gradient routine based on simulation type.
        """
        if(self.simulation_type==1): #Full 3D
            return grad(u)
        elif(self.simulation_type==2):
            return self.axi_grad_vector(u)
        elif(self.simulation_type==3):
            return self.pe_grad_vector(u)

    # SCALAR FIELD GRADIENT

    def axi_grad_scalar(y):
        """Gradient of a scalar field in axisymmetric coordinates."""
        grad_y = grad(y)
        return as_vector([grad_y[0], grad_y[1], 0.])

    def pe_grad_scalar(y):
        """Gradient of a scalar field in plane strain coordinates."""
        grad_y = grad(y)
        return as_vector([grad_y[0], grad_y[1], 0.])

    def grad_scalar(self, y):
        """
        Dispatch to the correct scalar gradient routine based on simulation type.
        """
        if(self.simulation_type==1): #Full 3D
            return grad(y)
        elif(self.simulation_type==2):
            return self.axi_grad_scalar(y)
        elif(self.simulation_type==3):
            return self.pe_grad_scalar(y)

    def F_3D_calc(self, u):
        """
        Compute the 3x3 deformation gradient tensor F = I + Grad(u).
        """
        F  = Identity(3) + self.grad_vector(u)
        return F

    '''''''''''''''''''''''''''''
    Constitutive equations
    '''''''''''''''''''''''''''''
    def setUpWeakform(self):
        """
        Define kinematics, constitutive relations, and the weak form residual.

        Constitutive model: incompressible neo-Hookean bulk + isotropic surface tension.
        The incompressibility constraint J=1 is enforced weakly via a perturbed Lagrangian
        approach: the constraint (J-1) - p/Kbulk = 0 is added to the weak form, where
        p is the Lagrange multiplier (hydrostatic pressure) and Kbulk is a large penalty
        parameter interpreted as a bulk modulus.
        """
        # deformation gradient F and its determinant J (= 1 for incompressible material)
        self.F_3D    = variable(self.F_3D_calc(self.u))
        self.J = det(self.F_3D)
        F_inv_T = variable((inv(self.F_3D.T)))
        # I1 = tr(F^T F): first invariant of the right Cauchy-Green tensor
        I1 = variable(inner(self.F_3D, self.F_3D))

        # C = F^T F: right Cauchy-Green deformation tensor
        C_3D    = variable(self.F_3D.T*self.F_3D)
        # distortional part of C (volume-preserving): removes volumetric contribution
        Cdis_3D = variable(det(C_3D)**(-1/3)*C_3D)
        # I1_dis: first invariant of the distortional Cauchy-Green tensor
        I1_dis = variable(tr(Cdis_3D))

        # reference outward normal on boundary facets (extended to 3D for axisymmetric case)
        if(self.simulation_type==1):
            self.nR_3D = self.Facet_normal
        else:
            self.nR_3D = as_tensor([self.Facet_normal[0], self.Facet_normal[1], 0])

        # deformed normal vector (Nanson's formula): maps reference normal to deformed configuration
        self.deformed_normal = variable(dot(self.J*F_inv_T, self.nR_3D))
        # areal stretch: ratio of deformed to reference surface area element
        self.ArealStretch = variable(sqrt(dot(self.deformed_normal, self.deformed_normal)))
        # unit normal in the deformed configuration
        self.deformed_unit_normal = self.deformed_normal/self.ArealStretch

        # --- Free energy densities ---
        if(self.bulk_elastic_model==1):
            if(self.elastic_incompressible==1):
                # neo-Hookean strain energy density (distortional part only)
                self.psi_bulk = self.Gshear/2*(I1_dis-3)
                # perturbed Lagrangian: adds pressure work term and penalty for volumetric deviation
                # the term -(1/2/Kbulk)*p^2 is the compliance correction for the penalty approach
                self.psi = variable(self.psi_bulk + self.p*(self.J-1) -(1/2.0/self.Kbulk)*self.p**2.0)
            else:
                # compressible neo-Hookean: standard form with logarithmic volumetric term
                self.psi_bulk = (self.Gshear/2)*(I1-3-2*ln(self.J)) + (self.lambdaparam/2)*(ln(self.J)*ln(self.J))
                self.psi = variable(self.psi_bulk)

        if(self.surface_elastic_model==1):
            # surface energy density: proportional to areal stretch (constant surface tension gambar)
            self.psi_gamma = self.gambar*self.ArealStretch

        # first Piola-Kirchhoff stress via automatic differentiation of the free energy w.r.t. F
        Piola = diff(self.psi, self.F_3D)
        # Cauchy stress: T = (1/J) * P * F^T
        self.T    = Piola*self.F_3D.T/self.J

        # surface Piola stress via automatic differentiation of surface energy w.r.t. F
        Surface_S = diff(self.psi_gamma, self.F_3D)

        # --- Weak form residual ---
        # for axisymmetric problems, integrals include the factor 2*pi*r (here just x[0] since
        # 2*pi cancels in the dimensionless formulation); the surface term is restricted to
        # the cavity boundary only
        if(self.simulation_type==2):
            self.L = inner(Piola, self.grad_vector(self.u_test))*self.x[0]*self.dxm \
                   + inner(Surface_S, self.grad_vector(self.u_test))*self.x[0]*self.dsm(self.surface_energy_boundaries)
            if(self.elastic_incompressible==1):
                # weak form of the perturbed incompressibility constraint
                self.L = self.L + inner(((self.J-1.0) - self.p/self.Kbulk), self.p_test)*self.x[0]*self.dxm
        else:
            self.L = inner(Piola, self.grad_vector(self.u_test))*self.dxm \
                   + inner(Surface_S, self.grad_vector(self.u_test))*self.dsm(self.surface_energy_boundaries)
            if(self.elastic_incompressible==1):
                self.L = self.L + inner(((self.J-1.0) - self.p/self.Kbulk), self.p_test)*self.dxm

        # tangent stiffness matrix (Jacobian of the residual): computed by automatic differentiation
        self.a = derivative(self.L, self.w, self.dw)


def writeResults(BVP, file_results, t):
    """
    Write displacement, pressure, and stress fields to the XDMF results file at time t.

    Parameters
    ----------
    BVP : Bulk_SurfaceElasticity
        The boundary value problem object containing solution fields
    file_results : XDMFFile
        The open XDMF file to write results into
    t : float
        Current pseudo-time (load step) used as the time stamp in the output file
    """
    # split the mixed solution into displacement and pressure components
    u_Vis, p_Vis = BVP.w.split()
    u_Vis.rename("displacement", " ")
    file_results.write(u_Vis, t)

    if BVP.elastic_incompressible == 1:
        p_Vis.rename("p", " ")
        file_results.write(p_Vis, t)

    # scalar function space for projecting derived quantities (J, stress components)
    W = FunctionSpace(BVP.mesh, BVP.P1)

    # J = det(F): volume ratio (should remain close to 1 for incompressible material)
    J_Vis = project(BVP.J, W)
    J_Vis.rename("J", " ")
    file_results.write(J_Vis, t)

    # project all Cauchy stress tensor components onto the scalar space for visualization
    T11_Vis = project(BVP.T[0,0], W);  T11_Vis.rename("T11", " ")
    T22_Vis = project(BVP.T[1,1], W);  T22_Vis.rename("T22", " ")
    T33_Vis = project(BVP.T[2,2], W);  T33_Vis.rename("T33", " ")

    T12_Vis = project(BVP.T[0,1], W);  T12_Vis.rename("T12", " ")
    T13_Vis = project(BVP.T[0,2], W);  T13_Vis.rename("T13", " ")
    T23_Vis = project(BVP.T[1,2], W);  T23_Vis.rename("T23", " ")

    T21_Vis = project(BVP.T[1,0], W);  T21_Vis.rename("T21", " ")
    T31_Vis = project(BVP.T[2,0], W);  T31_Vis.rename("T31", " ")
    T32_Vis = project(BVP.T[2,1], W);  T32_Vis.rename("T32", " ")

    # trace of Cauchy stress: tr(T) = T11 + T22 + T33 (= 3 * mean normal stress)
    tr_T_Vis   = project(tr(BVP.T),     W);  tr_T_Vis.rename("tr_T",  " ")
    # effective (Frobenius norm) stress: scalar measure of overall stress magnitude
    T_eff_Vis  = project(sqrt(inner(BVP.T, BVP.T)), W)
    T_eff_Vis.rename("Teff", " ")

    # write all stress fields to the XDMF file
    file_results.write(T11_Vis,  t)
    file_results.write(T22_Vis,  t)
    file_results.write(T33_Vis,  t)
    file_results.write(T12_Vis,  t)
    file_results.write(T13_Vis,  t)
    file_results.write(T23_Vis,  t)
    file_results.write(T21_Vis,  t)
    file_results.write(T31_Vis,  t)
    file_results.write(T32_Vis,  t)
    file_results.write(tr_T_Vis, t)
    file_results.write(T_eff_Vis, t)