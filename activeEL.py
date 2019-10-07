from firedrake import *
import numpy as np
import matplotlib.pyplot as plt

def ActiveEL(eg=0, nStep=8, degx=2, nex=-1, degt=-1,
       T=1.0, ndim=2,O_eta = 1):
    "Solve Active Stress Erickesn Leslie equations"
    6
    if nex < 1:
        nex = nStep

    if degt < 1:
        degt = degx-1

    Lx = 1
    Ly = 6
    
    
    ney = nex * Ly
    tau = T*(1.0/ 64)
    nStep = int(1.0/tau)
    print(nStep)

    psi  = FiniteElement("CG", interval, degt)
    psi0 = FiniteElement("DG", interval, 0)

    print("Ericksen Leslie Code")

    if ndim == 2:
        #from examples2d import Examples, EgMax, mkWofRel
        from examples2dCopy1 import Examples, EgMax, mkWofRel
        assert (degx > 1)

        phiv = FiniteElement("CG", quadrilateral, degx)
        phip = FiniteElement("CG", quadrilateral, degx-1)
        phid = FiniteElement("CG", quadrilateral, degx)

        m0 = RectangleMesh(nex, ney, Lx,Ly, quadrilateral=True)

        print("Solution on [0,1]*[0,9]")
        print("Mesh Boundary Codes: [left,right,bottom,top] = [1,2,3,4]")
    elif ndim == 3:
        from examples3d import Examples, EgMax, mkWofRel

        assert (degx > 1)

        phiv = FiniteElement("CG", tetrahedron, degx)
        phip = FiniteElement("CG", tetrahedron, degx-1)
        phid = FiniteElement("CG", tetrahedron, degx)

        m0 = UnitCubeMesh(nex, nex, nex)

        print("Solution on [0,1]^{:d}".format(ndim))
        print("Mesh Boundary Codes: [x=0,x=1,y=0,y=1,z=0,z=1] = [1,2,3,4,5,6]")

        bdy0 = [1,3,5]   # Dirichlet data at x=0 and y = 0 and z=0
        bdy1 = [2,4,6]   # Neumann   data at x=1 and y = 1 and z=1
    else:
        print("Dimension ndim = {:d} must be 2 or 3".format(ndim))
        exit(1)

    if (eg < 0) or (eg >= EgMax):
        print("eg = {:d} not in [0,{:d}) -- terminating".format(eg,EgMax))
        exit(1)
    print(eg)
    exact = Examples[eg].exact
    frhs  = Examples[eg].rhs
    rho, bdy0, bdy1, penalty, kk, alpha = Examples[eg].initialize()
         # active 
    eta = Examples[eg].eta
    eta = O_eta 
    plotFile="Output/eta_{:f}_{:f}/eta_{:f}_{:f}_Ael.pvd".format(eta, eg,eta,eg)

    print("Example {:d}, {:s}".format(eg, Examples[eg].desc))

    #rho, bdy0, bdy1 = Examples[eg].initialize()

    print("Polynomial Degree in space degx ..... ..... {:d}".format(degx))
    print("Polynomial Degree in time  degt ........ .. {:d}".format(degt))
    print("Number of time steps nStep ................ {:d}".format(nStep))
    print("Number elements in each direction nex ..... {:d}".format(nex))
    print("Dirichlet boundary codes bdy0 ............. {}".format(bdy0))
    print("Neumann   boundary codes bdy1 ............. {}".format(bdy1))
    print("Final time T .............................. {:f}".format(T))
    print("Time step tau ............................. {:f}".format(tau))
    print("Plot to file plotFile ..................... {:s}".format(plotFile))
    print("")

    set_log_level(ERROR)            # Only print errors, eliminate warnings

    # penalty = 0.1             # Penatly for |d|^2 = 1
    # kk = [1,1,1,0,0]          # Frank + Chiral Constants
    # alpha = [0,0,1,1,0,0]     # Leslie Constants

    Wof, Rel = mkWofRel(penalty, kk, alpha)

    # Space-time mesh

    mm = ExtrudedMesh(m0, 1, tau)
    
    nn = FacetNormal(mm)         # Boundary normals for Neuman BC
    n0 = np.array(nn)[0:ndim]    # nn[0:ndim] not allowed
    
    xxt = SpatialCoordinate(mm)

    tt = xxt[ndim]
    xx = np.array(xxt)[0:ndim]

    VV = VectorFunctionSpace(mm, TensorProductElement(phiv, psi), dim=ndim)
    PP = FunctionSpace(mm, TensorProductElement(phip, psi))
    DD = VectorFunctionSpace(mm, TensorProductElement(phid, psi), dim=ndim)

    VPD = MixedFunctionSpace([VV,PP,DD])

    vpd = Function(VPD)            # vv,pp,dd = Function(VPD) not allowed
    vv,pp,dd = split(vpd)          # vpd.split() compiles but is WRONG!!!
    ww,qq,ee = TestFunctions(VPD)  # TestFunctions() (plural!)

    dvdt = grad(vv)[:, ndim]
    dddt = grad(dd)[:, ndim]
    dedt = grad(ee)[:, ndim]

    rng = range(0,ndim)
    
    dvdx = np.array([[grad( vv )[i,j] for j in rng] for i in rng])
    dwdx = np.array([[grad( ww )[i,j] for j in rng] for i in rng])
    dddx = np.array([[grad( dd )[i,j] for j in rng] for i in rng])
    detx = np.array([[grad(dedt)[i,j] for j in rng] for i in rng])

    ddot = dddt + np.dot(dddx,vv)

    # Solution at the top of previous time step = initial data for next step
    # Should really live on mesh m0; however, the two meshes don't talk

    V0  = VectorFunctionSpace(mm, TensorProductElement(phiv, psi0), dim=ndim)
    P0  =       FunctionSpace(mm, TensorProductElement(phip, psi0))
    D0  = VectorFunctionSpace(mm, TensorProductElement(phid, psi0), dim=ndim)

    VPD0 = MixedFunctionSpace([V0,P0,D0])
    v0,p0,d0 = split(TrialFunction(VPD0))
    w0,q0,e0 = TestFunctions(VPD0)          # TestFunctions() (plural!)

    a0 = (np.dot(v0,w0) + p0*q0 + np.dot(d0,e0)) * ds_t

    vpd0 = Function(VPD0, name=['Velocity','Pressure','Director'])
    v0,p0,d0 = vpd0.split()

    dvdx0 = np.array([[grad( v0 )[i,j] for j in rng] for i in rng]) 
    dddx0 = np.array([[grad( d0 )[i,j] for j in rng] for i in rng])
    
    


    fv = Function(VV)      # right hand side
    fd = Function(DD)      # right hand side
    wof, dWdd, dWddx = Wof(dd, dddx)
    rel, dRddt, dRdvx = Rel(dd, ddot, dvdx)

    


    wp = (np.dot(rho*dvdt + (rho/2)*np.dot(dvdx,vv)
                 + np.dot(np.transpose(dddx), dRddt-fd) - fv, ww)
          - (rho/2) * np.dot(np.dot(dwdx,vv), vv) 
          + np.tensordot(dRdvx+eta*np.outer(dd,dd), dwdx) 
          - pp * np.trace(dwdx)   
          + np.trace(dvdx) * qq 
          + np.dot(dRddt + dWdd - fd, dedt) 
          + np.tensordot(dWddx, detx)) * dx \
    + (np.dot(vv-v0,ww) + np.dot(dd-d0,ee)) *ds_b




    sparms = {'mat_type': 'aij',
              'snes_type': 'newtonls',
              'ksp_type': 'preonly',
              'pc_type': 'lu',
              'pc_factor_mat_solver_type': 'mumps',
              'snes_max_it': 50,
              'ksp_atol': 1e-12,
              'ksp_stol': 1e-12,
              'ksp_rtol': 1e-12}

    # Initial Data

    t = 0.0

    vpde, dvpdex, dvpdet = exact(t,xx);
    wofe, dWdde, dWddxe = Wof(vpde[2], dvpdex[2])

    solve(a0 == (np.dot(vpde[0],w0)+(vpde[1]+wofe)*q0+np.dot(vpde[2],e0))*ds_t, vpd0, solver_parameters = sparms)

    output = File(plotFile);
    output.write(v0,p0,d0, time=t)
    print("\nTime T = {:f}".format(t))

    print("L2  norms at T: |v| = {:.6e}, |vh| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(np.dot(vpde[0],vpde[0]) * ds_t)),
        sqrt(assemble(np.dot(v0,v0) * ds_t)),
        sqrt(assemble(np.dot(vpde[0]-v0, vpde[0]-v0) * ds_t))))

    print("H1 s-norm at T: |v| = {:.6e}, |vh| = {:.6e}, error = {:.6e}".format(
        sqrt(assemble(np.tensordot(dvpdex[0],dvpdex[0]) * ds_t)),
        sqrt(assemble(np.tensordot(dvdx0,dvdx0) * ds_t)),
        sqrt(assemble(np.tensordot(dvpdex[0]-dvdx0, dvpdex[0]-dvdx0)* ds_t))))

    phate = vpde[1] + wofe
    
    print("L2  norms at T: |p| = {:.6e}, |ph| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(phate*phate * ds_t)),
        sqrt(assemble(p0*p0 * ds_t)),
        sqrt(assemble((p0-phate)*(p0-phate) * ds_t))))

    print("L2  norms at T: |d| = {:.6e}, |dh| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(np.dot(vpde[2],vpde[2]) * ds_t)),
        sqrt(assemble(np.dot(d0,d0) * ds_t)),
        sqrt(assemble(np.dot(vpde[2]-d0, vpde[2]-d0) * ds_t))))

    print("H1 s-norm at T: |d| = {:.6e}, |dh| = {:.6e}, error = {:.6e}".format(
        sqrt(assemble(np.tensordot(dvpdex[2],dvpdex[2]) * ds_t)),
        sqrt(assemble(np.tensordot(dddx0,dddx0) * ds_t)),
        sqrt(assemble(np.tensordot(dvpdex[2]-dddx0, dvpdex[2]-dddx0)* ds_t))))

    print("")

    error_array = np.zeros(nStep)
    for iStep in range(0,nStep):
        vpde, dvpdex, dvpdet = exact(tt+t,xx);
        d_old_L2 = sqrt(assemble(np.dot(d0,d0) * ds_t))

        # Note: Can not use VV in place of VPD.sub(0) or DD for VPD.sub(2)
        # or [DirichletBC(VPD, Expression(vpde), i) for i in bdy0]
        
        bcV = [DirichletBC(VPD.sub(0), vpde[0], i) for i in bdy0]
        bcD = [DirichletBC(VPD.sub(2), vpde[2], i) for i in bdy0]
        
        fve, fde = frhs(tt+t,xx)
        
        fv.interpolate(as_vector(fve))
        fd.interpolate(as_vector(fde))

        ddote = dvpdet[2] + np.dot(dvpdex[2], vpde[0])

        wofe, dWdde, dWddxe = Wof(vpde[2], dvpdex[2])
        rel, dRddte, dRdvxe = Rel(vpde[2], ddote, dvpdex[0])

        gg = [(
            np.dot(np.dot(dRdvxe+eta*np.outer(vpde[2],vpde[2]),n0)-(vpde[1]+wofe)*n0
                  -(rho/2)*np.dot(vpde[0],n0)*vpde[0], ww) 
            +np.dot(np.dot(dWddxe,n0), dedt))
             * ds_v(i) for i in bdy1]
        gg =  []# Do nothing bc

        eqn = wp - sum(gg)
        
        solve(eqn == 0, vpd, bcs=bcV+bcD, solver_parameters=sparms)

        solve(a0 == (np.dot(vv,w0)+pp*q0+np.dot(dd,e0))*ds_t, vpd0)

        t += tau
        

        output.write(v0, p0, d0, time=t)
        print("\nL2  norms at T = {:4f}".format(t), "|dh| = {:.4e}, error = {:.4e}".format(\
        sqrt(assemble(np.dot(d0,d0) * ds_t)),
        np.absolute(d_old_L2 - sqrt(assemble(np.dot(d0,d0) * ds_t)))))
        error_array[iStep] = np.absolute(d_old_L2 - sqrt(assemble(np.dot(d0,d0) * ds_t)))
        if error_array[iStep] > 10e-2:
            break
        
        

        # END for iStep in range(0,nStep):

    print("\nTime T = {:f}".format(t))
    
    

    vpde, dvpdex, dvpdet = exact(tt+t,xx);

    print("L2  norms at T: |v| = {:.6e}, |vh| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(np.dot(vpde[0],vpde[0]) * ds_t)),
        sqrt(assemble(np.dot(v0,v0) * ds_t)),
        sqrt(assemble(np.dot(vpde[0]-v0, vpde[0]-v0) * ds_t))))
 


    print("H1 s-norm at T: |v| = {:.6e}, |vh| = {:.6e}, error = {:.6e}".format(
        sqrt(assemble(np.tensordot(dvpdex[0],dvpdex[0]) * ds_t)),
        sqrt(assemble(np.tensordot(dvdx0,dvdx0) * ds_t)),
        sqrt(assemble(np.tensordot(dvpdex[0]-dvdx0, dvpdex[0]-dvdx0)* ds_t))))

    wofe, dWdde, dWddxe = Wof(vpde[2], dvpdex[2])
    phate = vpde[1] + wofe
    
    print("L2  norms at T: |p| = {:.6e}, |ph| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(phate*phate * ds_t)),
        sqrt(assemble(p0*p0 * ds_t)),
        sqrt(assemble((p0-phate)*(p0-phate) * ds_t))))

    print("L2  norms at T: |d| = {:.6e}, |dh| = {:.6e}, error = {:.6e}".format(\
        sqrt(assemble(np.dot(vpde[2],vpde[2]) * ds_t)),
        sqrt(assemble(np.dot(d0,d0) * ds_t)),
        sqrt(assemble(np.dot(vpde[2]-d0, vpde[2]-d0) * ds_t))))

    print("H1 s-norm at T: |d| = {:.6e}, |dh| = {:.6e}, error = {:.6e}".format(
        sqrt(assemble(np.tensordot(dvpdex[2],dvpdex[2]) * ds_t)),
        sqrt(assemble(np.tensordot(dddx0,dddx0) * ds_t)),
        sqrt(assemble(np.tensordot(dvpdex[2]-dddx0, dvpdex[2]-dddx0)* ds_t))))

    print("")
    


    return error_array
