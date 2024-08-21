# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Structural FEA on a L shaped domain fixed at top 
with adaptive mesh refinement using a posteriori error estimator 

@author: bhttchr6

"""
import numpy as np
import scipy as S
import scipy.linalg as la
import scipy.sparse as sparse
import scipy.sparse.linalg as sla

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import meshpy.triangle as triangle

plt.close('all')

def f(xvec):
    x, y = xvec
    if (x==1) & (y==0.5):
        return 0.5
    else:
        return 0.0

class MatrixBuilder:
    def __init__(self):
        self.rows = []
        self.cols = []
        self.vals = []
        
    def add(self, rows, cols, submat):
        for i, ri in enumerate(rows):
            for j, cj in enumerate(cols):
                self.rows.append(ri)
                self.cols.append(cj)
                self.vals.append(submat[i, j])
                
    def coo_matrix(self):
        return sparse.coo_matrix((self.vals, (self.rows, self.cols)))
 
def coo_submatrix_pull(matr, rows, cols):
    """
    Pulls out an arbitrary i.e. non-contiguous submatrix out of
    a sparse.coo_matrix. 
    """
    if type(matr) != S.sparse.coo_matrix:
        raise TypeError('Matrix must be sparse COOrdinate format')
    
    gr = -1 * np.ones(matr.shape[0])
    gc = -1 * np.ones(matr.shape[1])
    
    lr = len(rows)
    lc = len(cols)
    
    ar = np.arange(0, lr)
    ac = np.arange(0, lc)
    gr[rows[ar]] = ar
    gc[cols[ac]] = ac
    mrow = matr.row
    mcol = matr.col
    newelem = (gr[mrow] > -1) & (gc[mcol] > -1)
    newrows = mrow[newelem]
    newcols = mcol[newelem]
    return S.sparse.coo_matrix((matr.data[newelem], np.array([gr[newrows],
        gc[newcols]])),(lr, lc))    
    
def round_trip_connect(start, end):
    return [(i, i+1) for i in range(start, end)] + [(end, start)]

def make_mesh(flag, marked_elem,itr):
    points = [(0, 0), (1, 0), (1, 0.5), (0.5, 0.5), (0.5,1), (0,1)]
    #points = [(0, 0), (1, 0), (1,1), (0,1)]
    facets = round_trip_connect(0, len(points)-1)

    '''circ_start = len(points)
    points.extend(
            (0.25 * np.cos(angle), 0.25 * np.sin(angle))
            for angle in np.linspace(0, 2*np.pi, 30, endpoint=False))

    facets.extend(round_trip_connect(circ_start, len(points)-1))'''

    def needs_refinement(vertices, area):
        bary = np.sum(np.array(vertices), axis=0)/3
        max_area = 0.01 + la.norm(bary, np.inf)*0.001
        #max_area = 1 + la.norm(bary, np.inf)*1
        return bool(area > max_area)

    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_facets(facets)

    built_mesh = triangle.build(info, refinement_func=needs_refinement)
    if flag==1:
        built_mesh.element_volumes.setup()
        for i in range(len(built_mesh.elements)):
            built_mesh.element_volumes[i] = -1
        for i in range(0, len(built_mesh.elements)):
            if i in marked_elem:
             built_mesh.element_volumes[i] = 0.001
    
        built_mesh = triangle.refine(built_mesh)
        
    return np.array(built_mesh.points), np.array(built_mesh.elements)

def max_edge_len(x0,x1,x2):
    
    eps = 0.0001
    
    
    
    v_1_x = x1[0]-x0[0]+eps
    v_1_y = x1[1]-x0[1]+eps
    e_1 = (v_1_x**2+v_1_y**2)**(1/2)
    b1 = (1/(1+(v_1_y**2/v_1_x**2)))**(1/2)
    a1 = -b1*(v_1_y/v_1_x)
    
    v_1 = np.array([a1,b1])
    
    # edge 2
    v_2_x = x2[0]-x1[0]+eps
    v_2_y = x2[1]-x1[1]+eps
    e_2 = (v_2_x**2+v_2_y**2)**(1/2)
    b2 = (1/(1+(v_2_y**2/v_2_x**2)))**(1/2)
    a2 = -b2*(v_2_y/v_2_x)
    v_2 = np.array([a2,b2])
    
    # edge 1
    v_3_x = x2[0]-x0[0]+eps
    v_3_y = x2[1]-x0[1]+eps
    e_3 = (v_3_x**2+v_3_y**2)**(1/2)
    b3 = (1/(1+(v_3_y**2/v_3_x**2)))**(1/2)
    a3 = -b3*(v_3_y/v_3_x)
    v_3 = np.array([a3,b3])
        
    h = np.amax([e_1, e_2, e_3])

    return h,v_1,v_2,v_3
# Form Stiffness matrix and Internal stress vectors
def FEM_Ktan_Fint(V,E, x_load, y_load, U,load_dir):
    ne = len(E)
    nv = len(V)
    dbasis = np.array([
        [-1, 1, 0],
        [-1, 0, 1]])
    
    
    
    dN_ds = np.array([
        [-1, 0, 1, 0, 0, 0],
        [0, -1, 0, 0, 0, 1]])
    
    YM = 10
    mu = 0.3
    
    cond=1
    if cond==1:  # plane stress
        sc = YM/(1-mu**2)
        C = np.array([
            [sc, mu*sc, 0],
            [mu*sc, sc, 0],
            [0, 0, (1-mu)*sc]])
    
    
    if cond == 0: # plane strain
        sc = YM/((1-2*mu)*(1+mu))
        C = np.array([
            [(1-mu)*sc, mu*sc, 0],
            [mu*sc, (1-mu)*sc, 0],
            [0, 0, ((1-2*mu)/2)*sc]])
    
    
    a_builder = MatrixBuilder()
    Fint = np.zeros(2*nv)
    ext_dof = []
    
    for ei in range(0, ne):
        vert_indices = E[ei, :]
        
        
        x0, x1, x2 = el_verts = V[vert_indices]
        
        x_dofs = vert_indices*2
        y_dofs = vert_indices*2+1
        #print(x_dofs)
        
        if (x0[0]==x_load) & (x0[1]==y_load):
            ext_dof = vert_indices[0]*2+load_dir
        if (x1[0]==x_load) & (x1[1]==y_load):
            ext_dof = vert_indices[1]*2+load_dir    
        if (x2[0]==x_load) & (x2[1]==y_load):
            ext_dof = vert_indices[2]*2+load_dir    
        el_indices = np.array([
            x_dofs[0],
            y_dofs[0],
            x_dofs[1],
            y_dofs[1],
            x_dofs[2],
            y_dofs[2]])
        #print(el_indices)
        el_disp = np.array([
            U[x_dofs[0]],
            U[y_dofs[0]],
            U[x_dofs[1]],
            U[y_dofs[1]],
            U[x_dofs[2]],
            U[y_dofs[2]]
            ])
        centroid = np.mean(el_verts, axis=0)
    
        J = np.array([x1-x0, x2-x0])
        
        invJT = la.inv(J).T
        detJ = la.det(J.T)
        dN_dx = invJT @ dbasis
        #print(dN_dx)
        B = np.array([[dN_dx[0,0], 0 , dN_dx[0,1], 0, dN_dx[0,2], 0],
                     [0, dN_dx[1,0], 0, dN_dx[1,1], 0, dN_dx[1,2]],
                      [dN_dx[1,0],dN_dx[0,0],dN_dx[1,1],dN_dx[0,1],dN_dx[1,2],dN_dx[0,2]]])
        
        #print(B)
        eps= B@el_disp
        sigma = C@eps
        fint = (detJ / 2.0)*B.T@sigma
        Aelem = (detJ / 2.0) * B.T @C@ B
    
        #print(detJ)
        gl_idx = np.array([
            x_dofs[1],
            y_dofs[1],
            x_dofs[2],
            y_dofs[2],
            x_dofs[0],
            y_dofs[0]])
        #print(gl_idx)
        a_builder.add(el_indices, el_indices, Aelem)
        for i, vi in enumerate(el_indices):
            Fint[vi] += fint[i]
        
    Ktan = a_builder.coo_matrix().tocsr().tocoo()
    return ext_dof,Ktan, Fint

def int_stress(E,V,U):
    ne = len(E)
    nv = len(V)
    dbasis = np.array([
        [-1, 1, 0],
        [-1, 0, 1]])
    
    
    
    dN_ds = np.array([
        [-1, 0, 1, 0, 0, 0],
        [0, -1, 0, 0, 0, 1]])
    
    YM = 10
    mu = 0.3
    
    cond=1
    if cond==1:  # plane stress
        sc = YM/(1-mu**2)
        C = np.array([
            [sc, mu*sc, 0],
            [mu*sc, sc, 0],
            [0, 0, (1-mu)*sc]])
    
    
    if cond == 0: # plane strain
        sc = YM/((1-2*mu)*(1+mu))
        C = np.array([
            [(1-mu)*sc, mu*sc, 0],
            [mu*sc, (1-mu)*sc, 0],
            [0, 0, ((1-2*mu)/2)*sc]])
    
    
    a_builder = MatrixBuilder()
    Fint = np.zeros(2*nv)
    ext_dof = []
    
    for ei in range(0, ne):
        vert_indices = E[ei, :]
        
        
        x0, x1, x2 = el_verts = V[vert_indices]
        
        x_dofs = vert_indices*2
        y_dofs = vert_indices*2+1
        #print(x_dofs)
        
            
        el_indices = np.array([
            x_dofs[0],
            y_dofs[0],
            x_dofs[1],
            y_dofs[1],
            x_dofs[2],
            y_dofs[2]])
        #print(el_indices)
        el_disp = np.array([
            U[x_dofs[0]],
            U[y_dofs[0]],
            U[x_dofs[1]],
            U[y_dofs[1]],
            U[x_dofs[2]],
            U[y_dofs[2]]
            ])
        centroid = np.mean(el_verts, axis=0)
    
        J = np.array([x1-x0, x2-x0])
        
        invJT = la.inv(J).T
        detJ = la.det(J.T)
        dN_dx = invJT @ dbasis
        #print(dN_dx)
        B = np.array([[dN_dx[0,0], 0 , dN_dx[0,1], 0, dN_dx[0,2], 0],
                     [0, dN_dx[1,0], 0, dN_dx[1,1], 0, dN_dx[1,2]],
                      [dN_dx[1,0],dN_dx[0,0],dN_dx[1,1],dN_dx[0,1],dN_dx[1,2],dN_dx[0,2]]])
        
        #print(B)
        eps= B@el_disp
        sigma = C@eps
        fint = (detJ / 2.0)*B.T@sigma
        Aelem = (detJ / 2.0) * B.T @C@ B
    
        #print(detJ)
        gl_idx = np.array([
            x_dofs[1],
            y_dofs[1],
            x_dofs[2],
            y_dofs[2],
            x_dofs[0],
            y_dofs[0]])
        #print(gl_idx)
        
        for i, vi in enumerate(el_indices):
            Fint[vi] += fint[i]
        norm_fint = la.norm(Fint, np.inf)    
    return Fint, norm_fint


def body_force (E, V):
    ne = len(E)
    nv = len(V)
    Fb = np.zeros(2*nv)
    rho = 0.1
    g = 10
    for ei in range(0, ne):
        vert_indices = E[ei, :]
        x0, x1, x2 = el_verts = V[vert_indices]
        centroid = np.mean(el_verts, axis=0)
        x_dofs = vert_indices*2
        y_dofs = vert_indices*2+1
        
        J = np.array([x1-x0, x2-x0]).T
        detJ = la.det(J)
    
        el_indices = np.array([
            x_dofs[0],
            y_dofs[0],
            x_dofs[1],
            y_dofs[1],
            x_dofs[2],
            y_dofs[2]])
    
        belem = (detJ / 6.0) * np.array([
            0,
            rho*g,
            0,
            rho*g,
            0,
            rho*g
            ])
    
        for i, vi in enumerate(el_indices):
            Fb[vi] += belem[i] 
    return Fb  

def error_estimator(V,E, u):
    X, Y = V[:, 0], V[:, 1]
    nv = len(V)
    ne = len(E)
    dbasis = np.array([
    [-1, 1, 0],
    [-1, 0, 1]])
    rho = 0.1
    g = 10
    YM = 10
    mu = 0.3
    K = YM/(1-mu)
    cond=1
    if cond==1:  # plane stress
        sc = YM/(1-mu**2)
        C = np.array([
            [sc, mu*sc, 0],
            [mu*sc, sc, 0],
            [0, 0, (1-mu)*sc]])
    
    
    if cond == 0: # plane strain
        sc = YM/((1-2*mu)*(1+mu))
        C = np.array([
            [(1-mu)*sc, mu*sc, 0],
            [mu*sc, (1-mu)*sc, 0],
            [0, 0, ((1-2*mu)/2)*sc]])
   
    
    tol =1e-12
    
    is_boundary = ((np.abs(Y-1) < tol))
    
    eta_K = np.zeros(ne)
    e_rel = np.zeros(ne)
    eta = 0
    
    mark =[]
    ele_size = 0
    for ei in range(0, ne):
            vert_indices = E[ei, :]
            x0, x1, x2 = el_verts = V[vert_indices]
            x_dofs = vert_indices*2
            y_dofs = vert_indices*2+1
            val0,val1,val2 = is_boundary[vert_indices]
            
            
            #u_el = u[vert_indices]
            el_disp = np.array([
            u[x_dofs[0]],
            u[y_dofs[0]],
            u[x_dofs[1]],
            u[y_dofs[1]],
            u[x_dofs[2]],
            u[y_dofs[2]],
            ])
            
            h,v_1,v_2,v_3 =  max_edge_len(x0,x1,x2)
            
            
            
            
            ele_size = ele_size+h
            #error_el = u_ex_el-u_el
            J1 = la.norm(v_1@dbasis,2) 
            J2 = la.norm(v_2@dbasis,2)
            J3 = la.norm(v_3@dbasis,2)
            
            
            Jacob = np.array([x1-x0, x2-x0]).T
            invJT = la.inv(Jacob.T)
            detJ = la.det(Jacob)
            dN_dx = invJT @ dbasis
            #print(dN_dx)
            B = np.array([[dN_dx[0,0], 0 , dN_dx[0,1], 0, dN_dx[0,2], 0],
                         [0, dN_dx[1,0], 0, dN_dx[1,1], 0, dN_dx[1,2]],
                          [dN_dx[1,0],dN_dx[0,0],dN_dx[1,1],dN_dx[0,1],dN_dx[1,2],dN_dx[0,2]]])
            eps = B@el_disp.T
            sigma = C@eps
            centroid = np.mean(el_verts, axis=0)
            belem = (detJ / 6.0) * np.array([
            0,
            rho*g,
            0,
            rho*g,
            0,
            rho*g
            ])
            
            R_val = belem+B.T@sigma
            R = la.norm(R_val,2)
            
            if ((val0 ==True) & (val1 ==True)):
                J1 ==0
            if ((val1 ==True) & (val2 ==True)):
                J2 ==0
            if ((val2 ==True) & (val0 ==True)):
                J3 ==0    
            J =J1**2+J2**2+J3**2
            c1 = h**2/(24*K)
            c2 = h/(24*K)
            eta_K[ei] = (c1*R**2+ c2*(J))
            eta = eta +eta_K[ei]
            if ((val0==False) | (val1==False)|(val2==False)):
                
                e_rel[ei] = eta_K[ei]/la.norm(el_disp,2)
            if (e_rel[ei]>0.1):
                mark.append(ei)
    ele_size=ele_size/ne
    return mark, eta_K,eta, e_rel, ele_size


def FEM_sol(V,E):
    
    nv = len(V)
    ne = len(E)
    
    #print(V.shape)
    #print(E.shape)
    #print(E.max())
    # Extract the XY coordinates
    X, Y = V[:, 0], V[:, 1]    
    
    
    # Specify load location
    x_load = 1
    y_load = 0.5
    load_dir = 1
    
    #x_load = 1
    #y_load = 1
    
    # Specify external load value
    Fext_val = -0.09
    
    
    val=0
    '''
    for ei in range(0, ne):
        vert_indices = E[ei, :]
        
        
        x0, x1, x2 = el_verts = V[vert_indices]
        
        x_dofs = vert_indices*2
        y_dofs = vert_indices*2+1
        
        
        if (x0[0]>=0.8) & (x0[1]==0.5):
            #ext_dofs[ctr] = vert_indices[0]*2+1
            val=val+1
            #print(vert_indices[0])
        if (x1[0]>=0.8) & (x1[1]==0.5):
            #ext_dofs[ctr] = vert_indices[1]*2+1
            val=val+1
        if (x2[0]>=0.8) & (x2[1]==0.5):
            #ext_dofs[ctr] = vert_indices[2]*2+1 
            val=val+1
    
    #print(val)
    
    ext_dofs = np.zeros((val,),dtype=int)
    ctr=0
    for ei in range(0, ne):
        vert_indices = E[ei, :]
        
        
        x0, x1, x2 = el_verts = V[vert_indices]
        
        x_dofs = vert_indices*2
        y_dofs = vert_indices*2+1
        
        
        if (x0[0]>=0.8) & (x0[1]==0.5):
            ext_dofs[ctr] = vert_indices[0]*2+1
            ctr=ctr+1
            #print(vert_indices[0])
        if (x1[0]>=0.8) & (x1[1]==0.5):
            ext_dofs[ctr] = vert_indices[1]*2+1
            ctr=ctr+1
        if (x2[0]>=0.8) & (x2[1]==0.5):
            ext_dofs[ctr] = vert_indices[2]*2+1 
            ctr=ctr+1
    
    ext_dofs = np.unique(ext_dofs)
    #print(ext_dofs)
    '''
    U = np.zeros(2*nv)
    
    # Get Ktan and Fint
    ext_dof, Ktan, Fint = FEM_Ktan_Fint(V,E, x_load, y_load, U, load_dir)
    
    
    
    #print(ext_dof)
    
    # Form Fext vector
    Fext = np.zeros(2*nv)
    Fext[ext_dof] = Fext_val
    
    
    
    # Body force
    Fb = body_force (E, V)   
    
    tol = 1e-12
    #is_boundary = ((np.abs(X) < tol))
    
    
    # BC for re-entrant corner
    is_boundary = ((np.abs(Y-1) < tol))
    
    bc_num = np.count_nonzero(is_boundary==True)
    
    bc_dofs_x = np.zeros((bc_num,),dtype=int)
    bc_dofs_y = np.zeros((bc_num,),dtype=int)
    ctr=0
    for i, vi in enumerate(is_boundary):
        if(is_boundary[i]==True):
            
            bc_dofs_x[ctr]= i*2
            bc_dofs_y[ctr]= i*2+1
            ctr=ctr+1
    
    
    all_dofs = np.arange(0,nv*2)
    bc_dofs_p = np.sort(np.concatenate((bc_dofs_x,bc_dofs_y)))
    bc_dofs_f = np.setdiff1d(all_dofs, bc_dofs_p)
    #print(bc_dofs_f)
    
    Fint_f = np.take(Fint, bc_dofs_f)
    Fext_f = np.take(Fext, bc_dofs_f)
    Fb_f = np.take(Fb, bc_dofs_f)
    #Ktan_1 = Ktan.tocsr()
    #Ktan_f = Ktan_1[bc_dofs_f,bc_dofs_f]
    Ktan_f = coo_submatrix_pull(Ktan, bc_dofs_f, bc_dofs_f)
    #print(Ktan_f)
    
    F_eq = Fext_f-Fb_f
    
    
    # Obtain displacement vector
    uhat = sla.spsolve(Ktan_f.tocsr(), F_eq)
    
    
    U[bc_dofs_f]=uhat

    return U 


 
if __name__ == '__main__':
    
    #from tempfile import TemporaryFile
    #u_exact_val = TemporaryFile()
    
    V, E = make_mesh(0,[],0)
    nv = len(V)
    ne = len(E)
    
    #print(V.shape)
    #print(E.shape)
    #print(E.max())
    # Extract the XY coordinates
    X, Y = V[:, 0], V[:, 1] 
    
    
    
    
    # Solve FEM problem to get structural displacements
    U = FEM_sol(V,E)
    print('original mesh d.o.f=',len(U))
    
    #np.save(u_exact_val, U)
    # Get the distribution of internal stress           
    Fint, norm_f_old = int_stress(E,V,U)
    
    # Get deformed coordinates to plot displaced mesh
    U_mat = U.reshape((nv,2))
    X_new = X+ U_mat[:,0]
    Y_new = Y+ U_mat[:,1]
    #print(U_mat)
    
    # Plot original mesh and deformed mesh
    plt.figure(1,figsize=(7,7))
    plt.gca().set_aspect("equal")
    plt.triplot(X, Y, E)
    plt.triplot(X_new, Y_new, E)
    plt.title('Undeformed and deformed FEA meshes=',y=1.05, fontsize=10)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('Figs/Undeformed and deformed FEA meshes.png')

    
    #Plot Internal stress distribution
   
    fig = plt.figure(2,figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(X, Y, Fint, triangles=E, cmap=plt.cm.jet, linewidth=0.2)
    plt.title('Internal stress distribution over domain',y=1.05, fontsize=10)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('Figs/Internal stress distribution over domain.png')
    
    
    
    # Call error estimator and mark the triangles for refinement
    mark,_,eta,_,esz=error_estimator(V,E, U)
    

    
    # Generate the refined mesh
    t=0
    #for t in range(0,5):
    V_new, E_new = make_mesh(1,mark,t)

    # Solve FEM on new mesh
    U_new = FEM_sol(V_new,E_new)
    print('refined mesh d.o.f=',len(U_new))
    # Get the new internal stress distribution
    Fint_new,norm_f_new = int_stress(E_new,V_new,U_new)
    X_new, Y_new = V_new[:, 0], V_new[:, 1]
    mark,_,eta,_,esz=error_estimator(V_new,E_new, U_new)

    # Plot new mesh
    plt.figure(t*3,figsize=(7,7))
    plt.gca().set_aspect("equal")
    plt.triplot(X_new, Y_new, E_new)
    plt.title('Adaptive mesh refinement=',y=1.05, fontsize=10)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('Figs/Adaptive mesh refinement.png')
    
    
    #Plot new internal stress distribution
    fig = plt.figure(4,figsize=(8,8))
    ax = fig.add_subplot(projection='3d')
    ax.plot_trisurf(X_new, Y_new, Fint_new, triangles=E_new, cmap=plt.cm.jet, linewidth=0.2)
    plt.title('Internal stress distribution over refined mesh',y=1.05, fontsize=10)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('Figs/Internal stress distribution over refined mesh.png')
    
    print('L_inf norm for original mesh', norm_f_old)
    print('L_inf norm for refined mesh', norm_f_new)
    
    '''
    U_mat_new = U_new.reshape((len(V_new),2))
    
    X_new_bar = X_new+ U_mat_new[:,0]
    Y_new_bar = Y_new+ U_mat_new[:,1]
    
    # Plot original mesh
    plt.figure(1,figsize=(7,7))
    plt.gca().set_aspect("equal")
    plt.triplot(X, Y, E)
    
    # Plot deformed mesh
    plt.figure(3,figsize=(7,7))
    plt.gca().set_aspect("equal")
    plt.triplot(X_new_bar, Y_new_bar, E_new)
    '''
    
    
    
    
    
    '''
    U[is_g_boundary] = 0.0
    rhs = b - A @ U
    
    rhs[is_boundary] = 0.0
    print(is_boundary)'''
