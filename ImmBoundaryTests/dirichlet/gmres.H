#include "AMReX_Reduce.H"
#define NMAX (AMREX_SPACEDIM*AMREX_SPACEDIM*AMREX_SPACEDIM+1)
#define NMAX2 (NMAX* NMAX)
#define KSPSIZE NMAX
#define NOPC 0
#define GJPC 1
#define SGSPC 2
#define NEARZERO 1e-15

//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
noprecond(
    amrex::Real MinvX[NMAX],
    amrex::Real /*A*/[NMAX][NMAX],
    amrex::Real X[NMAX],int n)
{
    for (int i = 0; i < n; i++) 
    {
        MinvX[i] = X[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
GJprecond(
    amrex::Real MinvX[NMAX], amrex::Real A[NMAX][NMAX], 
    amrex::Real X[NMAX], int n)
{
    // solve Dy = X
    for (int i = 0; i < n; i++) 
    {
        MinvX[i] = X[i] / A[i][i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
SGSprecond(
    amrex::Real MinvX[NMAX], 
    amrex::Real A[NMAX][NMAX], amrex::Real X[NMAX],int n)
{
    amrex::Real sum;
    amrex::Real y_d[NMAX], y_dd[NMAX];
    amrex::Real diag_inv[NMAX];
    amrex::Real diag[NMAX];

    // solve (D+L)D^-1(D+U)y = X
    // MinvX = y

    for (int i = 0; i < n; i++) 
    {
        diag[i] = A[i][i];
        diag_inv[i] = 1.0 / diag[i];
    }

    // solve (D+L) y' = X
    for (int i = 0; i < n; i++) 
    {
        sum = 0.0;
        for (int j = 0; j < i; j++) 
        {
            sum += A[i][j] * y_d[j];
        }
        y_d[i] = (X[i] - sum) * diag_inv[i];
    }

    // solve D^(-1)y'' = y' = (D+L)^(-1) X
    for (int i = 0; i < n; i++) 
    {
        y_dd[i] = diag[i] * y_d[i];
    }

    // solve (D+U) y = y'' = D (D+L)^(-1) X
    for (int i = n - 1; i >= 0; i--) 
    {
        sum = 0.0;
        for (int j = n - 1; j > i; j--) 
        {
            sum += A[i][j] * MinvX[j];
        }

        MinvX[i] = (y_dd[i] - sum) * diag_inv[i];
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
precond(
    amrex::Real MinvX[NMAX],
    amrex::Real A[NMAX][NMAX],
    amrex::Real X[NMAX],
    int n,
    int type)
{
    if (type == NOPC) 
    {
        noprecond(MinvX, A, X, n);
    } 
    else if (type == GJPC) 
    {
        GJprecond(MinvX, A, X, n);
    } 
    else if (type == SGSPC) 
    {
        SGSprecond(MinvX, A, X, n);
    } 
    else 
    {
        amrex::Abort("Unknown preconditioner type in BDF");
    }
}
//==============================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
findAX(amrex::Real AX[NMAX], 
       amrex::Real A[NMAX][NMAX], 
       amrex::Real X[NMAX], int n)
{
    for (int i = 0; i < n; i++) 
    {
        AX[i] = 0.0;
        for (int j = 0; j < n; j++) 
        {
            AX[i] += A[i][j] * X[j];
        }
    }
}
//========================================================================
#ifndef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
printvec(std::string str, amrex::Real v[NMAX], int n)
{
    amrex::Print() << "\n" << str << "\t";
    for (int i = 0; i < n; i++) 
    {
        amrex::Print() << v[i] << "\t";
    }
    amrex::Print() << "\n";
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
printmat(std::string str, amrex::Real mat[NMAX][NMAX], int n)
{
    amrex::Print() << "\n" << str << "\n";

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            amrex::Print() << mat[i][j] << "  ";
        }
        amrex::Print() << "\n";
    }
    amrex::Print() << "\n";
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
print_hessmat(std::string str, amrex::Real hmat[KSPSIZE + 1][KSPSIZE])
{
    amrex::Print() << "\n" << str << "\n";

    for (int i = 0; i < KSPSIZE + 1; i++) 
    {
        for (int j = 0; j < KSPSIZE; j++) 
        {
            amrex::Print() << hmat[i][j] << "  ";
        }
        amrex::Print() << "\n";
    }
    amrex::Print() << "\n";
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
print_arnoldivecs(std::string str, amrex::Real arnoldivecs[KSPSIZE + 1][NMAX], int n)
{
    amrex::Print() << "\n" << str << "\n";

    for (int i = 0; i < (KSPSIZE + 1); i++) 
    {
        for (int j = 0; j < n; j++) 
        {
            amrex::Print() << arnoldivecs[i][j] << "  ";
        }
        amrex::Print() << "\n";
    }
    amrex::Print() << "\n";
}
#endif
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
scalevector(amrex::Real v1[NMAX], int n, amrex::Real factor)
{
    for (int i = 0; i < n; i++) 
    {
        v1[i] = v1[i] * factor;
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
multiplyvectors(
    amrex::Real v1[NMAX], amrex::Real v2[NMAX], amrex::Real v12[NMAX], int n)
{
    for (int i = 0; i < n; i++) 
    {
        v12[i] = v1[i] * v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
findnorm(amrex::Real v1[NMAX], int n)
{
    amrex::Real norm = 0.0;
    for (int i = 0; i < n; i++) 
    {
        norm = norm + v1[i] * v1[i];
    }

    return (sqrt(norm));
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
innerproduct(amrex::Real v1[NMAX], amrex::Real v2[NMAX], int n)
{
    amrex::Real innerprod = 0.0;
    for (int i = 0; i < n; i++) 
    {
        innerprod = innerprod + v1[i] * v2[i];
    }

    return (innerprod);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
addvectors(
    amrex::Real v1[NMAX],
    amrex::Real v2[NMAX],
    amrex::Real v12[NMAX],
    int n,
    amrex::Real a,
    amrex::Real b)
{
    for (int i = 0; i < n; i++) 
    {
        v12[i] = a * v1[i] + b * v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
copyvector(
    amrex::Real v1[NMAX],
    amrex::Real v2[NMAX],int n) //(dest,source,size)
{
    for (int i = 0; i < n; i++) 
    {
        v1[i] = v2[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
getkspvector(
    amrex::Real v1[NMAX], amrex::Real kspvectors[KSPSIZE + 1][NMAX], int n,
    int vecnum)
{
    for (int i = 0; i < n; i++) 
    {
        v1[i] = kspvectors[vecnum][i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
setkspvector(
    amrex::Real vec[NMAX],
    amrex::Real kspvectors[KSPSIZE + 1][NMAX],
    int n,
    int vecnum)
{
    for (int i = 0; i < n; i++) 
    {
        kspvectors[vecnum][i] = vec[i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
isit_NaN(amrex::Real val)
{
    amrex::Real valbyval;
    amrex::Real valplus3;

    valbyval = val / val;
    valplus3 = val + 3.0;

    if (val != val) 
    {
        return (true);
    } 
    else if ((valbyval != valbyval) && (valplus3 == val)) 
    {
        return (true);
    } 
    else 
    {
        return (false);
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
arnoldi(
    int j,
    amrex::Real mat[NMAX][NMAX],
    amrex::Real kspvectors[KSPSIZE + 1][NMAX],
    amrex::Real Hessmat[KSPSIZE + 1][KSPSIZE],
    int n,
    int precond_type,
    amrex::Real tol,
    bool& nanflag)
{
    int i;
    amrex::Real Avj[NMAX], vj[NMAX], vi[NMAX];
    amrex::Real wj[NMAX], tempvec[NMAX];
    amrex::Real MinvAvj[NMAX];

    bool lucky; // when norm becomes 0,
                // KSP vectors is no longer linearly independent.
                // we would have gotten the best solution.

    lucky = false;
    nanflag = false;

    getkspvector(vj, kspvectors, n, j);
    findAX(Avj, mat, vj, n);
    precond(MinvAvj, mat, Avj, n, precond_type);
    copyvector(Avj, MinvAvj, n);
    // Avj is now M^-1 A vj
    // remember we are solving M^-1 A X = M^-1 b

    for (i = 0; i <= j; i++) {
        getkspvector(vi, kspvectors, n, i);
        Hessmat[i][j] = innerproduct(Avj, vi, n);
    }

    copyvector(wj, Avj, n);

    for (i = 0; i <= j; i++) 
    {
        getkspvector(vi, kspvectors, n, i);
        addvectors(wj, vi, tempvec, n, 1.0, -Hessmat[i][j]);
        copyvector(wj, tempvec, n);
    }

    Hessmat[j + 1][j] = findnorm(wj, n);

    if (Hessmat[j + 1][j] > tol) 
    {
        for (i = 0; i < n; i++) 
        {
            wj[i] = wj[i] / Hessmat[j + 1][j];
        }
        setkspvector(wj, kspvectors, n, j + 1);
    } 
    else 
    {
        if (!isit_NaN(Hessmat[j + 1][j])) 
        {
            lucky = true;
        } 
        else 
        {
            amrex::Abort("NaN or Inf detected in Hessenberg matrix\n");
            nanflag = true;
        }
    }

    return (lucky);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
givens_rotate(amrex::Real r, amrex::Real s, amrex::Real& cos, amrex::Real& sin)
{
    if (r == 0) 
    {
        cos = 0.0;
        sin = 1.0;
    } 
    else 
    {
        cos = r / std::sqrt(r * r + s * s);
        sin = s / std::sqrt(r * r + s * s);
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
triangularize(
    amrex::Real Hessmat[KSPSIZE + 1][KSPSIZE],
    amrex::Real cos_arr[KSPSIZE],
    amrex::Real sin_arr[KSPSIZE],
    amrex::Real beta_e1[KSPSIZE + 1],
    int k)
{

    amrex::Real temp, cos, sin;

    // apply for ith column
    for (int i = 0; i <= (k - 1); i++) 
    {
        temp = cos_arr[i] * Hessmat[i][k] + sin_arr[i] * Hessmat[i + 1][k];
        Hessmat[i + 1][k] =
        -sin_arr[i] * Hessmat[i][k] + cos_arr[i] * Hessmat[i + 1][k];
        Hessmat[i][k] = temp;
    }

    // update the next sin cos values for rotation
    givens_rotate(Hessmat[k][k], Hessmat[k + 1][k], cos, sin);

    // eliminate H(i + 1, i)
    Hessmat[k][k] = cos * Hessmat[k][k] + sin * Hessmat[k + 1][k];
    Hessmat[k + 1][k] = 0.0;

    cos_arr[k] = cos;
    sin_arr[k] = sin;

    beta_e1[k + 1] = -sin * beta_e1[k];
    beta_e1[k] = cos * beta_e1[k];

    amrex::Real error = std::abs(beta_e1[k + 1]);

    return (error);
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
uppertrisolve(
    amrex::Real Hessmat[KSPSIZE + 1][KSPSIZE],
    amrex::Real y[KSPSIZE],
    amrex::Real beta_e1[KSPSIZE + 1],
    int k)
{
    for (int i = k - 1; i >= 0; i--) 
    {
        y[i] = beta_e1[i];
        for (int j = i + 1; j <= k; j++) 
        {
            y[i] = y[i] - Hessmat[i][j] * y[j];
        }
        y[i] = y[i] / Hessmat[i][i];
    }
}
//========================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
performgmres(
    amrex::Real Mat[NMAX][NMAX],
    amrex::Real b[NMAX],
    amrex::Real x0[NMAX], // initial soln for GMRES
    amrex::Real x[NMAX],  // current solution for GMRES
    int n,
    int precond_type,
    int restart_it,
    int kspiters,
    amrex::Real tol,
    int printflag)
{
    amrex::Real beta{0.0};
    bool arnoldistopped{false};
    bool success{false};
    bool nanflag{false};

    amrex::Real v1[NMAX] = {0.0};
    amrex::Real tempvec[NMAX] = {0.0};
    amrex::Real v[NMAX] = {0.0};
    amrex::Real y[KSPSIZE] = {0.0};
    amrex::Real r[NMAX] = {0.0};
    amrex::Real r0[NMAX] = {0.0};
    amrex::Real Ax0[NMAX] = {0.0};
    amrex::Real Minvr[NMAX] = {0.0};
    amrex::Real kspvectors[KSPSIZE + 1][NMAX] = {{0.0}};
    amrex::Real Hessmat[KSPSIZE + 1][KSPSIZE] = {{0.0}};
    amrex::Real cos_arr[KSPSIZE] = {0.0};
    amrex::Real sin_arr[KSPSIZE] = {0.0};
    amrex::Real beta_e1[KSPSIZE + 1] = {0.0};
    amrex::Real residnorm{0.0};

    success = true;

    // finding r0
    findAX(Ax0, Mat, x0, n);
    addvectors(b, Ax0, r0, n, 1.0, -1.0);     // b-Ax0
    precond(Minvr, Mat, r0, n, precond_type); // Minv * r0
    copyvector(r0, Minvr, n);

    // initial residual is r0=M^-1(b-Ax0)
    // we are solving M^-1 A X = M^-1 b
    copyvector(r, r0, n);
    copyvector(x, x0, n);

#ifndef AMREX_USE_GPU
    if (printflag)
        amrex::Print() << "initial norm of residual:" << findnorm(r0, n) << "\n";
#endif

    for (int it = 0; it < restart_it; it++) 
    {
#ifndef AMREX_USE_GPU
        if (printflag)
            amrex::Print() << "restart iteration:" << it << "\n";
#endif

        copyvector(x0, x, n);
        beta = findnorm(r, n);
        if(beta < tol)
        {
           break;
        }

        for (int i = 0; i < n; i++) {
            v1[i] = r[i] / beta;
        }
        beta_e1[0] = beta;

        setkspvector(v1, kspvectors, n, 0);

        int optkspsize = kspiters;
        for (int kspdim = 0; kspdim < kspiters; kspdim++) 
        {
            copyvector(x, x0, n);
            // finds the ksp vector at kspdim+1
            arnoldistopped = arnoldi(
                kspdim, Mat, kspvectors, Hessmat, n, precond_type, tol * 0.01, nanflag);

            residnorm = triangularize(Hessmat, cos_arr, sin_arr, beta_e1, kspdim);

            if (nanflag) 
            {
                success = !nanflag;
                amrex::Abort("nan detected in arnoldi iteration\n");
                break;
            }

#ifndef AMREX_USE_GPU
            if (printflag)
                amrex::Print() << "norm of residual:" << residnorm << "\t" << tol
                << "\t" << kspdim << "\n";
#endif

            if (arnoldistopped) 
            {
                optkspsize = kspdim + 1;
#ifndef AMREX_USE_GPU
                if (printflag)
                    amrex::Print() << "lucky condition:" << kspdim << "\n";
#endif
                break;
            }

            if (residnorm <= tol) 
            {
                optkspsize = kspdim + 1;
                break;
            }
        }

        uppertrisolve(Hessmat, y, beta_e1, optkspsize);
        for (int i = 0; i < optkspsize; i++) 
        {
            //amrex::Print()<<"y:"<<y[i]<<"\n";
            getkspvector(v, kspvectors, n, i);
            addvectors(x, v, tempvec, n, 1.0, y[i]);
            copyvector(x, tempvec, n);
        }

        if ((residnorm <= tol) || (!success) || (arnoldistopped)) 
        {
            break;
        }
    }

    return (success);
}
//========================================================================
