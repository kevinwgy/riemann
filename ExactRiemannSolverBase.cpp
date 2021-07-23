#include<ExactRiemannSolverBase.h>
#include<array>
#include<utility> //std::pair
#include<bits/stdc++.h> //std::swap

#ifndef WITHOUT_BOOST
#include<boost/math/tools/roots.hpp>
using namespace boost::math::tools;
#endif

using std::pair;
using std::cout;
using std::endl;

//-----------------------------------------------------

ExactRiemannSolverBase::ExactRiemannSolverBase(std::vector<VarFcnBase*> &vf_, 
                            ExactRiemannSolverData &iod_riemann) : vf(vf_)
{
  maxIts_main          = iod_riemann.maxIts_main;
  maxIts_shock         = iod_riemann.maxIts_shock;
  numSteps_rarefaction = iod_riemann.numSteps_rarefaction;
  tol_main             = iod_riemann.tol_main;
  tol_shock            = iod_riemann.tol_shock;
  tol_rarefaction      = iod_riemann.tol_rarefaction;
  min_pressure         = iod_riemann.min_pressure;
  failure_threshold    = iod_riemann.failure_threshold;
  pressure_at_failure  = iod_riemann.pressure_at_failure;
}

//-----------------------------------------------------
/** Solves the one-dimensional Riemann problem. Extension of Kamm 2015 
 * to Two Materials. See KW's notes for details
 */
void
ExactRiemannSolverBase::ComputeRiemannSolution(int dir/*0~x,1~y,2~z*/, 
                            double *Vm, int idl /*"left" state*/, 
                            double *Vp, int idr /*"right" state*/, 
                            double *Vs, int &id /*solution at xi = 0 (i.e. x=0) */,
                            double *Vsm /*left 'star' solution*/,
                            double *Vsp /*right 'star' solution*/)
{
  // Convert to a 1D problem (i.e. One-Dimensional Riemann)
  double rhol  = Vm[0];
  double ul    = Vm[dir+1];
  double pl    = Vm[4];
  double rhor  = Vp[0];
  double ur    = Vp[dir+1];
  double pr    = Vp[4];
  //fprintf(stderr,"1DRiemann: left = %e %e %e (%d) : right = %e %e %e (%d)\n", rhol, ul, pl, idl, rhor, ur, pr, idr);

  double el = vf[idl]->GetInternalEnergyPerUnitMass(rhol, pl);
  double cl = vf[idl]->ComputeSoundSpeedSquare(rhol, el);

  if(rhol<=0 || cl<0) {
    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l)." 
            " rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
            rhol, ul, pl, el, cl, idl);
    exit(-1);
  } else
    cl = sqrt(cl);

  double er = vf[idr]->GetInternalEnergyPerUnitMass(rhor, pr);
  double cr = vf[idr]->ComputeSoundSpeedSquare(rhor, er);

  if(rhor<=0 || cr<0) {
    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(r)."
            " rho = %e, u = %e, p = %e, e = %e, c^2 = %e, ID = %d.\n",
            rhor, ur, pr, er, cr, idr);
    exit(-1);
  } else
    cr = sqrt(cr);


  // Declare variables in the "star region"
  double p0, ul0, ur0, rhol0, rhor0;
  double p1, ul1, ur1, rhol1, rhor1; //Secand Method ("k-1","k" in Kamm, (19))
  double p2, ul2, ur2, rhol2, rhor2; // "k+1"
  double f0, f1, f2; //difference between ul and ur


  // -------------------------------
  // Now, Solve The Riemann Problem
  // -------------------------------
  bool success = true;

  // monitor if the solution involves a transonic rarefaction. This is special as the solution 
  // at xi = x = 0 is within the rarefaction fan.
  bool trans_rare = false;
  double Vrare_x0[3]; //rho, u, and p at x = 0, in the case of a transonic rarefaction


  // A Trivial Case
  if(ul == ur && pl == pr) {
    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol, rhor, ul, pl, 
                   trans_rare, Vrare_x0, //inputs
                   Vs, id, Vsm, Vsp/*outputs*/);
    return;
  }

  // -------------------------------
  // Step 1: Initialization
  //         (find initial interval [p0, p1])
  // -------------------------------
  success = FindInitialInterval(rhol, ul, pl, el, cl, idl, rhor, ur, pr, er, cr, idr, /*inputs*/
                                p0, rhol0, rhor0, ul0, ur0,
                                p1, rhol1, rhor1, ul1, ur1/*outputs*/);
  /* our convention is that p0 < p1 */

  if(!success) { //failed to find a bracketing interval. Output the state corresponding smallest "f"

    // get sol1d, trans_rare and Vrare_x0
#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif
    success = ComputeRhoUStar(1, rhol, ul, pl, p1, idl, rhol0, rhol0*1.1, rhol2, ul2,
                  &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);
    success = success && ComputeRhoUStar(3, rhor, ur, pr, p1, idr, rhor0, rhol0*1.1, rhor2, ur2,
                    &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      cout << "*** Error: Riemann solver failed (really strange. software bug)." << endl;
      exit(-1);
    }

    FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, 0.5*(ul2+ur2), p1,
                     trans_rare, Vrare_x0, /*inputs*/
                     Vs, id, Vsm, Vsp /*outputs*/);

    return;
  }

  f0 = ul0 - ur0;
  f1 = ul1 - ur1;

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Found initial interval: p0 = " << p0 << ", f0 = " << f0 << ", p1 = " << p1 << ", f1 = " << f1 << endl;
#endif
                                          
  // -------------------------------
  // Step 2: Main Loop (Secant Method, Safeguarded) 
  // -------------------------------
  double denom = 0; 

  int iter = 0;
  double err_p = 1.0, err_u = 1.0;

  //p2 (and f2) is always the latest one
  p2 = p1; 
  f2 = f1;

#if PRINT_RIEMANN_SOLUTION == 1
  sol1d.clear();
#endif

  for(iter=0; iter<maxIts_main; iter++) {

    // 2.1: Update p using the Brent method (safeguarded secant method)
    denom = f1 - f0;
    if(denom == 0) {
      cout << "*** Error: Division-by-zero while using the secant method to solve the Riemann problem." << endl;
      cout << "           left state: " << rhol << ", " << ul << ", " << pl << ", " << idl << " | right: " 
           << rhor << ", " << ur << ", " << pr << ", " << idr << endl;
      cout << "           dir = " << dir << ", f0 = " << f0 << ", f1 = " << f1 << endl;
      exit(-1);
    }

    p2 = p2 - f2*(p1-p0)/denom;  // update p2

    if(p2<=p0 || p2>=p1) //discard and switch to bisection
      p2 = 0.5*(p0+p1);

    //fprintf(stderr,"iter = %d, p0 = %e, p1 = %e, p2 = %e, f0 = %e, f1 = %e.\n", iter, p0, p1, p2, f0, f1);



    // 2.2: Calculate ul2, ur2 
    success = ComputeRhoUStar(1, rhol, ul, pl, p2, idl/*inputs*/, 
                  rhol0, rhol1/*initial guesses for Hugo. eq.*/,
                  rhol2, ul2/*outputs*/, 
                  &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      fprintf(stderr,"*** Error: Exact Riemann solver failed. left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
              rhol, ul, pl, idl, rhor, ur, pr, idr);
      exit(-1);
    }

    success = ComputeRhoUStar(3, rhor, ur, pr,  p2, idr/*inputs*/, 
                    rhor0, rhor1/*initial guesses for Hugo. erq.*/,
                    rhor2, ur2/*outputs*/,
                    &trans_rare, Vrare_x0/*filled only if found a trans. rarefaction*/);

    if(!success) {
      fprintf(stderr,"*** Error: Exact Riemann solver failed (2). left: %e %e %e (%d) | right: %e %e %e (%d).\n", 
              rhol, ul, pl, idl, rhor, ur, pr, idr);
      exit(-1);
    }

    f2 = ul2 - ur2;
    
 
    // 2.3: Update for the next iteration
    if(f0*f2<0.0) {
      p1 = p2;
      f1 = f2;
      rhol1 = rhol2;
      rhor1 = rhor2;
    } else {
      p0 = p2;
      f0 = f2;
      rhol0 = rhol2;
      rhor0 = rhor2;
    }


    // 2.4: Check stopping criterion
    err_p = fabs(p1 - p0)/std::max(fabs(pl + 0.5*rhol*ul*ul), fabs(pr + 0.5*rhor*ur*ur));
    err_u = fabs(f2)/std::max(cl, cr);

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "Iter " << iter << ": err_p = " << err_p << ", err_u = " << err_u << "." << endl;
#endif

    if( (err_p < tol_main && err_u < tol_main) || (err_p < tol_main*1e-3) || (err_u < tol_main*1e-3) )
      break; // converged

    trans_rare = false; //reset

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.clear();
#endif

  }



  if(iter == maxIts_main) {
    cout << "*** Error: Exact Riemann solver failed to converge. err_p = " << err_p
         << ", err_u = " << err_u << "." << endl;
    cout << "    Vm = [" << Vm[0] << ", " << Vm[1] << ", " << Vm[2] << ", " << Vm[3] << ", " << Vm[4] << "] (" << idl
         << "),  Vp = [" << Vp[0] << ", " << Vp[1] << ", " << Vp[2] << ", " << Vp[3] << ", " << Vp[4] << "] (" << idr
         << ")" << endl;
      
    exit(-1);
  }
  

  // -------------------------------
  // Step 3: Find state at xi = x = 0 (for output)
  // -------------------------------
  double u2 = 0.5*(ul2 + ur2);
  FinalizeSolution(dir, Vm, Vp, rhol, ul, pl, idl, rhor, ur, pr, idr, rhol2, rhor2, u2, p2, 
                   trans_rare, Vrare_x0, //inputs
                   Vs, id, Vsm, Vsp/*outputs*/);

}

//-----------------------------------------------------

void
ExactRiemannSolverBase::FinalizeSolution(int dir, double *Vm, double *Vp,
                                         double rhol, double ul, double pl, int idl, 
                                         double rhor, double ur, double pr, int idr, 
                                         double rhol2, double rhor2, double u2, double p2,
                                         bool trans_rare, double Vrare_x0[3], /*inputs*/
                                         double *Vs, int &id, double *Vsm, double *Vsp /*outputs*/)
{
  // find material id at xi = x = 0
  if(u2>=0)
    id = idl;
  else
    id = idr;

#if PRINT_RIEMANN_SOLUTION == 1
  // the 2-wave
  sol1d.push_back(vector<double>{u2 - std::max(1e-6, 0.001*fabs(u2)), rhol2, u2, p2, (double)idl});
  sol1d.push_back(vector<double>{u2, rhor2, u2, p2, (double)idr});
#endif

  Vs[0] = Vs[1] = Vs[2] = Vs[3] = Vs[4] = 0.0;

  if(trans_rare) {
    Vs[0]     = Vrare_x0[0];
    Vs[dir+1] = Vrare_x0[1];
    Vs[4]     = Vrare_x0[2];
  }
  else { 
    //find state variables at xi = x = 0

    if(u2>=0) { //either Vl or Vlstar --- check the 1-wave

      bool is_star_state = false;

      if(pl >= p2) {//1-wave is rarefaction
        double el2 = vf[idl]->GetInternalEnergyPerUnitMass(rhol2, p2);
        double cl2 = vf[idl]->ComputeSoundSpeedSquare(rhol2, el2);

        if(rhol2<=0 || cl2<0) {
          fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(l2)."
                  " rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
                  rhol2, pl, el2, cl2, idl);
          exit(-1);
        } else
          cl2 = sqrt(cl2);

        if(u2 - cl2 <= 0) //rarefaction tail speed
          is_star_state = true;
      } 
      else {//1-wave is shock
        double us = (rhol2*u2 - rhol*ul)/(rhol2 - rhol); //shock speed
        if(us <= 0)      
          is_star_state = true;
      }

      if(is_star_state) {
        Vs[0]     = rhol2;
        Vs[dir+1] = u2;
        Vs[4]     = p2;
      } else {
        Vs[0]     = rhol;
        Vs[dir+1] = ul;
        Vs[4]     = pl;
      }

    } else { //either Vr or Vrstar --- check the 3-wave

      bool is_star_state = false;

      if(pr >= p2) {//3-wave is rarefaction
        double er2 = vf[idr]->GetInternalEnergyPerUnitMass(rhor2, p2);
        double cr2 = vf[idr]->ComputeSoundSpeedSquare(rhor2, er2);

        if(rhor2<=0 || cr2<0) {
          fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed) in ComputeRiemannSolution(r2)." 
                  " rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
                  rhor2, p2, er2, cr2, idr);
          exit(-1);
        } else
          cr2 = sqrt(cr2);

        if(u2 - cr2 >= 0)
          is_star_state = true;
      }
      else {//3-wave is shock
        double us = (rhor2*u2 - rhor*ur)/(rhor2 - rhor); //shock speed
        if(us >= 0)
          is_star_state = true;
      }

      if(is_star_state) {
        Vs[0]     = rhor2;
        Vs[dir+1] = u2;
        Vs[4]     = p2;
      } else {
        Vs[0]     = rhor;
        Vs[dir+1] = ur;
        Vs[4]     = pr;
      }

    }
  }

  // determine the tangential components of velocity -- upwinding
  int k;
  if(u2>0) {
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      Vs[k] = Vm[k];
    }
  } else if(u2<0) {
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      Vs[k] = Vp[k];
    }
  } else {//u2 == 0
    for(int i=1; i<=2; i++) {
      k = (dir + i)%3 + 1;
      Vs[k] = 0.5*(Vm[k]+Vp[k]);
    }
  }


  // determine Vsm and Vsp, i.e. the star states on the minus and plus sides of the contact discontinuity
  Vsm[0]     = rhol2;
  Vsm[dir+1] = u2;
  Vsm[4]     = p2;
  Vsp[0]     = rhor2;
  Vsp[dir+1] = u2;
  Vsp[4]     = p2;
  for(int i=1; i<=2; i++) { //tangential velocity
    k = (dir + i)%3 + 1;
    Vsm[k] = Vm[k];
    Vsp[k] = Vp[k];
  }


#if PRINT_RIEMANN_SOLUTION == 1
  std::sort(sol1d.begin(), sol1d.end(), 
            [](vector<double> v1, vector<double> v2){return v1[0]<v2[0];});
  int last = sol1d.size()-1;
  double xi_span = sol1d[last][0] - sol1d[0][0];
  sol1d.insert(sol1d.begin(), vector<double>{sol1d[0][0]-xi_span, sol1d[0][1], sol1d[0][2], sol1d[0][3], sol1d[0][4]});
  last++;
  sol1d.push_back(vector<double>{sol1d[last][0]+xi_span, sol1d[last][1], sol1d[last][2], sol1d[last][3], sol1d[last][4]});

  FILE* solFile = fopen("RiemannSolution.txt", "w");
  print(solFile, "## One-Dimensional Riemann Problem.\n");
  print(solFile, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
        rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(solFile, "## xi(x/t) | density | velocity | pressure | internal energy per mass | material id\n");

  for(auto it = sol1d.begin(); it != sol1d.end(); it++) 
    print(solFile,"% e    % e    % e    % e    % e    % d\n", (*it)[0], (*it)[1], (*it)[2], (*it)[3], 
                   vf[(int)(*it)[4]]->GetInternalEnergyPerUnitMass((*it)[1], (*it)[3]), (int)(*it)[4]);

  fclose(solFile);
#endif

}

//----------------------------------------------------------------------------------
//! find a bracketing interval [p0, p1] (f0*f1<=0)
bool
ExactRiemannSolverBase::FindInitialInterval(double rhol, double ul, double pl, double el, double cl, int idl,
                            double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
                            double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0,
                            double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  /*convention: p0 < p1*/

  bool success = true;

  // Step 1: Find two feasible points (This step should never fail)
  success = FindInitialFeasiblePoints(rhol, ul, pl, el, cl, idl, rhor, ur, pr, er, cr, idr, /*inputs*/
                                      p0, rhol0, rhor0, ul0, ur0, p1, rhol1, rhor1, ul1, ur1/*outputs*/);

#if PRINT_RIEMANN_SOLUTION == 1
  fprintf(stderr, "Found two initial points: p0 = %e, f0 = %e, p1 = %e, f1 = %e.\n", p0, ul0-ur0, p1, ul1-ur1);
  fprintf(stderr, "Searching for a bracketing interval...\n");
#endif

  if(!success) //This should never happen (unless user's inputs have errors)!
    exit(-1);
  
  // Step 2: Starting from the two feasible points, find a bracketing interval
  //         This step may fail, which indicates a solution may not exist for arbitrary left & right states
  //         If this happens, return the point with smallest absolute value of "f"

  double p2, rhol2, rhor2, ul2, ur2;
  double f0, f1;
  int i;

  // These are variables corresponding to the smallest magnitude of "f" --- used only if a bracketing interval
  // cannot be found. This is to just to minimize the chance of code crashing...
  double fmin, p_fmin, rhol_fmin, rhor_fmin, ul_fmin, ur_fmin;
  if(fabs(ul0-ur0) < fabs(ul1-ur1)) {
    fmin = fabs(ul0 - ur0);
    p_fmin = p0;  rhol_fmin = rhol0;  rhor_fmin = rhor0;  ul_fmin = ul0;  ur_fmin = ur0;
  } else {
    fmin = fabs(ul1 - ur1);
    p_fmin = p1;  rhol_fmin = rhol1;  rhor_fmin = rhor1;  ul_fmin = ul1;  ur_fmin = ur1;
  }

  for(i=0; i<maxIts_main; i++) {
    
    f0 = ul0 - ur0;
    f1 = ul1 - ur1;

    if(f0*f1<=0.0)
      return true;

    // find a physical p2 that has the opposite sign
    if(fabs(f0-f1)>1e-9) {
      p2 = p1 - f1*(p1-p0)/(f1-f0); //the Secant method
      if(p2<p0)
        p2 -= 0.1*(p1-p0);
      else //p2 cannot be between p0 and p1, so p2>p1
        p2 += 0.1*(p1-p0);
    } else {//f0 == f1
      p2 = 1.1*p1; 
      //p2 = p1 + *(p1-p0); 
    }

    if(p2<min_pressure || i==int(maxIts_main/2) ) {//does not look right. reset to a small non-negative pressure
      p2 = 1.0e-8; 
    }

    success = ComputeRhoUStar(1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
    // compute the 3-wave only if the 1-wave is succeeded
    success = success && ComputeRhoUStar(3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);

    if(!success) {

#if PRINT_RIEMANN_SOLUTION == 1
      fprintf(stderr, "  -- p2 = %e (failed)\n", p2);
#endif
      //move closer to [p0, p1]
      for(int j=0; j<maxIts_main; j++) {
        if(p2<p0)     
          p2 = p0 - 0.5*(p0-p2);
        else //p2>p1
          p2 = p1 + 0.5*(p2-p1);

        success = ComputeRhoUStar(1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
        // compute the 3-wave only if the 1-wave is succeeded
        success = success && ComputeRhoUStar(3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);


        if(success)
          break;
      }
    }

    if(!success)
      break;

    // check & update fmin (only used if the Riemann solver fails -- a "failsafe" feature)
    if(fabs(ul2-ur2)<fmin) {
      fmin = fabs(ul2-ur2);
      p_fmin = p2;  rhol_fmin = rhol2;  rhor_fmin = rhor2;  ul_fmin = ul2;  ur_fmin = ur2;
    }

    // update p0 or p1
    if(p2<p0) {
      p1 = p0;  rhol1 = rhol0;  rhor1 = rhor0;  ul1 = ul0;  ur1 = ur0;
      p0 = p2;  rhol0 = rhol2;  rhor0 = rhor2;  ul0 = ul2;  ur0 = ur2;
    } else {//p2>p1
      p0 = p1;  rhol0 = rhol1;  rhor0 = rhor1;  ul0 = ul1;  ur0 = ur1;
      p1 = p2;  rhol1 = rhol2;  rhor1 = rhor2;  ul1 = ul2;  ur1 = ur2;
    }

#if PRINT_RIEMANN_SOLUTION == 1
    fprintf(stderr, "  -- p0 = %e, f0 = %e, p1 = %e, f1 = %e (success)\n", p0, ul0-ur0, p1, ul1-ur1);
#endif

  }

  if(!success || i==maxIts_main) {
    cout << "Warning: Exact Riemann solver failed. (Unable to find a bracketing interval) " << endl;
    cout << "   left: " << std::setprecision(10) << rhol << ", " << std::setprecision(10) << ul << ", " 
                           << std::setprecision(10) << pl << " (" << idl << "); right: "
                           << std::setprecision(10) << rhor << ", " << std::setprecision(10) << ur << ", " 
                           << std::setprecision(10) << pr << " (" << idr << ")."
                           << " Residual (|ulstar-urstar|): " << std::setprecision(10) << fmin << endl;
    if(fmin<failure_threshold*fabs(ul-ur)) {
      cout << "*** Best approximate solution: rhols = " << rhol_fmin << ", ps = " << p_fmin << ", us = ("
           << ul_fmin << "(l) + " << ur_fmin << "(r))/2, rhors = " << rhor_fmin << "." << endl;
      p0    = p1    = p_fmin;
      rhol0 = rhol1 = rhol_fmin;
      rhor0 = rhor1 = rhor_fmin;
      ul0   = ul1   = ul_fmin; 
      ur0   = ur1   = ur_fmin;
    } else { //it could be that the Riemann problem has no solution!
      p2 = pressure_at_failure;
      success = ComputeRhoUStar(1, rhol, ul, pl, p2, idl, rhol0, rhol1, rhol2, ul2);
      // compute the 3-wave only if the 1-wave is succeeded
      success = success && ComputeRhoUStar(3, rhor, ur, pr,  p2, idr, rhor0, rhor1, rhor2, ur2);
      if(success) {
        cout << "*** Prescribed solution: rhols = " << rhol2 << ", ps = " << p2 << ", us = ("
             << ul2 << "(l) + " << ur2 << "(r))/2, rhors = " << rhor2 << "." << endl;
        p0    = p1    = p2;
        rhol0 = rhol1 = rhol2; 
        rhor0 = rhor1 = rhor2;
        ul0   = ul1   = ul2; 
        ur0   = ur1   = ur2;
      } else {
        cout << "*** Best approximation: rhols = " << rhol_fmin << ", ps = " << p_fmin << ", us = ("
             << ul_fmin << "(l) + " << ur_fmin << "(r))/2, rhors = " << rhor_fmin << "." << endl;
        p0    = p1    = p_fmin;
        rhol0 = rhol1 = rhol_fmin;
        rhor0 = rhor1 = rhor_fmin;
        ul0   = ul1   = ul_fmin; 
        ur0   = ur1   = ur_fmin;
      }
    }

    return false;
  }

  return true;
}

//----------------------------------------------------------------------------------

bool
ExactRiemannSolverBase::FindInitialFeasiblePoints(double rhol, double ul, double pl, double el, double cl, int idl,
                            double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
                            double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
                            double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  double dp;
  int found = 0;
  bool success = true;

  // Method 1: Use the acoustic theory (Eqs. (20)-(22) of Kamm) to find p0, p1
  found = FindInitialFeasiblePointsByAcousticTheory(rhol, ul, pl, el, cl, idl, 
              rhor, ur, pr, er, cr, idr, /*inputs*/
              p0, rhol0, rhor0, ul0, ur0, p1, rhol1, rhor1, ul1, ur1/*outputs*/);

  if(found==2)
    return true; //yeah

  if(found==1) //the first one (p0) is good --> only need to find p1
    goto myLabel;

  // Method 2: based on dp (fixed width search)

  // 2.1. find the first one (p0)
  dp = (pl!=pr) ? fabs(pl-pr) : 0.5*pl;
  for(int i=0; i<maxIts_main; i++) {
    p0 = std::min(pl,pr) + 0.01*(i+1)*(i+1)*dp;
    if(p0<min_pressure)
      p0 = pressure_at_failure; 
    success = ComputeRhoUStar(1, rhol, ul, pl, p0, idl, 
                              rhol, (p0>pl) ? rhol*1.1 : rhol*0.9,
                              rhol0, ul0);
    success = success && ComputeRhoUStar(3, rhor, ur, pr, p0, idr, 
                                         rhor, (p0>pr) ? rhor*1.1 : rhor*0.9,
                                         rhor0, ur0);
    if(success)
      break;
  }
  if(!success) {//search in the opposite direction
    for(int i=0; i<maxIts_main; i++) {
      p0 = std::min(pl,pr) - 0.01*(i+1)*(i+1)*dp;
      if(p0<min_pressure)
        p0 = pressure_at_failure; 
      success = ComputeRhoUStar(1, rhol, ul, pl, p0, idl, 
                                rhol, (p0>pl) ? rhol*1.1 : rhol*0.9,
                                rhol0, ul0);
      success = success && ComputeRhoUStar(3, rhor, ur, pr, p0, idr, 
                                           rhor, (p0>pr) ? rhor*1.1 : rhor*0.9,
                                           rhor0, ur0);
      if(success)
        break;
    } 
  }

  if(!success) {
    fprintf(stderr,"*** Error: Failed to find the first initial guess (p0) in the 1D Riemann solver (it = %d).\n",
                    maxIts_main);
    return false;
  }
      
myLabel:
  // 2.2. find the second one (p1)
  dp = std::min(fabs(p0-pl), fabs(p0-pr));
  for(int i=0; i<maxIts_main; i++) {
    p1 = p0 + 0.01*(i+1)*(i+1)*dp;
    success = ComputeRhoUStar(1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
    success = success && ComputeRhoUStar(3, rhor, ur, pr, p1, idr, rhor, rhor0, rhor1, ur1);
    if(success)
      break;
  }
  if(!success) //search in the opposite direction
    for(int i=0; i<maxIts_main; i++) {
      p1 = p0 - 0.01*(i+1)*(i+1)*dp;
      if(p1<min_pressure)
        p1 = pressure_at_failure*1000.0; //so it is not the same as p0
      success = ComputeRhoUStar(1, rhol, ul, pl, p1, idl, rhol, rhol0, rhol1, ul1);
      success = success && ComputeRhoUStar(3, rhor, ur, pr, p1, idr, rhor, rhor0, rhor1, ur1);
      if(success)
        break;
    } 
  if(!success) {
    fprintf(stderr,"*** Error: Failed to find the second initial guess (p1) in the 1D Riemann solver (it = %d).\n",
                    maxIts_main);
    return false;
  }
 
  // Make sure p0<p1
  if(p0>p1) {
    std::swap(p0,p1);
    std::swap(rhol0, rhol1);
    std::swap(rhor0, rhor1);
    std::swap(ul0, ul1);
    std::swap(ur0, ur1);
  } 

  return true;
}

//----------------------------------------------------------------------------------

int
ExactRiemannSolverBase::FindInitialFeasiblePointsByAcousticTheory(double rhol, double ul, 
                            double pl, double el, double cl, int idl,
                            double rhor, double ur, double pr, double er, double cr, int idr, /*inputs*/
                            double &p0, double &rhol0, double &rhor0, double &ul0, double &ur0, 
                            double &p1, double &rhol1, double &rhor1, double &ul1, double &ur1/*outputs*/)
{
  int found = 0;
  bool success = true;

  // 1.1: Initialize p0 using acoustic theory ((20) of Kamm)
  double Cl = rhol*cl; //acoustic impedance
  double Cr = rhor*cr; //acoustic impedance
  p0 = (Cr*pl + Cl*pr + Cl*Cr*(ul - ur))/(Cl + Cr);

  success = ComputeRhoUStar(1, rhol, ul, pl, p0, idl/*inputs*/,
                rhol, (p0>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
                rhol0, ul0/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, rhor, ur, pr, p0, idr/*inputs*/,
                rhor, (p0>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
                rhor0, ur0/*outputs*/);
  if(!success)
    return found;

  found = 1; //found p0!

  // 1.2. Initialize p1 ((21)(22) of Kamm)
  double Clbar = (ul0 == ul) ? Cl : fabs(p0 - pl)/fabs(ul0 - ul);
  double Crbar = (ur0 == ur) ? Cr : fabs(p0 - pr)/fabs(ur0 - ur);
  p1 = (Crbar*pl + Clbar*pr + Clbar*Crbar*(ul - ur))/(Clbar + Crbar);
  double tmp = std::max(fabs(p0), fabs(p1));
  if(fabs(p1 - p0)/tmp<1.0e-8)
    p1 = p0 + 1.0e-8*tmp; //to avoid f0 = f1 (divide-by-zero)

  success = ComputeRhoUStar(1, rhol, ul, pl, p1, idl/*inputs*/,
                rhol, rhol0/*initial guesses for Hugo. eq.*/,
                rhol1, ul1/*outputs*/);
  if(!success)
    return found;

  success = ComputeRhoUStar(3, rhor, ur, pr, p1, idr/*inputs*/,
                rhor, rhor0/*initial guesses for Hugo. eq.*/,
                rhor1, ur1/*outputs*/);
  if(!success)
    return found;

  found = 2; //found p0 and p1!
  return found;
}

//----------------------------------------------------------------------------------

//! Connect the left/right initial state with the left/right star state (the 1-wave or 3-wave)
bool  //true: success  | false: failure
ExactRiemannSolverBase::ComputeRhoUStar(int wavenumber /*1 or 3*/,
                                        double rho, double u, double p, double ps, int id/*inputs*/,
                                        double rhos0, double rhos1/*initial guesses for Hugo. eq.*/,
                                        double &rhos, double &us/*outputs*/, 
                                        bool *trans_rare, double *Vrare_x0/*filled only if found tran rf*/)
{
  // default
  rhos = rho;
  us   = u;

  if(p > ps) {//rarefaction --- numerical integration

    double dp_max = 1.25*(p-ps)/numSteps_rarefaction;
    double dp_target = dp_max/1.25;

    // initialize drho based on linear approximation
    double e = vf[id]->GetInternalEnergyPerUnitMass(rho,p);
    double dpdrho = vf[id]->GetDpdrho(rho,e);
    double drho_op1 = (p-ps)/dpdrho; //use EOS
    double drho_op2 = rho/(numSteps_rarefaction*2.5); //initial step size
    double drho = std::min(drho_op1, drho_op2);

    // prepare for numerical integration
    double rhos_0 = rho, us_0 = u, ps_0 = p, xi_0; //start point of each step
    double rhos_1 = rho, us_1 = u, ps_1 = p, xi_1; //end point of each step
    double dp;

    double c = vf[id]->ComputeSoundSpeedSquare(rho, e);

    if(rho<=0 || c<0) {
      fprintf(stderr,"Warning: Negative density or c^2 (square of sound speed) in ComputeRhoUStar." 
                     "rho = %e, p = %e, e = %e, c^2 = %e, id = %d.\n",
                     rho, p, e, c, id);
      return false; //failure
    } else
      c = sqrt(c);

    double xi = (wavenumber == 1) ? u - c : u + c; // xi = u -/+ c
    xi_0 = xi;

#if PRINT_RIEMANN_SOLUTION == 1
    sol1d.push_back(vector<double>{xi, rho, u, p, (double)id});
#endif

    //fprintf(stderr,"rho = %e, p = %e, ps = %e\n", rho, p, ps);
    // integration by Runge-Kutta 4
    bool done = false;
    for(int i=0; i<numSteps_rarefaction*5; i++) {

      bool success = Rarefaction_OneStepRK4(wavenumber/*1 or 3*/, id,
                             rhos_0, us_0, ps_0 /*start state*/, drho /*step size*/,
                             rhos_1, us_1, ps_1, xi_1 /*output: end state*/); 
//      fprintf(stderr,"RK4 step: rhos_0 = %e, us_0 = %e, ps_0 = %e, drho = %e, rhos_1 = %e, us_1 = %e, ps_1 = %e | ps = %e | success = %d.\n",
//              rhos_0, us_0, ps_0, drho, rhos_1, us_1, ps_1, ps, success);
      if(!success) {
        drho /= 2.0;
        continue;
      }

      dp = ps_0 - ps_1;

      //Check if it went too far. If so, rewind and reduce step size
//      fprintf(stderr,"RK4 step: -- dp = %e, dp_target = %e, dp_max = %e, ps_1 - ps = %e, tol_rarefaction = %e\n",
//              dp, dp_target, dp_max, ps_1 - ps, tol_rarefaction);

      if(dp > dp_max) {
        drho = drho/dp*dp_target;
        continue;
      }
      if(ps_1 - ps < -tol_rarefaction) { //went over the hill...
        if(dp!=0)
          drho = drho/dp*(ps_0 - ps);
        else
          drho /= 2.0;
        continue;
      } 

#if PRINT_RIEMANN_SOLUTION == 1
      sol1d.push_back(vector<double>{xi_1, rhos_1, us_1, ps_1, (double)id});
#endif

      if(trans_rare && Vrare_x0 && xi_0*xi_1<=0) {//transonic rarefaction, crossing x = xi = 0
        *trans_rare = true;
        double w0 = fabs(xi_1), w1 = fabs(xi_0);
        double ww = w0 + w1;
        w0 /= ww;
        w1 /= ww;
        Vrare_x0[0] = w0*rhos_0 + w1*rhos_1;
        Vrare_x0[1] = w0*us_0   + w1*us_1;
        Vrare_x0[2] = w0*ps_0   + w1*ps_1;

#if PRINT_RIEMANN_SOLUTION == 1
        sol1d.push_back(vector<double>{0.0, Vrare_x0[0], Vrare_x0[1], Vrare_x0[2], (double)id});
#endif

      }

      //fprintf(stderr,"drho = %e, rho: %e -> %e,  u: %e -> %e,  p: %e -> %e\n", drho, rhos_0, rhos_1, us_0, us_1, ps_0, ps_1);

      // Check if we have reached the final pressure ps
      if(fabs(ps_1 - ps) <= tol_rarefaction) {
        rhos = rhos_1;
        us   = us_1;
      
        if(vf[id]->CheckState(rhos,ps)) {
#if PRINT_RIEMANN_SOLUTION == 1
          cout << "Rarefaction solver reached a nonphysical state!" << endl;         
#endif
          return false;
        }

#if PRINT_RIEMANN_SOLUTION == 1
        cout << "  " << wavenumber << "-wave: rarefaction, integration completed in " << i << " steps" << endl;
#endif
        done = true;
        break; //done!
      }

      // If the solver gets here, it means it hasn't reached the final pressure ps.
      // Adjust step size, then update state
      //
//      fprintf(stderr,"RK4 step: adjusting drho. old drho: %e, rhos_0 - rhos_1 = %e, a = %e, b = %e.\n", drho, rhos_0 - rhos_1, (rhos_0-rhos_1)/dp*std::min(dp_target,ps_1-ps), drho*4.0);

      drho = std::min( drho/dp*std::min(dp_target,ps_1-ps), //don't go beyond ps
                       drho*4.0); //don't increase too much in one step
       
      rhos_0 = rhos_1;
      us_0   = us_1;
      ps_0   = ps_1;
      xi_0   = xi_1;

    }

    if(!done) {
      if(vf[id]->CheckState(rhos_1,ps_1)) {
#if PRINT_RIEMANN_SOLUTION == 1
        cout << "  " << wavenumber << "-wave: rarefaction, solver failed (unphysical state: rhos = "
             << rhos_1 << ", ps = " << ps_1 << "!)" << endl;
#endif
        return false; //failed
      } else {
#if PRINT_RIEMANN_SOLUTION == 1
        cout << "  " << wavenumber << "-wave: rarefaction, solver did not converge (final sol.: rhos_1 = "
             << rhos_1 << ", ps_1 = " << ps_1 << "; inputs: rho = " << rho << ", p = " << p << ", ps = " << ps << ")" << endl;
#endif
        return true; 
      }
    }
  }

  else {// shock (p<=ps, rho<=rhos)

    HugoniotEquation hugo(vf[id],rho,p,ps);

    //find a bracketing interval
    double f0, f1;
    double drho = std::max(fabs(rhos0 - rhos1), 0.001*rhos0); 
    bool found_rhos0 = false, found_rhos1 = false;
    if(std::min(rhos0,rhos1)>=rho) {//both rhos0 and rhos1 are physically admissible
      f0 = hugo(rhos0);
      f1 = hugo(rhos1);
      if(f0*f1<=0) {
        /*found bracketing interval [rhos0, rhos1]*/
        if(rhos0>rhos1) {
          std::swap(rhos0,rhos1);
          std::swap(f0,f1);
        }
        found_rhos0 = found_rhos1 = true;
      } else {
        rhos0 = rhos1; //this is our starting point (presumably rhos1 is closer to sol'n)
        f0    = f1; 
      }
    } else { //at least, the smaller one among rhos0, rhos1 is non-physical
      if(rhos1>rhos0) {
        rhos0 = rhos1; 
        f0    = hugo(rhos0);
      }
      if(rhos0<rho) {//this one is also non-physical
        rhos0 = rho;
        f0    = hugo(rhos0);
        found_rhos0 = true;
      } else {
        /*rhos0 = rhos0;*/
        f0 = hugo(rhos0);
      }
    } 

    if(!found_rhos0 || !found_rhos1) {
      int i = 0;
      double factor = 1.5; 
      double tmp, ftmp;
      // before the search, rhos0 = rhos1 = an adimissible point > rho
      rhos1 = rhos0;
      f1    = f0;
      while(!found_rhos0) {
        if(++i>=maxIts_shock) {
//          cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
//               << " iterations (in the solution of the Hugoniot equation)." << endl;
          return false;
        }
        tmp = rhos1;
        ftmp = f1;
        //move to the left (towards rho)
        rhos1 = rhos0;
        f1    = f0;
        rhos0 = rhos1 - factor*drho;
        if(rhos0<=rho) {
          rhos0 = rho;
          found_rhos0 = true;
        }
        f0 = hugo(rhos0);

        if(f0*f1<=0) {
          found_rhos0 = found_rhos1 = true;
        } else {
          //move to the right
          rhos1 = tmp;
          f1    = ftmp;
          tmp   = rhos0; //don't forget the smallest point
          ftmp  = f0;
          rhos0 = rhos1;
          f0    = f1;
          rhos1 = rhos0 + factor*drho;
          f1 = hugo(rhos1);
          if(f0*f1<=0) {
            found_rhos0 = found_rhos1 = true;
          } else {
            rhos0 = tmp;
            f0    = ftmp;
            drho  = rhos1 - rhos0; //update drho
          }
        }
      }

      if(!found_rhos1) {//keep moving to the right
        i = 0;
        double factor = 2.5;
        while(!found_rhos1) {
          if(++i>=maxIts_shock) {
//            cout << "*** Error: Unable to find a bracketing interval after " << maxIts_shock 
//                 << " iterations (in the solution of the Hugoniot equation (2))." << endl;
            return false; //failure
          }
          rhos0 = rhos1;
          f0    = f1;
          rhos1 = rhos0 + factor*drho;
          f1    = hugo(rhos1);
          if(f0*f1<=0) {
            found_rhos1 = true;
          } else
            drho = rhos1 - rhos0;
        }
      }

    }

    pair<double,double> sol;

#ifndef WITHOUT_BOOST
    //*******************************************************************
    // Calling boost function for root-finding
    // Warning: "maxit" is BOTH AN INPUT AND AN OUTPUT
    boost::uintmax_t maxit = maxIts_shock;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      sol = toms748_solve(hugo, rhos0, rhos1, f0, f1,
                          [=](double r0, double r1){return r1-r0<std::min(tol_shock,0.001*(rhos1-rhos0));}, 
                          maxit);
    }
    //*******************************************************************
#else
    //*******************************************************************
    // Using a hybrid (Brent) method for root-finding
    int maxit = 0;
    if(f0==0.0)
      sol.first = sol.second = rhos0;
    else if(f1==0.0)
      sol.first = sol.second = rhos1;
    else {
      double rhos2 = rhos1; //rhos2 is always the latest one
      double f2    = f1;
      int it;
      for(it = 0; it<maxIts_shock; it++) {
        drho = rhos1 - rhos0;
        rhos2 = rhos2 - f2*(rhos1 - rhos0)/(f1 - f0); //secant method
        if(rhos2 >= rhos1 || rhos2 <= rhos0) //discard and switch to bisection
          rhos2 = 0.5*(rhos0+rhos1);
        f2 = hugo(rhos2);
        if(f2==0.0) {
          sol.first = sol.second = rhos2;
          break;
        }
        if(f2*f0<0) {
          rhos1 = rhos2;
          f1    = f2;
        } else {
          rhos0 = rhos2;
          f0    = f2;
        }
        if(rhos1-rhos0<tol_shock) {
          sol.first  = rhos0;
          sol.second = rhos1;
          break;
        }
      }
      if(it==maxIts_shock) {
        fprintf(stderr,"*** Error: Root-finding method failed to converge after %d iterations.\n", it);
        return false;
      }
      maxit = it;
    }
    //*******************************************************************
#endif
    

#if PRINT_RIEMANN_SOLUTION == 1
    cout << "  " << wavenumber << "-wave: shock, converged in " << maxit << " iterations. fun = " 
         << hugo(0.5*(sol.first+sol.second)) << "." << endl;
#endif

    rhos = 0.5*(sol.first+sol.second);

    double du = -(ps-p)*(1.0/rhos-1.0/rho);
    if(du<0) {
      //cout << "Warning: Violation of hyperbolicitiy when enforcing the Rankine-Hugoniot jump conditions (du = "
      //     << du << ")." << endl;
      return false;
    }

    if(vf[id]->CheckState(rhos,ps))
      return false; //nonphysical...

    us = (wavenumber==1) ? u - sqrt(du) : u + sqrt(du);
    

#if PRINT_RIEMANN_SOLUTION == 1
    double xi = (rhos*us - rho*u)/(rhos-rho);
    if(wavenumber==1) {
      sol1d.push_back(vector<double>{xi-0.0001*fabs(xi), rho, u, p, (double)id});
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
    } else {
      sol1d.push_back(vector<double>{xi, rhos, us, ps, (double)id});
      sol1d.push_back(vector<double>{xi+0.0001*fabs(xi), rho, u, p, (double)id});
    }
#endif

  }

  return true; //yeah!
}

//----------------------------------------------------------------------------------

bool
ExactRiemannSolverBase::Rarefaction_OneStepRK4(int wavenumber/*1 or 3*/, int id,
                            double rho_0, double u_0, double p_0 /*start state*/, 
                            double drho /*step size*/,
                            double &rho, double &u, double &p, double &xi /*output*/)
{
  drho = -drho; //drho is positive when passed in. It is actually negative if we follow Kamm's paper
                // Equations (36 - 42)
                
  double e_0 = vf[id]->GetInternalEnergyPerUnitMass(rho_0, p_0);
  double c_0_square = vf[id]->ComputeSoundSpeedSquare(rho_0, e_0);

  if(rho_0<=0 || c_0_square<0) {
//    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(0)." 
//            " rho = %e, p = %e, e = %e, id = %d.\n",
//            c_0_square, rho_0, p_0, e_0, id);
    return false;
  } 

  double c_0 = sqrt(c_0_square);

  double p_1 = p_0 + 0.5*drho*c_0_square;
  double rho_1 = rho_0 + 0.5*drho;
  double e_1 = vf[id]->GetInternalEnergyPerUnitMass(rho_1, p_1);
  double c_1_square = vf[id]->ComputeSoundSpeedSquare(rho_1, e_1);

  if(rho_1<=0 || c_1_square<0) {
//    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(1)." 
//                   " rho = %e, p = %e, e = %e, id = %d.\n",
//                   c_1_square, rho_1, p_1, e_1, id);
    return false;
  } 

  double c_1 = sqrt(c_1_square);

  double p_2 = p_0 + 0.5*drho*c_1_square;
  double rho_2 = rho_1;
  double e_2 = vf[id]->GetInternalEnergyPerUnitMass(rho_2, p_2);
  double c_2_square = vf[id]->ComputeSoundSpeedSquare(rho_2, e_2);

  if(rho_2<=0 || c_2_square<0) {
//    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(2)." 
//            " rho = %e, p = %e, e = %e, id = %d.\n",
//            c_2_square, rho_2, p_2, e_2, id);
    return false;
  } 

  double c_2 = sqrt(c_2_square);
 
  double p_3 = p_0 + drho*c_2_square;
  double rho_3 = rho_0 + drho;
  double e_3 = vf[id]->GetInternalEnergyPerUnitMass(rho_3, p_3);
  double c_3_square = vf[id]->ComputeSoundSpeedSquare(rho_3, e_3);

  if(rho_3<=0 || c_3_square<0) {
//    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(3)." 
//            " rho = %e, p = %e, e = %e, id = %d.\n",
//            c_3_square, rho_3, p_3, e_3, id);
    return false;
  } 

  double c_3 = sqrt(c_3_square);

  // now, calculate the outputs
  //
  p = p_0 + 1.0/6.0*drho*(c_0_square + 2.0*(c_1_square+c_2_square) + c_3_square);

  double du = 1.0/6.0*drho*( c_0/rho_0 + 2.0*(c_1/rho_1 + c_2/rho_2) + c_3/rho_3);
  u = (wavenumber == 1) ? u_0 - du : u_0 + du;

  rho = rho_0 + drho;  //KW: Note that if drho is tiny compared to rho_0, then rho = rho_0 due to roundoff (I have observed this!)

  double e = vf[id]->GetInternalEnergyPerUnitMass(rho, p);
  double c = vf[id]->ComputeSoundSpeedSquare(rho, e);

  if(rho<=0 || c<0) {
//    fprintf(stderr,"*** Error: Negative density or c^2 (square of sound speed, %e) in Rarefaction_OneStepRK4(final)." 
//            " rho = %e, p = %e, e = %e, id = %d.\n",
//            c, rho, p, e, id);
    return false;
  } else
    c = sqrt(c);

  xi = (wavenumber == 1) ? u - c : u + c;

  return true;
}

//----------------------------------------------------------------------------------

void
ExactRiemannSolverBase::PrintStarRelations(double rhol, double ul, double pl, int idl,
                                           double rhor, double ur, double pr, int idr,
                                           double pmin, double pmax, double dp)
{

  vector<std::array<double,3> > left; //(p*, rhol*, ul*)
  vector<std::array<double,3> > right; //(p*, rhor*, ur*)

  double ps, rhols, rhors, uls,  urs;
  bool success;
  ps = pmin;

  while(true) {

    success = ComputeRhoUStar(1, rhol, ul, pl, ps, idl/*inputs*/,
                  rhol, (ps>pl) ? rhol*1.1 : rhol*0.9/*initial guesses for Hugo. eq.*/,
                  rhols, uls/*outputs*/);
    if(success)
      left.push_back(std::array<double,3>{{ps,rhols,uls}});
    else 
      fprintf(stderr," -- ComputeRhoUStar(1) failed. left state: %e %e %e (%d), ps = %e.\n",
              rhol, ul, pl, idl, ps);
    
    success = ComputeRhoUStar(3, rhor, ur, pr, ps, idr/*inputs*/,
                  rhor, (ps>pr) ? rhor*1.1 : rhor*0.9/*initial guesses for Hugo. eq.*/,
                  rhors, urs/*outputs*/);
    if(success)
      right.push_back(std::array<double,3>{{ps,rhors,urs}});
    else
      fprintf(stderr," -- ComputeRhoUStar(3) failed. right state: %e %e %e (%d), ps = %e.\n",
              rhor, ur, pr, idr, ps);

    if(ps>=pmax)
      break;

    ps = std::min(ps+dp, pmax);

  }
   
  
  FILE* file = fopen("LeftStarState.txt", "w");
  print(file, "## One-Dimensional Riemann Problem.\n");
  print(file, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
        rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(file, "## pmin = %e, pmax = %e, dp = %e.\n", pmin, pmax, dp);

  for(auto it = left.begin(); it != left.end(); it++) 
    print(file, "%e    %e    %e\n", (*it)[0], (*it)[1], (*it)[2]);
  fclose(file);

 
  file = fopen("RightStarState.txt", "w");
  print(file, "## One-Dimensional Riemann Problem.\n");
  print(file, "## Initial State: %e %e %e, id %d (left) | (right) %e %e %e, id %d.\n", 
        rhol, ul, pl, idl, rhor, ur, pr, idr);
  print(file, "## pmin = %e, pmax = %e, dp = %e.\n", pmin, pmax, dp);

  for(auto it = right.begin(); it != right.end(); it++) 
    print(file, "%e    %e    %e\n", (*it)[0], (*it)[1], (*it)[2]);
  fclose(file);


}

//----------------------------------------------------------------------------------


//----------------------------------------------------------------------------------
