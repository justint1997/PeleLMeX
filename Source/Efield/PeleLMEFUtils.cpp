#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMEF_K.H>

using namespace amrex;

Vector<Array<MultiFab*,AMREX_SPACEDIM>>
PeleLM::getNLgradPhiVVect() {
   Vector<Array<MultiFab*,AMREX_SPACEDIM>> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(GetArrOfPtrs(m_leveldatanlsolve[lev]->gPhiVOld));
   }
   return r;
}

Vector<Array<MultiFab*,AMREX_SPACEDIM>>
PeleLM::getUeffVect() {
   Vector<Array<MultiFab*,AMREX_SPACEDIM>> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(GetArrOfPtrs(m_leveldatanlsolve[lev]->uEffnE));
   }
   return r;
}

Vector<MultiFab*>
PeleLM::getNLresidVect() {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldatanlsolve[lev]->nlResid));
   }
   return r;
}

Vector<MultiFab*>
PeleLM::getNLstateVect() {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldatanlsolve[lev]->nlState));
   }
   return r;
}

void PeleLM::getNLStateScaling(Real &nEScale, Real &phiVScale)
{
   Array<Real,2> r = {0.0,0.0};
   for (int comp = 0; comp < 2; comp++) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (lev != finest_level) {
            r[comp] = std::max(r[comp],
                               m_leveldatanlsolve[lev]->nlState.norm0(*m_coveredMask[lev],comp,0,true));
         } else {
            r[comp] = std::max(r[comp],
                               m_leveldatanlsolve[lev]->nlState.norm0(comp,0,true,true));
         }
      }
      ParallelDescriptor::ReduceRealMax(r[comp]);
   }
   nEScale = r[0];
   phiVScale = r[1];
}

void PeleLM::getNLResidScaling(Real &nEScale, Real &phiVScale)
{
   Array<Real,2> r = {0.0,0.0};
   for (int comp = 0; comp < 2; comp++) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (lev != finest_level) {
            r[comp] = std::max(r[comp], 
                               m_leveldatanlsolve[lev]->nlResid.norm0(*m_coveredMask[lev],comp,0,true));
         } else {
            r[comp] = std::max(r[comp], 
                               m_leveldatanlsolve[lev]->nlResid.norm0(comp,0,true));
         }
      }
      ParallelDescriptor::ReduceRealMax(r[comp]);
   }
   nEScale = r[0];
   phiVScale = r[1];
}

void PeleLM::scaleNLState(const Real &nEScale, const Real &phiVScale)
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      m_leveldatanlsolve[lev]->nlState.mult(1.0/nE_scale,0,1,1);
      m_leveldatanlsolve[lev]->nlState.mult(1.0/phiV_scale,1,1,1);
   }
}

void PeleLM::scaleNLResid(const Vector<MultiFab*> &a_resid, const Real &nEScale, const Real &phiVScale)
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      a_resid[lev]->mult(1.0/FnE_scale,0,1,1);
      a_resid[lev]->mult(1.0/FphiV_scale,1,1,1);
   }
}

BCRec
PeleLM::hackBCChargedParticle(const Real &charge,
                              const BCRec &bc_in) {

   BCRec bc_hacked;

   const int* lo_bc = bc_in.lo();
   const int* hi_bc = bc_in.hi();

   for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
   {

      int lo = lo_bc[idim];
      int hi = hi_bc[idim];
      // Spec is In/Out and it's cathode (neg electrode)
      if ( ( lo_bc[idim] == amrex::BCType::ext_dir ||
             lo_bc[idim] == amrex::BCType::foextrap ) &&
           ( m_phiV_bcpol.lo(idim) == 2 ) ) {
         if ( charge > 0.0 ) { // Outflow for cation
            lo = amrex::BCType::foextrap;
         } else {            // Dirich = 0 for anion
            lo = amrex::BCType::ext_dir;
         }
      } else if ( ( lo_bc[idim] == amrex::BCType::ext_dir ||
                    lo_bc[idim] == amrex::BCType::foextrap ) &&
                  ( m_phiV_bcpol.lo(idim) == 1 ) ) {
         if ( charge > 0.0 ) { // Dirich = 0 for cation
            lo = amrex::BCType::ext_dir;
         } else {            // Outflow for anion
            lo = amrex::BCType::foextrap;
         }
      }
      if ( ( hi_bc[idim] == amrex::BCType::ext_dir ||
             hi_bc[idim] == amrex::BCType::foextrap ) &&
           ( m_phiV_bcpol.hi(idim) == 2 ) ) {
         if ( charge > 0.0 ) { // Outflow for cation
            hi = amrex::BCType::foextrap;
         } else {            // Dirich = 0 for anion
            hi = amrex::BCType::ext_dir;
         }
      } else if ( ( hi_bc[idim] == amrex::BCType::ext_dir ||
                    hi_bc[idim] == amrex::BCType::foextrap ) &&
                  ( m_phiV_bcpol.hi(idim) == 1 ) ) {
         if ( charge > 0.0 ) { // Dirich = 0 for cation
            hi = amrex::BCType::ext_dir;
         } else {            // Outflow for anion
            hi = amrex::BCType::foextrap;
         }
      }
      bc_hacked.setLo(idim,lo);
      bc_hacked.setHi(idim,hi);
   }
   return bc_hacked;
}


void PeleLM::addLorentzVelForces(int lev,
                                 const Box&       bx,
                                 const Real&      a_time,
                                 Array4<      Real> const& force,
                                 Array4<const Real> const& rhoY,
                                 Array4<const Real> const& phiV,
                                 Array4<const Real> const& nE)
{
   const auto  dx       = geom[lev].CellSizeArray();
   GpuArray<int,3> blo = bx.loVect3d();
   GpuArray<int,3> bhi = bx.hiVect3d();

   amrex::ParallelFor(bx, [force, rhoY, phiV, nE, a_time, dx, blo, bhi, zk=zk]
   AMREX_GPU_DEVICE(int i, int j, int k) noexcept
   {
      addLorentzForce(i,j,k, blo, bhi, a_time, dx, zk, rhoY, nE, phiV, force);
   });
}

void PeleLM::initializeElectronNeutral()
{
   // Prob/PMF datas
   ProbParm const* lprobparm = prob_parm.get();

   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get level data new time pointer
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->species, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho      = ldata_p->density.array(mfi);
         auto const& rhoY     = ldata_p->species.array(mfi);
         auto const& rhoH     = ldata_p->rhoh.array(mfi);
         auto const& temp     = ldata_p->temp.array(mfi);
         auto const& nE       = ldata_p->nE.array(mfi);
         amrex::ParallelFor(bx, [rho, rhoY, rhoH, temp, nE, lprobparm]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            initElecNeutral(i,j,k,rho,rhoY,rhoH,temp,nE,*lprobparm);
         });
      }

      // Convert I_R(Y_nE) into I_R(nE) and set I_R(Y_nE) to zero
      auto ldataR_p   = getLevelDataReactPtr(lev);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataR_p->I_R, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& YnEdot = ldataR_p->I_R.array(mfi,E_ID);
         auto const& nEdot  = ldataR_p->I_R.array(mfi,NUM_SPECIES);
         auto eos = pele::physics::PhysicsType::eos();
         Real invmwt[NUM_SPECIES] = {0.0};
         eos.inv_molecular_weight(invmwt);
         ParallelFor(bx, [YnEdot,nEdot,invmwt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            nEdot(i,j,k) = YnEdot(i,j,k) * Na * invmwt[E_ID] * 1.0e3;
            YnEdot(i,j,k) = 0.0;
         });
      }
   }
}

void PeleLM::initializeElectronFromMassFraction()
{
}
