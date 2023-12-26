#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_molwt[NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_charge[NUM_ALL_SPECIES]={0.0};

    void init()
    {
        specnames[AR_ID]="Ar";
        specnames[ARp_ID]="Ar+";
        specnames[ARm_ID]="Arm";
        
        spec_charge[AR_ID]   =  0.0;
        spec_charge[ARp_ID]  =  1.0;
        spec_charge[ARm_ID]  =  0.0;
        spec_charge[EDN_ID]  = -1.0;

        spec_molwt[AR_ID]   = 39.948*M_AMU;
        spec_molwt[ARp_ID]  = 39.948*M_AMU;
        spec_molwt[ARm_ID]  = 39.948*M_AMU;
        spec_molwt[EDN_ID]  = ME;
    }    
    
    void close()
    {
        specnames.clear();
    }
    
    int find_id(std::string specname)
    {
        int loc=-1;
        auto it=std::find(specnames.begin(),specnames.end(),specname);
        if(it != specnames.end())
        {
            loc=it-specnames.begin();
        }
        return(loc);
    }
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_charge(int specid)
    {
        return(spec_charge[specid]);
    }
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_molwt(int specid)
    {
        return(spec_molwt[specid]);
    }
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_bg_molwt(Real specden[NUM_ALL_SPECIES])
    {
        return(spec_molwt[AR_ID]);
    }
   
    AMREX_GPU_HOST_DEVICE void get_wdot(amrex::Real Te,amrex::Real Tg, amrex::Real Pg,
                                        amrex::Real efield, amrex::Real spec[NUM_ALL_SPECIES],
                                        amrex::Real spec_wdot[NUM_ALL_SPECIES+1])
    {
        const amrex::Real cm_to_m=0.01;        
        amrex::Real rate;
        
        for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
        {
            spec_wdot[sp]=0.0;
        }

        //E + AR --> ARm + E
        rate = 1.1748e-8*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.046639)*std::exp(-1.3856e5/Te) * spec[EDN_ID] * spec[AR_ID];
        spec_wdot[ARm_ID] += rate;
        spec_wdot[EEN_ID] += -11.56*ECHARGE*rate;


        //E + AR --> AR+ + 2E
        rate = 7.0708e-11*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.60983)*std::exp(-1.8712e5/Te) * spec[EDN_ID] * spec[AR_ID];
        spec_wdot[EDN_ID] += rate;
        spec_wdot[ARp_ID] += rate;
        spec_wdot[EEN_ID] += -16.0*ECHARGE*rate;

        //E + ARm --> AR+ + 2E
        rate = 1.2456e-7*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.050382)*std::exp(-6.0524e4/Te) * spec[EDN_ID] * spec[ARm_ID];
        spec_wdot[EDN_ID] += rate;
        spec_wdot[ARm_ID] -= rate;
        spec_wdot[ARp_ID] += rate;
        spec_wdot[EEN_ID] += -4.43*ECHARGE*rate;
        
        //2 ARm --> AR+ + AR + E
        rate = 6.2e-10*(std::pow(cm_to_m,3.0)) * std::pow(spec[ARm_ID],2.0);
        spec_wdot[ARm_ID]  -= 2.0*rate;
        spec_wdot[ARp_ID]  += rate;
        spec_wdot[EDN_ID]  += rate;
        spec_wdot[EEN_ID]  += 7.0*ECHARGE*rate;
        
        //E + ARm --> AR + E
        rate = 2.0e-7*(std::pow(cm_to_m,3.0)) * spec[EDN_ID] * spec[ARm_ID];
        spec_wdot[ARm_ID] -= rate;
        spec_wdot[EEN_ID] += 11.56*ECHARGE*rate;

        //ARm + AR --> 2AR
        rate = 2.5e-15*(std::pow(cm_to_m,3.0)) * spec[ARm_ID] * spec[AR_ID];
        spec_wdot[ARm_ID]  -= rate;
        spec_wdot[EEN_ID]  += 0.0;

    }
    
    AMREX_GPU_DEVICE 
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg)

    {
        amrex::Real Arp_cs = 120e-20; //A^2 (cross section)
        amrex::Real E_cs = 40e-20; //A^2 (cross section)
        amrex::Real Ng = Pg/K_B/Tg;
        amrex::Real nu_elec = std::sqrt(8.0*K_B*Te/PI/ME)*Ng*E_cs;
        amrex::Real elecmob = -ECHARGE/ME/nu_elec;

        amrex::Real mob=0.0;
        if(specid==EDN_ID)
        {
            mob=elecmob;
        }
        if(specid==ARp_ID)
        {
            amrex::Real nu_Arp=std::sqrt(8.0*K_B*Tg/PI/spec_molwt[ARp_ID])*Ng*Arp_cs;
            mob=ECHARGE/spec_molwt[ARp_ID]/nu_Arp;
        }
        if(specid==EEN_ID)
        {
            mob=fivebythree*elecmob;
        }

        return(mob);
    }
    
    AMREX_GPU_DEVICE 
    amrex::Real diffusion_coeff(int specid,
                                amrex::Real Te,
                                amrex::Real efield,
                                amrex::Real Tg, amrex::Real Pg)

    {

        //mu=e/(m nu) , D=kT/(m nu)
        // m nu = e/mu  D=kT/(e/mu) = mu * T/(e/k) 
        // = mu * T_in_ev

        amrex::Real elecmob=mobility(EDN_ID,Te,efield,Tg,Pg);
        amrex::Real ionmob=mobility(ARp_ID,Te,efield,Tg,Pg);
        amrex::Real Te_in_ev=Te/eV;
        amrex::Real Tg_in_ev=Tg/eV;
        amrex::Real dcoeff=0.0;

        if(specid==EDN_ID)
        {
            dcoeff=Te_in_ev*amrex::Math::abs(elecmob);
        }
        if(specid==ARp_ID)
        {
            dcoeff=Tg_in_ev*ionmob;
        }
	if(specid==ARm_ID)
        {
           dcoeff=Tg_in_ev*ionmob;
        }
        if(specid==AR_ID)
        {
           dcoeff=Tg_in_ev*ionmob;
        }
        if(specid==EEN_ID)
        {
            dcoeff=fivebythree*Te_in_ev*amrex::Math::abs(elecmob);
        }

        return(dcoeff);
    }

    AMREX_GPU_DEVICE 
    amrex::Real electron_collision_freq(
                               amrex::Real Te,
                               amrex::Real efield,
                               Real Tg, Real Pg)

    {
        amrex::Real mu = mobility(EDN_ID, Te, efield, Tg, Pg);
        amrex::Real collfreq = -ECHARGE/ME/mu;
        return(collfreq);
    }
}
