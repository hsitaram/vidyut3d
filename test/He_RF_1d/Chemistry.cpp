#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_molwt[NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_charge[NUM_ALL_SPECIES]={0.0};

    void init()
    {
        specnames[HE_ID]="He";
        specnames[HEp_ID]="He+";
        specnames[HE2p_ID]="He2+";
        specnames[HEm_ID]="Hem";
        specnames[HE2m_ID]="He2m";
        
        spec_charge[HE_ID]   =  0.0;
        spec_charge[HEp_ID]  =  1.0;
        spec_charge[HE2p_ID] =  1.0;
        spec_charge[HEm_ID]  =  0.0;
        spec_charge[HE2m_ID] =  0.0;
        spec_charge[EDN_ID]  = -1.0;

        spec_molwt[HE_ID]   = 4.0*M_AMU;
        spec_molwt[HEp_ID]  = 4.0*M_AMU;
        spec_molwt[HE2p_ID] = 8.0*M_AMU;
        spec_molwt[HEm_ID]  = 4.0*M_AMU;
        spec_molwt[HE2m_ID] = 8.0*M_AMU;
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
        return(spec_molwt[HE_ID]);
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

        //E + HE --> HEm + E
        rate = 2.308e-10*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.31)*std::exp(-2.297e5/Te) * spec[EDN_ID] * spec[HE_ID];
        spec_wdot[HEm_ID] += rate;
        spec_wdot[EEN_ID] += -19.8*ECHARGE*rate;

        //E + HEm --> HE + E
        rate = 1.099e-11*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.31) * spec[EDN_ID] * spec[HEm_ID];
        spec_wdot[HEm_ID] -= rate;
        spec_wdot[EEN_ID] += 19.8*ECHARGE*rate;

        //E + HE --> HE+ + 2E
        rate = 2.584e-12*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.68)*std::exp(-2.854092e5/Te) * spec[EDN_ID] * spec[HE_ID];
        spec_wdot[EDN_ID] += rate;
        spec_wdot[HEp_ID] += rate;
        spec_wdot[EEN_ID] += -24.6*ECHARGE*rate;

        //E + HEm --> HE+ + 2E
        rate = 4.661e-10*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.6)*std::exp(-5.546e4/Te) * spec[EDN_ID] * spec[HEm_ID];
        spec_wdot[EDN_ID] += rate;
        spec_wdot[HEm_ID] -= rate;
        spec_wdot[HEp_ID] += rate;
        spec_wdot[EEN_ID] += -4.78*ECHARGE*rate;

        //E + HE2m --> HE2p + 2E
        rate = 1.268e-12*(std::pow(cm_to_m,3.0))
        *std::pow(Te,0.71)*std::exp(-3.945e4/Te) * spec[EDN_ID] * spec[HE2m_ID];
        spec_wdot[EDN_ID]  += rate;
        spec_wdot[HE2m_ID] -= rate;
        spec_wdot[HE2p_ID] += rate;
        spec_wdot[EEN_ID]  += -3.4*ECHARGE*rate; 

        //E + HE2p --> HEm + HE
        rate = 5.386e-7*(std::pow(cm_to_m,3.0))
        *std::pow(Te,-0.5) * spec[EDN_ID] *  spec[HE2p_ID];
        spec_wdot[EDN_ID]  -= rate;
        spec_wdot[HE2p_ID] -= rate;
        spec_wdot[HEm_ID]  += rate;
        spec_wdot[EEN_ID]  += 0.0; 

        //2 HEm --> HE+ + HE + E
        rate = 2.7e-10*(std::pow(cm_to_m,3.0)) * std::pow(spec[HEm_ID],2.0);
        spec_wdot[HEm_ID]  -= 2.0*rate;
        spec_wdot[HEp_ID]  += rate;
        spec_wdot[EDN_ID]  += rate;
        spec_wdot[EEN_ID]  += 15.0*ECHARGE*rate;

        //HEm + 2HE --> HE2m + HE
        rate = 1.3e-33*(std::pow(cm_to_m,3.0)) * spec[HEm_ID] * std::pow(spec[HE_ID],2.0);
        spec_wdot[HEm_ID]  -= rate;
        spec_wdot[HE2m_ID] += rate;
        spec_wdot[EEN_ID]  += 0.0;

        //HE+ + 2HE --> HE2+ + HE
        rate = 1.e-31*(std::pow(cm_to_m,3.0)) * spec[HEp_ID] * std::pow(spec[HE_ID],2.0);
        spec_wdot[HEp_ID]  -= rate;
        spec_wdot[HE2p_ID] += rate;
        spec_wdot[EEN_ID]  += 0.0;

    }
    
    AMREX_GPU_DEVICE 
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg)

    {
        amrex::Real mob=0.0;
        //mobility is scaled by charge
        
        //reduced mobility obtained from Yuan and Raja,
        //IEEE. Trans. Plasma Sci.,31,4,2003
        //they say these are 393K, may be temperature
        //scaling is necessary
        //amrex::Real Pres_ratio = (P_NTP/Pg);
        //amrex::Real elecmob=-0.1132*Pres_ratio;

        amrex::Real N=Pg/K_B/Tg;
        amrex::Real EbyN=efield/N/TDUNIT;
        amrex::Real meanE = 1.5*Te/eV;

        //mobility from 
        //Turner, Miles M., et al. "Simulation benchmarks for low-pressure
        //plasmas: Capacitive discharges." 
        //Physics of Plasmas 20.1 (2013): 013507.
        amrex::Real Hepmob=2.69*std::pow((1.0+1.2e-3*std::pow(EbyN,2.0)+4.2e-8*std::pow(EbyN,4.0)),-0.125);
        amrex::Real He2pmob=Hepmob;
        
        //computed by Taaresh Taneja (U Minnesota) using Turner's cross sections
        amrex::Real elecmob = (-1.0)*std::exp(55.0 + 0.3942*std::log(meanE) + 2.134/meanE 
                   -0.6433/std::pow(meanE,2.0) + (0.7112e-1)/std::pow(meanE,3.0)) / N;
        
        if(specid==EDN_ID)
        {
            mob=elecmob;
        }
        if(specid==HEp_ID)
        {
            //mob=0.001482*Pres_ratio;
            mob=Hepmob;
        }
        if(specid==HE2p_ID)
        {
            //mob=0.002403*Pres_ratio;
            mob=He2pmob;
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
        amrex::Real Pres_ratio = (P_NTP/Pg);
        amrex::Real elecdcoeff=0.1737*(Te/17406.0)*Pres_ratio;
        amrex::Real dcoeff=0.0;
        if(specid==EDN_ID)
        {
            dcoeff=elecdcoeff;
        }
        if(specid==HEp_ID)
        {
            dcoeff=5.026e-5*Pres_ratio;
            //dcoeff=0.0;
        }
	if(specid==HE2p_ID)
        {
           dcoeff=8.148e-5*Pres_ratio;
        }
	if(specid==HEm_ID)
        {
           dcoeff=4.116e-4*Pres_ratio;
        }
        if(specid==HE2m_ID)
        {
	  dcoeff=2.029e-4*Pres_ratio;
        }
        if(specid==EEN_ID)
        {
            dcoeff=fivebythree*elecdcoeff;
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
