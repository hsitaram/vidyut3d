#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_molwt[NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_charge[NUM_ALL_SPECIES]={0.0};

    void init()
    {
        specnames[N_ID]="Neutral";
        specnames[I_ID]="Ion";

        spec_charge[N_ID]   =  0.0;
        spec_charge[I_ID]   =  1.0;
        spec_charge[EDN_ID] = -1.0;

        spec_molwt[N_ID]   = 4.0*M_AMU;
        spec_molwt[I_ID]   = 4.0*M_AMU;
        spec_molwt[EDN_ID] = ME;
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
    amrex::Real get_bg_molwt(Real specden[NUM_ALL_SPECIES])
    {
        return(spec_molwt[N_ID]);
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
    void get_wdot(amrex::Real Te,
                  amrex::Real Tg, amrex::Real Pg,
                  amrex::Real efield, 
                  amrex::Real spec[NUM_ALL_SPECIES],
                  amrex::Real spec_wdot[NUM_ALL_SPECIES+1])
    {
        for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
        {
            spec_wdot[sp]=0.0;
        }

        amrex::Real rateconst;
        
        //E + N --> I + 2E
        rateconst = 2.5e-12
        *std::exp(-278508.0/Te);
        spec_wdot[EDN_ID] += rateconst*spec[EDN_ID]*spec[N_ID];
        spec_wdot[I_ID]   += rateconst*spec[EDN_ID]*spec[N_ID];
        spec_wdot[EEN_ID] += -24.0*ECHARGE*rateconst*spec[EDN_ID]*spec[N_ID];
    }

    AMREX_GPU_DEVICE 
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg)

    {
        //mobility is scaled by charge
        amrex::Real mob=0.0;
        amrex::Real elecmob=-20.0;
        if(specid==EDN_ID)
        {
            mob=elecmob;
        }
        if(specid==I_ID)
        {
            mob=0.2;
            //mob=0.0;
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
        amrex::Real dcoeff=0.0;
        amrex::Real elecdcoeff=100.0;
        if(specid==EDN_ID)
        {
            dcoeff=elecdcoeff;
        }
        if(specid==I_ID)
        {
            dcoeff=0.01;
            //dcoeff=0.0;
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
