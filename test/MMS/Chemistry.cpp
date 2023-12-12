#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_molwt[NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real spec_charge[NUM_ALL_SPECIES]={0.0};

    void init()
    {
        specnames[S1_ID]="S1";
        
        spec_charge[S1_ID]   =  0.0;
        spec_charge[EDN_ID]  = -1.0;

        spec_molwt[S1_ID]   = M_AMU;
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
    amrex::Real get_bg_molwt(amrex::Real specden[NUM_ALL_SPECIES])
    {
        amrex::Real molwt=M_AMU;    
        return(molwt);
    }
    
    AMREX_GPU_HOST_DEVICE  
    void get_wdot(amrex::Real Te,
                  amrex::Real Tg, amrex::Real Pg,
                  amrex::Real efield_mag, 
                  amrex::Real spec[NUM_ALL_SPECIES],
                  amrex::Real spec_wdot[NUM_ALL_SPECIES+1])
    {
        for(int sp=0;sp<(NUM_ALL_SPECIES+1);sp++)
        {
            spec_wdot[sp]=0.0;
        }
    }
    
    AMREX_GPU_DEVICE 
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg)

    {
        //mobility is scaled by charge
        amrex::Real mob=-1.0;
        return(mob);
    }
    
    AMREX_GPU_DEVICE 
    amrex::Real diffusion_coeff(int specid,
                                amrex::Real Te,
                                amrex::Real efield,
                                amrex::Real Tg, amrex::Real Pg)

    {
        amrex::Real dcoeff=0.0;
        if(specid==S1_ID)
        {
           dcoeff = 1.0;
        }
        return(dcoeff);
    }

    AMREX_GPU_DEVICE 
    amrex::Real electron_collision_freq(
                               amrex::Real Te,
                               amrex::Real efield,
                               Real Tg, Real Pg)

    {
        return(1.0);
    }
}
