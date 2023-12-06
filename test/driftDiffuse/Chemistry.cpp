#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);
    AMREX_GPU_DEVICE_MANAGED amrex::Real rct_rxnarray[NUM_REACTIONS][NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real pdt_rxnarray[NUM_REACTIONS][NUM_ALL_SPECIES]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real reaction_elecenergy[NUM_REACTIONS]={0.0};
    AMREX_GPU_DEVICE_MANAGED amrex::Real reaction_gasenergy[NUM_REACTIONS]={0.0};

    void init()
    {
        specnames[S1_ID]="S1";
        specnames[S2_ID]="S2";
        
        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_ALL_SPECIES;j++)
            {
                rct_rxnarray[i][j]=0.0;
                pdt_rxnarray[i][j]=0.0;
            }
        }

        int rnum=0;
    
        //just a dummy reaction    
        //E + S1 --> S2 + E
        rct_rxnarray[rnum][EDN_ID] = 1.0;
        rct_rxnarray[rnum][S1_ID]  = 1.0;

        pdt_rxnarray[rnum][S2_ID]  = 1.0;
        pdt_rxnarray[rnum][EDN_ID] = 1.0;
        
        reaction_elecenergy[rnum]  = 0.0;
        reaction_gasenergy[rnum]   = 0.0;
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
        amrex::Real charge=0.0;
        switch(specid)
        {
            case S1_ID:
                charge=0.0;
                break;               
            case S2_ID:
                charge=0.0;
                break;               
            case EDN_ID:
                charge=-1.0;
                break;               
            default:
                charge=0.0;
        }
        return(charge);
    }
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_molwt(int specid)
    {
        amrex::Real molwt=M_AMU;
        return(molwt);
    }
    
    AMREX_GPU_HOST_DEVICE 
    amrex::Real get_bg_molwt(amrex::Real specden[NUM_ALL_SPECIES])
    {
        amrex::Real molwt=M_AMU;    
        return(molwt);
    }
    
    AMREX_GPU_HOST_DEVICE  
    void get_reaction_rateconstants(amrex::Real Te,
                                    amrex::Real Tg, amrex::Real Pg,
                                     amrex::Real efield_mag, 
                                     amrex::Real spec[NUM_ALL_SPECIES],
                                    amrex::Real rateconst[NUM_REACTIONS])
   {
        int rnum=0;
        rateconst[rnum]=0.0; 
   }
    
    AMREX_GPU_DEVICE 
    amrex::Real mobility(int specid,
                         amrex::Real Te, 
                         amrex::Real efield,
                         Real Tg, Real Pg)

    {
        return(0.0);
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
           dcoeff = 0.1;
        }
        else
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
