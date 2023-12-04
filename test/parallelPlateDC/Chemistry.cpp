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
        specnames[I_ID]="Ion";
        
        for(int i=0;i<NUM_REACTIONS;i++)
        {
            for(int j=0;j<NUM_ALL_SPECIES;j++)
            {
                rct_rxnarray[i][j]=0.0;
                pdt_rxnarray[i][j]=0.0;
            }
        }

        int rnum=0;
        
        //E + N --> I + 2E
        rct_rxnarray[rnum][EDN_ID] = 1.0;
        rct_rxnarray[rnum][N_ID]   = 1.0;

        pdt_rxnarray[rnum][I_ID]   = 1.0;
        pdt_rxnarray[rnum][EDN_ID] = 2.0;
        
        reaction_elecenergy[rnum]  = -24.0*ECHARGE;
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
            case N_ID:
                charge=0.0;
                break;               
            case I_ID:
                charge=1.0;
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
        switch(specid)
        {
            case N_ID:
                molwt=4.0*M_AMU;
                break;               
            case I_ID:
                molwt=4.0*M_AMU;
                break;               
            case EDN_ID:
                molwt=ME;
                break;               
            default:
                molwt=M_AMU;
        }
        return(molwt);
    }
    
    AMREX_GPU_HOST_DEVICE void get_reaction_rateconstants(amrex::Real Te,amrex::Real Tg, amrex::Real Pg,
                                                   amrex::Real efield, amrex::Real spec[NUM_ALL_SPECIES],
                                                   amrex::Real rateconst[NUM_REACTIONS])
   {
        int rnum=0;
        //E + N --> I + 2E
        rateconst[rnum] = 2.5e-12
            *std::exp(-278508.0/Te);
   }
}
