#include<Chemistry.H>
#include<VarDefines.H>

namespace plasmachem
{
    amrex::Vector<std::string> specnames(NUM_SPECIES);

    void init()
    {
        specnames[I_ID]="Ion";
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
}