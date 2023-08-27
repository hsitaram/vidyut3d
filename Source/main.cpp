#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <Vidyut.H>
#include <Chemistry.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        // timer for profiling
        BL_PROFILE_VAR("main()", pmain);

        // wallclock time
        const Real strt_total = amrex::second();

        plasmachem::init();
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        Vidyut vidyut_obj;

        // initialize AMR data
        vidyut_obj.InitData();
        
        // wallclock time
        const Real strt_evolve = amrex::second();

        // advance solution to final time
        vidyut_obj.Evolve();
        
        // wallclock time
        Real end_evolve = amrex::second() - strt_evolve;

        // wallclock time
        Real end_total = amrex::second() - strt_total;

        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total, ParallelDescriptor::IOProcessorNumber());
        if (vidyut_obj.Verbose())
        {
            amrex::Print() << "\nEvolve_Time: " <<  ParallelDescriptor::NProcs()<<"\t"<<end_total << '\n';
            amrex::Print() << "\nTotal_Time: " <<  ParallelDescriptor::NProcs()<<"\t"<<end_total << '\n';
        }

        // destroy timer for profiling
        BL_PROFILE_VAR_STOP(pmain);

        plasmachem::close();
    }

    amrex::Finalize();
}
