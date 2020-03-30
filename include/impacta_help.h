/*
**********************************************************
Help for IMPACT - not complete

Version 1.20
1
AGRT

12/4/07
**********************************************************

*/
int IMPACT_Help(int argnum, char **argstr)
{
  int result =0;
  result +=IMPACT_Cmp_CmdLine(argnum, argstr, "-help");
  result +=IMPACT_Cmp_CmdLine(argnum, argstr, "-h");
  if (result>0)
    {
      std::cout<<
	"Options for IMPACT\n\nDiagnostics:\n-echo_input: Prints information obtained from input deck\n";
      std::cout<<"-echo_infuncs: Prints the terms of the functions in the input deck for diagnostics\n";
      std::cout<<"-show_matrix: Prints the full matrix to screen if on a single processor\n";
      std::cout<<"-view_matrix: Graphically shows the full matrix (needs X-windows)\n";
      std::cout<<"-show_vector: Prints the full RHS vector to screen if on a single processor\n";
      std::cout<<"-show_all_initial: Prints all the initial conditions to screen\n";
      std::cout<<"-probe_input: Prints the variable names extracted from the input deck sequentially\n";
      std::cout<<"\nCommand line control:\n";
      std::cout<<"-rd <Root Dir>: Set root directory for IMPACT\n";
      std::cout<<"-read <imstdin>: Use input deck filename\n";
      std::cout<<"-iterate_matrix <value>: Iterate matrix solve with coefficient <value>\n";
      std::cout<<"-init_DLM_order <order>: Initialize f0 with DLM with order <order>\n";
      std::cout<<"-Fix_dE/dt <value>: Fix dE/dt=0 term for matrix solve (0-1)\n";
      std::cout<<"-B_implicit_frac <value>: Balance between B and f1 in f1 equation to be implicit\n";
      std::cout<<"-B_explicit: Calculate B explicitly\n";
      std::cout<<"-time_centering: E,j and B taken at t=n+1/2\n";
      std::cout<<"-no_calc_ni: Don't calculate ni from Z and ne, get from input deck\n";
      std::cout<<"-max_picard_its <value>: Specify maximum Picard iterations, 0 is no maximum (default = "<<zerotolerance::max_picard_its<<")\n";
      std::cout<<"-static_ions: Ion motion switched off\n";
      std::cout<<"-LPC: Use linear solution as start solution for matrix\n";
      std::cout<<"-EPC: Use explicit solution as start solution for matrix\n";
      std::cout<<"-KSP_lower_initial_rtol: On first lagged loop, lower rtol to "<<zerotolerance::low_rtol*100<<" %\n";
      std::cout<<"-no_dump: Do not write data to file\n";
      std::cout<<"-using_xnbody: Use Xnbody to view data\n";
      std::cout<<"-no_div_check: Don't check for divergence in Picard iteration\n";
      std::cout<<"-renorm_temp <value>: Renormalize temperature to value. Smaller than one should help converge the matrix\n";
      // Note the following are now in the unput deck
      /* std::cout<<"-fixed: all boundaries fixed\n";
      std::cout<<"-fixed_i: x boundaries fixed\n";
      std::cout<<"-fixed_j: y boundaries fixed\n";
      std::cout<<"-fixed_iu: x upper boundary fixed\n";
      std::cout<<"-fixed_il: x lower boundary fixed\n";
     std::cout<<"-fixed_ju: y upper boundary fixed\n";
     std::cout<<"-fixed_jl: y lower boundary fixed\n";*/
      exit(0);
    }
  return result;
}
