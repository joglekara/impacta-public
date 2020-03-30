/*

Main program tracing a ray through a density map then mapping the intensity onto a grid.
Bilinear Interpolation used for interpolating density matrix as well as mapping intensity.
Timestep is declared in main()
ReadinGrids() retrieves x, y and n. 
Traceonestep traces one step through x,y,n,b

Archis Joglekar
02/01/2013

03/07 - Incoroporated into IMPACTA through impacta_filehandling.h and impacta_headers.h

08/2013 - Added X wave and calculation of dielectric constant

*/
void intensity_diffusion(IMPACT_Config *conf)
	//int Nx, int Ny, double sx, double sy)
{	
	int iminus,iplus,jminus,jplus;
	double sx, sy;
	IMPACT_Matrix tempint;
	tempint.ChangeSize(conf->Nx(),conf->Ny());
	tempint.setall(0.0);

	for (int t=1;t<IMPACT_Heating::diffsteps;++t)
	{	
		for (int i=1;i<=conf->Nx();++i)
		{
			for (int j=1;j<=conf->Ny();++j)
			{	
				sx=1e-1;sy=1e-1;
				if (i==1)
				{
					iminus = 1;
					iplus = 2;
				}
				else if (i==conf->Nx())
				{
					iplus = conf->Nx();
					iminus = conf->Nx()-1;
				}
				else
				{
					iminus = i-1;
					iplus = i+1;
				}
				if (j==1)
				{
					jminus = 1;
					jplus = 2;
				}
				else if (j==conf->Ny())
				{
					jplus = conf->Ny();
					jminus = conf->Ny()-1;
				}
				else
				{
					jminus = j-1;
					jplus = j+1;
				}

				// sy *= (conf->dy(&jminus)*conf->dy(&jminus))/(conf->dx(&iminus)*conf->dx(&iminus)+conf->dy(&jminus)*conf->dy(&jminus));
				// sx *= (conf->dx(&iminus)*conf->dx(&iminus))/(conf->dx(&iminus)*conf->dx(&iminus)+conf->dy(&jminus)*conf->dy(&jminus));

				//sx *= 1.0/(1.0+pow(conf->dy(&jminus)/conf->dx(&iminus),4));
				//sy *= 1.0/(1.0+pow(conf->dx(&iminus)/conf->dy(&jminus),4));

				tempint.set(i,j,IMPACT_Heating::i_mat.get(i,j)
					+sx*(IMPACT_Heating::i_mat.get(iplus,j)-2*IMPACT_Heating::i_mat.get(i,j)+IMPACT_Heating::i_mat.get(iminus,j))
					+sy*(IMPACT_Heating::i_mat.get(i,jplus)-2*IMPACT_Heating::i_mat.get(i,j)+IMPACT_Heating::i_mat.get(i,jminus)));
			}
		}
		for (int i=1;i<=conf->Nx();++i)
		{
			for (int j=1;j<=conf->Ny();++j)
			{
				IMPACT_Heating::i_mat.set(i,j,tempint.get(i,j));
			}
		}
	}
}
void updatelocalintensity(IMPACT_Matrix &tempint,IMPACT_Matrix &localint)
{
	int Nx, Ny;
	localint.size(&Nx,&Ny);

	for (int i=1;i<=Nx;++i)
	{
		for (int j=1;j<=Ny;++j)
		{
			localint.set(i,j,localint.get(i,j)+tempint.get(i,j));
		}
	}
}
void intensitycalc(IMPACT_Config *conf, IMPACT_Matrix &n_mat, IMPACT_Vector &ray_vec, 
					IMPACT_Matrix &localint, double mult)
{
	
	int tempi;
	double tempv;
	bool check;

	// int lowx=0, uppx=conf->Nx()+1, lowy=0, uppy=conf->Ny()+1;
	int lowx=-1, uppx=conf->Nx()+1, lowy=-1, uppy=conf->Ny()+1;
	double xpos = ray_vec.Get(1);
	double ypos = ray_vec.Get(2);


	// if (std::isnan(tempvG))
	// {
	// 	tempvG=1e12;
	// }

	// Find indices for x
	
	check = 0;
	
	// tempi = (uppx-lowx)/2;
	// while (!check)
	// {
	// 	tempv = conf->xpos(tempi);
	// 	if (xpos < tempv)
	// 	{
	// 		uppx=tempi;
	// 	}
	// 	else if (xpos > tempv)
	// 	{
	// 		lowx=tempi;
	// 	}
	// 	else if (xpos == tempv)
	// 	{
	// 		lowx=tempi;
	// 		uppx=tempi+1;
	// 	}
	// 	check = ((uppx-lowx)==1);
	// 	tempi = lowx+(uppx-lowx)/2;	
	// }

	// tempxuppV = fabs(xpos-conf->xpos(lowx))/(conf->xpos(uppx)-conf->xpos(lowx));
	// tempxlowV = 1-tempxuppV;

	// // Find indices for y
	// check = 0;
	
	// tempi = (uppy-lowy)/2;
	// while (!check)
	// {
	// 	tempv = conf->ypos(tempi);
	// 	if (ypos < tempv)
	// 	{
	// 		uppy=tempi;
	// 	}
	// 	else if (ypos > tempv)
	// 	{
	// 		lowy=tempi;
	// 	}
	// 	else if (ypos == tempv)
	// 	{
	// 		lowy=tempi;
	// 		uppy=tempi+1;
	// 	}
	// 	check = ((uppy-lowy)==1);
	// 	tempi = lowy+(uppy-lowy)/2;	
	// }

	// tempyuppV = fabs(ypos-conf->ypos(lowy))/(conf->ypos(uppy)-conf->ypos(lowy));
	// tempylowV = 1-tempyuppV;

	// mult /= ((conf->ypos(uppy)-conf->ypos(lowy))*(conf->xpos(uppx)-conf->xpos(lowx)));


	// if (uppx == conf->Nx()+1) uppx=0;
	// if (uppy == conf->Ny()+1) uppy=0;

	// if ((lowx*lowy*uppx*uppy) == 0)
	// {
	// 	if (lowx==0 && (lowy*uppy) !=0 )		// Left
	// 	{
	// 		localint.set(uppx,lowy,localint.get(uppx,lowy)+mult*(tempxuppV*tempylowV));
	// 		localint.set(uppx,uppy,localint.get(uppx,uppy)+mult*(tempxuppV*tempyuppV));
	// 	}
	// 	else if (lowy == 0 && (lowx * uppx)!=0 ) // Bottom
	// 	{
	// 		localint.set(lowx,uppy,localint.get(lowx,uppy)+mult*(tempxlowV*tempyuppV));
	// 		localint.set(uppx,uppy,localint.get(uppx,uppy)+mult*(tempxuppV*tempyuppV));		
	// 	}
	// 	else if (uppx == 0 && (lowy * uppy)!=0) // Right
	// 	{
	// 		localint.set(lowx,lowy,localint.get(lowx,lowy)+mult*(tempxlowV*tempylowV));
	// 		localint.set(lowx,uppy,localint.get(lowx,uppy)+mult*(tempxlowV*tempyuppV));
	// 	}
	// 	else if (uppy == 0 && (lowx * uppx)!=0)	// Top
	// 	{
	// 		localint.set(lowx,lowy,localint.get(lowx,lowy)+mult*(tempxlowV*tempylowV));
	// 		localint.set(uppx,lowy,localint.get(uppx,lowy)+mult*(tempxuppV*tempylowV));
	// 	}
	// 	else if (lowx == 0 && lowy == 0) // Lower Left
	// 	{
	// 		localint.set(uppx,uppy,localint.get(uppx,uppy)+mult*(tempxuppV*tempyuppV));
	// 	}
	// 	else if (lowx == 0 && uppy == 0) // Upper Left
	// 	{
	// 		localint.set(uppx,lowy,localint.get(uppx,lowy)+mult*(tempxuppV*tempylowV));
	// 	}
	// 	else if (uppx == 0 && lowy == 0) // Lower Right
	// 	{
	// 		localint.set(lowx,uppy,localint.get(lowx,uppy)+mult*(tempxlowV*tempyuppV));
	// 	}
	// 	else if (uppx == 0 && uppy == 0) // Upper Right
	// 	{
	// 		localint.set(lowx,lowy,localint.get(lowx,lowy)+mult*(tempxlowV*tempylowV));
	// 	}
	// }
	// else	// Rest
	// {
	// 	localint.set(lowx,lowy,localint.get(lowx,lowy)+mult*(tempxlowV*tempylowV));
	// 	localint.set(lowx,uppy,localint.get(lowx,uppy)+mult*(tempxlowV*tempyuppV));
	// 	localint.set(uppx,lowy,localint.get(uppx,lowy)+mult*(tempxuppV*tempylowV));
	// 	localint.set(uppx,uppy,localint.get(uppx,uppy)+mult*(tempxuppV*tempyuppV));
	// }

	tempi = (uppx-lowx)/2;
	while (!check)
	{
		tempv = conf->xb(tempi);
		if (xpos < tempv)
		{
			uppx=tempi;
		}
		else if (xpos > tempv)
		{
			lowx=tempi;
		}
		else if (xpos == tempv)
		{
			lowx=tempi;
			uppx=tempi+1;
		}
		check = ((uppx-lowx)==1);
		tempi = lowx+(uppx-lowx)/2;	
	}



	// Find indices for y
	check = 0;
	
	tempi = (uppy-lowy)/2;
	while (!check)
	{
		tempv = conf->yb(tempi);
		if (ypos < tempv)
		{
			uppy=tempi;
		}
		else if (ypos > tempv)
		{
			lowy=tempi;
		}
		else if (ypos == tempv)
		{
			lowy=tempi;
			uppy=tempi+1;
		}
		check = ((uppy-lowy)==1);
		tempi = lowy+(uppy-lowy)/2;	
	}
	// 		std::cout << "xpos=" << xpos << "ypos=" << ypos <<std::endl;
	// std::cout << "lowx=" << lowx << ", uppx=" << uppx << std::endl;
	if (xpos == conf->xb(lowx) && ypos != conf->yb(lowy))
	{
		if (lowx==0)
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		else if (lowx==conf->Nx())
		{
			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
		}
		else
		{
			// std::cout << "hello2=" << std::endl;
			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
	}
	else if (xpos != conf->xb(lowx) && ypos == conf->yb(lowy))
	{
		if (lowy==0)
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		else if (lowy==conf->Ny())
		{
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.5/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else
		{
			// std::cout << "hello3=" << std::endl;
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.5/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.5/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
	}
	else if (xpos == conf->xb(lowx) && ypos == conf->yb(lowy))
	{
		if (lowx==0 && lowy==0)
		{				
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else if (lowx==conf->Nx() && lowy==conf->Ny())
		{
			localint.set(lowx,lowy,localint.get(lowx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
		}
		else if (lowx==0 && lowy==conf->Ny())
		{
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else if (lowx==conf->Nx() && lowy==0)
		{
			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
		}
		else if (lowx==0)
		{
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else if (lowy==0)
		{
			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else if (lowx==conf->Nx())
		{
			localint.set(lowx,lowy,localint.get(lowx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
		}
		else if (lowy==conf->Ny())
		{
			localint.set(lowx,lowy,localint.get(lowx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}
		else
		{
			// std::cout << "hello4=" << std::endl;
			localint.set(lowx,lowy,localint.get(lowx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(uppx,lowy,localint.get(uppx,lowy)+0.25/((conf->yb(lowy)-conf->yb(lowy-1))*(conf->xb(uppx)-conf->xb(lowx)))*mult);

			localint.set(lowx,uppy,localint.get(lowx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(lowx)-conf->xb(lowx-1)))*mult);
			localint.set(uppx,uppy,localint.get(uppx,uppy)+0.25/((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)))*mult);
		}	
	}
	else
	{	
		mult /= ((conf->yb(uppy)-conf->yb(lowy))*(conf->xb(uppx)-conf->xb(lowx)));
		// mult /= ((conf->dy(&uppy))*(conf->dx(&uppx)));
		localint.set(uppx,uppy,localint.get(uppx,uppy)+mult);
	}

	// }

}

void traceonestep(IMPACT_Config *conf,IMPACT_Vector &ray_vec,
		IMPACT_Matrix &eps_mat,IMPACT_Matrix &gradepsx, IMPACT_Matrix &gradepsy,
	double eps, double ds, bool &flag, bool &firststep)
{
	IMPACT_Vector r_h(3), r0(3), r1(3);
	IMPACT_Vector r_temp(3);
	IMPACT_Vector v0(3), v1(3);
	IMPACT_Vector gradn_h(3),gradb_h(3),gradI_h(3), omega_0(3), omega_h(3), grad_eps0(3);

	IMPACT_Vector temp(3), temp2(3), tempVstep(3);
	
	double n_h, tempdouble,tempdouble2,tempindex;

	
	// double xmx = x_vec.Get(x_vec.length())+dx/2;
	// double ymx = y_vec.Get(y_vec.length())+dy/2;
	// double xmn = x_vec.Get(1)-dx/2;
	// double ymn = y_vec.Get(1)-dy/2;

	double vals [3];

	r0.Set(1,ray_vec.Get(1));
	r0.Set(2,ray_vec.Get(2));
	r0.Set(3,0.0);

	v0.Set(1,ray_vec.Get(3));
	v0.Set(2,ray_vec.Get(4));
	v0.Set(3,0.0);

	// First Half Step  
	// STEP 1 
	temp.Duplicate(&v0);  	 // DO I WANT GETVEC OR DUPLICATE?
	temp.MultiplyAll(ds/2.0);  // temp is now v0*ds/2
	addVecs(r0,temp,r_h);

	//std::cout << "\n\n";gradepsy.Print();std::cout << "\n\n";

	interp23(conf,gradepsx,gradepsy,eps_mat,r_h.Get(1),r_h.Get(2),vals);
	gradn_h.Set(1,vals[0]); //INTERPOLATIONS!
	gradn_h.Set(2,vals[1]);
	gradn_h.Set(3,0.0);

	n_h = vals[2]-eps;  
	gradn_h.MultiplyAll(0.5/n_h);  // Is now grad index / index
	//	std::cout << "gradnh is: "; gradn_h.PrintArray(); std::cout<< "and n_h is: " << n_h << "\n\n";


	// STEP 2
	crossVecs(gradn_h,temp,omega_0); // With ds
	//	std::cout << "Omega_0 is: "; omega_0.PrintArray();

	// Step 3
	temp.Duplicate(&omega_0);
	// temp.MultiplyAll(ds/2.0);

	crossVecs(v0,omega_0,temp);

	addVecs(v0,temp,temp);
	temp.MultiplyAll(ds/2.0);
	crossVecs(gradn_h,temp,omega_h);

	//	std::cout << "Omega_(1/2) is: "; omega_h.PrintArray();

	// Step 4
	// double tempVstep;
	// tempVstep = omega_h.VecSum();
	// tempVstep = 2.0/(1.0+tempVstep*ds*ds/4.0);

	specialVstep(omega_h,ds,tempVstep);

	crossVecs(v0,omega_h,temp);
	// temp.MultiplyAll(ds/2.0);
	addVecs(v0,temp,temp2);
	crossVecs(temp2,omega_h,temp2);
	// temp2.MultiplyAll(ds/2.0);

	// temp2.MultiplyAll(tempVstep);
	elementMultVecs(temp2,tempVstep,temp2);
	addVecs(v0,temp2,v1);

	//	std::cout << "v_1 is: "; v1.PrintArray();	

	v1.MultiplyAll(ds/2.0);

	// Check for correctness 
	addVecs(r_h,v1,r_temp);
	//	std::cout << "r_temp is: "; r_temp.PrintArray(); std::cout<< "\n";

	// INBOUNDS
	// if (r_temp.Get(1) <= conf->xpos(conf->Nx()+1) && r_temp.Get(2) <= conf->ypos(conf->Ny()+1) 
	// 	&& r_temp.Get(1) >= conf->xpos(0) && r_temp.Get(2) >= conf->ypos(0))
	if (r_temp.Get(1) < conf->xb(conf->Nx()) && r_temp.Get(2) < conf->yb(conf->Ny()) 
		&& r_temp.Get(1) > conf->xb(0) && r_temp.Get(2) > conf->yb(0))
	{
		tempindex = (interp2(conf,eps_mat,r_temp.Get(1),r_temp.Get(2))-eps);

		//std::cout << "tempindex = " << tempindex << "\n";

		if (tempindex > 0) // Interpolate
		{   // Step 5 If allowed to enter 
			ray_vec.Set(1,r_temp.Get(1));
			ray_vec.Set(2,r_temp.Get(2));
			ray_vec.Set(3,2.0/ds*v1.Get(1));
			ray_vec.Set(4,2.0/ds*v1.Get(2));       	  
		}
		else if (firststep)
		{
			flag = 0;
		}
		else
		{

			// Reflection if not allowed
			grad_eps0.reset();
			
			// Gradient calculation @ r0
			grad_eps0.Set(1,interp2(conf,gradepsx,r0.Get(1),r0.Get(2))); // Interpolate
			grad_eps0.Set(2,interp2(conf,gradepsy,r0.Get(1),r0.Get(2))); // Interpolate
			grad_eps0.Set(3,0.0);

			tempdouble=dotVecs(v0,grad_eps0);
			tempdouble2=grad_eps0.VecSum();
			
			double LcosA= -1.0*tempdouble/tempdouble2;
			
			
			grad_eps0.MultiplyAll(LcosA);
			// std::cout << "\n\n";grad_eps0.PrintArray();std::cout << "\n\n";
			//std::cout << "\n\n"<< LcosA << "\n\n";
			
			IMPACT_Vector v_par(3);
			
			v_par.Duplicate(&grad_eps0);	        
			v_par.MultiplyAll(2.0);

			addVecs(v0,v_par,v1);
			ray_vec.Set(3,v1.Get(1));
			ray_vec.Set(4,v1.Get(2));

			// std::cout << "vx = " << v1.Get(1) << "\n";
			// std::cout << "vy = " << v1.Get(2) << "\n";
			
			// v1.MultiplyAll(ds/2);
			tempdouble = interp2(conf,eps_mat,r0.Get(1),r0.Get(2));

			// std::cout << "eps = " << tempdouble << "\n";

			// std::cout << "\n\n";v0.PrintArray();std::cout << "\n\n";
			// std::cout << "\n\n";grad_eps0.PrintArray();std::cout << "\n\n";

			addVecs(v0,grad_eps0,v1);


			v1.MultiplyAll(4*LcosA*tempdouble);
			//std::cout << "\n\n";v1.PrintArray();std::cout << "\n\n";


			ray_vec.Set(1,r0.Get(1)+v1.Get(1));
			ray_vec.Set(2,r0.Get(2)+v1.Get(2));


			 // std::cout << "dx = " << v1.Get(1) << "\n";
	   //       std::cout << "dy = " << v1.Get(2) << "\n";

			// addVecs(r_h,v1,r1);
			// ray_vec.Set(1,r1.Get(1));
			// ray_vec.Set(2,r1.Get(2));
		} 
	}
	else // OUT OF BOUNDS
	{

		flag=0;
	} 

}
void tracerunner(IMPACT_Matrix &n_mat, IMPACT_Matrix &Bz_mat,  IMPACT_StenOps *sten, IMPACT_Config *conf, IMPACT_MPI_Config *M)
{
	// For the Main Loop
	bool flag=1, firststep=1, discard = 0;
	IMPACT_Matrix tempint, localint;

	localint.ChangeSize(conf->Nx(),conf->Ny());
	localint.setall(0.0);

	tempint.ChangeSize(conf->Nx(),conf->Ny());
	tempint.setall(0.0);
	
	// for (int zs=0;zs<conf->Nx()+1;++zs)
	// std::cout<<"xb[0]: " << conf->xb(zs)<< "\n";

	
///////////////------------------------------------------------------------///////////////	
	// Initialization of plasma properties

	// IMPACT_Matrix gradnx, gradny;
	// IMPACT_Matrix gradxBz, gradyBz;
	IMPACT_Matrix eps_mat,gradepsx,gradepsy;

	IMPACT_Heating::i_mat.ChangeSize(conf->Nx(),conf->Ny());
	IMPACT_Heating::i_mat.setall(0.0);

	// double xmin=x_vec.Get(1)-dx_vec.Get(1);
	// double ymin=y_vec.Get(1)-dy_vec.Get(1);
	// double xmax=x_vec.Get(conf->Nx)+dx_vec.Get(conf->Nx);
	// double ymax=y_vec.Get(conf->Ny)+dy_vec.Get(conf->Ny);
			// double xmin=x_vec.Get(1)-dx/2;
	// double ymin=y_vec.Get(1)-dy/2;
	// double xmax=x_vec.Get(Nx)+dx/2;
	// double ymax=y_vec.Get(conf->Ny)+dy/2;


	/////////////////////////////////////////////////////////////
	// Determine STEP SIZE //////////////////////////////////////
	int n_max;
	
	int i =1;
	double dxmin = conf->dx(&i);
	double dymin = conf->dy(&i);
	double ds = IMPACT_Heating::ds;
	
	for (i=2;i<=conf->Nx();++i)
	{
		if (dxmin > conf->dx(&i))
		{
			dxmin = conf->dx(&i);
		}
	}

	for (i=2;i<=conf->Ny();++i)
	{
		if (dymin > conf->dy(&i))
		{
			dymin = conf->dy(&i);
		}
	}


	if (dxmin > dymin)
	{
		ds *= dymin;
		n_max = (conf->get_xmax()-conf->get_xmin())/ds*4;
	}
	else 
	{
		//globalconsts::pi*
		ds *= dxmin; 
		n_max = (conf->get_ymax()-conf->get_ymin())/ds*4;
	}
	// std::cout<< "\n\n ds: " << ds << "\n\n";
	/////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////
	
	epscalc(eps_mat,n_mat,Bz_mat);
	gradcalc(sten,conf,eps_mat,gradepsx,gradepsy);

///////////////------------------------------------------------------------///////////////	
	
	// Ray Trace Quantities
	IMPACT_Vector rayvec(4);double xc;
	//IMPACT_Vector rayfull(4);
	int r; // Position
	int n=1; // Time Step
	
	double theta_live;
	double mult;
	
///////////////------------------------------------------------------------///////////////	
	// VARIATION IN CRITICAL DENSITY ACCORDING TO AIRY FUNCTION
	double tempeps; // For Airy noise
	double deltac=globalconsts::c_L/globalconsts::omega_p_nuei_n; // Airy

	double max_loc=findindex(conf,eps_mat,0,IMPACT_Heating::direction);

	// bool airyflag = max_loc>1;
	bool airyflag=0;

	double scln;
	// = IMPACT_Heating::direction*fabs(1.0/interp2(x_vec,y_vec,gradepsy,floor(Nx/2.0),max_loc))+
	// (-1.0*IMPACT_Heating::direction+1)*
	// 	fabs(1.0/interp2(x_vec,y_vec,gradepsx,max_loc,floor(Ny/2.0)));
	if (airyflag)
	{
		if (IMPACT_Heating::direction)
		{
			scln = fabs(1.0/interp2(conf,gradepsy,conf->xpos(conf->Nx()/2),max_loc));
	
		}
		else
		{
			scln = fabs(1.0/interp2(conf,gradepsx,max_loc,conf->ypos(conf->Ny()/2)));
	
		}	
	}
	// // FIND ETA = 1,  LOCATION OF MAX INTENSITY
	max_loc-=pow(deltac*deltac*scln/IMPACT_Heating::n_c,1.0/3.0);

	// For Airy Loop	
	bool ncflg;
	// bool airyflag = 1;
	
///////////////-------------------------------------------------------------///////////////

	boost::posix_time::ptime time_t_epoch(boost::gregorian::date(2014,1,1));
	boost::posix_time::ptime tick = boost::posix_time::microsec_clock::local_time();
	
	boost::posix_time::time_duration diff = tick - time_t_epoch;
	
	long x = diff.total_microseconds();

	// Random Number Generation
	typedef boost::mt19937 RNGType;
	long seed = 1e4*(M->rank())+x;

//	printf("Seed: %.17g\n", seed);

	 // std::cout<< "\n\n microseconds: " << x << "\n\n";
	// std::cout<< "\n\n seed: " << seed << "\n\n";


	RNGType rng(seed);
	//boost::uniform_real<> degreesample(0,1.0);
	boost::uniform_real<> Xsample(-1.0,2.0);
	boost::uniform_real<> Ysample(0.0,1.0);
	boost::uniform_real<> degsample(-1.0,1.0);

	//boost::variate_generator< RNGType, boost::uniform_real<> > deg_gen(rng, degreesample);
	boost::variate_generator< RNGType, boost::uniform_real<> > x_gen(rng, Xsample);
	boost::variate_generator< RNGType, boost::uniform_real<> > y_gen(rng, Ysample);	
	boost::variate_generator< RNGType, boost::uniform_real<> > deg_gen(rng, degsample);	

///////////////------------------------------------------------------------///////////////	
	// Live Loop Quantity containers
	double livemult;
	double pertcoeff;
	double temploc;
	double dx_orig = conf->xb(1)-conf->xb(0); 
	double dy_orig = conf->yb(1)-conf->yb(0);  

	// Group Velocity
	double vG0;
	double ntemp;
	double btemp;
	// Airy Loop
	double Xtemp,Ytemp;
	
	// each step
	double tempvG;
	bool check;
	int tempi;
	double tempv;
	int upp,low;

	// Ray Diagnostics!!!!!!!!!
	// 
	std::string dir = "MNT/Te/";
	dir = IMPACT_Messages::Data_Directory + dir;
	std::string rayfile = "ray_";
	std::string raydmpname = IMPACT_Out(dir,rayfile,0);
	
	const int ifnotdump = 0; //((*tstep+1) % ndump)+1;

	//std::cout<< "\n ifnotdump : " << ifnotdump << "\n";
	int rayndump=1;//IMPACT_Heating::beam_res/10;
	int ifraydump;

	std::ofstream raystream;

	if (!ifnotdump && !(M->rank())) raystream.open(raydmpname.c_str());

	bool bool_u = 0,bool_s = 0,bool_g = 0; 

	
	if (IMPACT_Heating::shape==0) 	{bool_u = 1; bool_s = 0; bool_g = 0;}
	if (IMPACT_Heating::shape==1) 	{bool_u = 0; bool_s = 1; bool_g = 0;}
	if (IMPACT_Heating::shape==2) 	{bool_u = 0; bool_s = 0; bool_g = 1;}

	int local_ray_ind=1+M->rank()*IMPACT_Heating::beam_res;
	int local_ray_total=(M->rank()+1)*IMPACT_Heating::beam_res;
	int totalbeamres = M->size()*IMPACT_Heating::beam_res;
	double sigma = IMPACT_Heating::beam_width/sqrt(8*log(2));
	double totalGaussianbeamwidth = 2.0*sigma;
	bool currentraydir;

	double x0,y0;
	double turnx,turny;
	bool turn;

	int beam;
	int totalbeams = 4;

	mult = (bool_g*2.0*sigma+(!bool_g)*IMPACT_Heating::beam_width)/IMPACT_Heating::beam_res*ds/dx_orig/dy_orig;
	
	// mult = IMPACT_Heating::beam_width/totalbeamres*ds/dx/dy;

	// double sx = 0.01*(dymin*dymin)/(dxmin*dxmin+dymin*dymin);
	// double sy = 0.01*(dxmin*dxmin)/(dxmin*dxmin+dymin*dymin);

	// for (r=1;r<=IMPACT_Heating::beam_res;++r)
	for (r=local_ray_ind;r<=local_ray_total;++r)
	{	
		for (beam=1;beam<=totalbeams;++beam)
		{
			flag = 1;
			turn = 0;
			firststep=1;
			n = 1;
			ncflg=1;


			if (beam==1)
			{
				// INNER BEAM 1 (23.5 degrees)
				currentraydir = IMPACT_Heating::direction;
				theta_live=IMPACT_Heating::theta;
				
				sigma = IMPACT_Heating::beam_width/sqrt(8*log(2))/sin(theta_live);
				// y0 = 6.5*IMPACT_Heating::ray_y0;
				x0 = IMPACT_Heating::ray_x0;
				y0 = IMPACT_Heating::ray_y0;
				
				livemult = mult;//livemult = mult;
			}
			else if (beam==2)
			{
				// INNER BEAM 2 (30 degrees)
				currentraydir = !IMPACT_Heating::direction;
				theta_live=30.0*globalconsts::pi/180.0;
				
				sigma = IMPACT_Heating::beam_width/sqrt(8*log(2))/sin(theta_live);
				// y0 = 12.0*IMPACT_Heating::ray_y0;
				x0 = IMPACT_Heating::ray_x0;
				y0 = 1300+IMPACT_Heating::ray_y0;

				livemult = mult;//livemult = mult;
			}
			else if (beam==3)
			{
				// OUTER BEAM 1 (44.5 degrees)
				currentraydir = IMPACT_Heating::direction;
				theta_live=45.5*globalconsts::pi/180.0;
				
				sigma = IMPACT_Heating::beam_width/sqrt(8*log(2))/sin(theta_live);
				sigma = 1.5*sigma;
				// totalGaussianbeamwidth = 1.5*sigma;
				
				x0 = 100+IMPACT_Heating::ray_x0;
				y0 = IMPACT_Heating::ray_y0;
				livemult = mult;
			}
			else if (beam==4)
			{
				// OUTER BEAM 1 (50 degrees)
				currentraydir = !IMPACT_Heating::direction;
				theta_live=50.0*globalconsts::pi/180.0;

				sigma = IMPACT_Heating::beam_width/sqrt(8*log(2))/sin(theta_live);
				sigma = 1.5*sigma;
				// totalGaussianbeamwidth = 1.5*sigma;

				// y0 = 4.0*IMPACT_Heating::ray_y0;
				x0 = IMPACT_Heating::ray_x0;
				y0 = 150+IMPACT_Heating::ray_y0;

				livemult = mult;

			}
			
			// dTheta goes here
			// theta_live=IMPACT_Heating::theta+IMPACT_Heating::dtheta*deg_gen();
			theta_live=theta_live+IMPACT_Heating::dtheta*(1.0-2.0*(r-1)/totalbeamres);
	
			// Perturbation to intensity profile goes here
			// pertcoeff = 1+0.0*cos(4*2.0*globalconsts::pi*r/IMPACT_Heating::beam_res);
	
			// Initialize Ray positions and directions. 
			// Also initialize intensity
			// if (!IMPACT_Heating::direction)
			if (!currentraydir)
			{	
				xc = x0+(2.0*sigma)/2.0;
				low = 0;upp=conf->Nx()+1;
				// rayvec.Set(1,IMPACT_Heating::ray_x0+(r-0.5)/IMPACT_Heating::beam_res*IMPACT_Heating::beam_width);
				rayvec.Set(1,x0+(r-0.5)/totalbeamres*
					(bool_g*(2.0*sigma)+(!bool_g)*IMPACT_Heating::beam_width));
				rayvec.Set(2,conf->yb(0)+0.1*ds);
				rayvec.Set(3,cos(theta_live));
				rayvec.Set(4,sin(theta_live));
	
	
				check = 0;
		
				tempi = (upp-low)/2;
				while (!check)
				{
					tempv = conf->xpos(tempi);
					if (rayvec.Get(1) < tempv)
					{
						upp=tempi;
					}
					else if (rayvec.Get(1) > tempv)
					{
						low=tempi;
					}
					else if (rayvec.Get(1) == tempv)
					{
						low=tempi;
						upp=tempi+1;
					}
					check = ((upp-low)==1);
					tempi = low+(upp-low)/2;	
				}
	
	
				livemult = livemult*(conf->xpos(upp)-conf->xpos(low))*dy_orig*
					fabs(rayvec.Get(4))*((bool_u)+
					(bool_s*(pow(sin((rayvec.Get(1)-x0)/IMPACT_Heating::beam_width*globalconsts::pi),2)))+
					(bool_g*exp(-1.0/2.0*pow((rayvec.Get(1)-xc)/sigma,4.0))));
				
	
				// (bool_g*(exp(-1.0*pow((rayvec.Get(1)-xc)/(IMPACT_Heating::beam_width/6.0),2.0)/2))));
				// livemult*= pertcoeff;
	
				// Ray Diag
				//rayndump = IMPACT_Heating::beam_res/10;//(IMPACT_Heating::beam_width/dx);
	
			}
			else
			{
				xc = y0+(2.0*sigma)/2.0;
				low = 0;upp=conf->Ny()+1;
	
				rayvec.Set(1,conf->xb(0)+0.1*ds);		
				// rayvec.Set(2,IMPACT_Heating::ray_y0+(r-0.5)/IMPACT_Heating::beam_res*IMPACT_Heating::beam_width);
				// rayvec.Set(2,y0+(r-0.5)/totalbeamres*
				rayvec.Set(2,y0+(r-0.5)/totalbeamres*
					(bool_g*(2.0*sigma)+(!bool_g)*IMPACT_Heating::beam_width));
				rayvec.Set(3,sin(theta_live));
				rayvec.Set(4,cos(theta_live));
	
	
				check = 0;
		
				tempi = (upp-low)/2;
				while (!check)
				{
					tempv = conf->ypos(tempi);
					if (rayvec.Get(2) < tempv)
					{
						upp=tempi;
					}
					else if (rayvec.Get(2) > tempv)
					{
						low=tempi;
					}
					else if (rayvec.Get(2) == tempv)
					{
						low=tempi;
						upp=tempi+1;
					}
					check = ((upp-low)==1);
					tempi = low+(upp-low)/2;	
				}
	
	
	
				livemult = livemult*(conf->ypos(upp)-conf->ypos(low))*dx_orig*fabs(rayvec.Get(3))*((bool_u)+
					(bool_s*(pow(sin((rayvec.Get(2)-y0)/IMPACT_Heating::beam_width*globalconsts::pi),2)))+
					// (bool_g*(exp(-1.0*pow((rayvec.Get(2)-xc)/(IMPACT_Heating::beam_width/6.0),2.0)/2))));
					(bool_g*exp(-1.0/2.0*pow((rayvec.Get(2)-xc)/sigma,4.0))));

				// livemult*= pertcoeff;
	
				// Ray Diag
				//rayndump = 4*IMPACT_Heating::beam_res/(IMPACT_Heating::beam_width/dy);
	
			}	
	
			// OBTAIN CRITICAL DENSITY VALUES ACCORDING TO AIRY FUNCTION ASYMPTOTIC EXPANSION USING RUDIMENTARY
			// MONTE CARLO
	
			if (airyflag)
			{		
				while (ncflg)
				{
					Xtemp = x_gen();
					Ytemp = y_gen();
	
					if (3.48979847492*pow(boost::math::airy_ai(Xtemp),2.0)>=Ytemp)
					{
						Xtemp+=1.0;
						Xtemp*=pow(deltac*deltac*scln/IMPACT_Heating::n_c,1.0/3.0);
	
						// std::cout<< "\n\n scln is " << Xtemp << "\n\n";
	
						// Xtemp*=pow(1.0*1.0*scln/IMPACT_Heating::n_c,1.0/3.0);
	
						temploc=max_loc+Xtemp;
	
						ncflg=0;
	
						// std::cout<< "\n\n temploc is " << temploc << "\n\n";
					}
				}
				if (IMPACT_Heating::direction)
				{			
					if (temploc > conf->ypos(conf->Ny())) {temploc = conf->ypos(conf->Ny());}
					tempeps = interp2(conf,eps_mat,conf->xpos(floor(conf->Nx()/2.0)),temploc);
				}
				else
				{	
					if (temploc > conf->xpos(conf->Nx())) {temploc = conf->xpos(conf->Nx());}	
					tempeps = interp2(conf,eps_mat,temploc,conf->xpos(floor(conf->Ny()/2.0)));
	
				}
			}
			else
			{
				tempeps = 0.0; //0.02*deg_gen();
			}
			
			//tempeps = 0;
			ntemp = interp2(conf,n_mat,rayvec.Get(1),rayvec.Get(2))/IMPACT_Heating::n_c;
			btemp = pow(interp2(conf,Bz_mat,rayvec.Get(1),rayvec.Get(2))/globalconsts::omega_p_nuei_n,2.0)/IMPACT_Heating::n_c;
	
			vG0 = sqrt(interp2(conf,eps_mat,rayvec.Get(1),rayvec.Get(2))-tempeps); // sqrt(epsilon)
	
			vG0 = (vG0)/(vG0*vG0+
				ntemp*((1-ntemp)*(1-ntemp)+ntemp*btemp)/(1-ntemp-btemp)/(1-ntemp-btemp));
			
			discard=0;
			livemult=vG0*livemult;
			// Ray Diagnostics
			// Ray Diag
			
			ifraydump = 0; //(r % rayndump);
	
			// std::cout << "xpos=" << rayvec.Get(1) << "ypos=" << rayvec.Get(2) <<std::endl;
	
			// MAIN LOOP
			while (flag)
			{	
				if ((!ifraydump) && !(M->rank()) && !ifnotdump)
				{ 
					rayvec.PrintArraytoFile(raystream);
				}
				 
				// if (!(M->rank())) std::cout << "\n\n xpos=" << rayvec.Get(1) << ", ypos=" << rayvec.Get(2) << "\n\n";		
				ntemp = interp2(conf,n_mat,rayvec.Get(1),rayvec.Get(2))/IMPACT_Heating::n_c;
				btemp = pow(interp2(conf,Bz_mat,rayvec.Get(1),rayvec.Get(2))/globalconsts::omega_p_nuei_n,2.0)/IMPACT_Heating::n_c;
	
				tempvG = sqrt(interp2(conf,eps_mat,rayvec.Get(1),rayvec.Get(2))-tempeps);
				
				tempvG = (tempvG)/(tempvG*tempvG+
					ntemp*((1-ntemp)*(1-ntemp)+ntemp*btemp)/(1-ntemp-btemp)/(1-ntemp-btemp));
	
	
				// std::cout << "xpos=" << rayvec.Get(1) << " ypos=" << rayvec.Get(2) <<std::endl;//
				// std::cout << "ntemp=" << ntemp << "\n\n";
	
				intensitycalc(conf,n_mat,rayvec,tempint,livemult/tempvG);//
				traceonestep(conf,rayvec,eps_mat,gradepsx,gradepsy,tempeps,ds,flag,firststep);
				
				++n;
	
				if (rayvec.Get(3)<0.0)
				{
					if (!turn)
					{
						turnx=rayvec.Get(1);
						turny=rayvec.Get(2);
						turn=1;
					}
					else
					{
						livemult*=pow((1.0-1e-3),sqrt((rayvec.Get(1)-turnx)*(rayvec.Get(1)-turnx)+(rayvec.Get(2)-turny)*(rayvec.Get(2)-turny)));
						turnx=rayvec.Get(1);
						turny=rayvec.Get(2);
					}
					// flag = 0;
					//discard = 1; std::cout << "uh oh \n";
				}			

				
	
				firststep=0;
			}
	
			if (!discard) updatelocalintensity(tempint,localint);
			tempint.setall(0.0);
	
			if ((!ifraydump) && !(M->rank()) && !ifnotdump)
			{
				raystream << "\r\n" << "\r\n";
				
			}
		}
	}

		// Share Intensity here:
	double localval;
	double gval;
	double numprocs = M->size();
	
	for (int i=1;i<=conf->Nx();++i)
	{
		for (int j=1;j<=conf->Ny();++j)
		{

			localval = localint.get(i,j);

			MPI_Reduce(&localval,&gval,1, MPI_DOUBLE, MPI_SUM,0,MPI::COMM_WORLD);
			MPI_Bcast(&gval,1,MPI_DOUBLE,0,MPI::COMM_WORLD);

			//std::cout<<"\n From rank: "<< M->rank() << " ----- localval = " << localval <<"------ gval = " << gval << "\n";
	//         if (!(M->rank())) {
			IMPACT_Heating::i_mat.set(i,j,gval/numprocs);
	 //       }
		}
			// if (!(M->rank())) {
	}
   //          	std::cout<<"\n From rank: "<< M->rank() << " --\n"; 
   //          	IMPACT_Heating::i_mat.Print();
   //          } 
   //          else
   //          {
   //          	std::cout<<"\n From rank: "<< M->rank() << " --\n"; 
   //          	IMPACT_Heating::i_mat.Print();
   //          }
	intensity_diffusion(conf);



		if (!ifnotdump && !(M->rank())) raystream.close();


}