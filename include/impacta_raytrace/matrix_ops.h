void epscalc(IMPACT_Matrix &eps_mat, IMPACT_Matrix &n_mat, IMPACT_Matrix &Bz_mat)
{
	int Nx,Ny;
	double ntemp,btemp,temp;

	n_mat.size(&Nx,&Ny);
	eps_mat.ChangeSize(Nx,Ny);
	eps_mat.setall(0.0);

	for (int i=1;i<=Nx;++i)
	{
		for (int j=1;j<=Ny;++j)
		{
			ntemp = n_mat.get(i,j)/IMPACT_Heating::n_c;
			btemp = pow(Bz_mat.get(i,j)/globalconsts::omega_p_nuei_n,2.0)/IMPACT_Heating::n_c;
	
			temp = (1.0-ntemp*(1.0-ntemp)/(1.0-ntemp-btemp)); // sqrt(epsilon)
	
			//if (temp<0)	temp=0;	
			//temp = sqrt(temp);

			eps_mat.set(i,j,temp);
		}
	}	
}
void gradcalc(IMPACT_StenOps *sten, IMPACT_Config *conf, IMPACT_Matrix &n_mat, IMPACT_Matrix &gradnx, IMPACT_Matrix &gradny)
{

	double dx, dy, dxL, dyL, dxH, dyH;
	int i, j;

	gradnx.ChangeSize(conf->Nx(),conf->Ny());
	gradny.ChangeSize(conf->Nx(),conf->Ny());
	
	// for (i=1;i<=conf->Nx();i++)
	// 	for (j=1;j<=conf->Ny();j++)
	// 		{gradnx.set(i,j,-0.01);
	// 					gradny.set(i,j,0.0);}

	dxL = conf->xpos(2)-conf->xpos(1);
	dyL = conf->ypos(2)-conf->ypos(1);
	dxH = conf->xpos(conf->Nx())-conf->xpos(conf->Nx()-1);
	dyH = conf->ypos(conf->Ny())-conf->ypos(conf->Ny()-1);

	gradnx.set(1,1,(n_mat.get(2,1)-n_mat.get(1,1))/dxL );  			// Bottom Left
	gradny.set(1,1,(n_mat.get(1,2)-n_mat.get(1,1))/dyL );			// Bottom Left
	
	gradnx.set(conf->Nx(),1,(n_mat.get(conf->Nx(),1)-n_mat.get(conf->Nx()-1,1))/dxH );		// Bottom Right
	gradny.set(conf->Nx(),1,(n_mat.get(conf->Nx(),2)-n_mat.get(conf->Nx(),1))/dyL );		// Bottom Right
	
	gradnx.set(1,conf->Ny(),(n_mat.get(2,conf->Ny())-n_mat.get(1,conf->Ny()))/dxL );  		// Top Left
	gradny.set(1,conf->Ny(),(n_mat.get(1,conf->Ny())-n_mat.get(1,conf->Ny()-1))/dyH );		// Top Left

	gradnx.set(conf->Nx(),conf->Ny(),(n_mat.get(conf->Nx(),conf->Ny())-n_mat.get(conf->Nx()-1,conf->Ny()))/dxH );	// Top Right
	gradny.set(conf->Nx(),conf->Ny(),(n_mat.get(conf->Nx(),conf->Ny())-n_mat.get(conf->Nx(),conf->Ny()-1))/dyH );	// Top Right

	// Left and Right Boundaries
	for (j=2;j<=conf->Ny()-1;j++)
	{
		
		gradnx.set(1,j,(n_mat.get(2,j)-n_mat.get(1,j))/dxL );  			// Left
		gradnx.set(conf->Nx(),j,(n_mat.get(conf->Nx(),j)-n_mat.get(conf->Nx()-1,j))/dxH );		// Right
		
		dy = conf->ypos(j+1) - conf->ypos(j-1);

		gradny.set(1,j,(n_mat.get(1,j+1)-n_mat.get(1,j-1))/dy );		// Left
		gradny.set(conf->Nx(),j,(n_mat.get(conf->Nx(),j+1)-n_mat.get(conf->Nx(),j-1))/dy );	// Right

		// std::cout << "\n\n VAL: " << n_mat.get(conf->Nx(),j+1) << ", " << n_mat.get(conf->Nx(),j-1) << ", " << (n_mat.get(Nx,j+1)-n_mat.get(Nx,j-1))/dy << "\n\n";
		// std::cout << "\n\n j= " << j << " and, gradny= " << gradny.get(Nx,j) << "\n\n";
	}

	// Top and Bottom Boundaries
	for (i=2;i<=conf->Nx()-1;i++)
	{
		dx = conf->xpos(i+1) - conf->xpos(i-1);
		gradnx.set(i,1,(n_mat.get(i+1,1)-n_mat.get(i-1,1))/dx );		// Bottom
		gradnx.set(i,conf->Ny(),(n_mat.get(i+1,conf->Ny())-n_mat.get(i-1,conf->Ny()))/dx );	// Top
		
		gradny.set(i,1,(n_mat.get(i,2)-n_mat.get(i,1))/dyL );				// Bottom
		gradny.set(i,conf->Ny(),(n_mat.get(i,conf->Ny())-n_mat.get(i,conf->Ny()-1))/dyH );			// Top
	}
	
	// The Rest
	for (i=2;i<=conf->Nx()-1;++i)
	{
		for (j=2;j<=conf->Ny()-1;++j)
		{
			dx = conf->xpos(i+1) - conf->xpos(i-1);
			dy = conf->ypos(j+1) - conf->ypos(j-1);
			gradnx.set(i,j,(n_mat.get(i+1,j)-n_mat.get(i-1,j))/dx);
			gradny.set(i,j,(n_mat.get(i,j+1)-n_mat.get(i,j-1))/dy);
		}
	}	
 	// gradnx.Print();
	// gradny.Print();			
}
double interp2(IMPACT_Config *conf, IMPACT_Matrix &mat, double xpos, double ypos)
{
	// Interpolates mat(x,y) given mat(xvec,yvec) using bilinear interpolation
	int lowx=0, uppx=conf->Nx()+1, lowy=0, uppy=conf->Ny()+1;
	int tempi; 
	double tempv,answer;
	double dux, dlx,duy,dly;
	bool check;
	//mat.size(&Nx,&Ny);
	
	// IMPACT_Stencil ts = (*sten->ddxi(&i,&j,&x1));
	// Find indices for x
	check = 0;
	
	tempi = (uppx-lowx)/2;
	while (!check)
	{
		tempv = conf->xpos(tempi);
		
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

	dux = fabs(xpos-conf->xpos(lowx))/(conf->xpos(uppx)-conf->xpos(lowx));
	dlx = 1.0-dux;

	// Find indices for y
	check = 0;
	
	tempi = (uppy-lowy)/2;
	while (!check)
	{
		tempv = conf->ypos(tempi);

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

	duy = (ypos-conf->ypos(lowy))/(conf->ypos(uppy)-conf->ypos(lowy));
	dly = 1.0-duy;
			

	if (lowx==0){ lowx++;}
	if (uppx==conf->Nx()+1) {uppx--;}
	if (lowy==0) {lowy++;}
	if (uppy==conf->Ny()+1) {uppy--;}

	// std::cout << "xpos=" << xpos << std::endl;
	// std::cout << "lowx=" << lowx << ", uppx=" << uppx << std::endl;
	// std::cout << "mat.get(lowx)=" << mat.get(lowx,uppy) << ", uppx=" << mat.get(uppx,uppy) << std::endl;
	// if (uppx > Nx || lowx < 1 || uppy > Ny || lowy < 1) {
	// 	std::cout << "Interpolating out of bounds" << std::endl;
	// std::cout << "xpos=" << xpos << ", ypos=" << ypos << std::endl;
	// std::cout << "lowx=" << lowx << ", uppx=" << uppx << std::endl;
	// std::cout << "lowy=" << lowy << ", uppy=" << uppy << std::endl;
	
	// std::cout<<" \n\n hello : \n\n"; 
	// Bilinear Interpolation



	answer = mat.get(lowx,lowy)*(dlx)*(dly)
		+mat.get(uppx,lowy)*(dux)*(dly)
		+mat.get(lowx,uppy)*(dlx)*(duy)
		+mat.get(uppx,uppy)*(dux)*(duy);

	// std::cout << "xpos=" << xpos << ", ypos=" << ypos << std::endl;
	// std::cout << "a: " << answer << std::endl;

	return answer;
}
void interp23(IMPACT_Config *conf, IMPACT_Matrix &mat, IMPACT_Matrix &mat1, IMPACT_Matrix &mat2,  double xpos, double ypos, double vals[])
{
	int lowx=0, uppx=conf->Nx()+1, lowy=0, uppy=conf->Ny()+1;
	int tempi; 
	double tempv;
	double dux, dlx,duy,dly;
	bool check;
	//mat.size(&Nx,&Ny);
	
	// Find indices for x
	check = 0;
	
	tempi = (uppx-lowx)/2;
	while (!check)
	{
		tempv = conf->xpos(tempi);
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

	dux = fabs(xpos-conf->xpos(lowx))/(conf->xpos(uppx)-conf->xpos(lowx));
	dlx = 1.0-dux;

	// Find indices for y
	check = 0;
	
	tempi = (uppy-lowy)/2;
	while (!check)
	{
		tempv = conf->ypos(tempi);
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

	duy = (ypos-conf->ypos(lowy))/(conf->ypos(uppy)-conf->ypos(lowy));
	dly = 1.0-duy;
			

	if (lowx==0){ lowx++;}
	if (uppx==conf->Nx()+1) {uppx--;}
	if (lowy==0) {lowy++;}
	if (uppy==conf->Ny()+1) {uppy--;}
	
	// if (uppx > Nx || lowx < 1 || uppy > Ny || lowy < 1) {
	// 	std::cout << "Interpolating out of bounds" << std::endl;
		// std::cout << "xpos=" << xpos << ", ypos=" << ypos << std::endl;
		// std::cout << "lowx=" << conf->xpos(lowx) << ", uppx=" << conf->xpos(uppx) << std::endl;
		// std::cout << "lowy=" << conf->ypos(lowy) << ", uppy=" << conf->ypos(uppy) << std::endl << std::endl;
	
	// Bilinear Interpolation
	vals[0] = mat.get(lowx,lowy)*(dlx)*(dly)
		+mat.get(uppx,lowy)*(dux)*(dly)
		+mat.get(lowx,uppy)*(dlx)*(duy)
		+mat.get(uppx,uppy)*(dux)*(duy);

	vals[1] = mat1.get(lowx,lowy)*(dlx)*(dly)
		+mat1.get(uppx,lowy)*(dux)*(dly)
		+mat1.get(lowx,uppy)*(dlx)*(duy)
		+mat1.get(uppx,uppy)*(dux)*(duy);

	vals[2] = mat2.get(lowx,lowy)*(dlx)*(dly)
		+mat2.get(uppx,lowy)*(dux)*(dly)
		+mat2.get(lowx,uppy)*(dlx)*(duy)
		+mat2.get(uppx,uppy)*(dux)*(duy);

	// std::cout << "xpos=" << xpos << ", ypos=" << ypos << std::endl;
	// std::cout << "0 = " << vals[0] << ", 1 = " << vals[1] << ", 2 = " << vals[2] << std::endl;

	
}
int findindex(IMPACT_Config *conf, IMPACT_Matrix &mat, double val,bool direction)
// int findindex(IMPACT_Matrix &mat, double val,bool direction)
{
	// Interpolates mat(x,y) given mat(xvec,yvec) using bilinear interpolation
	double maxloc;
	int i=1;
	bool flg=1;
	maxloc =1;
	//std::cout<<"HERE: " << i << "and"<< j << "\n";
	// Find indices for x
	if (direction)
	{
		while(i<=conf->Ny() && flg)
		{
			if (mat.get(floor(conf->Nx()/2.0)+1,i)<val)
			{
				maxloc = (val-mat.get(floor(conf->Nx()/2.0)+1,i-1))/(mat.get(floor(conf->Nx()/2.0)+1,i)-mat.get(floor(conf->Nx()/2.0)+1,i-1))
				*(conf->ypos(i)-conf->ypos(i-1))+conf->ypos(i-1);

				// maxloc=i;

				flg=0;

			}
			++i;
		}	
	}
	else
	{
		while(i<=conf->Nx() && flg)
		{
			if (mat.get(i,floor(conf->Ny()/2.0)+1)<val)
			{
				maxloc = (val-mat.get(i-1,floor(conf->Ny()/2.0)+1))/(mat.get(i,floor(conf->Ny()/2.0)+1)-mat.get(i-1,floor(conf->Ny()/2.0)+1))
				*(conf->xpos(i)-conf->xpos(i-1))+conf->xpos(i-1);
				// maxloc=i;
				flg=0;
			}
			++i;
		}	
	}
	return (maxloc);
}