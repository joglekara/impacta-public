void crossVecs(IMPACT_Vector &vec1, IMPACT_Vector &vec2, IMPACT_Vector &vec_out)
{
//	double temp11 = vec1.Get(1);
//	double temp12 = vec1.Get(2);
//	double temp13 = vec1.Get(3);
//	
//	double temp21 = vec2.Get(1);
//	double temp22 = vec2.Get(2);
//	double temp23 = vec2.Get(3);
	
	double temp1, temp2, temp3;
	temp1 = vec1.Get(2)*vec2.Get(3)-vec1.Get(3)*vec2.Get(2);
	temp2 = vec1.Get(3)*vec2.Get(1)-vec1.Get(1)*vec2.Get(3);
	temp3 = vec1.Get(1)*vec2.Get(2)-vec1.Get(2)*vec2.Get(1); 
	
	// vec_out.Set(1,vec1.Get(2)*vec2.Get(3)-vec1.Get(3)*vec2.Get(2));
	// vec_out.Set(2,vec1.Get(3)*vec2.Get(1)-vec1.Get(1)*vec2.Get(3));
	// vec_out.Set(3,vec1.Get(1)*vec2.Get(2)-vec1.Get(2)*vec2.Get(1)); 
	
//	vec_out.Set(1,temp12*temp23-temp13*temp22);
//	vec_out.Set(2,temp13*temp21-temp11*temp23);
//	vec_out.Set(3,temp11*temp22-temp12*temp21);

	vec_out.Set(1,temp1);
	vec_out.Set(2,temp2);
	vec_out.Set(3,temp3);

}

double dotVecs(IMPACT_Vector &vec1, IMPACT_Vector &vec2)
{
	double temp=0.0;
	for (int i=1;i<=vec1.length();i++)
	{		
		temp+=(vec1.Get(i)*vec2.Get(i));
	}
	return temp;
}

void addVecs(IMPACT_Vector &vec1, IMPACT_Vector &vec2, IMPACT_Vector &vec_out)
{
	//IMPACT_Vector vec_out(vec1.length());
	for (int i=1;i<=vec1.length();i++)
	{		
		vec_out.Set(i,vec1.Get(i)+vec2.Get(i));
	}
	//return vec_out.GetVec();
}

void elementMultVecs(IMPACT_Vector &vec1, IMPACT_Vector &vec2, IMPACT_Vector &vec_out)
{	
	for (int i=1;i<=vec1.length();i++)
	{		
		vec_out.Set(i,vec1.Get(i)*vec2.Get(i));
	}	
}

void elementDivVecs(IMPACT_Vector &vec1, IMPACT_Vector &vec2, IMPACT_Vector &vec_out)
{
	for (int i=1;i<=vec1.length();i++)
	{		
		vec_out.Set(i,vec1.Get(i)/vec2.Get(i));
	}	
}

void specialVstep(IMPACT_Vector &omega_h, double ds, IMPACT_Vector &vec_out)
{
//	omega_h.MultiplyAll(ds/2);
	
	vec_out.Set(1,2/(1.0+ds/2.0*omega_h.Get(1))/(1.0+ds/2.0*omega_h.Get(1)));
	vec_out.Set(2,2/(1.0+ds/2.0*omega_h.Get(2))/(1.0+ds/2.0*omega_h.Get(2)));
	vec_out.Set(3,2/(1.0+ds/2.0*omega_h.Get(3))/(1.0+ds/2.0*omega_h.Get(3)));	
}

