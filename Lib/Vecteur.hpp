#ifndef VECTEUR_HPP
#define VECTEUR_HPP

#include <math.h>

class Vecteur  
{
	private :
	int dim;
	//stockage type matrice creuse...
	vector<double> v;
	
	public :
	Vecteur(int Dim, vector<double> V)
	{
		dim = Dim;
		for(int i=0; i<V.size();i++)
		{
			v.push_back(V[i]);
		}
	}
	
	Vecteur()
	{
		dim= 0;
	}
	
	~Vecteur()
	{
		//cout<<"destr V"<<endl;
	}
	
	Vecteur(int Dim)
	{
		dim = Dim;
		for(int i=0; i<dim;i++)
		{
			v.push_back(0.0);
		}
		
	}
	
	double operator() (int l)
	{
		return v[l];
	}
	
	bool operator== (Vecteur V)
	{
		if(V.getDim() != dim)
			return false;
		bool res = true;
		int i = 0;
		while(res && i < dim)
		{
			if(v[i] != V(i))
				res = false;
			i++;
		}
		return res;
	}
	
	Vecteur operator+(Vecteur V)
	{
		Vecteur Vres(dim);
		for(int i=0;i<dim;i++)
		{
			Vres.set(i, v[i]+V(i));
		}
		return Vres;
	}
	
    void operator=(Vecteur V)
    {
        for(int i=0;i<dim;i++)
            v[i] = V(i);
    }
    
	Vecteur operator-(Vecteur V)
	{
		Vecteur Vres(dim);
		for(int i=0;i<dim;i++)
		{
			Vres.set(i, v[i]-V(i));
		}
		return Vres;
	}
	
	double operator*(Vecteur V)
	{
		double res = 0;
		for (int i=0;i<dim;i++)
			res += v[i]*V(i);
		return res;
	}
	
    
    Vecteur prodTaT(Vecteur V)
    {
        Vecteur res(dim);
        for(int i=0;i<dim;i++)
            res.set(i,v[i]*V(i));
        return res;	
    }
	
	
	Vecteur operator*(double scal)
	{
		Vecteur res(dim);
		for(int i=0;i<dim;i++)
			res.set(i,v[i]*scal);
		return res;	
	}
	
	Vecteur operator/(double scal)
	{
		Vecteur res(dim);
		for(int i=0;i<dim;i++)
			res.set(i,v[i]/scal);
		return res;	
	}
	
	void set(int l,double val)
	{
		v[l] = val;
	}
	
	int getDim()
	{
		return dim;
	}
	
	double norme()
	{
		double res = 0 ;
		for (int i=0; i<dim; i++) {
			res += v[i]*v[i];
		}
		return sqrt(res);
	}
    
    vector<double> toVector() {
        return v;
    }
};
#endif