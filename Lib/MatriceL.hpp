#ifndef MATRICEL_HPP
#define MATRICEL_HPP

#include <math.h>

class MatriceL
{
	private :
		int dimX;
	    int dimY;
		vector< vector<int> > indexL;
		vector< vector<double> > val;
        bool is_factored;
		bool is_ILU;
		
    
	public :
		MatriceL(int DimX, int DimY) {
			dimX = DimX;
			dimY = DimY;
			for (int i=0; i<dimX; i++) {
				vector<double> tmpC;
				val.push_back(tmpC);
				vector<int> tmpI;
				indexL.push_back(tmpI);
			}
			is_ILU = false;
		};
	
		MatriceL(int Dim) {
			dimX = Dim;
			dimY = Dim;
			for (int i=0; i<dimX; i++) {
				vector<double> tmpC;
				val.push_back(tmpC);
				vector<int> tmpI;
				indexL.push_back(tmpI);
			}
			is_ILU = false;
		}
	
		MatriceL() {
		
		};
		
		MatriceL(const MatriceL &M)
		{
			dimX = M.dimX;
			dimY = M.dimY;
			for(int i=0; i < dimX; i++){
				vector<double> tmpC;
				val.push_back(tmpC);
				vector<int> tmpI;
				indexL.push_back(tmpI);
				for(int j=0;j<(M.indexL[i]).size();j++) {
					val[i].push_back(M.val[i][j]);
					indexL[i].push_back(M.indexL[i][j]);
				}
			}
		}
	
		~MatriceL() {
			//cout<<"Dest M"<<endl;
		};
	
		double operator() (int l, int c) {
            if(l>=dimX || c >= dimY) {
				cout<<"probleme !!"<<endl;
                cout<<"l c dimX dimY "<<l<<"    "<<c<<" "<<dimX<<"  "<<dimY<<endl;
            }
			for(int j=0;j<indexL[l].size();j++)
				if(indexL[l][j] == c)
					return val[l][j];
			return 0.0;
		};
	
		void clear() {
			indexL.clear();
			val.clear();
			for (int i=0; i<dimX; i++) {
				vector<double> tmpC;
				val.push_back(tmpC);
				vector<int> tmpI;
				indexL.push_back(tmpI);
			}
		}
	
		Vecteur extrC(int c) {
			Vecteur vRes(dimX);
			for(int i=0;i<dimX;i++) {
				for(int j=0;j<indexL[i].size();j++) {
					if(indexL[i][j] == c) {
						vRes.set(i, val[i][j]);
						break;
					}
				}
			}
			return vRes;
		};
	
		Vecteur extrL(int l) {
			Vecteur vRes(dimY);
			for(int i=0;i<indexL[l].size();i++) {
				vRes.set(indexL[l][i], val[l][i]);
			}
			return vRes;
		};	
	
		Vecteur operator*(Vecteur V) {
			Vecteur Vres(dimX);
			for(int i=0;i<dimX;i++) {
				//Vres.set(i,complex(0,0));
				double res = 0;
				for(int j=0;j<indexL[i].size();j++) {
					res = res+val[i][j]*V(indexL[i][j]);
				}
				Vres.set(i,res);
			}
			return Vres;
		};
	
		void set(int L, int C, double Val) {
            if(L<0 || L>dimX || C>dimY || C<0) {
				cout<<"Erreur indice L"<<endl;
                cout<<"L C : "<<L<<" "<<C<<endl;
                cout<<"dimX dimY : "<<dimX<<"   "<<dimY<<endl;
            }
			else {
				bool found = false;
				for(int i=0;i<indexL[L].size();i++) {
					if(indexL[L][i] == C) {
                        if(fabs(Val) == 0) {
                            val[L].erase(val[L].begin()+i);
                            indexL[L].erase(indexL[L].begin()+i);
                            is_factored = false;
							is_ILU = false;
                            return;
                        }
                        else {
                            val[L][i] = Val;
                            found = true;
                            is_factored = false;
							is_ILU = false;
                            return;
                        }
					}
				}
				if(!found && fabs(Val)>0) {
					val[L].push_back(Val);
					indexL[L].push_back(C);
                    is_factored = false;
					is_ILU = false;
				}
			}
		};
	
		
		int getNbTermL(int L) {
			return indexL[L].size();
		}
		int getIndLC(int L,int C) {
			return indexL[L][C];
		}
	
		int getNnz() {
			int nnz = 0;
			for (int l=0; l<dimX; l++) {
				nnz += indexL[l].size();
			}
			return nnz;
		}
	
	
		int getDimX() {
			return dimX;
		}
		int getDimY() {
			return dimY;
		}
		
		// void convertC(double** A)
// 		{
// 			for(int i=0;i<dimX;i++)
// 				for(int j=0;j<dimY;j++)
// 					A[i][j] = 0.0;
// 			for(int i=0;i<dimX;i++)
// 				for(int j=0;j<indexL[i].size();j++)
// 					A[i][indexL[i][j]] = val[i][j];
//
// 		};
	
        MatriceL transpose()
        {
            MatriceL res(dimX,dimY);
            for (int l=0; l<dimX; l++) {
                int nbEleL = getNbTermL(l);
                for (int j=0; j<nbEleL; j++) {
                    int c = getIndLC(l,j);
                    res.set(c,l,val[l][j]);
                }
            }
            return res;
        }
	
        MatriceL operator*(double scal)
        {
            MatriceL res(dimX,dimY);
			//#pragma omp parallel for num_threads(10)
            for (int l=0; l<dimX; l++) {
                int nbEleL = getNbTermL(l);
                for (int j=0; j<nbEleL; j++) {
                    int c = getIndLC(l,j);
                    res.set(l,c,val[l][j]*scal);
                }
            }
            return res;
        }
	
		MatriceL operator+(MatriceL mat)
		{
			MatriceL res(dimX);
			for (int l=0; l<dimX; l++) {
				int nbEleL = getNbTermL(l);
				for (int j=0; j<nbEleL; j++) {
					int c = getIndLC(l,j);
                    double tmpC = res(l,c);
					res.set(l,c,val[l][j]+tmpC);
				}
				nbEleL = mat.getNbTermL(l);
				for (int j=0; j<nbEleL; j++) {
					int c = mat.getIndLC(l,j);
                    double tmpC = res(l,c);
					res.set(l,c,tmpC+mat(l,c));
				}
			}
			return res;	
		}
	
		MatriceL operator-(MatriceL mat)
		{
			MatriceL res(dimX);
			for (int l=0; l<dimX; l++) {
				int nbEleL = getNbTermL(l);
				for (int j=0; j<nbEleL; j++) {
					int c = getIndLC(l,j);
                    double tmpC = res(l,c);
                    res.set(l,c,val[l][j]+tmpC);
				}
				nbEleL = mat.getNbTermL(l);
				for (int j=0; j<nbEleL; j++) {
					int c = mat.getIndLC(l,j);
                    double tmpC = res(l,c);
                    res.set(l,c,tmpC-mat(l,c));
				}
			}
			return res;	
		}
    
  
        Vecteur solveIte(Vecteur scmb, Vecteur x0, char* type, char* choix, double omega, int iteMax, int restart, double tol, int affich) {
			// Methode Iterative :
			int size = dimX;
			
			Vecteur r(size);
			Vecteur Sol(size);
			Sol = x0;
			
			if (scmb.norme() < 1e-15)
				return Sol*0;
			
			if (strcmp(type,"KrylovSym") == 0) {
				
				r = (this->operator*)(x0) - scmb;// A x0 - b
				
				double r0 = (r*r);
				
				Vecteur v1(size); Vecteur v2(size); Vecteur v3(size);
				Vecteur w1(size); Vecteur w2(size); Vecteur w3(size);
				Vecteur w(size);
				
				Vecteur Cr(dimX); Vecteur Cw(dimX);
				//Cr = r;//precond(r, choix, omega);
				Cr = precond(r, choix, omega);
				//v2 = Cr / sqrt((Cr*Cr));
				v2 = Cr / sqrt(Cr*r);
				//v2 = Cr;// Gv2Tilde
					
				v1 = v2*0.; v3 = v2*0.;
				w1 = v1*0.; w2 = v1*0.;
				
				//double theta1 = 0.; double theta2 = 0.; double theta3 = 0.;
				double sintheta1 = 0.; double sintheta2 = 0.; double sintheta3 = 0.;
				double costheta1 = 1.; double costheta2 = 1.; double costheta3 = 1.;
				
				int ite = 0;
				//double y1 = sqrt(r0); double y2 = 0.;
				double y1 = sqrt(r*Cr); double y2 = 0.;
				
				while((r*r) > tol*r0 and iteMax > ite ){
					if (affich == 1) {
						cout<<"Ite : "<<ite<<" / "<<iteMax<<endl;
						cout<<"Residu courant : "<<(r*r)/r0<<endl;
					}
					if (ite % restart == 0 and ite >0) {
						cout<<"Restart !"<<endl;
						Cr = precond(r, choix, omega);
						//v2 = Cr / sqrt((Cr*Cr));
						v2 = Cr / sqrt(Cr*r);
						v1 = v2*0.; v3 = v2*0.;
						w1 = v1*0.; w2 = v1*0.;
						
						sintheta1 = 0.; sintheta2 = 0.; sintheta3 = 0.;
						costheta1 = 1.; costheta2 = 1.; costheta3 = 1.;
						y1 = sqrt(r*Cr); y2 = 0.;
					}
					
					w = (this->operator*)(v2);
					Cw = precond(w, choix, omega);// GwTilde
					
					double h12 = (w*v1);
					double h22 = (w*v2);
					
			        Cw = Cw - v1*h12 - v2*h22;
			        double h23 = sqrt((Cw * w));
			        v3 = Cw / h23;
        
			        // Decomposition QR
			        //double r13 = sin(theta1)*h12;
			        double r13 = sintheta1*h12;
			        //double r23 = cos(theta1)*h12;
			        double r23 = costheta1*h12;
			        double r33 = h22;
			        //double tr23 = cos(theta2)*r23+sin(theta2)*r33;
			        double tr23 = costheta2*r23+sintheta2*r33;
			        //r33 = -sin(theta2)*r23+cos(theta2)*r33;
			        r33 = -sintheta2*r23+costheta2*r33;
			        r23 = tr23;
			        
					//theta3 = atan(h23 / r33);
					costheta3 = 1.0 / sqrt(1.0 + (h23 / r33) * (h23 / r33));
					sintheta3 = sqrt(1.0-costheta3*costheta3);
			        //double tr33 = cos(theta3)*r33 + sin(theta3)*h23;
			        double tr33 = costheta3*r33 + sintheta3*h23;
			        r33 = tr33;
     			   
			        //double ty1 = cos(theta3)*y1;
			        double ty1 = costheta3*y1;
			        //y2 = -sin(theta3)*y1;
			        y2 = -sintheta3*y1;
			        y1 = ty1;
        
			        w3 = (v2 - w1*r13 - w2*r23)/r33;
			        Sol = Sol-w3*y1;
        
			        r = (this->operator*)(Sol) - scmb;
        
			        // Iteration suivante :
			        v1 = v2;
			        v2 = v3;
        
			        w1 = w2;
			        w2 = w3;
        
			        //theta1 = theta2;
			        //theta2 = theta3;
			        sintheta1 = sintheta2;
			        sintheta2 = sintheta3;
					
			        costheta1 = costheta2;
			        costheta2 = costheta3;
						
			        y1 = y2;
					
					ite = ite + 1;
				}
				cout<<"Ite : "<<ite<<" / "<<iteMax<<endl;
				cout<<"Residu courant : "<<(r*r)/r0<<endl;
            }
			
			if (strcmp(type,"CG") == 0) {
				cout<<"ONK !!"<<endl;
				Vecteur w(size);
				Vecteur Aw(size);
				Vecteur Cr(size);
				int ite = 0;
				
				// Init. :
				
				r = (this->operator*)(x0) - scmb; // A*x0 - scmb
				w = precond(r, choix, omega); //Vecteur x, char* choix, double omega
				//w = r;
				double r0 = (r*r);
				while((r*r) > tol*r0 and iteMax > ite ){
		
					if (ite % restart == 0 and ite >0) {
						cout<<"Restart !"<<endl;
						w = r*1;
						w = precond(r, choix, omega);
					}
		
					if (affich == 1) {
						cout<<"Ite : "<<ite<<" / "<<iteMax<<endl;
						cout<<"Residu courant : "<<(r*r)/r0<<endl;
					}
		
					Aw = (this->operator*)(w);
		        	double z = ( (r*w) / (Aw*w) );
		        	Sol = Sol - w*z;
		        	r = r - Aw*z;
					Cr = precond(r, choix, omega);
		       	 	double l = ( (Cr*Aw) / (Aw*w) );
		        	w = Cr - w*l;
		
					ite = ite + 1;
				}
				cout<<"Ite : "<<ite<<" / "<<iteMax<<endl;
				cout<<"Residu courant : "<<(r*r)/r0<<endl;
            }
			
			if (strcmp(type,"GMRES") == 0) {
				
				Vecteur w(size); Vecteur y(size);
				vector<Vecteur> VV; 
				MatriceL H(size); MatriceL R(size);
				vector<double> list_theta;
				
				r = (this->operator*)(x0) - scmb;// A x0 - b
				r = precond(r, choix, omega);
				double r0 = sqrt(r*r);
				
				// Initialisation vecteur Vi
				VV.push_back(r / r0);
				y.set(0,r0);
				
				// Initialisation matrice H
				w = (this->operator*)(VV[0]);
				w = precond(w, choix, omega);
				
				H.set(0,0,w*VV[0]); w = w - VV[0]*H(0,0);
				double normw = sqrt(w*w);
				VV.push_back(w / normw); // On met un 2ème vecteur dans la base
				H.set(1,0,normw);
				
				// Initialisation matrice R
				R.set(0,0,H(0,0)); R.set(1,0,H(1,0));
				double theta = atan(R(1,0) / R(0,0));
				list_theta.push_back(theta);
				double tmp1, tmp2;
				tmp1 = cos(theta)*R(0,0) + sin(theta)*R(1,0);
				tmp2 = -sin(theta)*R(0,0) + cos(theta)*R(1,0);
				R.set(0,0,tmp1); R.set(1,0,tmp2);
				
				// Mise à jour 
				tmp1 = cos(theta)*y(0) + sin(theta)*y(1);
				tmp2 = -sin(theta)*y(0) + cos(theta)*y(1);
				y.set(0,tmp1); y.set(1,tmp2); 
				
				cout<<" Err ini : "<<fabs(y(1))<<endl;
				int p = 0;
				while(fabs(y(p+1)) > tol*r0 and iteMax > p) {
					p = p+1;
					if (affich == 1) {
						cout<<"Ite : "<<p<<" / "<<iteMax<<endl;
						cout<<"Residu courant : "<<fabs(y(p))/r0<<endl; // <<" "<<tol
					}
					w = (this->operator*)(VV[p]);
					w = precond(w, choix, omega);
					
					for (int i=0;i<=p; i++){
						double scal = w*VV[i];
						H.set(i,p,scal);
						R.set(i,p,scal);
						w = w - VV[i]*scal;
					}
					normw = sqrt(w*w);
					VV.push_back(w / normw);
					
					// Verif orthonormaux
					// for(int l=0;l<=p;l++){
// 						for (int c=0;c<=p;c++)
// 							cout<<l<<" "<<c<<" "<<VV[l]*VV[c]<<endl;
// 					}
// 					int PAUSE;
// 					cin>>PAUSE;
					H.set(p+1,p,normw);
					R.set(p+1,p,normw);
					
					
					for (int i=0;i<=p; i++){
						if (i==p){
							theta = atan(R(p+1,p) / R(p,p));
							list_theta.push_back(theta);
						}
						theta = list_theta[i];
						tmp1 = cos(theta)*R(i,p) + sin(theta)*R(i+1,p);
						tmp2 = -sin(theta)*R(i,p) + cos(theta)*R(i+1,p);
						R.set(i,p,tmp1); R.set(i+1,p,tmp2);
					}
					tmp1 = cos(theta)*y(p) + sin(theta)*y(p+1);
					tmp2 = -sin(theta)*y(p);
					y.set(p,tmp1); y.set(p+1,tmp2);
				}
				//cout<<"FIN "<<endl;
				Vecteur zz(p+1);
				// Remontee :
				for(int i=p;i>=0;i--){
					double tmp = y(i);
					for (int j=p;j>i;j--){
						tmp = tmp - R(i,j)*zz(j);
					}
					zz.set(i,tmp / R(i,i));
				}
				
				// Sol contient deja x0 !!!
				for (int i=0;i<=p;i++){
					Sol = Sol - VV[i]*zz(i);
				}
				
				cout<<"Ite : "<<p<<" / "<<iteMax<<endl;
				cout<<"Residu courant : "<<fabs(y(p+1))/r0<<endl;
            }
			
            return Sol;
        }
		
		Vecteur precond(Vecteur x, char* choix, double omega) {
			
			Vecteur Vres(dimX);
			Vres = x;
			
			if (strcmp(choix,"D") == 0) {
				for (int i=0;i<dimX;i++) {
					for (int j=0;j<indexL[i].size();j++) {
						if (indexL[i][j] == i and fabs(val[i][j]) != 0) {
							Vres.set(i, x(i)/val[i][j]);
						} 
					}
				}
			}
			if (strcmp(choix,"SSOR") == 0) {
				double tmpC = 0.0; 
				vector<double> Diag(dimX); 
				int ll;
				
		        for (int l=dimX-1;l>=0;l--){
					ll = -1;
					tmpC = Vres(l);
		        	for (int c=0;c<indexL[l].size();c++) {
		        		if(indexL[l][c] > l) {
							tmpC = tmpC - val[l][c]*Vres(indexL[l][c])*omega;
		        			//Vres.set(l,Vres(l) - val[l][c]*x(indexL[l][c])*omega );
		        		}
						if (indexL[l][c] == l) {
							ll = c;
						}
		        	}
					//if (ll != -1) {
					Diag[l] = val[l][ll];
						//Diag.push_back(val[l][ll]);
					Vres.set(l,tmpC / val[l][ll]);
						//Vres.set(l,Vres(l) / val[l][ll]);
						//}
					//else {
					//	cout<<"PRBL !!"<<endl; int PAUSE;
					//	cin>>PAUSE;
						//}
		        }
				
				for (int i=0;i<dimX;i++) {
					Vres.set(i, Vres(i)*Diag[i]);
					//Vres.set(i, Vres(i)*val[i][lDiag[i]]);
					/*
					for (int j=0;j<indexL[i].size();j++) {
						if (indexL[i][j] == i) {
							Vres.set(i, Vres(i)*val[i][j]);
						} 
					}
					*/
				}
		        for (int l=0;l<dimX;l++){
					//ll = -1;
					tmpC = Vres(l);
		        	for (int c=0;c<indexL[l].size();c++) {
		        		if(indexL[l][c] < l) {
		        			tmpC = tmpC - val[l][c]*Vres(indexL[l][c])*omega;
							//Vres.set(l,Vres(l) - val[l][c]*Vres(indexL[l][c])*omega );
		        		}
						//if (indexL[l][c] == l) {
						//	ll = c;
						//}
		        	}
					//if (lDiag[l] != -1) {
					//if (ll != -1) {
						//Vres.set(l,Vres(l) / val[l][ll]);
					Vres.set(l,tmpC / Diag[l]);
					//}
					//else {
					//	cout<<"PRBL !!"<<endl; int PAUSE;
					//	cin>>PAUSE;
					//}
		        }	
			}
			
			if (strcmp(choix,"ILU") == 0) {
				// if (!is_ILU) {
// 					for (int i=0; i<dimX;i++) {
// 						vector<double> tmpD;
// 						ILUval.push_back(tmpD);
// 						vector<int> tmpI;
// 						ILUindex.push_back(tmpD);
//
// 					}
// 				}
			}
			
			return Vres;
		}
    
		void setToIdentity() {
			this->clear();
			if (dimX != dimY) {
				cout<<"Matrice non carree !! Pb "<<endl;
				return;
			}
				
			for (int i=0;i<dimX;i++) 
				this->set(i,i,1.0);
		}
	
};
#endif
