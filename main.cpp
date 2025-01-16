using namespace std;

// Pour compilation
// g++ -O3 main.cpp

#include<iostream>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<cstring>

#include <time.h> 



// Include my lib :
#include "Lib/Vecteur.hpp"
#include "Lib/MatriceL.hpp"
#include "Lib/Meshes.hpp"

double donneeBord(Vecteur pos);

int main()
{
    
	int nThread = 16;
    	int PAUSE;
	double PI = atan(1.)*4;
	double epsilon = 1e-8; 
	cout << setprecision(20);    
	
	
			
//---------------------------------------------------------
	// Lecture maillage :
	Meshes myMeshes;
	myMeshes.readGMSH("Maillages/maillageSimple.msh");
	
	// Liste des points (avec les coordonnées) du maillages
	vector<Vecteur*> listeNoeuds = myMeshes.getNodes();
	
	// Nombre de points dans le maillage
	cout<<"Nombre de ddl : "<<listeNoeuds.size()<<endl;
	
	// Nombre de mailles (triangles)
	vector<int*> listeMaille = myMeshes.getMesh2D();
	cout<<"Nombre de mailles : "<<listeMaille.size()<<endl;
	
	// On récupère le nombre de noeuds (points)
	int nbFct = myMeshes.getNNodes();
	// On récupère le nombre de mailles 
	int nbMesh = myMeshes.getNMesh2D();
	
	// On récupère les mailes du bord du domaine :
	vector<int> lienBordVMesh;
	vector<int*> listMailleBord = myMeshes.getBoundary2D(1, 3, lienBordVMesh); 
//---------------------------------------------------------
	
	
//---------------------------------------------------------	
	// Déclaration des matrices :
	MatriceL Rig(nbFct);
	MatriceL Masse1DBord(nbFct);
//---------------------------------------------------------
	
		
//---------------------------------------------------------	
	// Assemblage :		
	cout<<"Assemblage matrices de masse et de rigidité"<<endl;
	// Pour numérotation locale :
	int num1, num2, num3;
	Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); // Coins triangle
	
	MatriceL T(2); // Matrice pour la transformation affine
	MatriceL InvT(2); // Matrice pour l'inverse de T
	
	// Matrices Elémentaires locales sur la maille K :
	MatriceL KEle(3);
	MatriceL MEle(3);
	
	Vecteur T2(3);
	Vecteur T3(3);
	
	// Boucle d'assemblage sur les mailles 2D :
	for (int m=0; m<nbMesh; m++) {
		if (m % 100 == 0 or true) 
			cout<<"Maille m= "<<m<<" / "<<listeMaille.size()<<endl;
		
		// On considere le triangle numero m :
		// Numero globale des coins du triangle :
		num1 = listeMaille[m][1];
		num2 = listeMaille[m][2];
		num3 = listeMaille[m][3];
		
		// p1,p2 et p3 contiennet les coordonnées des 3 sommets du triangle
		p1 = *listeNoeuds[num1];
		p2 = *listeNoeuds[num2];
		p3 = *listeNoeuds[num3];
		
		// Calcul de la matrice T formant la transformation :
		T2 = p2-p1;
		T3 = p3-p1;	
		
		// TODO
		T.set(0,0,T2(0));
		T.set(0,1,T3(0));
		T.set(1,0,T2(1));
		T.set(1,1,T3(1));
		
		// Calcul de déterminant de T :
		// TODO : 
		double detT = T(0,0) * T(1,1) - T(0,1) * T(1,0);
		
		
		// Calcul du volume de la maille K :
		double volume = detT*0.5;
		
		// Calcul de l'inverse de T :
		// TODO :
		InvT.set(0,0,T3(1) / detT);
		InvT.set(0,1,-T3(0) / detT); 
		InvT.set(1,0,-T2(1) / detT); 
		InvT.set(1,1,T2(0) / detT);
		
	
		// Matrice élémentaire :
		// > Rigidité : Nabla Phi_i . Nabla Phi_j
		double tmpx = InvT(0,0) + InvT(1,0);
		double tmpy = InvT(0,1) + InvT(1,1);
		// Completer ICI 
		KEle.set(0, 0, volume * (tmpx * tmpx + tmpy * tmpy)); // K(0,0)
		KEle.set(1, 1, volume * (InvT(0,0) * InvT(0,0) + InvT(0,1) * InvT(0,1))); // K(1,1)
		KEle.set(2, 2, volume * (InvT(1,0) * InvT(1,0) + InvT(1,1) * InvT(1,1))); // K(2,2) // Ok
//
		KEle.set(0, 1, volume * ((-1) * InvT(0,0) * tmpx + (-1) * InvT(0,1) * tmpy)); // K(0,1)
		KEle.set(0, 2, volume * ((-1) * InvT(1,0) * tmpx + (-1) * InvT(1,1) * tmpy)); // K(0,2) // Ok
//
		KEle.set(1, 0, volume * ((-1) * InvT(0,0) * tmpx + (-1) * InvT(0,1) * tmpy)); // K(1,0) = K(0,1)
		KEle.set(1, 2, volume * (InvT(0,0) * InvT(1,0) + InvT(0,1) * InvT(1,1)));      // K(1,2) // Ok
//
		KEle.set(2, 0, volume * ((-1) * InvT(1,0) * tmpx + (-1) * InvT(1,1) * tmpy)); // K(2,0) = K(0,2)
		KEle.set(2, 1, volume * (InvT(0,0) * InvT(1,0) + InvT(0,1) * InvT(1,1)));      // K(2,1)
		
		
		// Assemblage dans les matrices globales :
		double tmpC;
		for (int i=0;i<3;i++) {
			int indL = listeMaille[m][i+1];
			for (int j=0;j<3;j++) {
				int indC = listeMaille[m][j+1];
				
				// > Matrice Rigidite :
				tmpC = Rig(indL,indC);
				// Completer ICI :
				Rig.set(indL, indC, tmpC + KEle(i, j)); 
			}
		}
	}
	
	// // Boucle d'assemblage sur les mailles 2D :
	MatriceL MEle1D(2);
	for (int m=0; m<listMailleBord.size(); m++) {
		num1 = listMailleBord[m][1];
		num2 = listMailleBord[m][2];

		// p1,p2 contiennent les coordonnées des 2 extrémités du segment
		p1 = *listeNoeuds[num1];
		p2 = *listeNoeuds[num2];


		// Aide 
		double hh = (p1-p2).norme();
		MEle1D.set(0,0,hh/3); 
		MEle1D.set(0,1,hh/6);
		MEle1D.set(1,1,hh/3);
		MEle1D.set(1,0,hh/6);

		// Assemblage dans les matrices globales :
		double tmpC;
		for (int i=0;i<2;i++) {
			int indL = listMailleBord[m][i+1];
			for (int j=0;j<2;j++) {
				int indC = listMailleBord[m][j+1];

				// > Matrice Rigidite :
				tmpC = Masse1DBord(indL,indC);
				// Completer ICI :
				Masse1DBord.set(indL,indC,tmpC + MEle1D(i,j));
			}
		}
	}
	

	cout<<"Fin assemblage"<<endl;
//---------------------------------------------------------	
	

	
//---------------------------------------------------------		
	// On impose u = 0 à une position à l'intérieur du domaine 
	// Recherche d'un point dans l'intérieur du domaine
	int ind0 = 0;
	for(int i=0;i<nbFct;i++) {
		if(Masse1DBord(i,i) > 0) {
			ind0 = i;
			break;
		}
	}
	
	for(int i=0;i<nbFct;i++) {
		if(i != ind0) {
			Rig.set(ind0,i,0.0);
			Rig.set(i,ind0,0.0);
		}
	}
//---------------------------------------------------------
	
//---------------------------------------------------------
	// Calcul du terme source venant du bord:		
	// Déclaration du vecteur second membre :
	Vecteur scmb(nbFct);
	Vecteur F(nbFct);
	for(int i=0;i<nbFct;i++) {
		Vecteur pos = *listeNoeuds[i];
		F.set(i,donneeBord(pos));
	}
	scmb = Masse1DBord * F;	
//---------------------------------------------------------		

//---------------------------------------------------------		
// Résolution du probleme :
	Vecteur SOL(nbFct);
	FILE* fid;
	
	cout<<"Début résolution effective :"<<endl;
	SOL = Rig.solveIte(scmb,SOL*0.0,"CG","SSOR",1.47,5000, 4000, 1e-15,0);
	cout<<"Fin résolution"<<endl;
//---------------------------------------------------------	
	
//---------------------------------------------------------		
// Ecriture de la solution:
	cout<<"Ecriture du resultat :"<<endl;
	fid = fopen("resultat.vtk","w");
	
	fprintf(fid, "# vtk DataFile Version 3.6\n");
	fprintf(fid, "PolyDATA\n");
	fprintf(fid, "ASCII\n");
	fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

	fprintf(fid, "POINTS	%d	float\n",nbFct);
	// Liste Noeuds du maillage :
	for(int i=0;i<nbFct;i++) {
		fprintf(fid, "%lf	%lf	%lf\n",(*listeNoeuds[i])(0), (*listeNoeuds[i])(1), 0.0);
	}

	// Liste Triangles du maillage :
	fprintf(fid, "CELLS	%d	%d\n",int(listeMaille.size()), int(listeMaille.size())*4);

	for(int m=0;m<listeMaille.size();m++)
		fprintf(fid, "3	%d	%d	%d\n",listeMaille[m][1],listeMaille[m][2],listeMaille[m][3]);

	fprintf(fid, "CELL_TYPES	%d\n",int(listeMaille.size()));
	for(int m=0;m<listeMaille.size();m++) 
		fprintf(fid, "5\n");

	// Valeurs solution sur les différents noeuds du maillage :
	fprintf(fid, "POINT_DATA	%d\n",nbFct);
	fprintf(fid,"SCALARS solution float	1\n");
	fprintf(fid,"LOOKUP_TABLE default\n");

	for(int i=0;i<nbFct;i++) 
		fprintf(fid, "%lf\n",SOL(i)); 


	fclose(fid);
//---------------------------------------------------------			
		
    cout<<"FIN programme"<<endl;
	
}


double donneeBord(Vecteur pos) {
	double x = pos(0); double y = pos(1);
	double PI = 3.141592653;
	
	// Pour maillage simple :
	if(fabs(x) <= 1e-8 )
 		return -1.0;
 	if(fabs(x-20.0)<=1e-8)
 		return 1.0;
 	return 0.0;
}
