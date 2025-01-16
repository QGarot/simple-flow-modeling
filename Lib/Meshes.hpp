#ifndef Meshes_HPP
#define Meshes_HPP

class Meshes
{
	private :
	
	protected :
		vector<int*> Mesh3D;
		vector<int*> Mesh2D;
		vector<int*> Mesh1D;
		vector<int> Mesh0D;
		vector<int> Mesh3Dtag;
		vector<int> Mesh2Dtag;
		vector<int> Mesh1Dtag;
		vector<int> Mesh0Dtag;
		
		vector<int> lienM2DPkP1;
		
		vector<Vecteur*> Nodes;
		
		double eleSize2D;
		double eleSize3D;
		
		bool isConnec2D;
		MatriceL* mesh1Connect2D;
		MatriceL* mesh2Connect2D;
		
		
		bool isIn(Vecteur pt, Vecteur a, Vecteur b){
			Vecteur T = (b-a);
			Vecteur U = (pt-a);
			
			Vecteur cross(3);
			cross.set(0, T(1)*U(2)- T(2)*U(1));
			cross.set(1,-T(0)*U(2)+ T(2)*U(0));
			cross.set(2, T(0)*U(1)- T(1)*U(0));
			
			if(cross.norme() > 1e-10)
				return false;
			// if(T*U < 0)
// 				return false;
// 			if(T.norme() < U.norme())
// 				return false;
			return true;
		}
		
	public :
		Meshes() {
			eleSize2D = -1;
			isConnec2D = false;
		}
		
		~Meshes() {
			int size = Mesh3D.size();
			for (int i=0;i<size;i++){
				delete(Mesh3D[i]);
			}
			size = Mesh2D.size();
			for (int i=0;i<size;i++){
				delete(Mesh2D[i]);
			}
			size = Mesh1D.size();
			for (int i=0;i<size;i++){
				delete(Mesh1D[i]);
			}
			size = Nodes.size();
			for (int i=0;i<size;i++){
				delete(Nodes[i]);
			}
		}
		
		int getNNodes() {
			return Nodes.size();
		}
		
		int getNMesh3D() {
			return Mesh3D.size();
		}	
		
		int getNMesh2D() {
			return Mesh2D.size();
		}
		
		int getNMesh1D() {
			return Mesh1D.size();
		}
		
		int getNMesh0D() {
			return Mesh0D.size();
		}
		
		
		vector<int*> getMesh3D() {
			return Mesh3D;
		}	
		
		void  setMesh3D(vector<int*> newMesh3D) {
			Mesh3D = newMesh3D;
		}
		
		vector<int*>  getMesh2D() {
			return Mesh2D;
		}
		
		void  setMesh2D(vector<int*> newMesh2D) {
			Mesh2D = newMesh2D;
		}
		
		vector<int*>  getMesh1D() {
			return Mesh1D;
		}
		
		vector<int>  getMesh0D() {
			return Mesh0D;
		}
		
		vector<int> getMesh3Dtag() {
			return Mesh3Dtag;
		}	
		
		vector<int>  getMesh2Dtag() {
			return Mesh2Dtag;
		}
		
		vector<int>  getMesh1Dtag() {
			return Mesh1Dtag;
		}
		
		vector<int>  getMesh0Dtag() {
			return Mesh0Dtag;
		}
		
		vector<Vecteur*> getNodes() {
			return Nodes;
		}
		
		double getEleSize2D() {
			eleSize2D = -1;
			if(eleSize2D == -1){
				Vecteur p1(3); Vecteur p2(3); Vecteur p3(3);
				for(int m=0;m<Mesh2D.size();m++) {
					int num1 = Mesh2D[m][1];
					int num2 = Mesh2D[m][2];
					int num3 = Mesh2D[m][3];
					
					p1 = *(Nodes[num1]);
					p2 = *(Nodes[num2]);
					p3 = *(Nodes[num3]);
					
					Vecteur T1 = p2-p1;
					Vecteur T2 = p3-p1;
					Vecteur T3 = p3-p2;
					
					double vol = fabs(T1(0)*T2(1) - T1(1)*T2(0));
					if(vol < 1e-10){
						cout<<"PB vol null "<<vol<<endl;
						int PAUSE;
						cin>>PAUSE;
					}
					
					double tmp = (T1.norme() * T2.norme() * T3.norme())/vol;
					
					if(tmp > eleSize2D)
						eleSize2D = tmp;
					
				}
				
			}
			//cout<<eleSize2D<<endl;
			return eleSize2D;
		}
		
		double getEleSize3D() {
			eleSize3D = -1;
			MatriceL A(3);
			if(eleSize3D == -1){
				Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); Vecteur p4(3);
				for(int m=0;m<Mesh3D.size();m++) {
					int num1 = Mesh3D[m][1];
					int num2 = Mesh3D[m][2];
					int num3 = Mesh3D[m][3];
					int num4 = Mesh3D[m][4];
					
					p1 = *(Nodes[num1]);
					p2 = *(Nodes[num2]);
					p3 = *(Nodes[num3]);
					p4 = *(Nodes[num4]);
					
					Vecteur T1 = p2-p1;
					Vecteur T2 = p3-p1;
					Vecteur T3 = p4-p1;
					
					Vecteur T4 = p4-p2;
					Vecteur T5 = p4-p3;
					Vecteur T6 = p3-p2;
					
					A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
					A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
					A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
	
	
					double vol = fabs( ((A(0,0)*A(1,1)-A(0,1)*A(1,0))*A(2,2)-(A(0,0)*A(2,1)-A(0,1)*A(2,0))*A(1,2)+(A(1,0)*A(2,1)-A(1,1)*A(2,0))*A(0,2)));
					
					
					double tmp = vol / ( (T1.norme() * T2.norme()) +  (T1.norme()*T3.norme()) + (T2.norme()*T3.norme()));
					
					if(tmp > eleSize3D)
						eleSize3D = tmp;
					
				}
				
			}
			//cout<<eleSize2D<<endl;
			return eleSize3D;
		}
		
		void readGMSH(char* name){
			FILE* fid;
			fid = fopen (name,"r");
			
			cout<<"Starting reading the mesh: "<<endl;
			MatriceL A(3);
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); Vecteur p4(3);
			Vecteur T1(3); Vecteur T2(3); Vecteur T3(3); Vecteur T4(3);
			// Nodes:
			char s[1000];
			while (!(s[0] == '$' && s[1] =='N' && s[2] == 'o' && s[3] == 'd' && s[4] == 'e' && s[5] == 's')) {
				fscanf(fid,"%s",s);
			}
			int nNodes =0;
			int nBlock =0;
			int minTag = 0;
			int maxTag = 0;
			//fscanf(fid,"%d	%d	%d	%d",&nBlock,&nNodes,&minTag,&maxTag);
			fscanf(fid,"%d",&nBlock);
			fscanf(fid,"%d",&nNodes);
			fscanf(fid,"%d",&minTag);
			fscanf(fid,"%d",&maxTag);
			cout<<"Number of nodes: "<<nNodes<<endl;
			
			int cpt = 0;
			for(int b=0;b<nBlock;b++) {
				int entDim = 0; int tagEnt = 0; int param = 0; int nNodeB = 0;
				//fscanf(fid,"%d	%d	%d	%d",&entDim,&tagEnt,&param,&nNodeB);
				fscanf(fid,"%d",&entDim);
				fscanf(fid,"%d",&tagEnt);
				fscanf(fid,"%d",&param);
				fscanf(fid,"%d",&nNodeB);
				
				for(int i=0;i<nNodeB;i++){
					int num;
					fscanf(fid,"%d",&num);
				}
				for(int i=0;i<nNodeB;i++){
					double x,y,z;
					fscanf(fid,"%lf",&x);
					fscanf(fid,"%lf",&y);
					fscanf(fid,"%lf",&z);
					
					Vecteur* p = new Vecteur(3);
					p->set(0,x);
					p->set(1,y);
					p->set(2,z);
					Nodes.push_back(p);
				}
				
			}
			
			/*
			for (int i=0; i<nNodes; i++) {
				Vecteur* p = new Vecteur(3);
				double x,y,z;
				int num;
				fscanf(fid,"%d",&num);
				fscanf(fid,"%lf",&x);
				fscanf(fid,"%lf",&y);
				fscanf(fid,"%lf",&z);
				p->set(0,x);
				p->set(1,y);
				p->set(2,z);
				Nodes.push_back(p);
			}
			*/
			
			// Elements reading:
			fscanf(fid,"%s",s);
			while (!(s[0] == '$' && s[1] =='E' && s[2] == 'l' && s[3] == 'e' && s[4] == 'm' && s[5] == 'e')) {
				fscanf(fid,"%s",s);
			}
			int nElems = 0;
			int nElemB = 0;
			int minElemTag = 0;
			int maxElemTag = 0;
			//fscanf(fid,"%d",&nElems);
			//fscanf(fid,"%d	%d	%d	%d",&nElemB,&nElems,&minElemTag,&maxElemTag);
			fscanf(fid,"%d",&nElemB);
			fscanf(fid,"%d",&nElems);
			fscanf(fid,"%d",&minElemTag);
			fscanf(fid,"%d",&maxElemTag);
			
			cout<<"Number of Elements (Tot): "<<nElems<<endl;
			
			int numEle, typeEle, nbTag;
			int num1, num2, num3, num4;
			int num5, num6, num7, num8;
			int num9, num10, num11, num12;
			int num13, num14, num15, num16;
			int num17, num18, num19, num20;
			int tag;
			
			for(int b=0;b<nElemB;b++){
				int entDim; int entTag; int entTyp; int nbEleEnt;
				//fscanf(fid,"%d	%d	%d	%d",&entDim,&entTag,&entTyp,&nbEleEnt);
				fscanf(fid,"%d",&entDim);
				fscanf(fid,"%d",&entTag);
				fscanf(fid,"%d",&entTyp);
				fscanf(fid,"%d",&nbEleEnt);
				
				typeEle = entTyp;
				tag = entTag;
				for(int i=0;i<nbEleEnt;i++){
					fscanf(fid,"%d",&numEle);
					// Tetrahedron 3D
					if(typeEle == 4 or typeEle == 11 or typeEle == 29) {
			
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						fscanf(fid,"%d",&num4);
						num1--;
						num2--;
						num3--;
						num4--;
					
						// Check if null volume tetrahedron :
						p1 = *(Nodes[num1]);
						p2 = *(Nodes[num2]);
						p3 = *(Nodes[num3]);
						p4 = *(Nodes[num4]);
		
						T2 = p2-p1;
						T3 = p3-p1;	
						T4 = p4-p1;
		
						// Choice s.t. p1 <-> (0,0,0) && p2 <-> (1,0,0) & p3 <-> (0,1,0) & p4 <-> (0,0,1) 
						A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
						A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
						A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
		
		
						double volume = fabs( ((A(0,0)*A(1,1)-A(0,1)*A(1,0))*A(2,2)-(A(0,0)*A(2,1)-A(0,1)*A(2,0))*A(1,2)+(A(1,0)*A(2,1)-A(1,1)*A(2,0))*A(0,2)));
						if (volume > 1e-10) {
							int* tab ;
							if (typeEle == 4) {
								tab = new int[5];
								tab[0] = 4; // Element Type (3D -> tetrahedron)
								tab[1] = num1; tab[2] = num2;
								tab[3] = num3; tab[4] = num4;
							}
							if (typeEle == 11) {
								fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
								fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
								num5--; num6--; num7--; num8--; num9--; num10--;
							
								tab = new int[11];
								tab[0] = 11; // Element Type (3D -> tetrahedron o2)
								tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
								tab[5] = num5; tab[6] = num6; tab[7] = num7; 
								tab[8] = num8; tab[9] = num9; tab[10] = num10;
							}
							if (typeEle == 29) {
								fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
								fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
								fscanf(fid,"%d",&num11); fscanf(fid,"%d",&num12); fscanf(fid,"%d",&num13);
								fscanf(fid,"%d",&num14); fscanf(fid,"%d",&num15); fscanf(fid,"%d",&num16);
								fscanf(fid,"%d",&num17); fscanf(fid,"%d",&num18); fscanf(fid,"%d",&num19); fscanf(fid,"%d",&num20);
								num5--; num6--; num7--; num8--; num9--; num10--;
								num11--; num12--; num13--; num14--; num15--; num16--;
								num17--; num18--; num19--; num20--;
							
								tab = new int[21];
								tab[0] = 29; // Element Type (3D -> tetrahedron o3)
								tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
								tab[5] = num5; tab[6] = num6; tab[7] = num7; tab[8] = num8; 
								tab[9] = num9; tab[10] = num10; tab[11] = num11; tab[12] = num12;
								tab[13] = num13; tab[14] = num14; tab[15] = num15; tab[16] = num16;
								tab[17] = num17; tab[18] = num18; tab[19] = num19; tab[20] = num20;
							}
						
						
						
							Mesh3D.push_back(tab);
							Mesh3Dtag.push_back(tag);
						}
						else{
							int PAUSE;
							cout<<"Null tetrahedron volume !!!"<<endl;
							cout<<i<<"	"<<typeEle<<endl;
							cin>>PAUSE;
						}	
					}
				
				
					// square 2D
					if(typeEle == 3) {
						//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
						fscanf(fid,"%d",&num1);fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);fscanf(fid,"%d",&num4);
					
						num1--;
						num2--;
						num3--;
						num4--;
					
						int* tab = new int[5];
						tab[0] = 3; // Element Type (2D -> square o1)
						tab[1] = num1; tab[2] = num2;
						tab[3] = num3;  tab[4] = num4;
					
						Mesh2D.push_back(tab);
						Mesh2Dtag.push_back(tag);
					}
			
					// Triangle 2D
					if (typeEle == 2 or typeEle == 9 or typeEle == 21) {
						//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						num1--;
						num2--;
						num3--;
					
						int* tab;
						if (typeEle == 2) {
							 tab = new int[4];
							 tab[0] = 3; // Element Type (2D -> triangle)
							 tab[1] = num1; tab[2] = num2; tab[3] = num3; 
						}
						if (typeEle == 9) {
							fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
							num4--; num5--; num6--;
						
							tab = new int[7];
							tab[0] = 6; // Element Type (2D -> triangle o2)
							tab[1] = num1; tab[2] = num2; tab[3] = num3;  
							tab[4] = num4; tab[5] = num5; tab[6] = num6;
						}
						if (typeEle == 21) {
							fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
							fscanf(fid,"%d",&num7); fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
							num4--; num5--; num6--; num7--; num8--; num9--; num10--;
						
							tab = new int[11];
							tab[0] = 21; // Element Type (2D -> triangle o3)
							tab[1] = num1; tab[2] = num2; tab[3] = num3;  
							tab[4] = num4; tab[5] = num5; tab[6] = num6;
							tab[7] = num7; tab[8] = num8; tab[9] = num9; tab[10] = num10;
						}
					
						Mesh2D.push_back(tab);
						Mesh2Dtag.push_back(tag);
					}
				
				
					// Line 1D :
					if (typeEle == 1 or typeEle == 8 or typeEle == 26) {
						//fscanf(fid,"%s",s);fscanf(fid,"%s",s);
						int* tab;
						if (typeEle == 1) {
							fscanf(fid,"%d",&num1);
							fscanf(fid,"%d",&num2);
							num1--;
							num2--;
					
							tab = new int[3];
							tab[0] = 2; // Element Type (2D -> line)
							tab[1] = num1; tab[2] = num2;
						}
						// I do not know....
						if (typeEle == 8) { // Ligne ordre 2
							fscanf(fid,"%d",&num1);
							fscanf(fid,"%d",&num2);
							fscanf(fid,"%d",&num3);
							num1--;
							num2--;
							num3--;
						
							tab = new int[4];
							tab[0] = 8; // Element Type (1D -> line o2)
							tab[1] = num1; tab[2] = num2; tab[3] = num3;
					
						}
					
						if (typeEle == 26) { // Ligne ordre 3
							fscanf(fid,"%d",&num1);
							fscanf(fid,"%d",&num2);
							fscanf(fid,"%d",&num3);
							fscanf(fid,"%d",&num4);
							num1--;
							num2--;
							num3--;
							num4--;
						
							tab = new int[5];
							tab[0] = 8; // Element Type (1D -> line o3)
							tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
					
						}
					
						Mesh1D.push_back(tab);
						Mesh1Dtag.push_back(tag);
					
					}
				
					// Node 0D
					if (typeEle == 15) {
						fscanf(fid,"%d",&num1);
						num1--;
					
						Mesh0D.push_back(num1);
						Mesh0Dtag.push_back(tag);
					}
				}
				
			}
			
			/*
			for (int i=0; i<nElems; i++) {
				fscanf(fid,"%d",&numEle);
				fscanf(fid,"%d",&typeEle);
				fscanf(fid,"%d",&nbTag);
				fscanf(fid,"%d",&tag);
				
				// Tag reading: (usefull for domains informations)
				for (int j=1; j<nbTag; j++) {
					fscanf(fid,"%s",s);
				}
				
				// Tetrahedron 3D
				if(typeEle == 4 or typeEle == 11 or typeEle == 29) {
			
					fscanf(fid,"%d",&num1);
					fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);
					fscanf(fid,"%d",&num4);
					num1--;
					num2--;
					num3--;
					num4--;
					
					// Check if null volume tetrahedron :
					p1 = *(Nodes[num1]);
					p2 = *(Nodes[num2]);
					p3 = *(Nodes[num3]);
					p4 = *(Nodes[num4]);
		
					T2 = p2-p1;
					T3 = p3-p1;	
					T4 = p4-p1;
		
					// Choice s.t. p1 <-> (0,0,0) && p2 <-> (1,0,0) & p3 <-> (0,1,0) & p4 <-> (0,0,1) 
					A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
					A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
					A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
		
		
					double volume = fabs( ((A(0,0)*A(1,1)-A(0,1)*A(1,0))*A(2,2)-(A(0,0)*A(2,1)-A(0,1)*A(2,0))*A(1,2)+(A(1,0)*A(2,1)-A(1,1)*A(2,0))*A(0,2)));
					if (volume > 1e-10) {
						int* tab ;
						if (typeEle == 4) {
							tab = new int[5];
							tab[0] = 4; // Element Type (3D -> tetrahedron)
							tab[1] = num1; tab[2] = num2;
							tab[3] = num3; tab[4] = num4;
						}
						if (typeEle == 11) {
							fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
							fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
							num5--; num6--; num7--; num8--; num9--; num10--;
							
							tab = new int[11];
							tab[0] = 11; // Element Type (3D -> tetrahedron o2)
							tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
							tab[5] = num5; tab[6] = num6; tab[7] = num7; 
							tab[8] = num8; tab[9] = num9; tab[10] = num10;
						}
						if (typeEle == 29) {
							fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
							fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
							fscanf(fid,"%d",&num11); fscanf(fid,"%d",&num12); fscanf(fid,"%d",&num13);
							fscanf(fid,"%d",&num14); fscanf(fid,"%d",&num15); fscanf(fid,"%d",&num16);
							fscanf(fid,"%d",&num17); fscanf(fid,"%d",&num18); fscanf(fid,"%d",&num19); fscanf(fid,"%d",&num20);
							num5--; num6--; num7--; num8--; num9--; num10--;
							num11--; num12--; num13--; num14--; num15--; num16--;
							num17--; num18--; num19--; num20--;
							
							tab = new int[21];
							tab[0] = 29; // Element Type (3D -> tetrahedron o3)
							tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
							tab[5] = num5; tab[6] = num6; tab[7] = num7; tab[8] = num8; 
							tab[9] = num9; tab[10] = num10; tab[11] = num11; tab[12] = num12;
							tab[13] = num13; tab[14] = num14; tab[15] = num15; tab[16] = num16;
							tab[17] = num17; tab[18] = num18; tab[19] = num19; tab[20] = num20;
						}
						
						
						
						Mesh3D.push_back(tab);
						Mesh3Dtag.push_back(tag);
					}
					else{
						int PAUSE;
						cout<<"Null tetrahedron volume !!!"<<endl;
						cout<<i<<"	"<<typeEle<<endl;
						cin>>PAUSE;
					}	
				}
				
				
				// square 2D
				if(typeEle == 3) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
					fscanf(fid,"%d",&num1);fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);fscanf(fid,"%d",&num4);
					
					num1--;
					num2--;
					num3--;
					num4--;
					
					int* tab = new int[5];
					tab[0] = 3; // Element Type (2D -> square o1)
					tab[1] = num1; tab[2] = num2;
					tab[3] = num3;  tab[4] = num4;
					
					Mesh2D.push_back(tab);
					Mesh2Dtag.push_back(tag);
				}
			
				// Triangle 2D
				if (typeEle == 2 or typeEle == 9 or typeEle == 21) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
					fscanf(fid,"%d",&num1);
					fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);
					num1--;
					num2--;
					num3--;
					
					int* tab;
					if (typeEle == 2) {
						 tab = new int[4];
						 tab[0] = 3; // Element Type (2D -> triangle)
						 tab[1] = num1; tab[2] = num2; tab[3] = num3; 
					}
					if (typeEle == 9) {
						fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
						num4--; num5--; num6--;
						
						tab = new int[7];
						tab[0] = 6; // Element Type (2D -> triangle o2)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;  
						tab[4] = num4; tab[5] = num5; tab[6] = num6;
					}
					if (typeEle == 21) {
						fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
						fscanf(fid,"%d",&num7); fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
						num4--; num5--; num6--; num7--; num8--; num9--; num10--;
						
						tab = new int[11];
						tab[0] = 21; // Element Type (2D -> triangle o3)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;  
						tab[4] = num4; tab[5] = num5; tab[6] = num6;
						tab[7] = num7; tab[8] = num8; tab[9] = num9; tab[10] = num10;
					}
					
					Mesh2D.push_back(tab);
					Mesh2Dtag.push_back(tag);
				}
				
				
				// Line 1D :
				if (typeEle == 1 or typeEle == 8 or typeEle == 26) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					int* tab;
					if (typeEle == 1) {
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						num1--;
						num2--;
					
						tab = new int[3];
						tab[0] = 2; // Element Type (2D -> line)
						tab[1] = num1; tab[2] = num2;
					}
					// I do not know....
					if (typeEle == 8) { // Ligne ordre 2
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						num1--;
						num2--;
						num3--;
						
						tab = new int[4];
						tab[0] = 8; // Element Type (1D -> line o2)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;
					
					}
					
					if (typeEle == 26) { // Ligne ordre 2
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						fscanf(fid,"%d",&num4);
						num1--;
						num2--;
						num3--;
						num4--;
						
						tab = new int[5];
						tab[0] = 8; // Element Type (1D -> line o3)
						tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
					
					}
					
					Mesh1D.push_back(tab);
					Mesh1Dtag.push_back(tag);
					
				}
				
				// Node 0D
				if (typeEle == 15) {
					fscanf(fid,"%d",&num1);
					num1--;
					
					Mesh0D.push_back(num1);
					Mesh0Dtag.push_back(tag);
				}
			}
			*/
			fclose(fid);
		}
		
		
		
		
		void readGMSHold(char* name){
			FILE* fid;
			fid = fopen (name,"r");
			
			cout<<"Starting reading the mesh: "<<endl;
			MatriceL A(3);
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); Vecteur p4(3);
			Vecteur T1(3); Vecteur T2(3); Vecteur T3(3); Vecteur T4(3);
			// Nodes:
			char s[1000];
			while (!(s[0] == '$' && s[1] =='N' && s[2] == 'o' && s[3] == 'd' && s[4] == 'e' && s[5] == 's')) {
				fscanf(fid,"%s",s);
			}
			int nNodes =0;
			fscanf(fid,"%d",&nNodes);
			cout<<"Number of nodes: "<<nNodes<<endl;
			
			for (int i=0; i<nNodes; i++) {
				Vecteur* p = new Vecteur(3);
				double x,y,z;
				int num;
				fscanf(fid,"%d",&num);
				fscanf(fid,"%lf",&x);
				fscanf(fid,"%lf",&y);
				fscanf(fid,"%lf",&z);
				p->set(0,x);
				p->set(1,y);
				p->set(2,z);
				Nodes.push_back(p);
			}
			
			// Elements reading:
			fscanf(fid,"%s",s);
			while (!(s[0] == '$' && s[1] =='E' && s[2] == 'l' && s[3] == 'e' && s[4] == 'm' && s[5] == 'e')) {
				fscanf(fid,"%s",s);
			}
			int nElems = 0;
			fscanf(fid,"%d",&nElems);
			cout<<"Number of Elements: "<<nElems<<endl;
			
			int numEle, typeEle, nbTag;
			int num1, num2, num3, num4;
			int num5, num6, num7, num8;
			int num9, num10, num11, num12;
			int num13, num14, num15, num16;
			int num17, num18, num19, num20;
			int tag;
			
			for (int i=0; i<nElems; i++) {
				fscanf(fid,"%d",&numEle);
				fscanf(fid,"%d",&typeEle);
				fscanf(fid,"%d",&nbTag);
				fscanf(fid,"%d",&tag);
				
				// Tag reading: (usefull for domains informations)
				for (int j=1; j<nbTag; j++) {
					fscanf(fid,"%s",s);
				}
				
				// Tetrahedron 3D
				if(typeEle == 4 or typeEle == 11 or typeEle == 29) {
			
					fscanf(fid,"%d",&num1);
					fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);
					fscanf(fid,"%d",&num4);
					num1--;
					num2--;
					num3--;
					num4--;
					
					// Check if null volume tetrahedron :
					p1 = *(Nodes[num1]);
					p2 = *(Nodes[num2]);
					p3 = *(Nodes[num3]);
					p4 = *(Nodes[num4]);
		
					T2 = p2-p1;
					T3 = p3-p1;	
					T4 = p4-p1;
		
					// Choice s.t. p1 <-> (0,0,0) && p2 <-> (1,0,0) & p3 <-> (0,1,0) & p4 <-> (0,0,1) 
					A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
					A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
					A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
		
		
					double volume = fabs( ((A(0,0)*A(1,1)-A(0,1)*A(1,0))*A(2,2)-(A(0,0)*A(2,1)-A(0,1)*A(2,0))*A(1,2)+(A(1,0)*A(2,1)-A(1,1)*A(2,0))*A(0,2)));
					if (volume > 1e-10) {
						int* tab ;
						if (typeEle == 4) {
							tab = new int[5];
							tab[0] = 4; // Element Type (3D -> tetrahedron)
							tab[1] = num1; tab[2] = num2;
							tab[3] = num3; tab[4] = num4;
						}
						if (typeEle == 11) {
							fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
							fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
							num5--; num6--; num7--; num8--; num9--; num10--;
							
							tab = new int[11];
							tab[0] = 11; // Element Type (3D -> tetrahedron o2)
							tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
							tab[5] = num5; tab[6] = num6; tab[7] = num7; 
							tab[8] = num8; tab[9] = num9; tab[10] = num10;
						}
						if (typeEle == 29) {
							fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6); fscanf(fid,"%d",&num7);
							fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
							fscanf(fid,"%d",&num11); fscanf(fid,"%d",&num12); fscanf(fid,"%d",&num13);
							fscanf(fid,"%d",&num14); fscanf(fid,"%d",&num15); fscanf(fid,"%d",&num16);
							fscanf(fid,"%d",&num17); fscanf(fid,"%d",&num18); fscanf(fid,"%d",&num19); fscanf(fid,"%d",&num20);
							num5--; num6--; num7--; num8--; num9--; num10--;
							num11--; num12--; num13--; num14--; num15--; num16--;
							num17--; num18--; num19--; num20--;
							
							tab = new int[21];
							tab[0] = 29; // Element Type (3D -> tetrahedron o3)
							tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
							tab[5] = num5; tab[6] = num6; tab[7] = num7; tab[8] = num8; 
							tab[9] = num9; tab[10] = num10; tab[11] = num11; tab[12] = num12;
							tab[13] = num13; tab[14] = num14; tab[15] = num15; tab[16] = num16;
							tab[17] = num17; tab[18] = num18; tab[19] = num19; tab[20] = num20;
						}
						
						
						
						Mesh3D.push_back(tab);
						Mesh3Dtag.push_back(tag);
					}
					else{
						int PAUSE;
						cout<<"Null tetrahedron volume !!!"<<endl;
						cout<<i<<"	"<<typeEle<<endl;
						cin>>PAUSE;
					}	
				}
				
				
				// square 2D
				if(typeEle == 3) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
					fscanf(fid,"%d",&num1);fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);fscanf(fid,"%d",&num4);
					
					num1--;
					num2--;
					num3--;
					num4--;
					
					int* tab = new int[5];
					tab[0] = 3; // Element Type (2D -> square o1)
					tab[1] = num1; tab[2] = num2;
					tab[3] = num3;  tab[4] = num4;
					
					Mesh2D.push_back(tab);
					Mesh2Dtag.push_back(tag);
				}
			
				// Triangle 2D
				if (typeEle == 2 or typeEle == 9 or typeEle == 21) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					
					fscanf(fid,"%d",&num1);
					fscanf(fid,"%d",&num2);
					fscanf(fid,"%d",&num3);
					num1--;
					num2--;
					num3--;
					
					int* tab;
					if (typeEle == 2) {
						 tab = new int[4];
						 tab[0] = 3; // Element Type (2D -> triangle)
						 tab[1] = num1; tab[2] = num2; tab[3] = num3; 
					}
					if (typeEle == 9) {
						fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
						num4--; num5--; num6--;
						
						tab = new int[7];
						tab[0] = 6; // Element Type (2D -> triangle o2)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;  
						tab[4] = num4; tab[5] = num5; tab[6] = num6;
					}
					if (typeEle == 21) {
						fscanf(fid,"%d",&num4); fscanf(fid,"%d",&num5); fscanf(fid,"%d",&num6);
						fscanf(fid,"%d",&num7); fscanf(fid,"%d",&num8); fscanf(fid,"%d",&num9); fscanf(fid,"%d",&num10);
						num4--; num5--; num6--; num7--; num8--; num9--; num10--;
						
						tab = new int[11];
						tab[0] = 21; // Element Type (2D -> triangle o3)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;  
						tab[4] = num4; tab[5] = num5; tab[6] = num6;
						tab[7] = num7; tab[8] = num8; tab[9] = num9; tab[10] = num10;
					}
					
					Mesh2D.push_back(tab);
					Mesh2Dtag.push_back(tag);
				}
				
				
				// Line 1D :
				if (typeEle == 1 or typeEle == 8 or typeEle == 26) {
					//fscanf(fid,"%s",s);fscanf(fid,"%s",s);
					int* tab;
					if (typeEle == 1) {
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						num1--;
						num2--;
					
						tab = new int[3];
						tab[0] = 2; // Element Type (2D -> line)
						tab[1] = num1; tab[2] = num2;
					}
					// I do not know....
					if (typeEle == 8) { // Ligne ordre 2
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						num1--;
						num2--;
						num3--;
						
						tab = new int[4];
						tab[0] = 8; // Element Type (1D -> line o2)
						tab[1] = num1; tab[2] = num2; tab[3] = num3;
					
					}
					
					if (typeEle == 26) { // Ligne ordre 2
						fscanf(fid,"%d",&num1);
						fscanf(fid,"%d",&num2);
						fscanf(fid,"%d",&num3);
						fscanf(fid,"%d",&num4);
						num1--;
						num2--;
						num3--;
						num4--;
						
						tab = new int[5];
						tab[0] = 8; // Element Type (1D -> line o3)
						tab[1] = num1; tab[2] = num2; tab[3] = num3; tab[4] = num4;
					
					}
					
					Mesh1D.push_back(tab);
					Mesh1Dtag.push_back(tag);
					
				}
				
				// Node 0D
				if (typeEle == 15) {
					fscanf(fid,"%d",&num1);
					num1--;
					
					Mesh0D.push_back(num1);
					Mesh0Dtag.push_back(tag);
				}
			}
			fclose(fid);
		}
		
		
		vector<double> readTXT(char* name){
			vector<double> res;
			FILE* fid;
			fid = fopen (name,"r");
			int nbTerm;
			fscanf(fid,"%d",&nbTerm);
			for(int i=0;i<nbTerm;i++) {
				double val;
				fscanf(fid,"%lf",&val);
				res.push_back(val);
			}
			fclose(fid);
			return res;
		}
		
		
		void writeVTK3D(char* name, Vecteur U, vector<int*> meshEqP1,char* type){
		
			FILE* fid;
			fid = fopen(name,"w");

			fprintf(fid, "# vtk DataFile Version 3.6\n");
			fprintf(fid, "PolyDATA\n");
			fprintf(fid, "ASCII\n");
			fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

			int nbFct = getNNodes();
			
			fprintf(fid, "POINTS	%d	float\n",nbFct);

			for(int i=0;i<nbFct;i++) {
				fprintf(fid, "%lf	%lf	%lf\n",( (*(Nodes[i]))(0)), (*(Nodes[i]))(1), (*(Nodes[i]))(2));
			}
			
			/**/
			fprintf(fid, "CELLS	%d	%d\n",int(meshEqP1.size()), int(meshEqP1.size())*5);
			
			//for (int st=0; st<listeTertraDsMaille.size(); st++) {
			//	cout<<st<<"	"<<listeTertraDsMaille[st][0]<<" "<<listeTertraDsMaille[st][1]<<" "<<listeTertraDsMaille[st][2]<<" "<<listeTertraDsMaille[st][3]<<endl;
			//}
			//cin>>PAUSE;
			
			
			for(int m=0;m<meshEqP1.size();m++)
				fprintf(fid, "4	%d	%d	%d	%d\n",meshEqP1[m][1],meshEqP1[m][2],meshEqP1[m][3],meshEqP1[m][4]);
			
			fprintf(fid, "CELL_TYPES	%d\n",int(meshEqP1.size()));
			for(int m=0;m<meshEqP1.size();m++) {
				fprintf(fid, "10\n");
			}
			

			if (type == "scalar") {
				fprintf(fid, "POINT_DATA	%d\n",nbFct);
				fprintf(fid,"SCALARS scalars float	1\n");
				fprintf(fid,"LOOKUP_TABLE default\n");

				for(int i=0;i<nbFct;i++) 
					fprintf(fid, "%.15lf\n",U(i) );
			}
			if(type == "vectorial") {
				fprintf(fid, "POINT_DATA	%d\n",nbFct);
				//fprintf(fid,"SCALARS scalars float	1\n");
				fprintf(fid,"VECTORS deformation float\n");

				for(int i=0;i<nbFct;i++) 
					fprintf(fid, "%lf	%lf	%lf\n",U(3*i),U(3*i+1), U(3*i+2) );
			}
			if (type == "cell") {
				fprintf(fid, "CELL_DATA	%d\n",meshEqP1.size());
				fprintf(fid,"SCALARS cell_scalars float	1\n");
				fprintf(fid,"LOOKUP_TABLE default\n");
				for (int i=0;i<meshEqP1.size();i++){
					if(lienM2DPkP1.size() == 0)
						fprintf(fid, "%.15lf\n",U(i) );
					else
						fprintf(fid, "%.15lf\n",U(lienM2DPkP1[i]) );
					//fprintf(fid, "%.15lf\n",U(lienM2DPkP1[i]) );
				}
			}

			fclose(fid);
		}
		
		void writeTXT(char* name, vector<double> sig){
		
			FILE* fid;
			fid = fopen(name,"w");

			fprintf(fid, "%d\n",sig.size());
			for(int i=0;i<sig.size();i++) {
				fprintf(fid, "%.15lf\n",sig[i]);
			}
			fclose(fid);
		}
		
		
		void writeVTK2D(char* name, Vecteur U, vector<int*> meshEqP1,char* type) {
		
			FILE* fid;
			fid = fopen(name,"w");

			fprintf(fid, "# vtk DataFile Version 3.6\n");
			fprintf(fid, "PolyDATA\n");
			fprintf(fid, "ASCII\n");
			fprintf(fid, "DATASET UNSTRUCTURED_GRID\n");

			int nbFct = getNNodes();
			
			fprintf(fid, "POINTS	%d	float\n",nbFct);

			for(int i=0;i<nbFct;i++) {
				if (type == "scalar3D")
					fprintf(fid, "%lf	%lf	%lf\n",( (*(Nodes[i]))(0)), (*(Nodes[i]))(1), U(i));
				else
					fprintf(fid, "%lf	%lf	%lf\n",( (*(Nodes[i]))(0)), (*(Nodes[i]))(1), (*(Nodes[i]))(2));
			}

			fprintf(fid, "CELLS	%d	%d\n",int(meshEqP1.size()), int(meshEqP1.size())*4);
			/*
			for (int st=0; st<listeTertraDsMaille.size(); st++) {
				cout<<st<<"	"<<listeTertraDsMaille[st][0]<<" "<<listeTertraDsMaille[st][1]<<" "<<listeTertraDsMaille[st][2]<<" "<<listeTertraDsMaille[st][3]<<endl;
			}
			cin>>PAUSE;
			*/
			
			for(int m=0;m<meshEqP1.size();m++)
				fprintf(fid, "3	%d	%d	%d\n",meshEqP1[m][1],meshEqP1[m][2],meshEqP1[m][3]);
			
			fprintf(fid, "CELL_TYPES	%d\n",int(meshEqP1.size()));
			for(int m=0;m<meshEqP1.size();m++) {
				fprintf(fid, "5\n");
			}
			
			if (type == "scalar" or type == "scalar3D") {
				fprintf(fid, "POINT_DATA	%d\n",nbFct);
				fprintf(fid,"SCALARS scalars float	1\n");
				fprintf(fid,"LOOKUP_TABLE default\n");

				for(int i=0;i<nbFct;i++) 
					fprintf(fid, "%.15lf\n",U(i) );
			}
			if(type == "vectorial") {
				fprintf(fid, "POINT_DATA	%d\n",nbFct);
				fprintf(fid,"VECTORS deformation float\n");

				for(int i=0;i<nbFct;i++) 
					fprintf(fid, "%lf	%lf	0.0\n",U(2*i),U(2*i+1) );
			}
			if (type == "cell") {
				fprintf(fid, "CELL_DATA	%d\n",meshEqP1.size());
				fprintf(fid,"SCALARS cell_scalars float	1\n");
				fprintf(fid,"LOOKUP_TABLE default\n");
				for (int i=0;i<meshEqP1.size();i++){
					if(lienM2DPkP1.size() == 0)
						fprintf(fid, "%.15lf\n",U(i) );
					else
						fprintf(fid, "%.15lf\n",U(lienM2DPkP1[i]) );
					//fprintf(fid, "%.15lf\n",U(lienM2DPkP1[i]) );
				}
			}
			fclose(fid);
		}
		
		void convertP1ToPk2D(int nbPtPMaille, double** treilli) {
			MatriceL lienAdj(Nodes.size());
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3);
			Vecteur T2(3); Vecteur T3(3);
			Vecteur pRef(3); Vecteur ImpRef(3);
			
			int num1, num2, num3;
			MatriceL A(3);
			vector<int*>* newMesh = new vector<int*>;
			vector<int>* newMeshtag = new vector<int>;
			
			// Triangles:
			for (int m=0; m<Mesh2D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh2D.size()<<endl;
				
				num1 = Mesh2D[m][1]; num2 = Mesh2D[m][2]; num3 = Mesh2D[m][3];
				int tag = Mesh2Dtag[m];
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1;
				
				A.set(0,0,T2(0)); A.set(1,0,T2(1));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,2,1.0);
				
				
				int* tab = new int[nbPtPMaille+1];
				tab[0] = Mesh2D[m][0]; // Mesh Type
				tab[1] = Mesh2D[m][1]; // Coins triangle
				tab[2] = Mesh2D[m][2]; // Coins triangle
				tab[3] = Mesh2D[m][3]; // Coins triangle
				
				int cpt = 4;
				for(int p=0; p<nbPtPMaille; p++) {
					pRef.set(0,treilli[p][0]); pRef.set(1,treilli[p][1]); pRef.set(2,0.0);
					ImpRef = A*pRef + p1;
					if ((ImpRef-p1).norme()>1e-10 && (ImpRef-p2).norme()>1e-10 && (ImpRef-p3).norme()>1e-10) { 
						if (fabs(pRef(0)*pRef(1))>1e-10 && fabs(pRef(0)+pRef(1))<(1.- 1e-10) ) { 
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
						} // Interior Point 
					
						if (lienAdj(num1,num2) == 0 && fabs(pRef(1)) < 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
						}	
						if (lienAdj(num1,num2) != 0 && fabs(pRef(1)) < 1e-10) {
							int mailleAdj = lienAdj(num1,num2)-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
						}
					
						if (lienAdj(num1,num3) == 0 && fabs(pRef(0)) <1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
						}	
						if (lienAdj(num1,num3) != 0 && fabs(pRef(0)) <1e-10) {
							int mailleAdj = lienAdj(num1,num3)-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) {
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
						}
					
						if (lienAdj(num2,num3) == 0 && fabs(pRef(0)+pRef(1)-1.0)<1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
						}	
						if (lienAdj(num2,num3) != 0 && fabs(pRef(0)+pRef(1)-1.0)<1e-10) {
							int mailleAdj = lienAdj(num2,num3)-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) {
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
						}
						cpt++;
					}
				}
				(*newMesh).push_back(tab);
				(*newMeshtag).push_back(tag);
				lienAdj.set(num1,num2,(*newMesh).size()); lienAdj.set(num2,num1,(*newMesh).size());
				lienAdj.set(num1,num3,(*newMesh).size()); lienAdj.set(num3,num1,(*newMesh).size());
				lienAdj.set(num2,num3,(*newMesh).size()); lienAdj.set(num3,num2,(*newMesh).size());
			}
			
			Mesh2D = (*newMesh);
			Mesh2Dtag = (*newMeshtag);
		}
		
		
		
		vector<int*> convertPkToP12D(int nbPtPMaille, double** treilli, vector<int*> decoup) {
			
			int ordre = int( ( -1 + sqrt(1+8*nbPtPMaille))*0.5 -1 );
			//cout<<"Ordre : "<<ordre<<endl;
			
			Vecteur pRef(3); Vecteur ImpRef(3);
			Vecteur p1Ref(3); Vecteur p2Ref(3); Vecteur p3Ref(3); Vecteur p4Ref(3);
			Vecteur Imp1Ref(3); Vecteur Imp2Ref(3); Vecteur Imp3Ref(3); Vecteur Imp4Ref(3);
			
			Vecteur p(3);
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3);
			Vecteur T2(3); Vecteur T3(3);
			//Vecteur pRef(3); Vecteur ImpRef(3);
			
			int num1, num2, num3;
			int ind1, ind2, ind3, ind4;
			MatriceL A(3);
			vector<int*> meshP1;
			
			
			int* correspond = new int[nbPtPMaille];
			for(int i=0;i<nbPtPMaille;i++)
				correspond[i] = -1;
			
			for (int ii=0; ii< nbPtPMaille; ii++){
				num1 = Mesh2D[0][1]; num2 = Mesh2D[0][2]; num3 = Mesh2D[0][3];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1;
				A.set(0,0,T2(0)); A.set(1,0,T2(1));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,2,1.0);
				
				pRef.set(0,treilli[ii][0]); pRef.set(1,treilli[ii][1]); pRef.set(2,treilli[ii][2]);
				ImpRef = p1 + A * pRef;
				for(int jj=0; jj< nbPtPMaille; jj++) {
					int num = Mesh2D[0][jj+1];
					p = *(Nodes[num]);
					
					if ((p - ImpRef).norme() < 1e-10) {
						correspond[ii] = jj;
						break;
					}
				}
			}
			
			
			// Triangles:
			lienM2DPkP1.clear();
			for (int m=0; m<Mesh2D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh2D.size()<<endl;
				
				
				num1 = Mesh2D[m][1]; num2 = Mesh2D[m][2]; num3 = Mesh2D[m][3];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1;
				A.set(0,0,T2(0)); A.set(1,0,T2(1));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,2,1.0);
				
				for(int ii=0;ii<decoup.size();ii++) {
					int ind1 = correspond[decoup[ii][0]]+1;
					int ind2 = correspond[decoup[ii][1]]+1;
					int ind3 = correspond[decoup[ii][2]]+1;
					
					int* tab = new int[4];
					tab[0] = 2; //Mesh2D[m][0];
					tab[1] = Mesh2D[m][ind1]; tab[2] = Mesh2D[m][ind2]; tab[3] = Mesh2D[m][ind3];
					meshP1.push_back(tab);
					
					lienM2DPkP1.push_back(m);
				}
				
			}
			
			// Triangles:
			/*
			for (int m=0; m<Mesh2D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh2D.size()<<endl;
				
				num1 = Mesh2D[m][1]; num2 = Mesh2D[m][2]; num3 = Mesh2D[m][3];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1;
				A.set(0,0,T2(0)); A.set(1,0,T2(1));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,2,1.0);
				
				for(int i=0; i<ordre; i++) {
					for (int j=0; j<ordre-i; j++) {
						p1Ref.set(0,j*1.0/ordre); p1Ref.set(1,i*1.0/ordre);
						p2Ref.set(0,(j+1)*1.0/ordre); p2Ref.set(1,i*1.0/ordre);
						p3Ref.set(0,j*1.0/ordre); p3Ref.set(1,(i+1)*1.0/ordre);
						p4Ref.set(0,(j+1)*1.0/ordre); p4Ref.set(1,(i+1)*1.0/ordre);
						
						Imp1Ref = A*p1Ref + p1; Imp2Ref = A*p2Ref + p1; 
						Imp3Ref = A*p3Ref + p1; Imp4Ref = A*p4Ref + p1;
						
 						for (int ii=0; ii< nbPtPMaille; ii++){
 							p = *(Nodes[Mesh2D[m][ii+1]]);
							if ((p - Imp1Ref).norme() < 1e-10)
								ind1 = Mesh2D[m][ii+1];
							if ((p - Imp2Ref).norme() < 1e-10)
								ind2 = Mesh2D[m][ii+1];
							if ((p - Imp3Ref).norme() < 1e-10)
								ind3 = Mesh2D[m][ii+1];
							if ((p - Imp4Ref).norme() < 1e-10)
								ind4 = Mesh2D[m][ii+1];
 						}
						int* tab = new int[4];
						tab[0] = 2; //Mesh2D[m][0];
						tab[1] = ind1; tab[2] = ind2; tab[3] = ind3;
						meshP1.push_back(tab);
						
						if (j < ordre-i-1) { //  Correction 21 12 19
							int* tabbis = new int[4];
							tabbis[0] = 2; //Mesh2D[m][0];
							tabbis[1] = ind2; tabbis[2] = ind4; tabbis[3] = ind3;
							meshP1.push_back(tabbis);
						}
					}
				}
			}
			*/
			return meshP1;
		}
		
		
		vector<int*> convertPkToP12DwTag(int nbPtPMaille, double** treilli, vector<int>* m2DTag) {
			
			int ordre = int( ( -1 + sqrt(1+8*nbPtPMaille))*0.5 -1 );
			//cout<<"Ordre : "<<ordre<<endl;
			
			m2DTag->clear();
			
			Vecteur p1Ref(3); Vecteur p2Ref(3); Vecteur p3Ref(3); Vecteur p4Ref(3);
			Vecteur Imp1Ref(3); Vecteur Imp2Ref(3); Vecteur Imp3Ref(3); Vecteur Imp4Ref(3);
			
			Vecteur p(3);
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3);
			Vecteur T2(3); Vecteur T3(3); Vecteur T4(3); // Artificielle
			Vecteur Tortho(3);
			Vecteur pRef(3); Vecteur ImpRef(3);
			
			int num1, num2, num3;
			int ind1, ind2, ind3, ind4;
			MatriceL A(3);
			vector<int*> meshP1;
			
			// Triangles:
			for (int m=0; m<Mesh2D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh2D.size()<<endl;
				
				num1 = Mesh2D[m][1]; num2 = Mesh2D[m][2]; num3 = Mesh2D[m][3];
				int tag = Mesh2Dtag[m];
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1;
				
				T4.set(0,T2(1)*T3(2) - T2(2)*T3(1));
				T4.set(1,-T2(0)*T3(2) + T2(2)*T3(0));
				T4.set(2,T2(0)*T3(1) - T2(1)*T3(0)); 
				
				A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
				A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2)); 
				
				for(int i=0; i<ordre; i++) {
					for (int j=0; j<ordre-i; j++) {
						p1Ref.set(0,j*1.0/ordre); p1Ref.set(1,i*1.0/ordre);
						p2Ref.set(0,(j+1)*1.0/ordre); p2Ref.set(1,i*1.0/ordre);
						p3Ref.set(0,j*1.0/ordre); p3Ref.set(1,(i+1)*1.0/ordre);
						p4Ref.set(0,(j+1)*1.0/ordre); p4Ref.set(1,(i+1)*1.0/ordre);
						
						Imp1Ref = A*p1Ref + p1; Imp2Ref = A*p2Ref + p1; 
						Imp3Ref = A*p3Ref + p1; Imp4Ref = A*p4Ref + p1;
						
 						for (int ii=0; ii< nbPtPMaille; ii++){
 							p = *(Nodes[Mesh2D[m][ii+1]]);
							if ((p - Imp1Ref).norme() < 1e-10)
								ind1 = Mesh2D[m][ii+1];
							if ((p - Imp2Ref).norme() < 1e-10)
								ind2 = Mesh2D[m][ii+1];
							if ((p - Imp3Ref).norme() < 1e-10)
								ind3 = Mesh2D[m][ii+1];
							if ((p - Imp4Ref).norme() < 1e-10)
								ind4 = Mesh2D[m][ii+1];
 						}
						int* tab = new int[4];
						tab[0] = 2; //Mesh2D[m][0];
						tab[1] = ind1; tab[2] = ind2; tab[3] = ind3;
						meshP1.push_back(tab);
						m2DTag->push_back(tag);
						
						if (j < ordre-i-1) { //  Correction 21 12 19
							int* tabbis = new int[4];
							tabbis[0] = 2; //Mesh2D[m][0];
							tabbis[1] = ind2; tabbis[2] = ind4; tabbis[3] = ind3;
							meshP1.push_back(tabbis);
							m2DTag->push_back(tag);
						}
					}
				}
			}
			
			return meshP1;
		}
		
		
		
		void convertP1ToPk3D(int nbPtPMaille, double** treilli) {
			MatriceL connec1D(Nodes.size());
			vector< int* > connec2D;
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); Vecteur p4(3);
			Vecteur T2(3); Vecteur T3(3); Vecteur T4(3);
			Vecteur pRef(3); Vecteur ImpRef(3);
			
			int num1, num2, num3, num4;
			MatriceL A(3);
			vector<int*>* newMesh = new vector<int*>;
			vector<int>* newMeshtag = new vector<int>;
			
			// Triangles:
			for (int m=0; m<Mesh3D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh3D.size()<<endl;
				
				num1 = Mesh3D[m][1]; num2 = Mesh3D[m][2]; num3 = Mesh3D[m][3]; num4 = Mesh3D[m][4];
				int tag = Mesh3Dtag[m];
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				p4 = *(Nodes[num4]);
				
				T2 = p2 - p1; T3 = p3 - p1; T4 = p4 - p1;
				
				A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2));
				A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
				
				int* tab = new int[nbPtPMaille+1];
				tab[0] = Mesh3D[m][0]; // Mesh Type
				tab[1] = Mesh3D[m][1]; // Coins tetra
				tab[2] = Mesh3D[m][2]; // Coins tetra
				tab[3] = Mesh3D[m][3]; // Coins tetra
				tab[4] = Mesh3D[m][4]; // Coins tetra
				
				int cpt = 5;
				int i123 = -1;
				int i124 = -1;
				int i134 = -1;
				int i234 = -1;
				
				// num1 <-> (0,0,0)
				// num2 <-> (1,0,0)
				// num3 <-> (0,1,0)
				// num4 <-> (0,0,1)
				
				int num123 = -1;
				int num124 = -1;
				int num134 = -1;
				int num234 = -1;
				for(int i=0;i<connec2D.size();i++) {
					//cout<<"i : "<<i<<endl;
					int n1 = connec2D[i][0];
					int n2 = connec2D[i][1];
					int n3 = connec2D[i][2];
					int mm = connec2D[i][3];
					if( (n1 == num1 and n2 == num2 and n3 == num3) or (n1 == num1 and n3 == num2 and n2 == num3) or (n2 == num1 and n1 == num2 and n3 == num3) or (n2 == num1 and n3 == num2 and n1 == num3) or (n3 == num1 and n1 == num2 and n2 == num3) or (n3 == num1 and n2 == num2 and n1 == num3) ) {
						num123 = mm;
						i123 = i;
					}
					if( (n1 == num1 and n2 == num2 and n3 == num4) or (n1 == num1 and n3 == num2 and n2 == num4) or (n2 == num1 and n1 == num2 and n3 == num4) or (n2 == num1 and n3 == num2 and n1 == num4) or (n3 == num1 and n1 == num2 and n2 == num4) or (n3 == num1 and n2 == num2 and n1 == num4) ) {
						num124 = mm;
						i124 = i;
					}
					if( (n1 == num1 and n2 == num3 and n3 == num4) or (n1 == num1 and n3 == num3 and n2 == num4) or (n2 == num1 and n1 == num3 and n3 == num4) or (n2 == num1 and n3 == num3 and n1 == num4) or (n3 == num1 and n1 == num3 and n2 == num4) or (n3 == num1 and n2 == num3 and n1 == num4) ) {
						num134 = mm;
						i134 = i;
					}
					if( (n1 == num2 and n2 == num3 and n3 == num4) or (n1 == num2 and n3 == num3 and n2 == num4) or (n2 == num2 and n1 == num3 and n3 == num4) or (n2 == num2 and n3 == num3 and n1 == num4) or (n3 == num2 and n1 == num3 and n2 == num4) or (n3 == num2 and n2 == num3 and n1 == num4) ) {
						num234 = mm;
						i234 = i;
					}
				}
				int num12 = connec1D(num1,num2);
				int num13 = connec1D(num1,num3);
				int num23 = connec1D(num2,num3);
				int num14 = connec1D(num1,num4);
				int num24 = connec1D(num2,num4);
				int num34 = connec1D(num3,num4);
				
				//cout<<"onk : "<<i123<<"	"<<i124<<"	"<<i134<<"	"<<i234<<endl;
				//cout<<"onk : "<<num123<<"	"<<num124<<"	"<<num134<<"	"<<num234<<endl;
				for(int p=0; p<nbPtPMaille; p++) {
					//pRef.set(0,treilli[nbPtPMaille-1-p][0]); pRef.set(1,treilli[nbPtPMaille-1-p][1]); pRef.set(2,treilli[nbPtPMaille-1-p][2]);
					//cout<<treilli[nbPtPMaille-1-p][0]<<" "<<treilli[nbPtPMaille-1-p][1]<<" "<<treilli[nbPtPMaille-1-p][2]<<endl;
					pRef.set(0,treilli[p][0]); pRef.set(1,treilli[p][1]); pRef.set(2,treilli[p][2]);
					//cout<<treilli[p][0]<<" "<<treilli[p][1]<<" "<<treilli[p][2]<<endl;
					ImpRef = A*pRef + p1;
					// Verifie qu'on est pas sur un coin du tetra :
					if ((ImpRef-p1).norme()>1e-10 && (ImpRef-p2).norme()>1e-10 && (ImpRef-p3).norme()>1e-10 && (ImpRef-p4).norme()>1e-10) { 
						
						bool found = false;
						/**/
						if (fabs(pRef(0)*pRef(1)*pRef(2))>1e-10 && fabs(pRef(0)+pRef(1)+pRef(2) - 1.0) > 1e-10 ) { 
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
							
						} // Interior Point  
						//cout<<"Cpt : "<<cpt<<"	/	"<<nbPtPMaille<<endl;
						
						// Interieur face
						if (!found && num123 == -1 && fabs(pRef(2)) < 1e-10 && fabs(pRef(0)*pRef(1)*(pRef(0)+pRef(1)-1.0)) > 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
							
						}	
						if (!found && num123 != -1 && fabs(pRef(2)) < 1e-10 && fabs(pRef(0)*pRef(1)*(pRef(0)+pRef(1)-1.0)) > 1e-10) {
							int mailleAdj = num123;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
							
						}
						
						if (!found && num124 ==-1 && fabs(pRef(1)) < 1e-10 && fabs(pRef(0)*pRef(2)*(pRef(0)+pRef(2)-1.0)) > 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
							
						}	
						if (!found && num124 != -1 && fabs(pRef(1)) < 1e-10 && fabs(pRef(0)*pRef(2)*(pRef(0)+pRef(2)-1.0)) > 1e-10) {
							int mailleAdj = num124;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
							
						}
						
						if (!found && num134 == -1 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1)*pRef(2)*(pRef(1)+pRef(2)-1.0)) > 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
							
						}	
						if (!found && num134 != -1 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1)*pRef(2)*(pRef(1)+pRef(2)-1.0)) > 1e-10) {
							int mailleAdj = num134;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
							
						}
						
						if (!found && num234 == -1 && fabs(pRef(0)+pRef(1)+pRef(2)-1.0) < 1e-10 && fabs(pRef(2)*pRef(1)*pRef(0)) > 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
							
						}	
						if (!found && num234 != -1 && fabs(pRef(0)+pRef(1)+pRef(2)-1.0) < 1e-10  && fabs(pRef(2)*pRef(1)*pRef(0)) > 1e-10) {
							int mailleAdj = num234;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
							
						}
						/**/
						if(!found && num12 == 0 && fabs(pRef(1)) < 1e-10 && fabs(pRef(2)) < 1e-10 ) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num12 != 0 && fabs(pRef(1)) < 1e-10 &&  fabs(pRef(2)) < 1e-10) {
							int mailleAdj = num12-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(!found && num13 == 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(2)) < 1e-10 ) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num13 != 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(2)) < 1e-10) {
							int mailleAdj = num13-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(!found && num23 == 0 && fabs((pRef(0) + pRef(1) -1.0)) < 1e-10 && fabs(pRef(2)) < 1e-10 ) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num23 != 0 && fabs((pRef(0) + pRef(1) -1.0)) < 1e-10 && fabs(pRef(2)) < 1e-10) {
							int mailleAdj = num23-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(!found && num14 == 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1)) < 1e-10 ) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num14 != 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1)) < 1e-10) {
							int mailleAdj = num14-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(!found && num24 == 0 && fabs(pRef(1)) < 1e-10 && fabs(pRef(0) + pRef(2) - 1.0) < 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num24 != 0 && fabs(pRef(1)) < 1e-10 && fabs(pRef(0) + pRef(2) - 1.0) < 1e-10) {
							int mailleAdj = num24-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(!found && num34 == 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1) + pRef(2) - 1.0) < 1e-10) {
							Vecteur* newImpRef = new Vecteur(3);
							newImpRef->set(0,ImpRef(0)); newImpRef->set(1,ImpRef(1)); newImpRef->set(2,ImpRef(2));
							Nodes.push_back(newImpRef);
							tab[cpt] = Nodes.size()-1;
							found = true;
						}
						if(!found && num34 != 0 && fabs(pRef(0)) < 1e-10 && fabs(pRef(1) + pRef(2) - 1.0) < 1e-10) {
							int mailleAdj = num34-1;
							for	(int ele=0; ele<nbPtPMaille;ele++) {
								if ((*(Nodes[(*newMesh)[mailleAdj][ele+1]])-ImpRef).norme() < 1e-10) { 
									tab[cpt] = (*newMesh)[mailleAdj][ele+1];
									break;
								}
							}
							found = true;
						}
						
						if(found == false){
							cout<<"Probleme !!"<<endl;
							int PAUSE;
							cin>>PAUSE;
						}
						
						cpt++;
					}
				}
				//cout<<"AICI "<<endl;
				//cout<<num1<<" "<<num2<<" "<<num3<<" "<<num4<<endl;
				//for(int i=0;i<nbPtPMaille;i++)
				//	cout<<tab[i+1]<<" ";
				//cout<<endl;
				(*newMesh).push_back(tab);
				(*newMeshtag).push_back(tag);
				int nbTtoPop = 0;
				int nbT = connec2D.size();
				
				if(i123 != -1){
					int* tmp = connec2D[i123];
					connec2D[i123] = connec2D[nbT-1];
					connec2D[nbT-1] = tmp;
					if(nbT-1-nbTtoPop == i124)
						i124 = i123;
					if(nbT-1-nbTtoPop == i134)
						i134 = i123;
					if(nbT-1-nbTtoPop == i234)
						i234 = i123;
					nbTtoPop = nbTtoPop+1;
					
					
				}
				
				if(i124 != -1){
					int* tmp = connec2D[i124];
					connec2D[i124] = connec2D[nbT-1-nbTtoPop];
					connec2D[nbT-1-nbTtoPop] = tmp;
					if(nbT-1-nbTtoPop == i134)
						i134 = i124;
					if(nbT-1-nbTtoPop == i234)
						i234 = i124;
					nbTtoPop = nbTtoPop+1;
				}
				
				if(i134 != -1){
					int* tmp = connec2D[i134];
					connec2D[i134] = connec2D[nbT-1-nbTtoPop];
					connec2D[nbT-1-nbTtoPop] = tmp;
					if(nbT-1-nbTtoPop == i234)
						i234 = i134;
					nbTtoPop = nbTtoPop+1;
				}
				
				if(i234 != -1){
					int* tmp = connec2D[i234];
					connec2D[i234] = connec2D[nbT-1-nbTtoPop];
					connec2D[nbT-1-nbTtoPop] = tmp;
					nbTtoPop = nbTtoPop+1;
				}
				/**/
				
				for(int i=0;i<nbTtoPop;i++) {
					
					//nbT = connec2D.size();
					//cout<<"POP ! "<<nbT<<endl;
					//delete(connec[nbT-1]);
					connec2D.pop_back();
				}
				/**/
				
				if(i123 == -1) {
					int* tabC = new int[4]; 
					tabC[0] = num1;
					tabC[1] = num2;
					tabC[2] = num3;
					tabC[3] = (*newMesh).size()-1;
					connec2D.push_back(tabC);
				}
				if(i124 == -1) {
					int* tabC = new int[4]; 
					tabC[0] = num1;
					tabC[1] = num2;
					tabC[2] = num4;
					tabC[3] = (*newMesh).size()-1;
					connec2D.push_back(tabC);
				}
				if(i134 == -1) {
					int* tabC = new int[4]; 
					tabC[0] = num1;
					tabC[1] = num3;
					tabC[2] = num4;
					tabC[3] = (*newMesh).size()-1;
					connec2D.push_back(tabC);
				}
				if(i234 == -1) {
					int* tabC = new int[4]; 
					tabC[0] = num2;
					tabC[1] = num3;
					tabC[2] = num4;
					tabC[3] = (*newMesh).size()-1;
					connec2D.push_back(tabC);
				}
				if (num12 == 0) {
					connec1D.set(num1,num2,(*newMesh).size());
					connec1D.set(num2,num1,(*newMesh).size());
				}
				if(num13 == 0) {
					connec1D.set(num1,num3,(*newMesh).size());
					connec1D.set(num3,num1,(*newMesh).size());
				}
				if(num23 == 0) {
					connec1D.set(num2,num3,(*newMesh).size());
					connec1D.set(num3,num2,(*newMesh).size());
				}
				if(num14 == 0) {
					connec1D.set(num1,num4,(*newMesh).size());
					connec1D.set(num4,num1,(*newMesh).size());
				}
				if(num24 == 0) {
					connec1D.set(num2,num4,(*newMesh).size());	
					connec1D.set(num4,num2,(*newMesh).size());
				}
				if(num34 == 0) {
					connec1D.set(num3,num4,(*newMesh).size());
					connec1D.set(num4,num3,(*newMesh).size());
				}
				
				
				//cout<<"len connec2D : "<<connec2D.size()<<endl;
			}
			//cout<<"ON K FIN "<<endl;
			//for(int i=0;i<connec.size();i++)
			//	delete(connec[i]);
			
			Mesh3D = (*newMesh);
			Mesh3Dtag = (*newMeshtag);
		}
		
		vector<int*> convertPkToP13D(int nbPtPMaille, double** treilli, int ordre, vector<int*> decoup) {
			
			//cout<<"Convert Pk to P1 3D :"<<endl;
			//cout<<"Ordre : "<<ordre<<endl;
			
			Vecteur p1Ref(3); Vecteur p2Ref(3); Vecteur p3Ref(3); Vecteur p4Ref(3);
			Vecteur p5Ref(3); Vecteur p6Ref(3); Vecteur p7Ref(3); Vecteur p8Ref(3);
			Vecteur Imp1Ref(3); Vecteur Imp2Ref(3); Vecteur Imp3Ref(3); Vecteur Imp4Ref(3);
			Vecteur Imp5Ref(3); Vecteur Imp6Ref(3); Vecteur Imp7Ref(3); Vecteur Imp8Ref(3);
			
			Vecteur p(3);
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3); Vecteur p4(3);
			Vecteur T2(3); Vecteur T3(3); Vecteur T4(3);
			Vecteur pRef(3); Vecteur ImpRef(3);
			
			int num1, num2, num3, num4;
			int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
			MatriceL A(3);
			vector<int*> meshP1;
			
			
			
			int* correspond = new int[nbPtPMaille];
			for(int i=0;i<nbPtPMaille;i++)
				correspond[i] = -1;
			
			num1 = Mesh3D[0][1]; num2 = Mesh3D[0][2]; num3 = Mesh3D[0][3]; num4 = Mesh3D[0][4];
			
			p1 = *(Nodes[num1]);
			p2 = *(Nodes[num2]);
			p3 = *(Nodes[num3]);
			p4 = *(Nodes[num4]);
			
			T2 = p2 - p1; T3 = p3 - p1; T4 = p4-p1;
			A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
			A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2)); 
			A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
			
			for (int ii=0; ii< nbPtPMaille; ii++){
				pRef.set(0,treilli[ii][0]); pRef.set(1,treilli[ii][1]); pRef.set(2,treilli[ii][2]);
				ImpRef = p1 + A * pRef;
				for(int jj=0; jj< nbPtPMaille; jj++) {
					int num = Mesh3D[0][jj+1];
					p = *(Nodes[num]);
					
					if ((p - ImpRef).norme() < 1e-10) {
						correspond[ii] = jj;
						break;
					}
				}
			}
			
			
			// Triangles:
			for (int m=0; m<Mesh3D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh3D.size()<<endl;
				
				
				num1 = Mesh3D[m][1]; num2 = Mesh3D[m][2]; num3 = Mesh3D[m][3]; num4 = Mesh3D[m][4];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				p4 = *(Nodes[num3]);
				
				T2 = p2 - p1; T3 = p3 - p1; T4 = p4-p1;
				A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2)); 
				A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
				
				for(int ii=0;ii<decoup.size();ii++) {
					int ind1 = correspond[decoup[ii][0]]+1;
					int ind2 = correspond[decoup[ii][1]]+1;
					int ind3 = correspond[decoup[ii][2]]+1;
					int ind4 = correspond[decoup[ii][3]]+1;
					
					
					int* tab = new int[5];
					tab[0] = 4; //Mesh3D[m][0];
					tab[1] = Mesh3D[m][ind1]; tab[2] = Mesh3D[m][ind2]; tab[3] = Mesh3D[m][ind3]; tab[4] = Mesh3D[m][ind4];
					meshP1.push_back(tab);
				}
				
			}
			
			// Tetra:
			/*
			for (int m=0; m<Mesh3D.size();m++){
				if (m % 1000 == 0)
					cout<<"Mesh "<<m<<" / "<<Mesh3D.size()<<endl;
				
				num1 = Mesh3D[m][1]; num2 = Mesh3D[m][2]; num3 = Mesh3D[m][3]; num4 = Mesh3D[m][4];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				p4 = *(Nodes[num4]);
				
				T2 = p2 - p1; T3 = p3 - p1; T4 = p4 - p1;
				
				A.set(0,0,T2(0)); A.set(1,0,T2(1)); A.set(2,0,T2(2));
				A.set(0,1,T3(0)); A.set(1,1,T3(1)); A.set(2,1,T3(2)); 
				A.set(0,2,T4(0)); A.set(1,2,T4(1)); A.set(2,2,T4(2));
				
				for(int i=0; i<ordre; i++) {
					for (int j=0; j<ordre-i; j++) {
						for(int k=0; k<ordre-i-j; k++) {
							p1Ref.set(0,k*1.0/ordre); p1Ref.set(1,j*1.0/ordre); p1Ref.set(2,i*1.0/ordre);
							p2Ref.set(0,(k+1)*1.0/ordre); p2Ref.set(1,j*1.0/ordre); p2Ref.set(2,i*1.0/ordre);
							p3Ref.set(0,k*1.0/ordre); p3Ref.set(1,(j+1)*1.0/ordre); p3Ref.set(2,i*1.0/ordre);
							p4Ref.set(0,(k+1)*1.0/ordre); p4Ref.set(1,(j+1)*1.0/ordre); p4Ref.set(2,i*1.0/ordre);
							
							p5Ref.set(0,k*1.0/ordre); p5Ref.set(1,j*1.0/ordre); p5Ref.set(2,(i+1)*1.0/ordre);
							p6Ref.set(0,(k+1)*1.0/ordre); p6Ref.set(1,j*1.0/ordre); p6Ref.set(2,(i+1)*1.0/ordre);
							p7Ref.set(0,k*1.0/ordre); p7Ref.set(1,(j+1)*1.0/ordre); p7Ref.set(2,(i+1)*1.0/ordre);
							p8Ref.set(0,(k+1)*1.0/ordre); p8Ref.set(1,(j+1)*1.0/ordre); p8Ref.set(2,(i+1)*1.0/ordre);
							
							Imp1Ref = A*p1Ref + p1; Imp2Ref = A*p2Ref + p1; 
							Imp3Ref = A*p3Ref + p1; Imp4Ref = A*p4Ref + p1;
							Imp5Ref = A*p5Ref + p1; Imp6Ref = A*p6Ref + p1; 
							Imp7Ref = A*p7Ref + p1; Imp8Ref = A*p8Ref + p1;
							
							
							ind1 = -1; ind2 = -1; ind3 = -1; ind4 = -1; 
							ind5 = -1; ind6 = -1; ind7 = -1; ind8 = -1;
							
	 						for (int ii=0; ii< nbPtPMaille; ii++){
	 							p = *(Nodes[Mesh3D[m][ii+1]]);
								if ((p - Imp1Ref).norme() < 1e-10)
									ind1 = Mesh3D[m][ii+1];
								if ((p - Imp2Ref).norme() < 1e-10)
									ind2 = Mesh3D[m][ii+1];
								if ((p - Imp3Ref).norme() < 1e-10)
									ind3 = Mesh3D[m][ii+1];
								if ((p - Imp4Ref).norme() < 1e-10)
									ind4 = Mesh3D[m][ii+1];
								if ((p - Imp5Ref).norme() < 1e-10)
									ind5 = Mesh3D[m][ii+1];
								if ((p - Imp6Ref).norme() < 1e-10)
									ind6 = Mesh3D[m][ii+1];
								if ((p - Imp7Ref).norme() < 1e-10)
									ind7 = Mesh3D[m][ii+1];
								if ((p - Imp8Ref).norme() < 1e-10)
									ind8 = Mesh3D[m][ii+1];
	 						}
							
							if ((ind1 != -1) and (ind2 != -1) and (ind3 != -1) and (ind5 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind1; tab[2] = ind2; tab[3] = ind3; tab[4] = ind5;
								meshP1.push_back(tab);
							}
							if ((ind5 != -1) and (ind6 != -1) and (ind7 != -1) and (ind3 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind5; tab[2] = ind6; tab[3] = ind7; tab[4] = ind3;
								meshP1.push_back(tab);
								
							}
							if ((ind2 != -1) and (ind5 != -1) and (ind6 != -1) and (ind3 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind2; tab[2] = ind5; tab[3] = ind6; tab[4] = ind3;
								meshP1.push_back(tab);
								
							}
							
							
							if ((ind2 != -1) and (ind3 != -1) and (ind4 != -1) and (ind6 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind2; tab[2] = ind3; tab[3] = ind4; tab[4] = ind6;
								meshP1.push_back(tab);
							}
							if ((ind4 != -1) and (ind6 != -1) and (ind7 != -1) and (ind8 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind4; tab[2] = ind6; tab[3] = ind7; tab[4] = ind8;
								meshP1.push_back(tab);
							}
							if ((ind3 != -1) and (ind6 != -1) and (ind7 != -1) and (ind4 != -1) ) {
								int* tab = new int[5];
								tab[0] = 4; 
								tab[1] = ind3; tab[2] = ind6; tab[3] = ind7; tab[4] = ind4;
								meshP1.push_back(tab);
							}
						}
					}
				}
			}
			*/
			return meshP1;
		}
		
		vector<int*> getMesh1Din(Vecteur a, Vecteur b, int ordre, int nbPtMaille) {
			vector<int*> res;
			Vecteur p1(3); Vecteur p2(3);
			Vecteur T(3);
			T = (b-a) / (b-a).norme();
			double coordb = (b-a)*T;
			
			if(T.norme() <1e-8) {
				cout<<"Problemme T !!!"<<endl;
				int PAUSE;
				cin>>PAUSE;
			}

			vector< vector<int> > MatAdj;
			for(int i=0;i<getNNodes();i++) {
				vector<int> tp;
				MatAdj.push_back(tp);
			}

			for (int m=0; m<Mesh2D.size();m++){
				// int num1 = Mesh2D[m][0];
// 				int num2 = Mesh2D[m][1];
// 				int num3 = Mesh2D[m][2];

				for(int a1=0;a1<3;a1++){
					for(int a2=a1+1;a2<3;a2++){
						int num1 = Mesh2D[m][a1+1];
						int num2 = Mesh2D[m][a2+1];

						if(num1 > num2) {
							int num = num2;
							num2 = num1;
							num1 = num;
						}
						if(num1 == num2){
							int PAUSE;
							cout<<a1+1<<"	"<<a2+1<<endl;
							cout<<Mesh2D[m][a1+1]<<"	"<<Mesh2D[m][a2+1]<<endl;
							cout<<"Pb !!! "<<num1<<"	"<<num2<<endl;
							cin>>PAUSE;
						}
						// On verifie qu'on a pas déjà parcouru l'arete
						int found = -1;
						for(int i=0;i<MatAdj[num1].size();i++){
							if(MatAdj[num1][i] == num2){
								found = i;
								break;
							}
						}
						if (found == -1) {
							MatAdj[num1].push_back(num2);
							MatAdj[num2].push_back(num1);
						}
						if(found==-1){
							p1 = *(Nodes[num1]);
							p2 = *(Nodes[num2]);
							
							if( isIn(p1,a,b) and isIn(p2,a,b) ) {
								double coordp1 = (p1-a)*T;
								double coordp2 = (p2-a)*T;
								if (coordp1 >=0 and coordp1 <=coordb and coordp2 >= 0 and coordp2 <=coordb) {
									// les deux points p1 et p2 sont dans le segment [a,b]
									int* tab = new int[ordre+2];
									tab[0] = 1;
									tab[1] = num1;
									tab[2] = num2;
									int cpt = 3;
									for(int j=4;j<=nbPtMaille;j++) {
										int ind = Mesh2D[m][j];
										Vecteur p(3);
										p = *(Nodes[ind]);
										Vecteur T1(3);
										T1 = (p2-p1) / (p2-p1).norme();
										double lenn = (p2-p1).norme();
										if (isIn(p,p1,p2) ){ // and cpt<ordre+2
											double coordp = (p-p1)*T1;
											if(coordp > 0 and coordp < lenn) {
												tab[cpt] = ind;
												cpt++;
											}
											
										}
										// if(cpt>=ordre+2){
// 											cout<<"Probleme !!! "<<cpt<<" "<<ordre+2<<endl;
// 											int PAUSE;
// 											cin>>PAUSE;
// 										}
									}
									res.push_back(tab);
								}
							}
						}
					}
				}
			}
			return res;
		}
		
		
		vector<int> getNodesin(Vecteur a, Vecteur b) {
			vector<int> res;
			for(int i=0; i<getNNodes();i++) {
				Vecteur p(3);
				p = *(Nodes[i]);
				if(isIn(a,b,p)) 
					res.push_back(i);
			}
			return res;
		}
		
		// Recupere le numero des mailles dans la boule de centre cc et de rayon R
		vector<int> getMesh2Din(Vecteur cc, double R, int norm) {
			vector<int> res;
			Vecteur p1(3); Vecteur p2(3); Vecteur p3(3);
			Vecteur Bary(3);

			for (int m=0; m<Mesh2D.size();m++){
				int num1 = Mesh2D[m][1];
				int num2 = Mesh2D[m][2];
				int num3 = Mesh2D[m][3];
				
				p1 = *(Nodes[num1]);
				p2 = *(Nodes[num2]);
				p3 = *(Nodes[num3]);
				
				Bary = (p1 + p2 + p3)/3.0;
				
				double distBC = ((Bary - cc) * (Bary - cc));
				if(norm == 1)
					distBC = fabs(Bary(0)-cc(0)) + fabs(Bary(1)-cc(1)); 
				if(distBC < R*R)
					res.push_back(m);
			}
			return res;
		}
		
		
		vector<int> triangConn(int l, int c) {
			
			if(!isConnec2D) {
				isConnec2D = true;
				
				mesh1Connect2D = new MatriceL(Nodes.size());
				mesh2Connect2D = new MatriceL(Nodes.size());
				
				
				for(int m=0;m<Mesh2D.size();m++) {
					int i1 = Mesh2D[m][1];
					int i2 = Mesh2D[m][2];
					int i3 = Mesh2D[m][3];
					
					if((*mesh1Connect2D)(i1,i2) != 0) {
						mesh2Connect2D->set(i1,i2,m+1);
						mesh2Connect2D->set(i2,i1,m+1);
					}
					else {
						mesh1Connect2D->set(i1,i2,m+1);
						mesh1Connect2D->set(i2,i1,m+1);
					}
					
					if((*mesh1Connect2D)(i1,i3) != 0) {
						mesh2Connect2D->set(i1,i3,m+1);
						mesh2Connect2D->set(i3,i1,m+1);
					}
					else {
						mesh1Connect2D->set(i1,i3,m+1);
						mesh1Connect2D->set(i3,i1,m+1);
					}
					
					if((*mesh1Connect2D)(i2,i3) != 0) {
						mesh2Connect2D->set(i2,i3,m+1);
						mesh2Connect2D->set(i3,i2,m+1);
					}
					else {
						mesh1Connect2D->set(i2,i3,m+1);
						mesh1Connect2D->set(i3,i2,m+1);
					}
				}
			}
			
			vector<int> res;
			
			if((*mesh1Connect2D)(l,c) > 0){
				res.push_back(int((*mesh1Connect2D)(l,c)-1));
			}
			if((*mesh2Connect2D)(l,c) > 0){
				res.push_back(int((*mesh2Connect2D)(l,c)-1));
			}
			
			return res;
		}
		
		
		vector<int*> getBoundary2D(int ordre, int nbPtMaille, vector<int>& lienBordVMesh) { // , vector<int>* resMeshBord
			vector<int*> res;
			// resMeshBord -> Maille qui longent le bord
			
			Vecteur p1(3); Vecteur p2(3);
			Vecteur T(3);
			
			MatriceL lien(getNNodes());

			//vector< vector<int> > MatAdj;
			// for(int i=0;i<getNNodes();i++) {
// 				vector<int> tp;
// 				MatAdj.push_back(tp);
// 			}

			for (int m=0; m<Mesh2D.size();m++){
				// int num1 = Mesh2D[m][0];
// 				int num2 = Mesh2D[m][1];
// 				int num3 = Mesh2D[m][2];

				for(int a1=0;a1<3;a1++){
					for(int a2=a1+1;a2<3;a2++){
						int num1 = Mesh2D[m][a1+1];
						int num2 = Mesh2D[m][a2+1];

						if(num1 > num2) {
							int num = num2;
							num2 = num1;
							num1 = num;
						}
						if(num1 == num2 or num1 < 0 or num2 < 0){
							int PAUSE;
							cout<<a1+1<<"	"<<a2+1<<endl;
							cout<<Mesh2D[m][a1+1]<<"	"<<Mesh2D[m][a2+1]<<endl;
							cout<<"Pb !!! "<<num1<<"	"<<num2<<endl;
							cin>>PAUSE;
						}
						
						// Compte le nombre de fois qu'on passe par une arete
						lien.set(num1,num2,lien(num1,num2)+1);
						lien.set(num2,num1,lien(num2,num1)+1);
					}
				}
			}
			
			
			lienBordVMesh.clear();
			for(int m=0;m<Mesh2D.size();m++) {
				for(int a1=0;a1<3;a1++){
					for(int a2=a1+1;a2<3;a2++){
						int num1 = Mesh2D[m][a1+1];
						int num2 = Mesh2D[m][a2+1];
						
						if(num1 == num2){
							int PAUSE;
							cout<<a1+1<<"	"<<a2+1<<endl;
							cout<<Mesh2D[m][a1+1]<<"	"<<Mesh2D[m][a2+1]<<endl;
							cout<<"Pb !!! "<<num1<<"	"<<num2<<endl;
							cin>>PAUSE;
						}
						
						int nbMailleAdj = lien(num1,num2);
						if(nbMailleAdj == 1) { // On a une maille ayant un coté sur le bord
							/*
							int alreadyfound = -1;
							for(int mm=0;mm<resMeshBord->size();mm++) {
								if( (*resMeshBord)[mm] == m) {
									alreadyfound = 1;
									break;
								}
							}
							if(alreadyfound == -1)
								resMeshBord->push_back(m);
							*/
							Vecteur aa = *Nodes[num1];
							Vecteur bb = *Nodes[num2];
							int* tab = new int[ordre+2];
							tab[0] = 1.0;
							tab[1] = num1; 
							tab[2] = num2;
							int cpt = 3;
							for(int ipt=0;ipt<nbPtMaille;ipt++) {
								int ind = Mesh2D[m][ipt+1];
								if(ind != num1 and ind != num2) {
									Vecteur tmp = *Nodes[ind];
									
									if(isIn(tmp,aa,bb)) {
										if(cpt > ordre+1) {
											cout<<"Probleme bords !!!"<<endl;
											cout<<ind<<" "<<num1<<" "<<num2<<endl;
											int PAUSE;
											cin>>PAUSE;
										}
										tab[cpt] = ind;
										cpt++;
										//cout<<"AIIICIIII"<<endl;
									}	
								}	
							}
							if(cpt<= ordre){
								cout<<"Probleme bords pas assez de terme !!!"<<endl;
								int PAUSE;
								cin>>PAUSE;
							}
							res.push_back(tab);	
							lienBordVMesh.push_back(m);
						}
					}
				}	
			}
			
			//std::cout<<lienBordVMesh.size()<<"	"<<res.size()<<endl;
			return res;
		}
		
	
};
#endif
