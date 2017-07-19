//Functions for the BPA algorithm





#include "triangulations/MeshGrowing.h"
#include <iostream>


    //Some  macros

    //quando un triangolo viene rimosso bisogno aggiornare la struttura
    //dati in modo che gli edge "lasciati soli" abbiano un orientazione
    //coerente con l'unico triangolo che li contiene. E' questo un problema
    //che si verifica nellegenerazione di mesh per avanzamento di fronte quando
    //si va a cancellare un triangolo. L'edge solo deve essere un edge al fronte
    #define Reset_e2t(ide)\
            if (e[ide].t2!=idt)\
            {e[ide].t1=e[ide].t2;}\
                    e[ide].t2=-1;\
                    if(e[ide].t1>=0)\
                    {T1=e[ide].t1;\
                             if (t[T1].p1==e[ide].p2 && t[T1].p2==e[ide].p1){Swap(ide)}\
                             else if (t[T1].p2==e[ide].p2 &&t[T1].p3==e[ide].p1){Swap(ide)}\
                             else if (t[T1].p3==e[ide].p2 && t[T1].p1==e[ide].p1){Swap(ide)}\
                    }

//Swap the edge to preserve orientation
                    #define Swap(ide)\
                            temp=e[ide].p1;\
                            e[ide].p1=e[ide].p2;\
                                    e[ide].p2=temp;



////boundary edges not deleted
#define  IsBound(e1) (e[e1].t1>=0 && e[e1].t2<0)

//look for over boundary triagnles
#define LookForNMV_OverBound_Triangle(e1)\
  if (IsBound(e1))\
   {    if(nbe[e[e1].p1]>2 || nbe[e[e1].p2]>2){PostP_Delete_NMV_OverBound_Triangle(e[e1].t1);}\
    }
//looks for boundary triagnles belonging to the point idp
#define LookForNMV_BoundTriangle_idp(e1)\
if (IsBound(e1))\
{T2=e[e1].t1;\
 if (t[T2].p1==idp || t[T2].p2==idp || t[T2].p3==idp)\
 {PostP_Delete_NMV_BoundTriangle_idp(T2, idp);}\
}




//constructor
MESHGROWING::MESHGROWING() {

    p=NULL;
    N=0;
    R=0;
    Ndel=0;
    Ntdel=0;

    //entities counters
    countt=0;
    counte=0;


    //entities
    t=NULL;
    e=NULL;
    queue=NULL;
    NE=NULL;
    NT=NULL;
    nbe=NULL;
    FlatToll=NULL;

    //maximumnumbers of entities

    MAXE=N*3.2;
    MAXT=N*2.1;
    Ntdel=0;

    InitParameters();
}

//destructor
MESHGROWING::~MESHGROWING() {
    Memory_Deallocate(Mem_All);
}

void MESHGROWING::MakePublic()
//copy some private to public memebrs
{
    nt=countt;//passing to public member
    ne=counte;
}
void MESHGROWING::InitParameters(){//set default values for internalparameters

    //Internal P (Deafault values in the constructor)




    //BPA
    P.BPA_TooSharp_Torus=-.5;
    P.BPA_TooSharp_SR=0;
    P.BPA_Max_BCAngle=1.1;//Maximum BC angle
    P.BPA_FlatnessToll=-.7;
	P.BPA_ReverseNormal=false;//false =default normal outside


    //SCB
    P.SCB_TooSharp=-.5;//flatnnes test
    P.SCB_SR_PriorityLevel=2;//Priority levels for search radius
    P.SCB_Dotp_PriorityLevel=3;// 3 per il fltness test
    //Attenzione la modifica del parametro sopra comporta il cambiamento dell'array in Get_SCB_Triangle


    //ALL

    P.NN_Filter_Cutdist=20;//maeandist/cutdist will be used in the filter. <=0 means disabled
    P.Max_Edge_Length=HUGE_VAL;;//maximum allowed edge length
    P.EnablePreProcessing=true;
    P.EnablePostProcessing=true;

}


void MESHGROWING::PrintParameters(int mode){//set default values for internalparameters
    //mode==0 Prints BPA Parameters
    //mode==1 Prints SCB Parameters

    cout<<"Algorithm started with the following parameters:"<<endl;
    if (mode==0){//BPA Parameters
        cout<<"-TooSharp_Torus: "<<P.BPA_TooSharp_Torus<<endl;
        cout<<"-TooSharp_SR: "<<P.BPA_TooSharp_SR<<endl;
        cout<<"-Max_BCAngle: "<<P.BPA_Max_BCAngle<<endl;
        cout<<"-FlatnessToll: "<<P.BPA_FlatnessToll<<endl;

    }
    else if(mode==1){//SCB PARAMETERS
        cout<<"-TooSharp: "<<P.SCB_TooSharp<<endl;
        cout<<"-SR_PriorityLevel: "<<P.SCB_SR_PriorityLevel<<endl;
        cout<<"-Dotp_PriorityLevel: "<<P.SCB_Dotp_PriorityLevel<<endl;
    }
    else{
        Error("Invalid Print Parameters Mode");
    }

    //ALL PARAMTERS
    cout<<"-NN_Filter_Cutdist: "<<P.NN_Filter_Cutdist<<endl;
    cout<<"-Max_Edge_Length: "<<P.Max_Edge_Length<<endl;
    cout<<"-EnablePreProcessing: "<<P.EnablePreProcessing<<endl;
    cout<<"-EnablePostProcessing: "<<P.EnablePostProcessing<<endl;



}



void MESHGROWING::Memory_Deallocate(DeallocateMode Mode) {
    //Default Mode is Mem_All

    //Free Memory (no need of errors checks)

    switch(Mode)
    {case Mem_All:

         Deallocate(&t);
         Deallocate(&e);

         break;

        case Mem_Save_t:

            Deallocate(&e);

        case Mem_Save_t_e:

            break;

            //these are always deallocated
            Deallocate(&NT);
            Deallocate(&NE);
            Deallocate(&queue);
            SDS.Deallocate();
            EPMap.Deallocate();

    }


}




void MESHGROWING::Memory_Allocate() {

    MAXE=N*3.2;
    MAXT=N*2.1;
    //Oversize the triangles arrays and call the tesselletion routine
    Allocate(&t, MAXT);
    Allocate(&e, MAXE);
    Allocate(&queue, MAXE);
    AllocateAndInit(&NE, N, int16_t(0));
    AllocateAndInit(&NT, N, int16_t(0));


    EPMap.BuildPointEdgeMap(N, MAXE); //EPmap

    //inizilizza alcuni valori
    int i;
    for (i=0;i<MAXE;i++)
    {e[i].t2=-1;}//flag all edges as boundary

}



void MESHGROWING::SeedTriangle(bool BPAtest=false) { //aggiungere l'error check nel caso il fronte non parte

    //Attenzione questa funzione funziona solo per il primo triangolo non per il front refresh

    //Questa funzione ha lo scopo di creare il primo triangolo in modo automatico.
    //Viene adottata la seguente strategia:
    // IL primo triangolo deve avere la normale orientata all'esterno. Il volume formato dal cilindro:
    //-il cui raggio è uguale a quelllo della sfera D2.5D
    //-il cui asse è parallelo alla normale del triangolo e passante per il circocentro
    //-il cilindro è tagliato dal piano passante per il triagnolo, la parte considerata è quella in cui punta la normale
    //QUESTO VOLUME DEVE ESSERE PRIVO DI PUNTI
    // -nel BPA questo triangolo deve anche essere BE (Ball Empty) dal lato della normale


    //COSTANT DEFINITION


    //MEMORY ALLOCATION

      int i,j;
    double mindist;//minum distance
    Coord3D pm, cc, BC;//Punto medio

    double  r2;//circumcenter and circum radius (squared and normal)
    int dtest, betest;
    bool GoodTriangle;
    bool up0;


    Coord3D tnorm;

    cout<<"Seed Triangle...";

    for (i=0;i<N;i++)//loop trough all points until a good start triangle has been found
    {
        P1=i;
        //Find the closet to i point
        SDS.SearchClosestExclusive( &p[P1], &P2, &mindist, P1);


        //Find all the points in k range of edge midpoint
        MidPoint(p[P1]., p[P2]., pm.)
        SDS.SearchRadius( &pm, K*mindist);


        //Stampa i punti trovati
        /*for (j=0;j<SDS.npts;j++)
         * {
         * cout<<j<<" "<< SDS.idStore[j]<< endl;
         * }*/

        GoodTriangle=false;//inizializzo nel caso non ci siano punti nel range

        //per ogni punto verifica se si può troavre un delaunay2_5D triangle
        for (j=0;j<SDS.npts;j++) {
            P4=SDS.idStore[j];

            if (P4==P1 || P4==P2)//non accettare il terzo punto uguale al aprimo o al secondo
            {
                continue;
            }
            //Get circumcenters and circumradius


            // la normale è calcolata come: cross(v(idc-P1),(v(id3-P1))
            //gli input della funzione sono P1 P2 P4 per calcolare cross(v(P4-P1),(v(P2-P1))
            // Il triangolo dovrà essere [P1 P4 P2]=[P1 idc id3]
            CircumCenter(&p[P1],  &p[P4], &p[P2], &tnorm, &r2, &cc);
            dtest=SDS.EmptyBallTest( &cc, .9999999*r2);//run Dtest

            GoodTriangle=false;//Prima del Dtest il Triangolo è presunto cattivo


            if(dtest==-1 && CylinderTest(cc, r2, tnorm, &up0)) //se supera il Dtest e il cylinder test
            {
                if (BPAtest) {//in the ball pivoting we also need to test if the triagnle is ball empty
                    if(up0){Reverse(tnorm.);}
                    if (!BallCenter(&cc, &r2, &tnorm, &BC)){continue;}//continuee if traignle is too big
                    //test if the traignle is BPA
                    betest=SDS.EmptyBallTest(&BC, R*R*.999999);
                    if (betest<0 ) {
                        GoodTriangle=true;break;}
                }
                else
                {GoodTriangle=true;break;}
            }

        }//loop trough all points in k-range


        if(GoodTriangle) {
            //aggiungi il nuovo triangolo

			if (P.BPA_ReverseNormal)
			{  cout<<"Reversing Normal: Inner Ball"<<endl;
				up0=!up0;}//reverse normal make the algorithm behave at the contrary
           
            //check wheter reverse normal
            if(up0)//all points were over the triangle
            {
                // Reverse(t[0].tnorm.)//la normale punta verso l'alto

                //per cambiare l'orientazione del trianglolo scambio p1 e p2
                UpdateFront_Seed(P2, P1, P4);

            }
			else
			{ UpdateFront_Seed(P1, P2, P4);}

            cout<<"Found: "<<t[0].p1<<" "<<t[0].p2<<" "<<t[0].p3<<endl;
            return;//found a good triangle return
        }//good Triangle Found


    }//loop trough al points until a good triangle has been found

    Error("Unable to find Seed Triangle");
    return;//no seed triangle found
}

bool MESHGROWING::CylinderTest(Coord3D cc, double r2, Coord3D tnorm, bool* up0) {//testa se tutti ipunti sono sopra o sotto il cilindro definito da cc,r2 e tnorm
    //se sono tutti sopra up0=true altrimenti falso

    //calcolo la distanza dei punti dall'asse del cilindro nel seguente modo:
    //-calcolo la distanza punto circocentro (ipotenusa)
    //-calcolo la distanza del primo cateto come dot(v(p-cc),n), con n normale del triangolo
    //-calcolo la distanza con il teorema di pitagora  sqrt(ipotenusa^2-cateto1^2)
    int  q;
    bool* inside=NULL;
    double r, dist, leng;
    Coord3D v21;//vettori per il calcolo della normale
    bool up, GoodTriangle;
    double d;//coefficeinte del piano secondo l'eq ax+by+cz+d;


    Allocate(&inside, N);//non lo inizializzo lo scrivo poi come vero e falso
    GoodTriangle=false;

    r=sqrt(r2);//raggio SCB
    for (q=0;q<N;q++) {

        Distance(cc., p[q]., dist);//distanza punto circocentro (mindist)
        DiffPoints(p[q]., cc. , v21.);//vettore (p-cc) salvato in v21.
        DotProduct(v21., tnorm., leng);//lunghezza cateto (leng)
        leng=sqrt(dist*dist-leng*leng);//distanza dall'asse (leng)

        if (leng<r)// il punto è nel cilindro?
        {
            inside[q]=true;
        }
        else
        {inside[q]=false;}
    }
    //i punti del triangolo giaciono sul piano ma per evitare errori numerici li consdiero fuori
    inside[P1]=false;
    inside[P2]=false;
    inside[P4]=false;

    //ora dobbiamo verificare che i punti stiano tutti dalla stessa parte del triangolo,
    //se ciò è vero la direzione della normale sarà quella che non "guarda" neesuno punto


    DotProduct(cc., tnorm., d);//coefficiente del piano

    //ora vedo se il primo punto non inside si trova sotto o sopra

    q=0;
    while (q<N) {//faccio il loop fino ad N ma in realtà esco una volta trovato il primo
        if (inside[q]) {
            DotProduct(p[q]., tnorm., leng);// Dot(p,n)
            *up0=(leng-d)>0; // up >0 punto sopra il piano
            q++;
            break;//bail out
        }
        q++;
    }


    //Se non ci sono punti nel cilindro dobbiamo passare al prossimo punto
    if (q>=N)//no point in the cylinder
    {return GoodTriangle;}//try next point

    //ora tutti i punti devono avere la stessa orientazione
    while (q<N) {
        if (inside[q]) {
            DotProduct(p[q]., tnorm., leng);// Dot(p,n)
            up=(leng-d)>0; // up >0 punto sopra il piano

            if (up!=*up0)//orientazione diversa
            {
                //GoodTriangle=false;
                break;//bail out
            }
        }
        q++;
    }


    if (q>=N)//We reached the end the triagnle is good
    {GoodTriangle=true;}

//deallocate
    Deallocate(&inside);

    return GoodTriangle;
}

void MESHGROWING::UpdateFront_Seed(int i, int idc, int id3) {
    //Update the data structure with the first triangle


    t[0].p1=i;
    t[0].p2=idc;
    t[0].p3=id3;

    //New edges (preserving orientation)
    e[0].p1=t[0].p1;
    e[0].p2=t[0].p2;
    EPMap.AddEdge(t[0].p1, t[0].p2, 0);//add edge to the map

    e[1].p1=t[0].p2;
    e[1].p2=t[0].p3;
    EPMap.AddEdge(t[0].p3, t[0].p2, 1);//add edge to the map


    e[2].p1=t[0].p3;
    e[2].p2=t[0].p1;
    EPMap.AddEdge(t[0].p3, t[0].p1, 2);//add edge to the map


    //new triangles edges connectivity
    e[0].t1=0;
    e[1].t1=0;
    e[2].t1=0;

    //add edges to queue
    Pqueue.PushHigh(0);
    Pqueue.PushHigh(1);
    Pqueue.PushHigh(2);


    //t[0].e1=0;
    //t[0].e2=1;
    //t[0].e3=2;

    //counters
    countt=1;
    counte=3;


    //for not manifold topology
    NE[i]=NE[idc]=NE[id3]=2;
    NT[i]=NT[idc]=NT[id3]=1;



}



//The tesselation routine
void MESHGROWING::RunFront_BPA() {

    int i=0;//queue iterator
    int idedge1, idedge2;//id dei presunti nuovi edge
    int NumSR=0;
    int NumTorus=0;
    //ADVANCING FRONT
    cout<<"Running Front... ";
    //pop an elemnt from the priority queue at each iteration
    while (counte<MAXE && Pqueue.Pop(&ide)) {


        //cout<<"i: "<<i<<" ide: "<<ide<<endl; //il bug is verifica per i=9024 ide=9023;

        if (e[ide].t2>=0){continue;}//edge is no more on front skip
#ifdef _DEBUG
        if(e[ide].p1==21116 && e[ide].p2==18740)
        {cout<<"edge under debug"<<endl;}
#endif
        //GetTriangle_BPA_Torus();
        //GetTriangle_BPA_SR();
        if (e[ide].t2==-1)
        {GetTriangle_BPA_SR();
         if (P4==-1)
         {Pqueue.Push(ide, 1);// add to queue with lower priority
          e[ide].t2--;}
         else{NumSR++;}
        }
        else
        {GetTriangle_BPA_Torus();if(P4>=0){NumTorus++;}}

        if (P4<0) {continue;}//no candidate

        //CHECK EDGES CONNECTIVITY
        idedge1=EPMap.GetEdge(P1, P4);
        if(!CheckEdgeConformity(idedge1, P1)){continue;}//edge not conform
        idedge2=EPMap.GetEdge(P4, P2);
        if(!CheckEdgeConformity(idedge2, P4)){continue;}//edge not conform

        //Se è stato trovato un buon triangolo dobbiamo aggiornare la struttura dati
        UpdateFront(idedge1, idedge2);


    }

    cout<<countt<<" Triangles generated: "<<"NumSr :"<<NumSR<<" NumTorus: "<<NumTorus<<endl;
}

// Tessellation strategies
void MESHGROWING::BallPivoting(double inputr) {
    //Classical Ball Pivoting without points normals
    int flag;
    double meandist;

    PrintParameters(0);//Print BPA Parameters

    R=inputr;//set the radius value

    Pqueue.Set(queue, 2);//setting the queue with 2 priority leveles (SR or Torus triagnle)

    SDS.step=R;//setting the SDS step to R



    cout<<"Building Search Data Structure...";
    flag=SDS.BuildSDS(p, N);if(flag<0)Error("Error duting SDS construction");//SDS uses -1 for out of memory
     cout<<"Done"<<endl;


   cout<<endl<<"PREPROCESSING"<<endl;
     PreP_NN_MeanDist(&meandist);//conpute NN graph

    if(P.NN_Filter_Cutdist>0)PreP_NN_Filter(meandist/P.NN_Filter_Cutdist);//removing points too close too each other

    //auto computation of the ball radius
    if(R==0){R=meandist*4;}


    cout<<"Radius set to: "<<R<<endl;

    cout<<endl<<"MESHING"<<endl;
    SeedTriangle(true);//seed triagnle with BPATest

 

    RunFront_BPA();

    //PostP
    if (P.EnablePostProcessing)
    { cout<<endl<<"POSTPROCESSING"<<endl;
     PostP_PostProcessing();//the gateway routine that inlcudes postprocessing operations
    }

    MakePublic();
    Pqueue.Reset();



};

void MESHGROWING::SCBMesher() {
    //SCB MEsher algorithm
    int flag;
	double meandist;

    cout<<endl<<"PRE-PROCESSING"<<endl;
    //setting SDS
	cout<<"Building search data structure"<<endl;
    SDS.density=.5;
    flag=SDS.BuildSDS(p, N);if(flag<0)Error("Error during SDS Construction");////SDS uses -1 for out of memory



     PreP_NN_MeanDist(&meandist);//conpute NN graph
     PreP_NN_Filter(meandist/P.NN_Filter_Cutdist);//removing points too close too each other

    Pqueue.Set(queue, P.SCB_Dotp_PriorityLevel+P.SCB_SR_PriorityLevel);//setting the queue


	 cout<<endl<<"MESHING"<<endl;
    SeedTriangle();

    RunFront_SCB();

	 cout<<endl<<"POST-PROCESSING"<<endl;
	PostP_PostProcessing();

    MakePublic();//pass to public members

	Pqueue.Reset();//Reset the queue

};

void MESHGROWING::PostP_PostProcessing()//the gataway routine that inlcudes postprocessing operations
{

	 PostP_FillTriplets();

     PostP_HealNMV();

	 //EdgeQueue_Fill_With_Boundary(2);
	 //PostP_HoleFiller();
     //Here we could delete NT and NE but that's just a little memory and we wont do it!
     if(Ntdel>0){PostP_Restore_t();}//if there are deleted triagnles restore connections
}

void MESHGROWING::EdgeQueue_Fill_With_Boundary(int sr_value)//fill the queue with boundary edges with the given search radius value
{   int j;
    Pqueue.Empty();//empty the priority queue
	 
			for (j=0;j<counte;j++) {
                if IsBound(j)//edge di boundary non cancellato
				{ 
					Pqueue.PushHigh(j);//add edge to the queue
					e[j].t2=-sr_value;//resetto il search radius
                }
            }
}
	 void MESHGROWING::BuildTorus(Torus3D* Torus) {//Cosntruisce il toro di pivoting partendo dall'edge P1-P2
    //memory allocation


    double leng;



    //trova il punto medio dell'edge (centro del toro)
    MidPoint(p[P1]., p[P2]., Torus->Center.);

//Asse del toro
    DiffPoints(p[P2]., p[P1]., Torus->Axis.);
    Normalize(Torus->Axis., leng);

    //Cooeff d del toro
    DotProduct(Torus->Center., Torus->Axis., Torus->d);

    //lunghezza dell'edge
    Distance(p[P2]., p[P1]., leng)
    //Calcolare la lunghezza di proiezione in base al raggio della palla
    Torus->R=sqrt(R*R-leng*leng*0.25);//occhio R da solo è il raggio della palla (il raggio del toro è torus.R
    Torus->r=R;



}



bool MESHGROWING::InsideTorus(Torus3D* Torus, Coord3D* point) {//Detects if a point lies inside a torus

    double a, b;
    double dist, dista;

    //Distanza punto centro toro
    Distance(point->, Torus->Center., dista);

    //distanza punto piano toro
    DotProduct(point->, Torus->Axis., dist);
    dist=dist-Torus->d;


    a=sqrt(dista*dista-dist*dist);
    b=Torus->R-a;

    return sqrt(dist*dist+b*b)<=Torus->r;

}

void MESHGROWING::GetCandidates_Torus(Torus3D* Torus) {
    double Cuboid[6];
    double sr=Torus->R+Torus->r;

    //building the Cuboid
    Cuboid[0]=Torus->Center.x-sr;Cuboid[1]=Torus->Center.x+sr;
    Cuboid[2]=Torus->Center.y-sr;Cuboid[3]=Torus->Center.y+sr;
    Cuboid[4]=Torus->Center.z-sr;Cuboid[5]=Torus->Center.z+sr;

    SDS.GetPointsInRange(&Cuboid[0]);


    //Standard search with search radius
    // Coord3D sp;
    // SearchPoint(&p[P1],&p[P2],&t[T1].tnorm,&sp,&sr,1);
    // SDS.SearchRadius(&sp,.999*sr);



    //SearchPointBPA(&p[P1],&p[P2],&t[T1].tnorm,&sp);
    // SDS.SearchRadius(&sp,R);



}


void MESHGROWING::SelectCandidate_BPA_Torus(Torus3D* Torus)
//Select Canditae Point among a set of candidate points for the front edge ide
{
    int j, P3, Ptemp;
    double dist;

    Coord3D BC;//Ball center
    Coord3D m;//Midpoint
    double Angle, Dotp;//we must choose the minimum angle ball center
    //ANgle is a modified dotprduct so we must choose the maximum "Angle"

    double  MaxAngle=-2;
    double StartAngle;//to aovid the ball to rotate at the contrary
    Coord3D T1Norm, T2Norm, tempnorm;
    double cr2;//squared circumradius
    Coord3D cc;//circumcenter


#ifdef _DEBUG
Distance(p[P1]., p[P2]., dist);
if (dist>2*R)
{printf("Edge too long found\n");return;}
#endif
Setdiff(t[T1]., e[ide]., P3);
MidPoint(p[P1]., p[P2]., m.);//get the midpoint



//Get the starting position of the ball
CircumCenter(&p[P1], &p[P3], &p[P2], &T1Norm, &cr2, &cc);//normal computed as cross(v(p2-p1),v(p4-p1));
BallCenter(&cc, &cr2, &T1Norm, &BC);
StartAngle=BCAngle(&m, &BC, &T1Norm);//compute angle (gamma) with T1


//Start angle can not be greater thean the set tolerance this should avoid the front inversion in presence of holes
//if(StartAngle>P.BPA_Max_BCAngle){StartAngle=P.BPA_Max_BCAngle;}


//Get traingle data
DiffPoints(p[P2]., p[P1]., v21.);



//loop trough candidates
for (j=0;j<SDS.npts;j++)
{ Ptemp=SDS.idStore[j];


  if(P1==Ptemp || P2==Ptemp || P3==Ptemp) {continue;}//verte belongs to T1
  if(!InsideTorus(Torus, &p[Ptemp])){continue;}//point outside torus

  //getting the normal of the new triagnle
  DiffPoints(p[Ptemp]., p[P1]., v41.);
  CrossProduct(v41., v21., tempnorm.);
  Normalize(tempnorm., dist);


  CircumCenter_Fast(&p[P1], &p[P2], &p[Ptemp], &tempnorm, &cr2, &cc);

  if (!BallCenter(&cc, &cr2, &tempnorm, &BC)){continue;}//continuee if traignle is too big


  //Distance(p[P3].,BC.,dist);
  //if (dist<R*.999999)
  //    {continue;}//unnatural pivoting



  Angle=BCAngle(&m, &BC, &T1Norm);//compute angle (gamma) with T1
  if (Angle>MaxAngle && Angle<StartAngle && Angle>-StartAngle) {
      MaxAngle=Angle;
      P4=Ptemp;
      T2Norm=tempnorm;//compute normal
  }

}

//checking if the ball is truly empty




if(P4>=0){
    DotProduct( T2Norm., T1Norm., Dotp);
    if (Dotp<P.BPA_TooSharp_Torus)
    {{P4=-3;return;}}
    if(NE[P4]==NT[P4] && NE[P4]!=0)
    {P4=-4;return;}//not manifold vertex
}
#ifdef  _DEBUG //check for not BPA triagnles
if(P4>=0) {


    CircumCenter(&p[P1], &p[P2], &p[P4], &tempnorm, &cr2, &cc);
    BallCenter(&cc, &cr2, &T2Norm, &BC);
    int test;
    test=SDS.EmptyBallTest(&BC, R*R*.999);
    if (test>=0 )//|| MaxAngle>1.2
    {

        double dist;
        Distance(BC., p[test]., dist);
        cout<<" P1="<<P1<<" P2="<<P2<<" P3="<<P3<<" P4="<<P4<<" test="<<test<<" Distance(BC-test)="<<dist<<endl;
        P4=-1;//exclude this triagnle
    }
}
#endif








}



bool MESHGROWING::CheckEdgeConformity(int idedge, int point) {
    if (idedge>=0)//solo se l'edge esiste già
    {if(e[idedge].t2>=0){return false;}//posto occupato
     else {

         //controlla l'orientamento
         if (e[idedge].p1==point){return false;}//degenere

         //Checking the angle
         //double Angle;
         //DotProduct( t[countt].tnorm.,t[e[idedge].t1].tnorm.,Angle);
         // if (Angle<-0.5){return false;}//{continue;}

     }
    }
    return true;


}

void MESHGROWING::UpdateFront(int idedge1, int idedge2) {
    //Check the first edge
    if  (idedge1<0)// %if edge1 is new
    {
        e[counte].p1=P1;//add new edge
        e[counte].p2=P4;//add new edge
        e[counte].t1=countt;//add new triangle
        Pqueue.PushHigh(counte);//add to queue


        NE[P1]++;
        NE[P4]++;

        EPMap.AddEdge(P1, P4, counte);//add edge to the map

        //t[countt].e1=counte;

        counte++;//increase counter edge

    }
    else//edge is old just remove it from front
    {
        //t[countt].e1=idedge1;
        e[idedge1].t2=countt;//add new triangle
    }


    //Check the second edge
    if  (idedge2<0)// %if edge1 is new
    {
        e[counte].p1=P4;//add new edge
        e[counte].p2=P2;//add new edge
        e[counte].t1=countt;//add new triangle
        Pqueue.PushHigh(counte);//add to queue

        NE[P2]++;
        NE[P4]++;

        EPMap.AddEdge(P2, P4, counte);//add edge to the map

        //t[countt].e2=counte;

        counte++;//increase counter edge

    }
    else//edge is old just remove it from front
    {
        //t[countt].e2=idedge2;
        e[idedge2].t2=countt;//add new triangle
    }


    //Edge is no more on front
    e[ide].t2=countt;

    //Edge on front belongs to the current triangle
    //t[countt].e3=ide;


    //Add new triangle
    t[countt].p1=P1;
    t[countt].p2=P4;
    t[countt].p3=P2;


    //for not manifold check

    NT[P4]++;NT[P1]++;NT[P2]++;

    // cout<<P1<<" "<<P2<<" "<<P4<<endl;//display triangles
    countt++;//a new triagnle has been created
}


bool MESHGROWING::BallCenter(Coord3D* cc, double* cr2, Coord3D* tnorm, Coord3D* BC)
//gets the Ball center of a triangle
{


    double offset;



    if (*cr2>R*R) return false;//points outside the torus

    offset=sqrt(R*R-*cr2);

    BC->x=cc->x+offset*tnorm->x;
    BC->y=cc->y+offset*tnorm->y;
    BC->z=cc->z+offset*tnorm->z;



    return true;


}






double MESHGROWING::BCAngle(Coord3D* m, Coord3D* BC, Coord3D* tnorm) {

    double leng;
    Coord3D vBCm, v21, vp;
    double Cos, Sin;

    //build the vector BB-m
    DiffPoints(BC->, m->, vBCm.);
    Normalize(vBCm., leng);

    //compute the cosine of angle btween BC-m and tnorm
    DotProduct(vBCm., tnorm->, Cos);

    //Computing sin
    DiffPoints(p[P2]., p[P1]., v21.);//find the vector p2-p1
    Normalize(v21., leng);
    //calcola la direzione del search point
    CrossProduct(v21., tnorm->, vp.);//vettore perpendicolare al lato

    //normalizza
    Normalize(vp., leng);
    DotProduct(vp., vBCm., Sin);

    if (Sin<0) //sharp corners
    {
        if (Cos>0)Cos=2-Cos;//spigolo tagliente
        else Cos=-2-Cos;//bordo acuto

    }
    return Cos;

}






void MESHGROWING::CircumCenter(Coord3D* p1, Coord3D* p2, Coord3D* p4, Coord3D* n, double *r, Coord3D* cc)

{ //Ritorna il centro ed il quadrato del raggio della sfera SCB
    //calcola inoltre il vettore v21 v41 (non normalizzati) e la normale del triangolo


    Coord3D b, m;
    double leng;
    double rtemp;

//Getting the normal

    DiffPoints(p2->, p1->, v21.);
    DiffPoints(p4->, p1->, v41.);
    CrossProduct(v41., v21., n->);
    Normalize(n->, leng);


    //Assembling the system


    //vettore dei termini noti
    DotProduct(n->, p1->, b.x);//primo elemento vettore termini noti

    MidPoint(p1->, p2->, m.);//salva in m il puntomedio

    DotProduct(v21., m., b.y);//secondo elemento vettore termini noti

    MidPoint(p1->, p4->, m.);//salva in m il puntomedio

    DotProduct(v41., m., b.z);//terzo elemento vettore termini noti


//Solve the system
    FastCramer(&n->x, &v21.x, &v41.x, &b.x, &cc->x);

//CircumRadius
    SquaredDistance(cc->, p1->, *r);
//take the minum distance
//avoids the triangle to self destruct during Dtest
    SquaredDistance(cc->, p2->, rtemp);
    if (rtemp<*r)
    {*r=rtemp;}
    SquaredDistance(cc->, p4->, rtemp);
    if (rtemp<*r)
    {*r=rtemp;}

    #ifdef _DEBUG
            double d1, d2 , d3;
    Distance(cc->, p1->, d1);
    Distance(cc->, p2->, d2);
    Distance(cc->,  p4->, d3);
#endif


}
void MESHGROWING::CircumCenter_Fast(Coord3D* p1, Coord3D* p2, Coord3D* p4, Coord3D* n, double *r, Coord3D* cc)

{ //Ritorna il centro ed il quadrato del raggio della sfera SCB
    //NON calcola inoltre il vettore v21 v41 (non normalizzati) e la normale del triangolo
    // DEVONO ESSERE GIà CALCOLATI

    Coord3D b, m;
    double rtemp;


    //Assembling the system


    //vettore dei termini noti
    DotProduct(n->, p1->, b.x);//primo elemento vettore termini noti

    MidPoint(p1->, p2->, m.);//salva in m il puntomedio

    DotProduct(v21., m., b.y);//secondo elemento vettore termini noti

    MidPoint(p1->, p4->, m.);//salva in m il puntomedio

    DotProduct(v41., m., b.z);//terzo elemento vettore termini noti


//Solve the system
    FastCramer(&n->x, &v21.x, &v41.x, &b.x, &cc->x);

//CircumRadius
    SquaredDistance(cc->, p1->, *r);
//take the minum distance
//avoids the triangle to self destruct during Dtest
    SquaredDistance(cc->, p2->, rtemp);
    if (rtemp<*r)
    {*r=rtemp;}
    SquaredDistance(cc->, p4->, rtemp);
    if (rtemp<*r)
    {*r=rtemp;}

    #ifdef _DEBUG
            double d1, d2 , d3;
    Distance(cc->, p1->, d1);
    Distance(cc->, p2->, d2);
    Distance(cc->,  p4->, d3);
#endif
}







void MESHGROWING::GetSearchPoint(Coord3D* P1, Coord3D* P2, Coord3D* tnorm, Coord3D* sp, double* sr, double kr)//gets the search point for SCBMesher
{



    //memory allocation

    Coord3D cosdir;
    double leng;


    //trova i vettori
    DiffPoints(P2->, P1->, v21.);//find the vector p2-p1


    //calcola la direzione del search point
    CrossProduct(v21., tnorm->, cosdir.);

    //normalizza

    Normalize(cosdir., leng);

    //trova il punto medio dell'edge (ricicla v21)
    MidPoint(P1->, P2->, v21.);

    //sqrt(3)/2=0.866025403784439

    //Calcola il search radius come lunghezza dell'edge
    Distance(P1->, P2->, *sr);

    *sr=(*sr+*sr*kr)*.5;
    //calcola il search point in grandendo un pò sqrt(3)/2=0.866025403784439
    sp->x=v21.x+cosdir.x**sr*0.866025403785;
    sp->y=v21.y+cosdir.y**sr*0.866025403785;
    sp->z=v21.z+cosdir.z**sr*0.866025403785;


    Distance(P1->, sp->, *sr);//ricalcola il search radius per il caso maxr>1
    *sr=*sr*.9999999;

}






void MESHGROWING::GetTriangle_BPA_Torus() {// Gets ball pivoting traignle using the torus method

    //retunrs flag
    // P4>= candidate found
    // P4=-1 no points in range or no good candidate
    // P4=-2 edge is too long
    // P4=-3 sharp triagnle
    // P4=-4 not manifold vertex

    double leng;
    Torus3D Torus;


    P4=-1;//set to minus one in case triagnle can not be found

    P1=e[ide].p1;
    P2=e[ide].p2;
    T1=e[ide].t1;

    Distance(p[P1]., p[P2]., leng);
    if (leng>P.Max_Edge_Length)
    { P4=-2;return;}//edge is too long


    //BUILD THE PIVOTING TORUS
    BuildTorus(&Torus);
    //GET CANDIDATE POINTS (inside the torus

    GetCandidates_Torus(&Torus);

    //For all points in range che check if there is a Delaunay2_5D triangle
    if (SDS.npts==0){P4=-1;return;}//nopoints in range

    //CANDIDATE SELECTION

    SelectCandidate_BPA_Torus(&Torus);

}


void MESHGROWING::GetTriangle_BPA_SR() {// Gets ball pivoting traignle using Search radius

    //retunrs flag
    // P4>= candidate found
    // P4=-1 no points in range
    // P4=-2 edge is too long
    // P4=-3 sharp triagnle
    // P4=-4 not manifold vertex

    double leng, sr;
    Coord3D T1norm, sp;


    P4=-1;//set to minus one in case triagnle can not be found

    P1=e[ide].p1;
    P2=e[ide].p2;
    T1=e[ide].t1;


    Distance(p[P1]., p[P2]., leng);
    if (leng>P.Max_Edge_Length)
    { P4=-2;return;}//edge is too long

    //Get the search point
    TNormTriangle(T1, &T1norm);
    GetSearchPoint(&p[P1], &p[P2], &T1norm, &sp, &sr, .5);

    //Get all neighbors within the search radius
    SDS.SearchRadius( &sp, sr);



    //For all points in range che check if there is a Delaunay2_5D triangle
    if (SDS.npts==0)
    {P4=-1;return;}//nopoints in range


    //CANDIDATE SELECTION

    SelectCandidate_BPA_SR();

}

void MESHGROWING::  SelectCandidate_BPA_SR() {// Gets th e points that forms a BPA triagnle
    int j, Ptemp;

    Coord3D BC;//Ball center
    double Dotp;//we must choose the minimum angle ball center
    Coord3D T1Norm, T2Norm, tempnorm;
    double cr2;//squared circumradius
    Coord3D cc;//circumcenter
    int test;


//initiate the ball center
    //loop trough candidates
    for (j=0;j<SDS.npts;j++)
    { Ptemp=SDS.idStore[j];

      //NOTA! : We can not use the fast circumcenter since we don't have v21 already computed
      //Ball center
      CircumCenter(&p[P1], &p[P2], &p[Ptemp], &tempnorm, &cr2, &cc);
      if (!BallCenter(&cc, &cr2, &tempnorm, &BC)){continue;}//continuee if traignle is too big

      //test if the traignle is BPA
      test=SDS.EmptyBallTest(&BC, R*R*.999999);
      if (test<0 ) {
          P4=Ptemp;
          T2Norm=tempnorm;
          break;//we found one no need of more
      }
    }


    if(P4>=0){
        TNormTriangle(T1, &T1Norm);
        DotProduct( T2Norm., T1Norm., Dotp);
        if (Dotp<P.BPA_TooSharp_SR){{P4=-3;return;}}//sharp triagnle
        if(NE[P4]==NT[P4] && NE[P4]!=0){P4=-4;return;}//not manifold vertex
    }

}

void MESHGROWING:: TNormPoints(Coord3D* p1, Coord3D* p2, Coord3D* p4, Coord3D* tnorm) {
//gets the normal of the trianle formed by the points p1 p2 p4
// it writes data into the globals v21 v14
    //Coord3D v21, v41; uncomment to use local variables
    double leng;
    //computing the normal of the triagnle
    //Get normal of new triangle
    DiffPoints(p2->, p1->, v21.);
    DiffPoints(p4->, p1->, v41.);
    CrossProduct(v21., v41., tnorm->);
    Normalize(tnorm->, leng);//la normale potrebbe essere normalizzata  afine dtest pensaci!!!
}


//The tesselation routine
void MESHGROWING::RunFront_SCB() {

    //Viene adottata la seguente strategia
    //vengono settati alcuni livelli di priorità, prima per il flatness toll poi per il search radius
    //I valori di priorità vengono salvati in e2t. Il flatnnes test viene salvato nelle unità il search radius nelle centinaia

    int i=0;//queue iterator
    int idedge1, idedge2;//id dei presunti nuovi edge
    int plevel;

    Allocate(&FlatToll,P.SCB_Dotp_PriorityLevel);

    //initializing flatoll
    for(i=0;i<P.SCB_Dotp_PriorityLevel;i++)FlatToll[i]=cos((M_PI/2)/P.SCB_Dotp_PriorityLevel*(i+1));



    //ADVANCING FRONT
    cout<<"Running Front... ";
    while (counte<MAXE && Pqueue.Pop(&ide)) {



        //cout<<"i: "<<i<<" ide: "<<ide<<endl; //il bug is verifica per i=9024 ide=9023;

        if (e[ide].t2>=0){continue;}//edge is no more on front skip





        GetTriangle_SCB();

        if (P4<0) {

             if(P4==-1)//no points in search radius
                 {e[ide].t2-=100;
              plevel=P.SCB_Dotp_PriorityLevel-e[ide].t2/100-1;
              Pqueue.Push(ide, plevel);}
             if(P4==-3)//triagnle no flat enough
                 { plevel=-e[ide].t2%100;
                 if(plevel<P.SCB_Dotp_PriorityLevel)
                 {Pqueue.Push(ide,plevel);}//Triangle is not flat enough push the edge in the queue with a different priority}
                  e[ide].t2--;//priority is save into e[i].t2, decrease it
                 }


            continue;
            }//no candidate





        //CHECK EDGES CONNECTIVITY
        idedge1=EPMap.GetEdge(P1, P4);
        if(!CheckEdgeConformity(idedge1, P1)){continue;}//edge not conform
        idedge2=EPMap.GetEdge(P4, P2);
        if(!CheckEdgeConformity(idedge2, P4)){continue;}//edge not conform

        //Se è stato trovato un buon triangolo dobbiamo aggiornare la struttura dati
        UpdateFront(idedge1, idedge2);


    }
    cout<<countt<< " triangles generated"<<endl;

    //Deallocate memory
    Deallocate(&FlatToll);
}

void MESHGROWING::GetTriangle_SCB() {
//Set the value of P4 attempting to find an SCB trianle
       //retunrs flag
    // P4>= candidate found
    // P4=-1 no points in range
    // P4=-2 edge is too long
    // P4=-3 sharp triagnle
    // P4=-4 not manifold vertex

    int j;
    double leng, cr2, Dotp, dist;
    int Ptemp;
    double SearchRadius;
    int plevel;//priority level

    Coord3D sp, T1norm, T2norm;//search point
    double sr;//search radius


    Coord3D cc;//circumcenter

    P1=e[ide].p1;
    P2=e[ide].p2;
    T1=e[ide].t1;


    P4=-1;//setting to default value (triangle can not be found)
    //controllo la lunghezza dell'edge (se troppo lungo salta)
    Distance(p[P1]., p[P2]., leng);
    if (leng>P.Max_Edge_Length)
    { P4=-2;return;}//edge is too long



    //Get the search point
    SearchRadius=pow(2,(double)(-e[ide].t2/100))*.5;//.5*2^i
    //SearchRadius=(-e[ide].t2/100)*.5;
    TNormTriangle(T1, &T1norm);
    GetSearchPoint(&p[P1], &p[P2], &T1norm, &sp, &sr, SearchRadius);

    //Get all neighbors within the search radius
    SDS.SearchRadius( &sp, sr);



    //For all points in range che check if there is a Delaunay2_5D triangle
    if (SDS.npts==0)
    { P4=-1;return;}//nopoints in range

    //RUN SCB TEST

    P4=SDS.idStore[0];
    // if(P4==P1 || P4==P2)//Controlla che il punto non appartenga già al triangolo
    //Questo non dovrebbe accadere ma nelle prime versioni lo controlliamo per sicurezza
    //{continue;}//trova un punto diverso

    //Get normal of new triangle
    DiffPoints(p[P2]., p[P1]., v21.);
    DiffPoints(p[P4]., p[P1]., v41.);
    CrossProduct(v41., v21., T2norm.);
    //we normalize in the end

    //Scan all points to check if we can find an SCB triagnle
    for (j=1;j<SDS.npts;j++) //esegui il dtest solo se nel search radius ci sono più di punti
    {

        // trova il circocentro del triangolo per effettuare il test SCB
        CircumCenter_Fast(&p[P1], &p[P2], &p[P4], &T2norm, &cr2, &cc);

        //get a point form the search region
        Ptemp=SDS.idStore[j];//punto su cui fare il dtest

        SquaredDistance(cc., p[Ptemp]., dist);//calcola la distanza dal punto sotto test

        if (dist<cr2)//se il punto è dentro il test fallisce *0.9999
        { P4=Ptemp;//ricomincia da questo punto
          //non delaunay recompute CC with new point, but first we need to get the noirmal of the new triangle

          //Get normal of new triangle
          //DiffPoints(p[P2]., p[P1]., v21.); //we already have it
          DiffPoints(p[P4]., p[P1]., v41.);
          CrossProduct(v41., v21., T2norm.);
          // Normalize(t[countt].tnorm., cr);//la normale potrebbe essere normalizzata  afine dtest pensaci!!!

        }
    }

    //Then check if it breaks some manifolds

    if(NE[P4]==NT[P4] && NE[P4]!=0)
    {P4=-4;return;}//not manifold vertex

    Normalize(T2norm., leng);//normalize the normal

    //Flatness test
    DotProduct(T2norm., T1norm., Dotp);//coseno dell'angolo fra 2 triangoli
    plevel=-e[ide].t2%100-1;
    if(Dotp<=FlatToll[plevel] )
        {P4=-3;return;}








}
void MESHGROWING::PreP_NN_MeanDist(double*meandist) {
    //trova la distanza media tra i punti calcolando il NN Graph
    // Elimina i punti doppi per cui dist=0;

    int i, idc;

    double dist;


    //Memory allocation


    *meandist=0;//init

    //get the NN for each point
    for(i=0;i<N;i++) {
        SDS.SearchClosestExclusive(&p[i], &idc, &dist, i);

        if (dist>0)
        {*meandist+=dist;}
        else //duplicate point
        {SDS.RemovePoint(i);Ndel++;}//remove duplicate point

    }

    //Getting the mean value of the distances
    *meandist/=(N-Ndel);//mean distance
}


void MESHGROWING::PreP_NN_Filter(double cutdist) {  //filtro che rimuove tutti i punti nel raggio di cutdist
    // at the end no points will be closer than cutdist

    int i, j;
    bool* removed=NULL;AllocateAndInit(&removed, N, false);//flag if points has been removed

    //remove bad points
    for(i=0;i<N;i++) {

        if (!removed[i]) {
            SDS.SearchRadiusExclusive(&p[i], cutdist, i);
            for(j=0;j<SDS.npts;j++)
            {SDS.RemovePoint(SDS.idStore[j]);Ndel++;
             removed[SDS.idStore[j]]=true;
            }//remove points in range


        }
    }

    cout<<"Deleted "<<Ndel<<" points"<<endl;

    if ((N-Ndel)<3)Error("Too many points deleted by the filter, reduce filtering stregth");

}



void MESHGROWING::PostP_FillTriplets() { //funtion to fill empty triplets
    int i;
    int P1, P2, P4;
    int e14, e24;
	int startt=countt;

	cout<<"Filling triplets: ";

    int* Map=NULL; AllocateAndInit(&Map, N, -1);//initializza to -1
    int* StoredEdge=NULL; Allocate(&StoredEdge, N); //initializza to -1

    bool found=false;

    //loop trough all the edges

    for (i=0;i<counte;i++) {
        if(e[i].t2<0 && e[i].t1>=0)//boundary edge not deleted
        {
            P1=e[i].p1;P2=e[i].p2;
            if (Map[P1]>=0)//il punto era già allocato nella mappa
            {	//edge p1-p4 was in the map
                P4=Map[P1];//trova il terzo punto chiamato P4 secondo le precedenti convenzioni

                e24=EPMap.GetEdge( P2, P4);
                if (e24>=0 && e[e24].t2<0 && e[e24].t1>=0)//se l'edge esiste ed è boundary abbiamo trovato un triangolo
                {
                    e14=StoredEdge[P1];
                    if (e[e14].p1==P1 || e[e24].p1==P4) {//cout<<P1+1<<" "<< P2+1<<" "<< P4+1<<" "<<endl;
                        continue;}//degenere}
                    found=true;
                }
            }
            else if(Map[P2]>=0) {	//edge p2-p4 was in the map
                P4=Map[P2];//trova il terzo punto chiamato P4 secondo le precedenti convenzioni

                e14=EPMap.GetEdge( P1, P4);
                if (e14>=0 && e[e14].t2<0 && e[e14].t1>=0)//se l'edge esiste, è boundary e non è cancellato abbiamo trovato un triangolo
                {
                    e24=StoredEdge[P2];
                    //Check the orientation
                    if (e[e14].p1==P1 || e[e24].p1==P4) {//cout<<P1+1<<" "<< P2+1<<" "<< P4+1<<" "<<endl;
                        continue;}//degenere}

                    found=true;}
            }
            else//aggiungi i punti alla mappa
            {
                Map[P1]=P2;
                Map[P2]=P1;
                StoredEdge[P1]=StoredEdge[P2]=i;
                continue;//go to next edge//no triangle found}
            }


            // in relatà il found non serve perhè nell'else sopra c'è continue
            if (found)	//formto da P1 P2 P4 e con gli edge i,idedge
            {
                NT[P1]++;NT[P2]++; NT[P4]++; //
                t[countt].p1=P1;
                t[countt].p2=P4;
                t[countt].p3=P2;



                //update e2t
                e[i].t2=countt;
                e[e14].t2=countt;
                e[e24].t2=countt;






                countt++;//MANCA IL CONTROLLO DI OVERFLOW
                found=false;//reset to false;

                //ResetMap
                Map[P1]=Map[P2]=Map[P4]=-1;



            }


        }
    }



	cout<<countt-startt<<"Triangles "<<endl;

    Deallocate(&Map);
    Deallocate(&StoredEdge);
}







void MESHGROWING::ImportPoints_Pointers(double* pointer, int inputN) {
    //copy the pointer to the internal data structure and allocate memory
    N=inputN;


    p=(Coord3D*)pointer;

    //We ahve N now we can allocate memory
    Memory_Allocate();//allocate all the memory in the class constructor
}


inline void MESHGROWING::TNormTriangle(int idt, Coord3D* tnorm) {//gets the normal tnorm of the triangle idt
    Coord3D v21, v31;

    double leng;
    int p1, p2, p3;
    p1=t[idt].p1;p2=t[idt].p2;p3=t[idt].p3;
    //computing the normal of the triagnle
    //Get normal of new triangle
    DiffPoints(p[p2]., p[p1]., v21.);
    DiffPoints(p[p3]., p[p1]., v31.);
    CrossProduct(v21., v31., tnorm->);
    Normalize(tnorm->, leng);//la normale potrebbe essere normalizzata  afine dtest pensaci!!!

}

void MESHGROWING::PostP_HealNMV()//delete triagnle related to not manifold vertxes
{//if a manifold vertex has boundary edges, delete recursively all triangles with boundary edges linked with the NMV
    ////if a manifold vertex has no boundary edge, delete all triangles related to the NMV and remove it from SDS
	// then delete all traingles related to overboundary points (more than 2 boundary edges)
    int i, nNMV;
    int T1;
    char* manifold=NULL;Allocate(&manifold, N);

    AllocateAndInit(&nbe, N,(char)0);//number of boundary edges for each point


    //counting the number of boundary edges for each point
for (i=0;i<counte;i++)
   {
     if (IsBound(i))//boundary edges not deleted
         {nbe[e[i].p1]++;nbe[e[i].p2]++;}
}



 nNMV=PostP_FindNMV(manifold);



  cout<<"Found "<<nNMV<<" starting not manifold vertxes"<<endl;

/*
#ifdef _DEBUG
	cout<<"Strong Not manifold points"<<endl;
for (i=0;i<N;i++) {
    if( manifold[i]==0){cout<<i<<endl; }
}
#endif
*/


//delete closed NMV
for (i=0;i<N;i++)
   {
        if(manifold[i]==0 && nbe[i]==0)
            {PostP_Delete_NMV_Triangles_All(i);}
}


//loop trough all edges and delete recursively boundary triangle related to NMV

for (i=0;i<counte;i++) {

    if (IsBound(i))//boundary edges not deleted
    {T1=e[i].t1;

     //if the triangle has one not manifold vetex
     if (manifold[t[T1].p1]==0){PostP_Delete_NMV_BoundTriangle_idp(T1, t[T1].p1); }
     if(manifold[t[T1].p2]==0){PostP_Delete_NMV_BoundTriangle_idp(T1, t[T1].p2); }
     if(manifold[t[T1].p3]==0){PostP_Delete_NMV_BoundTriangle_idp(T1, t[T1].p3);}

    }

}



//Ora dobbiamo eliminare i NMV che hanno più di 2 boundary edges

for (i=0;i<counte;i++)
   {
     LookForNMV_OverBound_Triangle(i)
}


#ifdef _DEBUG
//check the nbe status: if it gets odds number  it is corrupted.
for (i=0;i<N;i++)
   {
	   if ((nbe[i]%2)==1){cout<<"point"<< i<< " nbe corrupted"<<endl;}
	   //if (nbe[i]>2){cout<<"point"<< i<< " is still not manifold"<<endl; }
	     if (NT[i]>0 && manifold[i]<0){cout<<"point"<< i<< " is wrong marked"<<endl; }//uncheked points can not have no triagnle
}
#endif


//deallocate
Deallocate(&manifold);
Deallocate(&nbe);
}

void MESHGROWING::PostP_Delete_NMV_Triangles_All(int idp)//Deletes all triagles related with NMV idp
{
  int idt;
  int e1, e2, e3;

  for (idt=0;idt<countt;idt++)
   {
     if (t[idt].p1==idp || t[idt].p2==idp || t[idt].p3==idp)//boundary edges not deleted
         {

                          if (t[idt].p1<0)continue;//triangle already deleted

                            //Update triangle edges

                            e1=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p2,idt);
                            e2=PostP_UpdateEdge_DeletedTriangle(t[idt].p2, t[idt].p3,idt);
                            e3=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p3,idt);

                            t[idt].p1=-1;Ntdel++;//delete triangle
         }

     }


}

void MESHGROWING::PostP_Delete_NMV_OverBound_Triangle(int idt)//Deletes the over boundary triangles idt
{
      int e1, e2, e3;


                            if (t[idt].p1<0)return;//triangle already deleted


                            //Update triangle edges
                            e1=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p2,idt);
                            e2=PostP_UpdateEdge_DeletedTriangle(t[idt].p2, t[idt].p3,idt);
                            e3=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p3,idt);

                            t[idt].p1=-1;Ntdel++;//delete triangle

                        //call recursively
                            LookForNMV_OverBound_Triangle(e1)
                            LookForNMV_OverBound_Triangle(e2)
                            LookForNMV_OverBound_Triangle(e3)

}
void MESHGROWING::PostP_Delete_NMV_BoundTriangle_idp(int idt, int idp)//Deletes triangle idt related with NMV idp
{



//function body
                            int e1, e2, e3, T2;

                            if (t[idt].p1<0)return;//triangle already deleted



                          //Update triangle edges
                            e1=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p2,idt);
                            e2=PostP_UpdateEdge_DeletedTriangle(t[idt].p2, t[idt].p3,idt);
                            e3=PostP_UpdateEdge_DeletedTriangle(t[idt].p1, t[idt].p3,idt);


                            t[idt].p1=-1;Ntdel++;//delete triangle


#ifdef _DEBUG
cout<<"Triangle: "<<idt<<" is not manifold due to point: "<<idp<<endl;

#endif

//call recursively trough neighbour trianlges
LookForNMV_BoundTriangle_idp(e1)
LookForNMV_BoundTriangle_idp(e2)
LookForNMV_BoundTriangle_idp(e3)

}

int MESHGROWING::PostP_UpdateEdge_DeletedTriangle(int p1,int p2,int idt)//update the edge status before a triagnle idt that contains the edge is deleted
    {
        int ide,temp;
        bool bound;
        ide=EPMap.GetEdge(p1,p2);
         IsBound(ide)? bound=true : bound=false;
       Reset_e2t(ide)//update e2t connectivy
      if(!bound){nbe[p1]++;nbe[p2]++;}
      else {nbe[p1]--;nbe[p2]--;}//boundar edge no more exists

      return ide;
    }


void MESHGROWING::PostP_Restore_t()//restore the t array due to deleted triagnles
{
    int i, c;

    c=0;
    for (i=0;i<countt;i++) {


        if (t[i].p1>=0) {
			t[c]=t[i];
		     c++;
		}

    }

    countt=c;//reset counters
}

int MESHGROWING::PostP_FindNMV(char* manifold)//Finds  not manifold vertices
{//richiede che nbe sia calcolato

//Strategia:
//
// esamina e2t connctivity fino a che non si hiude il cerchio attorno ad un punto. Il cerchio viene identificato da un
// Back ed un Front
// Per gestire il controllo delgi edge viene usata una coda a doppia entrata
// Gli edge hanno un orientazione, in e2t T1 va a T2 in senso antiorario

//char*manifold -1=untested 0=not manifold 1=manifold

//WARNING: global(queue) will be reset

    int  i, P1, T1, T2, iq, last, lastold, nNMV;
    bool checked;//counts how many edges have correspondence
    int* Back=NULL;AllocateAndInit(&Back, N, -2);//init to -2 per non confonderlo con l'edge di boundary marcato con -1
    int* Front=NULL;AllocateAndInit(&Front, N, -2);
    int16_t* Count=NULL;AllocateAndInit(&Count, N, int16_t(0));


    nNMV=0;

//Presume all points are unchecked
    Alls(manifold, N, (char)-1);

//Flag overboundary points as not manifold
	for (i=0;i<N;i++)
	{if (nbe[i]>2){manifold[i]=0;nNMV++;}}

//push all edges inside the queue
    for(i=0;i<counte;i++) {
        queue[i+1]=i;//i+1 the array position 0 is occupied by the sentinel
        if(e[i].t2<0){e[i].t2=-1;}//to equal all boundary triangles
    }


    last=0;
    lastold=counte;
    while (1)//loop until we have analyzed all edges
    {  iq=counte;//initiate front queue iterator
       queue[last]=-1;//setting sentinel
       last=counte;//initiate back front iterator
       while(1)//loop until sentinel
       {
           ide=queue[iq];if (ide<0)break;//break when meet sentinel



           if(ide<counte)
           {i=ide;P1=e[i].p1;T1=e[i].t1;T2=e[i].t2;}//counter clockwise  orientation
           else
           {i=ide-counte; P1=e[i].p2;T2=e[i].t1;T1=e[i].t2;}//clockwise orientation

            if(e[i].t1<0)
			   { iq--;continue;}//edge deleted


#ifdef _DEBUG
if(P1==0)
{cout<<"Point under debug"<<endl;}
#endif

//trying to close the loop

checked=false;//correspondence not found yet;
if (manifold[P1]>=0){  checked=true;}//correspondence found points already falgged as manifold or not manifold
else if(Back[P1]==-2)//Point circle is still empty
{
    checked=true;//correspondence found
    Back[P1]=T1;Front[P1]=T2;//copy
    if(T1>=0){Count[P1]++;}//one more triangle
    if(T2>=0){Count[P1]++;}//one more triangle
}
else if( Front[P1]==T1 )//advance front of P1 in counterclockwise order
{
    checked=true;//correspondence found
    Front[P1]=T2;
    if(T2>=0 ){Count[P1]++;}//one more triangle
    if (Front[P1]==Back[P1])//loop is closed check manifold
    {
        if(T2>=0){Count[P1]--;}//remove joining triagnle
        if(NT[P1]==Count[P1] || NT[P1]==1)//the number of traignles to close the loop must equals the total number of neighbours triangles
        { manifold[P1]=1;}
        else
        { manifold[P1]=0;nNMV++;}

    }

}
else if(Back[P1]==T2)//advance back of P1 in counterclockwise order
{checked=true;//correspondence found
 Back[P1]=T1;
 if(T1>=0){Count[P1]++;}//one more triangle
 if (Front[P1]==Back[P1])//loop is closed check manifold
 {
     if(T1>=0){Count[P1]--;}//remove joining triagnle
     if(NT[P1]==Count[P1] || NT[P1]==1)//the number of traignles to close the loop must equals the total number of neighbours triangles
     { manifold[P1]=1;}
     else
     { manifold[P1]=0;nNMV++;}

 }
}



if(!checked)//no corresponde found, re-push the edge in the back queue
{queue[last]=ide;//push in the back queue
 last--;
}
else if(ide<counte)//if we have just analyzed the edge with counterclockkwise orientation, analyze it with clockwise
{queue[iq]=ide+counte;
 iq++;//to compensate the minus 1
}

iq--;//decrease queue iterator
       }

       if (lastold==last)
       {break;}
       lastold=last;
    }


//Deallocate memory
    Deallocate(&Back);
    Deallocate(&Front);
    Deallocate(&Count);
    Pqueue.Reset();//clean the queue before exting

    return nNMV;
}



void MESHGROWING::PostP_HoleFiller()
//Questa funzione può essere ottimizzata nel seguente modo:
//-Flaggare i punti di boundary e considerare solo quelli come candidati
//-In alternativa costruire un secondo GLTree solo con ipunti di boundary

{
    //FUNCTIONS DELARATIONS



    //MEMORY ALLOCATION




    int i=0;//queue iterator
    int ide;

    int  j;//counter



    Coord3D sp;//Search point
    double sr;//Searc radius


    int idedge1, idedge2;//id dei presunti nuovi edge

//    Coord3D v21;//vettori per le operazioni vettoriali e midpoint per il sistema

    //point id map

//          P2
//        / | \
// 	T1  P3  |  P4(newpoint) T2
//        \ | /
//		    P1

    //int flag;
    double dist;//double per operazioni varie
    Coord3D T1norm;
    bool found;//per vedere se un triangolo è stato trovato


    bool*bound=NULL;AllocateAndInit(&bound, N, true);

    //flagginf boundary points
    for (i=0;i<counte;i++) {
        if (e[i].t2<0 && e[i].t1>=0)//se l'edge è di boundary e non è stato cancellato
        {bound[e[i].p1]=bound[e[i].p2]=true;}//flag as bnoundary point
    }


    //ADVANCING FRONT
    i=0;
    while (counte<MAXE && Pqueue.Pop(&ide)) {

        if (e[ide].t2>=0 || e[ide].t1<0){continue;}//edge is no more on front or has been deleted

        P1=e[ide].p1;
        P2=e[ide].p2;

        //controllo la lunghezza dell'edge (se troppo lungo salta)
        Distance(p[P1]., p[P2]., dist);
        if (dist>P.Max_Edge_Length)
        {continue;}//edge is too long


        T1=e[ide].t1;

        //Get the search point
        TNormTriangle(T1, &T1norm);
        GetSearchPoint(&p[P1], &p[P2], &T1norm, &sp, &sr, -e[ide].t2);

        //Get all neighbors within the search radius


        SDS.SearchRadius(&sp, sr);

        //L'HOle filler fa una sola iterazione non c'è bisogno di incrementare il search radius
        //if (SDS.npts==0)
        //{ e[ide].t2--;//flag as boundary
        // continue;}//nopoints in range





        //Analizza tutti i punti per vedere se riusciamo a trovare un  edge noto e confome
        found=false;
        for (j=0;j<SDS.npts;j++) {
            P4=SDS.idStore[j];

            if (!bound[P4]){continue;}//point is not boundary

            idedge1=EPMap.GetEdge(P1, P4);


            if (idedge1>=0)//solo se l'edge esiste già
            {if(e[idedge1].t2>=0 || e[idedge1].t1<0 ){continue;}//posto occupato o edge cancellato
             else {

                 //controlla l'orientamento
                 if (e[idedge1].p1==P1){continue;}//orientazione non conforme

             }
            }

            idedge2=EPMap.GetEdge(P4, P2);



            if (idedge2>=0)//solo se l'edge esiste già
            {
                if(e[idedge2].t2>=0 || e[idedge2].t1<0 ){continue;}//posto occupato o edge cancellato
                else{

                    //controlla l'orientamento
                    if (e[idedge2].p1==P4){continue;}//orientazione non conforme

                }

            }

            if(idedge1>=0 || idedge2>=0)//almeno uno dei 2 esiste già
            {found=true;break;}//trovato un triangolo conforme esci

        }


        if (!found){continue;}




        //Se è stato trovato un buon triangolo dobbiamo aggiornare la struttura dati

        //Inizializza la "libertà" dei punti a falso
        // FreePoint[P1]=FreePoint[P2]=FreePoint[P4]=false;

        //Check the first edge
        if  (idedge1<0)// %if edge1 is new
        {
            e[counte].p1=P1;//add new edge
            e[counte].p2=P4;//add new edge
            e[counte].t1=countt;//add new triangle
            Pqueue.PushHigh(counte);//add to queue
            // FreePoint[P1]=FreePoint[P4]=true;//liberi

            EPMap.AddEdge(P1, P4, counte);//add edge to the map



            counte++;//increase counter edge

        }
        else//edge is old just remove it from front
        {

            e[idedge1].t2=countt;//add new triangle
        }


        //Check the second edge
        if  (idedge2<0)// %if edge1 is new
        {
            e[counte].p1=P4;//add new edge
            e[counte].p2=P2;//add new edge
            e[counte].t1=countt;//add new triangle
            Pqueue.PushHigh(counte);//add to queue
            //FreePoint[P2]=FreePoint[P4]=true;//liberi

            EPMap.AddEdge(P2, P4, counte);//add edge to the map



            counte++;//increase counter edge

        }
        else//edge is old just remove it from front
        {
            //t[countt].e2=idedge2;
            e[idedge2].t2=countt;//add new triangle
        }


        //Edge is no more on front
        e[ide].t2=countt;

        //Edge on front belongs to the current triangle


        //Add new triangle
        t[countt].p1=P2;
        t[countt].p2=P1;
        t[countt].p3=P4;


        // cout<<P1<<" "<<P2<<" "<<P4<<endl;//display triangles
        countt++;//a new triagnle has been created

//         Debug
// Vedi se riesce a trovare i veritici non manifold



    }

    //deallocate
    Deallocate(&bound);

}


void MESHGROWING::PostP_FillQueue(int value) {//Fills the queue with all boundary edges
    int i;
    for(i=0;i<counte;i++)
    {if (e[i].t2<0){
         e[i].t2=value;
         Pqueue.PushHigh(i);
     }
    }//add to queue

}




