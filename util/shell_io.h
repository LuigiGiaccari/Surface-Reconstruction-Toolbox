#ifndef _shell_io_h
#define _shell_io_h

//Header file that manages input outputs for shell commands


#include <time.h>
#include "util/filemanager.h"
#include <fstream>
//GLOBALS
    char inputfile[256]="NONE";//name of the input file
    char   outputfile[256]="NONE";//name of output file


    double R=0;
    const char NOFILE[256]="NONE";
 FILE_MANAGER FileManager;//group of routines to handle files








void PrintHeader(char*algoname,char*version,char*author,char* user)
    {


 char* license=" The Software provided has to be considered ""as is"" and it is without any kind of warranty. The authors deny any kind of warranty concerning the Software  as well as any kind of responsibility for problems and damages which may be caused by the use of the Software  itself.";
//Getting current time
time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

cout<<"INFOS"<<endl;
 cout<<algoname<<" "<<version<<endl;
 cout<<"Author: "<<author<<endl;
 cout<<"User: "<<user<<endl;
  cout<<"Build Time: "<< asctime (timeinfo)<<endl<<endl;;

  cout<<"LICENSE TERMS"<<endl;
  cout<<license<<endl;



    }

void AskForInputs(int mode=0)//asks for command line when is not provided
{
//ScbMesher and RC are fine with mode=0 (default value)
//BPA require to have mode=1 to ask the radius value

char inputstring[256]="__________";//a file that can not be found (I hope so:-))
while(1){
cout<<"Please insert the Input File path"<<endl;
cin>>inputstring;
if (!ifstream(inputstring))
{Warning("file not found!");}
else
{strcpy(inputfile,inputstring);break;}
}
cout<<"Please insert the Output File path"<<endl;
cin>>outputfile;
if (mode==1)
{cout<<"Please insert the Ball Radius (set zero for auto-mode) "<<endl;
cin>>R;}

}

void ReadInputs(int argc,char *argv[])
{

    //legge gli input dell'exe
	// ritorna:
	// 0 se tutto ok
	// 1 per un errore negli input
	int temp;
	int i;



	//loop trough all the inputs parameters with step=2
	//jump the first input which is the path
	for(i=1;i<argc;i=i+2){


    //INPUTFILE
    if (strcmp ("-in",argv[i])==0)//user selct input file
	{ strcpy (inputfile,argv[i+1]);
	cout<<"-Input File= "<<inputfile<<endl;}

	   //OUTPUTFILE
	else   if (strcmp ("-out",argv[i])==0)//user selct output file
	{ strcpy (outputfile,argv[i+1]);
	cout<<"-Output File= "<<outputfile<<endl;}

	 //PAUSE at end execution
	else	if (strcmp ("-pa",argv[i])==0)//user selct input file
	{temp=atoi(argv[i+1]);if(temp==0)PAUSE_BEFORE_EXIT=false ;}//set pause =false if value ==0

	//RADIUS for BPA
    else if (strcmp ("-r",argv[i])==0)//user select Radius only for BPA
	{R=atof(argv[i+1]);
	if (R<0) Error("Negative Radius Value");
	cout<<"-Ball radius= "<<R<<endl;
	}

	else{cout<<"Unknow input Type: '"<<argv[i]<<"' for input number: "<< (i+1)/2<<endl;
	Error("Wrong command line");
	}

	}




 if	(strcmp (NOFILE,inputfile)==0){Error("No Input File");}
 if	(strcmp (NOFILE,outputfile)==0){Error("No Output File");}

    cout<<endl;//crea spazion per le operazioni successive


}





#endif
