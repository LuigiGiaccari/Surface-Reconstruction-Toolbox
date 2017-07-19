
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "timers/hr_time.h"
#include "util/ExceptionHandler.h"


using namespace std;

void ValidateTestFile(char* filename)//checks if all tools are present
{
 //Dichiarazioni
    FILE * pFile ;

    pFile = fopen(filename, "rb");//Open the input binary file
	if (pFile  == NULL)
	{ cout<<filename<<" not found"<<endl;
	Error("Unable to find test file");}
	fclose(pFile);


}

void ValidateAlgo(char* command,char* algoname)//checks if all tools are present
{
 //Dichiarazioni
	int flag;

	cout<<"Validating " <<algoname<<endl;
    flag=system(command);
	if(flag!=0){Error("Error during test");}
	cout<<algoname<<" Validated!"<<endl;
}




int main()

{
 PAUSE_BEFORE_EXIT=true;//artificially avoids exit without pause
 cout<<"TEST STARTED"<<endl;

 ValidateTestFile("bunny.dat");
  ValidateTestFile("bunny.cgo");
  ValidateTestFile("kili_490k.cgo");

#ifdef WIN32
ValidateTestFile("scbmesher.exe ");
ValidateTestFile("ballpivoting.exe");
 ValidateTestFile("robustcrust.exe");
  ValidateTestFile("quickterrain.exe");
 ValidateAlgo("scbmesher.exe -in bunny.cgo -out bunny_SCB_cgo.stl -pa 0", "SCBMesher");
 ValidateAlgo("scbmesher.exe -in bunny.dat -out bunny_SCB_dat.stl -pa 0", "SCBMesher");
 ValidateAlgo("ballpivoting.exe -in bunny.cgo -out bunny_BPA_cgo.stl -r 0 -pa 0", "BallPivoting");
 ValidateAlgo("ballpivoting.exe -in bunny.dat -out bunny_BPA_dat.stl -r 0 -pa 0", "BallPivoting");
 ValidateAlgo("robustcrust.exe -in bunny.cgo -out bunny_RC_cgo.stl -pa 0", "RobustCrust");
 ValidateAlgo("robustcrust.exe -in bunny.dat -out bunny_RC_dat.stl -pa 0", "RobustCrust");
   ValidateAlgo("quickterrain.exe -in kili_490k.cgo -out kili_cgo.stl -pa 0", "QuickTerrain");
#else
ValidateTestFile("./scbmesher");
ValidateTestFile("./ballpivoting");
 ValidateTestFile("./robustcrust");
  ValidateTestFile("./quickterrain");
   ValidateAlgo("./robustcrust -in bunny.cgo -out bunny_RC_cgo.stl -pa 0", "RobustCrust");
 ValidateAlgo("./robustcrust -in bunny.dat -out bunny_RC_dat.stl -pa 0", "RobustCrust");
 ValidateAlgo("./scbmesher -in bunny.cgo -out bunny_SCB_cgo.stl -pa 0", "SCBMesher");
 ValidateAlgo("./scbmesher -in bunny.dat -out bunny_SCB_dat.stl -pa 0", "SCBMesher");
 ValidateAlgo("./ballpivoting -in bunny.cgo -out bunny_BPA_cgo.stl -r 0 -pa 0", "BallPivoting");
 ValidateAlgo("./ballpivoting -in bunny.dat -out bunny_BPA_dat.stl -r 0 -pa 0", "BallPivoting");

  ValidateAlgo("./quickterrain -in kili_490k.cgo -out kili_cgo.stl -pa 0", "QuickTerrain");
#endif

 cout<<"TEST ENDED"<<endl;

 cout<<endl<<"UNMESHABLE IS NOTHING!"<<endl<<endl;

 ExitProgram(0);
}
