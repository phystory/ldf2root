// Universial Root Tree builder for Oak Ridge Data(LDF) at ORNL
// Based on Jeff's code for Yale data and Milan's prototype code 

// This reads PAC file to get ID information of the data, and build trees based on the file.

// This can deal with multiple runfiles in one ldf.
// See the ints start and stop defined below.

// To use provide the name of the ldf file, the name of the output root file and optionally the block range desired.

//Latest Update: February 11th, 2015, version 5.1 for including the case of no cid in the pac/ldf file.
//               Thanks to Andrew Ratkiewicz.
//Latest Update: February 10th, 2015, version 5.0 for generating a readroot.C code to provide a simple root analysis code to users.
//               The code generates TH1D histograms and fill the adc channel and value.
//               Thanks to Andrew Ratkiewicz.
//Latest Update: March 31st, 2014, version 4.0 for changing structure of branches due to the TTree.Fill() update all branches for every events.
//               Now the program reads PAC file only for finding max pacid.
//               Thanks to Justin Browne.
//Latest Update: March 14th, 2014, version 3.0 for changing size of some arrays.
//Latest Update: November 22th, 2011, version 2.0 for adding run number and run title.

#include <iostream>
#include <string>
#include <fstream>		// File reading stuff
#include <iomanip>		// Needed for hex values
#include <cstdlib>		// Contains exit() function
#include <unistd.h>
#include <limits.h>
#include "TFile.h"		// Root headers
#include "TTree.h"		// Ditto
#include "TRandom.h"
#include "TH1.h"
using namespace std;

// To compile you must include ROOT libraries:
//LIBS=-L/usr/lib/root -lCore  -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRIO -lRint
//GCC=g++
//ldf2root.o :  ldf2root.cpp
//    $(GCC) -g -O -Wall -fPIC -pthread -I/usr/include/root -c ldf2root.cpp
//ldf2root :  ldf2root.o
//    $(GCC) -O ldf2root.o $(LIBS)  -o ldf2root


// The binary will expect three arguments -
// pacFilename.pac, InputFilename.ldf, and OutputFilename.root.
// Files do not require byteswapping

// argc (number of arguments) should be 4 or 6 (binary name, 1st input, 2nd input, 3rd input, [startblock, endblock])
// argv holds arguments - argv[1] is input, argv[2] is input, argv[3] is output

// To give a fancy progressbar
void progressbar(int percent);

// To generate a root code reading the root file
int makerootcode(string rootname, string treename, int maxpacid);

int main(int argc, const char* argv[]){	

	int filelength,nblock=0,nevt=0,nword=0;		// Counters for data blocks, events and words	
	int i=0,block=0,tevt=0;				// Current block, total events
	//char buf[32776];				// Buffer to hold each block
	unsigned short word[16388];		// 16388 2-byte words
	int adc,value;					// Tempory values before wrting to disk
	long curdata;					//
	int header=65552;				// Number of bytes in the header
	char buffer[300];
	string bstr;
	int hbegin=0,hend=0,runno=0;
	char filename[30],head[5],rundate[17],runtitle[81];
	int blocksize=32776;				// Number of bytes in each block
	string pacname,ldfname,rootname;
	
	// There should be either 3 or 6 arguments passed
	// If only 3 arguments then all blocks are assumed to be part of the same run
	//if (argc<=2 && argc==5){
	if (argc<=2){
		cout << "Missing arguments" << endl;
		cout << "Usage: " << argv[0] << " PacFile EventFile [RootFileName StartBlock EndBlock]" << endl;
		exit(1);
	}
    else if (argc==3){
		cout << "RootFileName will be copied from EventFileName" << endl;
        pacname = argv[1];
        ldfname = argv[2];
		rootname = ldfname.substr(0,ldfname.length()-3)+"root";
		cout << "root Filename is " << rootname << endl;
    }
    else{
        pacname = argv[1];
        ldfname = argv[2];
        rootname = argv[3];
    }

	bool partial=false;
	int start=0,stop=0;
	if (argc==5){
		partial=true;
		}

	// Read in from the first argument on the command line
	string treename,adcname;
	ifstream fs(ldfname.c_str(),ios::binary);

    // Read PAC file and get treename, branchname, and leafname
    ifstream pac(pacname.c_str(),ios::in);
    string line;
    char cline[50];
    char adcs[50][10][30];
    int result,maxcmdline;
    int chbegin=0,chend=0,idbegin=0,idstep=1,pacid=0,maxpacid=0;
    int iscid=0;
    if(pac.is_open()){
      i=0;
      while(getline(pac,line)){
        if(line[0]!=';'){
          if(line.empty() != 1){
            strcpy(cline, line.c_str());
            result=sscanf(cline, "%s %s %s id%s %s\n",&adcs[i][0][0],&adcs[i][1][0],
                                       &adcs[i][2][0],&adcs[i][3][0],&adcs[i][4][0]);
            i++;
          }
        }
      }
      maxcmdline = i;
      pac.close();
    }
    else {
      cout << "Unable to open " << pacname << endl;
      return 0;
    }

	// Create a root file with the name of the second argument
	TFile *ff = new TFile(rootname.c_str(),"UPDATE");
	// Give names and descriptions
	treename =ldfname.substr(0,ldfname.length()-4);
	TTree *tt = new TTree(treename.c_str(),"Root Data generated by ldf2root");
	// Find max pacid from PAC file
	for(i=0;i<maxcmdline;i++){
      line = adcs[i][0];
      if(line == "$vme"){
	    adcname = adcs[i][1];
        line = adcs[i][2];
        chbegin = atoi((line.substr(1,line.find("-")-1)).c_str());
        chend = atoi((line.substr(line.find("-")+1,line.find("-")-1)).c_str());
        line = adcs[i][3];
        idbegin = atoi((line.substr(0,line.find(","))).c_str());
        idstep = atoi((line.substr(line.find(",")+1,line.find(","))).c_str());
        pacid = idbegin + (chbegin-chbegin)*idstep;
		if(maxpacid<pacid) maxpacid = pacid;
        pacid = idbegin + (chend-chbegin)*idstep;
		if(maxpacid<pacid) maxpacid = pacid;
      }else if(line == "$cid"){
	    adcname = "cid";
            iscid = 1;
        pacid = atoi(adcs[i][1]);
		if(maxpacid<pacid) maxpacid = pacid;
        pacid = atoi(adcs[i][2]);
		if(maxpacid<pacid) maxpacid = pacid;
      }
	} 

	// Define three branches for multiplicity, pacid and value
	int multiplicity,bpacid[maxpacid+1],bvalue[maxpacid+1];					// Arrays for pacid and value
	tt->Branch("multiplicity",&multiplicity,"multiplicity/I");
	tt->Branch("pacid",&bpacid,"bpacid[multiplicity]/I");
	tt->Branch("value",&bvalue,"bvalue[multiplicity]/I");
	
        int counter[maxpacid+1];
        for(int i=0;i<maxpacid+1;i++) counter[i]=0;

	//get size of file and number of blocks
	// Each block is 'blocksize' bytes
	// The preamble is 'header' bytes
	fs.seekg(0,ios::end);
	filelength = fs.tellg();
	nblock = (filelength-header)/(blocksize);
	cout << "File Length: "<< filelength << " bytes; " << nblock << " blocks" << endl;
	
	//return get pointer to beginning of file
  	fs.seekg(0,ios::beg);
	
	//make four branches for header information
	sprintf(filename,"%s",ldfname.c_str());
	//move pointer past first 32776 bytes (the begining point of header)
	fs.seekg(32776,ios::beg);
  	fs.read((char *)buffer,sizeof(buffer));
  	curdata = fs.tellg();
	// read header information and save some of them into each vars.
	// Total size of header is 264 bytes
	bstr=""; hbegin=0; hend=4; // HEAD
	for(int i=hbegin;i<hend;i++) bstr+=buffer[i];
    strcpy(head, bstr.c_str());
	if(!strcmp(head,"HEAD")){
	  bstr=""; hbegin=8; hend=16; // HHIRF
	  bstr=""; hbegin=16; hend=24; // L003
	  bstr=""; hbegin=24; hend=40; // LIST DATA
	  bstr=""; hbegin=40; hend=56; // 04/01/11 09:25  
	  for(int i=hbegin;i<hend;i++) bstr+=buffer[i];
      strcpy(rundate, bstr.c_str());
	  bstr=""; hbegin=56; hend=136; // htit
	  for(int i=hbegin;i<hend;i++) bstr+=buffer[i];
      strcpy(runtitle, bstr.c_str());
	  bstr=""; hbegin=136; hend=140; // hnum (Run No.)
  	  runno = buffer[hbegin]&0xff;
	  cout << "Run # is " << runno << "." <<endl;
	  char *ptr;
      for (ptr=runtitle+strlen(runtitle)-1;(ptr >= runtitle)&&isspace(*ptr);--ptr);
      ptr[1] = '\0';

	}else{
	  cout << "Header format is incorrect. Run No. was not stored!" << endl;
	}

	//return get pointer to beginning of file
	fs.seekg(0,ios::beg);
	//move pointer past first 65552 bytes (the file header) to start of first block
	fs.seekg(header,ios::beg);
	curdata = fs.tellg();	
	
	
	// Blocks to read between
	if (partial){
		start = atoi(argv[4]);
		stop = atoi(argv[5]);
		}
	
	//********** Start reading data from file **********// 
	
	
	// Enter loop over each 32776-byte block
	while (block < nblock){	
		//read in buffer: 32776 8-bit words; 16388 16-bit (2-byte) words 
		fs.read((char *)word,sizeof(word));
		curdata = fs.tellg();
		
		// If only certain blocks are requested skip over others
		if (partial){
			if (block < start-1 || block > stop-1){
				block+=1;
				continue;
				}
			}
		
		nword = 4; // Skip first four words (8 bytes)
		
		// Enter loop over words in block
		while (nword < (blocksize/2)){		
			// If there are two 'ffff' blocks in a row, have hit padding -> exit block
			if (word[nword]==0xffff){
				if (word[nword+2]==0xffff){
					break;
					}
				else {
					nword+=2;
					}
				}			
				
			// First word can be any detector, must start with 8
			if (int(word[nword]&0xf000)!=32768){    //32768 = 0x8000
				break;
				}
			
			// Enter loop over event -> end when buffer is hit
			// Data in detector-value pairs, one word (two bytes) each
			multiplicity = 0;
			while (word[nword]!=0xffff){
				adc = word[nword]&0xfff;
				value = word[nword+1]&0x3fff;
				bpacid[multiplicity]=adc;
				bvalue[multiplicity]=value;
                                counter[adc]++;
				multiplicity++;
				nword+=2;
				}
			
			// Add the data to tree
			if(multiplicity>0) {
                           if(iscid) multiplicity-=3;
                           tt->Fill();
                        }
            //progressbar(100*block/nblock);	
			// Record event
			nevt +=1;
			// Skip to next event 
			nword+=2;
			}
		
		tevt = tevt+nevt;
		nevt=0;
		block +=1;
        progressbar(100*block/nblock);	
		}
	
	tt->Branch("RunNo",&runno,"RunNo/I");
	tt->Branch("FileName",&filename,"FileName/C");
	tt->Branch("RunTitle",&runtitle,"RunTitle/C");
	tt->Branch("RunDate",&rundate,"RunDate/C");
	tt->Fill();
	ff->cd();
	//tt->Print();
	tt->Write();
	ff->Close();
	fs.close();
	
        cout << endl;
	cout << "# of events in ADC channel:" << endl;
        for(int i=0;i<maxpacid+1;i++) cout << "adc_" << i << " : " << counter[i] << endl;
	cout << tevt << " events in total" << endl;
        result=makerootcode(rootname, treename, maxpacid);

	return 0;
}

//  this function prints progress bar
//  pre: percentage complete
//  post: progress bar
void progressbar(int percent){
    static int x = 0;
    string slash[4];
    slash[0] = "\\";
    slash[1] = "-";
    slash[2] = "/";
    slash[3] = "|";
    cout << "\r"; // carriage return back to beginning of line
    cout << slash[x] << " " << percent << " %"; // print the bars and percentage
    x++; // increment to make the slash appear to rotate
    if(x == 4)
    x = 0; // reset slash animation
}

int makerootcode(string rootname, string treename, int maxpacid){
  ifstream InData;
  string line;
  ofstream OutData;
  OutData.open("readroot.C");
  char* cwd;
  char buff[PATH_MAX + 1];
  cwd = getcwd( buff, PATH_MAX + 1 );
  if( cwd != NULL ) {
    printf("Make sure readroot.1-6 files should be placed in %s.\n", cwd);
    printf("Or, the readroot.C will be mostly empty.\n");
  }
  InData.open("readroot.1");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData << Form("  TH1I *hSingle[%d];",maxpacid+1)<<endl;
  InData.open("readroot.2");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData << Form("  int multiplicity,bpacid[%d],bvalue[%d];",maxpacid+1,maxpacid+1)<<endl;
  InData.open("readroot.3");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData << Form("  rootname=\"%s\";",rootname.c_str())<<endl;
  OutData << Form("  treename=\"%s\";",treename.c_str())<<endl;
  InData.open("readroot.4");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData << Form("  for(int i=0;i<%d;i++){",maxpacid+1)<<endl;
  InData.open("readroot.5");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData << Form("      if(adc>=0 && adc<%d) {",maxpacid+1)<<endl;
  InData.open("readroot.6");
  while(InData.good()){
    getline(InData,line);
    OutData << line << endl;
  }
  InData.close();
  OutData.close();
 
  return 0;
}
